"""Geometric Verifier Agent — semantic-connectivity rewrite.

Design (L3GO-style, see plan): instead of all-pairs AABB guessing (which
produced false positives that misled the Fixer), this verifier only reports
ZERO-false-positive HARD ERRORS, all at severity "error":

  E1. Empty scene / empty mesh         — nothing was generated.
  E2. Degenerate or explosive scale    — a part ≈0 on an axis, or the whole
                                          model outside [0.05 m, 15 m].
  E3. Declared connectivity broken     — the model code declares a module-level
                                          CONNECTIONS list of (part_a, part_b)
                                          pairs that MUST touch; we check only
                                          those pairs. A declared pair separated
                                          by > 1 cm is a real assembly bug.

If the model code does not define CONNECTIONS, E3 is skipped entirely (graceful
degradation — never guess). The old noisy checks (manifold, all-pairs gap/
overlap, centre-of-mass, area-heuristic leg detection) are intentionally gone.

Architecture (unchanged):
  - build_verify_snippet() returns a bpy snippet injected BETWEEN model code
    and render code, so CONNECTIONS (defined in the model code) is visible.
  - parse_from_stdout() extracts the JSON issue list from captured stdout.
"""

import json
import re

_SENTINEL_START = "GEOM_VERIFY_START"
_SENTINEL_END = "GEOM_VERIFY_END"

# -------------------------------------------------------------------
# The snippet runs INSIDE Blender, so it may use bpy / mathutils.
# All temporary variable names are prefixed with _G_ to avoid collisions
# with the model code above it.  CONNECTIONS, if defined by the model code,
# lives in the same module globals and is read via globals().get(...).
# -------------------------------------------------------------------
_GEOM_VERIFY_SNIPPET = """\
# ── Programmatic Geometry Verification (injected by GeomVerifier) ───────────
import bpy as _bpy, json as _json
import mathutils as _mu

# CRITICAL: the model code sets obj.scale/obj.location by direct assignment,
# which does NOT refresh obj.matrix_world until the depsgraph updates. Without
# this call every bbox below would be read from the UNSCALED unit cube, making
# all geometry checks wrong. Force the transforms to settle first.
_bpy.context.view_layer.update()

_G_ISSUES = []


def _G_bbox_world(obj):
    c = [obj.matrix_world @ _mu.Vector(v) for v in obj.bound_box]
    return (min(v.x for v in c), min(v.y for v in c), min(v.z for v in c),
            max(v.x for v in c), max(v.y for v in c), max(v.z for v in c))


def _G_aabb_gap(b1, b2):
    # Smallest 3-D separation between two AABBs (0 if they touch/overlap).
    _dx = max(0.0, max(b1[0], b2[0]) - min(b1[3], b2[3]))
    _dy = max(0.0, max(b1[1], b2[1]) - min(b1[4], b2[4]))
    _dz = max(0.0, max(b1[2], b2[2]) - min(b1[5], b2[5]))
    return (_dx * _dx + _dy * _dy + _dz * _dz) ** 0.5


_G_mesh_objs = [o for o in _bpy.data.objects if o.type == 'MESH']

# ── E1: Empty scene / empty mesh ─────────────────────────────────────────
if not _G_mesh_objs:
    _G_ISSUES.append({
        "source": "geom", "severity": "error",
        "message": "No mesh objects found in scene — nothing was generated.",
        "component": "scene",
    })
else:
    for _G_obj in _G_mesh_objs:
        if len(_G_obj.data.vertices) == 0:
            _G_ISSUES.append({
                "source": "geom", "severity": "error",
                "message": f"'{_G_obj.name}' is an empty mesh (zero vertices).",
                "component": _G_obj.name,
            })

    _G_bboxes = {o.name: _G_bbox_world(o) for o in _G_mesh_objs}

    # ── E2: Degenerate / explosive scale ────────────────────────────────
    # Per-part: any axis essentially zero (≈ flat/degenerate part).
    for _G_name, _G_b in _G_bboxes.items():
        _G_dims = (_G_b[3] - _G_b[0], _G_b[4] - _G_b[1], _G_b[5] - _G_b[2])
        if min(_G_dims) < 1e-4:
            _G_ISSUES.append({
                "source": "geom", "severity": "error",
                "message": (
                    f"'{_G_name}' is degenerate (size "
                    f"{_G_dims[0]:.4f}x{_G_dims[1]:.4f}x{_G_dims[2]:.4f} m) — "
                    "one axis is ~0. Check its scale."
                ),
                "component": _G_name,
            })

    # Whole-model: overall bounding box far outside furniture range.
    _G_Xmin = min(b[0] for b in _G_bboxes.values())
    _G_Ymin = min(b[1] for b in _G_bboxes.values())
    _G_Zmin = min(b[2] for b in _G_bboxes.values())
    _G_Xmax = max(b[3] for b in _G_bboxes.values())
    _G_Ymax = max(b[4] for b in _G_bboxes.values())
    _G_Zmax = max(b[5] for b in _G_bboxes.values())
    _G_dims_all = {"width": _G_Xmax - _G_Xmin,
                   "depth": _G_Ymax - _G_Ymin,
                   "height": _G_Zmax - _G_Zmin}
    for _G_dname, _G_dval in _G_dims_all.items():
        if _G_dval > 0 and (_G_dval < 0.05 or _G_dval > 15.0):
            _G_ISSUES.append({
                "source": "geom", "severity": "error",
                "message": (
                    f"Overall {_G_dname} is {_G_dval:.3f} m — outside the valid "
                    "furniture range [0.05, 15.0] m. The whole model is mis-scaled."
                ),
                "component": "scene",
            })

    # ── E3: Declared connectivity ───────────────────────────────────────
    # Read CONNECTIONS from the model code's globals (if it declared any).
    _G_conns = globals().get("CONNECTIONS", None)
    if isinstance(_G_conns, (list, tuple)):
        def _G_match(spec):
            # Resolve a CONNECTIONS name to actual object bboxes:
            # exact name first, else prefix match (e.g. "leg" -> leg_1..leg_4).
            if spec in _G_bboxes:
                return [spec]
            _hits = [n for n in _G_bboxes if n.startswith(spec)]
            return _hits

        _G_TOL = 0.01  # 1 cm: declared-to-touch parts farther than this = bug
        _G_reported = set()
        for _G_pair in _G_conns:
            if not (isinstance(_G_pair, (list, tuple)) and len(_G_pair) == 2):
                continue
            _G_a_names = _G_match(str(_G_pair[0]))
            _G_b_names = _G_match(str(_G_pair[1]))
            if not _G_a_names or not _G_b_names:
                continue  # name not found — skip rather than guess
            for _G_an in _G_a_names:
                # A declared part is satisfied if it touches ANY resolved
                # counterpart (e.g. leg_1 only needs to touch the tabletop).
                _G_min_gap = min(_G_aabb_gap(_G_bboxes[_G_an], _G_bboxes[_G_bn])
                                 for _G_bn in _G_b_names if _G_bn != _G_an)
                if _G_min_gap > _G_TOL:
                    _G_key = tuple(sorted((_G_an, str(_G_pair[1]))))
                    if _G_key in _G_reported:
                        continue
                    _G_reported.add(_G_key)
                    _G_ISSUES.append({
                        "source": "geom", "severity": "error",
                        "message": (
                            f"Declared connection '{_G_an}' <-> '{_G_pair[1]}' is "
                            f"broken: parts are {_G_min_gap * 100:.1f} cm apart but "
                            "were declared to touch. Reposition/resize so they meet."
                        ),
                        "component": f"{_G_an},{_G_pair[1]}",
                    })

print("GEOM_VERIFY_START")
print(_json.dumps(_G_ISSUES, ensure_ascii=False))
print("GEOM_VERIFY_END")
"""


class GeomVerifier:
    """Programmatic geometry checker (runs inside Blender as injected code)."""

    @staticmethod
    def build_verify_snippet() -> str:
        """Return bpy snippet to inject BEFORE the render snippet.

        Checks run after the model code so all mesh objects are in the scene.
        Results are printed between GEOM_VERIFY_START / GEOM_VERIFY_END sentinels.
        """
        return _GEOM_VERIFY_SNIPPET

    @staticmethod
    def parse_from_stdout(stdout: str) -> list[dict]:
        """Extract and parse the geometry-issue JSON from Blender's captured stdout."""
        m = re.search(
            r"GEOM_VERIFY_START\n(.*?)\nGEOM_VERIFY_END",
            stdout,
            re.DOTALL,
        )
        if not m:
            return []
        try:
            data = json.loads(m.group(1).strip())
            if isinstance(data, list):
                return data
        except json.JSONDecodeError:
            pass
        return []

    def verify(self, blend_file: str | None = None, mesh_data: dict | None = None) -> list[dict]:
        """Not used in the integrated pipeline.

        Geom verification is injected into every Blender run via
        build_verify_snippet() + parse_from_stdout().  This method is kept
        for the interface contract only.
        """
        raise NotImplementedError(
            "Use build_verify_snippet() / parse_from_stdout() instead. "
            "Geom checks are injected into the main Blender pipeline run."
        )
