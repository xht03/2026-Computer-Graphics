"""Geometric Verifier Agent — Day 4 implementation.

Architecture:
  - build_verify_snippet() returns a bpy code snippet injected into the
    Blender run (BETWEEN model code and render code).
  - parse_from_stdout() extracts the JSON issue list from captured stdout.

Five checks performed inside Blender:
  1. Manifold integrity + isolated vertices
  2. Bounding-box overlap / gap detection between parts
  3. Overall proportion constraints (scale, H/W ratio)
  4. Empty mesh guard (zero-vertex objects)
  5. Centre-of-mass stability (vertex-centroid proxy)
"""

import json
import re

_SENTINEL_START = "GEOM_VERIFY_START"
_SENTINEL_END = "GEOM_VERIFY_END"

# -------------------------------------------------------------------
# The snippet runs INSIDE Blender, so it may use bpy / bmesh / mathutils.
# All temporary variable names are prefixed with _G_ to avoid collisions
# with the model code above it.
# -------------------------------------------------------------------
_GEOM_VERIFY_SNIPPET = """\
# ── Programmatic Geometry Verification (injected by GeomVerifier) ───────────
import bpy as _bpy, bmesh as _bmesh, json as _json, math as _math
import mathutils as _mu

_G_ISSUES = []


def _G_bbox_world(obj):
    c = [obj.matrix_world @ _mu.Vector(v) for v in obj.bound_box]
    return (min(v.x for v in c), min(v.y for v in c), min(v.z for v in c),
            max(v.x for v in c), max(v.y for v in c), max(v.z for v in c))


_G_mesh_objs = [o for o in _bpy.data.objects if o.type == 'MESH']

if not _G_mesh_objs:
    _G_ISSUES.append({
        "source": "geom", "severity": "error",
        "message": "No mesh objects found in scene — nothing was generated.",
        "component": "scene",
    })
else:
    # ── Check 1: Manifold integrity + isolated vertices ──────────────────
    for _G_obj in _G_mesh_objs:
        _G_bm = _bmesh.new()
        _G_bm.from_mesh(_G_obj.data)
        _G_nm = sum(1 for e in _G_bm.edges if not e.is_manifold)
        _G_iv = sum(1 for v in _G_bm.verts if not v.link_edges)
        _G_bm.free()

        if _G_nm > 0:
            _G_ISSUES.append({
                "source": "geom", "severity": "warning",
                "message": (f"'{_G_obj.name}' has {_G_nm} non-manifold edge(s). "
                            "This prevents clean 3D-printing or physics simulation."),
                "component": _G_obj.name,
            })
        if _G_iv > 0:
            _G_ISSUES.append({
                "source": "geom", "severity": "warning",
                "message": (f"'{_G_obj.name}' has {_G_iv} isolated vertex/vertices "
                            "(dangling geometry with no edges)."),
                "component": _G_obj.name,
            })

        # ── Check 4: Empty mesh guard ────────────────────────────────────
        if len(_G_obj.data.vertices) == 0:
            _G_ISSUES.append({
                "source": "geom", "severity": "error",
                "message": f"'{_G_obj.name}' is an empty mesh (zero vertices).",
                "component": _G_obj.name,
            })

    # ── Check 2: Bounding-box overlap / gap between parts ────────────────
    _G_bboxes = {o.name: _G_bbox_world(o) for o in _G_mesh_objs}
    _G_names  = list(_G_bboxes.keys())
    _G_gap_candidates = []  # (gap_m, name1, name2) for deferred reporting
    for _G_i in range(len(_G_names)):
        for _G_j in range(_G_i + 1, len(_G_names)):
            _G_n1, _G_n2 = _G_names[_G_i], _G_names[_G_j]
            _G_b1, _G_b2 = _G_bboxes[_G_n1], _G_bboxes[_G_n2]

            # Overlap: positive intersection volume in all three axes
            _G_ox = max(0.0, min(_G_b1[3], _G_b2[3]) - max(_G_b1[0], _G_b2[0]))
            _G_oy = max(0.0, min(_G_b1[4], _G_b2[4]) - max(_G_b1[1], _G_b2[1]))
            _G_oz = max(0.0, min(_G_b1[5], _G_b2[5]) - max(_G_b1[2], _G_b2[2]))
            if _G_ox > 0.01 and _G_oy > 0.01 and _G_oz > 0.01:
                _G_vol = _G_ox * _G_oy * _G_oz
                if _G_vol > 5e-4:  # ignore < 0.5 cm³ floating-point tolerance
                    _G_ISSUES.append({
                        "source": "geom", "severity": "warning",
                        "message": (
                            f"'{_G_n1}' and '{_G_n2}' interpenetrate "
                            f"(overlap ≈ {_G_vol * 1e6:.1f} cm³). "
                            "Parts are colliding — check positions/scales."
                        ),
                        "component": f"{_G_n1},{_G_n2}",
                    })
            else:
                # Gap: smallest 3-D separation between the two AABBs
                _G_dx = max(0.0, max(_G_b1[0], _G_b2[0]) - min(_G_b1[3], _G_b2[3]))
                _G_dy = max(0.0, max(_G_b1[1], _G_b2[1]) - min(_G_b1[4], _G_b2[4]))
                _G_dz = max(0.0, max(_G_b1[2], _G_b2[2]) - min(_G_b1[5], _G_b2[5]))
                _G_gap = _math.sqrt(_G_dx * _G_dx + _G_dy * _G_dy + _G_dz * _G_dz)
                # Flag suspicious gaps: 1 cm – 12 cm.
                # Upper limit 12 cm prevents false positives: parts > 12 cm apart
                # were never meant to be joined (e.g. legs across a table's width).
                if 0.01 < _G_gap < 0.12:
                    _G_gap_candidates.append((_G_gap, _G_n1, _G_n2))

    # Report only the top-3 smallest gaps (most suspicious) to keep issue count low.
    _G_gap_candidates.sort()
    for _G_gap, _G_n1, _G_n2 in _G_gap_candidates[:3]:
        _G_ISSUES.append({
            "source": "geom", "severity": "info",
            "message": (
                f"Gap between '{_G_n1}' and '{_G_n2}': "
                f"{_G_gap * 100:.1f} cm. "
                "Parts may not be properly joined — VLM cannot detect this."
            ),
            "component": f"{_G_n1},{_G_n2}",
        })

    # ── Check 3: Overall proportion constraints ───────────────────────────
    _G_all_c = []
    for _G_obj in _G_mesh_objs:
        _G_all_c.extend([_G_obj.matrix_world @ _mu.Vector(v) for v in _G_obj.bound_box])
    _G_Xmin = min(c.x for c in _G_all_c); _G_Xmax = max(c.x for c in _G_all_c)
    _G_Ymin = min(c.y for c in _G_all_c); _G_Ymax = max(c.y for c in _G_all_c)
    _G_Zmin = min(c.z for c in _G_all_c); _G_Zmax = max(c.z for c in _G_all_c)
    _G_W = _G_Xmax - _G_Xmin
    _G_D = _G_Ymax - _G_Ymin
    _G_H = _G_Zmax - _G_Zmin

    for _G_dname, _G_dval in [("width", _G_W), ("depth", _G_D), ("height", _G_H)]:
        if _G_dval > 0 and (_G_dval < 0.05 or _G_dval > 15.0):
            _G_ISSUES.append({
                "source": "geom", "severity": "warning",
                "message": (
                    f"Overall {_G_dname} is {_G_dval:.3f} m — outside typical "
                    "furniture range [0.05, 15.0] m. Check scale."
                ),
                "component": "scene",
            })

    if _G_H > 0 and _G_W > 0 and _G_H / _G_W > 6.0:
        _G_ISSUES.append({
            "source": "geom", "severity": "warning",
            "message": (
                f"Height/width ratio {_G_H / _G_W:.2f} is very large. "
                "Object may be unrealistically tall or narrow."
            ),
            "component": "scene",
        })

    # ── Check 5: Centre-of-mass stability (vertex centroid proxy) ────────
    _G_tx = _G_ty = 0.0
    _G_n = 0
    for _G_obj in _G_mesh_objs:
        for _G_v in _G_obj.data.vertices:
            _G_wv = _G_obj.matrix_world @ _G_v.co
            _G_tx += _G_wv.x
            _G_ty += _G_wv.y
            _G_n  += 1
    if _G_n > 0:
        _G_cx = _G_tx / _G_n
        _G_cy = _G_ty / _G_n
        _G_margin = max(_G_W, _G_D) * 0.15
        if not (_G_Xmin - _G_margin <= _G_cx <= _G_Xmax + _G_margin and
                _G_Ymin - _G_margin <= _G_cy <= _G_Ymax + _G_margin):
            _G_ISSUES.append({
                "source": "geom", "severity": "info",
                "message": (
                    f"Estimated centre-of-mass ({_G_cx:.2f}, {_G_cy:.2f}) "
                    "is outside the footprint. Object may be structurally unstable."
                ),
                "component": "scene",
            })

    # ── Check 6: Leg-to-tabletop connectivity ────────────────────────────
    # Find the candidate "top surface" = object with largest XY footprint.
    # Then verify that all narrow objects (legs) below it actually touch it.
    _G_xy_areas = {
        o.name: (_G_bboxes[o.name][3] - _G_bboxes[o.name][0]) *
                (_G_bboxes[o.name][4] - _G_bboxes[o.name][1])
        for o in _G_mesh_objs
    }
    if len(_G_xy_areas) >= 2:
        _G_top_name = max(_G_xy_areas, key=lambda n: _G_xy_areas[n])
        _G_top_bbox = _G_bboxes[_G_top_name]
        _G_top_area = _G_xy_areas[_G_top_name]
        _G_top_z_bottom = _G_top_bbox[2]   # min-z of the top surface = its bottom face

        for _G_n, _G_b in _G_bboxes.items():
            if _G_n == _G_top_name:
                continue
            _G_part_area = _G_xy_areas[_G_n]
            # Only examine narrow objects (legs/posts): area < 15% of top surface
            if _G_part_area >= _G_top_area * 0.15:
                continue
            # Must overlap in XY with the top surface (i.e. be "under" it)
            _G_xy_ox = max(0.0, min(_G_b[3], _G_top_bbox[3]) - max(_G_b[0], _G_top_bbox[0]))
            _G_xy_oy = max(0.0, min(_G_b[4], _G_top_bbox[4]) - max(_G_b[1], _G_top_bbox[1]))
            if _G_xy_ox <= 0 or _G_xy_oy <= 0:
                continue
            # Must be below the top surface
            _G_part_z_top = _G_b[5]
            if _G_part_z_top >= _G_top_z_bottom + 0.001:
                continue  # object is above the top surface (not a leg)
            _G_v_gap = _G_top_z_bottom - _G_part_z_top
            if _G_v_gap > 0.005:  # gap > 5 mm is a bug
                _G_target_leg_h = _G_top_z_bottom
                _G_ISSUES.append({
                    "source": "geom", "severity": "warning",
                    "message": (
                        f"Leg/post '{_G_n}' top is at z={_G_part_z_top:.4f} m but "
                        f"'{_G_top_name}' bottom is at z={_G_top_z_bottom:.4f} m "
                        f"(gap = {_G_v_gap * 100:.1f} cm). "
                        f"Fix: set leg_h = {_G_target_leg_h:.4f} and "
                        f"leg.location.z = {_G_target_leg_h / 2:.4f}, "
                        f"leg.scale.z = {_G_target_leg_h / 2:.4f}."
                    ),
                    "component": f"{_G_n},{_G_top_name}",
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
