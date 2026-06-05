#!/usr/bin/env python3
"""
Day 4 milestone test: Programmatic Geometry Verifier.

Tests both the integrated pipeline AND an ablation comparison:
  - VLM-only critique (use_geom=False)
  - Geometry-verifier-only critique (use_vlm=False, use_geom=True)
  - Combined (use_vlm=True, use_geom=True)

The key goal for Day 4 is proving that the geometry verifier catches issues
that the VLM cannot see from rendered images (floating parts, interpenetration,
scale anomalies).

Usage:
    cd pj3_blender_agent
    python test_day4.py                    # run all standard cases
    python test_day4.py --case 方桌 --geom-only
    python test_day4.py --ablation         # VLM-only vs geom vs combined
"""
import argparse
import json
import os
import sys
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()
load_dotenv(".env.example")

from main import pipeline, _blender_run, _save_code, _make_run_dir
from agents.geom_verifier import GeomVerifier

# ── Standard test cases ────────────────────────────────────────────────────
TEST_CASES = ["长凳", "方桌", "屏风"]

# ── Intentionally broken scripts for ablation ─────────────────────────────
# These scripts have geometry bugs that LOOK fine in renders but are caught
# by the verifier.  Use them to demonstrate "VLM misses, α catches."

_BROKEN_FLOATING_LEG = """\
import bpy, math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("Wood")
mat.use_nodes = True
bpy.data.materials["Wood"].node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)

# COMPONENT: tabletop
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.77))
top = bpy.context.active_object
top.name = "tabletop"
top.scale = (0.45, 0.45, 0.025)
top.data.materials.append(mat)

# COMPONENT: leg_1  (deliberately 5 cm above the ground — floating!)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.35, 0.35, 0.42))
leg1 = bpy.context.active_object
leg1.name = "leg_1"
leg1.scale = (0.025, 0.025, 0.22)
leg1.data.materials.append(mat)

# COMPONENT: leg_2
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.35, -0.35, 0.22))
leg2 = bpy.context.active_object
leg2.name = "leg_2"
leg2.scale = (0.025, 0.025, 0.22)
leg2.data.materials.append(mat)

# COMPONENT: leg_3
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.35, 0.35, 0.22))
leg3 = bpy.context.active_object
leg3.name = "leg_3"
leg3.scale = (0.025, 0.025, 0.22)
leg3.data.materials.append(mat)

# COMPONENT: leg_4
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.35, -0.35, 0.22))
leg4 = bpy.context.active_object
leg4.name = "leg_4"
leg4.scale = (0.025, 0.025, 0.22)
leg4.data.materials.append(mat)
"""

_BROKEN_INTERPENETRATING_LEGS = """\
import bpy

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("Wood")
mat.use_nodes = True
bpy.data.materials["Wood"].node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)

# COMPONENT: tabletop
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.77))
top = bpy.context.active_object
top.name = "tabletop"
top.scale = (0.45, 0.45, 0.025)
top.data.materials.append(mat)

# COMPONENT: leg_1 + leg_2 deliberately overlap (same position!)
for i, loc in enumerate([(0.35, 0.35, 0.365), (0.33, 0.33, 0.365),
                          (-0.35, 0.35, 0.365), (-0.35, -0.35, 0.365)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=loc)
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (0.025, 0.025, 0.365)
    leg.data.materials.append(mat)
"""

_BROKEN_WRONG_SCALE = """\
import bpy

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("Wood")
mat.use_nodes = True
bpy.data.materials["Wood"].node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)

# COMPONENT: tabletop — scale in mm instead of m (100x too large!)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 77))
top = bpy.context.active_object
top.name = "tabletop"
top.scale = (45, 45, 2.5)
top.data.materials.append(mat)

for i, (lx, ly) in enumerate([(35, 35), (35, -35), (-35, 35), (-35, -35)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(lx, ly, 36.5))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (2.5, 2.5, 36.5)
    leg.data.materials.append(mat)
"""

_ABLATION_CASES = [
    ("floating_leg",      _BROKEN_FLOATING_LEG,          "方桌（悬浮腿）"),
    ("interpenetrating",  _BROKEN_INTERPENETRATING_LEGS,  "方桌（腿部穿模）"),
    ("wrong_scale",       _BROKEN_WRONG_SCALE,            "方桌（尺度错误）"),
]


# ── Helpers ────────────────────────────────────────────────────────────────

def _run_geom_only(code: str, run_dir: str, label: str) -> list[dict]:
    """Run code in Blender with only geom verifier (no VLM)."""
    import time
    iter_dir = os.path.join(run_dir, "geom_check")
    success, result, geom_issues = _blender_run(code, iter_dir, use_geom=True)
    return geom_issues


def _run_vlm_only(code: str, run_dir: str, description: str) -> list[dict]:
    """Run code in Blender then run VLM critic on renders (no geom)."""
    iter_dir = os.path.join(run_dir, "vlm_check")
    success, result, _ = _blender_run(code, iter_dir, use_geom=False)
    if not success:
        return [{"source": "vlm", "severity": "error", "message": "Blender run failed."}]
    from agents.vlm_critic import VLMCritic
    critic = VLMCritic()
    renders = [os.path.join(iter_dir, f"{v}.png") for v in ("front", "side", "top", "iso")]
    return critic.critique(renders, description)


# ── Standard pipeline test ─────────────────────────────────────────────────

def run_standard(cases: list[str], iterations: int, use_vlm: bool) -> None:
    results = []
    for desc in cases:
        print(f"\n{'#'*60}\n# {desc}  (geom=True, vlm={use_vlm}, iters={iterations})\n{'#'*60}")
        result = pipeline(
            description=desc,
            use_planner=True,
            use_rag=False,
            max_iterations=iterations,
            use_vlm=use_vlm,
            use_geom=True,
        )
        total_geom = sum(
            sum(1 for i in rnd if i.get("source") == "geom")
            for rnd in result.get("issues_per_round", [])
        )
        results.append((desc, result["success"], result["run_dir"], total_geom))

    print(f"\n{'='*60}\n  Day 4 Standard Test Summary\n{'='*60}")
    for desc, ok, run_dir, n_geom in results:
        status = "✓ PASS" if ok else "✗ FAIL"
        print(f"  {status}  {desc:8s}  geom_issues_total={n_geom}  → {run_dir}")
    passed = sum(1 for _, ok, _, _ in results if ok)
    print(f"\n  {passed}/{len(results)} passed\n{'='*60}\n")
    sys.exit(0 if passed == len(results) else 1)


# ── Ablation: VLM-only vs Geom-only vs Combined ───────────────────────────

def run_ablation() -> None:
    """Demonstrate that geom verifier catches issues VLM cannot see."""
    print(f"\n{'='*60}")
    print("  Day 4 Ablation: VLM-only vs Geometry Verifier")
    print(f"{'='*60}\n")

    rows = []
    for case_id, code, label in _ABLATION_CASES:
        run_dir = _make_run_dir(f"ablation_{case_id}")
        os.makedirs(run_dir, exist_ok=True)

        print(f"\n── Case: {label} ──")
        _save_code(code, os.path.join(run_dir, "code_v0.py"))

        # Geom-only check
        print("  [Geom verifier]...")
        geom_issues = _run_geom_only(code, run_dir, label)
        print(f"    → {len(geom_issues)} issue(s) found")
        for iss in geom_issues:
            print(f"      [{iss['severity']}] {iss['message'][:80]}")

        # VLM-only check
        print("  [VLM critic]...")
        try:
            vlm_issues = _run_vlm_only(code, run_dir, label)
            print(f"    → {len(vlm_issues)} issue(s) found")
            for iss in vlm_issues[:3]:
                print(f"      [{iss.get('severity','?')}] {iss.get('message','')[:80]}")
        except Exception as e:
            vlm_issues = []
            print(f"    → VLM skipped ({e})")

        rows.append({
            "case": label,
            "geom_issues": len(geom_issues),
            "vlm_issues":  len(vlm_issues),
            "geom_catches_vlm_misses": len(geom_issues) > 0 and len(vlm_issues) == 0,
        })

    print(f"\n{'='*60}\n  Ablation Summary\n{'='*60}")
    print(f"  {'Case':<25} {'Geom':>6} {'VLM':>6} {'Geom-only catch?':>18}")
    print(f"  {'-'*55}")
    for r in rows:
        flag = "★ YES" if r["geom_catches_vlm_misses"] else "no"
        print(f"  {r['case']:<25} {r['geom_issues']:>6} {r['vlm_issues']:>6} {flag:>18}")
    print(f"\n  ★ marks cases where geom catches issues VLM misses.")
    print(f"{'='*60}\n")

    # Save ablation results for the paper
    ablation_path = os.path.join("outputs", "ablation_day4.json")
    with open(ablation_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print(f"  Ablation results saved → {ablation_path}\n")


# ── Entry point ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Day 4 Geometry Verifier Tests")
    parser.add_argument("--case", help="Single standard case to test")
    parser.add_argument("--iterations", type=int, default=2)
    parser.add_argument("--vlm", action="store_true", help="Also enable VLM critic")
    parser.add_argument("--geom-only", action="store_true",
                        help="Run standard cases with geom only (no VLM)")
    parser.add_argument("--ablation", action="store_true",
                        help="Run ablation comparison (VLM-only vs geom-only)")
    args = parser.parse_args()

    if args.ablation:
        run_ablation()
    else:
        cases = [args.case] if args.case else TEST_CASES
        use_vlm = args.vlm and not args.geom_only
        run_standard(cases, args.iterations, use_vlm)
