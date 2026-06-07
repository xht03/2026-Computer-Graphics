#!/usr/bin/env python3
"""
Main pipeline entry point.

Usage:
    python main.py "a red cube"
    python main.py "做一张明式案几" --iterations 3 --vlm --geom
    python main.py "长凳" --iterations 2
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()
load_dotenv(".env.example")

from agents.coder import CoderAgent
from agents.planner import PlannerAgent
from agents.fixer import FixerAgent
from agents.api_rag import ApiDocRetriever
from blender.render import build_render_script
from blender.runner import run_script

OUTPUT_ROOT = os.getenv("OUTPUT_ROOT", "outputs")


def _make_run_dir(description: str) -> str:
    ts = time.strftime("%Y%m%d_%H%M%S")
    slug = description[:20].replace(" ", "_").replace("/", "-")
    run_dir = os.path.join(OUTPUT_ROOT, f"{ts}_{slug}")
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    return run_dir


def _save_code(code: str, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(code)


def _blender_run(
    code: str,
    output_dir: str,
    use_geom: bool = False,
    verbose: bool = False,
) -> tuple[bool, dict, list]:
    """Append geometry-check + render snippets and run Blender headlessly.

    Returns:
        (success, result_dict, geom_issues)

    Geometry issues are detected by injecting a bpy verification snippet
    BETWEEN the model code and the render snippet.  The snippet prints JSON
    to stdout; we parse it after the process exits.  This catches problems
    (floating parts, interpenetration, scale errors) that VLM critique misses.
    """
    abs_output_dir = str(Path(output_dir).resolve())

    # Inject geometry verification BEFORE render so the clean scene is checked
    geom_snippet = ""
    if use_geom:
        from agents.geom_verifier import GeomVerifier
        geom_snippet = "\n\n" + GeomVerifier.build_verify_snippet()

    render_snippet = build_render_script(abs_output_dir)
    full_script = code + geom_snippet + "\n\n" + render_snippet
    result = run_script(full_script)

    # Parse geometry issues from stdout
    geom_issues: list[dict] = []
    if use_geom:
        from agents.geom_verifier import GeomVerifier
        geom_issues = GeomVerifier.parse_from_stdout(result["stdout"])
        sev_counts = {}
        for iss in geom_issues:
            s = iss.get("severity", "info")
            sev_counts[s] = sev_counts.get(s, 0) + 1
        summary = ", ".join(f"{v} {k}" for k, v in sorted(sev_counts.items())) if sev_counts else "none"
        print(f"  [Geom] {len(geom_issues)} issue(s): {summary}")
        if verbose:
            if geom_issues:
                for iss in geom_issues:
                    sev = iss.get('severity','?').upper()
                    msg = iss.get('message','')
                    comp = iss.get('component','')
                    comp_str = f" [{comp}]" if comp else ""
                    print(f"    │  [geom] {sev}{comp_str}: {msg}")
            else:
                print(f"    │  (geom verifier found no issues — geometry looks structurally correct)")

    # Show key output lines
    for line in result["stdout"].splitlines():
        if any(kw in line for kw in ("RENDER_OK", "ALL_RENDERS", "Error", "Traceback")):
            print(f"    [blender] {line}")
    if result.get("error_detail"):
        print(f"  [error]\n{result['error_detail']}")

    return result["success"], result, geom_issues


def pipeline(
    description: str,
    use_planner: bool = True,
    max_iterations: int = 1,
    use_vlm: bool = False,
    use_geom: bool = False,
    verbose: bool = False,
) -> dict:
    """Run the full agent pipeline.

    Returns a summary dict with keys:
        success, run_dir, iterations, issues_per_round, plan
    """
    run_dir = _make_run_dir(description)
    print(f"\n{'='*60}")
    print(f"  PJ3 Blender Agent")
    print(f"  Description : {description}")
    print(f"  Output dir  : {run_dir}")
    print(f"  VLM critic  : {use_vlm}  |  Geom verifier: {use_geom}")
    print(f"  Max iters   : {max_iterations}")
    print(f"{'='*60}\n")

    # ── Step 1: Plan ──────────────────────────────────────────────
    coder_desc = description
    plan = {}
    if use_planner:
        print("[1/N] Planning...")
        planner = PlannerAgent()
        plan = planner.plan(description)
        coder_desc = planner.plan_to_coder_description(plan)
        plan_path = os.path.join(run_dir, "plan.json")
        with open(plan_path, "w", encoding="utf-8") as f:
            json.dump(plan, f, ensure_ascii=False, indent=2)
        print(f"  → Plan: {plan.get('object_name', description)}")
        print(f"     Components: {[c['name'] for c in plan.get('components', [])]}")

    # ── Step 2: API-doc retrieval (plan-driven, targeted) ─────────
    api_docs = []
    print("[2/N] API-doc retrieval...")
    api_rag = ApiDocRetriever()
    api_docs = api_rag.retrieve_for_plan(plan)
    print(f"  → Retrieved {len(api_docs)} API ref(s)")

    # ── Step 3: Initial code generation ──────────────────────────
    print("[3/N] Generating bpy code...")
    coder = CoderAgent()
    code = coder.generate(coder_desc, api_docs=api_docs or None)
    _save_code(code, os.path.join(run_dir, "code_v0.py"))
    print(f"  → Code generated ({len(code)} chars)")

    # ── Step 4: Execute + render (+ optional geom check) ─────────
    print("[4/N] Running Blender (iteration 1)...")
    iter_dir = os.path.join(run_dir, "iter_0")
    success, result, geom_issues = _blender_run(code, iter_dir, use_geom=use_geom)

    if not success:
        # Auto-retry: feed the runtime traceback to the fixer
        error_msg = result.get("error_detail") or result["stderr"][-600:]
        print("  ✗ Blender error — auto-fixing...")
        fixer = FixerAgent()
        code = fixer.fix(code, [{
            "source": "blender_runtime",
            "severity": "error",
            "message": f"Script crashed: {error_msg}",
        }])
        _save_code(code, os.path.join(run_dir, "code_v0_fixed.py"))
        iter_dir = os.path.join(run_dir, "iter_0_fixed")
        success, result, geom_issues = _blender_run(code, iter_dir, use_geom=use_geom)
        if not success:
            print("  ✗ Still failing after auto-fix.")
            error_msg2 = result.get("error_detail") or result["stderr"][-400:]
            print(f"  Error: {error_msg2}")
            return {
                "success": False,
                "run_dir": run_dir,
                "iterations": 1,
                "issues_per_round": [],
                "plan": plan,
            }
    print(f"  ✓ Render complete → {iter_dir}/")

    issues_per_round: list[list] = []

    # ── Step 5: Critic + Fixer loop ───────────────────────────────
    if max_iterations > 1 and (use_vlm or use_geom):
        fixer = FixerAgent(verbose=verbose)
        # Geom issues from the *previous* run seed the next fix round
        prev_geom_issues = geom_issues

        for iteration in range(1, max_iterations):
            print(f"\n[5/N] Critic round {iteration}...")
            # Only error-level geom issues reach the Fixer. The geom verifier
            # now emits exclusively zero-false-positive hard errors, but we
            # filter explicitly so any future soft checks can't leak noise in.
            geom_for_fixer = [i for i in prev_geom_issues if i.get("severity") == "error"]
            all_issues: list[dict] = list(geom_for_fixer)

            if use_vlm:
                from agents.vlm_critic import VLMCritic
                critic = VLMCritic(verbose=verbose)
                renders = [
                    os.path.join(iter_dir, f"{v}.png")
                    for v in ("front", "side", "top", "iso")
                ]
                vlm_issues = critic.critique(renders, description)
                all_issues.extend(vlm_issues)
                print(f"  VLM found {len(vlm_issues)} issue(s)")

            if geom_for_fixer:
                print(f"  Geom carried {len(geom_for_fixer)} error(s) from previous run")

            issues_per_round.append(all_issues)

            if not all_issues:
                print("  ✓ No issues found — stopping early")
                break

            print(f"  Fixing {len(all_issues)} issue(s) total...")
            code = fixer.fix(code, all_issues)
            version = f"v{iteration}"
            _save_code(code, os.path.join(run_dir, f"code_{version}.py"))

            iter_dir = os.path.join(run_dir, f"iter_{iteration}")
            success, result, prev_geom_issues = _blender_run(
                code, iter_dir, use_geom=use_geom
            )
            if not success:
                print(f"  ✗ Blender error after fix:\n{result['stderr'][-400:]}")
                break
            print(f"  ✓ Render complete → {iter_dir}/")

    print(f"\n{'='*60}")
    print(f"  Done. Output: {run_dir}")
    print(f"{'='*60}\n")

    return {
        "success": success,
        "run_dir": run_dir,
        "iterations": 1 + len(issues_per_round),
        "issues_per_round": issues_per_round,
        "plan": plan,
    }


def main():
    parser = argparse.ArgumentParser(description="PJ3 Blender Agent")
    parser.add_argument("description", help="Object description (Chinese or English)")
    parser.add_argument("--no-planner", action="store_true", help="Skip planner step")
    parser.add_argument("--iterations", type=int, default=1, help="Max critic-fixer iterations")
    parser.add_argument("--vlm", action="store_true", help="Enable VLM critic (Day 3)")
    parser.add_argument("--geom", action="store_true", help="Enable geometry verifier (Day 4)")
    parser.add_argument("--verbose", action="store_true", help="Print full agent dialogue (issues + fixer I/O)")
    args = parser.parse_args()

    result = pipeline(
        description=args.description,
        use_planner=not args.no_planner,
        max_iterations=args.iterations,
        use_vlm=args.vlm,
        use_geom=args.geom,
        verbose=args.verbose,
    )

    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
