#!/usr/bin/env python3
"""
Main pipeline entry point.

Usage:
    # Interactive (default): pause each round to review renders and steer
    python main.py "做一张明式案几" --vlm --geom
    python main.py "a red cube"

    # Auto / batch (non-interactive): iterate up to --iterations rounds
    python main.py "长凳" --auto --iterations 2 --vlm --geom

Interactive menu each round:
    [a] accept & finish   [f] fix detected issues
    [c] type your own feedback (auto-routed: structural -> replan, else -> fixer)
    [q] quit
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
from agents.planner import PlannerAgent, summarize_plan_change
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


# ── Issue helpers ─────────────────────────────────────────────────────────────

_VIEWS = ("front", "side", "top", "iso")


def _print_renders(iter_dir: str) -> bool:
    """Print clickable paths to the 4-view renders. Returns True if any exist."""
    found = False
    print("  Renders:")
    for v in _VIEWS:
        p = os.path.join(iter_dir, f"{v}.png")
        if os.path.exists(p):
            print(f"    {v:5s} → {p}")
            found = True
    if not found:
        print("    (no render images — the script likely crashed)")
    return found


def _print_issues(issues: list[dict]) -> None:
    if not issues:
        print("  Auto-detected issues: none")
        return
    print(f"  Auto-detected issues ({len(issues)}):")
    order = {"error": 0, "warning": 1, "info": 2}
    for i in sorted(issues, key=lambda x: order.get(x.get("severity", "info"), 2)):
        src = i.get("source", "?")
        sev = i.get("severity", "?").upper()
        kind = i.get("kind", "parameter")
        print(f"    [{src}/{kind}] {sev}: {i.get('message', '')}")


def _detect_issues(
    iter_dir: str,
    description: str,
    geom_issues: list[dict],
    use_vlm: bool,
    use_geom: bool,
    verbose: bool,
) -> list[dict]:
    """Run the enabled critics on a finished render and return their issues."""
    issues: list[dict] = []
    if use_geom:
        issues.extend(i for i in geom_issues if i.get("severity") == "error")
    if use_vlm:
        from agents.vlm_critic import VLMCritic
        renders = [os.path.join(iter_dir, f"{v}.png") for v in _VIEWS]
        issues.extend(VLMCritic(verbose=verbose).critique(renders, description))
    return issues


def _route_issues(issues: list[dict]) -> tuple[str, list[dict]]:
    """Split issues by kind. If any structural issue is present, route the whole
    round to replan; otherwise patch the parameter issues with the Fixer."""
    structural = [i for i in issues if i.get("kind") == "structural"]
    if structural:
        return "replan", structural
    return "fix", issues


def _interactive_review(
    iter_dir: str, issues: list[dict], planner
) -> tuple[str, list[dict]]:
    """Show the round's result and ask the user what to do next.

    Returns (action, issues_to_apply) where action is one of:
        "accept"  — finish, keep current result
        "fix"     — Fixer patches issues_to_apply
        "replan"  — Planner revises the plan from issues_to_apply
        "quit"    — abort
    """
    _print_renders(iter_dir)
    _print_issues(issues)
    print(
        "\n  What next?\n"
        "    [a] 接受并结束\n"
        "    [f] 按自动检测到的问题修改\n"
        "    [c] 我来提意见（自由输入）\n"
        "    [q] 退出"
    )
    while True:
        try:
            choice = input("  > ").strip().lower()
        except EOFError:
            return "accept", []
        if choice in ("a", "accept", ""):
            return "accept", []
        if choice in ("q", "quit"):
            return "quit", []
        if choice in ("f", "fix"):
            if not issues:
                print("  (没有自动检测到的问题；请用 [c] 自己提意见。)")
                continue
            return _route_issues(issues)
        if choice in ("c", "comment"):
            try:
                feedback = input("  请输入你的修改意见：").strip()
            except EOFError:
                feedback = ""
            if not feedback:
                print("  (意见为空，已忽略。)")
                continue
            kind = planner.classify_feedback(feedback)
            user_issue = {
                "source": "user",
                "severity": "error",
                "kind": kind,
                "message": feedback,
            }
            print(f"  → 已归类为 {kind} 问题。")
            return ("replan" if kind == "structural" else "fix"), [user_issue]
        print("  (无效输入，请输入 a / f / c / q。)")


def pipeline(
    description: str,
    use_planner: bool = True,
    max_iterations: int = 1,
    use_vlm: bool = False,
    use_geom: bool = False,
    interactive: bool = True,
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
    print(f"  Mode        : {'interactive' if interactive else f'auto (max {max_iterations} rounds)'}")
    print(f"{'='*60}\n")

    # ── Step 1: Plan ──────────────────────────────────────────────
    coder_desc = description
    plan = {}
    planner = PlannerAgent()  # always available for replan / feedback classification
    if use_planner:
        print("[1/N] Planning...")
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
        fixer = FixerAgent(api_rag=api_rag)
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

    # ── Step 5: Review → (Fix | Replan) loop ──────────────────────
    # Interactive mode: pause each round, show result, let the user steer.
    # Auto mode: run up to max_iterations rounds, routing each batch of issues
    # to the Fixer (parameter problems) or the Planner (structural problems).
    fixer = FixerAgent(verbose=verbose, api_rag=api_rag)
    prev_geom_issues = geom_issues
    # Only run the critics up front if we might actually act on them.
    if interactive or max_iterations > 1:
        issues = _detect_issues(iter_dir, description, prev_geom_issues, use_vlm, use_geom, verbose)
    else:
        issues = []
    iteration = 0

    while True:
        # ── Decide what to do this round ──
        if interactive:
            print(f"\n[review] 第 {iteration} 轮结果 → {iter_dir}/")
            action, to_apply = _interactive_review(iter_dir, issues, planner)
        else:
            if iteration >= max_iterations - 1:
                break
            if not issues:
                print("  ✓ No issues found — stopping early")
                break
            action, to_apply = _route_issues(issues)

        if action in ("accept", "quit"):
            if action == "quit":
                print("  已退出。")
            break

        issues_per_round.append(issues)
        iteration += 1
        version = f"v{iteration}"

        # ── Apply: structural → replan + regenerate; parameter → fixer patch ──
        if action == "replan" and not plan:
            print("  (无可修订的计划，改用 Fixer 处理。)")
            action = "fix"

        if action == "replan":
            print(f"  ↻ Replanning from {len(to_apply)} structural issue(s)...")
            old_plan = plan
            plan = planner.replan(plan, to_apply)
            with open(os.path.join(run_dir, f"plan_{version}.json"), "w", encoding="utf-8") as f:
                json.dump(plan, f, ensure_ascii=False, indent=2)
            print(f"     Components: {[c['name'] for c in plan.get('components', [])]}")
            # Plan-only revise → Fixer applies the diff INCREMENTALLY, so already-tuned
            # parts are preserved instead of being regenerated from scratch.
            changes = summarize_plan_change(old_plan, plan)
            if changes.strip():
                print("  Plan diff → applying incrementally:")
                for ln in changes.splitlines():
                    print(f"    {ln}")
                new_api_docs = api_rag.retrieve_for_plan(plan)
                code = fixer.apply_plan_change(code, changes, api_docs=new_api_docs or None)
            else:
                print("  (计划无结构性变化 — 退回普通 Fixer 微调。)")
                code = fixer.fix(code, to_apply)
        else:  # fix
            print(f"  ✎ Fixing {len(to_apply)} issue(s)...")
            code = fixer.fix(code, to_apply)

        _save_code(code, os.path.join(run_dir, f"code_{version}.py"))

        # ── Re-run Blender ──
        iter_dir = os.path.join(run_dir, f"iter_{iteration}")
        success, result, prev_geom_issues = _blender_run(
            code, iter_dir, use_geom=use_geom, verbose=verbose
        )
        if not success:
            error_msg = result.get("error_detail") or result["stderr"][-600:]
            print("  ✗ Blender error after change.")
            issues = [{
                "source": "blender_runtime", "severity": "error",
                "kind": "parameter", "message": f"Script crashed: {error_msg}",
            }]
            if interactive:
                # Surface the crash in the next review so the user can decide.
                continue
            # Auto mode: one fixer pass on the traceback, then give up if still broken.
            code = fixer.fix(code, issues)
            _save_code(code, os.path.join(run_dir, f"code_{version}_fixed.py"))
            iter_dir = os.path.join(run_dir, f"iter_{iteration}_fixed")
            success, result, prev_geom_issues = _blender_run(
                code, iter_dir, use_geom=use_geom, verbose=verbose
            )
            if not success:
                print("  ✗ Still failing — stopping.")
                break
        print(f"  ✓ Render complete → {iter_dir}/")
        issues = _detect_issues(iter_dir, description, prev_geom_issues, use_vlm, use_geom, verbose)

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
    parser.add_argument("--iterations", type=int, default=1, help="Max rounds in --auto mode")
    parser.add_argument("--vlm", action="store_true", help="Enable VLM critic (Day 3)")
    parser.add_argument("--geom", action="store_true", help="Enable geometry verifier (Day 4)")
    parser.add_argument("--auto", action="store_true",
                        help="Non-interactive batch mode: auto-iterate up to --iterations rounds")
    parser.add_argument("--verbose", action="store_true", help="Print full agent dialogue (issues + fixer I/O)")
    args = parser.parse_args()

    result = pipeline(
        description=args.description,
        use_planner=not args.no_planner,
        max_iterations=args.iterations,
        use_vlm=args.vlm,
        use_geom=args.geom,
        interactive=not args.auto,
        verbose=args.verbose,
    )

    sys.exit(0 if result["success"] else 1)


if __name__ == "__main__":
    main()
