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

from agents import ui
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
    connections: list | None = None,
) -> tuple[bool, dict, list]:
    """Append geometry-check + render snippets and run Blender headlessly.

    Returns:
        (success, result_dict, geom_issues)

    Geometry issues are detected by injecting a bpy verification snippet
    BETWEEN the model code and the render snippet.  The snippet prints JSON
    to stdout; we parse it after the process exits.  This catches problems
    (floating parts, interpenetration, scale errors) that VLM critique misses.

    The plan's *connections* are injected as a module-level CONNECTIONS literal
    so the in-Blender verifier (a separate process that cannot see the plan dict)
    can check declared joints and flag unintended interpenetration. We inject it
    deterministically here rather than trusting the Coder to copy it into code.
    """
    abs_output_dir = str(Path(output_dir).resolve())

    # Inject geometry verification BEFORE render so the clean scene is checked
    geom_snippet = ""
    if use_geom:
        from agents.geom_verifier import GeomVerifier
        conns = [
            [str(p[0]), str(p[1])]
            for p in (connections or [])
            if isinstance(p, (list, tuple)) and len(p) == 2
        ]
        conn_prelude = f"\n\nCONNECTIONS = {conns!r}\n"
        geom_snippet = conn_prelude + "\n" + GeomVerifier.build_verify_snippet()

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
        ui.log(f"geom · {len(geom_issues)} issue(s): {summary}")
        if verbose:
            if geom_issues:
                for iss in geom_issues:
                    sev = iss.get('severity', '?').upper()
                    msg = iss.get('message', '')
                    comp = iss.get('component', '')
                    comp_str = f" [{comp}]" if comp else ""
                    ui.log(f"  [geom] {sev}{comp_str}: {msg}")
            else:
                ui.log("  (geom verifier found no issues)")

    # Surface errors/tracebacks only; success lines are summarised by the caller.
    if result.get("error_detail"):
        ui.log(f"[red]{result['error_detail']}[/]")

    return result["success"], result, geom_issues


# ── Issue helpers ─────────────────────────────────────────────────────────────

_VIEWS = ("front", "side", "top", "iso")


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


def _missing_component_names(code: str, plan: dict) -> list[str]:
    """Plan component names that never appear in the generated code — the tell-tale
    of a truncated/omitted generation (the code does not implement the full plan)."""
    return [
        c["name"] for c in plan.get("components", [])
        if c.get("name") and c["name"] not in code
    ]


def _route_issues(issues: list[dict]) -> tuple[str, list[dict]]:
    """Pick the round's path but ALWAYS carry every issue. If any issue is structural
    the round goes through replan+fix_structure (which fixes structure AND the rest in
    one pass); otherwise fix_parameters handles them. Parameter issues are never dropped."""
    if any(i.get("kind") == "structural" for i in issues):
        return "replan", issues
    return "fix", issues


def _interactive_review(
    iteration: int, iter_dir: str, issues: list[dict], planner
) -> tuple[str, list[dict]]:
    """Show the round's result and ask the user what to do next.

    Returns (action, issues_to_apply) where action is one of:
        "accept"  — finish, keep current result
        "fix"     — Fixer patches issues_to_apply
        "replan"  — Planner revises the plan from issues_to_apply
        "quit"    — abort
    """
    ui.review(iteration, iter_dir, issues)
    ui.console.print()
    ui.menu()
    while True:
        choice = ui.ask().lower()
        if choice in ("a", "accept", ""):
            return "accept", []
        if choice in ("q", "quit"):
            return "quit", []
        if choice in ("f", "fix"):
            if not issues:
                ui.note("没有自动检测到的问题；请用 [c] 自己提意见。")
                continue
            return _route_issues(issues)
        if choice in ("c", "comment"):
            feedback = ui.ask("意见 › ")
            if not feedback:
                ui.note("意见为空，已忽略。")
                continue
            with ui.working("Classifying feedback"):
                kind = planner.classify_feedback(feedback)
            user_issue = {
                "source": "user",
                "severity": "error",
                "kind": kind,
                "message": feedback,
            }
            ui.note(f"→ 已归类为 {kind} 问题")
            return ("replan" if kind == "structural" else "fix"), [user_issue]
        ui.note("无效输入，请输入 a / f / c / q。")


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
    spin = not verbose
    ui.banner(description, vlm=use_vlm, geom=use_geom,
              interactive=interactive, max_iterations=max_iterations)

    # ── Step 1: Plan ──────────────────────────────────────────────
    coder_desc = description
    plan = {}
    planner = PlannerAgent()  # always available for replan / feedback classification
    if use_planner:
        with ui.working("Planning", spinner=spin):
            plan = planner.plan(description)
            coder_desc = planner.plan_to_coder_description(plan)
        with open(os.path.join(run_dir, "plan.json"), "w", encoding="utf-8") as f:
            json.dump(plan, f, ensure_ascii=False, indent=2)
        ui.ok("Plan", f"{plan.get('object_name', description)} · "
                       f"{len(plan.get('components', []))} components")

    # ── Step 2: API-doc retrieval (plan-driven, targeted) ─────────
    with ui.working("Retrieving API docs", spinner=spin):
        api_rag = ApiDocRetriever()
        api_docs = api_rag.retrieve_for_plan(plan)
    ui.ok("API refs", str(len(api_docs)))

    # One shared Fixer for the whole run (crash repair, completeness, the loop).
    fixer = FixerAgent(verbose=verbose, api_rag=api_rag)

    # ── Step 3: Initial code generation ──────────────────────────
    with ui.working("Generating bpy code", spinner=spin):
        coder = CoderAgent()
        code = coder.generate(coder_desc, api_docs=api_docs or None)
    _save_code(code, os.path.join(run_dir, "code_v0.py"))
    ui.ok("Code", f"{len(code)} chars")

    # ── Step 4: Execute + render (+ optional geom check) ─────────
    iter_dir = os.path.join(run_dir, "iter_0")
    with ui.working("Rendering iter_0", spinner=spin):
        success, result, geom_issues = _blender_run(
            code, iter_dir, use_geom=use_geom, connections=plan.get("connections"))

    if not success:
        # Auto-retry: feed the runtime traceback to the fixer
        error_msg = result.get("error_detail") or result["stderr"][-600:]
        ui.err("Blender error", "auto-fixing")
        with ui.working("Fixing crash", spinner=spin):
            code = fixer.fix_parameters(code, [{
                "source": "blender_runtime",
                "severity": "error",
                "message": f"Script crashed: {error_msg}",
            }])
        _save_code(code, os.path.join(run_dir, "code_v0_fixed.py"))
        iter_dir = os.path.join(run_dir, "iter_0_fixed")
        with ui.working("Rendering iter_0_fixed", spinner=spin):
            success, result, geom_issues = _blender_run(
                code, iter_dir, use_geom=use_geom, connections=plan.get("connections"))
        if not success:
            error_msg2 = result.get("error_detail") or result["stderr"][-400:]
            ui.err("Still failing after auto-fix", error_msg2.splitlines()[-1] if error_msg2 else "")
            ui.footer(run_dir, success=False)
            return {
                "success": False,
                "run_dir": run_dir,
                "iterations": 1,
                "issues_per_round": [],
                "plan": plan,
            }
    ui.ok("Render", iter_dir)

    # ── Step 4.5: Completeness — recover components dropped by truncation/omission ──
    # replan can't catch this (it diffs plan-vs-plan); we diff plan-vs-CODE and add
    # whatever the generator failed to emit, so the model actually implements the plan.
    if plan.get("components"):
        missing = _missing_component_names(code, plan)
        if missing:
            ui.warn("Incomplete code", f"{len(missing)} planned component(s) missing: {', '.join(missing)}")
            present = {**plan, "components": [
                c for c in plan["components"] if c.get("name") not in missing]}
            add_diff = summarize_plan_change(present, plan)
            if add_diff.strip():
                with ui.working("Adding missing components", spinner=spin):
                    code = fixer.fix_structure(
                        code, add_diff, api_docs=api_rag.retrieve_for_plan(plan) or None,
                        issues=[{"source": "selfcheck", "severity": "error", "kind": "structural",
                                 "message": "These planned components are missing from the code: "
                                            + ", ".join(missing)}],
                        plan_spec=planner.plan_to_coder_description(plan))
                _save_code(code, os.path.join(run_dir, "code_v0_complete.py"))
                iter_dir = os.path.join(run_dir, "iter_0_complete")
                with ui.working("Rendering iter_0_complete", spinner=spin):
                    success, result, geom_issues = _blender_run(
                        code, iter_dir, use_geom=use_geom, connections=plan.get("connections"))
                still = _missing_component_names(code, plan)
                ui.ok("Components", f"added {len(missing) - len(still)}/{len(missing)}") if success \
                    else ui.err("Render failed after adding components")

    issues_per_round: list[list] = []

    # ── Step 5: Review → (Fix | Replan) loop ──────────────────────
    # Interactive mode: pause each round, show result, let the user steer.
    # Auto mode: run up to max_iterations rounds, routing each batch of issues
    # to the Fixer (parameter problems) or the Planner (structural problems).
    prev_geom_issues = geom_issues
    # Only run the critics up front if we might actually act on them.
    if interactive or max_iterations > 1:
        with ui.working("Critiquing", spinner=spin):
            issues = _detect_issues(iter_dir, description, prev_geom_issues, use_vlm, use_geom, verbose)
    else:
        issues = []
    iteration = 0

    while True:
        # ── Decide what to do this round ──
        if interactive:
            action, to_apply = _interactive_review(iteration, iter_dir, issues, planner)
        else:
            if iteration >= max_iterations - 1:
                break
            if not issues:
                ui.ok("No issues — stopping early")
                break
            action, to_apply = _route_issues(issues)

        if action in ("accept", "quit"):
            break

        issues_per_round.append(issues)
        iteration += 1
        version = f"v{iteration}"

        # ── Apply. A structural round revises the plan, then fixes structure AND every
        # other issue in ONE Fixer pass with the new plan as the spec. A parameter round
        # adjusts numbers (also plan-aware). Either way NO issue is dropped. ──
        n_struct = sum(1 for i in to_apply if i.get("kind") == "structural")
        if action == "replan" and not plan.get("components"):
            action = "fix"  # nothing to replan against

        if action == "replan":
            with ui.working(f"Replanning ({n_struct} structural, {len(to_apply)} total)", spinner=spin):
                old_plan = plan
                plan = planner.replan(plan, [i for i in to_apply if i.get("kind") == "structural"])
            with open(os.path.join(run_dir, f"plan_{version}.json"), "w", encoding="utf-8") as f:
                json.dump(plan, f, ensure_ascii=False, indent=2)
            changes = summarize_plan_change(old_plan, plan)
            if changes.strip():
                ui.plan_diff(changes)
            plan_spec = planner.plan_to_coder_description(plan)
            with ui.working("Restructuring + fixing", spinner=spin):
                code = fixer.fix_structure(
                    code, changes, api_docs=api_rag.retrieve_for_plan(plan) or None,
                    issues=to_apply, plan_spec=plan_spec)
        else:  # fix — parameter-only round
            plan_spec = planner.plan_to_coder_description(plan) if plan.get("components") else ""
            with ui.working(f"Fixing ({len(to_apply)} issue(s))", spinner=spin):
                code = fixer.fix_parameters(code, to_apply, plan_spec=plan_spec)

        _save_code(code, os.path.join(run_dir, f"code_{version}.py"))

        # ── Re-run Blender ──
        iter_dir = os.path.join(run_dir, f"iter_{iteration}")
        with ui.working(f"Rendering iter_{iteration}", spinner=spin):
            success, result, prev_geom_issues = _blender_run(
                code, iter_dir, use_geom=use_geom, verbose=verbose,
                connections=plan.get("connections"),
            )
        if not success:
            error_msg = result.get("error_detail") or result["stderr"][-600:]
            ui.err("Blender error after change")
            issues = [{
                "source": "blender_runtime", "severity": "error",
                "kind": "parameter", "message": f"Script crashed: {error_msg}",
            }]
            if interactive:
                # Surface the crash in the next review so the user can decide.
                continue
            # Auto mode: one fixer pass on the traceback, then give up if still broken.
            with ui.working("Fixing crash", spinner=spin):
                code = fixer.fix_parameters(code, issues)
            _save_code(code, os.path.join(run_dir, f"code_{version}_fixed.py"))
            iter_dir = os.path.join(run_dir, f"iter_{iteration}_fixed")
            with ui.working(f"Rendering iter_{iteration}_fixed", spinner=spin):
                success, result, prev_geom_issues = _blender_run(
                    code, iter_dir, use_geom=use_geom, verbose=verbose,
                    connections=plan.get("connections"),
                )
            if not success:
                ui.err("Still failing — stopping")
                break
        ui.ok("Render", iter_dir)
        with ui.working("Critiquing", spinner=spin):
            issues = _detect_issues(iter_dir, description, prev_geom_issues, use_vlm, use_geom, verbose)

    ui.footer(run_dir, success=success)

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
