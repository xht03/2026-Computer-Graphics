#!/usr/bin/env python3
"""
Day 3 milestone test: VLM Critic + Fixer loop.

Usage:
    python test_day3.py            # run all cases with 2 iterations
    python test_day3.py --case 长凳 --iterations 2
"""
import argparse
import sys
from dotenv import load_dotenv

load_dotenv()
load_dotenv(".env.example")

from main import pipeline

TEST_CASES = ["长凳", "方桌", "屏风"]


def run(cases: list[str], iterations: int) -> None:
    results = []
    for desc in cases:
        print(f"\n{'#'*60}\n# {desc}  (max {iterations} iterations)\n{'#'*60}")
        result = pipeline(
            description=desc,
            use_planner=True,
            use_rag=False,
            max_iterations=iterations,
            use_vlm=True,
            use_geom=False,
        )
        results.append((desc, result["success"], result["run_dir"],
                        result.get("issues_per_round", [])))

    print(f"\n{'='*60}\n  Day 3 Test Summary\n{'='*60}")
    for desc, ok, run_dir, issues in results:
        status = "✓ PASS" if ok else "✗ FAIL"
        total_issues = sum(len(r) for r in issues)
        print(f"  {status}  {desc:8s}  issues={total_issues}  → {run_dir}")
    passed = sum(1 for _, ok, _, _ in results if ok)
    print(f"\n  {passed}/{len(results)} passed\n{'='*60}\n")
    sys.exit(0 if passed == len(results) else 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", help="Single case to test")
    parser.add_argument("--iterations", type=int, default=2)
    args = parser.parse_args()
    cases = [args.case] if args.case else TEST_CASES
    run(cases, args.iterations)
