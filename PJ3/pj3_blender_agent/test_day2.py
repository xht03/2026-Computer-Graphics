#!/usr/bin/env python3
"""
Day 2 milestone test: run pipeline on 3 Chinese furniture descriptions.

Usage:
    python test_day2.py
    python test_day2.py --case 长凳   # run single case
"""
import argparse
import sys
from dotenv import load_dotenv

load_dotenv()
load_dotenv(".env.example")

from main import pipeline

TEST_CASES = [
    "长凳",
    "方桌",
    "屏风",
]


def run_all(cases: list[str]) -> None:
    results = []
    for desc in cases:
        print(f"\n{'#'*60}")
        print(f"# Test case: {desc}")
        print(f"{'#'*60}")
        result = pipeline(
            description=desc,
            use_planner=True,
            use_rag=False,
            max_iterations=1,
        )
        results.append((desc, result["success"], result["run_dir"]))

    print(f"\n{'='*60}")
    print("  Day 2 Test Summary")
    print(f"{'='*60}")
    passed = 0
    for desc, ok, run_dir in results:
        status = "✓ PASS" if ok else "✗ FAIL"
        print(f"  {status}  {desc:10s}  →  {run_dir}")
        if ok:
            passed += 1
    print(f"\n  {passed}/{len(results)} cases passed")
    print(f"{'='*60}\n")

    sys.exit(0 if passed == len(results) else 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", help="Run a single test case")
    args = parser.parse_args()

    cases = [args.case] if args.case else TEST_CASES
    run_all(cases)
