#!/usr/bin/env python3
"""
Day 5 milestone test: Domain RAG + Chinese furniture cases.

Tests:
  1. RAG retrieval accuracy: verify top-1 result is semantically relevant
  2. RAG-augmented pipeline: 5-6 Chinese furniture cases with --rag flag
  3. RAG ablation: with-RAG vs without-RAG (optional, slow)

Usage:
    cd pj3_blender_agent
    python test_day5.py                        # quick retrieval + 3 pipeline cases
    python test_day5.py --retrieval-only       # only test RAG retrieval (no Blender)
    python test_day5.py --cases 长凳 方桌 案几  # custom cases
    python test_day5.py --full                 # all 6 cases + ablation
"""
import argparse
import json
import os
import sys
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()
load_dotenv(".env.example")


# ── RAG retrieval unit test ─────────────────────────────────────────────────

def test_retrieval() -> bool:
    from agents.rag import RAGRetriever

    print(f"\n{'='*60}\n  RAG Retrieval Test\n{'='*60}")
    rag = RAGRetriever()

    # Expected: query → at least one snippet whose tags/title match
    test_pairs = [
        ("椅子的腿",          ["椅腿", "方腿", "锥形腿", "腿"],    ),
        ("弧形靠背圈椅",       ["椅圈", "弧形", "圈椅"],           ),
        ("桌子回纹装饰材质",   ["回纹", "雷纹", "装饰", "纹样"],   ),
        ("八仙桌方桌整体",     ["方桌", "八仙桌", "桌子"],         ),
        ("屏风三折框架",       ["屏风", "折屏", "框架"],           ),
        ("案几条案牙板",       ["案几", "条案", "牙板"],           ),
    ]

    passed = 0
    for query, expected_keywords in test_pairs:
        results = rag.retrieve_with_meta(query, top_k=3)
        if not results:
            print(f"  ✗ '{query}' → no results")
            continue

        # Check if any expected keyword appears in any top result
        all_text = " ".join(
            r["title"] + " " + r["description"] + " " + " ".join(r["tags"])
            for r in results
        )
        found = [kw for kw in expected_keywords if kw in all_text]
        ok = len(found) > 0

        top_titles = [r["title"] for r in results[:2]]
        score_str  = f"score={results[0].get('score', 0):.3f}" if results else ""
        status     = "✓" if ok else "✗"
        print(f"  {status} '{query}'")
        print(f"      top: {top_titles}  {score_str}")
        if not ok:
            print(f"      expected any of {expected_keywords}")
        if ok:
            passed += 1

    total = len(test_pairs)
    print(f"\n  {passed}/{total} retrieval tests passed")
    print(f"{'='*60}\n")
    return passed >= total * 0.7   # 70% pass threshold (fallback mode ok)


# ── Pipeline tests with RAG ─────────────────────────────────────────────────

ALL_CASES = ["长凳", "案几", "方桌", "屏风", "太师椅", "八仙桌"]
DEFAULT_CASES = ["长凳", "案几", "方桌"]


def run_pipeline_cases(cases: list[str], iterations: int, use_vlm: bool, use_geom: bool) -> bool:
    from main import pipeline

    print(f"\n{'='*60}\n  Day 5 Pipeline Test (RAG enabled)\n{'='*60}")
    results = []
    for desc in cases:
        print(f"\n{'#'*60}\n# {desc}  (rag=True, vlm={use_vlm}, geom={use_geom})\n{'#'*60}")
        result = pipeline(
            description=desc,
            use_planner=True,
            use_rag=True,
            max_iterations=iterations,
            use_vlm=use_vlm,
            use_geom=use_geom,
        )
        results.append((desc, result["success"], result["run_dir"]))

    print(f"\n{'='*60}\n  Day 5 Pipeline Summary\n{'='*60}")
    for desc, ok, run_dir in results:
        status = "✓ PASS" if ok else "✗ FAIL"
        print(f"  {status}  {desc:8s}  → {run_dir}")
    passed = sum(1 for _, ok, _ in results if ok)
    print(f"\n  {passed}/{len(results)} cases passed")
    print(f"{'='*60}\n")
    return passed == len(results)


# ── RAG ablation: with vs without ─────────────────────────────────────────

def run_ablation(cases: list[str] | None = None) -> None:
    from main import pipeline

    test_cases = cases or ["长凳", "案几", "太师椅"]
    print(f"\n{'='*60}\n  Day 5 Ablation: With-RAG vs Without-RAG\n{'='*60}")
    rows = []
    for desc in test_cases:
        for use_rag, label in [(False, "no_rag"), (True, "rag")]:
            print(f"\n── {desc} ({label}) ──")
            result = pipeline(
                description=desc,
                use_planner=True,
                use_rag=use_rag,
                max_iterations=1,
                use_vlm=False,
                use_geom=False,
            )
            rows.append({
                "case": desc,
                "rag": use_rag,
                "success": result["success"],
                "run_dir": result["run_dir"],
            })

    print(f"\n{'='*60}\n  Ablation Results\n{'='*60}")
    print(f"  {'Case':<12} {'RAG':>5} {'Success':>8}")
    print(f"  {'-'*28}")
    for r in rows:
        status = "✓" if r["success"] else "✗"
        print(f"  {r['case']:<12} {str(r['rag']):>5} {status:>8}")

    # Save for paper
    out_path = os.path.join("outputs", "ablation_day5_rag.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, ensure_ascii=False, indent=2)
    print(f"\n  Results saved → {out_path}")
    print(f"{'='*60}\n")


# ── Entry point ────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Day 5 RAG + Chinese Furniture Tests")
    parser.add_argument("--retrieval-only", action="store_true",
                        help="Only test RAG retrieval (no Blender run)")
    parser.add_argument("--full", action="store_true",
                        help="Run all 6 cases + ablation")
    parser.add_argument("--ablation", action="store_true",
                        help="Run RAG ablation (with vs without)")
    parser.add_argument("--cases", nargs="+", help="Custom case list")
    parser.add_argument("--iterations", type=int, default=1)
    parser.add_argument("--vlm", action="store_true")
    parser.add_argument("--geom", action="store_true")
    args = parser.parse_args()

    all_ok = True

    if args.retrieval_only:
        ok = test_retrieval()
        sys.exit(0 if ok else 1)

    # Always run retrieval test first
    ok = test_retrieval()
    all_ok = all_ok and ok

    if args.ablation:
        run_ablation(args.cases)
    elif args.full:
        ok = run_pipeline_cases(ALL_CASES, args.iterations, args.vlm, args.geom)
        all_ok = all_ok and ok
        run_ablation()
    else:
        cases = args.cases or DEFAULT_CASES
        ok = run_pipeline_cases(cases, args.iterations, args.vlm, args.geom)
        all_ok = all_ok and ok

    sys.exit(0 if all_ok else 1)
