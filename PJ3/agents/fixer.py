import os
import re
from pathlib import Path
from openai import OpenAI

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"

_FIXER_SYSTEM = (_PROMPT_DIR / "fixer_system.txt").read_text(encoding="utf-8")
_FIXER_USER   = (_PROMPT_DIR / "fixer_user.txt").read_text(encoding="utf-8")
_REPLAN_SYSTEM = (_PROMPT_DIR / "fixer_replan_system.txt").read_text(encoding="utf-8")
_REPLAN_USER   = (_PROMPT_DIR / "fixer_replan_user.txt").read_text(encoding="utf-8")

KIMI_BASE_URL = "https://api.moonshot.cn/v1"
KIMI_MODEL    = "kimi-k2.5"


class FixerAgent:
    def __init__(self, model: str = KIMI_MODEL, verbose: bool = False,
                 use_api_rag: bool = True, api_rag=None):
        self.client = OpenAI(
            api_key=os.getenv("KIMI_API_KEY"),
            base_url=KIMI_BASE_URL,
        )
        self.model = model
        self.verbose = verbose
        # Prefer a shared retriever so the embedding model is loaded only once per
        # run; fall back to building our own only if none was passed in.
        self._api_rag = api_rag
        if self._api_rag is None and use_api_rag:
            try:
                from agents.api_rag import ApiDocRetriever
                self._api_rag = ApiDocRetriever()
            except Exception:
                self._api_rag = None

    def _retrieve_api_docs(self, issues: list[dict]) -> str:
        """For runtime crashes, fetch exact API signatures/sockets for the error."""
        if self._api_rag is None:
            return ""
        # Only error-severity messages (tracebacks) carry useful API signal.
        err_text = "\n".join(
            i.get("message", "") for i in issues if i.get("severity") == "error"
        )
        if not err_text.strip():
            return ""
        try:
            chunks = self._api_rag.retrieve_for_error(err_text, top_k=4)
        except Exception:
            return ""
        block = self._api_rag.as_prompt_block(
            chunks,
            "\nAPI REFERENCE (exact signatures from THIS Blender — trust over memory):",
        )
        if block and self.verbose:
            print(f"\n  ┌── Fixer API-RAG ({len(chunks)} chunk(s)) ─────────────────")
            for c in chunks:
                print(f"  │  {c.splitlines()[0]}")
            print(f"  └─────────────────────────────────────────────────────────────\n")
        return block

    def fix_parameters(self, code: str, issues: list[dict], plan_spec: str = "") -> str:
        """PARAMETER-level repair: adjust numbers (positions, sizes, colours) to fix
        the given issues. Forbidden from adding/removing parts or changing a part's
        construction — that is fix_structure()'s job. Issues come from the runtime,
        geom verifier, VLM critic, or the user. *plan_spec* (the plan's target sizes/
        positions) is given so numeric fixes aim at the intended values."""
        if not issues:
            return code

        # Group by severity so errors appear first
        sorted_issues = sorted(issues, key=lambda x: {"error": 0, "warning": 1, "info": 2}.get(x.get("severity", "info"), 2))
        issues_text = "\n".join(
            f"- [{i.get('source','?')}] {i.get('severity','warning').upper()}: {i['message']}"
            for i in sorted_issues
        )
        plan_block = (
            f"\nPLAN SPEC (the authoritative target — fix numbers toward these):\n{plan_spec}\n"
            if plan_spec else ""
        )

        if self.verbose:
            print(f"\n  ┌── Fixer INPUT ({len(sorted_issues)} issues) ─────────────────────────")
            for iss in sorted_issues:
                src = iss.get('source', '?')
                sev = iss.get('severity', '?').upper()
                msg = iss.get('message', '')
                comp = iss.get('component', '')
                comp_str = f" [{comp}]" if comp else ""
                print(f"  │  [{src}] {sev}{comp_str}: {msg}")
            print(f"  └─────────────────────────────────────────────────────────────\n")

        api_docs = self._retrieve_api_docs(sorted_issues)

        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": _FIXER_SYSTEM},
                {"role": "user", "content": _FIXER_USER.format(
                    code=code, issues=issues_text, api_docs=api_docs, plan_spec=plan_block)},
            ],
            temperature=1,
            max_tokens=64000,
        )
        raw = response.choices[0].message.content.strip()

        if self.verbose:
            print(f"\n  ┌── Fixer OUTPUT (raw LLM response, first 60 lines) ──────────")
            for line in raw.splitlines()[:60]:
                print(f"  │  {line}")
            if len(raw.splitlines()) > 60:
                print(f"  │  ... ({len(raw.splitlines())} lines total)")
            print(f"  └─────────────────────────────────────────────────────────────\n")

        return self._postprocess(raw, code)

    def fix_structure(
        self, code: str, changes: str, api_docs: list[str] | None = None,
        issues: list[dict] | None = None, plan_spec: str = "",
    ) -> str:
        """STRUCTURE-level repair for a round that has structural problems. Applies the
        ADD/REMOVE/RESHAPE diff in *changes*, AND fixes every other problem in *issues*
        (parameter ones too — proportions, broken joints) in the SAME pass, using
        *plan_spec* (the revised plan's target sizes/positions) as the authoritative
        spec. This is one round addressing all feedback at once, not just structure.
        *changes* is the diff from planner.summarize_plan_change.
        """
        if not changes.strip() and not (issues and plan_spec):
            return code

        if api_docs:
            api_block = (
                "\nAPI REFERENCE (exact signatures from THIS Blender — trust over memory):\n"
                + "\n\n".join(api_docs) + "\n"
            )
        else:
            api_block = ""

        plan_block = (
            f"\nPLAN SPEC (the authoritative target — every part should match this):\n{plan_spec}\n"
            if plan_spec else ""
        )

        if issues:
            context = (
                "ALL problems to fix this round (apply the structural changes below AND "
                "correct these — proportions, broken joints, colours, placement):\n"
                + "\n".join(
                    f"- [{i.get('source','?')}/{i.get('kind','parameter')}] {i.get('message', '')}"
                    for i in issues if i.get("message"))
                + "\n\n"
            )
        else:
            context = ""

        if self.verbose:
            print(f"\n  ┌── Fixer(structure) changes ─────────────────────────────────")
            for line in changes.splitlines():
                print(f"  │  {line}")
            print(f"  └─────────────────────────────────────────────────────────────\n")

        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": _REPLAN_SYSTEM},
                {"role": "user", "content": _REPLAN_USER.format(
                    code=code, changes=changes, api_docs=api_block,
                    context=context, plan_spec=plan_block)},
            ],
            temperature=1,
            max_tokens=64000,
        )
        raw = response.choices[0].message.content.strip()
        return self._postprocess(raw, code)

    def _postprocess(self, raw: str, original_code: str) -> str:
        """Strip markdown fences, sanitize, and (verbose) print a diff summary."""
        raw = re.sub(r"^```(?:python)?\n?", "", raw, flags=re.MULTILINE)
        raw = re.sub(r"\n?```$", "", raw, flags=re.MULTILINE)
        from agents.coder import sanitize_bpy_code
        fixed = sanitize_bpy_code(raw.strip())

        if self.verbose:
            import difflib
            diff = list(difflib.unified_diff(
                original_code.splitlines(), fixed.splitlines(),
                fromfile="before", tofile="after", lineterm=""
            ))
            added   = [l[1:] for l in diff if l.startswith("+") and not l.startswith("+++")]
            removed = [l[1:] for l in diff if l.startswith("-") and not l.startswith("---")]
            print(f"  ── Diff summary: +{len(added)} lines, -{len(removed)} lines")
            for l in removed[:6]:
                print(f"  \033[91m-  {l}\033[0m")
            for l in added[:6]:
                print(f"  \033[92m+  {l}\033[0m")
            if not diff:
                print(f"  !! WARNING: Fixer output is IDENTICAL to input — no changes made!")
            print()

        return fixed
