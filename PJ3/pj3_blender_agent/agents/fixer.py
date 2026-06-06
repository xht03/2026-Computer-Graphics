import os
import re
from pathlib import Path
from openai import OpenAI

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"

_FIXER_SYSTEM = (_PROMPT_DIR / "fixer_system.txt").read_text(encoding="utf-8")
_FIXER_USER   = (_PROMPT_DIR / "fixer_user.txt").read_text(encoding="utf-8")

KIMI_BASE_URL = "https://api.moonshot.cn/v1"
KIMI_MODEL    = "kimi-k2.5"


class FixerAgent:
    def __init__(self, model: str = KIMI_MODEL, verbose: bool = False):
        self.client = OpenAI(
            api_key=os.getenv("KIMI_API_KEY"),
            base_url=KIMI_BASE_URL,
        )
        self.model = model
        self.verbose = verbose

    def fix(self, code: str, issues: list[dict]) -> str:
        """Return corrected bpy code. Issues come from VLM critic and/or geom verifier."""
        if not issues:
            return code

        # Group by severity so errors appear first
        sorted_issues = sorted(issues, key=lambda x: {"error": 0, "warning": 1, "info": 2}.get(x.get("severity", "info"), 2))
        issues_text = "\n".join(
            f"- [{i.get('source','?')}] {i.get('severity','warning').upper()}: {i['message']}"
            for i in sorted_issues
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

        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": _FIXER_SYSTEM},
                {"role": "user", "content": _FIXER_USER.format(code=code, issues=issues_text)},
            ],
            temperature=1,
            max_tokens=16000,
        )
        raw = response.choices[0].message.content.strip()

        if self.verbose:
            print(f"\n  ┌── Fixer OUTPUT (raw LLM response, first 60 lines) ──────────")
            for line in raw.splitlines()[:60]:
                print(f"  │  {line}")
            if len(raw.splitlines()) > 60:
                print(f"  │  ... ({len(raw.splitlines())} lines total)")
            print(f"  └─────────────────────────────────────────────────────────────\n")

        raw = re.sub(r"^```(?:python)?\n?", "", raw, flags=re.MULTILINE)
        raw = re.sub(r"\n?```$", "", raw, flags=re.MULTILINE)
        from agents.coder import sanitize_bpy_code
        fixed = sanitize_bpy_code(raw.strip())

        if self.verbose:
            import difflib
            diff = list(difflib.unified_diff(
                code.splitlines(), fixed.splitlines(),
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
