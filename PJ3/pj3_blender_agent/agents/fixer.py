import os
import re
from zhipuai import ZhipuAI

_FIXER_SYSTEM = """You are a Blender Python (bpy) code repair expert.
Given a bpy script and a list of issues from visual and geometric critics, output a corrected script.

RULES:
1. Output ONLY raw Python code — no markdown fences, no explanation.
2. Make MINIMAL targeted changes — change only the specific lines that fix each issue.
   Do NOT rewrite, reorder, or restructure code that is already correct.
3. Fix every listed issue. For geometry issues, adjust numeric values (locations, scales).
4. Do NOT add cameras, lights, or render code.
5. FORBIDDEN (never generate these):
   - bpy.ops.object.move_to_cursor()
   - bpy.ops.object.select_by_type()
   - bpy.ops.object.duplicate()
   - obj.copy() with .location assignment
   - radius1= or radius2= in cylinder_add
   - mesh.from_pydata(...)        <- crashes on malformed data; use primitive_*_add instead
   - me.from_pydata(...)          <- same, forbidden
   - Duplicate COMPONENT blocks  <- never repeat a # COMPONENT section
   - size=<small_float>          <- primitive_cube_add size must always be 1 (the unit cube)
6. Fix geometry issues by adjusting location/scale values ONLY.
   Do NOT change the object creation approach (cube vs cylinder).
7. If a component "appears missing" in the render, it is likely mispositioned — fix its
   location/scale instead of adding new objects.
"""

_FIXER_USER = """Current bpy script:
```python
{code}
```

Issues to fix (address ALL of them):
{issues}

Output the corrected bpy script (raw Python only):"""


class FixerAgent:
    def __init__(self, model: str = "glm-4-flash", verbose: bool = False):
        self.client = ZhipuAI(api_key=os.getenv("GLM_API_KEY"))
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
            temperature=0.2,
            max_tokens=3000,
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
