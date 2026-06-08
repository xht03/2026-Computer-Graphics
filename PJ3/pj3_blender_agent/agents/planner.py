import json
import os
import re
from pathlib import Path
from openai import OpenAI

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"

KIMI_BASE_URL = "https://api.moonshot.cn/v1"
KIMI_MODEL = "kimi-k2.5"


def summarize_plan_change(old_plan: dict, new_plan: dict) -> str:
    """Diff two plans by component name and render the change as Fixer instructions.

    Returns an ADD / REMOVE / CHANGE block the structural Fixer can apply
    incrementally, or "" if the component set is structurally identical (only
    style notes / wording changed). Connections of added parts are surfaced so the
    Fixer keeps them attached.
    """
    def by_name(p: dict) -> dict:
        return {c.get("name"): c for c in p.get("components", []) if c.get("name")}

    old_c, new_c = by_name(old_plan), by_name(new_plan)
    added   = [new_c[n] for n in new_c if n not in old_c]
    removed = [old_c[n] for n in old_c if n not in new_c]
    changed = []
    for n in new_c:
        if n not in old_c:
            continue
        o, nw = old_c[n], new_c[n]
        diffs = [
            f"{f}: {o.get(f, '?')} -> {nw.get(f, '?')}"
            for f in ("shape", "approx_size_m", "count", "position_note")
            if str(o.get(f, "")) != str(nw.get(f, ""))
        ]
        if diffs:
            changed.append((nw, diffs))

    conns = new_plan.get("connections") or []

    def partners_of(name: str) -> list[str]:
        out = []
        for pair in conns:
            if isinstance(pair, (list, tuple)) and len(pair) == 2:
                if pair[0] == name:
                    out.append(pair[1])
                elif pair[1] == name:
                    out.append(pair[0])
        return out

    lines: list[str] = []
    if added:
        lines.append("ADD these NEW components (write a new '# COMPONENT' block for each):")
        for c in added:
            ln = f"  - {c.get('name')} (x{c.get('count', 1)}, {c.get('shape', 'cube')}): {c.get('description', '')}"
            if c.get("approx_size_m"):
                ln += f", size={c['approx_size_m']}"
            if c.get("position_note"):
                ln += f", position: {c['position_note']}"
            partners = partners_of(c.get("name"))
            if partners:
                ln += f", must physically touch: {', '.join(partners)}"
            lines.append(ln)
    if removed:
        lines.append("REMOVE these components (delete their entire '# COMPONENT' block):")
        for c in removed:
            lines.append(f"  - {c.get('name')}")
    if changed:
        lines.append("CHANGE these existing components (edit their block in place):")
        for c, diffs in changed:
            lines.append(f"  - {c.get('name')}: " + "; ".join(diffs))
    return "\n".join(lines)


def _extract_json_object(raw: str) -> str:
    """Strip markdown fences and isolate the outermost JSON object."""
    raw = raw.strip()
    raw = re.sub(r"^```(?:json)?\n?", "", raw, flags=re.MULTILINE)
    raw = re.sub(r"\n?```$", "", raw, flags=re.MULTILINE)
    # Reasoning models may wrap the JSON in prose; isolate the outermost object.
    start, end = raw.find("{"), raw.rfind("}")
    if start != -1 and end > start:
        raw = raw[start : end + 1]
    return raw


class PlannerAgent:
    """Decomposes a user description into a structured JSON task plan."""

    def __init__(self, model: str = KIMI_MODEL):
        self.client = OpenAI(
            api_key=os.getenv("KIMI_API_KEY"),
            base_url=KIMI_BASE_URL,
        )
        self.model = model
        self._system = (_PROMPT_DIR / "planner_system.txt").read_text(encoding="utf-8")
        self._replan_user = (_PROMPT_DIR / "planner_replan_user.txt").read_text(encoding="utf-8")

    def plan(self, description: str) -> dict:
        """Return a parsed plan dict from a Chinese or English description."""
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": self._system},
                {"role": "user", "content": f"Plan the 3D modeling task for: {description}"},
            ],
            temperature=1,
            max_tokens=16000,
        )
        raw = _extract_json_object(response.choices[0].message.content)
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            # Fallback: return minimal plan so pipeline can continue
            return {
                "object_name": "object",
                "object_name_zh": description,
                "description_en": description,
                "components": [],
                "style_notes": "",
                "estimated_difficulty": "simple",
            }

    def replan(self, plan: dict, issues: list[dict]) -> dict:
        """Revise an existing plan to fix STRUCTURAL problems (missing/extra parts,
        wrong shape). Returns the revised plan, or the original if parsing fails."""
        if not issues:
            return plan
        issues_text = "\n".join(
            f"- {i.get('message', '')}" for i in issues if i.get("message")
        )
        user_msg = self._replan_user.format(
            plan=json.dumps(plan, ensure_ascii=False, indent=2),
            issues=issues_text,
        )
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": self._system},
                {"role": "user", "content": user_msg},
            ],
            temperature=1,
            max_tokens=16000,
        )
        raw = _extract_json_object(response.choices[0].message.content)
        try:
            return json.loads(raw)
        except json.JSONDecodeError:
            # Keep the old plan rather than crashing the pipeline.
            return plan

    def classify_feedback(self, feedback: str) -> str:
        """Classify free-form user feedback as 'structural' or 'parameter'.

        structural -> needs the plan to change (add/remove a part, change a shape).
        parameter  -> only adjusts existing parts (position/size/colour/proportion).
        Defaults to 'parameter' on any ambiguity or error.
        """
        system = (
            "Classify a user's feedback about a 3D model into exactly ONE word: "
            "'structural' or 'parameter'. "
            "structural = the fix requires ADDING or REMOVING a part, or changing a "
            "part's basic shape/primitive, or reworking the overall form. "
            "parameter = the fix only adjusts EXISTING parts: position, size, "
            "proportion, rotation, colour, or material. "
            "When unsure, answer 'parameter'. Output ONLY the one word, nothing else."
        )
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": system},
                    {"role": "user", "content": feedback},
                ],
                temperature=0,
                max_tokens=8,
            )
            word = (response.choices[0].message.content or "").strip().lower()
            return "structural" if "structural" in word else "parameter"
        except Exception:
            return "parameter"

    def plan_to_coder_description(self, plan: dict) -> str:
        """Convert a plan dict to a rich description string for CoderAgent."""
        lines = []
        name = plan.get("object_name", "furniture")
        name_zh = plan.get("object_name_zh", "")
        lines.append(f"{name} ({name_zh}): {plan.get('description_en', name)}")

        if plan.get("overall_dimensions_m"):
            lines.append(f"Overall size: {plan['overall_dimensions_m']}")

        if plan.get("style_notes"):
            lines.append(f"Style: {plan['style_notes']}")

        if plan.get("components"):
            lines.append("Components to model:")
            for c in plan["components"]:
                count = c.get("count", 1)
                size = c.get("approx_size_m", "")
                pos = c.get("position_note", "")
                shape = c.get("shape", "cube")
                desc = c.get("description", "")
                part_of = c.get("part_of", "")
                line = f"  - {c['name']} (x{count}, {shape}): {desc}"
                if size:
                    line += f", size={size}"
                if pos:
                    line += f", position: {pos}"
                if part_of:
                    line += f", part of: {part_of}"
                lines.append(line)

        connections = plan.get("connections")
        if connections:
            lines.append("Required connections (these parts MUST physically touch):")
            for pair in connections:
                if isinstance(pair, (list, tuple)) and len(pair) == 2:
                    lines.append(f"  - {pair[0]} <-> {pair[1]}")

        return "\n".join(lines)
