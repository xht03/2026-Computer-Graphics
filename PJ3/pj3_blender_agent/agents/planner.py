import json
import os
import re
from pathlib import Path
from zhipuai import ZhipuAI

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"


class PlannerAgent:
    """Decomposes a user description into a structured JSON task plan."""

    def __init__(self, model: str = "glm-4-flash"):
        self.client = ZhipuAI(api_key=os.getenv("GLM_API_KEY"))
        self.model = model
        self._system = (_PROMPT_DIR / "planner_system.txt").read_text(encoding="utf-8")

    def plan(self, description: str) -> dict:
        """Return a parsed plan dict from a Chinese or English description."""
        response = self.client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": self._system},
                {"role": "user", "content": f"Plan the 3D modeling task for: {description}"},
            ],
            temperature=0.1,
            max_tokens=1024,
        )
        raw = response.choices[0].message.content.strip()
        # Strip markdown fences if present
        raw = re.sub(r"^```(?:json)?\n?", "", raw, flags=re.MULTILINE)
        raw = re.sub(r"\n?```$", "", raw, flags=re.MULTILINE)
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
                line = f"  - {c['name']} (x{count}, {shape}): {desc}"
                if size:
                    line += f", size={size}"
                if pos:
                    line += f", position: {pos}"
                lines.append(line)

        return "\n".join(lines)
