"""VLM Critic Agent.

Primary:  Kimi k2.5  (Moonshot AI, vision-capable, OpenAI-compatible API)
Fallback: GLM-4V-Flash (ZhipuAI)
"""
import base64
import json
import os
import re
from pathlib import Path

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"

_SYSTEM = (_PROMPT_DIR / "vlm_critic_system.txt").read_text(encoding="utf-8")
_USER_TMPL = (_PROMPT_DIR / "vlm_critic_user.txt").read_text(encoding="utf-8")

_VIEW_LABELS = {
    "front": "Front view",
    "side":  "Side view",
    "top":   "Top view",
    "iso":   "Isometric (45°) view",
}

KIMI_BASE_URL = "https://api.moonshot.cn/v1"
KIMI_MODEL    = "kimi-k2.5"


def _b64(path: str) -> str:
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


class VLMCritic:
    """Critiques 4-view renders. Primary: Kimi k2.5. Fallback: GLM-4V-Flash."""

    def __init__(self, verbose: bool = False):
        self._kimi_client  = None
        self._glm_client   = None
        self.verbose       = verbose

    # ── clients (lazy) ────────────────────────────────────────────────────────

    def _get_kimi(self):
        if self._kimi_client is None:
            from openai import OpenAI
            self._kimi_client = OpenAI(
                api_key=os.getenv("KIMI_API_KEY"),
                base_url=KIMI_BASE_URL,
            )
        return self._kimi_client

    def _get_glm(self):
        if self._glm_client is None:
            from zhipuai import ZhipuAI
            self._glm_client = ZhipuAI(api_key=os.getenv("GLM_API_KEY"))
        return self._glm_client

    # ── main entry ─────────────────────────────────────────────────────────────

    def critique(self, image_paths: list[str], description: str) -> list[dict]:
        """Return list of issue dicts: {source, severity, message}."""
        existing = [(p, Path(p).stem) for p in image_paths if os.path.exists(p)]
        if not existing:
            print("  [VLM] No render images found — skipping.")
            return []

        try:
            return self._critique_kimi(existing, description)
        except Exception as e:
            print(f"  [VLM] Kimi failed ({e}), falling back to GLM-4V.")
            return self._critique_glm4v(existing, description)

    # ── Kimi path ──────────────────────────────────────────────────────────────

    def _critique_kimi(self, existing: list[tuple], description: str) -> list[dict]:
        content = []
        for path, stem in existing:
            label = _VIEW_LABELS.get(stem, stem)
            content.append({"type": "text", "text": f"[{label}]"})
            content.append({
                "type": "image_url",
                "image_url": {"url": f"data:image/png;base64,{_b64(path)}"},
            })
        content.append({"type": "text", "text": _USER_TMPL.format(description=description)})

        response = self._get_kimi().chat.completions.create(
            model=KIMI_MODEL,
            messages=[
                {"role": "system", "content": _SYSTEM},
                {"role": "user",   "content": content},
            ],
            temperature=1,
            max_tokens=16000,
        )
        choice = response.choices[0]
        raw    = choice.message.content or ""
        finish = choice.finish_reason
        usage  = getattr(response, "usage", None)

        if self.verbose or not raw.strip():
            print(f"\n  ┌── VLM (Kimi) raw response ───────────────────────────────────")
            print(f"  │  finish_reason : {finish}")
            if usage:
                print(f"  │  tokens: prompt={getattr(usage,'prompt_tokens','?')}  completion={getattr(usage,'completion_tokens','?')}")
            print(f"  │  content length: {len(raw)} chars")
            for line in raw.splitlines()[:20]:
                print(f"  │  {line}")
            if len(raw.splitlines()) > 20:
                print(f"  │  ... ({len(raw.splitlines())} lines total)")
            print(f"  └─────────────────────────────────────────────────────────────\n")

        if not raw.strip():
            raise ValueError(
                f"Kimi returned empty content (finish_reason={finish}). "
                "Check model vision support and API quota."
            )

        issues = self._parse_json(raw)
        print(f"  [VLM] Kimi found {len(issues)} issue(s).")
        if self.verbose:
            if issues:
                for iss in issues:
                    sev = iss.get('severity','?').upper()
                    msg = iss.get('message','')
                    print(f"    │  [vlm] {sev}: {msg}")
            else:
                print(f"    │  (VLM returned empty or could not parse issues from response)")
        return issues

    # ── GLM-4V fallback ────────────────────────────────────────────────────────

    def _critique_glm4v(self, existing: list[tuple], description: str) -> list[dict]:
        content = []
        for path, stem in existing:
            label = _VIEW_LABELS.get(stem, stem)
            content.append({"type": "text", "text": f"[{label}]"})
            content.append({
                "type": "image_url",
                "image_url": {"url": f"data:image/png;base64,{_b64(path)}"},
            })
        content.append({"type": "text", "text": _USER_TMPL.format(description=description)})

        response = self._get_glm().chat.completions.create(
            model="glm-4v-flash",
            messages=[
                {"role": "system", "content": _SYSTEM},
                {"role": "user",   "content": content},
            ],
            temperature=0.3,
            max_tokens=1024,
        )
        raw = response.choices[0].message.content
        if self.verbose:
            print(f"\n  ┌── VLM (GLM-4V) raw response ─────────────────────────────────")
            for line in raw.splitlines()[:20]:
                print(f"  │  {line}")
            print(f"  └─────────────────────────────────────────────────────────────\n")
        issues = self._parse_json(raw)
        print(f"  [VLM/GLM-4V] Found {len(issues)} issue(s).")
        if self.verbose:
            if issues:
                for iss in issues:
                    sev = iss.get('severity','?').upper()
                    msg = iss.get('message','')
                    print(f"    │  [vlm] {sev}: {msg}")
            else:
                print(f"    │  (VLM returned empty or unparseable response)")
        return issues

    # ── JSON parsing ───────────────────────────────────────────────────────────

    def _parse_json(self, raw: str) -> list[dict]:
        raw = raw.strip()
        raw = re.sub(r"^```(?:json)?\n?", "", raw, flags=re.MULTILINE)
        raw = re.sub(r"\n?```$", "", raw, flags=re.MULTILINE)
        try:
            data = json.loads(raw)
            if isinstance(data, list):
                for item in data:
                    item.setdefault("source", "vlm")
                    item.setdefault("severity", "warning")
                    item.setdefault("kind", "parameter")
                return data
        except json.JSONDecodeError:
            if raw.strip():
                return [{"source": "vlm", "severity": "warning",
                         "kind": "parameter", "message": raw.strip()}]
        return []
