import os
import re
from pathlib import Path
from openai import OpenAI

_PROMPT_DIR = Path(__file__).parent.parent / "prompts"


def _load_prompt(name: str) -> str:
    return (_PROMPT_DIR / name).read_text(encoding="utf-8")


def _extract_code(text: str) -> str:
    """若 LLM 将输出包在 markdown 代码块中，则剥除围栏标记。"""
    text = text.strip()
    text = re.sub(r"^```(?:python)?\n?", "", text, flags=re.MULTILINE)
    text = re.sub(r"\n?```$", "", text, flags=re.MULTILINE)
    return text.strip()


# bpy 中不存在或会导致报错的 API 调用模式 —— 匹配到的行将被移除。
_BAD_LINE_PATTERNS = [
    re.compile(r"bpy\.ops\.object\.move_to_cursor"),
    re.compile(r"bpy\.ops\.object\.select_by_type"),
    re.compile(r"bpy\.ops\.object\.duplicate\(\)"),
    re.compile(r"radius1\s*="),                        # primitive_cylinder_add 不支持 radius1 参数
    re.compile(r"radius2\s*="),
    re.compile(r"\.location\.rotate\("),               # Vector.rotate() 用法错误
    re.compile(r"\.rotate\(\s*\("),                    # rotate((轴), 角度) —— 签名错误
    re.compile(r"use_nodes\s*=\s*False"),              # 禁止关闭节点，否则材质失效
    re.compile(r"\.from_pydata\("),                    # 顶点数据格式错误时必然崩溃
    re.compile(r"primitive_cube_add\s*\(.*size\s*=\s*0\."), # size=0.x 会生成微型方块，应为 size=1
    # 循环内对创建后物体的 scale/location 累加操作（"弯腿"反模式）。
    # 会破坏 leg_h = TH - top_t 的不变量，导致桌腿与桌面脱节。
    re.compile(r"\.\s*scale\s*\.\s*[xyz]\s*[+\-\*\/]="),    # 例：leg.scale.z -= leg.scale.z/4
    re.compile(r"\.\s*location\s*\.\s*[xyz]\s*[+\-\*\/]="), # 例：leg.location.z += ...
]


def sanitize_bpy_code(code: str) -> str:
    """移除已知会导致 Blender 崩溃的 bpy API 调用行。"""
    out_lines = []
    for line in code.splitlines():
        stripped = line.strip()
        # 若匹配到危险模式，则注释掉该行
        if any(p.search(stripped) for p in _BAD_LINE_PATTERNS):
            out_lines.append(f"# [sanitized] {line}")
            continue
        out_lines.append(line)
    return "\n".join(out_lines)


KIMI_BASE_URL = "https://api.moonshot.cn/v1"
KIMI_MODEL    = "kimi-k2.5"


class CoderAgent:
    """调用 Kimi k2.5，根据给定描述生成对应的 bpy 脚本。"""

    def __init__(self, model: str = KIMI_MODEL):
        self.client = OpenAI(
            api_key=os.getenv("KIMI_API_KEY"),
            base_url=KIMI_BASE_URL,
        )
        self.model = model
        self._system = _load_prompt("coder_system.txt")
        self._user_tmpl = _load_prompt("coder_user.txt")

    def generate(self, description: str, rag_snippets: list[str] | None = None) -> str:
        """为给定的 *description* 生成 bpy 代码。

        Args:
            description: 自然语言物体描述（中文或英文均可）。
            rag_snippets: 可选，RAG 检索到的 bpy 代码示例列表。

        Returns:
            可直接执行的 bpy Python 代码字符串。
        """
        if rag_snippets:
            rag_block = "Here are relevant bpy code examples for reference:\n\n"
            for i, snippet in enumerate(rag_snippets, 1):
                rag_block += f"--- Example {i} ---\n{snippet}\n\n"
        else:
            rag_block = ""

        user_msg = self._user_tmpl.format(
            description=description,
            rag_context=rag_block,
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
        raw = response.choices[0].message.content
        code = _extract_code(raw)
        return sanitize_bpy_code(code)
