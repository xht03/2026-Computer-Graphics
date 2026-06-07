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
    # cylinder 不支持 radius1/radius2（只有 radius）；但 cone 合法使用它们做圆台，故仅拦 cylinder 行。
    re.compile(r"cylinder_add.*radius1\s*="),
    re.compile(r"cylinder_add.*radius2\s*="),
    re.compile(r"\.location\.rotate\("),               # Vector.rotate() 用法错误
    re.compile(r"\.rotate\(\s*\("),                    # rotate((轴), 角度) —— 签名错误
    re.compile(r"use_nodes\s*=\s*False"),              # 禁止关闭节点，否则材质失效
    # 注：from_pydata 曾被全面禁用，但 lathe/回转体 recipe 需要它建轮廓线。
    # 现仅靠 prompt 约束（faces 必须为 []，只用于 screw 轮廓）来防崩，不再正则拦截。
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

    def generate(self, description: str, api_docs: list[str] | None = None) -> str:
        """为给定的 *description* 生成 bpy 代码。

        Args:
            description: 自然语言物体描述（中文或英文均可）。
            api_docs: 可选，针对计划用到的高级构造检索到的 bpy API 文档片段
                      （来自 ApiDocRetriever，内省本机 Blender 的精确签名/socket）。

        Returns:
            可直接执行的 bpy Python 代码字符串。
        """
        if api_docs:
            api_block = (
                "Blender API reference for THIS version — use these EXACT "
                "signatures, parameters, and socket/property names:\n\n"
                + "\n\n".join(api_docs)
                + "\n"
            )
        else:
            api_block = ""

        user_msg = self._user_tmpl.format(
            description=description,
            api_docs=api_block,
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
