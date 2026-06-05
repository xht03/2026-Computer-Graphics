"""Domain RAG Agent — Day 5 implementation.

Retrieves relevant bpy code snippets from the Chinese furniture example library.

Architecture:
  - Primary:  sentence-transformers + BAAI/bge-small-zh-v1.5 embeddings,
              cosine similarity (numpy, no FAISS needed at ≤ 20 snippets).
  - Fallback: Chinese character-overlap keyword scoring when the model is
              unavailable (no GPU / not installed).

Each snippet file in examples/ carries header metadata:
    # TITLE: <title>
    # DESCRIPTION: <description>
    # TAGS: <tag1>,<tag2>,...
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

_EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
_BGE_MODEL    = "BAAI/bge-small-zh-v1.5"


# ── Snippet parser ──────────────────────────────────────────────────────────

def _parse_snippet(filepath: Path) -> dict | None:
    """Parse a snippet file with structured header comments.

    Returns a dict with keys: title, description, tags, code, file.
    Returns None if the file has no TITLE header.
    """
    text = filepath.read_text(encoding="utf-8")

    title_m = re.search(r"^#\s*TITLE:\s*(.+)$",       text, re.MULTILINE)
    desc_m  = re.search(r"^#\s*DESCRIPTION:\s*(.+)$", text, re.MULTILINE)
    tags_m  = re.search(r"^#\s*TAGS:\s*(.+)$",        text, re.MULTILINE)

    if not title_m:
        return None

    # Find the end of all header lines, then take the remaining text as code
    header_end = max(
        (m.end() for m in (title_m, desc_m, tags_m) if m),
        default=0,
    )
    # Skip any "# ──..." separator lines right after the headers
    code_raw = text[header_end:]
    code = re.sub(r"^\s*#[─\-=]{3,}.*\n", "", code_raw).lstrip("\n")

    return {
        "title":       title_m.group(1).strip(),
        "description": desc_m.group(1).strip() if desc_m else "",
        "tags":        [t.strip() for t in tags_m.group(1).split(",")] if tags_m else [],
        "code":        code,
        "file":        filepath.name,
    }


# ── RAGRetriever ────────────────────────────────────────────────────────────

class RAGRetriever:
    """Retrieves relevant bpy snippets using BGE embeddings + cosine similarity.

    Lazy-loads the model and builds the index on the first call to retrieve().
    Falls back to Chinese character-overlap scoring if sentence-transformers
    is unavailable.
    """

    def __init__(
        self,
        examples_dir: str | Path | None = None,
        model_name: str = _BGE_MODEL,
    ) -> None:
        self._examples_dir = Path(examples_dir) if examples_dir else _EXAMPLES_DIR
        self._model_name   = model_name
        self._snippets: list[dict] | None = None
        self._embeddings: "np.ndarray | None" = None
        self._model = None
        self._loaded = False

    # ── Internal helpers ────────────────────────────────────────────────────

    def _load(self) -> None:
        if self._loaded:
            return
        self._loaded = True

        # Gather snippet files
        self._snippets = []
        for fp in sorted(self._examples_dir.glob("*.py")):
            s = _parse_snippet(fp)
            if s:
                self._snippets.append(s)

        if not self._snippets:
            print("  [RAG] Warning: no snippets found in examples/")
            return

        print(f"  [RAG] Loaded {len(self._snippets)} snippet(s) from {self._examples_dir}")

        # Build BGE embeddings
        try:
            from sentence_transformers import SentenceTransformer
            self._model = SentenceTransformer(self._model_name)
            texts = [
                f"{s['title']} {s['description']} {' '.join(s['tags'])}"
                for s in self._snippets
            ]
            self._embeddings = self._model.encode(texts, normalize_embeddings=True)
            print(f"  [RAG] Embeddings built with {self._model_name}")
        except Exception as exc:
            print(f"  [RAG] sentence-transformers unavailable ({exc}); using keyword fallback.")
            self._model     = None
            self._embeddings = None

    def _keyword_score(self, query: str, snippet: dict) -> float:
        """Character-overlap score for Chinese keyword fallback."""
        text = " ".join([snippet["title"], snippet["description"]] + snippet["tags"])
        # Count unique Chinese character matches as proxy for relevance
        return sum(
            1.0 for ch in set(query)
            if "一" <= ch <= "鿿" and ch in text
        )

    # ── Public API ──────────────────────────────────────────────────────────

    def retrieve(self, query: str, top_k: int = 3) -> list[str]:
        """Return top_k relevant bpy code snippets for *query*.

        Each returned snippet is prefixed with a comment header so the coder
        agent can see the title/description as context.

        Returns an empty list if no snippets are available.
        """
        self._load()
        if not self._snippets:
            return []

        if self._model is not None and self._embeddings is not None:
            import numpy as np
            q_emb = self._model.encode([query], normalize_embeddings=True)
            scores: np.ndarray = (self._embeddings @ q_emb.T).flatten()
            top_idx = scores.argsort()[::-1][:top_k]
            top_snippets = [self._snippets[int(i)] for i in top_idx]
        else:
            # Keyword fallback
            scored = sorted(
                self._snippets,
                key=lambda s: self._keyword_score(query, s),
                reverse=True,
            )
            top_snippets = scored[:top_k]

        # Format: header comment + code
        return [
            f"# === {s['title']} ===\n# {s['description']}\n{s['code']}"
            for s in top_snippets
        ]

    def retrieve_with_meta(self, query: str, top_k: int = 3) -> list[dict]:
        """Like retrieve() but returns full snippet dicts (with scores if available)."""
        self._load()
        if not self._snippets:
            return []

        if self._model is not None and self._embeddings is not None:
            import numpy as np
            q_emb = self._model.encode([query], normalize_embeddings=True)
            scores: np.ndarray = (self._embeddings @ q_emb.T).flatten()
            top_idx = scores.argsort()[::-1][:top_k]
            return [
                {**self._snippets[int(i)], "score": float(scores[int(i)])}
                for i in top_idx
            ]
        else:
            scored = sorted(
                self._snippets,
                key=lambda s: self._keyword_score(query, s),
                reverse=True,
            )
            return [
                {**s, "score": self._keyword_score(query, s)}
                for s in scored[:top_k]
            ]
