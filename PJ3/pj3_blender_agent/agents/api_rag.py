"""bpy API-doc retriever (LL3M BlenderRAG style, sourced from RNA introspection).

Knowledge base: knowledge/bpy_api.jsonl, built by build_bpy_kb.py from the installed
Blender via introspection (exact, offline, version-matched). Each line is
    {"name", "kind", "tags", "text"}.

Retrieval: bge-small-en-v1.5 embeddings + cosine (cached to disk), with an English
keyword-overlap fallback when sentence-transformers is unavailable.

Used at two stages (mirroring LL3M):
  - Coder: retrieve(component/description text) before generating code.
  - Fixer: retrieve_for_error(traceback) to fetch the right API signatures on a crash.
"""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

_KB_PATH = Path(__file__).parent.parent / "knowledge" / "bpy_api.jsonl"
_EMB_PATH = Path(__file__).parent.parent / "knowledge" / "bpy_api_emb.npy"
_EMB_META = Path(__file__).parent.parent / "knowledge" / "bpy_api_emb.meta"
_EN_MODEL = "BAAI/bge-small-en-v1.5"

# Pull bpy identifiers out of a traceback / error string so the query matches the KB.
_API_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

# Generic words that appear in almost every traceback/entry and carry no signal about
# WHICH api is meant — excluded from keyword scoring so high-signal identifiers win.
_STOPWORDS = {
    "object", "type", "value", "default", "default_value", "self", "none", "true",
    "false", "name", "data", "bpy", "ops", "attribute", "error", "line", "module",
    "key", "has", "no", "not", "found", "got", "the", "for", "and", "with", "in",
    "is", "of", "to", "an", "argument", "keyword", "unexpected", "traceback",
    "most", "recent", "call", "last", "file", "context", "active", "new", "add",
}


class ApiDocRetriever:
    """Retrieves bpy API-doc chunks for a query or an error message."""

    def __init__(self, kb_path: str | Path | None = None, model_name: str = _EN_MODEL):
        self._kb_path = Path(kb_path) if kb_path else _KB_PATH
        self._model_name = model_name
        self._entries: list[dict] | None = None
        self._embeddings: "np.ndarray | None" = None
        self._model = None
        self._loaded = False

    # ── loading / indexing ───────────────────────────────────────────────────

    def _load(self) -> None:
        if self._loaded:
            return
        self._loaded = True

        if not self._kb_path.exists():
            print(f"  [API-RAG] knowledge base not found: {self._kb_path} "
                  f"(run build_bpy_kb.py); retrieval disabled.")
            self._entries = []
            return

        self._entries = [json.loads(l) for l in self._kb_path.read_text(encoding="utf-8").splitlines() if l.strip()]
        print(f"  [API-RAG] Loaded {len(self._entries)} API entries.")

        texts = [f"{e['name']} {e.get('tags','')} {e['text']}" for e in self._entries]
        digest = hashlib.md5(("\n".join(texts) + self._model_name).encode("utf-8")).hexdigest()

        # Try cached embeddings first.
        try:
            import numpy as np
            if _EMB_PATH.exists() and _EMB_META.exists() and _EMB_META.read_text().strip() == digest:
                self._embeddings = np.load(_EMB_PATH)
                self._model = self._lazy_model()  # still needed to embed queries
                if self._model is not None:
                    print(f"  [API-RAG] Loaded cached embeddings ({self._embeddings.shape[0]}).")
                    return
        except Exception:
            pass

        # Build embeddings.
        self._model = self._lazy_model()
        if self._model is not None:
            try:
                import numpy as np
                self._embeddings = self._model.encode(texts, normalize_embeddings=True, show_progress_bar=False)
                np.save(_EMB_PATH, self._embeddings)
                _EMB_META.write_text(digest)
                print(f"  [API-RAG] Built + cached embeddings with {self._model_name}.")
            except Exception as exc:
                print(f"  [API-RAG] embedding build failed ({exc}); using keyword fallback.")
                self._model = None
                self._embeddings = None

    def _lazy_model(self):
        try:
            from sentence_transformers import SentenceTransformer
            return SentenceTransformer(self._model_name)
        except Exception as exc:
            print(f"  [API-RAG] sentence-transformers unavailable ({exc}); keyword fallback.")
            return None

    # ── scoring ──────────────────────────────────────────────────────────────

    def _keyword_rank(self, query: str, top_k: int) -> list[int]:
        q_tokens = {t.lower() for t in _API_TOKEN.findall(query)} - _STOPWORDS
        scored = []
        for idx, e in enumerate(self._entries):
            hay = (e["name"] + " " + e.get("tags", "") + " " + e["text"]).lower()
            # name hits weigh more
            name_l = e["name"].lower()
            score = sum(2.0 for t in q_tokens if t in name_l) + sum(1.0 for t in q_tokens if t in hay)
            scored.append((score, idx))
        scored.sort(reverse=True)
        return [i for s, i in scored[:top_k] if s > 0]

    def _rank(self, query: str, top_k: int) -> list[int]:
        if self._model is not None and self._embeddings is not None:
            import numpy as np
            q = self._model.encode([query], normalize_embeddings=True)
            scores = (self._embeddings @ q.T).flatten()
            return [int(i) for i in scores.argsort()[::-1][:top_k]]
        return self._keyword_rank(query, top_k)

    # ── public API ─────────────────────────────────────────────────────────────

    def retrieve(self, query: str, top_k: int = 4) -> list[str]:
        """Return formatted API-doc chunks most relevant to *query*."""
        self._load()
        if not self._entries:
            return []
        return [self._entries[i]["text"] for i in self._rank(query, top_k)]

    @staticmethod
    def _focus_error(error_text: str) -> str:
        """Distil a noisy traceback into the signal that identifies the failing API.

        Keeps the offending source line(s) + the exception message, and amplifies the
        quoted socket names ("Specular"), bpy.* dotted paths, and keyword args (radius1=)
        that pin down the exact operator/modifier/socket. Drops 'File ...'/'Traceback'
        boilerplate that otherwise dilutes the match.
        """
        keep = []
        for ln in error_text.splitlines():
            s = ln.strip()
            if not s or s.startswith("Traceback") or s.startswith("File "):
                continue
            keep.append(s)
        focus = "\n".join(keep)
        quoted = re.findall(r"[\"']([A-Za-z][A-Za-z0-9 _]*)[\"']", error_text)
        dotted = re.findall(r"bpy\.[A-Za-z0-9_.]+", error_text)
        kwargs = re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*=", error_text)
        # Repeat the high-signal tokens so they dominate both keyword and semantic scores.
        amp = " ".join(quoted * 3 + dotted * 2 + kwargs * 2)
        return f"{focus}\n{amp}"

    def retrieve_for_error(self, error_text: str, top_k: int = 4) -> list[str]:
        """Return API-doc chunks relevant to a traceback / runtime error.

        Uses a focused query (see _focus_error) and a HYBRID rank: keyword/identifier
        matches first (precise for things like 'Specular' or 'radius1'), then semantic
        neighbours, merged and de-duplicated.
        """
        self._load()
        if not self._entries:
            return []
        query = self._focus_error(error_text)
        kw = self._keyword_rank(query, top_k)
        order = list(kw)
        if self._model is not None and self._embeddings is not None:
            for i in self._rank(query, top_k):
                if i not in order:
                    order.append(i)
        return [self._entries[i]["text"] for i in order[:top_k]]

    # Map a plan component's declared shape to one or more retrieval queries. A
    # "curve" part can be built two ways — a lathed Screw profile OR a hand-built
    # Bezier spline — so it pulls docs for BOTH (the bezier query targets the new
    # bpy.types.* data entries that carry the .co/handle arities).
    _SHAPE_QUERY = {
        "cone": ["primitive_cone_add cone frustum tapered radius1 radius2"],
        "torus": ["primitive_torus_add torus ring"],
        "sphere": ["primitive_uv_sphere_add sphere ball"],
        "curve": [
            "ScrewModifier lathe revolve profile around axis",
            "bpy.types.Curve Spline BezierSplinePoint splines.new bezier_points handle_left_type co bevel_depth",
        ],
        "custom": [
            "ScrewModifier lathe revolve profile around axis",
            "bpy.types.Curve Spline BezierSplinePoint splines.new bezier_points handle_left_type co bevel_depth",
        ],
    }

    def retrieve_for_plan(self, plan: dict, top_k_per: int = 1, max_total: int = 6) -> list[str]:
        """Targeted retrieval driven by the plan's declared shapes / style notes.

        Only the ADVANCED constructs the plan actually uses (cone/curve/torus/sphere,
        bevel in style notes) pull docs — trivial cube/cylinder add no noise. The
        Principled BSDF socket list is always included because every script sets
        materials and the 4.x→5.x socket renames are a top crash source.
        """
        self._load()
        if not self._entries:
            return []
        queries: list[str] = []
        shapes = {str(c.get("shape", "")).lower() for c in plan.get("components", [])}
        for s in shapes:
            queries.extend(self._SHAPE_QUERY.get(s, []))
        notes = (str(plan.get("style_notes", "")) + " " + str(plan.get("description_en", ""))).lower()
        if any(k in notes for k in ("bevel", "chamfer", "rounded edge", "round edge", "soften")):
            queries.append("BevelModifier bevel chamfer edges segments width")
        # Always reinforce the version-correct material sockets.
        queries.append("ShaderNodeBsdfPrincipled material Base Color Roughness Metallic sockets")

        seen: set[str] = set()
        out: list[str] = []
        for q in queries:
            for i in self._rank(q, top_k_per):
                name = self._entries[i]["name"]
                if name in seen:
                    continue
                seen.add(name)
                out.append(self._entries[i]["text"])
                if len(out) >= max_total:
                    return out
        return out

    def as_prompt_block(self, chunks: list[str], header: str) -> str:
        """Wrap retrieved chunks in a labelled block for prompt injection."""
        if not chunks:
            return ""
        body = "\n\n".join(chunks)
        return f"{header}\n{body}\n"
