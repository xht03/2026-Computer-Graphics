#!/usr/bin/env python3
"""Build the bpy API-doc knowledge base by INTROSPECTING the installed Blender.

Run inside Blender:
    blender --background --python build_bpy_kb.py

Writes knowledge/bpy_api.jsonl — one JSON object per line:
    {"name": <api path>, "kind": operator|modifier|shader_node|note, "tags": "...", "text": "..."}

This is the LL3M-BlenderRAG idea, but sourced from RNA introspection (exact, offline,
version-matched to THIS Blender) rather than scraped HTML. It captures the things that
actually cause crashes: valid operator parameters, modifier property names, and the
Blender-4.x shader socket names.
"""
import bpy
import json
import os

_OUT = os.path.join(os.path.dirname(__file__), "knowledge", "bpy_api.jsonl")

# Operator namespaces we actually use when modeling furniture/objects.
_OP_NAMESPACES = ["mesh", "object", "curve", "transform", "material"]

# Shader nodes whose dynamic socket names matter (Principled is the crash hotspot).
_SHADER_NODES = [
    "ShaderNodeBsdfPrincipled", "ShaderNodeOutputMaterial", "ShaderNodeTexNoise",
    "ShaderNodeTexChecker", "ShaderNodeTexWave", "ShaderNodeValToRGB",
    "ShaderNodeTexCoord", "ShaderNodeMapping", "ShaderNodeBump", "ShaderNodeMixRGB",
    "ShaderNodeEmission", "ShaderNodeBsdfDiffuse", "ShaderNodeBsdfGlass",
]


def _prop_str(p):
    """One-line description of an RNA property."""
    if p.type == "ENUM":
        vals = "|".join(e.identifier for e in p.enum_items)
        s = f"{p.identifier} (ENUM in {vals})"
    elif p.type in ("FLOAT", "INT") and getattr(p, "array_length", 0) > 1:
        s = f"{p.identifier} ({p.type}[{p.array_length}])"
    else:
        s = f"{p.identifier} ({p.type})"
    desc = getattr(p, "description", "") or ""
    return s + (f" - {desc}" if desc else "")


def build_operators(entries):
    for ns in _OP_NAMESPACES:
        mod = getattr(bpy.ops, ns, None)
        if mod is None:
            continue
        for op_name in dir(mod):
            if op_name.startswith("_"):
                continue
            try:
                op = getattr(mod, op_name)
                rna = op.get_rna_type()
            except Exception:
                continue
            params = [p for p in rna.properties if p.identifier != "rna_type"]
            sig = ", ".join(p.identifier for p in params)
            lines = [f"bpy.ops.{ns}.{op_name}({sig})"]
            opdesc = getattr(rna, "description", "") or ""
            if opdesc:
                lines.append(f"  {opdesc}")
            for p in params:
                lines.append(f"    {_prop_str(p)}")
            entries.append({
                "name": f"bpy.ops.{ns}.{op_name}",
                "kind": "operator",
                "tags": f"{ns} operator {op_name}",
                "text": "\n".join(lines),
            })


def build_modifiers(entries):
    # Base props common to every modifier — exclude so each entry shows only its specifics.
    base = {p.identifier for p in bpy.types.Modifier.bl_rna.properties}
    # A scratch mesh object to instantiate each modifier type on.
    me = bpy.data.meshes.new("_kb_mesh")
    obj = bpy.data.objects.new("_kb_obj", me)
    bpy.context.collection.objects.link(obj)
    type_items = bpy.types.Modifier.bl_rna.properties["type"].enum_items
    for it in type_items:
        tstr = it.identifier  # e.g. 'SCREW'
        try:
            mod = obj.modifiers.new(name="m", type=tstr)
        except Exception:
            continue
        if mod is None:  # some types aren't valid on a mesh object
            continue
        specific = [p for p in mod.bl_rna.properties if p.identifier not in base]
        lines = [
            f"{mod.bl_rna.identifier}  —  obj.modifiers.new(name, type='{tstr}')",
        ]
        if it.description:
            lines.append(f"  {it.description}")
        if specific:
            lines.append("  Properties:")
            for p in specific:
                lines.append(f"    {_prop_str(p)}")
        entries.append({
            "name": f"modifier:{tstr}",
            "kind": "modifier",
            "tags": f"modifier {tstr} {mod.bl_rna.identifier}",
            "text": "\n".join(lines),
        })
        obj.modifiers.remove(mod)


def build_shader_nodes(entries):
    mat = bpy.data.materials.new("_kb_mat")
    mat.use_nodes = True
    tree = mat.node_tree
    for node_type in _SHADER_NODES:
        try:
            node = tree.nodes.new(node_type)
        except Exception:
            continue
        ins = [s.name for s in node.inputs]
        outs = [s.name for s in node.outputs]
        lines = [
            f'{node_type}  —  tree.nodes.new("{node_type}")',
            f"  Input sockets (use these EXACT names): {', '.join(ins)}",
            f"  Output sockets: {', '.join(outs)}",
        ]
        entries.append({
            "name": f"shader_node:{node_type}",
            "kind": "shader_node",
            "tags": f"shader node material {node_type}",
            "text": "\n".join(lines),
        })
        tree.nodes.remove(node)


def main():
    entries = []
    build_operators(entries)
    build_modifiers(entries)
    build_shader_nodes(entries)

    os.makedirs(os.path.dirname(_OUT), exist_ok=True)
    with open(_OUT, "w", encoding="utf-8") as f:
        for e in entries:
            f.write(json.dumps(e, ensure_ascii=False) + "\n")
    print(f"KB_BUILT: {len(entries)} entries -> {_OUT}")
    # breakdown
    from collections import Counter
    c = Counter(e["kind"] for e in entries)
    print("KB_BREAKDOWN:", dict(c))


if __name__ == "__main__":
    main()
