# TITLE: 三扇屏风框架构造
# DESCRIPTION: 中式三折屏风，每扇0.50m宽1.80m高，框架+内填板双层结构
# TAGS: 屏风,折屏,框架,隔断,三折,装饰

import bpy

mat_frame = bpy.data.materials.new("ScreenFrame")
mat_frame.use_nodes = True
mat_frame.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.35, 0.18, 0.06, 1.0)
mat_frame.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.70

mat_fill = bpy.data.materials.new("ScreenFill")
mat_fill.use_nodes = True
mat_fill.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.90, 0.82, 0.65, 1.0)
mat_fill.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.80

n_panels = 3
panel_w = 0.50
panel_d = 0.05   # 框架厚度
panel_h = 1.80
gap = 0.02       # 扇间隙

for i in range(n_panels):
    px = (i - (n_panels - 1) / 2.0) * (panel_w + gap)

    # 外框
    bpy.ops.mesh.primitive_cube_add(size=1, location=(px, 0.0, panel_h / 2))
    frame = bpy.context.active_object
    frame.name = f"frame_{i + 1}"
    frame.scale = (panel_w / 2, panel_d / 2, panel_h / 2)
    frame.data.materials.append(mat_frame)

    # 内芯（薄板，比外框略小）
    bpy.ops.mesh.primitive_cube_add(size=1, location=(px, 0.0, panel_h / 2))
    inner = bpy.context.active_object
    inner.name = f"panel_{i + 1}"
    inner.scale = (panel_w * 0.85 / 2, panel_d * 0.25 / 2, panel_h * 0.88 / 2)
    inner.data.materials.append(mat_fill)
