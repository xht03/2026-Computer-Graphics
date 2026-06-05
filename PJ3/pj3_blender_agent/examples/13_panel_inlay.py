# TITLE: 嵌板装饰（框式嵌板结构）
# DESCRIPTION: 中式家具常见框式嵌板：外框+内嵌浅色面板，适用于椅背、桌腿间板、屏风
# TAGS: 嵌板,框架,装饰,椅背,面板,屏风,腿间,裙板

import bpy

mat_frame = bpy.data.materials.new("PanelFrame")
mat_frame.use_nodes = True
mat_frame.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.30, 0.14, 0.04, 1.0)
mat_frame.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.65

mat_inlay = bpy.data.materials.new("PanelInlay")
mat_inlay.use_nodes = True
mat_inlay.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.85, 0.74, 0.52, 1.0)
mat_inlay.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.75

# 单块嵌板参数：宽0.40m，高0.35m（可缩放复用）
p_w = 0.40   # 面板总宽
p_h = 0.35   # 面板总高
p_d = 0.04   # 面板总厚
frame_edge = 0.035  # 外框边宽

# 组合位置（可修改 loc_z 放置到椅背或桌裙位置）
loc_x, loc_y, loc_z = 0.0, 0.0, 0.20

# COMPONENT: panel_frame（外框）
bpy.ops.mesh.primitive_cube_add(size=1, location=(loc_x, loc_y, loc_z))
frame = bpy.context.active_object
frame.name = "panel_frame"
frame.scale = (p_w / 2, p_d / 2, p_h / 2)
frame.data.materials.append(mat_frame)

# COMPONENT: panel_inlay（内嵌板，稍薄且偏后）
bpy.ops.mesh.primitive_cube_add(size=1, location=(loc_x, loc_y - p_d * 0.15, loc_z))
inlay = bpy.context.active_object
inlay.name = "panel_inlay"
inlay.scale = ((p_w - frame_edge * 2) / 2,
               p_d * 0.3 / 2,
               (p_h - frame_edge * 2) / 2)
inlay.data.materials.append(mat_inlay)
