# TITLE: 椅背构造（靠背板+搭脑）
# DESCRIPTION: 椅背由搭脑（顶横梁）+竖背板+靠背立柱组成，背高0.60m
# TAGS: 椅背,靠背,搭脑,背板,太师椅,官帽椅,背柱

import bpy

mat = bpy.data.materials.new("BackMat")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.38, 0.19, 0.06, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.62

seat_h  = 0.50    # 座面高度
seat_w  = 0.60    # 座面宽
back_h  = 0.60    # 靠背高（座面以上）
back_z0 = seat_h  # 靠背底端高度
back_w  = 0.46    # 靠背内宽
back_d  = 0.04    # 靠背厚度

# COMPONENT: back_post_left（左靠背立柱，从座面延伸到搭脑）
post_h = back_h + 0.05
bpy.ops.mesh.primitive_cube_add(size=1, location=(back_w / 2, -0.02, back_z0 + post_h / 2))
post_l = bpy.context.active_object
post_l.name = "back_post_left"
post_l.scale = (0.025, 0.025, post_h / 2)
post_l.data.materials.append(mat)

# COMPONENT: back_post_right
bpy.ops.mesh.primitive_cube_add(size=1, location=(-back_w / 2, -0.02, back_z0 + post_h / 2))
post_r = bpy.context.active_object
post_r.name = "back_post_right"
post_r.scale = (0.025, 0.025, post_h / 2)
post_r.data.materials.append(mat)

# COMPONENT: back_panel（背板，稍内缩）
panel_z = back_z0 + back_h * 0.45
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.02, panel_z))
panel = bpy.context.active_object
panel.name = "back_panel"
panel.scale = (back_w * 0.85 / 2, 0.015, back_h * 0.7 / 2)
panel.data.materials.append(mat)

# COMPONENT: top_rail（搭脑，顶部横梁，略宽于靠背）
top_z = back_z0 + post_h - 0.025
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.02, top_z))
top_rail = bpy.context.active_object
top_rail.name = "top_rail"
top_rail.scale = ((back_w + 0.08) / 2, 0.035, 0.025)
top_rail.data.materials.append(mat)
