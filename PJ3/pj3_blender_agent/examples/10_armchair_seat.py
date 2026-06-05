# TITLE: 太师椅坐板与扶手
# DESCRIPTION: 太师椅坐板0.60m×0.55m，高0.50m，带厚实扶手（宽0.08m高0.20m）
# TAGS: 太师椅,扶手椅,坐板,座面,扶手,椅子,坐具

import bpy

mat = bpy.data.materials.new("ChairWood")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.35, 0.17, 0.05, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.60

seat_w = 0.60   # 坐板宽
seat_d = 0.55   # 坐板深
seat_h = 0.50   # 座面高
seat_t = 0.06   # 坐板厚

arm_w  = 0.08   # 扶手截面宽
arm_h  = 0.06   # 扶手截面高
arm_z  = seat_h + 0.22  # 扶手高度（座面以上0.22m）
arm_len = seat_d * 0.85  # 扶手长度

# COMPONENT: seat_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, seat_h - seat_t / 2))
seat = bpy.context.active_object
seat.name = "seat_board"
seat.scale = (seat_w / 2, seat_d / 2, seat_t / 2)
seat.data.materials.append(mat)

# COMPONENT: armrest_left
bpy.ops.mesh.primitive_cube_add(size=1, location=(seat_w / 2 - arm_w / 2, 0, arm_z))
arm_l = bpy.context.active_object
arm_l.name = "armrest_left"
arm_l.scale = (arm_w / 2, arm_len / 2, arm_h / 2)
arm_l.data.materials.append(mat)

# COMPONENT: armrest_right
bpy.ops.mesh.primitive_cube_add(size=1, location=(-(seat_w / 2 - arm_w / 2), 0, arm_z))
arm_r = bpy.context.active_object
arm_r.name = "armrest_right"
arm_r.scale = (arm_w / 2, arm_len / 2, arm_h / 2)
arm_r.data.materials.append(mat)

# COMPONENT: armrest_support_left（扶手立柱）
support_h = arm_z - seat_h
bpy.ops.mesh.primitive_cube_add(size=1, location=(seat_w / 2 - arm_w / 2, -seat_d / 2 + 0.06, seat_h + support_h / 2))
sup_l = bpy.context.active_object
sup_l.name = "arm_support_left"
sup_l.scale = (arm_w / 2, 0.03, support_h / 2)
sup_l.data.materials.append(mat)

# COMPONENT: armrest_support_right
bpy.ops.mesh.primitive_cube_add(size=1, location=(-(seat_w / 2 - arm_w / 2), -seat_d / 2 + 0.06, seat_h + support_h / 2))
sup_r = bpy.context.active_object
sup_r.name = "arm_support_right"
sup_r.scale = (arm_w / 2, 0.03, support_h / 2)
sup_r.data.materials.append(mat)
