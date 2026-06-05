# TITLE: 束腰牙板（桌裙/椅裙）
# DESCRIPTION: 桌面或坐板下方的束腰和牙板装饰构件，常见于明式桌椅，厚0.03m
# TAGS: 牙板,束腰,桌裙,椅裙,案几,方桌,八仙桌,装饰

import bpy

mat = bpy.data.materials.new("ApronMat")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.40, 0.20, 0.07, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.62

# 桌面参数（束腰紧贴桌面底部）
table_w = 0.90   # 桌面宽度
table_d = 0.90   # 桌面进深
top_h   = 0.78   # 桌面顶部高度
top_t   = 0.05   # 桌面板厚度
apron_t = 0.03   # 牙板厚度
apron_h = 0.08   # 牙板高度
apron_z = top_h - top_t - apron_h / 2   # 牙板中心高度

# 前牙板（沿 X 方向）
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, table_d / 2, apron_z))
apron_front = bpy.context.active_object
apron_front.name = "apron_front"
apron_front.scale = (table_w / 2, apron_t / 2, apron_h / 2)
apron_front.data.materials.append(mat)

# 后牙板
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -table_d / 2, apron_z))
apron_back = bpy.context.active_object
apron_back.name = "apron_back"
apron_back.scale = (table_w / 2, apron_t / 2, apron_h / 2)
apron_back.data.materials.append(mat)

# 左牙板（沿 Y 方向）
bpy.ops.mesh.primitive_cube_add(size=1, location=(table_w / 2, 0, apron_z))
apron_left = bpy.context.active_object
apron_left.name = "apron_left"
apron_left.scale = (apron_t / 2, table_d / 2, apron_h / 2)
apron_left.data.materials.append(mat)

# 右牙板
bpy.ops.mesh.primitive_cube_add(size=1, location=(-table_w / 2, 0, apron_z))
apron_right = bpy.context.active_object
apron_right.name = "apron_right"
apron_right.scale = (apron_t / 2, table_d / 2, apron_h / 2)
apron_right.data.materials.append(mat)
