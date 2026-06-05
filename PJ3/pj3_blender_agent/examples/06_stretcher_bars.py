# TITLE: 横枨（腿间连接杆）
# DESCRIPTION: 连接桌腿或椅腿之间的横向支撑杆，增强结构强度，典型高度在腿高1/3处
# TAGS: 横枨,连接杆,支撑,桌腿,椅腿,加固,长凳,案几

import bpy

mat = bpy.data.materials.new("StretcherMat")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.38, 0.20, 0.07, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.65

# 以长凳为例：腿间距0.80m（X方向），进深0.28m（Y方向），腿高0.44m
leg_x = 0.42   # 腿的 X 坐标（±）
leg_y = 0.12   # 腿的 Y 坐标（±）
bar_z = 0.15   # 横枨高度（约腿高1/3处）
bar_r = 0.015  # 横枨半径

# 前横枨（连接两前腿，沿 X 轴）
front_len = leg_x * 2 + 0.04  # 略超出腿中心
bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=front_len, location=(0, leg_y, bar_z))
bar_f = bpy.context.active_object
bar_f.name = "stretcher_front"
import math
bar_f.rotation_euler = (0, math.radians(90), 0)
bar_f.data.materials.append(mat)

# 后横枨
bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=front_len, location=(0, -leg_y, bar_z))
bar_b = bpy.context.active_object
bar_b.name = "stretcher_back"
bar_b.rotation_euler = (0, math.radians(90), 0)
bar_b.data.materials.append(mat)

# 左侧横枨（连接左前腿和左后腿，沿 Y 轴）
side_len = leg_y * 2 + 0.04
bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=side_len, location=(leg_x, 0, bar_z))
bar_l = bpy.context.active_object
bar_l.name = "stretcher_left"
bar_l.rotation_euler = (math.radians(90), 0, 0)
bar_l.data.materials.append(mat)

# 右侧横枨
bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=side_len, location=(-leg_x, 0, bar_z))
bar_r_obj = bpy.context.active_object
bar_r_obj.name = "stretcher_right"
bar_r_obj.rotation_euler = (math.radians(90), 0, 0)
bar_r_obj.data.materials.append(mat)
