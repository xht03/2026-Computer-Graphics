# TITLE: 长凳完整参考模型
# DESCRIPTION: 明式长凳完整构造：面板+四腿+横枨，尺寸1.20m×0.30m×0.45m
# TAGS: 长凳,凳子,坐具,明式,完整,参考

import bpy, math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("BenchWood")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.42, 0.23, 0.09, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.63

# 整体尺寸
L, W, H = 1.20, 0.30, 0.45   # 长×宽×高(m)
top_t   = 0.05   # 面板厚
leg_w   = 0.05   # 腿截面
leg_h   = H - top_t  # 腿高 = 0.40m

# COMPONENT: seat_top
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, H - top_t / 2))
top = bpy.context.active_object
top.name = "seat_top"
top.scale = (L / 2, W / 2, top_t / 2)
top.data.materials.append(mat)

# COMPONENT: legs (四腿，内缩3cm)
inset = 0.03
lx, ly = L / 2 - inset, W / 2 - inset
for i, (x, y) in enumerate([(lx, ly), (lx, -ly), (-lx, ly), (-lx, -ly)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i + 1}"
    leg.scale = (leg_w / 2, leg_w / 2, leg_h / 2)
    leg.data.materials.append(mat)

# COMPONENT: stretchers (横枨，腿高1/3处)
bar_z = leg_h / 3
bar_r = 0.015
# 前后横枨
for name, y0 in [("str_front", ly), ("str_back", -ly)]:
    bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=lx * 2, location=(0, y0, bar_z))
    b = bpy.context.active_object
    b.name = name
    b.rotation_euler = (0, math.radians(90), 0)
    b.data.materials.append(mat)
# 两侧横枨
for name, x0 in [("str_left", lx), ("str_right", -lx)]:
    bpy.ops.mesh.primitive_cylinder_add(radius=bar_r, depth=ly * 2, location=(x0, 0, bar_z))
    b = bpy.context.active_object
    b.name = name
    b.rotation_euler = (math.radians(90), 0, 0)
    b.data.materials.append(mat)
