# TITLE: 明式案几参考模型
# DESCRIPTION: 案几（条案），细长比例，1.0m×0.35m×0.85m，四腿内缩，有牙板
# TAGS: 案几,条案,画案,书案,供案,长案,明式,细腿

import bpy

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("AltarWood")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.36, 0.17, 0.05, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.58

# 案几比例特征：细长，面宽远大于面深
L, W, H = 1.00, 0.35, 0.85   # 长×宽×高
top_t   = 0.04   # 面板厚（比桌薄）
leg_w   = 0.04   # 腿截面（比桌腿细）
leg_h   = H - top_t  # 0.81m
apron_h = 0.07

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, H - top_t / 2))
top = bpy.context.active_object
top.name = "top_board"
top.scale = (L / 2, W / 2, top_t / 2)
top.data.materials.append(mat)

# COMPONENT: legs (内缩4cm)
inset = 0.04
lx = L / 2 - inset
ly = W / 2 - inset
for i, (x, y) in enumerate([(lx, ly), (lx, -ly), (-lx, ly), (-lx, -ly)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i + 1}"
    leg.scale = (leg_w / 2, leg_w / 2, leg_h / 2)
    leg.data.materials.append(mat)

# COMPONENT: aprons (仅长边两侧有牙板)
apron_z = H - top_t - apron_h / 2
for name, loc, sx, sy in [
    ("apron_left",  (0,  W / 2, apron_z), L / 2, 0.015),
    ("apron_right", (0, -W / 2, apron_z), L / 2, 0.015),
]:
    bpy.ops.mesh.primitive_cube_add(size=1, location=loc)
    a = bpy.context.active_object
    a.name = name
    a.scale = (sx, sy, apron_h / 2)
    a.data.materials.append(mat)
