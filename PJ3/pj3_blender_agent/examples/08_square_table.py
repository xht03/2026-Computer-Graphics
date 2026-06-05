# TITLE: 八仙桌/方桌完整参考
# DESCRIPTION: 明式方桌，0.90m×0.90m×0.78m，厚桌面+四腿+牙板，典型比例
# TAGS: 方桌,八仙桌,桌子,桌面,四方桌,完整,明式,参考

import bpy

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new("TableWood")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.40, 0.20, 0.07, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.60

# 尺寸
TW, TD, TH = 0.90, 0.90, 0.78   # 宽×深×高
top_t = 0.05    # 桌面厚
leg_w = 0.06    # 腿截面
leg_h = TH - top_t  # 0.73m
apron_t = 0.03  # 牙板厚
apron_h = 0.08  # 牙板高

# COMPONENT: tabletop
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, TH - top_t / 2))
top = bpy.context.active_object
top.name = "tabletop"
top.scale = (TW / 2, TD / 2, top_t / 2)
top.data.materials.append(mat)

# COMPONENT: legs (内缩5cm)
inset = 0.05
lx, ly = TW / 2 - inset, TD / 2 - inset
for i, (x, y) in enumerate([(lx, ly), (lx, -ly), (-lx, ly), (-lx, -ly)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i + 1}"
    leg.scale = (leg_w / 2, leg_w / 2, leg_h / 2)
    leg.data.materials.append(mat)

# COMPONENT: aprons (牙板，4面)
apron_z = TH - top_t - apron_h / 2
for name, loc, sx, sy in [
    ("apron_front",  (0,   TD / 2, apron_z), TW / 2, apron_t / 2),
    ("apron_back",   (0,  -TD / 2, apron_z), TW / 2, apron_t / 2),
    ("apron_left",   (TW / 2, 0,  apron_z), apron_t / 2, TD / 2),
    ("apron_right",  (-TW / 2, 0, apron_z), apron_t / 2, TD / 2),
]:
    bpy.ops.mesh.primitive_cube_add(size=1, location=loc)
    a = bpy.context.active_object
    a.name = name
    a.scale = (sx, sy, apron_h / 2)
    a.data.materials.append(mat)
