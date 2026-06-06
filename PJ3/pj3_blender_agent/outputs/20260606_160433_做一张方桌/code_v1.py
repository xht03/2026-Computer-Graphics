import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

mat = bpy.data.materials.new(name="WoodMat")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)
bsdf.inputs["Roughness"].default_value = 0.6

TH = 0.78
top_w, top_d, top_t = 0.9, 0.9, 0.06
apron_t = 0.02
leg_r = 0.09
leg_h = TH - top_t

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, TH - top_t / 2))
top = bpy.context.active_object
top.name = "top_board"
top.scale = (top_w / 2, top_d / 2, top_t / 2)
top.data.materials.append(mat)

# COMPONENT: apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, TH - top_t - apron_t / 2))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (top_w / 2, top_d / 2, apron_t / 2)
apron.data.materials.append(mat)

# COMPONENT: leg
inset = leg_r
lx_abs = top_w / 2 - inset
ly_abs = top_d / 2 - inset
for i, (lx, ly) in enumerate([(lx_abs, ly_abs), (lx_abs, -ly_abs), (-lx_abs, ly_abs), (-lx_abs, -ly_abs)]):
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_r, depth=leg_h, location=(lx, ly, leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.data.materials.append(mat)