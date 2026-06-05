# TITLE: 明式方腿构造
# DESCRIPTION: 典型明式家具方腿，细长比例1:10，四腿对称排布，坐板高0.50m
# TAGS: 椅子,凳子,太师椅,官帽椅,椅腿,方腿,明式

import bpy

mat = bpy.data.materials.new("WoodLeg")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.40, 0.22, 0.08, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.65

# 明式方腿标准参数：截面40mm×40mm，高450mm
leg_w, leg_d, leg_h = 0.04, 0.04, 0.45
# 座面宽0.60m，进深0.50m → 腿内缩10mm
for i, (lx, ly) in enumerate([(0.27, 0.22), (0.27, -0.22), (-0.27, 0.22), (-0.27, -0.22)]):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(lx, ly, leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i + 1}"
    leg.scale = (leg_w / 2, leg_d / 2, leg_h / 2)
    leg.data.materials.append(mat)
