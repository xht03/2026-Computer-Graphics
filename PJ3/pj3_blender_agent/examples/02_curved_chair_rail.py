# TITLE: 弧形椅圈（圈椅/官帽椅扶手背圈）
# DESCRIPTION: 用多段矩形截面近似弧形椅圈，从扶手延续到椅背，高出座面0.30~0.35m
# TAGS: 圈椅,官帽椅,椅圈,扶手,椅背,弧形,靠背

import bpy
import math

mat = bpy.data.materials.new("WoodRail")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.42, 0.23, 0.09, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.60

# 椅圈用12段折线近似半圆弧，圆心在座面中心上方
# 半径0.30m，圆心高度0.80m（从地面算）
seat_h = 0.50   # 座面高度
rail_r = 0.30   # 椅圈半径（水平投影）
rail_z = 0.80   # 椅圈顶端高度
n_segs = 12     # 分段数
rail_w = 0.04   # 截面宽度
rail_h = 0.04   # 截面高度

for i in range(n_segs):
    # 从正前方(angle=0)到正后方(angle=π)
    a0 = math.pi * i / n_segs
    a1 = math.pi * (i + 1) / n_segs
    # 两端点坐标（前：y负；后：y正）
    x0 = rail_r * math.sin(a0)
    y0 = -rail_r * math.cos(a0)
    x1 = rail_r * math.sin(a1)
    y1 = -rail_r * math.cos(a1)
    # 段中点
    mx = (x0 + x1) / 2
    my = (y0 + y1) / 2
    # 线段长度
    seg_len = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    # 高度：两端等高（椅圈在一个平面内，略前低后高）
    mz = seat_h + (rail_z - seat_h) * (i + 0.5) / n_segs

    bpy.ops.mesh.primitive_cube_add(size=1, location=(mx, my, mz))
    seg = bpy.context.active_object
    seg.name = f"rail_seg_{i}"
    seg.scale = (seg_len / 2 + 0.002, rail_w / 2, rail_h / 2)
    # 旋转使长轴对齐线段方向
    angle = math.atan2(y1 - y0, x1 - x0)
    seg.rotation_euler = (0, 0, angle)
    seg.data.materials.append(mat)
