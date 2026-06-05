# TITLE: 收腰锥形腿（用层叠薄片近似锥形）
# DESCRIPTION: 上宽下细的收腰腿，用多个尺寸渐变的方块叠加，适合明式案几和椅子
# TAGS: 锥形腿,收腰,渐变,案几,椅腿,桌腿,弧腿

import bpy

mat = bpy.data.materials.new("TaperedLeg")
mat.use_nodes = True
mat.node_tree.nodes["Principled BSDF"].inputs["Base Color"].default_value = (0.38, 0.19, 0.06, 1.0)
mat.node_tree.nodes["Principled BSDF"].inputs["Roughness"].default_value = 0.60


def make_tapered_leg(name: str, loc_x: float, loc_y: float,
                     top_w: float = 0.055, bot_w: float = 0.028,
                     height: float = 0.45, n_segs: int = 8) -> None:
    """创建收腰锥形腿，用 n 段立方体叠加模拟锥度。

    Args:
        name:    基础名称
        loc_x/y: 腿底部中心的 XY 坐标
        top_w:   腿顶部截面宽度（m）
        bot_w:   腿底部截面宽度（m）
        height:  腿高（m）
        n_segs:  分段数（越多越光滑，越慢）
    """
    seg_h = height / n_segs
    for k in range(n_segs):
        # 线性插值：底部最细，顶部最宽
        t = k / n_segs
        w = bot_w + (top_w - bot_w) * t
        z = k * seg_h + seg_h / 2
        bpy.ops.mesh.primitive_cube_add(size=1, location=(loc_x, loc_y, z))
        seg = bpy.context.active_object
        seg.name = f"{name}_seg{k}"
        seg.scale = (w / 2, w / 2, seg_h / 2 + 0.001)  # +1mm 避免缝隙
        if bpy.data.materials.get("TaperedLeg"):
            seg.data.materials.append(bpy.data.materials["TaperedLeg"])


# 示例：四腿配置（适合案几）
lx, ly = 0.44, 0.14
for i, (x, y) in enumerate([(lx, ly), (lx, -ly), (-lx, ly), (-lx, -ly)]):
    make_tapered_leg(f"leg_{i + 1}", x, y, top_w=0.05, bot_w=0.025, height=0.81)
