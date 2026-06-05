# TITLE: 明式家具标准比例常数
# DESCRIPTION: 明式各类家具的典型尺寸参数字典，可直接在代码中 import 使用
# TAGS: 比例,尺寸,参数,明式,标准,参考数据,所有家具

# 明式家具典型比例（单位：米）
# 参考《明式家具研究》（王世襄）及现代仿制品实测值

MING_DIMENSIONS = {
    # ── 椅凳类 ──────────────────────────────────────────────────────
    "官帽椅": {
        "width": 0.60, "depth": 0.52, "height": 1.20,
        "seat_height": 0.50, "seat_thickness": 0.05,
        "leg_section": 0.04,   # 腿截面尺寸
        "back_height": 0.70,   # 靠背高（座面以上）
        "armrest_height": 0.22,  # 扶手高（座面以上）
    },
    "太师椅": {
        "width": 0.68, "depth": 0.60, "height": 1.15,
        "seat_height": 0.52, "seat_thickness": 0.06,
        "leg_section": 0.06,
        "back_height": 0.63,
        "armrest_height": 0.24,
    },
    "圈椅": {
        "width": 0.64, "depth": 0.56, "height": 1.00,
        "seat_height": 0.48, "seat_thickness": 0.05,
        "leg_section": 0.038,
        "back_height": 0.52,  # 含椅圈
        "rail_radius": 0.30,  # 椅圈水平半径
    },
    "长凳": {
        "width": 1.20, "depth": 0.30, "height": 0.45,
        "seat_thickness": 0.05, "leg_section": 0.05,
    },
    # ── 桌案类 ──────────────────────────────────────────────────────
    "方桌": {
        "width": 0.90, "depth": 0.90, "height": 0.78,
        "top_thickness": 0.05, "leg_section": 0.06,
        "apron_height": 0.08,
    },
    "八仙桌": {
        "width": 0.95, "depth": 0.95, "height": 0.82,
        "top_thickness": 0.06, "leg_section": 0.065,
        "apron_height": 0.09,
    },
    "案几": {
        "width": 1.00, "depth": 0.35, "height": 0.85,
        "top_thickness": 0.04, "leg_section": 0.04,
        "apron_height": 0.07,
    },
    "书案": {
        "width": 1.50, "depth": 0.60, "height": 0.80,
        "top_thickness": 0.05, "leg_section": 0.05,
        "apron_height": 0.08,
    },
    # ── 床榻类 ──────────────────────────────────────────────────────
    "罗汉床": {
        "width": 2.00, "depth": 0.90, "height": 0.50,
        "bed_frame_h": 0.18, "railing_h": 0.55,
        "leg_section": 0.07,
    },
    # ── 杂项 ────────────────────────────────────────────────────────
    "屏风_单扇": {
        "width": 0.50, "depth": 0.05, "height": 1.80,
        "frame_thickness": 0.05, "panel_inset": 0.85,  # 内板占宽度比例
    },
}

# 用法示例：
# from examples.ming_proportions import MING_DIMENSIONS
# dims = MING_DIMENSIONS["官帽椅"]
# seat_h = dims["seat_height"]   # → 0.50
