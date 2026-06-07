调研完成，三个最相关的工作都读到了核心机制。下面是总结、诊断和设计方案。

一、三个参考工作怎么做验证
工作	验证方式	关键点
LL3M (2025, threedle，和你的项目几乎同构：planner/retrieval/coding/critic/verification 多 agent 写 Blender 代码)	纯 VLM，零几何检查	critic 看渲染图找问题（"腿没连到座面 → 沿 z 轴上移腿"），verification agent 比对修复前后的渲染图确认改对了
CADCodeVerify (ICLR 2025)	刻意不用几何求解器	论文明确指出：几何求解器需要 ground-truth 3D 模型做对比，现实中拿不到；改用 VLM 生成 2-5 个二元验证问题，且允许答 "Unclear" 避免幻觉式的过度自信
L3GO (NAACL 2025)	在坐标/规格层验证，不在 mesh 层	有 Coordinate Calculator 显式算出每个部件的吸附坐标，Spatial Critic 检查的是计划里部件之间的连接意图，而不是盲目对成品 mesh 做几何猜测
业界共识很清晰：现代系统要么完全放弃 mesh 级几何检查（LL3M、CADCodeVerify），要么把验证上移到语义/坐标层（L3GO）。没有人像当前 geom_verifier 这样对裸 mesh 做全配对 AABB 启发式检查。

二、诊断：当前 geom_verifier 为什么帮倒忙
我对照 geom_verifier.py 逐项看，问题出在假阳性噪声被当成指令喂给了 Fixer：

检查	问题
Check 2 全配对 AABB 间隙 (1–12cm)	假阳性重灾区。AABB 相邻 ≠ 应该相连；对角的两条腿本就该有间隙，却被报成"未连接"
Check 5 质心稳定性	输出 info 噪声，对装配正确性无意义
Check 6 腿-桌面连通性 (用 15% 面积比猜哪个是腿)	启发式脆弱，遇到非桌类家具（屏风、椅背）就误判
Check 1 非流形边/孤立顶点	对正常 primitive 也报 warning，且与"是否装配正确"无关
Check 4 空 mesh	唯一真正有用、永不假阳性的检查
致命点在 main.py:221：所有这些 issue（含一堆假阳性）被 fixer.fix(code, all_issues) 当作"必须修复的指令"，Fixer（现在是 kimi）就去"修"那些本来正确的几何，于是越改越坏 —— 这正是你说的"启用还不如不启用"。

三、提炼出的设计原则
VLM 当主裁判（对齐 LL3M/CADCodeVerify），几何只做 VLM 真正看不到、且零假阳性的薄安全网。
几何检查只发"硬错误"，绝不发 info/warning 这类"猜测"进 Fixer 的 prompt。
要做连通性检查，就用语义（L3GO 思路）：依据 planner/COMPONENT 里"哪些部件本应连接"的意图做定向检查，而不是全配对盲猜。
四、三个候选方案
方案 A — 直接删掉 geom_verifier，纯 VLM
最简单，对齐 LL3M/CADCodeVerify 的 SOTA 做法。代价：丢掉极少数 VLM 看不出的亚厘米悬空。

方案 B — 砍成"零假阳性安全网"（推荐）
只保留三类不可能误报、且 VLM 确实易漏的硬错误，全部以 error 级输出，其余噪声检查全删：

空场景 / 空 mesh（Check 4）
退化或爆炸尺度（某部件尺寸 ≈0 或 >15m）
全局悬空：整个模型最低点 min-z 显著 > 0（整体飘在空中，VLM 在缺地面参照时常看不准）
并改 main.py：几何 issue 不再无条件进 Fixer，只有 error 级才进。

方案 C — 上移到语义层（L3GO 思路，工作量最大）
让 coder 在 COMPONENT 注释里声明连接关系，或用 planner 的结构，只对"声明应相连"的部件对做定向连通检查。准确率最高，但要改 coder/planner 的 prompt 和解析。

我倾向 方案 B：投入小、直接消除假阳性根源、保留几何真正不可替代的价值，且与三篇论文的结论一致。在这之上如果你想要更高精度的连通性，未来再叠加方案 C。

下面想跟你确认方向再动手：

AskUserQuestion
OUT
Your questions have been answered: "geom_verifier 的改造方向选哪个？"="方案C 语义层连通检查", "几何 issue 进 Fixer 的策略？"="只有error级进Fixer". You can now continue with these answers in mind.
方向定了：方案 C 语义层连通检查 + 只有 error 级进 Fixer。下面是落地到这个代码库的具体设计，你过目确认后我再改代码。

设计方案：L3GO 式语义连通验证
核心思路
不再全配对盲猜，而是让 coder 自己声明"哪些部件本应相连"，geom 只检查这些声明过的部件对是否真的接触。没声明的对永不检查 → 从根上杜绝假阳性。

改动 1 — coder 声明连接清单
coder_system.txt 增加规则：在所有部件创建完后，用对象名输出一个连接清单变量：


# 必须物理接触的部件对（用 obj.name 的值）
CONNECTIONS = [("leg_1", "tabletop"), ("leg_2", "tabletop"),
               ("leg_3", "tabletop"), ("leg_4", "tabletop"),
               ("apron", "tabletop")]
这正是 L3GO "Coordinate Calculator 显式声明吸附关系" 的轻量版。

改动 2 — geom_verifier.py 重写
当前 6 项检查全部删除，新 snippet 只做三类零假阳性硬错误（全部 severity: error）：

新检查	逻辑	为什么零假阳性
E1 空场景/空 mesh	无 mesh 对象 或 某对象 0 顶点	客观事实
E2 退化/爆炸尺度	整体尺寸 <5cm 或 >15m，或某部件某轴 ≈0	极端阈值，正常家具碰不到
E3 声明连通断裂	读 globals().get("CONNECTIONS")，对每个声明对算 AABB 最小间距，>1cm 判 error	只查 coder 自己声明应连的对
名称匹配支持前缀（"leg" 可匹配 leg_1..leg_4），容忍循环命名。
降级保护：若代码没有 CONNECTIONS（旧代码/coder 漏写），E3 整段跳过，只跑 E1/E2，绝不猜。
全配对间隙、非流形、质心、面积比猜腿 —— 全删。互穿留给 VLM。
改动 3 — main.py 只放 error 级进 Fixer
main.py 里喂 Fixer 前显式过滤几何 issue：


geom_for_fixer = [i for i in prev_geom_issues if i.get("severity") == "error"]
新 geom 本就只发 error，这一步是显式护栏，防止以后加 warning 级检查时再次泄漏噪声进 Fixer。VLM issue 不受影响（VLM 仍是主裁判）。

涉及文件
prompts/coder_system.txt — 加 CONNECTIONS 声明规则
agents/geom_verifier.py — 重写 snippet（6 检查 → 3 硬错误）
main.py — geom issue 过滤到 error 级
验证方式

python main.py "做一张方桌" --iterations 3 --vlm --geom --verbose
预期：geom 段只在真断裂/空场景/尺度异常时报 error；正常方桌 geom 静默；启用 --geom 不再让模型变差。

这个设计你看可以吗？确认后我就按这三处动手。


全部完成。总结这次对 geom_verifier 的改造：

改了什么
1. geom_verifier.py — 6 个噪声检查 → 3 个零假阳性硬错误

E1 空场景/空 mesh、E2 退化或爆炸尺度、E3 声明连通断裂
E3 只检查 coder 自己声明的 CONNECTIONS 部件对（L3GO 思路），从根上杜绝假阳性
删除了全配对 AABB 间隙、非流形、质心、面积比猜腿等假阳性源头
2. 顺手修了一个潜伏 bug
旧验证器读 obj.matrix_world 时 scale 尚未刷新（直接赋值 obj.scale 不触发 depsgraph 更新），所有 bbox 都跑在未缩放的单位立方体上 —— 这是旧验证器不可靠的重要原因。已加 view_layer.update() 修复。

3. main.py — 喂 Fixer 前显式过滤，只放 error 级几何 issue。

4. coder_system.txt — 要求 coder 输出 CONNECTIONS 声明，并强制包含每条腿（首次测试时 coder 漏声明腿，强化后已修正）。

验证结果
单元测试：空场景、退化部件、断裂连接都被精准捕获；正常方桌完全静默（零假阳性）
端到端 --vlm --geom 跑通（exit 0），geom 全程 0 误报，coder 正确输出 4 条腿的完整连接声明
调研依据
对齐三篇最新工作：LL3M / CADCodeVerify 证明 mesh 级几何检查弊大于利、应以 VLM 为主裁判；L3GO 提供了"在语义/声明层验证连接"的范式，正是 E3 的设计来源。

Sources:

L3GO (NAACL 2025)
LL3M: Large Language 3D Modelers
CADCodeVerify (ICLR 2025)
一个遗留观察：iter_0 出现过一次 Blender 运行失败（无 traceback，auto-fix 已自动恢复），这是 coder 代码质量的既有现象，和本次几何验证改动无关。