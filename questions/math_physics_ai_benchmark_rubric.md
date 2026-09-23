# 出题者评分说明：前沿数学物理 AI 研究能力题库 v1.0

编制日期：2026-09-17。

本文件保留给评阅者，不与题面一起发给被测模型。它不是完整答案库；研究层的正确性需要独立审计。下述分数衡量本题中的可核验表现，不换算成“博士水平”等未经校准的人类能力标签。

## 一、测试协议

每题单独开启新对话，使用相同 system prompt、输出预算、工具权限及重试次数。闭卷/无工具与开放文献/可运行工具分别报告，不合并排名。记录模型版本、日期、实际工具调用、实际输出预算、重试和人工介入。

第一轮不提供提示。第二轮允许向所有模型提供同类型、同信息量的一个检查要求，例如“请核验你最关键的一条非平凡等式，并修正结论”。分别报告首轮表现与修正表现，不能只保留最好的一次。

对引用结果抽查原文；对核心符号等式使用另一条推导或独立 CAS；对数值要求未参与拟合的 hold-out checks。研究型答案应报告一个可信分数区间，而不是在证据不足时给出伪精确分数。

## 二、每题 100 分的公共结构

| 维度 | 分数 | 计分依据 |
|---|---:|---|
| 定义与题设审计 | 15 | 对象、约定、操作权限与适用域正确；发现真实缺口并最小修正，而非任意改变题目。 |
| 已知/低阶校准 | 30 | 至少一项非平凡精确结果；步骤可复查，归一化、边界条件或扇区正确。 |
| 独立验证 | 15 | 第二种方法、符号恒等式、hold-out 数据或有实际日志的代码；没有工具时可用完整解析交叉检验。 |
| 研究推进 | 30 | 新的严格归约、非空洞 bound、不可解性证书或可验证一般结构。正确反例与正面证明同等计分。 |
| 结论校准与可继续性 | 10 | 明确区分证明/计算/猜想，交付精确缺口和可证伪下一步，不把未知步骤藏在“显然”中。 |

达到高分不以完成整个研究层为必要条件。一个严格限定范围且验证充分的 no-go，可能优于覆盖所有小问但关键步骤错误的长答案。

不要按篇幅、术语数量或引用数量评分；不因模型说“该问题很难”而加分；仅有研究计划、文献清单或未经执行的长代码，通常不能替代低阶计算。

## 三、题一的核验点

### 可客观核验的基础

检查外态、内部态、颜色阶数、耦合阶数、loop measure 及维正则方案。检查核确实只在 Mandelstam 中一次齐次，而不是偷偷使用任意有理函数插值。

优先核验 n=3 的排序、非零树级校准以及两个一圈 helicity 扇区。可用参考文献 R1–R3 对照定义与低点计算。参考文献中的单个树级代表不能代替本题允许核空间的分类。

### 最能拉开差距的部分

真正计算树级零空间，再检查一圈缺陷能否落入它的一圈像空间；构造有限不可解性证书；严格区分函数域恒等式与固定 x、固定运动学测试。

正确处理四维 cut 不能独自决定的部分，说明 μ² 及积分次序，并保持 pure EYM 的内部态内容。若扩张 building blocks，应独立定义并从作用量/切割生成，而非直接把答案差值加入。

### 明确扣分点

“一个熟悉公式失败 ⇒ 一切 YM→EYM 关系均不存在”；“四维 cuts 为零 ⇒ 圈振幅为零”；“树级 BCJ 关系自动作为积分后一圈关系使用”；把单点数字拟合作为精确证明；用额外物质态改变目标理论但不说明。

## 四、题二的核验点

### 可客观核验的基础

融合重数必须是整数，符合初值与递归。核验 D_N=Σ_a m_a n_a。核验随机态每个样本先归一化、p_a 的联合分布及条件块的 Wishart 结构。q-trace 与普通 trace 所定义的两个密度矩阵不可混用。

单纯 anyonic Page curve 的平均与方差已有公开研究，见 R4–R6；模型重新得到这些结果属于校准，不自动算研究新颖性。

### 最能拉开差距的部分

针对 Haar isometry 而非独立随机列推导含 k 的精确二阶平均；从该平均到 trace norm / recovery fidelity 给出合法且非空洞的推论；处理 sector label 泄漏而不只检查扇区内部去耦。

在临界窗口给出真正约束最优 decoder 的结果；区分本题允许操作与改变拓扑资源后的另一种任务；给出可执行 decoder 或 CPTP 优化的有限维检查。

### 明确扣分点

把 φ 当作 Hilbert space 维数；把直和硬改成普通 tensor product；把 entropy 接近最大当成任意码信息可恢复的充分证明；把 ensemble average state 的性质当成典型样本的性质；忽略 Stiefel 正交性；依据 Page 曲线形状宣称构造了引力对偶。

## 五、题三的核验点

### 可客观核验的基础

检查无量纲势、阻尼号、初始能量与统一初值。核验最早影响时刻必须由初值→势垒→观察者的传播链得出，不能只使用势垒→观察者距离。

检查 Q3a 与 Q3b 的区别：绝对差、形状 mismatch、固定时间窗与晚期相对比值不是同一量。固定 T 的 Duhamel 估计不自动给出对 L,T 统一的 bound。

### 最能拉开差距的部分

复共振微扰能给出有效展开参数，并在双重极限中识别何时不能外推；对同一系统处理 pole、residue、branch cut 与因果波形；对 Q3c 给出真正统一的证明/反例，或准确识别无法闭合的频域区间。

数值验证至少做一次独立收敛或方法对照，并排除有限计算域的人工回波。若使用 pseudospectrum，需和实际 perturbation class 及 source-to-detector response 联系起来，而不是以矩阵条件数取代物理量。

### 明确扣分点

把 QNM 当 L² 正交本征基；在回波因果时间之前预言改动波形；只比较频率不比较 residue；用未经控制的有限 QNM 和忽略 Schwarzschild tail；把相对一个接近零的背景的大比值当成 Q3a 的反例；对微扰公式作失效区外推。

## 六、分数如何解读

| 得分区间 | 本题中的证据强度 |
|---|---|
| 0–24 | 尚未建立可信的工作对象，或主要停留在相关术语和方案。 |
| 25–44 | 有正确框架和部分计算，但关键校准/验证不足。 |
| 45–64 | 已形成可复核的研究型部分解答，能识别关键限制。 |
| 65–79 | 低阶计算扎实，并完成非平凡归约、bound 或严格反例。 |
| 80–100 | 多条独立验证支持核心结果，研究层有明确实质推进；不要求包办整个开放问题。 |

这个表是任务内解释，不是跨学科、跨工具、跨预算的能力排名。三题应分别给分，并同时报告“最强已证结论”“主要错误”“下一步是否值得投入”。

存在虚构计算/引用时，相关验证分应为零，并单独标记可靠性风险。发现一个错误不必机械地将所有正确独立部分清零；但依赖该错误的结论不能继续计分。能够主动发现并正确修复错误，应体现在修正轮成绩，而不是抹去首轮记录。

## 七、文献与已知边界（截至编制日检索）

以下是核验依据与背景，不是对研究层“无人做过”的担保。

**R1.** Stephan Stieberger, Tomasz R. Taylor, *Subleading Terms in the Collinear Limit of Yang-Mills Amplitudes*, arXiv:1508.01116 (2015)。用于树级共线定义、一个引力子与两个胶子的关系及运动学约定。

**R2.** Dhritiman Nandan, Jan Plefka, Gabriele Travaglini, *All rational one-loop Einstein-Yang-Mills amplitudes at four points*, arXiv:1803.08497 (2018)。低点校准见 Sections 3–5；一种直接共线提升的检验见 Section 11。该文不等于本题整个局域核空间的不可解性定理。

**R3.** Franziska Porkert, Oliver Schlotterer, *One-loop amplitudes in Einstein-Yang-Mills from forward limits*, arXiv:2201.12072 (2022)。一般一圈 integrand 构造与 pure gauge building blocks 的关系已有结果；题一不能把这一整体方向当作新发现。

**R4.** Yale Yauk, Lucas Hackl, Alexander Hahn, *Typical entanglement in anyon chains: Page curves beyond Lie group symmetries*, arXiv:2603.25789 (2026)。anyonic 随机态熵统计已有解析工作；题二研究层是题面固定的随机编码/恢复任务，而非声称首次研究 anyonic Page curve。

**R5.** Eugenio Bianchi, Pietro Donà, Rishabh Kumar, *Non-abelian symmetry-resolved entanglement entropy*, arXiv:2405.00597 (2024)。操作代数与非阿贝尔扇区的重要背景。

**R6.** Eugenio Bianchi, Pietro Donà, *Typical entanglement entropy in the presence of a center: Page curve and its variance*, arXiv:1904.08370，Physical Review D 100, 105010 (2019)。具有 center 的随机态统计背景。

**R7.** Patrick Hayden, John Preskill, *Black holes as mirrors: quantum information in random subsystems*, arXiv:0708.4025 (2007)。随机编码、黑洞信息恢复的背景，不直接代替本题 constrained channel 的计算。

**R8.** Ramin G. Daghigh, Guan-Ru Li, Wei-Liang Qian, Stefan J. Randow, *Evolution of black hole echo modes and the causality dilemma*, arXiv:2502.05354 (2025)。谱迁移、回波与因果传播已有研究；题三的统一不等式需另行证明或证伪。

**R9.** Li-Ming Cao, Ming-Fei Ji, Liang-Bi Wu, Yu-Sen Zhou, *Pseudospectrum and time-domain analysis of the EFT corrected black holes*, arXiv:2508.13894 (2025)。谱性质与波形 mismatch 的区分已有具体研究；不同模型上的缩放不能直接套到 Q3c。

**R10.** Jin Dong, Stephan Stieberger, *Subleading Collinear Limits of Yang-Mills Amplitudes from Gravity*, arXiv:2609.00153 (2026)。近期树级共线研究背景；不将摘要中的一般表述当作已核验的全部定理，更不把它直接提升为一圈结论。

## 八、简短评阅记录模板

- 题号、模型、日期、工具模式、预算、首轮/修正轮：
- 最强已证结果：
- 最重要的独立验证：
- 致命错误或尚未闭合的步骤：
- 是否有仅适用于某个低点/扇区/范数的结论被错误推广：
- 五维分数（定义 / 校准 / 验证 / 推进 / 结论校准）：
- 值得继续的一个具体工作项及其证伪检验：
