# 前沿数学物理 AI 研究能力题库 v1.0

编制日期：2026-09-17。

本题库包含三道独立研究型题目。每题都包含可核验的基础部分和研究层；研究层不预设候选结论成立。严格反例、限定适用范围的定理，以及可复现的非平凡部分解答，均属于有效成果。

“研究层”表示本题希望推进的具体问题，不表示已通过穷尽查新证明它是公开未解问题。题面所需的对象、允许操作及候选关系均在下文给定；不得以改变对象定义的方式制造解答。

## 给被测模型的统一说明

请直接开展研究，而不是只写研究计划或文献综述。优先完成能独立核验的低阶计算，再推进一般结论。

请分别标记：**已经证明**、**仅有精确低阶检验**、**仅有数值证据**、**猜想或尚未闭合的步骤**。引用外部结果时给出准确来源和适用条件；不得以相关文献的存在代替本题中的推导。

允许符号计算、数值计算和形式化证明。没有执行的代码必须标为 `NOT_RUN`；不得虚构运行日志、编译成功、数值输出或文献。不要求使用 Lean；使用 Lean 时必须说明哪些假设尚未消去。

若候选关系不成立，请给出最小反例或不可解性证书，并说明它排除了多大的候选类。若没有完成，请交付最强的已证命题、卡住的精确等式/引理，以及一个能使剩余猜想被证伪的下一步检验。


---

# 题一：树级 YM→EYM 共线重构能否提升到一圈？

## 研究目标

判断一个**同一套、与 helicity 无关、局域动量权重**的共线重构核，能否同时重构树级 EYM 振幅和指定有理 helicity 扇区的一圈 EYM 振幅。重点不是验证某个熟悉公式，而是区分“某个树级代表失败”与“整个允许的重构类不可能”。

## 1. 理论、正则化与振幅

考虑 minimal Einstein–Yang–Mills：

$$
S=\int d^D x\,\sqrt{-g}\left[
\frac{2}{\kappa^2}R-\frac14F^a_{\mu\nu}F^{a\mu\nu}
\right],\qquad
F^a_{\mu\nu}=\partial_\mu A^a_\nu-\partial_\nu A^a_\mu
+g f^{abc}A^b_\mu A^c_\nu .
$$

外态采用四维 spinor-helicity，圈动量采用 $D=4-2\epsilon$ 维正则化。采用 't Hooft–Veltman scheme，且两边使用一致的圈积分归一化与色生成元归一化。所有动量取 outgoing，允许一般复在壳运动学。

记 $M^{(\ell)}_{n;1}(1,\ldots,n;P)$ 为一个引力子、$n$ 个胶子的 leading-color single-trace 振幅，$A^{(\ell)}_{n+2}(\sigma)$ 为纯 YM 对应的 leading-color color-ordered 振幅。树级分别剥离 $\kappa g^{n-2}$ 与 $g^n$；一圈分别剥离 $\kappa g^n$ 与 $g^{n+2}$ 及一致的颜色/圈公共因子。

一圈只研究 $\kappa g^n$ 阶：内部物理态为胶子；协变规范计算应包含所需 ghost。**不包含**更高 $\kappa$ 阶引力圈，也不将额外 dilaton、antisymmetric tensor 或 matter states 混入目标理论。

## 2. 共线极限与允许的核

把引力子 $P^{+2}$ 替换为两个 $+1$ helicity 胶子 $a,b$，并固定

$$
\lambda_a=\sqrt{x}\lambda_P,\quad
\widetilde\lambda_a=\sqrt{x}\widetilde\lambda_P,\qquad
\lambda_b=\sqrt{1-x}\lambda_P,\quad
\widetilde\lambda_b=\sqrt{1-x}\widetilde\lambda_P,
\qquad 0<x<1.
$$

因此 $p_a=xP$、$p_b=(1-x)P$、$P^2=0$。负 helicity 引力子的树级校准使用两个负 helicity 胶子。对引力子取 $\varepsilon^{\pm2}_{\mu\nu}=\varepsilon^\pm_\mu\varepsilon^\pm_\nu$，其归一化应与以上剥离约定一致。

令 $\Pi_n$ 是保留硬胶子循环次序 $(1,\ldots,n)$、插入 $a,b$，并要求 $a,b$ 在循环意义下不相邻的全部排序，循环等价的排序只计一次。记

$$
\mathcal C_x A(\sigma)=
\left.A(\sigma)\right|_{p_a=xP,\,p_b=(1-x)P} .
$$

只在其余运动学一般、没有额外 factorization singularity 的区域定义该量。若使用其它表示，必须证明其中的表观奇异项已抵消。这里不允许把相邻两腿的 $1/\langle ab\rangle$ 项直接删去，冒充同一个定义。

允许的重构核为

$$
K_{n,\sigma}(x,s)=\sum_{i<j}c^{\sigma}_{ij}(x)\,s_{ij},
\qquad c^{\sigma}_{ij}(x)\in\mathbb Q(x),
\qquad s_{ij}=2p_i\!\cdot p_j,
$$

其中指标取遍 $1,\ldots,n,P$，并模掉在壳和动量守恒给出的关系。核在 Mandelstam 变量中是一次齐次多项式，不允许 Mandelstam 分母，不允许 helicity、polarization reference spinor 或圈动量依赖。

候选命题为：存在同一套 $K_{n,\sigma}$，使

$$
M^{(\ell)}_{n;1}
=\sum_{\sigma\in\Pi_n}
K_{n,\sigma}(x,s)\,\mathcal C_x A^{(\ell)}_{n+2}(\sigma).
\tag{Q1}
$$

在 $\ell=0$ 时要求对所有一般树级 helicity configurations 成立；在 $\ell=1$ 时先只要求在以下两个扇区成立，并取积分后的 $\epsilon^0$ 系数：

$$
(1^+,\ldots,n^+;P^{+2}),\qquad
(1^-,2^+,\ldots,n^+;P^{+2}),
$$

后一种情形也包括负 helicity 硬腿的所有位置。所有 $x$ 要求均指函数域恒等式，而不是仅在 $x=1/2$ 成立。

## A. 低点校准

从 $n=3$ 开始，明确写出排序集合、独立 Mandelstam 变量、至少一个非零树级校准，以及上面两个一圈 helicity 扇区的实际 EYM 结果。

请选择并明确写出一个标准树级重构代表，计算它在一圈的缺陷

$$
\Delta_{n}=M^{(1)}_{n;1}
-\sum_{\sigma}K^{\mathrm{tree}}_{n,\sigma}\,
\mathcal C_x A^{(1)}_{n+2}(\sigma).
$$

给出精确表达式及一个独立校验。只写“树级关系一般不能提升到一圈”不算完成本部分。

## B. 判定整个核空间，而不只判定一个代表

令 $\mathcal N_n$ 为所有允许核中在**全部树级 helicity configurations** 上给出零的核空间：

$$
\mathcal N_n=
\left\{N:\sum_\sigma N_{n,\sigma}\,
\mathcal C_x A^{(0)}_{n+2}(\sigma)=0\right\}.
$$

判定是否能通过 $K^{\mathrm{tree}}\mapsto K^{\mathrm{tree}}+N$、$N\in\mathcal N_n$ 消去两个一圈扇区的全部缺陷。优先完成 $n=3$，再处理 $n=4$。

若可解，给出核和剩余自由度；若不可解，给出与运动学、helicity 测试相对应的有限线性约束及不可解性证书，或同等强度的解析反证。随机点求解后必须在未参与求解的点验证，且不能将数值秩直接冒充函数域秩证明。

## C. 研究层：失败来自哪里，最小修正是什么？

用 $D$-dimensional unitarity 或等价的完整方法定位结果。须保留

$$
\ell^\mu=\bar\ell^\mu+\ell_\perp^\mu,
\qquad \mu^2=-\ell_\perp^2,
$$

并交代相关 $\mu^{2r}$ 项在积分后的作用。区分四维 cut、$D$ 维 cut、integrand identity 和 integrated identity。

若 (Q1) 失败，选择一种**事先限定自由度**的扩张，例如有限类 $\mu^{2r}$ 插入、明确的 dimension-shifted integral basis，或可由给定作用量独立计算的额外 building blocks；证明其至少能处理第一个障碍，并提出不依赖逐点拟合的因子化/归纳规则。不得把 $\Delta_n$ 本身重新命名为“修正项”作为完成。

若 (Q1) 成立，推进到 $n=5$ 或给出一般 $n$ 的证明结构，并明确接触项、无切割有理项及内部态闭合还剩哪些缺口。一个已经证明的低点 no-go 也是有效终点，不要求为了“给出关系”而改变原命题。

## 最低可评价交付

至少交付一个非平凡的精确低点结果，并说明它只约束某个核代表，还是已经约束整个允许核空间。
