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

# 题二：Fibonacci 融合约束下，Page 转折是否意味着量子信息可恢复？

## 研究目标

在一个完整给定的 random-code toy model 中，分别计算熵与信息恢复能力，判断 non-invertible fusion data 是否造成可观测的恢复阈值修正。该模型可用于讨论 black-hole information 与 holographic reconstruction，但**不预设它已经具有某个引力对偶**。

## 1. 融合空间与局域可观测量

采用 unitary Fibonacci fusion category：

$$
1\otimes a=a,\qquad \tau\otimes\tau=1\oplus\tau,
\qquad d_1=1,\quad d_\tau=\varphi=\frac{1+\sqrt5}{2}.
$$

在中间道基底 $(1,\tau)$ 中，取

$$
F^{\tau\tau\tau}_{\tau}=
\begin{pmatrix}
\varphi^{-1}&\varphi^{-1/2}\\
\varphi^{-1/2}&-\varphi^{-1}
\end{pmatrix}.
$$

本题不涉及 braiding。用与上述 unitary convention 一致的正交融合树基底。

定义

$$
\mathcal H_N^1=\operatorname{Hom}(1,\tau^{\otimes N}).
$$

将前 $N_R$ 个 anyons 记为辐射 $R$，其余 $N_B=N-N_R$ 个记为剩余系统 $B$，先取 $N_R,N_B\ge2$。定义 $\mathcal H_R^a=\operatorname{Hom}(a,\tau^{\otimes N_R})$、$\mathcal H_B^a=\operatorname{Hom}(a,\tau^{\otimes N_B})$，于是

$$
\mathcal H_N^1\cong
\bigoplus_{a\in\{1,\tau\}}
\mathcal H_R^a\otimes\mathcal H_B^a,
\qquad
m_a=\dim\mathcal H_R^a,\quad n_a=\dim\mathcal H_B^a,
\quad D_N=\sum_a m_an_a.
$$

$P_a$ 表示上述直和分解中的投影。只作用于辐射的可观测代数为

$$
\mathcal A_R=
\bigoplus_a\operatorname{End}(\mathcal H_R^a)\otimes I_{\mathcal H_B^a}.
$$

观察者可以测量 $a$ 并按扇区实施任意操作，但不能制造不同 $a$ 之间的相干。不得把非整数 $d_\tau$ 当作一个普通 Hilbert space 的维数。

## 2. 随机编码与允许的恢复

引入一个 topologically neutral 的 $k$ 维逻辑系统 $L$，$2\le k\le D_N$。令

$$
V:\mathbb C^k\longrightarrow\mathcal H_N^1,
\qquad V^\dagger V=I_k
$$

从 Haar/Stiefel measure 取样，即一个 $D_N\times D_N$ Haar unitary 的前 $k$ 列。

本题定义辐射与剩余系统的量子信道为

$$
\mathcal N_R(\rho)=
\bigoplus_a\operatorname{Tr}_{\mathcal H_B^a}
(P_aV\rho V^\dagger P_a),
$$

$$
\mathcal N_B(\rho)=
\bigoplus_a\operatorname{Tr}_{\mathcal H_R^a}
(P_aV\rho V^\dagger P_a).
$$

这些输出使用普通矩阵 trace 归一化。恢复器 $\mathcal D$ 是从直和输出代数到 $k$ 维 neutral system 的任意 CPTP map；它可以使用 neutral ancilla，但不额外提供带非平凡拓扑荷的资源。

取与 $L$ 最大纠缠的 neutral reference $Q$：

$$
|\Phi_k\rangle_{QL}=\frac1{\sqrt{k}}\sum_{j=1}^k|j\rangle_Q|j\rangle_L.
$$

定义最佳 entanglement recovery fidelity

$$
F_{\rm rec}(V)=\sup_{\mathcal D}
\langle\Phi_k|
\big[\operatorname{id}_Q\otimes(\mathcal D\circ\mathcal N_R)\big]
(|\Phi_k\rangle\langle\Phi_k|)
|\Phi_k\rangle.
\tag{Q2}
$$

## A. 熵的校准：先处理单个 Haar 随机纯态

暂时令 $k=1$，将随机态的系数按切分写成矩阵 $X_a\in\mathbb C^{m_a\times n_a}$，满足

$$
\sum_a\operatorname{Tr}(X_aX_a^\dagger)=1,
\qquad p_a=\operatorname{Tr}(X_aX_a^\dagger),
\qquad \rho_{R,a}=X_aX_a^\dagger/p_a.
$$

比较两种明确给定的熵：

$$
S_{\rm alg}=-\operatorname{Tr}(\rho_R^{\rm alg}\log\rho_R^{\rm alg}),
\qquad \rho_R^{\rm alg}=\bigoplus_aX_aX_a^\dagger,
$$

与

$$
\widetilde{\operatorname{Tr}}(Y)=\sum_a d_a\operatorname{Tr}(Y_a),
\qquad
\widetilde\rho_R=\bigoplus_a\frac{X_aX_a^\dagger}{d_a},
\qquad
S_{\rm qtr}=-\widetilde{\operatorname{Tr}}(
\widetilde\rho_R\log\widetilde\rho_R).
$$

从融合规则导出 $m_a,n_a$；推导 $p_a$ 的分布；给出 $\mathbb E S_{\rm alg}$ 与 $\mathbb E S_{\rm qtr}$ 的精确有限维公式，并推导 $N\to\infty$、$N_R-N_B=O(1)$ 时的首个非平凡修正。

至少计算一个二阶 replica moment，并核验上述 entropy conventions。必须区分 $\mathbb E[\log Z]$、$\log\mathbb E[Z]$ 和先对随机态归一化再平均。不改变物理切分的融合树基底变换不应改变可观测结论。

## B. 进入真正的信息恢复问题

恢复 $k\ge2$，定义

$$
\rho_{QB}=(\operatorname{id}_Q\otimes\mathcal N_B)
(|\Phi_k\rangle\langle\Phi_k|),\qquad
\rho_B=\operatorname{Tr}_Q\rho_{QB},
$$

$$
C_{QB}=\rho_{QB}-\frac{I_Q}{k}\otimes\rho_B,
\qquad \delta(V)=\frac12\|C_{QB}\|_1.
$$

从 Haar/Stiefel 平均出发，推导 $\mathbb E\operatorname{Tr}C_{QB}^2$ 的精确有限维表达式，或同等强度、带误差控制的结果。不能将 $V$ 的 $k$ 列当作互相独立且无需正交化的 Haar 随机态。

据此给出关于 $F_{\rm rec}$ 的非空洞 bounds，并处理 classical sector label 可能泄露到 $B$ 的信息。说明哪些维数条件足以保证典型的高保真恢复，哪些只能保证平均态看起来接近最大混合。

## C. 研究层：临界窗口与 island interpretation

先固定 $k=2$，再允许 $k$ 随 $N$ 增长。检验下列候选临界变量是否正确：

$$
\Xi=N_R-N_B-\log_{\varphi}k.
$$

在允许的整数子序列上，研究 $\Xi=O(1)$ 时 $F_{\rm rec}$ 的极限、上下界或极限分布。目标是以下任一种：可证明的 crossover formula；具有相同临界缩放的上下界；或一个推翻该变量/普适性的明确反例。

具体回答：两种熵的 Page 转折，能否单独决定 $F_{\rm rec}$ 的转折？$d_\tau$、扇区权重和融合重数分别在哪一步进入可操作的信息恢复？哪些影响只改变 entropy convention，哪些改变了本题已定义的量子信道？

最后提出与 quantum extremal surface / island prescription 比较所需的最小字典：什么对应 bulk entropy、什么可能对应 center/edge contribution，仍缺少哪些动力学和编码假设。不得仅凭熵曲线形状便宣称“证明了 island”或“构造了全息对偶”。

建议对至少两个小规模、包括不等切分的 $(N_R,N_B)$ 做独立验证。数值恢复可使用显式 decoder 或 CPTP 优化；只算纯态 entropy 不算完成恢复检验。

## 最低可评价交付

至少完成一个带扇区与正确归一化的精确平均量，并给出一条确实约束 (Q2) 的定理或定量 bound，而不只是画出 Page 曲线。
