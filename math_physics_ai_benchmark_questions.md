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


---

# 题三：黑洞 QNM 的谱不稳定性，究竟能否导致可观测振铃不稳定？

## 研究目标

对同一个明确给定的波动方程，同时研究复频共振极点和有限时间观测波形。定量判定“小扰动导致大的极点位移”与“小扰动导致大的波形变化”之间究竟有什么联系。

## 1. 模型与初值

采用 $G=c=1$，以 Schwarzschild 质量 $M$ 无量纲化：

$$
\tau=t/M,\qquad x=r_*/M,\qquad \rho=r/M,
\qquad x=\rho+2\log(\rho/2-1).
$$

考虑 $l=2$ odd-parity Regge–Wheeler master equation，加入一个外部小势垒：

$$
\left[\partial_\tau^2-\partial_x^2+V_{\varepsilon,L}(x)\right]
\psi(\tau,x)=0,
$$

$$
V_{\varepsilon,L}(x)=
\left(1-\frac2\rho\right)
\left(\frac6{\rho^2}-\frac6{\rho^3}\right)
+\varepsilon W(x-L),
\qquad 0<\varepsilon\ll1,
$$

$$
W(y)=
\begin{cases}
\exp\!\left(1-\dfrac1{1-y^2}\right),& |y|<1,\\
0,& |y|\ge1.
\end{cases}
$$

$
ho=\rho(x)>2$ 由 tortoise coordinate 关系确定。取 $L\ge L_0>12$，固定观察者 $x_o=10$，所以附加势垒位于观察者外侧。本题是指定 master equation 的扰动问题；不额外假设该势垒已经来自某个满足 Einstein equations 的物质模型。

初值固定为

$$
\psi(0,x)=C W(x),\qquad
\partial_\tau\psi(0,x)=-C W'(x),
$$

其中 $C>0$ 由未扰动背景上的初始能量

$$
E_0=\frac12\int_{\mathbb R}
\left(|\partial_\tau\psi|^2+|\partial_x\psi|^2+V_0|\psi|^2\right)dx=1
$$

确定，对所有 $\varepsilon,L$ 保持相同。定义观测信号

$$
h_{\varepsilon,L}(\tau)=\partial_\tau\psi_{\varepsilon,L}(\tau,x_o),
\qquad h_0=h_{0,L}.
$$

## 2. 共振与波形误差的定义

采用 $e^{-i\omega\tau}$ 时间约定。QNM 是频域 retarded Green function 解析延拓的极点，边界条件为

$$
\psi_\omega\sim e^{-i\omega x}\quad(x\to-\infty),
\qquad
\psi_\omega\sim e^{+i\omega x}\quad(x\to+\infty).
$$

阻尼态应有 $\operatorname{Im}\omega<0$。不要把这类 outgoing resonances 当作普通 $L^2(\mathbb R)$ 本征态。

定义绝对波形差的归一化观测量

$$
\mathcal E_{\varepsilon,L}(T)=
\left[
\frac{\int_0^T|h_{\varepsilon,L}-h_0|^2d\tau}
{\int_0^\infty|h_0|^2d\tau}
\right]^{1/2},
\tag{Q3a}
$$

以及在两信号范数非零时定义的白噪声形状 mismatch

$$
\mathcal M_{\varepsilon,L}(T)=1-
\frac{|\langle h_{\varepsilon,L},h_0\rangle_T|}
{\|h_{\varepsilon,L}\|_T\,\|h_0\|_T},
\qquad
\langle f,g\rangle_T=\int_0^T f(\tau)\overline{g(\tau)}d\tau.
\tag{Q3b}
$$

时间原点固定，不额外最大化 time shift。这是本题规定的数学度量，不等同于已经使用实际探测器 PSD 的检测结论。

## A. 先从因果传播计算波形

用 retarded Green function / Duhamel formula 导出 $h_{\varepsilon,L}-h_0$ 的首阶表达式，并说明余项控制范围。

根据初值支集、势垒支集和观察者位置，求出扰动最早能影响观察者的时间；在此之前证明两个波形精确相同。之后计算 (Q3a)、(Q3b) 的首个可能非零阶数，并区分 generic leading order、特殊抵消和归一化失效。

不能把 Gaussian 初值的“很小的尾巴”当作本题 compact support 初值；也不能把迟到回波相对于一个已经衰减至极小的背景信号的巨大比值，直接当作 (Q3a) 的巨大误差。

## B. 对同一系统研究极点迁移

对未扰动的一个 simple QNM $\omega_0$，从 Jost solutions、Wronskian 或 resolvent 出发求首阶共振位移，并给出真正控制微扰展开的参数。

研究

$$
L=c\log(1/\varepsilon),\qquad c>0
$$

的双重极限。判定何时可以出现不再随 $\varepsilon$ 消失的共振位移或新共振分支。若首阶展开已失效，不能继续线性外推，必须使用 resummation、精确共振方程或有误差控制的数值论证。

把此结果与同一 $(\varepsilon,L)$ 下的时域信号比较，并交代 pole residues、prompt response、Schwarzschild low-frequency branch cut / late-time tail 各自的作用。不给出适用时间窗，不得直接将一个有限 QNM sum 当作完整 retarded solution。

## C. 研究层：寻找统一的观测稳定性定理或反例

判定以下命题是否成立：存在 $\alpha>0$、$C_*<\infty$ 和 $\varepsilon_0>0$，使对给定初值和观察者

$$
\sup_{L\ge L_0}\sup_{T\ge0}
\mathcal E_{\varepsilon,L}(T)
\le C_*\varepsilon^\alpha,
\qquad 0<\varepsilon<\varepsilon_0.
\tag{Q3c}
$$

若成立，尽量求最佳幂次、可能需要的 logarithmic corrections，以及控制常数所依赖的条件；若不成立，构造一族 $\varepsilon_j\to0$、$L_j$、$T_j$ 和明确下界。固定 $T$ 的连续依赖估计不等于上述统一命题。

若统一命题尚无法闭合，至少在 $T=O(L)$ 且 $L=c\log(1/\varepsilon)$ 的窗口给出带误差控制的结果，并说明限制来自低频、共振寿命、远区 tail 还是估计工具。

最后把观测结果组织成一个带 source 和 detector 的 resolvent bound 或等价判据。若使用 pseudospectrum，必须给出函数空间、范数以及允许扰动类，并说明其与本题实际势垒扰动及单点观测之间的关系。目标是判断：大的谱位移能否与 (Q3c) 相容，而不是只画一张谱图。

建议数值部分使用两种不同的离散精度或算法核验，报告 outgoing boundary treatment、计算域、网格/阶数、误差及晚时域反射污染。未运行的代码照常可交，但必须标为 `NOT_RUN`。

## 最低可评价交付

至少完成一个严格因果结论和一个定量结果：共振位移的适用条件，或带明确时间窗的波形误差 bound。仅陈述“频域不稳定但时域稳定”不算完成。
