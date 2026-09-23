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
