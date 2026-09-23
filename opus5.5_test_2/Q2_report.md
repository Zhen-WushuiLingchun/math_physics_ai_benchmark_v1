# 题二研究报告：Fibonacci 融合约束下，Page 转折是否意味着量子信息可恢复？

所有代码在 `q2_work/` 中，**均已实际运行**，文中数值来自这些运行，见 §8 的清单。
状态标记：
**[证明]** 已经证明；
**[证明*]** 已证明，但依赖所引标准定理（出处和适用条件已注明）；
**[精确低阶]** 仅有精确低阶检验；
**[数值]** 仅有数值证据；
**[猜想]** 猜想或尚未闭合的步骤。

记号：$\varphi=d_\tau$，$\Delta=N_R-N_B$，$\Xi=\Delta-\log_\varphi k$，$c:=\varphi^{-\Xi}=k\varphi^{-\Delta}$，$\mathcal D^2=\sum_a d_a^2=1+\varphi^2=\sqrt5\,\varphi$，$d_R=\varphi^{N_R}$，$d_B=\varphi^{N_B}$。$F_n$ 为 Fibonacci 数（$F_0=0,F_1=1$），$H_n$ 为调和数，熵一律取自然对数。

---

## 0. 结论摘要

1. **维数与 F-符号 [证明]。** $m_1=F_{N_R-1}$，$m_\tau=F_{N_R}$，$n_1=F_{N_B-1}$，$n_\tau=F_{N_B}$，$D_N=F_{N-1}$。F-矩阵只用来具体实现 $\mathcal A_R$，已用数值精确核验。任何 Haar 平均量以及任一给定 $V$ 的可观测结论都不依赖 F-符号：不改变切分的基底变换恰好是 $\bigoplus_a W^R_a\otimes W^B_a$。

2. **A 部分 [证明]。**
   - $(p_1,p_\tau)\sim\mathrm{Dirichlet}(M_1,M_\tau)$，其中 $M_a=m_an_a$。
   - 逐点恒等式 $S_{\rm qtr}=S_{\rm alg}+p_\tau\ln\varphi$。
   - 两种熵均值都有精确有限维公式。
   - 在 $N\to\infty$、$\Delta$ 固定时：
   $$\mathbb E S_{\rm qtr}=\min(N_R,N_B)\ln\varphi-\tfrac12\varphi^{-|\Delta|}+\tfrac{\sqrt5\varphi}{2}\big(1+\tfrac13\varphi^{-|\Delta|}\big)\varphi^{-N}+O(\varphi^{-2N}),$$
   $$\mathbb E S_{\rm alg}=\mathbb E S_{\rm qtr}-\tfrac{\varphi}{\sqrt5}\ln\varphi+\tfrac{\varphi\ln\varphi}{\sqrt5}\big[(-1)^{N_B}\varphi^{\Delta}+(-1)^{N_R}\varphi^{-\Delta}\big]\varphi^{-N}+O(\varphi^{-2N}).$$
   - 两种约定的 Rényi-$n$ replica moment 在平面主阶上精确地等于"以量子维数 $d_R,d_B$ 代替 Hilbert 维数"的 Page moment：qtr 约定无前因子，alg 约定带常数因子 $\kappa_n=\sum_ap_ad_a^{n-1}$。

3. **B 部分 [证明]。** 精确有限维平均为
$$\boxed{\ \mathbb E\operatorname{Tr}C_{QB}^2=\Big(1-\frac1{k^2}\Big)\sum_a\frac{m_an_a\,(n_a-m_a/D)}{D^2-1}\ }$$
   并且各扇区分别成立。对每个 $V$ 都有 $(1-\delta)^2\le F_{\rm rec}\le1-\delta^2/4$。由精确二阶矩还得到两个有限 $N$ 下界：
   - 解耦下界：$\mathbb E F_{\rm rec}\ge(1-\bar\delta)^2$，其中 $\bar\delta\to\frac12\sqrt{1-k^{-2}}\,\varphi^{-\Xi/2}$；
   - 碰撞（Rényi-2）下界：$\mathbb E F_{\rm rec}\ge\sum_a\frac{(M_a/D)^3}{kn_a\mathbb E\operatorname{Tr}(\rho^{(a)}_{QB})^2}\to\dfrac{1}{1+\varphi^{-\Xi}}$。

   扇区标签的泄露有精确二阶矩：$\mathbb E\|Q_a-\tfrac{\operatorname{Tr}Q_a}{k}I\|_2^2=(k^2-1)M_a(D-M_a)/(D(D^2-1))$，随 $N$ 指数小。

4. **C 部分的结构定理 [证明]。**
   - 最优解码自动按扇区进行，$F_{\rm rec}=\sum_aq_aF_a=e^{-H_{\min}(Q|R)}/k$。
   - **扇区约化**：Fibonacci 码恰为"每个扇区一个普通 Haar 码"，前面接一个 $k\times k$ 滤波器，且二者独立。由此
   $$\Big|\mathbb EF^{\rm Fib}_{\rm rec}-\sum_ap_a\,\bar F^{\rm Haar}_k(m_a,n_a)\Big|\le\sum_a\sqrt{\tfrac{(k^2-1)M_a(D-M_a)}{kD(D^2-1)}}=O(\sqrt{k/D}).$$
   - 两个扇区的有效比值 $c_a=kn_a/m_a$ 都趋于同一个 $\varphi^{-\Xi}$，所以**$\Xi$ 是正确的临界变量，且没有依赖 $d_a$ 的常数偏移**，扇区权重也从极限中消失。
   - 对固定 $k$，**只依赖 $\Xi$ 的上下界 [证明*]**：
   $$\max\Big\{\tfrac1{k^2},\ \tfrac1{1+\varphi^{-\Xi}},\ \big(\mathbb E_{\rm MP(\varphi^{-\Xi})}\sqrt\lambda\big)^2\Big\}\ \le\ \liminf\mathbb EF_{\rm rec}\ \le\ \limsup\mathbb EF_{\rm rec}\ \le\ \min\Big\{1,\ \big(\varphi^{\Xi/2}+\tfrac1k\big)^2\Big\}.$$
   - $F_{\rm rec}$ 以 $\exp(-ckDt^2)$ 的速度集中，极限分布是点质量。

5. **数值 [数值]。** $k=2$ 沿 $\Delta$ 固定的整数子序列，$F_{\rm rec}$ 收敛到仅依赖 $\Xi$ 的值（§4.4 表），在 $c$ 精确匹配的点上与普通 Haar 码的差 $\le2.5\times10^{-4}$，与统计误差相当。另有两条渐近律：
$$1-F_{\rm rec}\simeq\tfrac14(1-k^{-2})\,\varphi^{-\Xi}\ (\Xi\gg1),\qquad F_{\rm rec}-\tfrac1{k^2}\simeq\tfrac{s_k}{k}\,\varphi^{\Xi/2}\ (\Xi\ll0),$$
   其中 $s_2\approx1.00$。第一条有启发式推导，第二条中 $s_k$ 被化为一个 GUE-SDP 常数。

6. **$k\to\infty$ [猜想，有强数值与证明概要]。**
$$F_{\rm rec}\to f_\infty(\varphi^{-\Xi}),\qquad f_\infty(c)=\Big(\mathbb E_{{\rm MP}(c)}\sqrt\lambda\Big)^2=\Big[\tfrac{\beta}{3\pi c}\big((\alpha^2+\beta^2)E(\kappa^2)-2\alpha^2K(\kappa^2)\big)\Big]^2,$$
   其中 $\alpha=|1-\sqrt c|$，$\beta=1+\sqrt c$，$\kappa^2=1-\alpha^2/\beta^2$。特别地 $F_{\rm rec}(\Xi=0)\to 64/(9\pi^2)\approx0.7205$。在 $\Xi\gtrsim-1$ 处偏差数值上约为 $O(k^{-2})$。

7. **对题目追问的回答。**
   - 两种熵的 Page 转折**不能**单独决定恢复转折。
   - $d_\tau$ 只通过融合重数的 Perron–Frobenius 增长进入（即 $\log_\varphi$）。
   - $\widetilde{\operatorname{Tr}}$ 中的 $\log d_a$ 是一个中心项，只改变熵约定，在一切条件熵以及 $H_{\min}(Q|R)$ 中精确抵消，因而不改变信道。
   - 扇区权重只以 $F=\sum q_aF_a$ 的混合权重进入，并在极限中消失。
   - 融合重数本身定义了信道。
   - 见 §5–§6。

---

## 1. 融合空间、可观测代数与基底

**1.1 维数 [证明]。** 设 $f_a(n)=\dim\operatorname{Hom}(a,\tau^{\otimes n})$。由 $\tau\otimes\tau=1\oplus\tau$ 得递推 $f_1(n+1)=f_\tau(n)$，$f_\tau(n+1)=f_1(n)+f_\tau(n)$，初值 $f_1(1)=0$，$f_\tau(1)=1$。于是 $f_\tau(n)=F_n$，$f_1(n)=F_{n-1}$，从而
$$m_1=F_{N_R-1},\quad m_\tau=F_{N_R},\quad n_1=F_{N_B-1},\quad n_\tau=F_{N_B},\quad D_N=F_{N_R-1}F_{N_B-1}+F_{N_R}F_{N_B}=F_{N-1},$$
最后一步用了恒等式 $F_{a+b}=F_aF_{b+1}+F_{a-1}F_b$。另有 $\sum_am_ad_a=\varphi^{N_R}$（即 $\tau^{\otimes N_R}$ 的量子维数），以及 $m_a=d_a\varphi^{N_R-1}/\sqrt5\,(1+O(\varphi^{-2N_R}))=\frac{d_ad_R}{\mathcal D^2}(1+\dots)$，$D=\frac{d_Rd_B}{\mathcal D^2}(1+\dots)$。

**1.2 $\mathcal A_R$ 由 R-局域算符生成 [精确低阶]**（`fib_basis.py`）。在标准左到右融合树基底中，用给定的 $F^{\tau\tau\tau}_\tau$（已核验 $F=F^T=F^{-1}$）构造相邻两 anyon 的融合道投影 $\Pi^{(i)}_1$。对 8 个切分 $(N_R,N_B)\in\{(2,2),(2,3),(3,2),(3,3),(4,2),(3,4),(4,4),(5,3)\}$ 核验了四件事：
- $\{\Pi^{(i)}_1\}_{i<N_R}\cup\{P_1\}$ 生成的 $*$-代数维数恰为 $\sum_am_a^2$；
- B 侧同理，维数为 $\sum_an_a^2$；
- R-局域算符与 B-局域算符对易到 $10^{-16}$；
- 跨越切分的 $\Pi^{(N_R)}_1$ 混合扇区（矩阵元 0.486）。

这说明题设的 $\mathcal A_R=\bigoplus_a\operatorname{End}(\mathcal H_R^a)\otimes I$ 正是 R 上局域可观测量生成的代数，其中心为 $\mathrm{span}\{P_a\}=Z(\mathcal A_B)$。

**1.3 基底无关性 [证明]。** 任何不改变物理切分的正交融合树基底变换，都由 R 内部的 F-move、B 内部的 F-move，以及把 R 的总荷 $a$ 作为旁观者的重耦合组成，因此都保持中心投影 $P_a$ 不变（$P_a$ 是 R 的总荷，与树的选择无关），形如 $W=\bigoplus_aW_a^R\otimes W_a^B$。由此：
- Haar 测度在 $V\mapsto WV$ 下不变；
- $S_{\rm alg}$、$S_{\rm qtr}$、所有谱、$\delta(V)$ 以及 $F_{\rm rec}(V)$ 在扇区内局域幺正下逐点不变。其中 $F_{\rm rec}(V)$ 的不变性是因为解码器可以吸收 $W_a^R$，而对 $B$ 求迹与 $W_a^B$ 无关。

所以一切结论与基底无关。附带结论：**F-符号不进入任何 Haar 平均量**，只通过"哪些算符是局域的"来实现 $\mathcal A_R$。

---

## 2. A 部分：熵的校准（$k=1$）

### 2.1 扇区概率 [证明]
单位球面上均匀分布的 $\psi\in\mathbb C^D$ 满足 $(|\psi_i|^2)_i\sim\mathrm{Dirichlet}(1,\dots,1)$。把分量按块求和，得到
$$(p_1,p_\tau)\sim\mathrm{Dirichlet}(M_1,M_\tau),\qquad p_\tau\sim\mathrm{Beta}(M_\tau,M_1),\qquad \mathbb Ep_a=\frac{M_a}{D},\qquad \operatorname{Var}p_a=\frac{M_1M_\tau}{D^2(D+1)}.$$
归一化的块 $X_a/\sqrt{p_a}$ 是 $\mathbb C^{M_a}$ 中相互独立的均匀单位向量，且与 $p$ 独立：在 Gaussian 表示 $\psi=g/|g|$ 中，块的范数与块的方向独立。当 $N\to\infty$ 时 $p_a\to d_a^2/\mathcal D^2$，即 $p_\tau^\infty=\varphi/\sqrt5\approx0.7236$，涨落为 $O(D^{-1/2})$。Monte Carlo 核验了 $\mathbb Ep_1$ 与 $\operatorname{Var}p_1$（`partA.py`，7 个切分）。

### 2.2 两种熵的关系 [证明]
$$S_{\rm alg}=H(p)+\sum_ap_aS(\rho_{R,a}),\qquad S_{\rm qtr}=-\sum_ad_a\operatorname{Tr}\Big[\tfrac{X_aX_a^\dagger}{d_a}\ln\tfrac{X_aX_a^\dagger}{d_a}\Big]=S_{\rm alg}+\sum_ap_a\ln d_a .$$
其中 $\widetilde{\operatorname{Tr}}\widetilde\rho_R=\sum_ap_a=1$，所以 $\widetilde\rho_R$ 已正确归一化。由于 $X_aX_a^\dagger$ 与 $X_a^TX_a^{*}$ 同谱，且 $R$、$B$ 共享同一个 $p$，有 $S_{\rm alg}(R)=S_{\rm alg}(B)$ 以及 $S_{\rm qtr}(R)=S_{\rm qtr}(B)$（数值上也逐样本成立）。两个上界为 $S_{\rm qtr}(R)\le\ln\sum_ad_am_a=N_R\ln\varphi$（**精确**），以及 $S_{\rm alg}(R)\le\ln F_{N_R+1}$。

### 2.3 精确有限维公式 [证明]
Page 公式（$m\le n$）为 $S_{\rm P}(m,n)=H_{mn}-H_n-\frac{m-1}{2n}$，出处 Page 1993；证明见 Foong–Kanno 1994、Sánchez-Ruiz 1995。Dirichlet 分布的熵平均为 $\mathbb E[-p_a\ln p_a]=\frac{M_a}{D}(H_D-H_{M_a})$。结合 2.1 的独立性：
$$\boxed{\mathbb E S_{\rm alg}=H_D-\sum_a\frac{M_a}{D}\Big[H_{\max(m_a,n_a)}+\frac{\min(m_a,n_a)-1}{2\max(m_a,n_a)}\Big],\qquad \mathbb E S_{\rm qtr}=\mathbb E S_{\rm alg}+\frac{M_\tau}{D}\ln\varphi.}$$
这里 $\mathbb EH(p)=H_D-\sum_a\frac{M_a}{D}H_{M_a}$ 与 Page 项中的 $H_{M_a}$ 相消了。

| $(N_R,N_B)$ | $D$ | $\mathbb ES_{\rm alg}$ 精确 | MC（$2\times10^4$ 样本） | $\mathbb ES_{\rm qtr}$ 精确 | MC |
|---|---|---|---|---|---|
| (3,2) | 3 | 0.50000 | 0.49907±0.00133 | 0.82081 | 0.82031±0.00103 |
| (2,4) | 5 | 0.58333 | 0.58185±0.00090 | 0.87206 | 0.86988±0.00079 |
| (4,4) | 13 | 1.14167 | 1.14170±0.00097 | 1.47482 | 1.47528±0.00090 |
| (5,3) | 13 | 0.92372 | 0.92377±0.00078 | 1.29389 | 1.29451±0.00064 |
| (3,6) | 21 | 0.98334 | 0.98386±0.00053 | 1.34998 | 1.35029±0.00043 |
| (6,6) | 89 | 2.04869 | 2.04843±0.00040 | 2.39473 | 2.39434±0.00037 |

### 2.4 $N\to\infty$、$\Delta$ 固定时的展开 [证明；系数已高精度数值核验]
推导步骤：
1. 将 $H_x=\ln x+\gamma+\frac1{2x}-\frac1{12x^2}+\dots$ 代入 2.3，得到（$\Delta\ge0$）$\mathbb ES_{\rm alg}=\ln D+\frac1{2D}-\sum_ap_a[\ln m_a+\frac{n_a}{2m_a}-\frac1{12m_a^2}]+O(D^{-2})$，其中 $p_a=M_a/D$ 取精确值。
2. 用 Binet 公式 $F_n=(\varphi^n-(-\varphi)^{-n})/\sqrt5$ 展开。
3. $O(\varphi^{-2N_B})$ 阶的奇偶项在 $S_{\rm qtr}$ 中精确相消，原因是 $p_\tau^\infty=\varphi^2p_1^\infty$，即 $p_a^\infty\propto d_a^2$。它在 $S_{\rm alg}$ 中残留为 $-\delta p_\tau\ln\varphi$。

结果即 §0 第 2 条中的两式（$\Delta<0$ 时由 $R\leftrightarrow B$ 对称得到）。

- **首个非平凡修正**是 $O(1)$ 的 $-\tfrac12\varphi^{-|\Delta|}=-\tfrac12\min(d_R,d_B)/\max(d_R,d_B)$。它正是把整数维换成量子维数后的 Page 修正。
- $S_{\rm alg}$ 还多一个常数中心项 $-p_\tau^\infty\ln\varphi=-0.348211$。
- 随后是 $\varphi^{-N}$ 阶的指数小修正。只有 $S_{\rm alg}$ 带 Fibonacci 奇偶振荡 $(-1)^{N_B},(-1)^{N_R}$。

核验（`partA_asym.py`，40 位精度）：$(\mathbb ES-{\rm leading})\varphi^N$ 在 $N=54$ 时，$\Delta=0,1,2,3$ 的值为 2.412023、2.181695、2.039345、1.951367。与 $\frac{\sqrt5\varphi}{2}(1+\varphi^{-\Delta}/3)$ 的差 $<10^{-6}$。$S_{\rm alg}$ 在偶/奇 $N_B$（$\Delta=0$）时的系数 3.108439 与 1.715606 也与公式一致。

### 2.5 Replica moments [证明]
对一般 $n$，精确有（$\gamma_n$ 为循环置换，$D^{\bar n}=D(D+1)\cdots(D+n-1)$）：
$$\mathbb E\operatorname{Tr}(\rho_R^{\rm alg})^n=\frac{1}{D^{\bar n}}\sum_a\sum_{\sigma\in S_n}m_a^{\#(\gamma_n^{-1}\sigma)}n_a^{\#(\sigma)},\qquad \mathbb E\,\widetilde{\operatorname{Tr}}\,\widetilde\rho_R^{\,n}=\frac{1}{D^{\bar n}}\sum_ad_a^{1-n}\sum_{\sigma\in S_n}(\cdots).$$
**二阶：**
$$\mathbb E\operatorname{Tr}(\rho^{\rm alg}_R)^2=\sum_a\frac{M_a(m_a+n_a)}{D(D+1)},\qquad \mathbb E\widetilde{\operatorname{Tr}}\widetilde\rho_R^2=\sum_a\frac{M_a(m_a+n_a)}{d_a\,D(D+1)}.$$
它们已由 MC 核验（`partA.py`）。例如 (6,6) 处精确值 0.15905 / 0.11022，MC 为 0.15914 / 0.11029。渐近地有
$$\mathbb E\widetilde{\operatorname{Tr}}\widetilde\rho^2=(d_R^{-1}+d_B^{-1})(1+O(\varphi^{-N})),\qquad \mathbb E\operatorname{Tr}(\rho^{\rm alg})^2=\tfrac{2\varphi}{\sqrt5}(d_R^{-1}+d_B^{-1})(1+\dots)$$
（`partA_asym.py`：(20,20) 处比值为 $1-1.6\times10^{-8}$）。

**平面结构：** 主阶只有非交叉置换贡献，满足 $\#(\gamma^{-1}\sigma)+\#\sigma=n+1$。代入 $m_a\simeq d_ad_R/\mathcal D^2$、$n_a\simeq d_ad_B/\mathcal D^2$、$D\simeq d_Rd_B/\mathcal D^2$，每一项的扇区和都给出同一个因子：
$$\mathbb E\widetilde{\operatorname{Tr}}\widetilde\rho^{\,n}=\sum_{\sigma\in NC(n)}d_R^{\#(\gamma^{-1}\sigma)-n}d_B^{\#\sigma-n}\,(1+O(\varphi^{-N})),\qquad \mathbb E\operatorname{Tr}(\rho^{\rm alg})^n=\kappa_n\times(\text{同一式}),\quad \kappa_n=\sum_ap_a^\infty d_a^{\,n-1}.$$
也就是说，**qtr 约定的所有 Rényi moment 都是以量子维数 $d_R,d_B$ 为"维数"的普通 Page moment**，alg 约定只差常数因子 $\kappa_n$。因此 $S_n^{\rm alg,ann}=S_n^{\rm qtr,ann}-\frac{\ln\kappa_n}{n-1}$，当 $n\to1$ 时回到 $-\sum_ap_a\ln d_a$，与 2.2 的精确恒等式一致。$\kappa_2=2\varphi/\sqrt5$。

### 2.6 $\mathbb E\ln Z$、$\ln\mathbb EZ$ 与归一化顺序的区分 [证明 + 数值]
- **归一化后再平均，von Neumann 情形。** 对归一化态 $\operatorname{Tr}\rho=1$，有 $-\partial_n\ln\mathbb E\operatorname{Tr}\rho^n|_{n=1}=\mathbb ES$，所以在 $n\to1$ 处 annealed 与 quenched 相等。
- **未归一化 Gaussian 态的 replica。** 对 $g\sim\mathcal{CN}(0,I/D)$，$|g|$ 与方向独立，$\mathbb E|g|^{2n}=\Gamma(D+n)/(\Gamma(D)D^n)$，从而
$$-\partial_n\ln\mathbb E\operatorname{Tr}\rho_g^n\big|_{n=1}=\mathbb ES-\big[\psi(D+1)-\ln D\big]=\mathbb ES-\tfrac1{2D}+O(D^{-2}).$$
  数值核验（`check_gauss_replica.py`，$2\times10^5$ 样本）：(4,4) 处 MC 1.10416，公式 1.10370。
- **Rényi-2。** $\mathbb E[-\ln Z]\ge-\ln\mathbb EZ$（Jensen）。数值上 (4,4)：quenched 0.97030，annealed 0.95551；(6,6)：1.84137 对 1.83853，差距随 $D$ 缩小。
- **归一化与未归一化的比值。** $\mathbb E\operatorname{Tr}\rho_\psi^2=\mathbb EZ_2(g)/\mathbb E[Z_1(g)^2]$ 是精确的（范数与方向独立）。但 $\mathbb EZ_2(g)/(\mathbb EZ_1(g))^2=(1+1/D)\,\mathbb E\operatorname{Tr}\rho_\psi^2$。

四种做法都在 $O(1/D)=O(\varphi^{-N})$ 阶上彼此不同，在有限 $N$ 下必须区分。

---

## 3. B 部分：信息恢复

### 3.1 信道结构 [证明]
令 $J:|a,r,b\rangle\mapsto|a,r\rangle_{R'}|a,b\rangle_{B'}$，其中 $R'=\bigoplus_a\mathcal H^a_R$，$B'=\bigoplus_a\mathcal H^a_B$。则 $\mathcal N_R=\operatorname{Tr}_{B'}(JV\cdot V^\dagger J^\dagger)$，$\mathcal N_B=\operatorname{Tr}_{R'}(\cdots)$，二者是**互补信道**，扇区标签作为共享中心被复制到两边。

$\mathcal N_R$ 的输出是块对角的，所以任何 CPTP 解码器都等价于"测量 $a$，再施加 $\mathcal D_a$"，不需要也无法利用扇区间的相干。由此逐 $V$ 地有：
$$F_{\rm rec}=\sum_aq_aF_a,\qquad q_a=\tfrac1k\operatorname{Tr}(V^\dagger P_aV),\qquad F_a=\max_{\mathcal D_a}\langle\Phi|(\mathrm{id}\otimes\mathcal D_a)(\hat\rho^{(a)}_{QR})|\Phi\rangle,$$
$$F_{\rm rec}=\max_{\sigma_B}F(\rho_{QB},\tfrac{I}{k}\otimes\sigma_B)=\tfrac1k\min\{\operatorname{Tr}\sigma_R:\ I_Q\otimes\sigma_R\ge\rho_{QR},\ \sigma_R\in\mathcal A_R\}=\tfrac1ke^{-H_{\min}(Q|R)}.$$
- 第二个等号用 Uhlmann 定理，两个方向都成立：
  - "≤"：解码后纯化态与 $|\Phi\rangle\otimes|\xi\rangle$ 的重叠；
  - "≥"：由 Uhlmann 等距构造解码器。
- 第三个等号是 König–Renner–Schaffner（IEEE TIT 55, 4337 (2009)）的操作意义定理，逐扇区应用。

**三路数值互验**（`check_frec.py`，10 个切分，含不等切分与 $k=2,3,4$）结果完全一致，到 $10^{-8}$：
1. 带对偶间隙证书的 $\max_\sigma$ 形式（L-BFGS 求解，证书为 $h(\sigma)\le h^\*\le h/2+\lambda_{\max}(\nabla h)$）；
2. 直接 CPTP SDP（上面的 $H_{\min}$ 形式，Clarabel 求解）；
3. **显式构造的扇区 Uhlmann 解码器作为 CPTP 映射作用后算出的保真度**。例如 (4,5)、$k=2$ 时三者都给出 0.66405867。

### 3.2 精确二阶矩 [证明；已精确低阶检验]
令 $W_e=\frac1{D^2-1}$，$W_s=\frac{-1}{D(D^2-1)}$（$n=2$ 的 Weingarten 函数，出处 Collins–Śniady, CMP 264, 773 (2006)），$\alpha_a=m_a^2n_a$，$\beta_a=m_an_a^2$。$V$ 的 $k$ 列按正交归一处理，没有当作独立向量。逐扇区有：
$$\mathbb E\operatorname{Tr}(\rho^{(a)}_{QB})^2=\tfrac1k\big[(\alpha_a+k\beta_a)W_e+(k\alpha_a+\beta_a)W_s\big],\qquad \mathbb E\operatorname{Tr}(\rho^{(a)}_{B})^2=\tfrac1k\big[(k\alpha_a+\beta_a)W_e+(\alpha_a+k\beta_a)W_s\big],$$
$$\mathbb E\operatorname{Tr}(C^{(a)})^2=\mathbb E\operatorname{Tr}(\rho^{(a)}_{QB})^2-\tfrac1k\mathbb E\operatorname{Tr}(\rho^{(a)}_B)^2=(1-k^{-2})(\beta_aW_e+\alpha_aW_s)=(1-k^{-2})\frac{m_an_a(n_a-m_a/D)}{D^2-1}.$$
$C_{QB}$ 是扇区块对角的，因此 $\mathbb E\operatorname{Tr}C_{QB}^2=\sum_a\mathbb E\operatorname{Tr}(C^{(a)})^2$。

核验（`partB_moments.py`）：
- 与独立实现的**暴力 Weingarten 张量收缩**比较，6 组 $(N_R,N_B,k)$、12 个扇区全部一致到 $10^{-10}$；
- 与 Monte Carlo（每组 4000 个 Haar 等距）一致。

几个检查点：
- 单扇区且 $n=1$ 时结果为 0（此时 $D=m$）；多扇区中 $n_a=1$ 时剩下 $(1-k^{-2})m_a(1-m_a/D)/(D^2-1)$，它恰好等于 §3.4 的标签泄露矩除以 $k^2$。此时 $C^{(a)}=(Q_a^T-q_aI)/k$，两个独立推导相互印证；
- $m=1$ 的单扇区情形给出 $1-k^{-2}$，与 $\rho_{QB}$ 为纯态时的直接计算一致；
- 渐近地 $\mathbb E\operatorname{Tr}C_{QB}^2\simeq(1-k^{-2})\frac{2\varphi}{\sqrt5}\varphi^{-N_R}$。

**若误把列当作独立向量**，推导会给出 $\frac{k-1}{k^2}\sum_a\big[\frac{k\beta_a-\alpha_a}{D^2}+\frac{\alpha_a+\beta_a}{D(D+1)}\big]$。它与正确结果主阶相同，但在 $O(1/D)$ 阶上不同。

### 3.3 关于 $F_{\rm rec}$ 的非空洞界 [证明]
- **(B1)** 对每个 $V$：$(1-\delta)^2\le F_{\rm rec}\le1-\delta^2/4$。
  - 下界：$F_{\rm rec}\ge F(\rho_{QB},I/k\otimes\rho_B)$，再用 Fuchs–van de Graaf 不等式。
  - 上界：取最优 $\sigma^\*$，则 $T(\rho_{QB},I/k\otimes\sigma^\*)\le\sqrt{1-F_{\rm rec}}$；再用三角不等式和 $\operatorname{Tr}_Q$ 下的单调性，得 $\delta\le2\sqrt{1-F_{\rm rec}}$。

  因此 $F_{\rm rec}\to1\iff\delta\to0$。
- **(B2) 解耦下界。** 由 $\|C^{(a)}\|_1\le\sqrt{kn_a}\|C^{(a)}\|_2$ 与 Jensen 不等式，
$$\mathbb E\delta\le\bar\delta:=\tfrac12\sum_a\sqrt{kn_a\,\mathbb E\operatorname{Tr}(C^{(a)})^2},\qquad \mathbb EF_{\rm rec}\ge(1-\bar\delta)^2\quad(\bar\delta\le1).$$
  渐近地 $\bar\delta\to\frac12\sqrt{1-k^{-2}}\,\varphi^{-\Xi/2}$。这里用了 $\sum_a\sqrt{m_an_a^3}/D\to\varphi^{-\Delta/2}$。
- **(B3) 碰撞下界。** 由 Hölder 不等式，$(\operatorname{Tr}\sqrt\rho)^2\ge1/\operatorname{Tr}\rho^2$，故 $F_a\ge F(\hat\rho^{(a)},I/(kn_a))\ge1/(kn_a\operatorname{Tr}\hat\rho^{(a)2})$。再对 $(x,y)\mapsto x^3/y$ 的凸性用 Jensen：
$$\mathbb EF_{\rm rec}\ \ge\ \sum_a\frac{(M_a/D)^3}{kn_a\,\mathbb E\operatorname{Tr}(\rho^{(a)}_{QB})^2}\ \xrightarrow{N\to\infty}\ \sum_a\frac{p_a}{1+c_a}=\frac1{1+\varphi^{-\Xi}}.$$
  这是精确的有限 $N$ 定理，只依赖 $\Xi$ 与精确二阶矩。它给出 $1-F\lesssim\varphi^{-\Xi}$，比解耦界的 $\varphi^{-\Xi/2}$ 更紧。
- **(B4) 逆向界。** 对每个 $V$，$q_aF_a\le\operatorname{rank}(\rho^{(a)}_{QB})\,\lambda_{\max}(\rho^{(a)}_B)/k$，由 $(\operatorname{Tr}\sqrt X)^2\le\operatorname{rank}X\cdot\operatorname{Tr}X$ 得到。等价的对偶形式是在 $H_{\min}$-SDP 中取 $\sigma_R=\lambda_{\max}I$。渐近形式见 §4.3。

### 3.4 扇区标签向 B 的泄露 [证明]
B 总能读出 $a$，泄露的是 $a$ 的分布对输入的依赖：$\Pr(a\,|\,\rho)=\operatorname{Tr}(Q_a\rho)$，$Q_a=V^\dagger P_aV$。精确二阶矩为
$$\mathbb E\big\|Q_a-\tfrac{\operatorname{Tr}Q_a}{k}I\big\|_2^2=\frac{(k^2-1)M_a(D-M_a)}{D(D^2-1)}\simeq\frac{(k^2-1)p_a(1-p_a)}{D},$$
已由 MC 核验。B 从标签上能获得的最大区分优势 $\le\sum_a\|Q_a-q_aI\|_\infty$，期望为 $O(k/\sqrt D)=O(k\varphi^{-N/2})$。

它对 $F_{\rm rec}$ 的影响被 §4.1 的约化定理精确控制。在极小系统中影响可见：(4,3)、$D=8$ 时 $F^{\rm Fib}=0.771\pm0.008$，扇区 Haar 混合为 $0.813\pm0.007$；到 $D\gtrsim10^2$ 时差异已不可分辨（§4.4）。

### 3.5 哪些维数条件保证什么
- **足以保证典型的高保真恢复：** 对每个权重不可忽略的扇区都有 $kn_a/m_a\to0$，也就是 $\Xi\to+\infty$。此时 $\mathbb E(1-F_{\rm rec})\le\varphi^{-\Xi}/(1+\varphi^{-\Xi})+o(1)$（B3）。再结合 §4.3 的集中性，几乎每个 $V$ 都成立。$D\gg k$ 自动满足，保证标签泄露可忽略。
- **只能保证某种"看起来接近最大混合"，不足以恢复：**
  1. $\mathbb E_V\rho_{QB}=\frac Ik\otimes\mathbb E_V\rho_B$ 对任何维数都成立（Haar 不变性），是空洞陈述。
  2. $\mathbb E\operatorname{Tr}C^2\simeq(1-k^{-2})\frac{2\varphi}{\sqrt5}\varphi^{-N_R}\to0$ 只要 $N_R\to\infty$ 就成立，即便 $N_B\gg N_R$、恢复完全失败。原因是 HS 范数小不等于迹范数小，缺了 $kn_a$ 的维数因子。
  3. $\rho_B$ 在每个扇区内接近 $q_aI/n_a$ 只需要 $n_a\ll km_a$，即 $\Xi>-2\log_\varphi k$。所以在窗口 $-2\log_\varphi k<\Xi<0$ 内，B 看起来是热的，但恢复失败。这就是单态 Page 条件 $n\ll m$ 与恢复条件 $kn\ll m$ 之间的差别。
  4. $q_a\approx p_a$ 只需要 $D\gg1$，与恢复无关。

---

## 4. C 部分：临界窗口

### 4.1 扇区约化定理 [证明]
取 $V=G(G^\dagger G)^{-1/2}$，其中 $G$ 是 $D\times k$ 复 Gaussian 矩阵。记 $G_a=P_aG=U_aH_a$，$U_a=G_a(G_a^\dagger G_a)^{-1/2}$，$H_a=(G_a^\dagger G_a)^{1/2}$。则
$$P_aV=U_aT_a,\qquad T_a=H_a(G^\dagger G)^{-1/2},$$
其中：
- $U_a$ 是 $\mathbb C^k\to\mathcal H^a_R\otimes\mathcal H^a_B$ 的 **Haar 等距**（要求 $M_a\ge k$）；
- $U_a$ 与 $T_a$ **独立**：由左不变性，$U_a$ 与 $G_a^\dagger G_a$ 独立，而 $G^\dagger G$ 只依赖 $H_a$ 与其他块；
- $T_a^\dagger T_a=Q_a$。

记 $u_a$ 为 $T_a$ 的极分解幺正部分。解码器能吸收 $u_a$，所以 $F_{\rm opt}(U_au_a)=F_{\rm opt}(U_a)$。又因为 $F_{\rm opt}$ 在迹距离下是 $\tfrac12$-Lipschitz 的，
$$q_a\,|F_a-F_{\rm opt}(U_a)|\le\|Q_a-q_aI\|_{2}/\sqrt k .$$
对 $V$ 取期望，并用 $q_a$ 与 $U_a$ 的独立性：
$$\boxed{\Big|\mathbb EF^{\rm Fib}_{\rm rec}-\sum_a\frac{M_a}{D}\,\bar F^{\rm Haar}_k(m_a,n_a)\Big|\le\sum_a\sqrt{\frac{(k^2-1)M_a(D-M_a)}{kD(D^2-1)}}}$$
其中 $\bar F^{\rm Haar}_k(m,n)$ 是普通 Haar 码 $\mathbb C^k\to\mathbb C^m\otimes\mathbb C^n$ 从 $m$ 因子恢复时的平均最优保真度。

也就是说，**Fibonacci 码在可操作层面上就是按 $p_a$ 加权的普通 Haar 码混合**，维数取融合重数。

### 4.2 比值普适性与条件极限定理 [证明]
由 1.1，
$$c_1=\frac{kF_{N_B-1}}{F_{N_R-1}},\qquad c_\tau=\frac{kF_{N_B}}{F_{N_R}},\qquad c_a=k\varphi^{-\Delta}\big(1+O(\varphi^{-2\min(N_R,N_B)})\big)\to\varphi^{-\Xi}.$$
修正项的符号随 Fibonacci 奇偶交替。**推论：** 若普通 Haar 码的交叉函数 $f_k(c)=\lim\bar F^{\rm Haar}_k(m,n)$（$m,n\to\infty$，$kn/m\to c$）存在、局部一致且连续，则 $\mathbb EF^{\rm Fib}_{\rm rec}\to f_k(\varphi^{-\Xi})$。这个极限：
- 与扇区权重 $p_a$ 无关，因为两个扇区的 $c_a$ 极限相同；
- 与 $\widetilde{\operatorname{Tr}}$ 中的 $d_a$ 无关；
- 与 F-符号无关。

$d_\tau$ 仅通过 $\log_\varphi$ 进入。$f_k$ 的存在性只有数值证据（§4.4），见 [猜想]。

**普适性的来源**有两种说法：
1. 从 moment 看（§2.5）：两个竞争的鞍点，例如 $m_a^2n_a$ 与 $km_an_a^2$，携带相同的扇区因子 $\sum_ad_a^3$，所以二者的比值与扇区无关。
2. 从 Perron–Frobenius 看：对任何单生成对象 $x$、融合图连通且非周期的融合范畴，都有 $n_a/m_a\to d_x^{-\Delta}$ 且与 $a$ 无关。因此结论应推广为 $\Xi_x=\Delta-\log_{d_x}k$ [猜想，论证概要]。

### 4.3 只依赖 $\Xi$ 的上下界与集中性 [证明*]
对固定 $k$、$N\to\infty$、$\Xi$ 固定：
- **下界：** $\tfrac1{1+\varphi^{-\Xi}}$ 来自 (B3)，**完全证明**。$\big(\mathbb E_{{\rm MP}(\varphi^{-\Xi})}\sqrt\lambda\big)^2$ 来自取 $\sigma=I/n_a$，即 $F\ge\sum_a(\operatorname{Tr}\sqrt{\rho^{(a)}_{QB}})^2/(kn_a)$ 对每个 $V$ 成立；再用 Marchenko–Pastur 定律（Math. USSR-Sb. 1, 457 (1967)），以及 $k\ll D$ 时 Stiefel 与 Gaussian 的差 $G^\dagger G/D\to I$。由 MP 律 $\mathbb E\lambda^2=1+c$ 与 Hölder 不等式，MP 下界 $\ge$ 碰撞下界。平凡下界 $1/k^2$ 由常值解码器给出。
- **上界：** 由 (B4) 与 Wishart 型矩阵最大特征值收敛到 MP 边缘 $(1+\sqrt y)^2$（Geman, Ann. Probab. 8, 252 (1980)；Yin–Bai–Krishnaiah, PTRF 78, 509 (1988)），其中 $y_a=n_a/(km_a)\to c/k^2$，得到
$$\limsup F_{\rm rec}\le\sum_ap_a\Big(\tfrac1{\sqrt{c_a}}+\tfrac1k\Big)^2=\big(\varphi^{\Xi/2}+\tfrac1k\big)^2 .$$
  上界在 $\Xi\to-\infty$ 时趋于 $1/k^2$（偏离量 $O(\varphi^{\Xi/2})$），下界在 $\Xi\to+\infty$ 时趋于 1（$1-F\le\varphi^{-\Xi}$）。两者都只依赖 $\Xi$，给出**相同的临界缩放**：转折发生在 $\Xi=O(1)$，且不含任何依赖 $d_a$ 的偏移。
- **集中性：** $|F(V)-F(V')|\le\|V-V'\|_{HS}/\sqrt k$，由迹距离 Lipschitz 性与 $\|\psi\psi^\dagger-\psi'\psi'^\dagger\|_1\le2\|\psi-\psi'\|$ 得到。由 $U(D)$ 上的 Gaussian 集中（E. Meckes, *The Random Matrix Theory of the Classical Compact Groups*, CUP 2019, §5.3，Thm 5.17），$\Pr(|F-\mathbb EF|\ge t)\le2e^{-kDt^2/12}$；只用到 $\exp(-c\,kDt^2)$ 这一形式。因此**极限分布是点质量**。数值上 sd$(F)$ 随 $D^{-1/2}$ 下降：$D=233,1597,4181$ 时为 0.012、0.0054、0.0028。

### 4.4 显式数值恢复：$k=2$，含不等切分 [数值]
`scan_k2.py`：对每个 $(N_R,N_B)$ 取 48–120 个 Haar 等距，每个样本都用带证书的最优解码（间隙 $<10^{-7}$）求出 $F_{\rm rec}$。下表取每个 $\Delta$ 的最大 $N$：

| $\Delta$ | $\Xi$ | $(N_R,N_B)$ | $D$ | $\mathbb EF_{\rm rec}$ | sd | $F_1$ / $F_\tau$ | 普通 Haar（同 $c$，$n=128$） | 碰撞下界（精确） | $(1-\bar\delta)^2$（精确） | 渐近上界 $(\varphi^{\Xi/2}+\frac12)^2$ | $1-\frac3{16}\varphi^{-\Xi}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| -4 | -5.44 | (8,12) | 4181 | **0.39328**±0.00040 | 0.0028 | 0.3939 / 0.3930 | — | 0.0680 | 0.0000 | 0.5930 | — |
| -3 | -4.44 | (9,12) | 6765 | **0.43580**±0.00031 | 0.0021 | 0.4358 / 0.4358 | 0.43362±0.00050（$c$=8.53） | 0.1056 | 0.0000 | 0.7116 | — |
| -2 | -3.44 | (10,12) | 10946 | **0.48945**±0.00025 | 0.0018 | 0.4900 / 0.4893 | 0.49073±0.00042（$c$=5.22） | 0.1604 | 0.0001 | 0.8780 | — |
| -1 | -2.44 | (11,12) | 17711 | **0.55888**±0.00026 | 0.0018 | 0.5587 / 0.5589 | 0.56086±0.00060（$c$=3.20） | 0.2361 | 0.0489 | 1.0000 | — |
| +0 | -1.44 | (12,12) | 28657 | **0.64667**±0.00019 | 0.0013 | 0.6467 / 0.6467 | 0.64692±0.00022 | 0.3333 | 0.1503 | 1.0000 | — |
| +1 | -0.44 | (13,12) | 46368 | **0.75301**±0.00015 | 0.0011 | 0.7528 / 0.7531 | 0.75320±0.00014 | 0.4472 | 0.2689 | 1.0000 | — |
| +2 | +0.56 | (14,12) | 75025 | **0.85198**±0.00008 | 0.0006 | 0.8519 / 0.8520 | 0.85203±0.00015 | 0.5669 | 0.3863 | 1.0000 | 0.85676 |
| +3 | +1.56 | (15,12) | 121393 | **0.91001**±0.00005 | 0.0004 | 0.9101 / 0.9100 | 0.90979±0.00009 | 0.6793 | 0.4935 | 1.0000 | 0.91147 |
| +4 | +2.56 | (16,12) | 196418 | **0.94472**±0.00004 | 0.0002 | 0.9447 / 0.9447 | 0.94477±0.00005 | 0.7741 | 0.5869 | 1.0000 | 0.94529 |
| +5 | +3.56 | (17,12) | 317811 | **0.96602**±0.00003 | 0.0002 | 0.9660 / 0.9660 | 0.96600±0.00003 | 0.8472 | 0.6660 | 1.0000 | 0.96619 |
| +6 | +4.56 | (17,11) | 196418 | **0.97901**±0.00003 | 0.0002 | 0.9790 / 0.9790 | — | 0.8997 | 0.7318 | 1.0000 | 0.97910 |
| +7 | +5.56 | (18,11) | 317811 | **0.98706**±0.00001 | 0.0001 | 0.9871 / 0.9870 | — | 0.9356 | 0.7856 | 1.0000 | 0.98708 |

要点：
1. 沿每条 $\Delta$ 固定的整数子序列，$F$ 随 $N$ 收敛，有限尺寸偏差在 $D\gtrsim10^2$ 后低于统计误差。
2. 大 $N$ 下两个扇区的 $F_1$ 与 $F_\tau$ 相等（小 $N$ 下因 $c_1\ne c_\tau$ 而不同），这直接验证了 §4.2。
3. 不等切分的数据（$\Delta=\pm1,\dots$）落在同一条仅依赖 $\Xi$ 的曲线上。

**与普通 Haar 码的普适性比较**（`universality.py`）：
- 在同一 $c$ 下，$n=128$ 时的普通 Haar 码给出：
  - $\Delta=0$（$c=2$）：0.64692±0.00022，Fibonacci $(12,12)$ 为 0.64667±0.00019；
  - $\Delta=2$：0.85203±0.00015，Fibonacci $(14,12)$ 为 0.85198±0.00008；
  - $\Delta=3$：0.90979±0.00009，Fibonacci $(15,12)$ 为 0.91001±0.00005；
  - $\Delta=5$：0.96600±0.00003，Fibonacci $(17,12)$ 为 0.96602±0.00003。
- 逐切分地，$F^{\rm Fib}$ 与 $\sum_ap_a\bar F^{\rm Haar}(m_a,n_a)$ 在 $(9,9)$、$(10,8)$、$(7,9)$、$(11,10)$ 上分别为 0.6469/0.6468、0.8520/0.8523、0.4902/0.4895、0.7533/0.7530。

**两端渐近律。**
- **$\Xi\gg1$ [数值 + 启发式推导]：** 在平坦点附近，Bures 度规的二阶展开给出 $1-F^\*\simeq\frac{kn}{4}\min_\sigma\|\rho-\frac Ik\otimes\sigma\|_2^2=\frac{kn}{4}\operatorname{Tr}C^2$。代入 3.2 的精确 $\mathbb E\operatorname{Tr}C^2$（含扇区求和）得
$$1-\mathbb EF_{\rm rec}\simeq\tfrac14(1-k^{-2})\,\varphi^{-\Xi}.$$
  $k=2$ 时：$\Xi=2.56$ 预测 0.94529，实测 0.94472；$\Xi=3.56$ 预测 0.96619，实测 0.96602；$\Xi=4.56$ 预测 0.97910，实测 0.97901；$\Xi=5.56$ 预测 0.98708，实测 0.98706。$k=3$ 时：$\Xi=4.72$ 预测 0.97702，实测 0.9769。
  同一展开给出 $I(Q\!:\!B)\simeq\frac{kn}{2}\operatorname{Tr}C^2$，因此在该区域 $1-F\simeq\tfrac12 I(Q\!:\!B)$。
- **$\Xi\ll0$ [数值 + 化约]：** 当 $c\gg k^2$ 时 $\rho_{QR}\simeq\frac1{km}(I+\frac{k}{\sqrt c}W)$，$W$ 为 GUE。代入 $H_{\min}$-SDP，在 $k/\sqrt c$ 的主阶上化约为
$$F_{\rm rec}-\tfrac1{k^2}\simeq\tfrac{s_k}{k}\varphi^{\Xi/2},\qquad s_k=\lim_m\tfrac1m\min\{\operatorname{Tr}S:\ I_k\otimes S\ge W\}\in(0,2].$$
  $s_k\le2$ 来自 $S=\lambda_{\max}(W)I$。数值上（`s_k.py`，Clarabel）$s_2=0.986\pm0.022\ (m=8)$、$1.0008\pm0.0074\ (m=16)$、$0.9946\pm0.0075\ (m=32)$，与 $s_2=1$ 相容。$k=3$ 及 $m=64$ 的 SDP 因内存不足未完成。$k=2$ 的扫描数据在 $\Xi=-5.44,-4.44,-3.44$ 给出 $2\sqrt c\,(F-\tfrac14)=1.06,1.08,1.10$，并随 $c$ 增大向 $s_2$ 减小，与 $O(1/c)$ 次级修正一致。

### 4.5 让 $k$ 随 $N$ 增长 [数值；猜想 + 证明概要]
$k\in\{3,5,8,13,21\}$，$N_B\in\{7,9\}$，每点 24 个样本（`scan_k.py`）。

| $k$ | $(N_R,N_B)$ | $\Xi$ | $\mathbb EF_{\rm rec}$ | MP 下界（$\sigma=I/n$） | 一阶证书上界 | $f_\infty(\varphi^{-\Xi})$ | $F-f_\infty$ |
|---|---|---|---|---|---|---|---|
| 3 | (8,9) | -3.28 | 0.3685±0.0008 | 0.1955 | 0.5962 | 0.1953 | +0.1732 |
| 3 | (9,9) | -2.28 | 0.4540±0.0007 | 0.3048 | 0.7090 | 0.3049 | +0.1491 |
| 3 | (10,9) | -1.28 | 0.5684±0.0007 | 0.4637 | 0.8356 | 0.4635 | +0.1049 |
| 3 | (11,9) | -0.28 | 0.7146±0.0004 | 0.6656 | 0.9903 | 0.6654 | +0.0491 |
| 3 | (12,9) | +0.72 | 0.8356±0.0003 | 0.8120 | 1.0714 | 0.8118 | +0.0238 |
| 3 | (13,9) | +1.72 | 0.9001±0.0002 | 0.8868 | 1.0930 | 0.8870 | +0.0132 |
| 3 | (14,9) | +2.72 | 0.9391±0.0001 | 0.9313 | 1.0935 | 0.9311 | +0.0080 |
| 3 | (15,9) | +3.72 | 0.9626±0.0001 | 0.9578 | 1.0838 | 0.9577 | +0.0048 |
| 3 | (16,9) | +4.72 | 0.9769±0.0000 | 0.9740 | 1.0734 | 0.9740 | +0.0029 |
| 5 | (9,9) | -3.34 | 0.2750±0.0005 | 0.1900 | 0.3815 | 0.1899 | +0.0851 |
| 5 | (10,9) | -2.34 | 0.3673±0.0004 | 0.2968 | 0.5025 | 0.2968 | +0.0705 |
| 5 | (11,9) | -1.34 | 0.4971±0.0004 | 0.4523 | 0.6590 | 0.4522 | +0.0449 |
| 5 | (12,9) | -0.34 | 0.6718±0.0003 | 0.6527 | 0.8408 | 0.6527 | +0.0191 |
| 5 | (13,9) | +0.66 | 0.8144±0.0002 | 0.8056 | 0.9573 | 0.8056 | +0.0088 |
| 5 | (14,9) | +1.66 | 0.8885±0.0001 | 0.8835 | 1.0058 | 0.8834 | +0.0051 |
| 5 | (15,9) | +2.66 | 0.9318±0.0001 | 0.9288 | 1.0251 | 0.9290 | +0.0028 |
| 5 | (16,9) | +3.66 | 0.9582±0.0001 | 0.9564 | 1.0327 | 0.9565 | +0.0018 |
| 5 | (17,9) | +4.66 | 0.9743±0.0000 | 0.9732 | 1.0327 | 0.9732 | +0.0011 |
| 8 | (10,9) | -3.32 | 0.2353±0.0002 | 0.1919 | 0.3019 | 0.1919 | +0.0434 |
| 8 | (11,9) | -2.32 | 0.3326±0.0002 | 0.2999 | 0.4209 | 0.2999 | +0.0327 |
| 8 | (12,9) | -1.32 | 0.4739±0.0001 | 0.4563 | 0.5776 | 0.4565 | +0.0174 |
| 8 | (13,9) | -0.32 | 0.6650±0.0001 | 0.6576 | 0.7714 | 0.6576 | +0.0075 |
| 8 | (14,9) | +0.68 | 0.8115±0.0001 | 0.8080 | 0.9017 | 0.8080 | +0.0035 |
| 8 | (15,9) | +1.68 | 0.8868±0.0001 | 0.8849 | 0.9602 | 0.8848 | +0.0021 |
| 8 | (16,9) | +2.68 | 0.9310±0.0000 | 0.9298 | 0.9897 | 0.9298 | +0.0012 |
| 8 | (17,9) | +3.68 | 0.9576±0.0000 | 0.9570 | 1.0034 | 0.9569 | +0.0007 |
| 8 | (18,9) | +4.68 | 0.9740±0.0000 | 0.9736 | 1.0106 | 0.9735 | +0.0005 |
| 13 | (11,9) | -3.33 | 0.2115±0.0001 | 0.1911 | 0.2526 | 0.1911 | +0.0204 |
| 13 | (12,9) | -2.33 | 0.3121±0.0001 | 0.2987 | 0.3686 | 0.2987 | +0.0134 |
| 13 | (13,9) | -1.33 | 0.4614±0.0001 | 0.4548 | 0.5276 | 0.4548 | +0.0065 |
| 13 | (14,9) | -0.33 | 0.6586±0.0001 | 0.6558 | 0.7245 | 0.6557 | +0.0029 |
| 13 | (15,9) | +0.67 | 0.8084±0.0001 | 0.8071 | 0.8640 | 0.8071 | +0.0014 |
| 13 | (16,9) | +1.67 | 0.8852±0.0001 | 0.8844 | 0.9301 | 0.8843 | +0.0009 |
| 13 | (17,9) | +2.67 | 0.9300±0.0000 | 0.9295 | 0.9658 | 0.9295 | +0.0005 |
| 13 | (18,9) | +3.67 | 0.9571±0.0000 | 0.9568 | 0.9855 | 0.9568 | +0.0003 |
| 13 | (19,9) | +4.67 | 0.9736±0.0000 | 0.9734 | 0.9962 | 0.9734 | +0.0002 |
| 21 | (12,9) | -3.33 | 0.2004±0.0001 | 0.1914 | 0.2282 | 0.1914 | +0.0090 |
| 21 | (13,9) | -2.33 | 0.3043±0.0000 | 0.2992 | 0.3408 | 0.2991 | +0.0051 |
| 21 | (14,9) | -1.33 | 0.4581±0.0000 | 0.4556 | 0.5007 | 0.4555 | +0.0027 |
| 21 | (15,9) | -0.33 | 0.6576±0.0001 | 0.6565 | 0.6981 | 0.6564 | +0.0012 |
| 21 | (16,9) | +0.67 | 0.8080±0.0000 | 0.8075 | 0.8420 | 0.8074 | +0.0006 |
| 21 | (17,9) | +1.67 | 0.8848±0.0000 | 0.8845 | 0.9124 | 0.8845 | +0.0004 |
| 21 | (18,9) | +2.67 | 0.9298±0.0000 | 0.9296 | 0.9519 | 0.9296 | +0.0002 |
| 21 | (19,9) | +3.67 | 0.9570±0.0000 | 0.9569 | 0.9744 | 0.9568 | +0.0001 |
| 21 | (20,9) | +4.67 | 0.9735±0.0000 | 0.9734 | 0.9874 | 0.9734 | +0.0001 |

（表中为 $N_B=9$ 的数据；$N_B=7$ 的数据在 `scan_kgrow.json` 中，与之最大相差 $2.7\times10^{-3}$，出现在 $k=3$、$\Xi\approx-2.3$ 处，其余更小。）

- 在固定 $\Xi$ 下，$F_k-f_\infty(\varphi^{-\Xi})$ 随 $k$ 单调减小。
  - 当 $\Xi\gtrsim-1.3$ 时，$k^2(F_k-f_\infty)$ 近似与 $k$ 无关：$\Xi\approx-0.3$ 时 $k=3,5,8,13,21$ 分别为 0.44、0.48、0.48、0.49、0.54；$\Xi\approx+0.7$ 时为 0.21、0.22、0.22、0.23、0.26。这与 $F_k=f_\infty+A(\Xi)/k^2+\dots$ 一致。
  - 当 $\Xi\approx-3.3$（$c\approx5$）时，$k^2(F_k-f_\infty)=1.6,2.1,2.8,3.4,3.9$ 仍随 $k$ 增长，说明在所考察的 $k$ 范围内收敛慢于 $k^{-2}$。这与修正项由 $c/k^2$ 控制、而 $c$ 较大时高阶项仍重要的图像一致，但**尚未被证实**。
- $\sigma=I/n$ 给出的 MP 下界在 $k\to\infty$ 时变紧。
- **猜想：** $\lim_{k\to\infty}\lim_NF_{\rm rec}=f_\infty(\varphi^{-\Xi})=(\mathbb E_{{\rm MP}(c)}\sqrt\lambda)^2$（闭式见 §0 第 6 条，已用 Mathematica 数值积分核验到 30 位）。性质：
  - 级数 $f_\infty(c)=1-\frac c4-\frac{c^2}{64}-\frac{3c^3}{512}-\cdots$；
  - $f_\infty(1)=64/(9\pi^2)$；
  - 对偶性 $f_\infty(c)=c^{-1}f_\infty(1/c)$（MP 律对 $XX^\dagger\leftrightarrow X^\dagger X$ 的对偶）。
- **证明概要（未闭合）：** 对任意 $\sigma_0=I/n$，凹性与 Euler 齐次性给出证书
$$\sqrt{kF}\le\tfrac12\Big[\tfrac{\operatorname{Tr}\sqrt\rho}{\sqrt n}+\sqrt n\,\lambda_{\max}(\operatorname{Tr}_Q\sqrt\rho)\Big].$$
  需要证明 $\lambda_{\max}(T)\le\frac{\operatorname{Tr}T}{n}(1+O(\sqrt{c/k}))$，其中 $T=\operatorname{Tr}_Q\sqrt\rho$。思路是：
  1. $X\mapsto|X|$ 在 HS 范数下 $\sqrt2$-Lipschitz（Araki–Yamagami, CMP 81, 89 (1981)）；
  2. 对固定单位向量 $v$，$\langle v|T|v\rangle$ 的 Gaussian 涨落为 $O(\sqrt{k/(mn)})$，均值为 $\sqrt{k/n}\,\mathbb E\sqrt\lambda$；
  3. 对 $v$ 取 $\varepsilon$-网并用并集界，相对偏差为 $O(\sqrt{n/m})=O(\sqrt{c/k})$。

  尚未写出的部分：常数、Stiefel→Gaussian 比较，以及扇区与 $q_a$ 涨落的逐项控制。表中"一阶证书上界"一列就是这个证书的数值，它随 $k$ 收紧，与概要一致。

---

## 5. 对 C 部分具体问题的回答

**(1) 两种熵的 Page 转折能否单独决定 $F_{\rm rec}$ 的转折？不能。**
- **单态 Page 曲线（$k=1$）** 在两种约定下都在 $\Delta=0$ 处转折。两者只差光滑的中心项 $p_\tau\ln\varphi\to$ 常数，再加指数小振荡，所以转折点位置相同。恢复转折却在 $\Delta=\log_\varphi k$，即 $\Xi=0$，两者相差 $\log_\varphi k$。
- **带参考系的码子空间 Page 曲线** $S(R)_{\rho_{\rm code}}\approx\min(N_R\ln\varphi,\ N_B\ln\varphi+\ln k)$ 在 $\Xi=0$ 处转折，两种约定给出同一位置。所以它只能把**位置**定到 $O(1)$ 精度。
- **窗口内的函数形状由单次熵决定，而不是由 von Neumann 熵决定。** 精确地有 $F_{\rm rec}=e^{-H_{\min}(Q|R)}/k$。在 $\Xi\to-\infty$ 端：
  - von Neumann 信号为 $I_c+\ln k\simeq\frac{k^2-1}{2}\varphi^{\Xi}$（Page 型修正）；
  - 而 $F-1/k^2\simeq\frac{s_k}k\varphi^{\Xi/2}$，二者**指数不同**：$F$ 的偏离量 $\propto\sqrt{I_c+\ln k}$，由谱顶（min-entropy）控制。数值核验（`check_ic_tail.py`，普通 Haar 码，$k=2$）：$c=16,32,64,100$ 时 $I_c+\ln k=0.0866,0.0454,0.0231,0.0148$，对应 $\frac{3}{2c}=0.0938,0.0469,0.0234,0.0150$；同时 $c=16,32$ 时 $F-\frac14=0.129,0.091$，对应 $\frac{1}{2\sqrt c}=0.125,0.088$。
  - 固定 $\Xi$ 时 $F$ 还依赖 $k$，例如 $\Xi=-0.28$ 时 $k=3$ 为 0.715，$k\to\infty$ 极限 $f_\infty=0.665$。
  - 只有在 $\Xi\gg1$ 的微扰区才有 $1-F\simeq\frac12I(Q\!:\!B)$。
- **数值对照**（`entropy_vs_recovery.py`；$N=20$，每点 24 个样本）：

  | $\Delta$ | $k=1$：$S_{\rm alg}(R)$ / $S_{\rm qtr}(R)$ | $k=2$：$S_{\rm alg}(R)$ / $S_{\rm alg}(B)$ | $k=2$：$\Xi$ | $H(Q\vert R)=-I_c$ | $H_{\min}(Q\vert R)$ | $F_{\rm rec}$ |
  |---|---|---|---|---|---|---|
  | −8 | 2.530 / 2.877 | 2.535 / 3.212 | −9.44 | +0.677 | +0.498 | 0.304 |
  | −4 | 3.427 / 3.776 | 3.466 / 4.049 | −5.44 | +0.583 | +0.238 | 0.394 |
  | −2 | 3.792 / 4.140 | 3.888 / 4.296 | −3.44 | +0.408 | +0.022 | 0.489 |
  | 0 | **3.964 / 4.313**（峰） | 4.214 / 4.214 | −1.44 | +0.001 | −0.256 | 0.646 |
  | +2 | 3.792 / 4.141 | **4.293**（峰附近）/ 3.887 | +0.56 | −0.407 | −0.533 | 0.852 |
  | +4 | 3.429 / 3.778 | 4.048 / 3.465 | +2.56 | −0.583 | −0.636 | 0.944 |
  | +8 | 2.530 / 2.877 | 3.212 / 2.534 | +6.56 | −0.678 | −0.685 | 0.992 |

  这张表显示了四点：
  1. 单态 Page 曲线在 $\Delta=0$ 处取峰，两种约定同时取峰，二者之差恒为 $0.348\approx p_\tau^\infty\ln\varphi$。
  2. 码态的 $S(R)$ 在 $\Delta=+2$ 处最大（本表步长为 2）。峰位于 $\Delta\in(0,4)$，与 $\Xi=0$（$\Delta=\log_\varphi2=1.44$）相容，而不在单态的 $\Delta=0$ 处。
  3. $I_c$ 在两种约定下逐样本相同，例如 $\Delta=+2$ 时 $4.6410-4.2345=4.2933-3.8868=0.4065$。
  4. 在窗口内 $H_{\min}(Q|R)$ 与 $H(Q|R)$ 相差 $O(1)$。例如 $\Xi=-3.44$ 时二者为 0.022 与 0.408；此处 $I_c<0$，但 $F=0.49>1/k^2$。

- **相同的 von Neumann 亏损对应不同的 $F$**（$N_B=8$）。$\ln k-I_c\approx0.9$ 时，$k=2,5,13$ 的 $F$ 分别为 0.552、0.497、0.462；$\ln k-I_c\approx0.57$ 时，$k=5,13$ 分别为 0.672、0.659，$k=2$ 插值约为 0.70。因此 $F_{\rm rec}$ 不是 von Neumann Page 数据的函数，还依赖 $k$ 以及谱的完整形状。
- **两种约定下的条件量相同 [证明]。** $I_c(Q\rangle R)=S(R)-S(B)$、$I(Q\!:\!R)$、$H(Q|R)$、$H_{\min}(Q|R)$ 在两种约定下严格相等，因为 $\sum_ap_a\ln d_a$ 在 $R$ 与 $QR$（或 $R$ 与 $B$）中相同，从而抵消。

**(2) $d_\tau$、扇区权重与融合重数各在哪一步进入可操作的恢复？**
- **融合重数 $(m_a,n_a)$** 定义了信道本身，也就是块的维数。它们通过 $\bar F^{\rm Haar}_k(m_a,n_a)$ 进入（§4.1）。
- **扇区权重** $q_a\to p_a=M_a/D\to d_a^2/\mathcal D^2$ 只作为混合权重进入 $F=\sum q_aF_a$。由于 $c_a$ 的极限相同，权重在极限中消失。在有限 $N$ 下，它们通过两条途径起作用：
  1. $c_1\ne c_\tau$ 带来奇偶交替的 $O(\varphi^{-2\min})$ 修正（表中小 $N$ 的 $F_1\neq F_\tau$）；
  2. 标签泄露 $\|Q_a-q_aI\|$ 带来 $O(\sqrt{k/D})$ 修正。
- **$d_\tau$** 只通过 PF 增长 $m_a\simeq d_ad_R/\mathcal D^2$ 进入，即 $c=kd_B/d_R$，所以临界变量以 $\log_\varphi$ 为底。$\widetilde{\operatorname{Tr}}$ 中显式出现的 $d_a$ 不改变 $\mathcal N_R$、$\mathcal N_B$，也不改变 $F_{\rm rec}$。

**(3) 哪些影响只改变熵约定，哪些改变了信道？**
- **只改变熵约定：** $\widetilde{\operatorname{Tr}}$ 中的 $d_a$ 权重，即中心项 $\sum p_a\ln d_a$，它在一切条件熵中抵消；此外还有 replica 的平均方式（§2.6）。
- **改变信道：**
  - 融合规则，通过 $m_a,n_a$；
  - 超选择，即输出是块对角的。它迫使恢复按扇区进行，并把 $H(p)$ 变成双方共享的经典中心。不过它不降低最优保真度，因为最优解码本来就按扇区进行。
- **两者都不改变：** F-符号与不改变切分的基底选择（§1.3）。
- **推广（超出本题定义）：** 若总荷为 $c\ne1$，同样的约化给出与扇区无关的比值 $kd_cd_B/d_R$，阈值平移 $\log_\varphi d_c$；$c=\tau$ 时恰为一个 anyon，因为 $\operatorname{Hom}(\tau,\tau^N)\cong\operatorname{Hom}(1,\tau^{N+1})$。这是本文找到的、量子维数以加性方式真正进入恢复阈值的情形 [证明概要]。

---

## 6. 与 QES/island 比较所需的最小字典，以及缺失的部分

| 本模型中的量 | 可能对应的对象 | 状态 |
|---|---|---|
| $\ln k$（$Q$–$L$ 纠缠，码子空间熵） | island 中的 bulk entropy $S_{\rm bulk}$ | 精确定义 |
| $N_B\ln\varphi=\ln d_B$（qtr 约定；alg 约定下为 $\ln\sum n_a$ 减中心项） | 剩余黑洞的粗粒化熵，即 "area/$4G$" | 仅为类比 |
| $\sum_ap_a\ln d_a$（中心上的算符 $L=\sum_a\ln d_aP_a$，属于 $Z(\mathcal A_R)=Z(\mathcal A_B)$） | Harlow（CMP 354, 865 (2017)）意义下的中心/"面积算符"，即 edge 贡献 | 结构上精确 |
| $H(p)$ | 中心上的经典 Shannon 项 | 结构上精确 |
| $\sum p_aS(\rho_{R,a})$ | 可蒸馏的 bulk 熵 | 结构上精确 |
| 两个"鞍点"，$\min(N_R\ln\varphi,\ N_B\ln\varphi+\ln k)$ | no-island / island 两个 QES 候选 | 类比：它们对应 Weingarten 和中 $\sigma=e$ 与 $\sigma=(12)$ 两类置换的主导 |
| $\Xi=0$ | 码子空间的 Page 时刻，即纠缠楔转换 | 已证明（位置） |
| $F_{\rm rec}=e^{-H_{\min}(Q\vert R)}/k$ 与窗口函数 $f_k$ | 纠缠楔重建保真度，以及 PSSY 中 Petz 重建的平滑转换 | 可操作量，已计算 |

**缺失的部分：**
1. **没有动力学。** "时间"是人为给定的 $N_R$，没有 Hamiltonian，也没有蒸发机制。
2. **没有几何或局域性。** $V$ 是全局 Haar，bulk 没有空间结构，面积项 $\ln d_a$ 是固定常数，不是动力学面积，也没有 $G_N\to0$ 的半经典极限。
3. **Haar 平均只是模拟了鞍点求和。** 它与 replica wormhole（PSSY, arXiv:1911.11977；AEMM, JHEP 12 (2019) 063；Penington, JHEP 09 (2020) 002）的关系仍是类比。
4. **编码假设。** $L$ 是中性的；恢复只针对最大纠缠输入；没有 Hayden–Preskill 式的先验纠缠；也没有 QES 的极值化条件（Engelhardt–Wall, JHEP 01 (2015) 073）。

因此本文**不**宣称证明了 island，也**不**宣称构造了全息对偶。可以陈述的是：
1. 本模型中码子空间熵的"最小化两鞍点"结构，与可操作恢复转折的位置一致，都在 $\Xi=0$；
2. 窗口形状由 $H_{\min}$ 给出，不由 von Neumann 型的 $S_{\rm gen}$ 给出；
3. 带 $\ln d_a$ 的中心项不影响可操作恢复。

---

## 7. 已证与未闭合，以及下一步可证伪的检验

**最强的已证命题。** §4.1 的约化定理，加上 §4.2 的比值普适性，再加上 §4.3 的只依赖 $\Xi$ 的双侧界与集中性。合起来，它们在 $f_k$ 存在的前提下证明了 $\Xi$ 是正确的临界变量，且没有扇区或 $d_a$ 依赖的常数偏移。在不需要该前提的情况下，它们也证明了转折发生在 $\Xi=O(1)$。其中碰撞下界是完全证明，MP 下界与 $\lambda_{\max}$ 上界依赖所引的 RMT 标准定理。

**卡住的精确位置：**
1. **$f_k(c)$ 的存在性与闭式**（固定 $k$）。它等价于大维数下随机 SDP
$$\max_{\sigma}\big(\operatorname{Tr}\sqrt{\textstyle\sum_jA_j^\dagger\sigma A_j}\big)^2/k^2$$
   的值，其中 $A_j$ 是 $n\times m$ Gaussian 块。这需要算子值自由概率。
2. $s_k$ 的闭式。数值上 $s_2\approx1.00$。
3. §4.5 概要中的算子范数集中还需要逐项写出。

**下一步可证伪的检验：**
- (a) 若在 $\Delta$ 固定时继续增大 $N$（例如 $N_B=14$、$k=2$），$F$ 仍以超过 $10^{-3}$ 的量偏离同一 $c$ 下的普通 Haar 码，则 §4.2 的推论被否定。
- (b) 若对 $k=3,5$ 测得 $\lim_{c\to0}(1-F)/c\ne(1-k^{-2})/4$，则 §4.4 的微扰律被否定。
- (c) 最直接的反例候选是总荷为 $\tau$ 的版本，它超出本题定义。本文预测阈值恰好平移一个 anyon；若数值上看到扇区依赖的分裂，则 PF 普适性论证有误。
- (d) 若 $k=8,13,21$ 在 $\Xi$ 固定时 $F_k-f_\infty$ 不按 $k^{-2}$ 收敛到 0，则 §4.5 的猜想被否定。

---

## 8. 可复现性（全部已运行，无 NOT_RUN 代码）

| 脚本 | 内容 |
|---|---|
| `fib_basis.py` | 融合树基底、F-move、$\mathcal A_R$ 代数维数与对易核验 |
| `partA.py` | A 部分精确公式与 MC、Rényi-2 的三种平均方式 |
| `partA_asym.py` | 40 位精度下的渐近系数核验 |
| `check_gauss_replica.py` | Gaussian replica 偏差 $\psi(D+1)-\ln D$ |
| `partB_moments.py` | 精确二阶矩、暴力 Weingarten 对照、MC、泄露矩 |
| `frec.py` | 最优恢复保真度：带证书的 $\max_\sigma$ 形式，L-BFGS |
| `check_frec.py` | 与 CPTP-SDP 以及显式 Uhlmann 解码器三路互验 |
| `bounds.py` | 解耦界、碰撞界 |
| `scan_k2.py` | $k=2$ 扫描，结果在 `scan_k2.json` |
| `universality.py` | 与普通 Haar 码比较 |
| `haar_code.py` | 普通 Haar 码 |
| `scan_k.py` | $k$ 增长扫描，结果在 `scan_kgrow.json` |
| `s_k.py` | GUE-SDP 常数 $s_k$。$k=2$ 完成了 $m=8,16,32$，$m=64$ 与 $k=3$ 因内存不足未完成 |
| `entropy_vs_recovery.py` | Page 曲线、$I_c$、$H_{\min}$ 与 $F$ 的对照。最后几个 $k=13$、$\Xi\ge1.67$ 的点因内存不足中断，已完成的结果见 `entropy_vs_recovery.log` |
| `check_ic_tail.py` | $\Xi\ll0$ 端 von Neumann 信号与 $F$ 的指数对照 |

$f_\infty$ 的闭式与级数用 Mathematica 核验。

**参考文献（仅限文中实际使用处）：**
- Page, PRL 71, 1291 (1993)；Foong–Kanno, PRL 72, 1148 (1994)；Sánchez-Ruiz, PRE 52, 5653 (1995)
- Collins–Śniady, CMP 264, 773 (2006)
- König–Renner–Schaffner, IEEE TIT 55, 4337 (2009)
- Fuchs–van de Graaf, IEEE TIT 45, 1216 (1999)
- Marchenko–Pastur (1967)；Geman (1980)；Yin–Bai–Krishnaiah (1988)
- Meckes (2019)
- Araki–Yamagami (1981)
- Harlow (2017)
- Hayden–Preskill, JHEP 0709:120 (2007)
- PSSY (2019)；AEMM (2019)；Penington (2020)；Engelhardt–Wall (2015)
- $\widetilde{\operatorname{Tr}}$ 约定可参照 Bonderson–Knapp–Patel, Ann. Phys. 385, 399 (2017)。本文的推导不依赖它。
