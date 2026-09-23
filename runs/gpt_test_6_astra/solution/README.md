# 题一：精确低点解答

正文：`q1_solution.tex`；编译文件：`../output/pdf/q1_solution.pdf`。

**结论**：以正文明确列出的一圈 EYM/YM 标准振幅为输入，n=3 的六个非相邻排序上，全部允许的局域一次 Mandelstam 核都无法同时重构全正和三个单负一圈扇区。显式函数域证书位于正文定理 1。树级零核空间维数为 8。

不是仅测试某个树核。也不把 n=3 结论声称为固定 n=4 的独立判定。最小扩张、全多重数闭合尚未完成。

## 复算

在原题目录运行：

```powershell
python solution/checks/verify_q1.py
& ./solution/build.ps1
```

依赖：Python、SymPy、XeLaTeX（含 ctex/Fandol 字体）。正式检查含 558 项精确断言和 12 组额外有理运动学点；没有浮点秩判断。`build/verification.json` 记录运行状态，`build/certificate.json` 给出矩阵、行标签、证书和样本。

## 来源和限制

`sources/manifest.json` 记录四篇原始文献的下载地址与 SHA-256。PDF 为研究来源，独立于解答存放。已核对 NPT 的纯胶子圈、最大 g 耦合阶和所用有限扇区；没有把标量计算代理当作额外物质。

证书是本次推导；EYM 的物理输入来自所引文献，其积分极限在本次脚本中复核。这不是另一套独立费曼图推导。

`checks/explore.py`、`certificate_explore.py`、`loop_certificate_explore.py`、`fast_certificate.py` 为探索记录。一个较慢的通用零空间求解被主动停止，随后正式脚本用 QQ(x) 上的算法与显式左零向量完成验证；以 `verify_q1.py` 的输出为正式证据。
