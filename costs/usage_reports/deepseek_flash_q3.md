> 公开副本：已省略会话标识和本机路径；原始会话日志未公开。

# Token 用量与美元计价报告（题三研究任务）

> 统计对象：本次对话（黑洞 QNM 谱不稳定性研究 + XeLaTeX 解答生成）
（含本机路径的行已省略）
> 模型：`deepseek/deepseek-flash`（DeepSeek V4.1 Flash，variant=max，provider=`deepseek`）
（含本机路径的行已省略）
（含本机路径的行已省略）
> 统计截止：2026-09-23 18:59:20（本回合收尾消息未计入，量级 < $0.005）

---

## 1. 核心数字（Token 总量）

| 类别 | Token 数 | 说明 |
|---|---:|---|
| input（未命中缓存） | 236,302 | 每轮新增的未缓存 prompt |
| **cache_read（缓存命中）** | **65,283,968** | 每轮重读的历史上下文、工具输出 |
| cache_write | 0 | 该 provider 无单独缓存写价 |
| output（可见输出） | 82,269 | 正文/代码/工具调用参数 |
| reasoning（思维链） | 326,671 | 按输出价计费 |
| **总计（计费口径）** | **65,929,210** | input + output + reasoning + cache_read |

- 助理回合数：189（另有 2 条 user 消息；总计 191 条消息）
- 会话时长：2026-09-23 17:45:50 → 18:59:20，约 **73.5 分钟**
- 单回合平均：input 1,250 · cache_read 345,418 · output 435 · reasoning 1,728
- 单回合平均 prompt 总量（input+cache_read）≈ **346,668 token**

## 2. 缓存 / 输出 比例

| 比例 | 数值 |
|---|---:|
| **缓存命中率** = cache_read / (cache_read + input) | **99.639 %** |
| 输出占比（output / 全部计费 token） | 0.125 % |
| 输出+推理占比（(output+reasoning) / 全部计费 token） | 0.620 % |
| (output+reasoning) / prompt 总量（input+cache_read） | 0.624 % |
| reasoning / output | **3.97** |

> 解读：这是一个"长上下文、多工具回合"型任务——prompt 侧 99.6% 由缓存命中承担；
> 真正的"生成"（output+reasoning）只占全部计费 token 的 0.62%。
> 若你的关注点是"模型实际写了多少字"，看 output（82k）+ reasoning（327k）。

## 3. 美元计价（USD）

opencode 对 `deepseek-flash` 的定价（models.dev 缓存，USD / 百万 token）：

| 项目 | 单价（USD / 1M tokens） |
|---|---:|
| input（cache miss） | 0.15 |
| output | 0.60 |
| reasoning | 0.60 |
| cache_read | 0.003 |

**费用分解（本次会话）：**

| 项目 | Token | 费用（USD） | 占比 |
|---|---:|---:|---:|
| input | 236,302 | $0.035445 | 7.4 % |
| output | 82,269 | $0.049361 | 10.4 % |
| reasoning | 326,671 | $0.196003 | 41.1 % |
| cache_read | 65,283,968 | $0.195852 | 41.1 % |
| **合计** | 65,929,210 | **$0.476661** | 100 % |

- 与 opencode 数据库逐条记录的 `cost` 求和**完全一致（差值 = 0.00）**，
  且与逐条消息用上式重算一致，可作为可信计价。
- 有效混合单价：$0.476661 / 65.93M ≈ **$0.00723 / 1M tokens**（含全部类别；
  因 99% 的计费 token 是单价仅 $0.003/M 的缓存读，混合价被显著拉低）。
- 折合单回合：**≈ $0.00252 / 回合**。
- 若**没有缓存**（64.9M cache_read 按 input 价 $0.15/M 计费）：约 **$10.07**，
  缓存节省 **$9.60（95.3 %）**。
- 关注"纯生成成本"：output+reasoning 合计 $0.245（占 51.5 %）；
  可理解为把全部开销摊到可见输出上 ≈ $5.79 / 1M output token。

## 4. 相关会话（对照参考，同一 ds_test 测试家族）

| 会话 | 目录 | 记录成本（USD） | input | output | reasoning | cache_read |
|---|---|---:|---:|---:|---:|---:|
| `ses_f50bb25f…vrDr` | ds_test | 0.5353 | 835,390 | 76,698 | 383,726 | 44,569,216 |
| `ses_f324dbcb…6NBr` | ds_test2 | 0.4350 | 252,666 | 79,474 | 294,155 | 57,624,320 |
| **本次 `ses_f3257c41…wh6a`** | **ds_test3** | **0.4767** | **236,302** | **82,269** | **326,671** | **65,283,968** |

> 注：表中数字为各会话在 opencode 会话表中的记录值（快照），
> 不同会话可能混用了不同 provider/模型（价格不同），仅作规模对照。

## 5. 复现方法

统计脚本：`qnm_study/db_report.py`（只读方式打开 opencode.db，不会影响正在运行的会话）。
会话结束后可重跑一次获得终值：

```powershell
（含本机路径的行已省略）
python qnm_study\db_report.py
```

计价公式（可直接复核）：

```
cost_USD = (input*0.15 + (output+reasoning)*0.60 + cache_read*0.003 + cache_write*0.15) / 1e6
```

## 6. 注意事项

1. **截止时点**：上面的数字取自 18:59:20 前已写入数据库的消息；本回合最后一条回复会再增加
   约 0.3–0.5M cache_read 与 1–3k output/reasoning（< $0.005），对比例与总价影响可忽略。
2. **reasoning 单独计费**：opencode 中 reasoning 与 output 同价（$0.60/M），
   本次 reasoning 是 output 的 3.97 倍，是最主要的"生成侧"成本项（41.1 %）。
3. **cache_read 是最大"流量"项**：65.3M（占全部计费 token 的 99.0 %），
   但因单价低（1/50 input 价），只占费用 41.1 %；这正是长会话依赖 prompt 缓存的原因。
4. 本任务在**同一对话内**完成（未派生 subagent 会话），故以上即该任务的全部用量；
   若把同一测试家族的三个会话合并计算，总记录成本约 **$1.45**。
5. 若需要按人民币计价：按汇率自行折算即可（例如 1 USD ≈ 7.1 CNY 时，本次约 ¥3.38）。
