"""Historical blind-batch assembler, not the adopted dataset builder.
Reinitializes identity/cost fields in the supplementary batch; do not use to
regenerate the published identity-integrated dataset.
"""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DIMENSIONS = {
    "definition": 15,
    "calibration": 30,
    "verification": 15,
    "research": 30,
    "calibration_of_claims": 10,
}
Q3_KEYS = {
    "definition_and_setup": "definition",
    "calibration": "calibration",
    "independent_verification": "verification",
    "research_progress": "research",
    "conclusion_calibration": "calibration_of_claims",
}
EXPECTED = {1: {"s_1a", "s_1b", "s_1c", "s_1d"},
            2: {"s_2a", "s_2b"}, 3: {"s_3a", "s_3b"}}
questions = []
for number, records_key in ((1, "answers"), (2, "records"), (3, "assessments")):
    data = json.loads((OUT / f"q{number}_scores.json").read_text(encoding="utf-8-sig"))
    answers = []
    for source in data[records_key]:
        scores = source["scores"]
        if number == 3:
            scores = {Q3_KEYS[key]: value for key, value in scores.items()}
        assert set(scores) == set(DIMENSIONS)
        assert all(isinstance(v, int) and 0 <= v <= DIMENSIONS[k] for k, v in scores.items())
        assert sum(scores.values()) == source["total"]
        assert source["interval"][0] <= source["total"] <= source["interval"][1]
        pdf = ROOT / "submissions" / "supplement" / f"{number:02d}" / f"{source['id']}.pdf"
        assert pdf.is_file()
        answer = {
            "id": source["id"],
            "question": number,
            "source_pdf": pdf.as_posix(),
            "model_name": None,
            "reasoning_effort": None,
            "cost_usd": None,
            "scores": scores,
            "total": source["total"],
            "interval": source["interval"],
            "summary": source["summary"],
            "strongest": source["strongest"],
            "errors": source["errors"],
            "next_step": source["next_step"],
            "source_detail_json": (OUT / f"q{number}_scores.json").as_posix(),
        }
        answers.append(answer)
    assert {a["id"] for a in answers} == EXPECTED[number]
    answers.sort(key=lambda a: (-a["total"], a["id"]))
    for answer in answers:
        answer["rank_within_supplement_question"] = 1 + sum(
            other["total"] > answer["total"] for other in answers)
    questions.append({
        "question": number,
        "assessment_report": (OUT / f"q{number}_assessment.md").as_posix(),
        "answers": answers,
    })

blinding_note = (
    "辅助评阅意外接触s_3a的PDF元数据中潜在身份线索；主评未接触该线索，"
    "并在查看拟评分前独立复核其主论证。未据此认定身份或调整分数。"
)
combined = {
    "date": "2026-09-23",
    "batch": "Supplement",
    "stage": "anonymous_scoring_complete_identity_and_cost_integration_pending",
    "rubric": (ROOT / "questions" / "math_physics_ai_benchmark_rubric.md").as_posix(),
    "dimensions": DIMENSIONS,
    "intervals_are_statistical_confidence_intervals": False,
    "cross_question_letter_mapping": False,
    "model_mapping_status": "pending_user_disclosure",
    "new_scores_finalized_before_model_disclosure": True,
    "global_workbook_updated": False,
    "blinding_note": blinding_note,
    "total_answers": 8,
    "questions": questions,
}
json_path = OUT / "补充评分明细.json"
json_path.write_text(json.dumps(combined, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

def link(label, path):
    return f"[{label}](<{path.as_posix()}>)"

lines = [
    "# Supplement 补充回答评分汇总",
    "",
    "评阅日期：2026-09-23。本轮 8 份补充回答已评分完成：第 1 题 4 份，第 2、3 题各 2 份。"
    "第 2、3 题由两名 Astra max 子智能体独立评阅，主评完成第 1 题并核对汇总与关键结论。"
    "沿用此前同一份评分细则、五维权重和匿名旧答卷的校准尺度。",
    "",
    "分数在模型身份公布前确定。不同题的相同字母不表示同一模型。当前只汇总补充评分，"
    "模型对应关系、实际成本及与原工作簿的合并留待身份公布后处理。"
    "区间是合理复评范围，不是统计置信区间；评分反映本次答卷证据，不能直接视作普遍能力排名。",
    "",
    "| 题目 | 回答 | 定义 /15 | 校准 /30 | 验证 /15 | 推进 /30 | 结论 /10 | 总分 | 评阅区间 |",
    "|---|---|---:|---:|---:|---:|---:|---:|---|",
]
for question in questions:
    for answer in question["answers"]:
        dims = " | ".join(str(answer["scores"][key]) for key in DIMENSIONS)
        lo, hi = answer["interval"]
        lines.append(f"| 第 {question['question']} 题 | {answer['id']} | {dims} | "
                     f"**{answer['total']}** | {lo}–{hi} |")
lines += [
    "",
    "本轮八份回答的总分互不相同。与此前同题回答比较时，s_2a 与旧 2g 均为 93 分，"
    "s_3a 与旧 3g 均为 81 分：它们的主要成立结果和核心缺口相近，因此保留同分。"
    "没有为了形成唯一名次而调整缺乏证据支持的分差。",
    "",
    "## 第 1 题：YM–EYM 一圈共线重构",
    "",
    "**s_1b：94 分。** 六个排序分别满足的树级—单负振幅关系，经精确符号复算全部成立；"
    "EYM 目标违反这一必要条件，因此其排除覆盖题面整个允许核类。树零空间、四个缺陷和 "
    "D 维积分后抵消均通过核验。它准确保留了跨 helicity 最小修复及高点推广的缺口。",
    "",
    "**s_1a：87 分。** 在明确有理运动学点上保留符号 x 的左零向量证书成立，"
    "三个 EYM 输入和四个缺陷全部复现。软方向留数归约与圈动量插入结构有实质价值；"
    "高精度有理重构缺少严格认证，退化排序的软展开及 n=4 数值秩仍需补证。",
    "",
    "**s_1c：53 分。** 树级校准、全正缺陷和维数移位积分分析保留得分。"
    "核心单负输入存在相对符号错误；它的整数证书确实证明了错误输入下的 12/13 秩，"
    "换成正确振幅后同一个全正加两个单负位置的系统变为 9/9，不能据原证书推出 no-go。",
    "",
    "**s_1d：47 分。** 树零空间和部分全正结果成立，但单负错号进入主要证书，"
    "负腿 2 的缺陷表又与自身公式不符。其显式全正修复有符号错误；"
    "从一列增广秩差为 1 推断最小物理扩张维数为 1，也没有成立。",
    "",
    link("第 1 题完整评阅、页码与精确证据", OUT / "q1_assessment.md"),
    "",
    "## 第 2 题：Fibonacci 任意子随机编码与恢复",
    "",
    "**s_2a：93 分。** 融合空间、熵与 Stiefel 矩完整，"
    "扇区极分解归约、碰撞/MP 下界和谱上界确实约束最优恢复。"
    "主要限制是集中性尚不能自动给出唯一均值极限，精确转折位置及增长 k 的最优函数仍未完成。",
    "",
    "**s_2b：75 分。** 核心 Stiefel 矩、Petz 恢复与最优 CPTP 的 SDP 公式正确，"
    "约 0.20 的经验常数有部分数据支持。"
    "平均约化态熵、HS 阈值和独立随机列数值对照有明确错误；"
    "只检查计算基上的标签分布不能排除量子叠加态的标签泄漏。窗口最优恢复的解析约束不足。",
    "",
    link("第 2 题完整评阅、反例与数值对照", OUT / "q2_assessment.md"),
    "",
    "## 第 3 题：QNM 谱迁移与观测波形稳定性",
    "",
    "**s_3a：81 分。** 精确 Jost 匹配、源—探测器表达式和显式能量界成立，"
    "在对数时间窗得到 O(ε log²(1/ε)) 的控制；多条数值比较也有实质证据。"
    "但全位置、全时间 α=1 的结论仍未证明，"
    "固定位置非零首阶只能排除 α>1，低频分母的 L⁻⁴ 主张还与静态变分矛盾。",
    "",
    "**s_3b：54 分。** 正确的因果时刻、固定位置首阶共振公式与辅助 delta 模型保留得分。"
    "mismatch 展开漏了投影项，临界参数处外推及全局 resolvent 对象有误；"
    "主定理所需的一致流量估计没有闭合，也未另给可靠的统一有限时间替代界。"
    "同参数数值表约 2% 的冲突降低验证强度，但不据此认定计算虚构。",
    "",
    link("第 3 题完整评阅、静态核验与两网格检查", OUT / "q3_assessment.md"),
    "",
    "## 后续合并所需信息与数据状态",
    "",
    "原有 19 份回答加上本轮 8 份，共有 27 份已评回答，三个题目各 9 份。"
    "这只是答卷数量；是否构成同一组模型的完整三题配对，需要用户提供模型身份对应关系后确认。"
    "本轮未改动原评分工作簿，也未因成本或猜测身份回调任何分数。",
    "",
    "收到 s_1a–s_1d、s_2a–s_2b、s_3a–s_3b 与实际模型及推理档位的对应关系后，"
    "再并入已有成本记录，更新总表与综合 sheet，重新输出实名逐题分析、逐模型评价、"
    "完整三题成本—得分曲线、Pareto 前沿与边际成本分析。成本是本测试中的实际美元支出；"
    "若没有独立算力测量，不把价格直接等同于 FLOPs 或计算时间。",
    "",
    link("可供后续合并的补充评分 JSON", json_path),
    "",
    "**评阅方法说明。** " + blinding_note + "所有评分均依据正文中的数学物理内容与可复核证据。",
    "",
]
report_path = OUT / "补充评分汇总.md"
report_path.write_text("\n".join(lines), encoding="utf-8")

reloaded = json.loads(json_path.read_text(encoding="utf-8"))
assert sum(len(q["answers"]) for q in reloaded["questions"]) == 8
assert {a["id"] for q in reloaded["questions"] for a in q["answers"]} == set.union(*EXPECTED.values())
assert all(Path(q["assessment_report"]).is_file() for q in reloaded["questions"])
print(json.dumps({
    "validated": True,
    "new_answers": 8,
    "question_counts": {q["question"]: len(q["answers"]) for q in reloaded["questions"]},
    "totals": {a["id"]: a["total"] for q in reloaded["questions"] for a in q["answers"]},
    "report": report_path.as_posix(),
    "json": json_path.as_posix(),
}, ensure_ascii=False, indent=2))
