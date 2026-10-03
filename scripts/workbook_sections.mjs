// Preserve the four-sheet workbook and extend its existing summary sections.
import fs from 'node:fs/promises';
import path from 'node:path';
import {SpreadsheetFile} from '@oai/artifact-tool';

export async function author({wb,data,root,work,scratch,out,file}){
 const {rows,full,questionData,fullFrontier,policies}=data;
 const byId=new Map(rows.map(r=>[r.id,r]));
 const n=full.length;
 const pos={base:6,compare:21,dims:36,policy:97,marginal:111,cost:133,detail:138,conclusion:172};
 data.workbookLayout=pos;
 for(let q=1;q<=3;q++){
  const sh=wb.worksheets.getItem(`Task ${q}`),rs=rows.filter(r=>r.question===q);
  // All previous source cells keep the same coordinates; append the new pair.
  const last=rs.length+1;
  for(let rr=10;rr<=last;rr++)sh.getRange(`A${rr}:D${rr}`).copyFrom(sh.getRange('A9:D9'),'all');
  sh.getRange(`A2:D${last}`).values=rs.map(r=>[r.id,r.model,r.cost,r.total]);
  sh.getRange(`A10:D${last}`).format.font={name:'Arial',size:11};
  sh.getRange(`A10:D${last}`).format.rowHeightPx=28;
  sh.getRange(`C2:C${last}`).setNumberFormat('0.000');sh.getRange(`D2:D${last}`).setNumberFormat('0');
 }
 const s=wb.worksheets.getItem('综合汇总');
 s.charts.deleteAll();s.getRange('A1:V195').unmerge();s.getRange('A1:V195').clear({applyTo:'all'});s.showGridLines=false;
 const purple='#7E0C6E',ink='#242126',gray='#66616B',light='#F4EFF4';
 for(const [c,w]of Object.entries({A:235,B:75,C:88,D:88,E:88,F:88,G:88,H:100,I:100,J:105,K:120,L:20}))s.getRange(`${c}:${c}`).format.columnWidthPx=w;
 s.getRange('A1:V190').format.font={name:'Arial',size:10,color:ink};s.getRange('A1:V190').format.rowHeightPx=28;s.getRange('A1:V190').format.verticalAlignment='center';
 s.getRange('A47:L92').format.rowHeightPx=23;
 const value=(a,v)=>s.getRange(a).values=[[v]];
 const formula=(a,v)=>s.getRange(a).formulas=[[v]];
 const ref=(id,col)=>{const r=byId.get(id);return `'${r.sheet}'!${col}${r.row}`;};
 const srcRow=id=>pos.detail+rows.findIndex(r=>r.id===id);
 function section(r,title){value(`A${r}`,title);s.getRange(`A${r}:K${r}`).format.font={size:12,bold:true,color:purple};s.getRange(`A${r}:K${r}`).format.borders={bottom:{style:'thin',color:'#C8BAC9'}};s.getRange(`A${r}:K${r}`).format.rowHeightPx=32;}
 function header(r,labels){const rg=s.getRangeByIndexes(r-1,0,1,labels.length);rg.values=[labels];rg.format={fill:purple,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:35};}
 function body(a,b,last='K'){s.getRange(`B${a}:${last}${b}`).format.horizontalAlignment='right';for(let i=a;i<=b;i++)if((i-a)%2)s.getRange(`A${i}:${last}${i}`).format.fill=light;s.getRange(`A${b}:${last}${b}`).format.borders={bottom:{style:'thin',color:'#C8BAC9'}};}
 function note(r,t){value(`A${r}`,t);s.getRange(`A${r}:K${r}`).format.font={size:10,color:gray};}
 value('A1','三道理论物理研究题：完整表现与记录成本');s.getRange('A1:K1').format.font={name:'Arial',size:17,bold:true,color:ink};s.getRange('A1:K1').format.rowHeightPx=37;
 note(2,`${data.date} · ${n} 个模型配置 × 3 题 = ${rows.length} 项采用结果 · 每题 100 分 · 美元参考成本 · 三题等权`);
 note(3,'分数沿用身份与用量揭晓前的评阅；新增 GPT-6.1 Sol xhigh / max，成本由参与者提供。');
 header(5,['模型配置','第 1 题','成本 USD','第 2 题','成本 USD','第 3 题','成本 USD','总分 /300','均分 /100','总成本 USD','最弱一题']);
 full.forEach((m,i)=>{const r=pos.base+i;value(`A${r}`,m.model);m.ids.forEach((id,j)=>{formula(`${['B','D','F'][j]}${r}`,`=${ref(id,'D')}`);formula(`${['C','E','G'][j]}${r}`,`=${ref(id,'C')}`);});formula(`H${r}`,`=SUM(B${r},D${r},F${r})`);formula(`I${r}`,`=H${r}/3`);formula(`J${r}`,`=SUM(C${r},E${r},G${r})`);formula(`K${r}`,`=MIN(B${r},D${r},F${r})`);});body(pos.base,pos.base+n-1);
 for(const c of ['C','E','G','J'])s.getRange(`${c}6:${c}${5+n}`).setNumberFormat('0.000');s.getRange(`I6:I${5+n}`).setNumberFormat('0.00');
 section(18,'三题综合比较与成本—得分前沿');
 header(20,['模型配置','排名','总分 /300','均分 /100','总成本 USD','总分 /USD','评阅下限','评阅上限','成本—得分','最弱一题','≥80 分题数']);
 const ca=pos.compare,cb=ca+n-1;
 full.forEach((m,i)=>{const r=ca+i,b=pos.base+i;value(`A${r}`,m.model);formula(`B${r}`,`=RANK(C${r},$C$${ca}:$C$${cb},0)`);formula(`C${r}`,`=H${b}`);formula(`D${r}`,`=I${b}`);formula(`E${r}`,`=J${b}`);formula(`F${r}`,`=C${r}/E${r}`);for(const [o,c]of [['G','I'],['H','J']])formula(`${o}${r}`,`=AVERAGE(${m.ids.map(id=>c+srcRow(id)).join(',')})`);formula(`I${r}`,`=IF(COUNTIFS($E$${ca}:$E$${cb},"<="&E${r},$C$${ca}:$C$${cb},">="&C${r})-COUNTIFS($E$${ca}:$E$${cb},E${r},$C$${ca}:$C$${cb},C${r})>0,"被支配","前沿")`);formula(`J${r}`,`=K${b}`);formula(`K${r}`,`=COUNTIF(B${b},">=80")+COUNTIF(D${b},">=80")+COUNTIF(F${b},">=80")`);});body(ca,cb);
 s.getRange(`D${ca}:D${cb}`).setNumberFormat('0.00');s.getRange(`E${ca}:E${cb}`).setNumberFormat('0.000');s.getRange(`F${ca}:H${cb}`).setNumberFormat('0.00');
 s.getRange(`I${ca}:I${cb}`).conditionalFormats.add('containsText',{text:'前沿',format:{fill:'#E5EEE3',font:{color:'#2D5828',bold:true}}});
 note(32,'评阅范围是三题判断范围端点的均值，非统计置信区间；前沿以中心分数判定。总分 / USD 仅为辅助指标。');
 s.getRange('A32:K32').merge();s.getRange('A32:K32').format.wrapText=true;
 section(34,'五个评分维度的三题均值');
 header(35,['模型配置','定义 /15','校准 /30','验证 /15','推进 /30','结论 /10','均分 /100','第 1 题','第 2 题','第 3 题','分差范围']);
 full.forEach((m,i)=>{const r=pos.dims+i,b=pos.base+i;value(`A${r}`,m.model);['C','D','E','F','G'].forEach((c,j)=>formula(`${['B','C','D','E','F'][j]}${r}`,`=AVERAGE(${m.ids.map(id=>c+srcRow(id)).join(',')})`));formula(`G${r}`,`=SUM(B${r}:F${r})`);['B','D','F'].forEach((c,j)=>formula(`${['H','I','J'][j]}${r}`,`=${c}${b}`));formula(`K${r}`,`=MAX(H${r}:J${r})-MIN(H${r}:J${r})`);});body(pos.dims,pos.dims+n-1);s.getRange(`B${pos.dims}:G${pos.dims+n-1}`).setNumberFormat('0.00');
 const chartBindings=[];
 function scatter(title,pts,top,a,b){
  s.getRange(`N${top}:P${top}`).values=[['Model','log10(USD)','Score /100']];
  const ch=s.charts.add('scatter');ch.title=title;ch.titleTextStyle.fontSize=12;ch.titleTextStyle.typeface='Arial';ch.hasLegend=true;ch.legend={position:'bottom',textStyle:{typeface:'Arial',fontSize:8}};
  ch.xAxis={title:{text:'log10(USD)',textStyle:{fontSize:10}},numberFormatCode:'0.0',numberFormatSourceLinked:false,textStyle:{fontSize:9}};
  ch.yAxis={title:{text:'Score / 100',textStyle:{fontSize:10}},numberFormatCode:'0',numberFormatSourceLinked:false,textStyle:{fontSize:9},majorGridlines:{fill:'#E4E1E6',style:'solid',width:.5}};
  pts.forEach((p,i)=>{const r=top+i+1;value(`N${r}`,p.model);formula(`O${r}`,`=LOG10(${p.costFormula})`);formula(`P${r}`,`=${p.scoreFormula}`);const sr=ch.series.add(p.short);sr.xValues=[Math.log10(p.cost)];sr.values=[p.score];sr.fill=p.color;sr.formula=`'综合汇总'!$P$${r}`;});ch.setPosition(a,b);
  chartBindings.push({title,mode:'xValues numeric snapshot; y cell formula',points:pts.map((p,i)=>({model:p.model,x:Math.log10(p.cost),y:p.score,sourceX:`O${top+i+1}`,sourceY:`P${top+i+1}`}))});
 }
 scatter('三题：均分与总成本',full.map((m,i)=>({...m,score:m.mean,costFormula:`J${pos.base+i}`,scoreFormula:`I${pos.base+i}`})),5,'A48','F69');
 for(let q=1;q<=3;q++)scatter(`第 ${q} 题：得分与记录成本`,questionData[q-1].map(r=>({...r,score:r.total,costFormula:ref(r.id,'C'),scoreFormula:ref(r.id,'D')})),18+(q-1)*12,...[['G48','L69'],['A71','F92'],['G71','L92']][q-1]);
 await fs.writeFile(path.join(work,'chart_bindings.json'),JSON.stringify(chartBindings,null,2));
 section(94,`事后按题选型：${data.portfolioCount} 种组合中的代表方案`);
 header(96,['选择方案','总分 /300','均分 /100','总成本 USD','题 1 分','题 2 分','题 3 分','题 1 USD','题 2 USD','题 3 USD','回答编号']);
 policies.forEach((p,i)=>{const r=pos.policy+i;value(`A${r}`,p.name);p.ids.forEach((id,j)=>{formula(`${['E','F','G'][j]}${r}`,`=${ref(id,'D')}`);formula(`${['H','I','J'][j]}${r}`,`=${ref(id,'C')}`);});formula(`B${r}`,`=SUM(E${r}:G${r})`);formula(`C${r}`,`=B${r}/3`);formula(`D${r}`,`=SUM(H${r}:J${r})`);value(`K${r}`,p.ids.join(' / '));});body(pos.policy,pos.policy+policies.length-1);s.getRange('A97:K104').format.rowHeightPx=44;s.getRange('A97:A104').format.wrapText=true;s.getRange('K97:K104').format.wrapText=true;s.getRange('C97:C104').setNumberFormat('0.00');s.getRange('D97:D104').setNumberFormat('0.000');s.getRange('H97:J104').setNumberFormat('0.000');
 note(106,'组合依据本次已知得分作事后选择，不能外推为未来新题的预期效果；完整模型名称见分析报告。');
 section(108,'统一模型配置的边际成本');
 header(110,['比较（全部三题）','总分增量','均分增量','美元增量','成本变化率','USD /分','前者总成本','后者总成本']);
 const changes=[['6.1 max 相对 6.1 xhigh','GPT-6.1 Sol max','GPT-6.1 Sol xhigh'],['6.1 xhigh 相对 6 Sol max','GPT-6.1 Sol xhigh','GPT-6 Sol max'],['6.1 max 相对 6 Sol max','GPT-6.1 Sol max','GPT-6 Sol max'],['6.1 xhigh 相对 Astra high','GPT-6.1 Sol xhigh','GPT-6 Astra high'],['6.1 xhigh 相对 Luna max','GPT-6.1 Sol xhigh','GPT-6 Luna max']];
 changes.forEach(([label,an,bn],i)=>{const r=pos.marginal+i,a=ca+full.findIndex(m=>m.model===an),b=ca+full.findIndex(m=>m.model===bn);value(`A${r}`,label);formula(`B${r}`,`=C${a}-C${b}`);formula(`C${r}`,`=D${a}-D${b}`);formula(`D${r}`,`=E${a}-E${b}`);formula(`E${r}`,`=D${r}/E${b}`);formula(`F${r}`,`=D${r}/B${r}`);formula(`G${r}`,`=E${a}`);formula(`H${r}`,`=E${b}`);});body(111,115,'H');s.getRange('A111:A115').format.wrapText=true;s.getRange('A111:H115').format.rowHeightPx=38;s.getRange('C111:C115').setNumberFormat('0.00');s.getRange('D111:D115').setNumberFormat('0.000');s.getRange('E111:E115').setNumberFormat('0.0%');s.getRange('F111:H115').setNumberFormat('0.000');
 section(118,'统计口径与数据来源');
 [
  '评分依据三批已完成评阅：定义 15、校准 30、验证 15、推进 30、结论校准 10；本次只合并身份和用量。',
  '新增映射：s2_1a / s2_2a / s2_3a = GPT-6.1 Sol xhigh；s2_1b / s2_2b / s2_3b = GPT-6.1 Sol max。',
  'xhigh 三题美元消耗：1.612、1.111、1.371；max：1.795、1.140、1.641。来源：用户 2026-09-30 披露。',
  '每个配置均有三题，均分等权。新增六项成本合计 8.670 USD；未提供 token 数，不倒推 token 或 FLOPs。',
  'Opus 三题采用补充报告评分，成本沿用原记录且各计一次；原始匿名答卷与评阅均保留。',
  '每个模型—题目只有一份采用结果，交互轮次和工具预算未证明严格一致；结论限于本次测试材料。',
  '评阅区间非统计置信区间；6.1 max 与 xhigh 的小分差不足以证明稳定的总体能力差异。',
  '评分按匿名编号开展；部分新增正文费用附录出现型号文字，未据此调分；同族隐性偏好仍可能存在。',
  '参考成本沿用参与者记录，不等于订阅实付、推理时间或算力；总分 / USD 不具有科研价值倍数含义。',
  '原生图横轴为 log10(USD) 且采用本次成本快照；编辑费用后应重新生成。PNG/SVG 图含评阅范围。',
  '来源：评分细则、三批 q1/q2/q3 评阅报告、第二轮补充评分明细，以及参与者确认的费用映射。',
  '项目与来源：https://github.com/Zhen-WushuiLingchun/math_physics_ai_benchmark_v1',
 ].forEach((t,i)=>note(120+i,t));
 value('A133',`本表 ${rows.length} 项采用结果的记录成本 USD`);formula('E133',"=SUM('Task 1'!C2:C11,'Task 2'!C2:C11,'Task 3'!C2:C11)");s.getRange('E133').setNumberFormat('0.000');
 section(135,'逐份回答：五维评分、评阅范围与溯源');
 header(137,['模型配置','题目','定义 /15','校准 /30','验证 /15','推进 /30','结论 /10','总分 /100','评阅下限','评阅上限','采用编号']);
 rows.forEach((r,i)=>{const rr=pos.detail+i;s.getRange(`A${rr}:G${rr}`).values=[[r.model,r.question,...['definition','calibration','verification','research','calibration_of_claims'].map(k=>r.scores[k])]];formula(`H${rr}`,`=${ref(r.id,'D')}`);s.getRange(`I${rr}:K${rr}`).values=[[r.interval[0],r.interval[1],r.id]];});body(pos.detail,pos.detail+rows.length-1);
 section(170,'本次完整样本的模型评价');
 const judgments={
  'GPT-6.1 Sol max':'均分最高；核空间证书、恢复界及最佳波形幂次成立，验证覆盖较完整。',
  'GPT-6.1 Sol xhigh':'高分与成本折中突出；比 max 少 3 个总分，费用少 0.482 USD。',
  'GPT-6 Astra high':'三题证据强；在本次中心分数与成本上被两个 6.1 配置支配。',
  'GPT-6 Sol max':'三题均衡且 ≥94；新 6.1 配置以更低成本取得更高分数。',
  'Opus 5.5 xhigh':'有价值的低点与解析归约；一般化论证和较高成本限制本次竞争力。',
  'GPT-6 Sol xhigh':'题 2、3 较强；题 1 的关键输入错误进入核心证书，跨题差异明显。',
  'GPT-5.6 Sol xhigh':'部分计算正确；目标错配和低频缺口限制结论，记录成本较高。',
  'GPT-6 Luna max':'最低成本的完整配置；题 2 基础与局部恢复结果较强，研究层仍有缺口。',
  'DeepSeek V4.1 Flash max':'部分结构成立，但核心输入与统一性推论不足；本次被 Luna 逐题支配。',
  'Grok 4.7 high':'有局部正确计算；树级相位与恢复基线等校准问题影响主要结论。',
 };
 full.forEach((m,i)=>note(pos.conclusion+i,`${m.model}：${m.score}/300，成本 ${m.cost.toFixed(3)} USD。${judgments[m.model]}`));
 for(const r of rows){const got=wb.worksheets.getItem(r.sheet).getRange(`A${r.row}:D${r.row}`).values[0];if(got.some((v,i)=>v!==[r.id,r.model,r.cost,r.total][i]))throw new Error(`Source mismatch ${r.id}`);}
 full.forEach((m,i)=>{const got=s.getRange(`H${pos.base+i}:K${pos.base+i}`).values[0];[m.score,m.mean,m.cost,m.minimum].forEach((v,j)=>{if(Math.abs(got[j]-v)>1e-8)throw new Error(`Summary mismatch ${m.model}`);});if((s.getRange(`I${ca+i}`).values[0][0]==='前沿')!==fullFrontier.some(x=>x.model===m.model))throw new Error(`Frontier mismatch ${m.model}`);});
 if(Math.abs(s.getRange('E133').values[0][0]-data.totalCost)>1e-8)throw new Error('Total cost mismatch');
 const check=await wb.inspect({kind:'table',range:'综合汇总!A5:K15',include:'values,formulas',tableMaxRows:11,tableMaxCols:11,maxChars:2500});await fs.writeFile(path.join(scratch,'export_inspect.ndjson'),check.ndjson);
 const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},summary:'final formula error scan'});await fs.writeFile(path.join(scratch,'formula_scan.ndjson'),errors.ndjson);console.log(errors.ndjson);
 for(const [name,sheetName,range,scale]of [
  ['summary_tables','综合汇总','A1:K45',1.15],['summary_charts','综合汇总','A48:L92',1.2],['summary_notes','综合汇总','A94:K133',1.15],['score_components','综合汇总','A135:K167',1.2],['summary_conclusions','综合汇总','A170:K181',1.2],
  ['after_q1','Task 1','A1:D11',1.5],['after_q2','Task 2','A1:D11',1.5],['after_q3','Task 3','A1:D11',1.5]]){
  const b=await wb.render({sheetName,range,scale});await fs.writeFile(path.join(scratch,`${name}.png`),new Uint8Array(await b.arrayBuffer()));
 }
 await (await SpreadsheetFile.exportXlsx(wb)).save(file);
 await fs.writeFile(path.join(work,'analysis_data.json'),JSON.stringify(data,null,2));
 console.log(JSON.stringify({output:file,answers:rows.length,models:full.length,charts:4,totalCost:data.totalCost}));
}
