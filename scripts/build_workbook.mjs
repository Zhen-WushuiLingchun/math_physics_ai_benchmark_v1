import fs from 'node:fs/promises';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import {FileBlob, SpreadsheetFile} from '@oai/artifact-tool';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const work=path.join(root,'data');
const scratch=path.join(root,'.work');
await fs.mkdir(scratch,{recursive:true});
const out=path.join(root,'results');
const file=path.join(out,'模型评分与成本汇总.xlsx');
const mode=process.argv[2]||'analyze';
const wb=await SpreadsheetFile.importXlsx(await FileBlob.load(file));
const original=JSON.parse(await fs.readFile(path.join(root,'evaluations/original/评分明细.json'),'utf8'));
const supplement=JSON.parse(await fs.readFile(path.join(root,'evaluations/supplement/补充评分明细.json'),'utf8'));
const grades=new Map([...original.questions,...supplement.questions].flatMap(q=>q.answers.map(a=>({...a,question:q.question}))).map(a=>[a.id,a]));
const originalRows=JSON.parse(await fs.readFile(path.join(work,'input_rows.json'),'utf8'));
const oldCosts=new Map(Object.values(originalRows).flatMap(rs=>rs.slice(1)).map(r=>[r[0],r[2]]));
const newCosts={'s_1b':4.118,'s_1c':.178,'s_1d':.414,'s_2b':.430,'s_3b':.477};
const replacements={'s_1a':'1e','s_2a':'2g','s_3a':'3g'};
const definitions=[
 {name:'GPT-6 Astra high',short:'Astra high',ids:['1b','2e','3e'],color:'#7E0C6E'},
 {name:'GPT-6 Sol max',short:'Sol max',ids:['s_1b','2c','3c'],color:'#167E9A'},
 {name:'GPT-6 Sol xhigh',short:'Sol xhigh',ids:['1c','2d','3d'],color:'#478B42'},
 {name:'GPT-6 Luna max',short:'Luna max',ids:['s_1c','2b','3b'],color:'#C88700'},
 {name:'GPT-5.6 Sol xhigh',short:'5.6 Sol xhigh',ids:['1a','2a','3a'],color:'#52677E'},
 {name:'Opus 5.5 xhigh',short:'Opus 5.5 xhigh',ids:['s_1a','s_2a','s_3a'],color:'#BE644B'},
 {name:'Grok 4.7 high',short:'Grok 4.7 high',ids:['1d','2f','3f'],color:'#888888'},
 {name:'DeepSeek V4.1 Flash max',short:'DeepSeek Flash max',ids:['s_1d','s_2b','s_3b'],color:'#6857B1'},
];
const taskIds=[['1a','1b','1c','1d','s_1a','s_1b','s_1c','s_1d'],['2a','2b','2c','2d','2e','2f','s_2a','s_2b'],['3a','3b','3c','3d','3e','3f','s_3a','s_3b']];
const round=x=>Math.round(x*1000)/1000;
const rows=taskIds.flatMap((ids,q)=>ids.map((id,i)=>{
 const g=grades.get(id),d=definitions.find(d=>d.ids.includes(id));
 const cost=newCosts[id]??oldCosts.get(replacements[id]??id);
 if(!g||!d||!(cost>0)||Object.values(g.scores).reduce((a,b)=>a+b,0)!==g.total)throw new Error(`Bad source: ${id}`);
 return {...g,model:d.name,short:d.short,color:d.color,cost,sheet:`Task ${q+1}`,row:i+2,replaces:replacements[id]??null,costSource:newCosts[id]?'用户补充用量':'原表用量'};
}));
if(rows.length!==24||new Set(rows.map(r=>r.id)).size!==24)throw new Error('Need 24 adopted answers');
const byId=new Map(rows.map(r=>[r.id,r]));
const frontier=pts=>pts.filter(a=>!pts.some(b=>b.cost<=a.cost&&b.score>=a.score&&(b.cost<a.cost||b.score>a.score))).sort((a,b)=>a.cost-b.cost);
const full=definitions.map(d=>{
 const rs=d.ids.map(id=>byId.get(id)),score=rs.reduce((a,r)=>a+r.total,0),cost=round(rs.reduce((a,r)=>a+r.cost,0));
 return {model:d.name,short:d.short,color:d.color,ids:d.ids,n:3,score,mean:score/3,cost,meanCost:cost/3,pointsPerDollar:score/cost,minimum:Math.min(...rs.map(r=>r.total)),at80:rs.filter(r=>r.total>=80).length,lower:rs.reduce((a,r)=>a+r.interval[0],0)/3,upper:rs.reduce((a,r)=>a+r.interval[1],0)/3,dims:Object.fromEntries(Object.keys(rs[0].scores).map(k=>[k,rs.reduce((a,r)=>a+r.scores[k],0)/3]))};
}).sort((a,b)=>b.score-a.score);
const questionData=[1,2,3].map(q=>rows.filter(r=>r.question===q).sort((a,b)=>b.total-a.total));
const portfolios=[];
for(const a of questionData[0])for(const b of questionData[1])for(const c of questionData[2])portfolios.push({ids:[a.id,b.id,c.id],names:[a.model,b.model,c.model],score:a.total+b.total+c.total,mean:(a.total+b.total+c.total)/3,cost:round(a.cost+b.cost+c.cost),minimum:Math.min(a.total,b.total,c.total)});
const cheapestAt=t=>portfolios.filter(p=>p.minimum>=t).sort((a,b)=>a.cost-b.cost||b.score-a.score)[0];
const makePolicy=(name,p)=>({name,...p});
const same=name=>{const m=full.find(d=>d.model===name);return portfolios.find(p=>p.ids.join()===m.ids.join());};
const policies=[makePolicy('三题均用 Luna max',same('GPT-6 Luna max')),makePolicy('每题至少 75 分：最低成本',cheapestAt(75)),makePolicy('每题至少 80 分：最低成本',cheapestAt(80)),makePolicy('每题至少 90 分：最低成本',cheapestAt(90)),makePolicy('三题均用 Sol max',same('GPT-6 Sol max')),makePolicy('逐题最高分组合',portfolios.toSorted((a,b)=>b.score-a.score||a.cost-b.cost)[0])];
const fullFrontier=frontier(full);
const data={date:'2026-09-23',scope:'8 configurations x 3 questions; Opus reformatted answers adopted; old cost reused once',definitions,rows,full,questionData,fullFrontier,questionFrontiers:questionData.map(rs=>frontier(rs.map(r=>({...r,score:r.total})))),portfolioFrontier:frontier(portfolios),portfolioCount:portfolios.length,policies,totalCost:round(rows.reduce((a,r)=>a+r.cost,0)),dimensions:original.dimensions};
await fs.writeFile(path.join(work,'analysis_data.json'),JSON.stringify(data,null,2));
if(mode==='analyze'){
 console.log(JSON.stringify({full:full.map(m=>({model:m.model,score:m.score,cost:m.cost,min:m.minimum,range:[m.lower,m.upper],dims:m.dims})),frontier:fullFrontier.map(m=>m.model),policies,portfolioFrontier:data.portfolioFrontier,totalCost:data.totalCost},null,2));process.exit(0);
}
if(mode!=='author')throw new Error('Use analyze or author');
// Update the three existing source sheets in place.
for(let q=1;q<=3;q++){
 const sh=wb.worksheets.getItem(`Task ${q}`),rs=rows.filter(r=>r.question===q);
 sh.getRange('A2:D9').values=rs.map(r=>[r.id,r.model,r.cost,r.total]);
 sh.getRange('A2:D9').format.font={name:'Arial',size:11};
 sh.getRange('A:A').format.columnWidthPx=75;sh.getRange('B:B').format.columnWidthPx=255;
 sh.getRange('C:D').format.columnWidthPx=125;sh.getRange('A1:D9').format.rowHeightPx=28;
 sh.getRange('C2:C9').setNumberFormat('0.000');sh.getRange('D2:D9').setNumberFormat('0');
}
const s=wb.worksheets.getItem('综合汇总');
s.charts.deleteAll();s.getRange('A1:V170').clear({applyTo:'all'});s.showGridLines=false;
const purple='#7E0C6E',ink='#242126',gray='#66616B',light='#F4EFF4';
for(const [c,w]of Object.entries({A:230,B:72,C:88,D:88,E:88,F:88,G:88,H:100,I:100,J:100,K:118,L:20}))s.getRange(`${c}:${c}`).format.columnWidthPx=w;
s.getRange('A1:V170').format.font={name:'Arial',size:10,color:ink};s.getRange('A1:V170').format.rowHeightPx=28;s.getRange('A1:V170').format.verticalAlignment='center';
s.getRange('A41:L87').format.rowHeightPx=23;
const value=(a,v)=>s.getRange(a).values=[[v]];
const formula=(a,v)=>s.getRange(a).formulas=[[v]];
const ref=(id,col)=>{const r=byId.get(id);return `'${r.sheet}'!${col}${r.row}`;};
const srcRow=id=>125+rows.findIndex(r=>r.id===id);
function section(r,title){value(`A${r}`,title);s.getRange(`A${r}:K${r}`).format.font={size:12,bold:true,color:purple};s.getRange(`A${r}:K${r}`).format.borders={bottom:{style:'thin',color:'#C8BAC9'}};s.getRange(`A${r}:K${r}`).format.rowHeightPx=32;}
function header(r,labels){const rg=s.getRangeByIndexes(r-1,0,1,labels.length);rg.values=[labels];rg.format={fill:purple,font:{name:'Arial',size:10,bold:true,color:'#FFFFFF'},wrapText:true,horizontalAlignment:'center',verticalAlignment:'center',rowHeight:35};}
function body(a,b,last='K'){s.getRange(`B${a}:${last}${b}`).format.horizontalAlignment='right';for(let i=a;i<=b;i++)if((i-a)%2)s.getRange(`A${i}:${last}${i}`).format.fill=light;s.getRange(`A${b}:${last}${b}`).format.borders={bottom:{style:'thin',color:'#C8BAC9'}};}
function note(r,t){value(`A${r}`,t);s.getRange(`A${r}:K${r}`).format.font={size:10,color:gray};}
value('A1','三道数学物理研究题：完整表现与记录成本');s.getRange('A1:K1').format.font={name:'Arial',size:17,bold:true,color:ink};s.getRange('A1:K1').format.rowHeightPx=37;
note(2,'2026-09-23 · 8 个模型配置 × 3 题 = 24 项有效结果 · 每题 100 分 · 美元参考成本 · 三题等权');
note(3,'Opus 5.5 xhigh 采用补充报告的评分；三题成本沿用原记录，各计一次。身份揭晓后沿用已完成的评分。');
header(5,['模型配置','第 1 题','成本 USD','第 2 题','成本 USD','第 3 题','成本 USD','总分 /300','均分 /100','总成本 USD','最弱一题']);
full.forEach((m,i)=>{const n=i+6;value(`A${n}`,m.model);m.ids.forEach((id,j)=>{formula(`${['B','D','F'][j]}${n}`,`=${ref(id,'D')}`);formula(`${['C','E','G'][j]}${n}`,`=${ref(id,'C')}`);});formula(`H${n}`,`=SUM(B${n},D${n},F${n})`);formula(`I${n}`,`=H${n}/3`);formula(`J${n}`,`=SUM(C${n},E${n},G${n})`);formula(`K${n}`,`=MIN(B${n},D${n},F${n})`);});body(6,13);
for(const c of ['C','E','G','J'])s.getRange(`${c}6:${c}13`).setNumberFormat('0.000');s.getRange('I6:I13').setNumberFormat('0.00');
section(16,'三题综合比较与成本—得分前沿');
header(18,['模型配置','排名','总分 /300','均分 /100','总成本 USD','总分 /USD','评阅下限','评阅上限','成本—得分','最弱一题','≥80 分题数']);
full.forEach((m,i)=>{const n=19+i,b=6+i;value(`A${n}`,m.model);formula(`B${n}`,`=RANK(C${n},$C$19:$C$26,0)`);formula(`C${n}`,`=H${b}`);formula(`D${n}`,`=I${b}`);formula(`E${n}`,`=J${b}`);formula(`F${n}`,`=C${n}/E${n}`);for(const [o,c]of [['G','I'],['H','J']])formula(`${o}${n}`,`=AVERAGE(${m.ids.map(id=>c+srcRow(id)).join(',')})`);formula(`I${n}`,`=IF(COUNTIFS($E$19:$E$26,"<="&E${n},$C$19:$C$26,">="&C${n})-COUNTIFS($E$19:$E$26,E${n},$C$19:$C$26,C${n})>0,"被支配","前沿")`);formula(`J${n}`,`=K${b}`);formula(`K${n}`,`=COUNTIF(B${b},">=80")+COUNTIF(D${b},">=80")+COUNTIF(F${b},">=80")`);});body(19,26);
s.getRange('D19:D26').setNumberFormat('0.00');s.getRange('E19:E26').setNumberFormat('0.000');s.getRange('F19:H26').setNumberFormat('0.00');
s.getRange('I19:I26').conditionalFormats.add('containsText',{text:'前沿',format:{fill:'#E5EEE3',font:{color:'#2D5828',bold:true}}});
note(28,'评阅上下限是三题判断范围端点的均值，不是统计置信区间；前沿以中心分数判定。总分 / USD 仅为辅助指标。');
section(30,'五个评分维度的三题均值');
header(31,['模型配置','定义 /15','校准 /30','验证 /15','推进 /30','结论 /10','均分 /100','第 1 题','第 2 题','第 3 题','分差范围']);
full.forEach((m,i)=>{const n=32+i,b=6+i;value(`A${n}`,m.model);['C','D','E','F','G'].forEach((c,j)=>formula(`${['B','C','D','E','F'][j]}${n}`,`=AVERAGE(${m.ids.map(id=>c+srcRow(id)).join(',')})`));formula(`G${n}`,`=SUM(B${n}:F${n})`);['B','D','F'].forEach((c,j)=>formula(`${['H','I','J'][j]}${n}`,`=${c}${b}`));formula(`K${n}`,`=MAX(H${n}:J${n})-MIN(H${n}:J${n})`);});body(32,39);s.getRange('B32:G39').setNumberFormat('0.00');
// The documented scatter xValues export numeric snapshots; y values remain cell-bound.
// Adjacent formula-backed helpers preserve an auditable cost and score source.
const chartBindings=[];
function scatter(title,pts,top,a,b){
 s.getRange(`N${top}:P${top}`).values=[['Model','log10(USD)','Score /100']];
 const series=pts.map((p,i)=>{const r=top+i+1;value(`N${r}`,p.model);formula(`O${r}`,`=LOG10(${p.costFormula})`);formula(`P${r}`,`=${p.scoreFormula}`);return {name:p.short,xValues:[Math.log10(p.cost)],values:[p.score],fill:p.color};});
 const ch=s.charts.add('scatter');
 ch.title=title;ch.titleTextStyle.fontSize=12;ch.titleTextStyle.typeface='Arial';
 ch.hasLegend=true;ch.legend={position:'bottom',textStyle:{typeface:'Arial',fontSize:8}};
 ch.xAxis={title:{text:'log10(USD)',textStyle:{fontSize:10}},numberFormatCode:'0.0',numberFormatSourceLinked:false,textStyle:{fontSize:9}};
 ch.yAxis={title:{text:'Score / 100',textStyle:{fontSize:10}},numberFormatCode:'0',numberFormatSourceLinked:false,textStyle:{fontSize:9},majorGridlines:{fill:'#E4E1E6',style:'solid',width:.5}};
 for(const p of series){const sr=ch.series.add(p.name);sr.xValues=p.xValues;sr.values=p.values;sr.fill=p.fill;}
 ch.setPosition(a,b);
 ch.series.items.forEach((sr,i)=>{sr.formula=`'综合汇总'!$P$${top+i+1}`;});
 if(ch.series.items.length!==pts.length)throw new Error('Missing chart series');
 chartBindings.push({title,mode:'xValues numeric snapshot; y cell formula',points:pts.map((p,i)=>({model:p.model,x:Math.log10(p.cost),y:p.score,sourceX:`O${top+i+1}`,sourceY:`P${top+i+1}`}))});
}
scatter('三题：均分与总成本',full.map((m,i)=>({...m,score:m.mean,costFormula:`J${i+6}`,scoreFormula:`I${i+6}`})),5,'A42','F63');
for(let q=1;q<=3;q++)scatter(`第 ${q} 题：得分与记录成本`,questionData[q-1].map(r=>({...r,score:r.total,costFormula:ref(r.id,'C'),scoreFormula:ref(r.id,'D')})),18+(q-1)*12,...[['G42','L63'],['A65','F86'],['G65','L86']][q-1]);
await fs.writeFile(path.join(work,'chart_bindings.json'),JSON.stringify(chartBindings,null,2));
section(88,'事后按题选型：512 种组合中的代表方案');
header(90,['选择方案','总分 /300','均分 /100','总成本 USD','题 1 分','题 2 分','题 3 分','题 1 USD','题 2 USD','题 3 USD','回答编号']);
policies.forEach((p,i)=>{const n=91+i;value(`A${n}`,p.name);p.ids.forEach((id,j)=>{formula(`${['E','F','G'][j]}${n}`,`=${ref(id,'D')}`);formula(`${['H','I','J'][j]}${n}`,`=${ref(id,'C')}`);});formula(`B${n}`,`=SUM(E${n}:G${n})`);formula(`C${n}`,`=B${n}/3`);formula(`D${n}`,`=SUM(H${n}:J${n})`);value(`K${n}`,p.ids.join(' / '));});body(91,96);s.getRange('A91:K96').format.rowHeightPx=42;s.getRange('A91:A96').format.wrapText=true;s.getRange('K91:K96').format.wrapText=true;s.getRange('C91:C96').setNumberFormat('0.00');s.getRange('D91:D96').setNumberFormat('0.000');s.getRange('H91:J96').setNumberFormat('0.000');
note(98,'组合利用本次已知得分选题，不代表未来新题可预先取得同样效果；完整选型名称见分析报告。');
section(100,'统一模型配置的边际成本');
header(102,['比较（全部三题）','总分增量','均分增量','美元增量','成本变化率','USD /分','前者总成本','后者总成本']);
const changes=[['Sol xhigh 相对 Luna max','GPT-6 Sol xhigh','GPT-6 Luna max'],['Sol max 相对 Sol xhigh','GPT-6 Sol max','GPT-6 Sol xhigh'],['Astra high 相对 Sol max','GPT-6 Astra high','GPT-6 Sol max']];
changes.forEach(([label,an,bn],i)=>{const n=103+i,a=19+full.findIndex(m=>m.model===an),b=19+full.findIndex(m=>m.model===bn);value(`A${n}`,label);formula(`B${n}`,`=C${a}-C${b}`);formula(`C${n}`,`=D${a}-D${b}`);formula(`D${n}`,`=E${a}-E${b}`);formula(`E${n}`,`=D${n}/E${b}`);formula(`F${n}`,`=D${n}/B${n}`);formula(`G${n}`,`=E${a}`);formula(`H${n}`,`=E${b}`);});body(103,105,'H');s.getRange('C103:C105').setNumberFormat('0.00');s.getRange('D103:D105').setNumberFormat('0.000');s.getRange('E103:E105').setNumberFormat('0.0%');s.getRange('F103:H105').setNumberFormat('0.000');
section(108,'统计口径与数据来源');
[
 '评分沿用原始及 Supplement 已完成评阅：定义 15、校准 30、验证 15、推进 30、结论校准 10。',
 '身份与新增成本来自用户本轮映射；原有费用保留。三份 Opus 补充报告替换其对应回答，成本不重复计入。',
 '每个配置均有三题，均分采用等权；模型排名、均分与合计成本由 Task 1–3 的单元格公式计算。',
 '参考：评分明细.json、补充评分明细.json，以及两批 q1/q2/q3_assessment.md；报告保留逐题依据。',
 '每个模型—题目仅一个样本。Opus 采用格式统一后的报告；这里比较最终提交版本，并非完全统一的首轮协议实验。',
 '成本是用户记录的美元参考成本；无 token、FLOPs、耗时及重复运行数据，不能推断算力缩放律或总体胜率。',
 '成本—得分前沿：不存在成本不高、得分不低且至少一项更优的另一配置；小分差应结合评阅范围理解。',
 '图中横轴为 log10(USD)：−1=$0.1，0=$1，1=$10；原生图横坐标为本次成本快照，改动成本后需重绘。',
 '误差棒与曲线的完整版本见随附 PNG/SVG；图上标注靠近数据点。前沿连线仅帮助阅读，不作连续性能拟合。',
].forEach((t,i)=>note(110+i,t));
value('A120','本表 24 项有效结果的记录成本 USD');formula('E120',"=SUM('Task 1'!C2:C9,'Task 2'!C2:C9,'Task 3'!C2:C9)");s.getRange('E120').setNumberFormat('0.000');
section(122,'逐份回答：五维评分、评阅范围与溯源');
header(124,['模型配置','题目','定义 /15','校准 /30','验证 /15','推进 /30','结论 /10','总分 /100','评阅下限','评阅上限','采用编号']);
rows.forEach((r,i)=>{const n=125+i;s.getRange(`A${n}:G${n}`).values=[[r.model,r.question,r.scores.definition,r.scores.calibration,r.scores.verification,r.scores.research,r.scores.calibration_of_claims]];formula(`H${n}`,`=${ref(r.id,'D')}`);s.getRange(`I${n}:K${n}`).values=[[r.interval[0],r.interval[1],r.id]];});body(125,148);
section(151,'本次完整样本的模型评价');
[
 'GPT-6 Astra high：总分 287、均分 95.67；三题均强，临界窗反向界与全时间统一稳定性界最突出。',
 'GPT-6 Sol max：总分 283、均分 94.33；题 1 最高，三题均 ≥94，成本比 Astra 低 42.6%。',
 'Opus 5.5 xhigh：总分 261、均分 87.00；有价值的解析归约，部分一般化与最优指数主张缺少充分证明。',
 'GPT-6 Sol xhigh：总分 236、均分 78.67；题 2/3 较强，题 1 的关键输入错号进入不可解性证书。',
 'GPT-5.6 Sol xhigh：总分 226、均分 75.33；基础计算与归约有价值，题 1 目标错配和题 3 低频缺口明显。',
 'GPT-6 Luna max：总分 215、均分 71.67，三题合计 $0.464；低成本优势明显，题 1 核心证书失效。',
 'DeepSeek V4.1 Flash max：总分 176、均分 58.67；部分结构正确，关键输入、推论与统一性验证不足。',
 'Grok 4.7 high：总分 169、均分 56.33；有局部正确结果，树级相位和恢复基线等校准错误限制结论。',
].forEach((t,i)=>note(153+i,t));
// Verify source cells and calculated summaries before export.
for(const r of rows){const got=wb.worksheets.getItem(r.sheet).getRange(`A${r.row}:D${r.row}`).values[0];const exp=[r.id,r.model,r.cost,r.total];if(got.some((v,i)=>v!==exp[i]))throw new Error(`Source mismatch ${r.id}`);}
full.forEach((m,i)=>{const got=s.getRange(`H${i+6}:K${i+6}`).values[0];[m.score,m.mean,m.cost,m.minimum].forEach((v,j)=>{if(Math.abs(got[j]-v)>1e-8)throw new Error(`Summary mismatch ${m.model} ${j}: ${got[j]}`);});const ff=s.getRange(`I${19+i}`).values[0][0];if((ff==='前沿')!==fullFrontier.some(x=>x.model===m.model))throw new Error(`Frontier mismatch ${m.model}`);});
if(Math.abs(s.getRange('E120').values[0][0]-data.totalCost)>1e-8)throw new Error('Total cost mismatch');
const errors=await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A|#NUM!|#NULL!|#SPILL!|#CALC!',options:{useRegex:true,maxResults:100},summary:'final formula error scan'});
await fs.writeFile(path.join(scratch,'formula_scan.ndjson'),errors.ndjson);console.log(errors.ndjson);
for(const [name,sheetName,range,scale]of [
 ['summary_tables','综合汇总','A1:K39',1.25],['summary_charts','综合汇总','A42:L86',1.15],['summary_notes','综合汇总','A88:K120',1.25],['score_components','综合汇总','A122:K148',1.2],['summary_conclusions','综合汇总','A151:K160',1.3],
 ['after_q1','Task 1','A1:D9',1.5],['after_q2','Task 2','A1:D9',1.5],['after_q3','Task 3','A1:D9',1.5]]){
 const p=await wb.render({sheetName,range,scale});await fs.writeFile(path.join(scratch,`${name}.png`),new Uint8Array(await p.arrayBuffer()));
}
await (await SpreadsheetFile.exportXlsx(wb)).save(file);
supplement.stage='identity_cost_integrated';supplement.model_mapping_status='user_confirmed';supplement.global_workbook_updated=true;
for(const q of supplement.questions)for(const a of q.answers){const r=byId.get(a.id);a.model_name=r.model;a.reasoning_effort=r.model.endsWith('xhigh')?'xhigh':'max';a.cost_usd=r.cost;a.cost_source=r.costSource;a.replaces_answer_id=r.replaces;}
await fs.writeFile(path.join(root,'evaluations/supplement/补充评分明细.json'),JSON.stringify(supplement,null,2));
console.log(JSON.stringify({output:file,answers:24,models:8,charts:4,totalCost:data.totalCost}));
