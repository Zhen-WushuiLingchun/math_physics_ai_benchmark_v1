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
if(mode==='inspect'){
 console.log((await wb.inspect({kind:'workbook,sheet,table',maxChars:3000,tableMaxRows:3,tableMaxCols:4})).ndjson);
 console.log((await wb.inspect({kind:'region',sheetId:'综合汇总',range:'A5:K13',maxChars:3000,tableMaxRows:10,tableMaxCols:11})).ndjson);
 for(const [name,sheetName,range]of [['before_summary','综合汇总','A1:K39'],['before_task','Task 1','A1:D9'],['before_charts','综合汇总','A42:L86']]){
  const b=await wb.render({sheetName,range,scale:1.2});await fs.writeFile(path.join(scratch,`${name}.png`),new Uint8Array(await b.arrayBuffer()));
 }
 process.exit(0);
}
const original=JSON.parse(await fs.readFile(path.join(root,'evaluations/original/评分明细.json'),'utf8'));
const supplement=JSON.parse(await fs.readFile(path.join(root,'evaluations/supplement/补充评分明细.json'),'utf8'));
const second=JSON.parse(await fs.readFile(path.join(root,'evaluations/supplement_2nd/第二轮补充评分明细.json'),'utf8'));
const secondCosts=JSON.parse(await fs.readFile(path.join(root,'costs/supplement_2nd_costs.json'),'utf8'));
const grades=new Map([...original.questions,...supplement.questions].flatMap(q=>q.answers.map(a=>({...a,question:q.question}))).concat(second.answers).map(a=>[a.id,a]));
const originalRows=JSON.parse(await fs.readFile(path.join(work,'input_rows.json'),'utf8'));
const oldCosts=new Map(Object.values(originalRows).flatMap(rs=>rs.slice(1)).map(r=>[r[0],r[2]]));
const newCosts={'s_1b':4.118,'s_1c':.178,'s_1d':.414,'s_2b':.430,'s_3b':.477};
for(const r of secondCosts.entries)newCosts[r.id]=r.cost_usd;
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
 {name:'GPT-6.1 Sol xhigh',short:'6.1 Sol xhigh',ids:['s2_1a','s2_2a','s2_3a'],color:'#15746B'},
 {name:'GPT-6.1 Sol max',short:'6.1 Sol max',ids:['s2_1b','s2_2b','s2_3b'],color:'#C43E66'},
];
const taskIds=[['1a','1b','1c','1d','s_1a','s_1b','s_1c','s_1d'],['2a','2b','2c','2d','2e','2f','s_2a','s_2b'],['3a','3b','3c','3d','3e','3f','s_3a','s_3b']];
taskIds.forEach((ids,q)=>ids.push(`s2_${q+1}a`,`s2_${q+1}b`));
const round=x=>Math.round(x*1000)/1000;
const rows=taskIds.flatMap((ids,q)=>ids.map((id,i)=>{
 const g=grades.get(id),d=definitions.find(d=>d.ids.includes(id));
 const cost=newCosts[id]??oldCosts.get(replacements[id]??id);
 if(!g||!d||!(cost>0)||Object.values(g.scores).reduce((a,b)=>a+b,0)!==g.total)throw new Error(`Bad source: ${id}`);
 const batch=id.startsWith('s2_')?'supplement_2nd':id.startsWith('s_')?'supplement':'original';
 return {...g,question:q+1,model:d.name,model_name:d.name,short:d.short,color:d.color,cost,cost_usd:cost,source_pdf:`submissions/${batch}/0${q+1}/${id}.pdf`,sheet:`Task ${q+1}`,row:i+2,replaces:replacements[id]??null,costSource:id.startsWith('s2_')?'用户披露：2026-09-30，USD':newCosts[id]?'用户补充用量':'原表用量'};
}));
if(rows.length!==30||new Set(rows.map(r=>r.id)).size!==30)throw new Error('Need 30 adopted answers');
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
const policies=[makePolicy('三题均用 Luna max',same('GPT-6 Luna max')),makePolicy('每题至少 75 分：最低成本',cheapestAt(75)),makePolicy('每题至少 80 分：最低成本',cheapestAt(80)),makePolicy('每题至少 90 分：最低成本',cheapestAt(90)),makePolicy('每题至少 95 分：最低成本',cheapestAt(95)),makePolicy('三题均用 6.1 Sol xhigh',same('GPT-6.1 Sol xhigh')),makePolicy('三题均用 6.1 Sol max',same('GPT-6.1 Sol max')),makePolicy('逐题最高分：最低成本',portfolios.toSorted((a,b)=>b.score-a.score||a.cost-b.cost)[0])];
const fullFrontier=frontier(full);
const data={date:'2026-09-30',scope:'10 configurations x 3 questions; scores fixed before user identity/cost mapping; Opus old cost reused once',definitions,rows,full,questionData,fullFrontier,questionFrontiers:questionData.map(rs=>frontier(rs.map(r=>({...r,score:r.total})))),portfolioFrontier:frontier(portfolios),portfolioCount:portfolios.length,policies,totalCost:round(rows.reduce((a,r)=>a+r.cost,0)),dimensions:original.dimensions};
await fs.writeFile(path.join(work,'analysis_data.json'),JSON.stringify(data,null,2));
if(mode==='analyze'){
 console.log(JSON.stringify({full:full.map(m=>({model:m.model,score:m.score,cost:m.cost,min:m.minimum,range:[m.lower,m.upper],dims:m.dims})),frontier:fullFrontier.map(m=>m.model),policies,portfolioFrontier:data.portfolioFrontier,totalCost:data.totalCost},null,2));process.exit(0);
}
if(mode!=='author')throw new Error('Use analyze or author');
await import('./workbook_sections.mjs').then(m=>m.author({wb,data,root,work,scratch,out,file}));
process.exit(0);
