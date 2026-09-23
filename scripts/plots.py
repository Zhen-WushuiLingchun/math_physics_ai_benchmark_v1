from pathlib import Path
import sys,json,math,zipfile
import xml.etree.ElementTree as ET
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'data';OUT=ROOT/'results'
FIG=OUT/'figures'
FIG.mkdir(parents=True,exist_ok=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,FixedLocator
from matplotlib import font_manager
import os
font_path=Path(os.environ.get('BENCHMARK_FONT',r'C:\Windows\Fonts\msyh.ttc'))
if font_path.is_file():font_manager.fontManager.addfont(str(font_path))
font_candidates=['Microsoft YaHei','Noto Sans CJK SC','Source Han Sans SC','SimHei','WenQuanYi Zen Hei']
available={f.name for f in font_manager.fontManager.ttflist}
font_family=next((f for f in font_candidates if f in available),'DejaVu Sans')
plt.rcParams.update({'font.family':font_family,'font.size':10,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#B8B3BC','axes.labelcolor':'#37313D','text.color':'#29232F','xtick.color':'#57505E','ytick.color':'#57505E','svg.fonttype':'none','pdf.fonttype':42})
data=json.loads((WORK/'analysis_data.json').read_text(encoding='utf-8'))

def save(fig,name,dpi=180):
 fig.savefig(FIG/(name+'.png'),dpi=dpi,bbox_inches='tight',facecolor='white')
 fig.savefig(FIG/(name+'.svg'),bbox_inches='tight',facecolor='white')
 svg=FIG/(name+'.svg')
 svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 plt.close(fig)

def cost_axis(ax,lo=.1,hi=65):
 ax.set_xscale('log');ax.set_xlim(lo,hi)
 ax.xaxis.set_major_locator(FixedLocator([x for x in [.1,.3,1,3,10,30,60] if lo<=x<=hi]))
 ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'${x:g}'))
 ax.tick_params(axis='x',which='minor',bottom=False)
 ax.grid(axis='y',color='#E8E4EB',lw=.7,zorder=0)
 ax.set_xlabel('参考成本（USD，对数坐标）')

def point(ax,r,mean=False):
 y=r['mean'] if mean else r['total'];lo=r['lower'] if mean else r['interval'][0];hi=r['upper'] if mean else r['interval'][1]
 ax.errorbar(r['cost'],y,yerr=[[y-lo],[hi-y]],fmt='o',color=r['color'],ecolor=r['color'],ms=6.5,capsize=3,lw=1.1,alpha=.95,zorder=3)

def label(ax,r,offset,ha,va,mean=False,size=8.2,short=False):
 y=r['mean'] if mean else r['total'];name=r['short'] if short else r['model']
 score=f'{y:.2f}' if mean else str(y)
 ax.annotate(f"{name}\n{score} 分 · ${r['cost']:.3f}",(r['cost'],y),xytext=offset,textcoords='offset points',fontsize=size,color=r['color'],ha=ha,va=va,linespacing=1.25,bbox={'facecolor':'white','edgecolor':'none','alpha':.93,'pad':1.0},zorder=4)

def conclusion_figures():
 import shutil
 positions={
  'GPT-6 Astra high':((9,1),'left','center'),
  'GPT-6 Sol max':((-9,-3),'right','top'),
  'Opus 5.5 xhigh':((-9,1),'right','center'),
  'GPT-6 Sol xhigh':((-9,1),'right','center'),
  'GPT-5.6 Sol xhigh':((9,1),'left','center'),
  'GPT-6 Luna max':((9,1),'left','center'),
  'DeepSeek V4.1 Flash max':((9,1),'left','center'),
  'Grok 4.7 high':((9,-2),'left','center'),
 }
 for scale,scale_name in [('log','对数坐标'),('linear','线性坐标')]:
  fig,ax=plt.subplots(figsize=(13.4,8.6))
  fs=data['fullFrontier']
  ax.plot([r['cost'] for r in fs],[r['mean'] for r in fs],color='#AAA0AF',lw=1.2,ls='--',zorder=1)
  for r in data['full']:
   point(ax,r,True)
   off,ha,va=positions[r['model']]
   if scale=='linear' and r['model']=='GPT-6 Sol xhigh':off,ha,va=(9,1),'left','center'
   label(ax,r,off,ha,va,True,size=8.2)
  if scale=='log':cost_axis(ax,.3,85)
  else:
   ax.set_xscale('linear');ax.set_xlim(0,56)
   ax.xaxis.set_major_locator(FixedLocator([0,10,20,30,40,50]))
   ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'${x:g}'))
   ax.grid(axis='y',color='#E8E4EB',lw=.7,zorder=0)
  ax.set_ylim(46,103);ax.set_yticks([50,60,70,80,90,100])
  ax.set_xlabel(f'完成全部三题的合计参考成本（USD，{scale_name}）',labelpad=10)
  ax.set_ylabel('三题等权均分 / 100',labelpad=10)
  ax.set_title('理论物理问题测试：模型得分与成本',loc='left',fontsize=18,pad=30)
  ax.text(0,1.025,'8 个模型配置 × 3 道理论物理问题 · 三题等权评分',transform=ax.transAxes,fontsize=10.5,color='#625A68')
  fig.text(.075,.172,'虚线连接中心分数下的 Pareto 前沿；竖线为三题评阅范围端点的均值，非统计置信区间。',fontsize=8.5,color='#625A68')
  fig.text(.075,.147,'成本采用本次记录的美元参考值，不代表实付金额、FLOPs 或推理时间；连线仅作视觉引导。',fontsize=8.5,color='#625A68')
  fig.add_artist(plt.Line2D([.075,.97],[.130,.130],transform=fig.transFigure,color='#DDD7E0',linewidth=.6))
  credits=[
   '题库与评分标准：真·无水零醇使用 GPT 6 Pro 辅助完成。',
   '问题作答：真·无水零醇和食司使用相关模型共同完成。',
   '评分：食司使用 GPT 6 Astra Max 完成。',
   '评分采用匿名方式，向评分模型隐去作答模型的型号；但不排除同一族模型在出题和评分时存在隐性偏好。',
  ]
  for y,t in zip([.105,.079,.053,.027],credits):fig.text(.075,y,t,fontsize=8.2,color='#655E6B')
  fig.subplots_adjust(left=.075,right=.97,top=.87,bottom=.265)
  save(fig,'结果'+('对数' if scale=='log' else '线性')+'版',dpi=240)

if '--conclusion-only' in sys.argv:
 conclusion_figures()
 print('Updated final logarithmic and linear cost-score figures (PNG/SVG).')
 sys.exit(0)

# All 24 individual results; common y scale and direct labels beside points.
fig,axs=plt.subplots(1,3,figsize=(19.2,6.8),sharey=True)
titles=['第 1 题：YM → EYM 一圈提升','第 2 题：Fibonacci 信息恢复','第 3 题：振铃稳定性']
positions={
 's_1b':((-7,4),'right','bottom'),'1b':((7,-3),'left','top'),'s_1a':((-7,-12),'right','top'),
 '1a':((7,0),'left','center'),'1c':((-7,-6),'right','top'),'s_1c':((7,2),'left','bottom'),
 's_1d':((7,-6),'left','top'),'1d':((7,0),'left','center'),
 '2e':((7,5),'left','bottom'),'2c':((-7,2),'right','bottom'),'s_2a':((7,-3),'left','top'),
 '2d':((7,-4),'left','top'),'2a':((-7,-7),'right','top'),'2b':((7,-4),'left','top'),
 's_2b':((7,-4),'left','top'),'2f':((7,-4),'left','top'),
 '3e':((7,3),'left','bottom'),'3c':((-7,3),'right','bottom'),'3d':((-7,-3),'right','top'),
 '3a':((-7,4),'right','bottom'),'s_3a':((7,-3),'left','top'),'3b':((7,-4),'left','top'),
 '3f':((7,-4),'left','top'),'s_3b':((7,-4),'left','top'),
}
for q,ax in enumerate(axs):
 fs=data['questionFrontiers'][q]
 ax.plot([r['cost'] for r in fs],[r['total'] for r in fs],color='#B3ABB8',lw=1.1,ls='--',zorder=1)
 for r in data['questionData'][q]:
  point(ax,r);off,ha,va=positions[r['id']];label(ax,r,off,ha,va,size=7.3,short=True)
 cost_axis(ax,.1,90);ax.set_ylim(20,108);ax.set_yticks([20,40,60,80,100]);ax.set_title(titles[q],loc='left',fontsize=12,pad=16)
axs[0].set_ylabel('本题得分 / 100')
fig.suptitle('逐题比较：8 个模型配置，各有完整三题结果',x=.055,ha='left',fontsize=17,y=.975)
fig.text(.055,.063,'标签中的 Astra、Sol、Luna 均为 GPT-6；5.6 Sol 为 GPT-5.6；DeepSeek Flash 为 DeepSeek V4.1 Flash。',fontsize=9,color='#625A68')
fig.text(.055,.026,'竖线为评阅判断范围，非统计置信区间；虚线连接本题成本—得分前沿，仅作视觉引导。每个点对应一份采用的回答。',fontsize=9,color='#625A68')
fig.subplots_adjust(left=.055,right=.99,top=.865,bottom=.18,wspace=.13)
save(fig,'逐题成本与得分')

# All three questions: final conclusion figures with attribution.
conclusion_figures()

fig,ax=plt.subplots(figsize=(12.4,7.0))
ps=data['portfolioFrontier'];xs=[r['cost'] for r in ps];ys=[r['mean'] for r in ps]
ax.step(xs+[15],ys+[ys[-1]],where='post',color='#7E0C6E',lw=1.8,zorder=2)
ax.scatter(xs,ys,c='#7E0C6E',s=24,zorder=3)
annotations={
 .464:('三题 Luna max\n$0.464 · 71.67 分',(8,6),'left','bottom'),
 4.404:('Sol max / Luna / Luna\n$4.404 · 85.33 分',(-8,9),'right','bottom'),
 5.904:('Sol max / Luna / Sol xhigh\n$5.904 · 90.00 分',(-8,15),'right','bottom'),
 7.575:('Sol max / Sol max / Sol xhigh\n$7.575 · 92.67 分',(8,-12),'left','top'),
 8.252:('三题 Sol max\n$8.252 · 94.33 分',(-8,16),'right','bottom'),
 12.874:('Sol max / Astra / Astra\n$12.874 · 96.33 分',(-8,13),'right','bottom'),
}
for r in ps:
 if r['cost'] in annotations:
  t,off,ha,va=annotations[r['cost']];ax.annotate(t,(r['cost'],r['mean']),xytext=off,textcoords='offset points',fontsize=8.0,ha=ha,va=va,linespacing=1.25,color='#4C354F',bbox={'facecolor':'white','edgecolor':'none','alpha':.95,'pad':1.1},zorder=4)
ax.set_xlim(0,15);ax.set_ylim(68,104);ax.set_xlabel('完成三题的参考成本预算上限（USD）');ax.set_ylabel('预算内观测到的最高三题均分 / 100');ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'${x:g}'));ax.grid(axis='y',color='#E8E4EB',lw=.7)
ax.set_title('按题选择模型：预算—得分阶梯',loc='left',fontsize=17,pad=29)
ax.text(0,1.025,'枚举 8 × 8 × 8 = 512 种组合；标签按第 1 / 2 / 3 题排列。',transform=ax.transAxes,fontsize=10,color='#625A68')
fig.text(.075,.05,'此图利用已知答案得分作事后选择，不能作为未来新题的预期表现；预算小于 $0.464 时无法覆盖全部三题。',fontsize=9,color='#625A68')
fig.subplots_adjust(left=.075,right=.98,top=.84,bottom=.17)
save(fig,'三题预算与得分前沿')

# Read-only validation of the finished Open XML export.
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
CN={'c':'http://schemas.openxmlformats.org/drawingml/2006/chart'}
with zipfile.ZipFile(OUT/'模型评分与成本汇总.xlsx') as z:
 shared=[]
 if 'xl/sharedStrings.xml' in z.namelist():shared=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
 wb=ET.fromstring(z.read('xl/workbook.xml'));names=[e.attrib['name'] for e in wb.find('s:sheets',NS)]
 assert names==['Task 1','Task 2','Task 3','综合汇总'],names
 cells=[];errors=[]
 for i in range(4):
  root=ET.fromstring(z.read(f'xl/worksheets/sheet{i+1}.xml'));cs={}
  for c in root.findall('.//s:sheetData/s:row/s:c',NS):
   v=c.find('s:v',NS);f=c.find('s:f',NS);t=c.attrib.get('t')
   if t=='s':value=shared[int(v.text)]
   elif t=='inlineStr':value=''.join(c.find('s:is',NS).itertext())
   elif v is None:value=None
   elif t in ['str','e']:value=v.text
   else:value=float(v.text) if v.text else None
   cs[c.attrib['r']]={'value':value,'formula':f.text if f is not None else None}
   if t=='e':errors.append((names[i],c.attrib['r'],value))
  cells.append(cs)
 for r in data['rows']:
  sh=cells[r['question']-1]
  for col,exp in zip('ABCD',[r['id'],r['model'],r['cost'],r['total']]):
   got=sh[f"{col}{r['row']}"]['value'];assert got==exp,(r['id'],col,got,exp)
 for i,m in enumerate(data['full']):
  for col,exp in [('H',m['score']),('I',m['mean']),('J',m['cost']),('K',m['minimum'])]:
   got=cells[3][f'{col}{i+6}'];assert got['formula'] and math.isclose(got['value'],exp,abs_tol=1e-8),(m['model'],col,got)
 assert not errors,errors
 assert math.isclose(cells[3]['E120']['value'],116.842,abs_tol=1e-8)
 chart_names=[n for n in z.namelist() if '/charts/chart' in n and n.endswith('.xml')]
 assert len(chart_names)==4
 binding=json.loads((WORK/'chart_bindings.json').read_text(encoding='utf-8'))
 for name,expected in zip(chart_names,binding):
  ch=ET.fromstring(z.read(name));series=ch.findall('.//c:ser',CN);assert len(series)==8
  for s,p in zip(series,expected['points']):
   xv=s.find('c:xVal/c:numLit/c:pt/c:v',CN);yv=s.find('c:yVal/c:numRef/c:numCache/c:pt/c:v',CN)
   assert math.isclose(float(xv.text),p['x'],abs_tol=1e-10)
   assert math.isclose(float(yv.text),p['y'],abs_tol=1e-10)
   assert s.find('c:yVal/c:numRef/c:f',CN) is not None
qa={'answers_verified':24,'models':8,'questions_per_model':3,'sheets':names,'native_scatter_charts':4,'native_scatter_points':32,'formula_errors':errors,'adopted_cost_sum':data['totalCost'],'opus_cost_counted_once':50.77,'supplemental_added_cost':5.617,'plot_pairs':4,'native_chart_x_values':'Current cost snapshot; worksheet note documents refresh requirement'}
(ROOT/'.work').mkdir(exist_ok=True)
(ROOT/'.work/plot_qa.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False,indent=2))
