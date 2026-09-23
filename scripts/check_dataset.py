"""Check the published benchmark records and artifact provenance (standard library only)."""
from pathlib import Path
from decimal import Decimal
from itertools import product
import csv,hashlib,json,math,zipfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
load=lambda p:json.loads((ROOT/p).read_text(encoding='utf-8-sig'))
data=load('data/analysis_data.json');manifest=load('data/submissions.json')
rows=data['rows'];byid={r['id']:r for r in rows}
assert len(rows)==len(byid)==24
assert len({r['model'] for r in rows})==8
assert len(manifest['submissions'])==27
assert sum(r['adopted'] for r in manifest['submissions'])==24
assert manifest['adopted_count']==24 and manifest['pdf_count']==27
rubric={'definition':15,'calibration':30,'verification':15,'research':30,'calibration_of_claims':10}
source_batches=[load('evaluations/original/评分明细.json'),load('evaluations/supplement/补充评分明细.json')]
source_scores={r['id']:r for batch in source_batches for q in batch['questions'] for r in q['answers']}
for r in rows:
 assert set(r['scores'])==set(rubric)
 assert all(isinstance(v,int) and 0<=v<=rubric[k] for k,v in r['scores'].items())
 assert sum(r['scores'].values())==r['total']==source_scores[r['id']]['total']
 assert r['scores']==source_scores[r['id']]['scores']
 assert r['interval']==source_scores[r['id']]['interval']
 assert r['interval'][0]<=r['total']<=r['interval'][1]
 assert r['cost']>0 and (ROOT/r['source_pdf']).is_file()
for r in manifest['submissions']:
 p=ROOT/r['path'];assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],r['id']
 if r['adopted']:
  adopted=byid[r['id']];assert r['path']==adopted['source_pdf'] and r['model']==adopted['model'] and r['question']==adopted['question']
 else:assert byid[r['replaced_by']]['replaces']==r['id']
for m in data['full']:
 group=[r for r in rows if r['model']==m['model']]
 assert len(group)==3 and {r['question'] for r in group}=={1,2,3}
 assert sum(r['total'] for r in group)==m['score']
 assert math.isclose(m['mean'],m['score']/3)
 assert sum(Decimal(str(r['cost'])) for r in group)==Decimal(str(m['cost']))
assert sum(Decimal(str(r['cost'])) for r in rows)==Decimal('116.842')
with (ROOT/'data/scores.csv').open(encoding='utf-8',newline='') as f:csvrows=list(csv.DictReader(f))
assert len(csvrows)==24
for r in csvrows:
 a=byid[r['id']]
 assert r['model']==a['model'] and int(r['question'])==a['question']
 assert int(r['total'])==a['total'] and Decimal(r['cost_usd'])==Decimal(str(a['cost']))
 assert r['source_pdf']==a['source_pdf']
 for k in rubric:assert int(r[k])==a['scores'][k]
def efficient(pts):
 return [a for a in pts if not any(b['cost']<=a['cost'] and b['score']>=a['score'] and (b['cost']<a['cost'] or b['score']>a['score']) for b in pts)]
assert {m['model'] for m in efficient(data['full'])}=={m['model'] for m in data['fullFrontier']}
combinations=[{'ids':[r['id'] for r in rs],'cost':sum(Decimal(str(r['cost'])) for r in rs),'score':sum(r['total'] for r in rs)} for rs in product(*[[r for r in rows if r['question']==q] for q in [1,2,3]])]
assert len(combinations)==data['portfolioCount']==512
assert {(tuple(r['ids']),r['cost'],r['score']) for r in efficient(combinations)}=={(tuple(r['ids']),Decimal(str(r['cost'])),r['score']) for r in data['portfolioFrontier']}
ns={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(ROOT/'results/模型评分与成本汇总.xlsx') as z:
 strings=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))] if 'xl/sharedStrings.xml' in z.namelist() else []
 sheetnames=[e.attrib['name'] for e in ET.fromstring(z.read('xl/workbook.xml')).find('s:sheets',ns)]
 assert sheetnames==['Task 1','Task 2','Task 3','综合汇总']
 for q in [1,2,3]:
  cells={}
  for c in ET.fromstring(z.read(f'xl/worksheets/sheet{q}.xml')).findall('.//s:sheetData/s:row/s:c',ns):
   v=c.find('s:v',ns);typ=c.attrib.get('t')
   if typ=='s':value=strings[int(v.text)]
   elif typ=='inlineStr':value=''.join(c.find('s:is',ns).itertext())
   elif v is None:value=None
   elif typ=='str':value=v.text
   else:value=float(v.text)
   cells[c.attrib['r']]=value
  for r in [r for r in rows if r['question']==q]:
   for col,expected in zip('ABCD',[r['id'],r['model'],r['cost'],r['total']]):assert cells[f"{col}{r['row']}"]==expected,(r['id'],col)
print('Verified: 27 PDF hashes; 24 adopted answers; 8 complete configurations; USD 116.842; 512 combinations; CSV and workbook match.')
