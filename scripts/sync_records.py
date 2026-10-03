"""Regenerate adopted CSV and PDF manifest without changing frozen scores."""
from pathlib import Path
import csv,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'data/analysis_data.json').read_text(encoding='utf-8'))
p=ROOT/'data/submissions.json'
manifest=json.loads(p.read_text(encoding='utf-8'))
byid={r['id']:r for r in d['rows']}
seen={r['id'] for r in manifest['submissions']}
for r in d['rows']:
 if r['id'] not in seen:
  manifest['submissions'].append({'id':r['id'],'question':r['question'],'model':r['model'],'adopted':True,'replaced_by':None,'path':r['source_pdf'],'sha256':hashlib.sha256((ROOT/r['source_pdf']).read_bytes()).hexdigest()})
for r in manifest['submissions']:
 if r['adopted']:
  a=byid[r['id']]
  assert r['path']==a['source_pdf'] and r['model']==a['model']
manifest['pdf_count']=len(manifest['submissions'])
manifest['adopted_count']=sum(r['adopted'] for r in manifest['submissions'])
p.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
dims=['definition','calibration','verification','research','calibration_of_claims']
fields=['id','question','model',*dims,'total','interval_lower','interval_upper','cost_usd','source_pdf','replaces']
with (ROOT/'data/scores.csv').open('w',encoding='utf-8',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=fields)
 writer.writeheader()
 for r in d['rows']:
  writer.writerow({'id':r['id'],'question':r['question'],'model':r['model'],**r['scores'],'total':r['total'],'interval_lower':r['interval'][0],'interval_upper':r['interval'][1],'cost_usd':f"{r['cost']:.3f}",'source_pdf':r['source_pdf'],'replaces':r['replaces'] or ''})
print(f"Synced {len(d['rows'])} adopted CSV rows and {manifest['pdf_count']} PDF manifest entries.")
