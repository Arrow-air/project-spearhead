import json
from pathlib import Path
from scipy.optimize import linear_sum_assignment
import numpy as np
base=Path(__file__).resolve().parents[1]
for slug in ['fuselage_v2','main_wing_structure','tail_wing_structure_v2_USER_EDITS_NOW','lg','motor_mount']:
 m=json.loads((base/slug/'manifest.json').read_text());bs=[]
 for c in m['components']:
  if c['id']==m['root_id']:bs += [dict(b,path='root') for b in c['bodies'] if b['visible'] and b['solid']]
 for o in m['occurrences']:
  if o['visible']:bs += [dict(b,path=o['path']) for b in o['bodies'] if b['visible'] and b['solid']]
 ss=json.loads((base/slug/'import-metrics.json').read_text())
 cost=np.array([[np.max(np.abs(np.array(b['bounds_mm'])-np.array([s['min'],s['max']]))) for b in bs] for s in ss])
 ri,ci=linear_sum_assignment(cost)
 rows=[dict(ss[i],source=bs[j],bounds_error=float(cost[i,j]),volume_error=ss[i]['volume']-bs[j]['volume_mm3']) for i,j in zip(ri,ci)]
 report={'matched':rows,'missing_native':[b for j,b in enumerate(bs) if j not in ci]}
 (base/slug/'body-map.json').write_text(json.dumps(report,indent=2))
 print(slug,'worst bbox',max(x['bounds_error'] for x in rows),'missing',report['missing_native'])
 for x in rows:
  if abs(x['volume_error'])>100 or x['bounds_error']>.1:print(x['index'],x['source']['path'],x['source']['name'],'bbox',x['bounds_error'],'vol',x['volume_error'])
