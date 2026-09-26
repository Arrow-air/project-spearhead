from pathlib import Path
import json
from build123d import import_step
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from spearhead import volume_mm3
b=Path(__file__).resolve().parents[1]
for slug in ['fuselage_v2','main_wing_structure','tail_wing_structure_v2_USER_EDITS_NOW','lg','motor_mount','conceptual_design']:
 s=import_step(b/slug/'source.step');ss=list(s.solids());print(slug,'solids',len(ss),'volume',sum(x.volume for x in ss),'valid',sum(x.is_valid for x in ss),'bounds',s.bounding_box())
 data=[]
 for i,x in enumerate(ss):
  bb=x.bounding_box();data.append({'index':i,'volume':volume_mm3(x),'default_nonadaptive_volume':x.volume,'min':list(bb.min),'max':list(bb.max)})
 (b/slug/'import-metrics.json').write_text(json.dumps(data,indent=2))
