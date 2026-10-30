"""Generate browser assets from the same checked Build123D aircraft assembly."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from spearhead import build_aircraft, part_numbers
from build123d import export_gltf
base=Path(__file__).resolve().parent
assembly=build_aircraft()
assert len(assembly.solids())==490
export_gltf(assembly,base/'assets/spearhead-v3.glb',binary=True,linear_deflection=.35,angular_deflection=.25)
parts={r['number']:{'description':r['description'],'material':r['material'],'qty':int(r['qty total'])} for r in part_numbers()}
(base/'assets/model.json').write_text(json.dumps({'solids':490,'revision':'pt1-v35','span_m':3.95,'source':'AR-SH-PT1-0000 PT1 v35','groups':[{'id':c.label,'solids':len(c.solids())} for c in assembly.children],'parts':parts},indent=2)+chr(10),newline=chr(10))
print('Generated', (base/'assets/spearhead-v3.glb').stat().st_size,'bytes')
