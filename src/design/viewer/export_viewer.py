"""Generate browser assets from the same checked Build123D aircraft assembly."""
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from spearhead import build_aircraft
from build123d import export_gltf
base=Path(__file__).resolve().parent
assembly=build_aircraft()
assert len(assembly.solids())==318
export_gltf(assembly,base/'assets/spearhead-v2.glb',binary=True,linear_deflection=.35,angular_deflection=.25)
(base/'assets/model.json').write_text(json.dumps({'solids':318,'revision':'center-ribs-rear-gear-v2','span_m':3.95,'snapshot':'2026-09-26','groups':[{'id':c.label,'solids':len(c.solids())} for c in assembly.children]},indent=2))
print('Generated', (base/'assets/spearhead-v2.glb').stat().st_size,'bytes')
