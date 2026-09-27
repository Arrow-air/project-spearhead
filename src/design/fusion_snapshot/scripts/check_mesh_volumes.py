"""Independent native mesh-volume check of anomalous STEP mass properties."""
import adsk.core,adsk.fusion,json,traceback
from pathlib import Path
out=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot')
a=adsk.core.Application.get();results=[]
checks=[('main_wing_structure','urn:adsk.wipprod:dm.lineage:3rj6o5Z1QYSj28d3Cyrr_w','root',['Body145','Body182']),('tail_wing_structure_v2_USER_EDITS_NOW','urn:adsk.wipprod:dm.lineage:O4Ca7m3FR7iOTq82--njbQ','tail:1+ribs:1',['Body3','Body6']),('tail_wing_structure_v2_USER_EDITS_NOW','urn:adsk.wipprod:dm.lineage:O4Ca7m3FR7iOTq82--njbQ','tail_connector:1',['Body1'])]
for slug,fid,path,names in checks:
 doc=a.documents.open(a.data.findFileById(fid).latestVersion,True);d=adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
 bodies=d.rootComponent.bRepBodies if path=='root' else next(o for o in d.rootComponent.allOccurrences if o.fullPathName==path).bRepBodies
 for b in bodies:
  if b.name not in names:continue
  rec={'slug':slug,'path':path,'body':b.name};results.append(rec)
  try:
   calc=b.meshManager.createMeshCalculator();calc.surfaceTolerance=.001
   mesh=calc.calculate();v=list(mesh.nodeCoordinatesAsDouble);indices=list(mesh.nodeIndices);vol=0
   for i in range(0,len(indices),3):
    q=[v[3*k:3*k+3] for k in indices[i:i+3]];x,y,z=q
    vol+=x[0]*(y[1]*z[2]-y[2]*z[1])+x[1]*(y[2]*z[0]-y[0]*z[2])+x[2]*(y[0]*z[1]-y[1]*z[0])
   rec.update(triangles=mesh.triangleCount,mesh_volume_mm3=abs(vol/6)*1000,source_modified=doc.isModified)
  except:rec['error']=traceback.format_exc()
  (out/'native-mesh-checks.json').write_text(json.dumps(results,indent=2))
# Inspect the placeholder without relying on activeProduct.
f=a.data.findFileById('urn:adsk.wipprod:dm.lineage:Pvg6BwuAQZSY2p7yAcUTTw').latestVersion
doc=a.documents.open(f,True);doc.activate();adsk.doEvents()
a.userInterface.workspaces.itemById('FusionSolidEnvironment').activate();adsk.doEvents()
design=adsk.fusion.Design.cast(a.activeProduct)
r={'name':doc.name,'design_available':bool(design),'modified':doc.isModified}
(out/'temp-inspection.json').write_text(json.dumps(r,indent=2))
