"""Archive native designs and high-accuracy body metrics, without saving sources."""
import adsk.core,adsk.fusion,json,traceback
from pathlib import Path
OUT=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot')
app=adsk.core.Application.get();index=json.loads((OUT/'index.json').read_text())
result={'status':'running','files':[]}
def save(): (OUT/'archives.json').write_text(json.dumps(result,indent=2))
for f in index['files']:
 r={'slug':f['slug']};result['files'].append(r);save()
 try:
  df=app.data.findFileById(f['id']).latestVersion;doc=app.documents.open(df,True);d=adsk.fusion.Design.cast(app.activeProduct)
  if not d:
   app.userInterface.workspaces.itemById('FusionSolidEnvironment').activate();adsk.doEvents();d=adsk.fusion.Design.cast(app.activeProduct)
  r['modified_before']=doc.isModified
  try:r['archive_ok']=d.exportManager.execute(d.exportManager.createFusionArchiveExportOptions(str(OUT/f['slug']/'source.f3d')))
  except:r['archive_error']=traceback.format_exc()
  if f['slug'] in ['fuselage_v2','main_wing_structure','tail_wing_structure_v2_USER_EDITS_NOW']:
   rows=[]
   def record(b,path):
    try:
     p=b.getPhysicalProperties(adsk.fusion.CalculationAccuracy.VeryHighCalculationAccuracy)
     rows.append({'path':path,'name':b.name,'solid':b.isSolid,'visible':b.isVisible,'volume_mm3':p.volume*1000,'area_mm2':p.area*100})
    except:rows.append({'path':path,'name':b.name,'error':traceback.format_exc()})
   for b in d.rootComponent.bRepBodies:record(b,'root')
   for o in d.rootComponent.allOccurrences:
    if o.isVisible:
     for b in o.bRepBodies:record(b,o.fullPathName)
   (OUT/f['slug']/'high-accuracy.json').write_text(json.dumps(rows,indent=2))
  r['modified_after']=doc.isModified
 except:r['error']=traceback.format_exc()
 save();adsk.doEvents()
result['status']='complete';save()
