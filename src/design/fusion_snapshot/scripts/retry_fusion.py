"""Run with exec(open(__file__).read(), {}) in Fusion's Python console.
Read-only snapshot of every f3d in the current Project Spearhead cloud project.
Does not update external references or save source documents.
"""
import adsk.core, adsk.fusion, json, re, traceback, datetime
from pathlib import Path
OUT=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot')
CATALOG=OUT/'catalog.json'
OUT.mkdir(parents=True,exist_ok=True)
app=adsk.core.Application.get()
files=[f for f in json.loads(CATALOG.read_text())['files'] if (f['path']=='Project Spearhead' or f['path'].startswith('Project Spearhead/')) and f['extension']=='f3d']
priority=['fuselage_v2','main_wing_structure','tail_wing_structure_v2_USER_EDITS_NOW','conceptual_design']
files.sort(key=lambda f: priority.index(f['name']) if f['name'] in priority else 99)
report=json.loads((OUT/'index.json').read_text());failed={f['id'] for f in report['files'] if f['status']=='failed'};files=[f for f in files if f['id'] in failed];report['files']=[f for f in report['files'] if f['id'] not in failed];report['status']='running'
def save(): (OUT/'index.json').write_text(json.dumps(report,indent=2))
def xyz(p):return [x*10 for x in p.asArray()]
def bounds(b):return [xyz(b.minPoint),xyz(b.maxPoint)]
def body(b):
 r={'name':b.name,'visible':b.isVisible,'solid':b.isSolid,'bounds_mm':bounds(b.preciseBoundingBox)}
 try:r['volume_mm3']=b.physicalProperties.volume*1000
 except:pass
 return r
for f in files:
 slug=re.sub('[^A-Za-z0-9_-]+','_',f['name']);out=OUT/slug;out.mkdir(exist_ok=True)
 r=dict(f,slug=slug,status='opening');report['files'].append(r);save()
 try:
  df=app.data.findFileById(f['id']).latestVersion
  doc=app.documents.open(df,True);d=adsk.fusion.Design.cast(app.activeProduct)
  if not d:
   app.userInterface.workspaces.itemById('FusionSolidEnvironment').activate();adsk.doEvents();d=adsk.fusion.Design.cast(app.activeProduct)
  r.update(version=df.versionNumber,modified_before=doc.isModified)
  meta={'file':dict(f),'document':doc.name,'root_id':d.rootComponent.id,'components':[],'occurrences':[],'parameters':[],'units':'mm; transform translations in cm'}
  if d.designType == adsk.fusion.DesignTypes.ParametricDesignType:
   for p in d.userParameters:meta['parameters'].append({'name':p.name,'expression':p.expression,'value_internal':p.value,'unit':p.unit})
  for c in d.allComponents:
   cr={'id':c.id,'name':c.name,'bodies':[],'mesh_count':c.meshBodies.count}
   for b in c.bRepBodies:
    try:cr['bodies'].append(body(b))
    except:cr.setdefault('errors',[]).append(traceback.format_exc())
   meta['components'].append(cr)
  for o in d.rootComponent.allOccurrences:
   oc={'path':o.fullPathName,'component_id':o.component.id,'visible':o.isVisible,'transform_cm':list(o.transform2.asArray()),'bodies':[]}
   for b in o.bRepBodies:
    try:oc['bodies'].append(body(b))
    except:oc.setdefault('errors',[]).append(traceback.format_exc())
   if o.isReferencedComponent:
    try:
     rf=o.documentReference.dataFile;oc['reference']={'id':rf.id,'name':rf.name,'version':rf.versionNumber,'latest':rf.latestVersionNumber}
    except:oc['reference_error']=traceback.format_exc()
   meta['occurrences'].append(oc)
  (out/'manifest.json').write_text(json.dumps(meta,indent=2))
  r.update(status='exporting',components=len(meta['components']),occurrences=len(meta['occurrences']));save()
  r['step_ok']=d.exportManager.execute(d.exportManager.createSTEPExportOptions(str(out/'source.step')))
  app.activeViewport.fit();app.activeViewport.saveAsImageFile(str(out/'fusion.png'),1400,1000)
  r.update(status='complete',modified_after=doc.isModified);save()
 except:r.update(status='failed',error=traceback.format_exc());save()
 adsk.doEvents()
report['status']='complete';save()
