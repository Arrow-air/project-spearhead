"""Recover the hidden fixed-wing ribs and rear gear from fuselage_v2 v41.

Source stays read-only. Only selected occurrence-context bodies are copied to
separate disposable direct-design documents; no duplicate center spar is added.
"""
import adsk.core,adsk.fusion,json,traceback
from pathlib import Path
OUT=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot/fuselage_v2')
a=adsk.core.Application.get();result={'groups':[]};scratch=None
try:
 doc=a.documents.open(a.data.findFileById('urn:adsk.wipprod:dm.lineage:IS_VnAkqRZqj7rEFbFcCzA').latestVersion,True)
 d=adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
 result.update(source_name=doc.name,source_version=doc.dataFile.versionNumber,source_modified_before=doc.isModified)
 assert doc.dataFile.versionNumber==41,'Reinspect new source version before mixing snapshots'
 selections=[('hidden_inner_ribs','fuselage:1+wing_inner:1',['Body13','Body14','Body15','Body16','Body17','Body18']),('hidden_rear_gear','lg:1',[f'Body{i} (1)' for i in range(1,12)])]
 manager=adsk.fusion.TemporaryBRepManager.get()
 for slug,path,names in selections:
  occ=next(o for o in d.rootComponent.allOccurrences if o.fullPathName==path)
  bodies=[b for b in occ.bRepBodies if b.name in names]
  assert len(bodies)==len(names),(slug,len(bodies))
  rec={'slug':slug,'occurrence':path,'transform_cm':list(occ.transform2.asArray()),'bodies':[]};result['groups'].append(rec)
  copies=[]
  for b in bodies:
   p=b.getPhysicalProperties(adsk.fusion.CalculationAccuracy.VeryHighCalculationAccuracy);bb=b.preciseBoundingBox
   rec['bodies'].append({'name':b.name,'visible':b.isVisible,'solid':b.isSolid,'volume_mm3':p.volume*1000,'bounds_mm':[[v*10 for v in pt.asArray()] for pt in [bb.minPoint,bb.maxPoint]]})
   copies.append((b.name,manager.copy(b)))
  scratch=a.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
  target=adsk.fusion.Design.cast(scratch.products.itemByProductType('DesignProductType'));target.designType=adsk.fusion.DesignTypes.DirectDesignType
  for name,body in copies:
   # One named component per body keeps picking meaningful in STEP/glTF.
   child=target.rootComponent.occurrences.addNewComponent(adsk.core.Matrix3D.create()).component;child.name=name
   created=child.bRepBodies.add(body);created.name=name
  rec['step_ok']=target.exportManager.execute(target.exportManager.createSTEPExportOptions(str(OUT/(slug+'.step'))))
  scratch.close(False);scratch=None
  (OUT/'center-gear-recovery.json').write_text(json.dumps(result,indent=2))
 result['source_modified_after']=doc.isModified
 result['status']='complete'
except:result.update(status='failed',error=traceback.format_exc())
finally:
 if scratch:scratch.close(False)
(OUT/'center-gear-recovery.json').write_text(json.dumps(result,indent=2))
