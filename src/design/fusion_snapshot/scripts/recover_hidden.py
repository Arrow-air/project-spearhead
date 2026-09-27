"""Export hidden fourth mount using a disposable copy; never changes source."""
import adsk.core,adsk.fusion,json,traceback
from pathlib import Path
out=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot/motor_mount')
a=adsk.core.Application.get();r={};scratch=None
try:
 source=a.documents.open(a.data.findFileById('urn:adsk.wipprod:dm.lineage:qZcgavhdTMSz-nhIpxyYcw').latestVersion,True)
 design=adsk.fusion.Design.cast(source.products.itemByProductType('DesignProductType'))
 r['source_modified_before']=source.isModified
 occ=next(o for o in design.rootComponent.allOccurrences if o.fullPathName=='motor_mounts:1')
 bodies=[b for b in occ.bRepBodies if b.name in ['Body6','Body7','Body8']]
 r['bodies']=[b.name for b in bodies]
 manager=adsk.fusion.TemporaryBRepManager.get();copies=[(b.name,manager.copy(b)) for b in bodies]
 scratch=a.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
 target=adsk.fusion.Design.cast(scratch.products.itemByProductType('DesignProductType'));target.designType=adsk.fusion.DesignTypes.DirectDesignType
 for name,body in copies:
  new=target.rootComponent.bRepBodies.add(body);new.name=name
 r['step_ok']=target.exportManager.execute(target.exportManager.createSTEPExportOptions(str(out/'hidden_fourth_mount.step')))
 r['source_modified_after']=source.isModified
except:r['error']=traceback.format_exc()
finally:
 if scratch:scratch.close(False)
(out/'hidden-recovery.json').write_text(json.dumps(r,indent=2))
