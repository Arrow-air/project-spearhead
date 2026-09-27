"""Open the Build123D-generated STEP as an unsaved review document in Fusion."""
import adsk.core,adsk.fusion,json,traceback
from pathlib import Path
repo=Path('/Users/hex/projects/project-spearhead');out=repo/'reviews/fusion-reconstruction'
a=adsk.core.Application.get()
try:
 opts=a.importManager.createSTEPImportOptions(str(repo/'src/design/output/spearhead.step'))
 doc=a.importManager.importToNewDocument(opts)
 d=adsk.fusion.Design.cast(doc.products.itemByProductType('DesignProductType'))
 cam=a.activeViewport.camera;cam.viewOrientation=adsk.core.ViewOrientations.IsoTopRightViewOrientation;cam.isFitView=True;a.activeViewport.camera=cam
 a.activeViewport.fit();a.activeViewport.saveAsImageFile(str(out/'fusion-roundtrip.png'),2000,1400)
 r={'document':doc.name,'saved':doc.isSaved,'components':d.allComponents.count,'occurrences':d.rootComponent.allOccurrences.count,'root_bodies':d.rootComponent.bRepBodies.count,'body_count':sum(c.bRepBodies.count for c in d.allComponents)}
 (out/'fusion-roundtrip.json').write_text(json.dumps(r,indent=2))
except:(out/'fusion-roundtrip-error.txt').write_text(traceback.format_exc())
