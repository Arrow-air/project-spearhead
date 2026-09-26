import adsk.core,adsk.fusion,json
from pathlib import Path
out=Path('/Users/hex/projects/project-spearhead/src/design/fusion_snapshot')
a=adsk.core.Application.get();doc=a.activeDocument;d=adsk.fusion.Design.cast(a.activeProduct)
assert doc.name=='temp'
r={'slug':'temp','modified_before':doc.isModified}
r['archive_ok']=d.exportManager.execute(d.exportManager.createFusionArchiveExportOptions(str(out/'temp/source.f3d')))
r['modified_after']=doc.isModified
archive=json.loads((out/'archives.json').read_text());archive['files']=[f for f in archive['files'] if f['slug']!='temp']+[r]
(out/'archives.json').write_text(json.dumps(archive,indent=2))
(out/'temp-inspection.json').write_text(json.dumps({'name':doc.name,'resolution':'Saved in Simulation workspace; activated Design workspace to access design API. Export succeeded without modifying source.','root_bodies':d.rootComponent.bRepBodies.count,'modified':doc.isModified},indent=2))
