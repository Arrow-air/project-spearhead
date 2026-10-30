"""Generate a repeatable review image of the Build123D assembly (no GUI needed)."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from OCP.BRep import BRep_Tool
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopAbs import TopAbs_Orientation
from OCP.TopLoc import TopLoc_Location
from spearhead import SYSTEMS, build_aircraft

def triangles(solid):
    """Triangles of a solid. OCP leaves a few small cone faces unmeshed; a preview can skip them."""
    BRepMesh_IncrementalMesh(solid.wrapped,1.0,False,0.3,True);result=[]
    for face in solid.faces():
        location=TopLoc_Location();poly=BRep_Tool.Triangulation_s(face.wrapped,location)
        if poly is None:continue
        trsf=location.Transformation();nodes=[poly.Node(i).Transformed(trsf) for i in range(1,poly.NbNodes()+1)]
        nodes=np.array([(p.X(),p.Y(),p.Z()) for p in nodes])
        for i in range(1,poly.NbTriangles()+1):
            a,b,c=poly.Triangle(i).Get()
            if face.wrapped.Orientation()==TopAbs_Orientation.TopAbs_REVERSED:b,c=c,b
            result.append(nodes[[a-1,b-1,c-1]])
    return np.array(result)

def render(output):
    aircraft=build_aircraft()
    fig=plt.figure(figsize=(16,10),facecolor='#101a24')
    ax=fig.add_axes([0,.13,1,.68],projection='3d',facecolor='#101a24');ax.set_proj_type('ortho')
    # One collection, so faces are depth-sorted across solids rather than per solid.
    faces,colors=[],[]
    for subsystem in aircraft.children:
        for solid in subsystem.solids():
            mesh=triangles(solid);faces.append(mesh);colors+=[SYSTEMS[subsystem.label]]*len(mesh)
    # Drop faces turned away from the camera, or the far side of the skins shows through.
    faces=np.concatenate(faces);normals=np.cross(faces[:,1]-faces[:,0],faces[:,2]-faces[:,0])
    elev,azim=np.radians(28),np.radians(-42);eye=np.array([np.cos(elev)*np.cos(azim),np.cos(elev)*np.sin(azim),np.sin(elev)])
    front=normals@eye>0
    mesh=Poly3DCollection(faces[front],facecolors=np.array(colors)[front],linewidth=0,shade=True,lightsource=matplotlib.colors.LightSource(azdeg=130,altdeg=50))
    mesh.set_edgecolor('none')
    ax.add_collection3d(mesh)
    bb=aircraft.bounding_box();padding=80
    ax.set(xlim=(bb.min.X-padding,bb.max.X+padding),ylim=(bb.min.Y-padding,bb.max.Y+padding),zlim=(bb.min.Z-padding,bb.max.Z+padding))
    ax.set_box_aspect(tuple(bb.size),zoom=1.35);ax.view_init(elev=28,azim=-42);ax.set_axis_off()
    fig.text(.05,.93,'SPEARHEAD',color='#e7ebee',fontsize=27,fontweight='bold')
    fig.text(.05,.891,'PT1 main assembly · Fusion v35 export',color='#b8c5cf',fontsize=12)
    fig.text(.05,.085,'3.95 m span  |  490 solids, 215 part numbers  |  Fusion coordinates retained',color='#b8c5cf',fontsize=11)
    fig.text(.05,.053,'Review model, not a manufacturing release. Battery, avionics and harness are not modelled.',color='#e7b98b',fontsize=10)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=150,facecolor=fig.get_facecolor(),bbox_inches='tight');plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).parent/'output/spearhead.png');render(p.parse_args().output)
