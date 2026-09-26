"""Generate a repeatable review image of the Build123D assembly (no GUI needed)."""
from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from spearhead import build_aircraft

def render(output):
    aircraft=build_aircraft()
    fig=plt.figure(figsize=(16,10),facecolor='#101a24')
    ax=fig.add_axes([0,.13,1,.68],projection='3d',facecolor='#101a24');ax.set_proj_type('ortho')
    palette=['#cead76','#d8c49b','#cfb27e','#657986','#ad987b','#87b4cf','#75c3ba']
    for subsystem,color in zip(aircraft.children,palette):
        for solid in subsystem.solids():
            vertices,triangles=solid.tessellate(1.0,0.3)
            vertices=np.array([tuple(v) for v in vertices])
            mesh=Poly3DCollection(vertices[np.array(triangles)],facecolors=color,linewidth=0,shade=True,lightsource=matplotlib.colors.LightSource(azdeg=130,altdeg=50))
            mesh.set_edgecolor('none')
            ax.add_collection3d(mesh)
    bb=aircraft.bounding_box();padding=80
    ax.set(xlim=(bb.min.X-padding,bb.max.X+padding),ylim=(bb.min.Y-padding,bb.max.Y+padding),zlim=(bb.min.Z-padding,bb.max.Z+padding))
    ax.set_box_aspect(tuple(bb.size),zoom=1.35);ax.view_init(elev=28,azim=-42);ax.set_axis_off()
    fig.text(.05,.93,'SPEARHEAD',color='#e7ebee',fontsize=27,fontweight='bold')
    fig.text(.05,.891,'Build123D structural reconstruction · Fusion snapshot 26 Sep 2026',color='#b8c5cf',fontsize=12)
    fig.text(.05,.085,'3.95 m span  |  Native coordinates retained  |  Blue: inferred mirrored wing',color='#b8c5cf',fontsize=11)
    fig.text(.05,.053,'Review model — not as-built. No fitted propulsion, electronics or covering skins.',color='#e7b98b',fontsize=10)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=150,facecolor=fig.get_facecolor(),bbox_inches='tight');plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).parent/'output/spearhead.png');render(p.parse_args().output)
