"""Spearhead structural reconstruction from the 2026-09-26 Fusion snapshot.

Coordinates: mm, X aft, Y span, Z up. Native source STEP assets are immutable.
Mirrored outer wing and restored hidden fourth motor mount are identified.
This is a review assembly, not an as-built or flight-ready definition.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from pathlib import Path
import json
from build123d import Color, Compound, Plane, export_step, import_step

HERE = Path(__file__).resolve().parent
SNAPSHOT = HERE / 'fusion_snapshot'
STRUCTURES = {
    'fuselage': ('fuselage_v2', '#cead76'),
    'main_wing': ('main_wing_structure', '#d8c49b'),
    'tail': ('tail_wing_structure_v2_USER_EDITS_NOW', '#cfb27e'),
    'motor_mounts_and_booms': ('motor_mount', '#657986'),
    'landing_gear': ('lg', '#ad987b'),
}


def source_structure(key: str) -> Compound:
    """Recover named solid children in aircraft coordinates, with import checks."""
    slug, color = STRUCTURES[key]
    raw = import_step(SNAPSHOT / slug / 'source.step')
    solids = list(raw.solids())
    mapping = json.loads((SNAPSHOT / slug / 'body-map.json').read_text())
    if len(solids) != len(mapping['matched']):
        raise ValueError(f'{slug}: source solid count changed; regenerate body map')
    groups = defaultdict(list)
    for rec in mapping['matched']:
        solid = solids[rec['index']]
        bb = solid.bounding_box()
        error = max(abs(a-b) for a,b in zip(list(bb.min)+list(bb.max),rec['min']+rec['max']))
        if error > 0.02:
            raise ValueError(f'{slug}: STEP body order/geometry changed at {rec["index"]}')
        source = rec['source']
        solid.label = source['name']
        solid.color = Color(color)
        groups[source['path']].append(solid)
    return Compound(label=key, children=[Compound(label=path, children=bodies)
                                          for path,bodies in groups.items()])


def build_aircraft(mirror_wing=True, restore_mount=True) -> Compound:
    """Assemble latest structures; keep the measured common Fusion coordinates."""
    fuselage = source_structure('fuselage')
    wing = source_structure('main_wing')
    tail = source_structure('tail')
    mounts = source_structure('motor_mounts_and_booms')
    landing_gear = source_structure('landing_gear')
    children = [fuselage, wing, tail, mounts, landing_gear]
    if mirror_wing:
        # Only the detailed outer panel is incomplete. Do not duplicate the
        # existing full-span spar, center geometry, or +/-550 mm boom fittings.
        mirrored = []
        for solid in (s for group in wing.children for s in group.children):
            if solid.bounding_box().max.Y <= -669.9:
                other = solid.mirror(Plane.XZ)
                other.label = f'{solid.label}_mirrored_inferred'
                other.color = Color('#87b4cf')
                mirrored.append(other)
        children.append(Compound(label='inferred_starboard_outer_wing', children=mirrored))
    if restore_mount:
        restored = import_step(SNAPSHOT / 'motor_mount' / 'hidden_fourth_mount.step')
        restored.label = 'restored_hidden_port_aft_motor_mount'
        restored.color = Color('#75c3ba')
        children.append(restored)
    return Compound(label='Spearhead_reconstruction_2026_09_26', children=children)


def volume_mm3(solid):
    """Adaptive integration avoids large default mass errors on spline faces."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    props = GProp_GProps()
    BRepGProp.VolumeProperties_s(solid.wrapped, props, 1e-9, True, False)
    return props.Mass()


def summary(aircraft):
    solids=list(aircraft.solids());bb=aircraft.bounding_box()
    return {'solids':len(solids),'valid_solids':sum(s.is_valid for s in solids),
            'summed_solid_volume_mm3':sum(volume_mm3(s) for s in solids),
            'bounds_mm':[list(bb.min),list(bb.max)],'size_mm':list(bb.size),
            'subsystems':[{ 'name':c.label,'solids':len(c.solids()) } for c in aircraft.children]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'output'/'spearhead.step')
    parser.add_argument('--as-exported',action='store_true',help='Only saved-visible geometry; no mirrored wing or restored hidden mount')
    parser.add_argument('--show',action='store_true',help='Show in an existing OCP CAD Viewer')
    args=parser.parse_args()
    aircraft=build_aircraft(not args.as_exported,not args.as_exported)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    export_step(aircraft,args.output)
    report=summary(aircraft)
    args.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if args.show:
        from ocp_vscode import show
        show(aircraft)

if __name__=='__main__':main()
