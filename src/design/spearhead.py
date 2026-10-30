"""Spearhead PT1 main assembly, read from the Fusion v35 STEP export.

Coordinates: mm, X aft, Y span, Z up. The STEP in pt1/ is the source of truth:
Fusion design AR-SH-PT1-0000 PT1 v35, exported without the hidden REF group.
This is a review assembly, not a manufacturing release.
"""
from __future__ import annotations
import argparse
import csv
import re
from pathlib import Path
import json
from build123d import Color, Compound, export_step, import_step

HERE = Path(__file__).resolve().parent
PT1 = HERE / 'pt1'
SOURCE = PT1 / 'AR-SH-PT1-0000_PT1_v35_noREF.step'
PART_NUMBERS = PT1 / 'pt1_part_numbers.csv'
SYSTEMS = {
    '10000-FuselageStructure': '#cead76',
    '20000-WingStructure': '#d8c49b',
    '30000-TailStructure': '#87b4cf',
    '40000-Propulsion': '#c8614d',
    '80000-BoomStructure': '#657986',
}
PART = re.compile(r'^(\d{3}(?!00)\d{2}[LR]?)-')


def _tree(node, color):
    """Rebuild a node with its solids as labeled children; export needs that for parts."""
    if node.children:
        children = [_tree(child, color) for child in node.children if child.solids()]
        return Compound(label=node.label, children=children)
    solids = list(node.solids())
    for solid in solids:
        solid.label = node.label
        solid.color = Color(color)
    return Compound(label=node.label, children=solids)


def build_aircraft():
    """Load the exported assembly and colour each system. Names come from Fusion."""
    raw = import_step(SOURCE)
    found = [c.label for c in raw.children]
    if found != list(SYSTEMS):
        raise ValueError(f'unexpected top-level systems in {SOURCE.name}: {found}')
    systems = [_tree(system, SYSTEMS[system.label]) for system in raw.children]
    return Compound(label='Spearhead_PT1_v35', children=systems)


def part_occurrences(node, parents=()):
    """Yield (part number, component name, grouping path, node) for each numbered part."""
    number = PART.match(node.label or '')
    if number:
        yield number.group(1), node.label, parents, node
        return
    for child in node.children:
        yield from part_occurrences(child, parents + (node.label,))


def part_numbers():
    with PART_NUMBERS.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


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
    parser.add_argument('--show',action='store_true',help='Show in an existing OCP CAD Viewer')
    args=parser.parse_args()
    aircraft=build_aircraft()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    export_step(aircraft,args.output)
    report=summary(aircraft)
    args.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    if args.show:
        from ocp_vscode import show
        show(aircraft)

if __name__=='__main__':main()
