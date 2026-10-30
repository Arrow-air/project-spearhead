"""Validate the PT1 export against the Fusion body dump and the part list; record the result."""
from collections import Counter
import csv, json, hashlib
from build123d import import_step, export_step
from spearhead import HERE, PT1, SOURCE, SYSTEMS, build_aircraft, part_numbers, part_occurrences, summary

REVIEW=HERE.parents[1]/'reviews'/'pt1-assembly'
NATIVE=PT1/'native-bodies.json'
PARTS,OCCURRENCES,SOLIDS,NUMBERED_BODIES=215,426,490,482
SPLIT_BODIES={'15105','15205'} # one Fusion body of 3 disjoint lumps; STEP writes 3 solids

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def measure(solids):
    """Adaptive volume and centroid of a set of solids (see spearhead.volume_mm3)."""
    from OCP.BRepGProp import BRepGProp
    from OCP.GProp import GProp_GProps
    total=0;moment=[0,0,0]
    for s in solids:
        props=GProp_GProps();BRepGProp.VolumeProperties_s(s.wrapped,props,1e-9,True,False)
        c=props.CentreOfMass();total+=props.Mass();moment=[m+props.Mass()*v for m,v in zip(moment,(c.X(),c.Y(),c.Z()))]
    return total,[m/total for m in moment]

def mesh_bounds(solids):
    """Bounds from a 0.02 mm mesh; OCP's analytic box is up to 59 mm loose on trimmed spline faces."""
    from OCP.BRep import BRep_Tool
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.TopLoc import TopLoc_Location
    low=[float('inf')]*3;high=[-float('inf')]*3
    for s in solids:
        BRepMesh_IncrementalMesh(s.wrapped,.02,False,.1,True)
        for f in s.faces():
            location=TopLoc_Location();poly=BRep_Tool.Triangulation_s(f.wrapped,location)
            if poly is None:continue
            trsf=location.Transformation()
            for i in range(1,poly.NbNodes()+1):
                p=poly.Node(i).Transformed(trsf);xyz=(p.X(),p.Y(),p.Z())
                low=[min(a,b) for a,b in zip(low,xyz)];high=[max(a,b) for a,b in zip(high,xyz)]
    return low+high

def side(y):return 'qty L' if y<-20 else 'qty R' if y>20 else 'qty centre'

def native_occurrences(native):
    """Group the Fusion bodies per occurrence: component name to volume and centroid."""
    grouped={}
    for b in native['bodies']:grouped.setdefault(b['path'],[]).append(b)
    result={}
    for path,bodies in grouped.items():
        total=sum(b['vol'] for b in bodies)
        bounds=[min(b['bbox'][k] for b in bodies) for k in range(3)]+[max(b['bbox'][k] for b in bodies) for k in range(3,6)]
        result.setdefault(bodies[0]['comp'],[]).append((total,[sum(b['c'][k]*b['vol'] for b in bodies)/total for k in range(3)],bounds))
    return result

def run():
    REVIEW.mkdir(parents=True,exist_ok=True)
    report={'source_files':[],'native_comparison':{},'checks':[],'warnings':[]}
    def check(name,condition):
        report['checks'].append({'name':name,'passed':bool(condition)})
    for path in [SOURCE,*sorted(PT1.glob('pt1_*.csv')),NATIVE]:
        report['source_files'].append({'name':path.name,'sha256':digest(path),'bytes':path.stat().st_size})
    native=json.loads(NATIVE.read_text())
    check('native dump is from Fusion v35 at very high accuracy with nothing unsaved',native['version']==35 and native['accuracy']=='very high' and not native['modified'])
    check(f'native dump has {NUMBERED_BODIES} numbered bodies',len(native['bodies'])==NUMBERED_BODIES)
    with (PT1/'pt1_materials_table.csv').open(encoding='utf-8',newline='') as stream:
        materials={r['part number']:set(r['material'].split('; ')) for r in csv.DictReader(stream)}
    check('every numbered body has the material of the materials table',all(b['mat'] in materials.get(b['pn'],()) for b in native['bodies']))
    aircraft=build_aircraft();report['assembly']=summary(aircraft)
    check('five systems in tree order',[c.label for c in aircraft.children]==list(SYSTEMS))
    check(f'assembly contains {SOLIDS} solids',report['assembly']['solids']==SOLIDS)
    check(f'all {SOLIDS} assembly solids pass OCP validity',report['assembly']['valid_solids']==SOLIDS)
    check('3950 mm span',abs(report['assembly']['size_mm'][1]-3950)<.001)
    occurrences=list(part_occurrences(aircraft))
    report['invalid_on_read']=[o[1] for o in occurrences if not all(x.is_valid for x in o[3].solids())]
    check(f'{OCCURRENCES} part occurrences',len(occurrences)==OCCURRENCES)
    rows=part_numbers();by_number={r['number']:r for r in rows}
    check(f'{PARTS} part numbers in the list',len(rows)==len(by_number)==PARTS)
    found=Counter(o[0] for o in occurrences)
    check('every part in the STEP is in the list and vice versa',set(found)==set(by_number))
    check('list quantities match the STEP',all(int(by_number[n]['qty total'])==q for n,q in found.items()))
    check('list quantities sum to the occurrences',sum(int(r['qty total']) for r in rows)==OCCURRENCES)
    pool=native_occurrences(native);sides=Counter();differences=[];all_differences=[];centroid_errors=[];bounds_errors=[];unmatched=[]
    for number,name,parents,node in occurrences:
        total,c=measure(node.solids());sides[number,side(c[1])]+=1
        candidates=pool.get(name,[])
        if not candidates:unmatched.append(name);continue
        best=min(candidates,key=lambda n:max(abs(a-b) for a,b in zip(n[1],c)));candidates.remove(best)
        delta=(total-best[0])/best[0]*100;all_differences.append(abs(delta));centroid_errors.append(max(abs(a-b) for a,b in zip(best[1],c)))
        centroid_error=centroid_errors[-1];bounds_errors.append(max(abs(a-b) for a,b in zip(best[2],mesh_bounds(node.solids()))));bounds_error=bounds_errors[-1]
        if abs(delta)>.1 or centroid_error>.05 or bounds_error>.02:differences.append({'name':name,'native_volume_mm3':best[0],'step_volume_mm3':total,'difference_percent':delta,'centroid_error_mm':centroid_error,'bounds_error_mm':bounds_error})
    leftovers=[k for k,v in pool.items() if v]
    report['native_comparison']={'matched_occurrences':len(all_differences),'unmatched_step':unmatched,'unmatched_native':leftovers,'max_adaptive_volume_difference_percent':max(all_differences),'max_centroid_error_mm':max(centroid_errors),'max_bounds_error_mm':max(bounds_errors),'over_limits':differences}
    check('every occurrence matched one native occurrence, both ways',not unmatched and not leftovers)
    check('adaptive volume within 0.1 percent of native',max(all_differences)<.1)
    check('centroid within 0.05 mm of native',max(centroid_errors)<.05)
    check('matched body bounds within 0.02 mm',max(bounds_errors)<=.02)
    check('list left, right and centre quantities match the STEP',all(int(r[k])==sides[r['number'],k] for r in rows for k in ('qty L','qty R','qty centre')))
    split=[o for o in occurrences if o[0] in SPLIT_BODIES]
    check('split bodies: 2 of each brace, 3 solids each',Counter(o[0] for o in split)==Counter({n:2 for n in SPLIT_BODIES}) and all(len(o[3].solids())==3 for o in split))
    check(f'{SOLIDS} solids less split lumps equals {NUMBERED_BODIES} bodies',report['assembly']['solids']-2*len(split)==NUMBERED_BODIES)
    check('four lift motors and four ESCs',found['41101']==4 and found['41301']==4)
    check('right ruddervator horn group present',all(found[n]==2 for n in ('33301','33302','33303')))
    output=HERE/'output';output.mkdir(exist_ok=True)
    export_step(aircraft,output/'spearhead.step')
    (output/'spearhead.json').write_text(json.dumps(report['assembly'],indent=2))
    roundtrip=import_step(output/'spearhead.step');rt=summary(roundtrip);report['roundtrip']=rt
    report['invalid_after_roundtrip']=[o[1] for o in part_occurrences(roundtrip) if not all(x.is_valid for x in o[3].solids())]
    check(f'STEP roundtrip retains {SOLIDS} solids',rt['solids']==SOLIDS)
    check('STEP roundtrip solids remain valid',rt['valid_solids']==SOLIDS)
    error=max(abs(a-b) for aa,bb in zip(rt['bounds_mm'],report['assembly']['bounds_mm']) for a,b in zip(aa,bb))
    check('STEP roundtrip bounds within 0.02 mm',error<.02)
    volume_error=abs(rt['summed_solid_volume_mm3']/report['assembly']['summed_solid_volume_mm3']-1)
    check('STEP roundtrip summed volume within 0.01 percent',volume_error<.0001)
    report['warnings']=[
        'Default non-adaptive OCP volume is inaccurate on some spline surfaces; use spearhead.volume_mm3. Validity, volume and centroid checks do not prove physical fit.',
        'Battery, avionics and harness are not modelled. Bought parts (motors, ESCs, aileron hinges) have no datasheet mass yet.',
        'The hidden REF group of Fusion (130 bodies) is not in this export.'
    ]
    report['passed']=all(r['passed'] for r in report['checks'])
    (REVIEW/'validation.json').write_text(json.dumps(report,indent=2)+'\n',newline='\n')
    print(json.dumps({'passed':report['passed'],'checks':len(report['checks']),'failed':[r for r in report['checks'] if not r['passed']],'native':report['native_comparison'],'assembly':report['assembly']},indent=2))
    return report['passed']

if __name__=='__main__':
    raise SystemExit(0 if run() else 1)
