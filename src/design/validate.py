"""Validate snapshot imports and reconstruction; record exchange discrepancies."""
from pathlib import Path
import json, hashlib
from build123d import import_step, export_step
from spearhead import HERE, SNAPSHOT, STRUCTURES, build_aircraft, summary, volume_mm3

REVIEW=HERE.parents[1]/'reviews'/'fusion-reconstruction'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def run(reuse_imports=False):
    REVIEW.mkdir(parents=True,exist_ok=True)
    previous=json.loads((REVIEW/'validation.json').read_text()) if reuse_imports and (REVIEW/'validation.json').exists() else {}
    verified={r['name']:r for r in previous.get('source_files',[])} if previous.get('passed') else {}
    report={'source_files':[],'native_comparison':{},'checks':[],'warnings':[]}
    def check(name,condition):
        report['checks'].append({'name':name,'passed':bool(condition)})
    index=json.loads((SNAPSHOT/'index.json').read_text())
    archives=json.loads((SNAPSHOT/'archives.json').read_text())
    check('26 native archives succeeded',sum(r.get('archive_ok') is True for r in archives['files'])==26)
    check('native source documents were not modified',all(not r.get('modified_before',False) and not r.get('modified_after',False) for r in archives['files']))
    for rec in index['files']:
        folder=SNAPSHOT/rec['slug']
        item={'name':rec['name'],'version':rec['version'],'status':rec['status'],'files':{}}
        for suffix in ['step','f3d']:
            path=folder/f'source.{suffix}'
            if path.exists():item['files'][suffix]={'sha256':digest(path),'bytes':path.stat().st_size}
        cached=verified.get(rec['name'],{})
        if rec.get('step_ok') and cached.get('files')==item['files'] and 'solids' in cached:
            item.update({k:cached[k] for k in ['solids','shells','faces','valid_solids']})
        elif rec.get('step_ok'):
            shape=import_step(folder/'source.step')
            item.update(solids=len(shape.solids()),shells=len(shape.shells()),faces=len(shape.faces()),valid_solids=sum(s.is_valid for s in shape.solids()))
        report['source_files'].append(item)
        (REVIEW/'validation.json').write_text(json.dumps(report,indent=2))
    for key,(slug,_) in STRUCTURES.items():
        bodymap=json.loads((SNAPSHOT/slug/'body-map.json').read_text())
        highfile=SNAPSHOT/slug/'high-accuracy.json'
        high={(r['path'],r['name']):r for r in json.loads(highfile.read_text())} if highfile.exists() else {}
        differences=[];all_differences=[]
        for r in bodymap['matched']:
            native=r['source'];vol=high.get((native['path'],native['name']),native)['volume_mm3']
            delta=(r['volume']-vol)/vol*100 if vol else None
            if delta is not None:all_differences.append(abs(delta))
            if delta is not None and abs(delta)>.25:
                differences.append({'path':native['path'],'name':native['name'],'native_volume_mm3':vol,'step_volume_mm3':r['volume'],'difference_percent':delta})
        report['native_comparison'][key]={'matched_solids':len(bodymap['matched']),'max_adaptive_volume_difference_percent':max(all_differences),'max_bounds_error_mm':max(r['bounds_error'] for r in bodymap['matched']),'missing_native':bodymap['missing_native'],'volume_differences_over_0_25_percent':differences}
        check(f'{key}: adaptive volume within 0.1 percent',max(all_differences)<.1)
        check(f'{key}: matched body bounds within 0.02 mm',all(r['bounds_error']<=.02 for r in bodymap['matched']))
        check(f'{key}: no missing nonzero-volume native solid',all(abs(r.get('volume_mm3',0))<1e-6 for r in bodymap['missing_native']))
    aircraft=build_aircraft();report['assembly']=summary(aircraft)
    check('assembly contains 301 solids',report['assembly']['solids']==301)
    check('all 301 assembly solids pass OCP validity',report['assembly']['valid_solids']==301)
    check('3950 mm span',abs(report['assembly']['size_mm'][1]-3950)<.001)
    mirrored=next(c for c in aircraft.children if c.label=='inferred_starboard_outer_wing')
    original=[s for s in aircraft.children[1].solids() if s.bounding_box().max.Y<=-669.9]
    check('47 mirrored wing solids only',len(mirrored.solids())==len(original)==47)
    check('mirror preserves summed volume',abs(sum(volume_mm3(s) for s in mirrored.solids())-sum(volume_mm3(s) for s in original))<.01)
    check('mirrored outer wing starts at Y=670 mm',abs(mirrored.bounding_box().min.Y-670)<.01)
    restored=aircraft.children[-1]
    check('three hidden mount solids restored',len(restored.solids())==3)
    output=HERE/'output';output.mkdir(exist_ok=True)
    export_step(aircraft,output/'spearhead.step')
    (output/'spearhead.json').write_text(json.dumps(report['assembly'],indent=2))
    roundtrip=import_step(output/'spearhead.step');rt=summary(roundtrip);report['roundtrip']=rt
    check('STEP roundtrip retains 301 solids',rt['solids']==301)
    check('STEP roundtrip solids remain valid',rt['valid_solids']==301)
    error=max(abs(a-b) for aa,bb in zip(rt['bounds_mm'],report['assembly']['bounds_mm']) for a,b in zip(aa,bb))
    check('STEP roundtrip bounds within 0.02 mm',error<.02)
    volume_error=abs(rt['summed_solid_volume_mm3']/report['assembly']['summed_solid_volume_mm3']-1)
    check('STEP roundtrip summed volume within 0.01 percent',volume_error<.0001)
    report['warnings']=[
        'Default non-adaptive OCP volume is inaccurate on some spline surfaces; use spearhead.volume_mm3. Adaptive volumes agree with native high-accuracy measurements. Validity, bounds and volume checks do not prove physical fit.',
        'One zero-volume tail placeholder omitted by STEP export.',
        'One outer wing is inferred by reflection. Restored fourth mount was hidden in source.',
        'Not as-built: current propulsion, avionics, battery, harness, coverings and fastening not reconciled.'
    ]
    report['passed']=all(r['passed'] for r in report['checks'])
    (REVIEW/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':len(report['checks']),'failed':[r for r in report['checks'] if not r['passed']],'assembly':report['assembly']},indent=2))
    return report['passed']

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--reuse-imports',action='store_true',help='Reuse prior passing import results only when SHA256 hashes match')
    raise SystemExit(0 if run(parser.parse_args().reuse_imports) else 1)
