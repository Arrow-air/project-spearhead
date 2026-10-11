"""Read the supplied DataFlash log without changing firmware or parameters."""
import collections, hashlib, json, pathlib
from pymavlink import DFReader
import numpy as np

root = pathlib.Path(__file__).parent
source = pathlib.Path('/Users/zeynep/Downloads/00000017.BIN')
reader = DFReader.DFReader_binary(str(source))
first, last, changes, counts, modes, arms, messages = {}, {}, [], collections.Counter(), [], [], []
samples = collections.defaultdict(list)
last_time = 0
while True:
    msg = reader.recv_msg()
    if msg is None:
        break
    kind = msg.get_type()
    counts[kind] += 1
    data = msg.to_dict()
    if 'TimeUS' in data:
        last_time = max(last_time, data['TimeUS'] / 1e6)
    if kind == 'PARM':
        name, value = str(msg.Name).rstrip('\x00'), float(msg.Value)
        first.setdefault(name, value)
        if name in last and last[name] != value:
            changes.append({'name': name, 'before': last[name], 'after': value, 'time_s': data.get('TimeUS', 0)/1e6})
        last[name] = value
    if kind == 'MODE':
        modes.append(data)
    if kind == 'ARM':
        arms.append(data)
    if kind == 'MSG':
        messages.append(data)
    for field in {'ARSP':['Airspeed','DiffPress','H'], 'BAT':['Volt','Curr'], 'QTUN':['ThH','ThO','Alt','DAlt'], 'VIBE':['VibeZ','Clip0','Clip1','Clip2'], 'RCOU':['C5','C6','C7','C8']}.get(kind, []):
        if field in data:
            samples[kind+'.'+field].append(float(data[field]))
stats = {k: {'n':len(v),'min':min(v),'median':float(np.median(v)), 'p95':float(np.percentile(v,95)), 'max':max(v)} for k,v in samples.items()}
summary = {'source':str(source), 'sha256':hashlib.sha256(source.read_bytes()).hexdigest(), 'duration_boot_s':last_time, 'parameter_count':len(last), 'modes':modes, 'arms':arms, 'messages':messages, 'statistics_whole_log_not_just_flight':stats, 'message_counts':counts}
for filename, value in [('log17-first-parameters.json',first),('log17-last-parameters.json',last),('log17-parameter-changes.json',changes),('log17-summary.json',summary)]:
    (root/filename).write_text(json.dumps(value, indent=2)+'\n')
(root/'log17-last-baseline.param').write_text('# Historical 4.7.0 log snapshot; preserve as evidence, not a hardware overlay.\n'+''.join(f'{k} {v:.9g}\n' for k,v in sorted(last.items())))
print(json.dumps({'parameter_count':len(last), 'duration_boot_s':last_time,'modes':modes,'arms':arms,'statistics':stats,'changed_parameters':len(changes)}, indent=2))
