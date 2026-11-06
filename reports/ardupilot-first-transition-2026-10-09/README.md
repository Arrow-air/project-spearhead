# ArduPilot first transition review — 9 October 2026

Start with **[report.md](report.md)**. It contains the full assessment, all 24 proposed changes and reasons, comparisons with the aircraft log and FDM/SITL, similar-aircraft references, unresolved measurements, and the step-by-step validation plan.

This is a saved review package. It does not change aircraft or simulator parameters.

**Later data update:** [control-data-update.md](control-data-update.md) records the integration and validation of the aero team's subsequently supplied `results.csv`. The original evidence files remain a snapshot from before that update.

## Recommendations and comparison

| File | Contents |
|---|---|
| [report.md](report.md) | Complete report and validation sequence |
| [proposed-changes.md](proposed-changes.md) | All 24 changes, current values and individual conditions |
| [parameter-comparison.csv](parameter-comparison.csv) | 291 rows comparing aircraft log, actual SITL, requested SITL overlay, proposal and peer values |
| [candidate.json](candidate.json) | Proposed values, rationale, retained settings and pending measurements |
| [candidate-review.param.txt](candidate-review.param.txt) | Commented parameter draft for review; all parameter lines are disabled |
| [geometry-calculations.json](geometry-calculations.json) | Wing loading, speed conversion, stall sensitivity and turning-radius calculations |

## Evidence from the supplied aircraft log

| File | Contents |
|---|---|
| [log17-first-parameters.json](log17-first-parameters.json) | First observed values in 00000017.BIN |
| [log17-last-parameters.json](log17-last-parameters.json) | Last observed values in that log |
| [log17-parameter-changes.json](log17-parameter-changes.json) | Parameter changes recorded during the log |
| [log17-last-baseline.param](log17-last-baseline.param) | Historical parameter snapshot; not an aircraft deployment overlay |
| [log17-summary.json](log17-summary.json) | Source hash, firmware messages, modes, arming events and whole-log statistics |
| [log17-last-hover-window.json](log17-last-hover-window.json) | Bounded hover-window tracking, motor, battery and vibration statistics |

## FDM and SITL evidence

| File | Contents |
|---|---|
| [sitl-last-parameters.json](sitl-last-parameters.json) | Actual stored parameters extracted from the 4.7.1 simulation log |
| [sitl-log-identity.json](sitl-log-identity.json) | Simulation-log source path, hash and firmware identity |
| [fdm-evidence.json](fdm-evidence.json) | Matching aircraft assumptions, trim, requested overlay and limitations |

## Sources and reproducibility

| File | Contents |
|---|---|
| [sources.json](sources.json) | Official documentation links and pinned peer-file provenance |
| [Foxtech_GreatShark.source.param](Foxtech_GreatShark.source.param) | Upstream GreatShark parameter reference |
| [Mugin_EV350.source.param](Mugin_EV350.source.param) | Upstream EV350 parameter reference |
| [extract_log.py](extract_log.py) | Original log-extraction script |
| [build_review.py](build_review.py) | Original comparison/calculation builder |
| [SHA256SUMS](SHA256SUMS) | Integrity hashes for every report and supporting file |

Raw binary logs and the ArduPilot source checkout are referenced by path and hash in the evidence files; they are not duplicated in this folder. The scripts require those original local inputs and the documented Python dependencies. The saved reports, CSV and JSON can be reviewed independently.
