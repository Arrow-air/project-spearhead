# FDM and Flight Deck checkpoint — 11 October 2026

Historical checkpoint: the active October moment reference was numerically corrected on 12 October. Earlier trim/flight results here do not validate the corrected model; see the [current validation summary](validation_summary.md).


Dated validation reports and recordings referenced below are retained locally and excluded from Git.

This checkpoint contains the October control-data update, independent PT1/Spearhead
session setup, QuadPlane mode commands and firmware warnings, the FDM horizon and
3D observer, PT1 propulsion estimates, and the ground-start takeoff/hover/land mission.

PT1 uses 21 kg, four KDE lift systems and the MAD V8013 150 KV / AMPX 80 A /
APC 20×12 pusher estimate. The preserved reference design uses 25 kg; the logged
aircraft comparison assumes approximately 18 kg. New SITL sessions use the October
AeroDB, rho 1.021 and the specified source moment origin (0,0,0).

Validation before committing: **299 FDM tests passed** with Python 3.10 on macOS.
The previously recorded real ArduPlane 4.7.1 ground missions at 5, 10 and 20 m
normally armed, took off, hovered, landed on the modeled gear and disarmed.
See the ground-start report (local report)
and pusher estimate (local report).
The simulator integration is validated within those scenarios; physical aircraft
behavior and forward transition remain unvalidated. Gear, CG/inertia, control
travel and propulsion predictions retain their documented assumptions.

The source, configuration, input data, reports and compact validation evidence
are versioned. Local SITL binaries, build environments and raw ground/pusher
run logs are excluded. Historical local artifacts outside this checkpoint remain
untouched. To run firmware-backed sessions on another machine, follow the
[native SITL setup](ardupilot_sitl_quickstart.md).

The browser is the interface; Python FDM and ArduPlane run on the GUI server's
machine. A server manages one active session. `/embed` supports website embedding;
Docker packaging, hosted independent sessions and forward-transition acceptance
are follow-up work.
