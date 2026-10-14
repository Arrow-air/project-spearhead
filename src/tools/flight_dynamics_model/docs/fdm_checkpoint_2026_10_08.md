# FDM checkpoint — 8 October 2026

Historical checkpoint: the active October moment reference was numerically corrected on 12 October. Earlier trim/flight results here do not validate the corrected model; see the [current validation summary](validation_summary.md).


Dated validation reports and recordings referenced below are retained locally and excluded from Git.

This checkpoint captures the aircraft-profile and October AeroDB work before
adding an ArduPilot SITL connection. It is an implementation and analysis
checkpoint, not a validated flight baseline.

## Included

- Independent reference-design and PT1 KDE hardware profiles, plus an
  independent aerodynamic dataset selector in the API, CLI and Flight Deck.
- KDE5215XF-220 / KDE-UAS75HVC+ / 21.5×7.3 lift configuration and static
  manufacturer performance data, with provenance and the original source file.
- October base database and 17 measured deflection cases; original input files
  are preserved. Runtime density is 1.021 kg/m³ and the source moment reference
  follows the user's confirmation of (0,0,0) m.
- Measured control increments, coverage diagnostics and bounded trim rejection
  when the force/moment balance cannot be achieved.
- Reproducible moment-comparison, Nondimit recentering and CG-search scripts.
- CG search, full-equilibrium checks and small-disturbance results (local report).

The latest CG calculations use only `tail_update_2026_10_07`, including its
deflection data. Original aerodynamics remain available for compatibility;
they are not mixed into the October calculations.

## Interpretation of the current results

For the preliminary 25 kg reference-design pusher configuration, physical
CG ≈ (0.01360,0,0.45) m gives locally stable full equilibria at the sampled
22,25,27,28,30,32,35 m/s conditions. The 25 and 28 m/s cases were recomputed
with the new dataset on 8 October. CG ≈ (0.01354,0,0.35) m also gives stable
sampled equilibria, with less elevator margin. These CG locations have not
been checked against the physical airframe's feasible mass distribution.

The original raw base CSV header specifies (-0.1,0,0.15) m, conflicting with
the user-confirmed origin. Treating that raw header as the actual export center
and transferring the loads to the origin reduces the pitch-moment RMS mismatch
against the three zero-control deflection cases from 30.77 to 0.99 N·m.
That is a diagnostic hypothesis requiring export-setting confirmation; it
does not silently change the selected runtime reference.

Rate derivatives, inertia, rotor locations and the pusher model remain
preliminary. The pusher currently acts through the selected CG. Updated
numeric tail chord and area are pending. Control coverage is partial and
beta measurements are limited to zero. PT1 currently has no installed pusher.

## Generated outputs

The source databases are versioned. CG assessments, per-center Nondimit
database copies and ZIP files are reproducible local outputs, excluded from Git:

```sh
PYTHONPATH=src/tools/flight_dynamics_model python \
  src/tools/flight_dynamics_model/scripts/export_nondimit_moment_centers.py \
  --output-dir artifacts/aerodb_tail_update_2026_10_07_moment_centers
```

Their manifests, references and validation results are retained locally. Firmware packages, unrelated flight-review data, Python
bytecode and Finder metadata are not part of this change.

## Next implementation boundary

Add an ArduPilot external-physics adapter around the shared dynamics engine:
PWM channel mapping, synchronized stepping and body-IMU/NED state messages.
Use the October dataset explicitly. Verify the protocol and sensor conventions
before running closed-loop fixed-wing tests, then add PT1 hover tests.
