# Case M8 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m8.csv`
- Moment/reference/CG: `(-0.080, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.8458 deg`
- theta: `1.8458 deg`
- elevator: `-3.3242 deg`
- throttle: `0.1630`
- residual norm: `3.985e-14`

## Static Stability Derivatives

- Cm_alpha: `+0.155594 /rad` (unstable longitudinal static sign)
- Cn_beta: `+0.028648 /rad` (stable directional static sign)
- Cl_beta: `-0.011018 /rad` (stable lateral static sign)
- CL_alpha: `+4.886141 /rad`
- CY_beta: `-0.227202 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.745163+0.728641j; -0.745163-0.728641j` | 1.0422 | 0.7150 | 1.3420 | stable |
| short_period | `-7.91482+0j; +0.574806+0j` | 0.5748 | -1.0000 | 1.7397 | unstable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2275+0j` | 30.2275 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00116757+0j` | 0.0012 | 1.0000 | 856.4773 | stable |

## Summary

- Static stability: longitudinal `unstable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `unstable`.
- Unstable modes: short_period.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M8 is close to the static pitch-stability transition but has not crossed into acceptable behavior. It remains statically unstable in pitch and dynamically unstable in the identified longitudinal response.
