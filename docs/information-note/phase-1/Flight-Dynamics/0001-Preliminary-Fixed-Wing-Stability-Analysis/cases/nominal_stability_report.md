# Nominal Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction.csv`
- Moment/reference/CG: `(-0.150, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.5005 deg`
- theta: `1.5005 deg`
- elevator: `1.5157 deg`
- throttle: `0.1563`
- residual norm: `1.279e-13`

## Static Stability Derivatives

- Cm_alpha: `+0.985487 /rad` (unstable longitudinal static sign)
- Cn_beta: `+0.024349 /rad` (stable directional static sign)
- Cl_beta: `-0.010028 /rad` (stable lateral static sign)
- CL_alpha: `+4.858682 /rad`
- CY_beta: `-0.229180 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.091489+0.764909j; -0.091489-0.764909j` | 0.7704 | 0.1188 | 10.9303 | stable |
| short_period | `-11.7041+0j; +3.06093+0j` | 3.0609 | -1.0000 | 0.3267 | unstable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.229+0j` | 30.2290 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00898588+0j` | 0.0090 | 1.0000 | 111.2858 | stable |

## Summary

- Static stability: longitudinal `unstable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `unstable`.
- Weakly damped modes: phugoid.
- Unstable modes: short_period.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

The nominal database is not an acceptable fixed-wing stability baseline in the current FDM workflow. It is statically unstable in pitch and has an unstable longitudinal dynamic mode, although it remains useful as the aft/reference starting point for the CG sensitivity sweep.
