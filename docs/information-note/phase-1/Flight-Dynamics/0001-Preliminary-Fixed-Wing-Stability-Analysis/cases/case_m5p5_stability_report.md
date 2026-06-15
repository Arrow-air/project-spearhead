# Case M5p5 Stability Report

> Preliminary FDM result: inertia, control derivatives, rate damping, and propulsion are placeholders; ADB v1.1 contains alpha-beta data only.

## Database Metadata

- File: `adb_v1_1_drag_correction_m5p5.csv`
- Moment/reference/CG: `(-0.055, 0.000, 0.150) m`
- rho: `1.090 kg/m^3`, q_ref: `340.625 Pa`, S: `1.805 m^2`, b: `3.950 m`, c: `0.457 m`

## Trim

- success: `True`
- alpha: `1.9787 deg`
- theta: `1.9787 deg`
- elevator: `-5.1878 deg`
- throttle: `0.1659`
- residual norm: `2.463e-14`

## Static Stability Derivatives

- Cm_alpha: `-0.214751 /rad` (stable longitudinal static sign)
- Cn_beta: `+0.028648 /rad` (stable directional static sign)
- Cl_beta: `-0.011398 /rad` (stable lateral static sign)
- CL_alpha: `+4.924211 /rad`
- CY_beta: `-0.226440 /rad`

## Dynamic Modes

| Mode | Eigenvalues | wn rad/s | zeta | time constant s | Classification |
|---|---:|---:|---:|---:|---|
| phugoid | `-0.04536+0.166775j; -0.04536-0.166775j` | 0.1728 | 0.2624 | 22.0459 | stable |
| short_period | `-4.37067+1.2534j; -4.37067-1.2534j` | 4.5468 | 0.9613 | 0.2288 | stable |
| dutch-roll-like lateral-directional | `No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization.` | n/a | n/a | n/a | not identified |
| roll | `-30.2268+0j` | 30.2268 | 1.0000 | 0.0331 | stable |
| spiral | `-0.00305176+0j` | 0.0031 | 1.0000 | 327.6796 | stable |

## Summary

- Static stability: longitudinal `stable`, directional `stable`, lateral `stable`.
- Dynamic stability among identified reduced-order modes: `stable`.
- No clearly identifiable Dutch-roll-like lateral-directional oscillatory mode was found in this linearization. Because the aircraft uses an inverse V-tail/ruddervator layout and the current lateral reduced eigenspectrum is entirely real, treat this as a model limitation and classification ambiguity to investigate before flight-test use.

## Conclusion

Case M5p5 is the first sampled case with stable static pitch behavior and stable identified longitudinal modes. It should be treated as the current sampled aft boundary of the stable region rather than as a high-margin baseline.
