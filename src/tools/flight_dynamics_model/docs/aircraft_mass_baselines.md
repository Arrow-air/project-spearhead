# Aircraft mass baselines

User-supplied update on 10 October 2026:

| Configuration | Mass | Use |
|---|---:|---|
| PT1.5 prototype (`pt15`) | 21 kg | KDE prototype with pusher; FDM, SITL and Flight Deck |
| PT1 as flown, log17 wings on (`pt1_as_built`) | 17.9 kg | Log comparison; hover thrust on the KDE curve gives 17.8 kg |
| Planned finished aircraft (`reference_design`) | 25 kg | Original planned configuration and reference fixed-wing model |

Changing prototype mass does not determine new CG or inertia. Rotor coordinates follow the 14 October drawing; CG and inertia are recorded per profile in [aircraft profiles](aircraft_profiles.md). The static hover scenario and SITL hover runner calculate weight support from the selected profile and voltage each time they start; the normalized PWM collective and firmware hover-thrust setting are not fixed constants.

The 9 October hover recording and the initial log comparison used a 25 kg PT1 assumption. Their saved manifests and numeric results remain historical evidence. New PT1.5 runs use 21 kg. See the mass-update and comparison deliverables under `artifacts/` for the revised interpretation and rerun.
