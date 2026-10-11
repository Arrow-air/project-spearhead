# Aircraft mass baselines

User-supplied update on 10 October 2026:

| Configuration | Mass | Use |
|---|---:|---|
| PT1 prototype (`pt1_as_built`) | 21 kg | Current KDE prototype FDM and QuadPlane hover |
| Aircraft in supplied log17 | Approximately 18 kg | Interpretation of the historical flight; user estimate, not measured by DataFlash |
| Planned finished aircraft (`reference_design`) | 25 kg | Original planned configuration and reference fixed-wing model |

Changing prototype mass does not determine new CG, inertia or rotor coordinates. Those existing provisional values remain pending measurements. The static hover scenario and SITL hover runner calculate weight support from the selected profile and voltage each time they start; the normalized PWM collective and firmware hover-thrust setting are not fixed constants.

The 9 October hover recording and the initial log comparison used a 25 kg PT1 assumption. Their saved manifests and numeric results remain historical evidence. New prototype runs use 21 kg. See the mass-update and comparison deliverables under `artifacts/` for the revised interpretation and rerun.
