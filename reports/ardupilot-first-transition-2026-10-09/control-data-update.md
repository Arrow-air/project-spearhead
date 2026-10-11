# Later control-data update — 9 October 2026

After the parameter-review snapshot was saved, the aero team supplied a valid `results.csv`. Its expanded controls are now integrated. The earlier report's discussion of waiting for corrected control data is historical; its saved FDM/SITL evidence and parameter comparisons still describe the pre-update snapshot.

Read the [integration and validation assessment](../../artifacts/control_results_2026_10_09/assessment.md) for the full intake audit, validation outputs and limitations. There are now 49 active control cases covering alpha −5…30 degrees for all three surfaces, with measured deflections still limited to ±5 degrees. All 17 overlapping original cases are unchanged. Five zero-control sideslip cases are retained for audit.

The full regression suite passed 210 tests and the model reproduced 245 source-increment comparisons at 18–28 m/s. The existing diagnostic setup now finds positive-alpha trims at 22 and 25 m/s without surface extrapolation; 18 and 20 m/s need elevator extrapolation. These are FDM checks with a diagnostic CG and provisional aircraft assumptions, not a PT1 transition validation.

No ArduPilot parameters were changed by this integration. The original candidate remains a review draft, with its airspeed, geometry, propulsion and physical-airframe conditions still applicable. CFD convergence is marked `not_assessed` in the new source. The historical evidence JSON and earlier simulation logs were preserved rather than overwritten.
