# Control Surface Case Matrix

Case matrix of the control surface deflection study, from [Updated Aerodynamic Database](../information-note.md). The case inputs and loads are in [deflection.csv](deflection.csv) and [deflection_cases_104_113.csv](deflection_cases_104_113.csv).

All cases are at beta = 0 deg unless stated. Neutral means all controls at 0 deg.

| Velocity | Group | Cases | Alpha, deg | Deflections, deg |
|----------|-------|------:|------------|------------------|
| 25 m/s | Neutral | 15 | -5, 0, 5, 10, 15, 18, 20, 25, 30; beta +/-5 at alpha 0, 10, 20 | - |
| 25 m/s | Aileron | 22 | -5, 0, 5, 10, 15, 20, 30 (+/-5); 0 (+/-15, +/-30); 15, 25 (+/-30) | +/-5, +/-15, +/-30 |
| 25 m/s | Elevator | 24 | -5, 0, 5, 10, 15, 20, 30 (+/-5); 0, 10, 18, 25, 30 (+/-15) | +/-5, +/-15 |
| 25 m/s | Rudder | 20 | -5, 0, 5, 10, 15, 20, 30 (+/-5); 0, 20, 25 (+/-15) | +/-5, +/-15 |
| 25 m/s | Elevator and rudder | 4 | 0 (elevator -10, rudder +/-15); 18 (elevator -15, rudder +/-15) | combined |
| 22 m/s | Neutral | 3 | 0, 2, 4 | - |
| 22 m/s | Elevator | 7 | 0, 2, 4 (-20, -30); 0 (-15) | -15, -20, -30 |
| 15 m/s | Neutral and elevator | 9 | 0, 10, 20 | 0, +/-5 |
| 10 m/s | Neutral and elevator | 9 | 0, 10, 20 | 0, +/-5 |
| | **Total** | **113** | | |

| Batch | Cases | Content |
|-------|-------|---------|
| 1 (October 6 to 9) | Case001 to Case073 | Neutral references and +/-5 deg single-surface deflections at 25 m/s from alpha -5 deg to 30 deg; neutral sideslip references at beta +/-5 deg; elevator +/-5 deg at 10 m/s and 15 m/s |
| 2 (October 11) | Case074 to Case103 | Full command travel: elevator and rudder +/-15 deg, aileron +/-15 deg and +/-30 deg; neutral references at alpha 18 deg and 25 deg; combined elevator and rudder commands |
| 3 (October 15) | Case104 to Case113 | Elevator -20 deg and -30 deg at 22 m/s and alpha 0, 2 and 4 deg, with matching neutral references and one -15 deg case |
