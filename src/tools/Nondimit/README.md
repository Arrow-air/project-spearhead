# Nondimit One Pager

## Purpose

Nondimit is a standalone desktop tool for turning dense aerodynamic force and moment tables into coefficient database files. It is intended for AeroDB and CFD export workflows where the source file contains dimensional forces and moments over alpha and beta, and the final file needs nondimensional body and wind-frame coefficients.

The tool can also reopen an existing Nondimit export and modify it without rerunning CFD. This can be used to move the moment reference center, change the nondimensionalization reference values, or do both in the same saved file.

![Nondimit reference frame notation](aero_frame_notation.png)

The illustration shows the force and moment directions. The table below defines the body and wind-frame notation used by the tool.

| Frame axis | Force name | Moment name | Force coefficient | Moment coefficient | Positive direction |
| --- | --- | --- | --- | --- | --- |
| Body `+x` | `Fx` | `Mx` | `cfx` | `cmx` | Forward toward the nose. |
| Body `+y` | `Fy` | `My` | `cfy` | `cmy` | Toward the aircraft right wing. |
| Body `+z` | `Fz` | `Mz` | `cfz` | `cmz` | Down. |
| Wind `+x` | `D` | `R` | `cd` | `cr` | Aft along the drag direction, opposite the velocity direction. |
| Wind `+y` | `S` | `M` | `cs` | `cm` | Side direction. At zero beta this matches body `+y`. |
| Wind `+z` | `L` | `N` | `cl` | `cn` | Up along the lift direction. |

## Launch

Run:

```bat
src\tools\Nondimit\Start_Nondimit.bat
```

The launcher starts `Nondimit.exe` when present. If the EXE is missing, it falls back to `nondimit.py`.

## Main Workflows

### Convert Raw Table

Use this tab when starting from a raw dense CSV with dimensional force and moment columns.

1. Select the input CSV.
2. Choose the output path.
3. Enter reference values:
   - density `rho`
   - reference area `S`
   - velocity `V`
   - span for roll and yaw moments
   - chord for pitch moment
   - output rounding accuracy
4. Enter old and new moment centers in body-axis coordinates.
5. Confirm column mapping for alpha, beta, `Fx`, `Fy`, `Fz`, and optionally `Mx`, `My`, `Mz`.
6. Set optional beta-zero cleanup and output group order.
7. Click `Convert and Save`.

If moments are mapped, the tool writes force and moment coefficients. If moment columns are not mapped, it writes force coefficients only.

### Modify Export

Use this tab when you already have a coefficient CSV from Nondimit and need to change the moment center, the nondimensionalization reference values, or the output formatting.

1. Select the existing Nondimit export.
2. The tool reads the metadata row and fills in `rho`, `S`, `V`, span, chord, output accuracy, and current moment center.
3. Edit the reference values if the output should use a different nondimensionalization.
4. Enter the new body-axis moment center if the moment reference should move. Leave it equal to the current center if only the reference values should change.
5. Confirm alpha and beta columns.
6. Set optional beta-zero cleanup and output group order.
7. Click `Modify and Save`.

The current center is read from `new_moment_center_body_xyz` in the input file metadata. The new output metadata records the previous center as `old_moment_center_body_xyz` and the selected center as `new_moment_center_body_xyz`.

The reference values in this tab are intentionally editable. When `rho`, `S`, `V`, span, or chord are changed, the tool recalculates all force and moment coefficients with the edited values:

```text
q = 0.5 rho V^2
force coefficients = force / qS
Mx, Mz coefficients = Mx,Mz / qS span
My coefficient = My / qS chord
```

Normal Nondimit exports include dimensional body force and moment columns, so the tool preserves the physical forces and moments, applies any requested moment-center shift, and then writes new coefficients from the edited reference values. If an export has been manually stripped down to coefficient-only columns, do not use this tab to change reference values, because the physical forces and moments cannot be recovered unambiguously.

## Optional Output Controls

### Beta = 0 Cleanup

The checkboxes `Fy`, `Mx`, and `Mz` force those components to zero only on rows where beta is numerically zero. This is useful for symmetric prop-off databases where small solver noise should not create side force, roll moment, or yaw moment at zero sideslip.

The cleanup is applied before coefficients and wind-frame values are written, so body and wind outputs remain internally consistent. Selected cleanup options are stored in metadata as `zero_at_beta0`.

### Output Group Order

The ordered list controls how the generated output columns are arranged. The four groups are:

- Body Forces: `Fx`, `Fy`, `Fz`
- Body Moments: `Mx`, `My`, `Mz`
- Wind Forces: `D`, `S`, `L`
- Wind Moments: `R`, `M`, `N`

Use `Up`, `Down`, and `Reset` to reorder them. Each group keeps its coefficient columns and dimensional columns together. The selected order is stored in metadata as `output_group_order`.

## Coordinate and Coefficient Conventions

The tool uses standard aircraft body axes and aerodynamic wind-axis notation.

Body frame directions are `x` forward, `y` right, and `z` downward. Wind-frame directions are `D` aft along the drag direction, `S` sideways, and `L` upward. The moment symbols `R`, `M`, and `N` are the wind-frame roll, pitch, and yaw moments about those same wind-frame axes.

The generated body coefficient columns are:

```text
cfx, cfy, cfz
cmx, cmy, cmz
```

The generated wind coefficient columns are:

```text
cd, cs, cl
cr, cm, cn
```

The generated dimensional body and wind columns are:

```text
Fx, Fy, Fz
Mx, My, Mz
D, S, L
R, M, N
```

Forces are nondimensionalized by:

```text
q S
```

Moments are nondimensionalized by:

```text
Mx, Mz: q S span
My:     q S chord
```

where:

```text
q = 0.5 rho V^2
```

## Moment Center Operation

Moment shifting is done in body axes before coefficient conversion:

```text
M_new = M_old + r_new_to_old x F
```

where:

```text
r_new_to_old = old_center - new_center
```

After the body-axis moment is shifted, the tool projects both force and moment vectors into the wind frame using the row's alpha and beta.

## Output File

The saved CSV starts with a metadata row, then the normal CSV header row, then table data. Generated columns include body and wind-frame coefficients plus dimensional body and wind-frame force and moment columns. Numeric table cells are rounded to the requested output accuracy.

Recommended sanity checks after export:

- Beta-zero rows have `Fy`, `Mx`, and `Mz` equal to zero if those cleanup options were selected.
- Moment signs change as expected when moving the reference center.
- Metadata values match the intended `rho`, `S`, `V`, span, chord, and moment centers.
