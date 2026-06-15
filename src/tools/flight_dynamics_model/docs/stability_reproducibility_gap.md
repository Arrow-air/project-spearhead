# Stability Reproducibility Gap

The repository currently tracks preliminary stability-analysis outputs under:

```text
src/tools/flight_dynamics_model/output/stability/
docs/information-note/phase-1/Flight-Dynamics/0001-Preliminary-Fixed-Wing-Stability-Analysis/
```

It also tracks:

```text
src/tools/flight_dynamics_model/spearhead/__pycache__/stability.cpython-310.pyc
```

In `codex/sim-core-platform`, the corresponding stability-analysis source module or script
was not tracked. That meant the published stability tables, figures, and reports could not
be regenerated from source without reconstructing the missing workflow.

This stacked branch restores a source implementation under:

```text
src/tools/flight_dynamics_model/spearhead/stability/
```

The restored implementation calls the shared simulation backend (`SimulationConfig`,
`run_simulation`, trim, force/moment, and dynamics) instead of introducing a separate
model path.

Remaining limitation: the tracked aerodynamic database is the nominal ADB only. The CG
sweep implementation varies runtime `cg_body_xyz` through the existing moment-shift path
and does not recreate missing historical Nondimit per-CG CSV cases.
