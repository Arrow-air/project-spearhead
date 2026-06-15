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

However, the corresponding stability-analysis source module or script is not tracked in the
current tree. That means the published stability tables, figures, and reports cannot be
regenerated from source without reconstructing the missing workflow.

TODO: restore stability-analysis source code as a follow-up PR and make it call the shared
simulation backend (`SimulationConfig`, `run_simulation`, trim, force/moment, and dynamics)
instead of introducing a separate model path.
