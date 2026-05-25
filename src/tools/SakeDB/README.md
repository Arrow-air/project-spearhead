# SakeDB v1.0 Quick Start

SakeDB is a desktop workflow tool for building aerodynamic databases from high-fidelity analysis cases. It manages the sampling plan, prepares each case folder, launches a user-defined external runner, tracks progress, and exports dense force and moment tables for flight dynamics and control work.

## Purpose

Aerodynamic databases need broad coverage across angle of attack, sideslip, and later control-surface deflections. Running every possible point directly is expensive. SakeDB reduces this cost by combining anchor cases, Latin-hypercube sampling, adaptive infill, and surrogate-surface export into one traceable workflow.

## What It Does

| Function | Role |
|----------|------|
| Study setup | Stores variables, bounds, outputs, budgets, anchors, symmetry settings, and runner settings in a **.sakedb** file |
| Case management | Creates one folder per analysis case and records artifacts, status, timing, and output paths |
| Sampling | Runs anchors first, then initial DOE points, then adaptive infill based on uncertainty, curvature, and spacing |
| Pause and resume | Stops or resumes studies without losing completed results |
| Discard control | Lets a bad current or completed case be excluded from surrogate fitting |
| Export | Writes raw results, symmetric results, and dense surrogate tables for downstream coefficient conversion |

## Runner Model

SakeDB is intentionally runner-neutral. A study defines an external runner, a model or mesh file, a setup script, and a launcher script. The launcher is responsible for executing the high-fidelity analysis and writing the expected output files. SakeDB then reads the requested force and moment columns and updates the study database.

The default template does not prefill a solver product or local program path. Each project supplies its own runner settings inside the study folder.

| Setting | Meaning |
|---------|---------|
| **sim_file** | The base model, mesh, or case file used by the external analysis. Keep it in the study folder for simple relative paths. SakeDB passes its path to the launcher as **SAKEDB_SIM_FILE** and **AERODB_SIM_FILE**. The launcher opens it directly and runs it with the setup file in the case folder |
| **setup_file** | A setup script or template copied into each case folder. SakeDB patches its **alpha** and **beta** definitions for the current sample point before the launcher runs |
| **run_script** | A launcher script copied into each case folder and executed from that folder. It should run the external analysis, write **output.out**, write **output.dat**, and optionally write **force.png** |

The setup script must expose angle definitions in this form so SakeDB can replace the contents for each case:

**double[] alpha = { 0.0 };**

**double[] beta = { 0.0 };**

The important part is the **alpha = { ... };** and **beta = { ... };** pattern. During a run, SakeDB rewrites those values with the current sample, for example **alpha = { -5.0 };** and **beta = { 10.0 };**.

For each sample, SakeDB creates a case folder such as **runs/cases/0007_infill_a-05.000_b+10.000**. That folder contains the patched setup script and copied launcher script. The launcher runs from inside this folder, using the environment variables supplied by SakeDB, and produces the case artifacts there. If a local model or mesh copy is required, the launcher should copy **SAKEDB_SIM_FILE** into the case folder before starting the analysis.

Required and expected case artifacts:

| File | Requirement |
|------|-------------|
| **output.out** | Text output from the external analysis. This is used for case inspection and troubleshooting |
| **output.dat** | Required force and moment table. It is also the default completion file, so the case is not accepted until it exists |
| **force.png** | Optional force-history image. If produced, it should be named exactly **force.png** so the Force Plot tab can find it |

**output.dat** must contain a **VARIABLES** header followed by numeric rows. The recommended format is:

**VARIABLES = "Beta" "Alpha" "Fx" "Fy" "Fz" "Mx" "My" "Mz"**

**10.0 -5.0 -12.3456 4.5678 -220.1234 0.1000 95.2000 -0.0300**

The force and moment columns **Fx**, **Fy**, **Fz**, **Mx**, **My**, and **Mz** are required. Column names are read from the quoted **VARIABLES** header, and numeric values may use standard decimal or scientific notation. A single final row is sufficient. If several complete rows are written, SakeDB uses the final rows according to the driver averaging setting.

## How To Use

Keep **SakeDB.exe**, **_internal**, **assets**, and **template.sakedb** together in the same folder, then launch **SakeDB.exe**. Use **New** to create a study from the bundled template or **Open** to load an existing **.sakedb** file. Put the model file, setup script, and launcher script beside the study file when you want simple relative paths.

You can also open a study directly from Windows by choosing **Open with** and selecting **SakeDB.exe** for a **.sakedb** file. Windows passes the study path to SakeDB, and SakeDB loads that study at startup. For this to work on another computer, copy the whole SakeDB folder, not only **SakeDB.exe**.

1. Define the study in **Setup**: set the study name, run folder, run budget, convergence stopping options, model file, setup script, launcher script, variables, export grid, and anchor points.
2. Save the study file, then press **Start**. SakeDB runs anchors first, then initial sampling, then adaptive infill.
3. Use **Pause**, **Resume**, and **Stop** from the sidebar to control a campaign. Use **Discard Current Case** when the active case is clearly bad.
4. Review sampled cases in **Samples**. Completed bad cases can be excluded with **Discard Selected Case**.
5. Press **Write Exports** when enough samples are available. SakeDB writes raw, symmetric, and dense surrogate outputs into the run folder.

| Area | Use |
|------|-----|
| Sidebar and status | Select study and run folders, start, pause, resume, stop runs, discard the current case, write exports, refresh the display, switch dark mode, and monitor state, sample count, budget, phase, current case, and stop reason |
| Setup tab | Edit the main study settings without touching YAML directly: budget, convergence threshold, model, setup, launcher files, variable bounds, export grid, anchor values, and extra anchor points |
| Samples tab | View the sample map, live current-case output, and the full results table. The table shows all samples while hiding long artifact paths and internal bookkeeping columns |
| Exports tab | Inspect generated dense export files and preview the dense force or moment table written for downstream flight-dynamics or coefficient-table work |
| Surrogate Surface tab | Select an output component and refresh a fitted surrogate surface to check trends, smoothness, and coverage |
| Force Plot tab | Select a completed case and inspect its per-case force plot when the runner produced one |
| Candidate Scores tab | Show the top current acquisition candidates from the fitted surrogate, including score, uncertainty, curvature, region, spacing, and distance terms |
| Study File tab | Edit or inspect the raw **.sakedb** YAML when an advanced setting is not exposed in the Setup tab |

## Why It Matters

SakeDB turns a large manual simulation campaign into a repeatable database process. The main benefit is not only fewer wasted cases, but cleaner traceability: every sampled point, discarded case, output artifact, symmetry operation, and dense export can be tied back to the study file and run folder.

## Near-Term Next Steps

- Add control-surface deflection dimensions or derivative tables.
- Keep the runner contract generic so future projects can connect different high-fidelity tools without changing the SakeDB workflow.
