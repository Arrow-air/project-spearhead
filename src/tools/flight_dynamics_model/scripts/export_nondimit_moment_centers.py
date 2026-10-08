"""Batch the real Nondimit backend at the previous study's moment centers.

The October input uses the user-confirmed origin; its contradictory original
header is preserved separately. Each case starts from the origin input, so
rounding does not accumulate across a chain of recentering operations.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import zipfile

import numpy as np

FDM_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FDM_ROOT))
sys.path.insert(0, str(FDM_ROOT.parent / "Nondimit"))
import nondimit as nd
from spearhead.aerodb import load_aerodynamic_database
from spearhead.deflections import load_deflection_database

SOURCE = FDM_ROOT / "data/aerodynamics/tail_update_2026_10_07"
CASES = {
    "origin": (0.0, 0.0, 0.0),
    "nominal": (-0.150, 0.0, 0.150),
    "m10": (-0.100, 0.0, 0.150),
    "m8": (-0.080, 0.0, 0.150),
    "m65": (-0.065, 0.0, 0.150),
    "m6": (-0.060, 0.0, 0.150),
    "m5p5": (-0.055, 0.0, 0.150),
    "m5": (-0.050, 0.0, 0.150),
    "x0_z015": (0.0, 0.0, 0.150),
}
FORCES = ("Fx", "Fy", "Fz")
MOMENTS = ("Mx", "My", "Mz")


def vectors(rows, columns):
    return np.array([[float(row[column]) for column in columns] for row in rows])


def validate_output(path, source_rows, center, settings, mapping):
    metadata, headers, rows = nd.read_csv_table_with_metadata(path)
    _, exported_settings = nd.read_export_settings(metadata)
    if len(rows) != len(source_rows) or exported_settings.new_center != center:
        raise AssertionError("Exported row count or moment center is incorrect")
    for key in ("rho", "area", "velocity", "span", "chord", "output_accuracy"):
        if getattr(exported_settings, key) != getattr(settings, key):
            raise AssertionError(f"Export changed reference value {key}")
    force = vectors(source_rows, FORCES)
    moments = vectors(source_rows, MOMENTS)
    expected = moments + np.cross(-np.asarray(center), force)
    saved_force = vectors(rows, FORCES)
    saved_moment = vectors(rows, MOMENTS)
    force_error = float(np.max(np.abs(saved_force - force)))
    moment_error = float(np.max(np.abs(saved_moment - expected)))
    tolerance = settings.output_accuracy / 2 + 1e-9
    if force_error > tolerance or moment_error > tolerance:
        raise AssertionError(f"Failed rigid-body transfer check: {force_error=}, {moment_error=}")
    # Validate dimensional/coefficient consistency within exported rounding.
    scales = np.array([settings.force_denominator] * 3 +
                      [settings.roll_yaw_denominator, settings.pitch_denominator, settings.roll_yaw_denominator])
    coeffs = vectors(rows, ("cfx", "cfy", "cfz", "cmx", "cmy", "cmz"))
    dimensional = np.column_stack([saved_force, saved_moment])
    coefficient_error = float(np.max(np.abs(coeffs - dimensional / scales)))
    # Nondimit rounds moments and coefficients independently. Check coefficients
    # against the unrounded transfer; rounded dimensional columns have their
    # own half-step error, which also contributes to displayed consistency.
    source_coefficient_error = float(np.max(np.abs(coeffs - np.column_stack([force, expected]) / scales)))
    if source_coefficient_error > tolerance:
        raise AssertionError("Exported coefficients do not match dimensional loads")
    alpha_column, beta_column = mapping["alpha"], mapping["beta"]
    for source, output in zip(source_rows, rows):
        for column in [alpha_column, beta_column, "case_id", "aileron", "elevator", "rudder", "velocity"]:
            if column in source and source[column].strip() != output[column].strip():
                try:
                    same = float(source[column]) == float(output[column])
                except ValueError:
                    same = False
                if not same:
                    raise AssertionError(f"Export changed case coordinate {column}")
    return {"rows": len(rows), "max_force_error_n": force_error,
            "max_moment_transfer_error_nm": moment_error,
            "max_coefficient_transfer_error": source_coefficient_error,
            "max_coefficient_consistency_error": coefficient_error}, headers, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    nd.self_test()

    metadata, headers, base_rows = nd.read_csv_table_with_metadata(SOURCE / "database.csv")
    _, reference = nd.read_export_settings(metadata)
    if reference.old_center != (0.0, 0.0, 0.0) or reference.rho != 1.021:
        raise ValueError("Source must use the user-confirmed origin and density")
    mapping = nd.default_mapping(headers)
    _, def_headers, def_rows = nd.read_csv_table_with_metadata(SOURCE / "sources/deflection.csv")
    def_headers = [column.strip() for column in def_headers]
    def_rows = [{key.strip(): value.strip() for key, value in row.items()} for row in def_rows]
    def_mapping = nd.default_mapping(def_headers)
    if any(float(row["velocity"]) != reference.velocity for row in def_rows):
        raise ValueError("Deflection velocities must match the export reference")

    validation = []
    provenance = {
        "tool": nd.APP_DISPLAY_NAME,
        "tool_source_sha256": hashlib.sha256(Path(nd.__file__).read_bytes()).hexdigest(),
        "source_database_sha256": hashlib.sha256((SOURCE / "database.csv").read_bytes()).hexdigest(),
        "source_deflection_sha256": hashlib.sha256((SOURCE / "sources/deflection.csv").read_bytes()).hexdigest(),
        "source_center_body_xyz_m": [0, 0, 0], "rho_kg_m3": reference.rho,
        "reference_speed_m_s": reference.velocity,
        "corrections_applied_this_operation": "none; preserve already-corrected physical loads",
        "beta_zero_cleanup_this_operation": "none; preserve supplied loads and control asymmetry",
        "source_note": "Origin explicitly confirmed by user; raw database header (-0.1,0,0.15) conflicts. Source database and deflection zero loads remain inconsistent; recentering does not resolve that discrepancy.",
        "validation_basis": "Original dimensional loads, not rounded source coefficients",
        "cases": [],
    }
    for name, center in CASES.items():
        case_dir = args.output_dir / name
        case_dir.mkdir(exist_ok=True)
        settings = replace(reference, new_center=center)
        output_headers, output_rows = nd.recenter_export_rows(
            base_rows, headers, mapping, settings, zero_at_beta=[], corrections=None,
        )
        database_path = case_dir / "database.csv"
        nd.write_csv_table(database_path, output_headers, output_rows,
                           metadata=settings.metadata_row(mapping, zero_at_beta=[], corrections=None))
        check, _, _ = validate_output(database_path, base_rows, center, settings, mapping)
        validation.append({"case": name, "table": "database", **check})

        output_headers, output_rows = nd.convert_rows(
            def_rows, def_headers, def_mapping, settings, zero_at_beta=[], corrections=None,
        )
        deflection_path = case_dir / "deflection_nondimit.csv"
        nd.write_csv_table(deflection_path, output_headers, output_rows,
                           metadata=settings.metadata_row(def_mapping, zero_at_beta=[], corrections=None))
        check, verified_headers, verified_rows = validate_output(deflection_path, def_rows, center, settings, def_mapping)
        validation.append({"case": name, "table": "deflection", **check})
        # Same Nondimit results, as a plain-header CSV compatible with the FDM
        # control loader. The moment reference is recorded in the paired file.
        nd.write_csv_table(case_dir / "deflection.csv", verified_headers, verified_rows)
        case = {"case": name, "moment_center_body_xyz_m": list(center),
                "rho_kg_m3": reference.rho, "V_ref_m_s": reference.velocity,
                "S_m2": reference.area, "b_m": reference.span, "c_m": reference.chord,
                "output_accuracy": reference.output_accuracy,
                "database": "database.csv", "deflections": "deflection.csv",
                "deflections_with_nondimit_metadata": "deflection_nondimit.csv"}
        (case_dir / "reference.json").write_text(json.dumps(case, indent=2) + "\n")
        provenance["cases"].append(case)
        print(f"{name}: {center}; base={len(base_rows)}, controls={len(def_rows)}, verified", flush=True)

    with (args.output_dir / "validation.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(validation[0]))
        writer.writeheader()
        writer.writerows(validation)
    provenance["max_force_error_n"] = max(row["max_force_error_n"] for row in validation)
    provenance["max_moment_transfer_error_nm"] = max(row["max_moment_transfer_error_nm"] for row in validation)
    provenance["validation_passed"] = True
    (args.output_dir / "manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    fdm_checks = []
    for case in provenance["cases"]:
        directory = args.output_dir / case["case"]
        database = load_aerodynamic_database(directory / "database.csv")
        controls = load_deflection_database(directory / "deflection.csv", database.reference)
        np.testing.assert_allclose(database.reference.moment_center_body_xyz, case["moment_center_body_xyz_m"], rtol=0, atol=0)
        if database.reference.rho != 1.021 or database.alpha_grid_deg.size * database.beta_grid_deg.size != 32761 or controls.case_count != 17:
            raise AssertionError("Generated database could not be loaded consistently by the FDM")
        fdm_checks.append({"case": case["case"], "base_grid_points": 32761,
                           "deflection_cases": 17, "fdm_load_passed": True})
    (args.output_dir / "fdm_validation.json").write_text(json.dumps(fdm_checks, indent=2) + "\n")
    lines = ["# October AeroDB — Nondimit moment-center exports", "",
             "Generated with the actual Nondimit v1.1 conversion backend. Each case is independently recentered from the user-confirmed source origin (0,0,0) m. Density remains 1.021 kg/m³; reference velocity 25 m/s, S=1.805 m², b=3.95 m, c=0.457 m. Existing drag correction is preserved and not reapplied. No beta-zero cleanup is applied during recentering.", "",
             "| Case | x (m) | y (m) | z (m) |", "|---|---:|---:|---:|"]
    lines.extend(f"| {name} | {center[0]:g} | {center[1]:g} | {center[2]:g} |" for name, center in CASES.items())
    lines += ["", "Each folder has `database.csv` (32,761 base points), `deflection_nondimit.csv` (17 cases, native Nondimit metadata and coefficient columns), `deflection.csv` (the same data with a plain CSV header for the FDM control loader), and `reference.json`. Surface angles, case IDs, source velocities and force values are preserved.", "",
              "`validation.csv` checks every row against M_new = M_origin + cross(-r_new, F). Exported density, center, geometry, row counts and coefficient consistency are also checked. Dimensional columns are used as the source, matching Nondimit Modify Export behavior; source coefficients are rounded and can differ slightly from coefficients regenerated from dimensional loads.", "",
              f"Maximum force error: {provenance['max_force_error_n']:.8g} N. Maximum moment-transfer error: {provenance['max_moment_transfer_error_nm']:.8g} N·m (output rounding accuracy 0.0001).", "",
              "The 8 non-origin centers match the previous CG study at z=0.15 m; the origin export is an additional baseline. A moment-reference change does not itself move the physical CG. When loading a case into the FDM, use its actual CSV reference and the desired CG; the FDM performs any remaining transfer once.", "",
              "`fdm_validation.json` confirms that the base and control loaders read all 9 cases with the correct references and grid/case counts. Reproduce the exports from the repository root: `python src/tools/flight_dynamics_model/scripts/export_nondimit_moment_centers.py --output-dir artifacts/aerodb_tail_update_2026_10_07_moment_centers` (FDM requirements plus Python tkinter required).", "",
              "These exports do not resolve the discrepancy between the source database header and the user-confirmed origin, or between the base database and deflection zero-case moments. They are algebraically verified exports, not validated stability/trim results. Numeric updated tail dimensions remain pending. No default aircraft profile or GUI selection was changed by this export task.", ""]
    (args.output_dir / "README.md").write_text("\n".join(lines))
    archive_path = args.output_dir.with_suffix(".zip")
    with zipfile.ZipFile(archive_path, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for path in sorted(args.output_dir.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(args.output_dir.parent))
    print(f"All cases verified; archive={archive_path}")


if __name__ == "__main__":
    main()
