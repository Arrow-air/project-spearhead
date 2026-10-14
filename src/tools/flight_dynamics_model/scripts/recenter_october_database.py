"""Recover origin-referenced loads from the preserved Nondimit base export.

The CFD origin and the coefficient export's final moment center are distinct.
Existing force corrections are already in the export and must not be repeated.
"""

import argparse
import hashlib
import json
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np

FDM_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = FDM_ROOT / "data/aerodynamics/tail_update_2026_10_07"
SOURCE_SHA256 = "6a9561724e7f3059e88e206d8793054ed357723af73b8e7aaec501f6eff83dc1"
sys.path.insert(0, str(FDM_ROOT.parent / "Nondimit"))
import nondimit as nd


def recenter_base(source, destination):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("Preserve the original export; use a different output path")
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ValueError("This converter is pinned to the preserved October export")
    metadata, headers, rows = nd.read_csv_table_with_metadata(source)
    _, original = nd.read_export_settings(metadata)
    if original.old_center != (-0.1, 0.0, 0.15) or original.rho != 1.021:
        raise ValueError("Unexpected source export reference")
    settings = replace(original, new_center=(0.0, 0.0, 0.0), output_accuracy=1e-10)
    mapping = nd.default_mapping(headers)
    output_headers, output_rows = nd.recenter_export_rows(
        rows, headers, mapping, settings,
    )
    force = np.array([[float(row[key]) for key in ("Fx", "Fy", "Fz")] for row in rows])
    moment = np.array([[float(row[key]) for key in ("Mx", "My", "Mz")] for row in rows])
    output_force = np.array([
        [float(row[key]) for key in ("Fx", "Fy", "Fz")] for row in output_rows
    ])
    output_moment = np.array([
        [float(row[key]) for key in ("Mx", "My", "Mz")] for row in output_rows
    ])
    expected = moment + np.cross(np.asarray(original.old_center), force)
    force_error = float(np.max(np.abs(output_force - force)))
    moment_error = float(np.max(np.abs(output_moment - expected)))
    if force_error > 1e-9 or moment_error > 1e-8:
        raise ValueError("Origin recovery failed the independent wrench calculation")
    output_metadata = settings.metadata_row(mapping)
    output_metadata.append([
        "# Origin recovery", "input_export_center_body_xyz=(-0.1,0,0.15)",
        "raw_cfd_origin_body_xyz=(0,0,0)", "method=Nondimit.recenter_export_rows",
        "existing_force_corrections_retained=true", "corrections_reapplied=false",
    ])
    output_metadata.append(["# Previously applied export options", *metadata[2][1:]])
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".csv.tmp")
    nd.write_csv_table(temporary, output_headers, output_rows, metadata=output_metadata)
    # Generated runtime data use LF; preserved deliveries keep their original bytes.
    temporary.write_text(temporary.read_text())
    temporary.replace(destination)
    return {
        "method": "Nondimit.recenter_export_rows",
        "source_sha256": SOURCE_SHA256,
        "source_export_center_body_xyz_m": list(original.old_center),
        "canonical_center_body_xyz_m": list(settings.new_center),
        "rho_kg_m3": settings.rho,
        "output_accuracy": settings.output_accuracy,
        "rows": len(rows),
        "max_force_change_n": force_error,
        "max_origin_transfer_error_nm": moment_error,
        "existing_corrections_reapplied": False,
        "active_sha256": hashlib.sha256(destination.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DATA_ROOT / "database.csv")
    args = parser.parse_args()
    report = recenter_base(DATA_ROOT / "sources/database.csv", args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
