"""Normalize a control-results export without modifying any aerodynamic loads.

Preserve the complete export; zero-control sideslip rows are audit data rather
than control increments. Nonzero controls at other beta planes require an
explicit extension of the current single-plane model and are rejected here.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

FIELDS = ("case_id", "alpha", "beta", "velocity", "aileron", "elevator", "rudder",
          "Fx", "Fy", "Fz", "Mx", "My", "Mz")
SURFACES = ("aileron", "elevator", "rudder")


def normalize(source, output, previous):
    import math
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not set(FIELDS) <= set(reader.fieldnames or []):
            raise ValueError("Results export is missing required case coordinates or loads")
        rows = list(reader)
    usable, audit_only, keys = [], [], set()
    for row in rows:
        if None in row or any(row[key] is None for key in FIELDS):
            raise ValueError("Malformed CSV row")
        numbers = {key: float(row[key]) for key in FIELDS[1:]}
        if not all(math.isfinite(v) for v in numbers.values()) or numbers["velocity"] <= 0:
            raise ValueError("Case coordinates and loads must be finite with positive speed")
        if row.get("status", "ok") != "ok":
            raise ValueError(f"Unsuccessful source case: {row['case_id']}")
        key = tuple(numbers[k] for k in FIELDS[1:7])
        if key in keys:
            raise ValueError("Duplicate aerodynamic case")
        keys.add(key)
        moved = [k for k in SURFACES if numbers[k] != 0]
        if len(moved) > 1:
            raise ValueError("Combined-surface cases need an explicit interaction model")
        if numbers["beta"] != 0:
            if moved:
                raise ValueError("Nonzero-beta control cases need a multi-plane model")
            audit_only.append({k: row[k] for k in FIELDS})
        else:
            usable.append({k: row[k] for k in FIELDS})
    with previous.open(newline="") as stream:
        reader = csv.DictReader(stream)
        reader.fieldnames = [s.strip() for s in reader.fieldnames]
        old = {row["case_id"]: row for row in reader}
    comparisons = []
    for row in usable:
        if row["case_id"] not in old:
            continue
        for key in FIELDS[1:]:
            before, after = float(old[row["case_id"]][key]), float(row[key])
            comparisons.append({"case_id":row["case_id"], "field":key,
                                "before":before, "after":after, "difference":after-before})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(usable)
    return {"source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
            "normalized_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
            "source_cases":len(rows), "active_cases":len(usable),
            "audit_only_sideslip_baselines":audit_only,
            "overlap_cases":len({r['case_id'] for r in comparisons}),
            "overlap_fields":len(comparisons),
            "changed_overlap_fields":[r for r in comparisons if r["difference"] != 0],
            "overlap_comparison":comparisons,
            "ignored_columns":sorted(set(rows[0])-set(FIELDS)) if rows else [],
            "source_convergence":sorted({r.get('convergence', 'unspecified') for r in rows})}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = normalize(args.source, args.output, args.previous)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k != "overlap_comparison"}, indent=2))
