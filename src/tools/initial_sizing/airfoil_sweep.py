"""
airfoil_sweep.py
----------------
Two-stage airfoil sweep through the aircraft sizing model.

Stage 1 — Wing sweep:  Sweep all wing airfoils from the database with a
          fixed symmetric tail. Rank and select the top N wing candidates.

Stage 2 — Tail sweep:  For each top wing candidate, sweep tail airfoils
          (filtered to low camber / near-symmetric). Rank all wing+tail
          combinations.

The Reynolds number selects which database to use.

Usage:
  python airfoil_sweep.py                              # Re=1M, defaults
  python airfoil_sweep.py --re 500000                  # Re=500k
  python airfoil_sweep.py --top-wing 30                # keep top 30 from stage 1
  python airfoil_sweep.py --tail-max-camber 1.0        # tail camber <= 1%
  python airfoil_sweep.py --min-cl-max 1.2 --top 20    # wing filter + show top 20
  python airfoil_sweep.py --detail 5                   # full sizing for top 5

Requires: aircraft_sizing.py, build_airfoil_db.py in the same directory.
"""

import math
import sys
import os
import csv
import argparse

# Allow importing from this directory
sys.path.insert(0, os.path.dirname(__file__))
from aircraft_sizing import size_aircraft
from build_airfoil_db import load_db_as_tuples, db_file_for_re
from fetch_airfoil_polars import closest_re

# ---------------------------------------------------------------------------
# Tee — write to console and file simultaneously
# ---------------------------------------------------------------------------
FILE_TOP_N = 30  # number of results to write to file


def _save_sweep_plot(r: dict, out_dir: str) -> None:
    """Save L and D vs AoA plots for a sizing result that includes sweep data."""
    if "sweep_alpha_deg" not in r:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    alpha = r["sweep_alpha_deg"]
    L = r["sweep_L"]
    D = r["sweep_D"]
    label = r.get("airfoil", "unknown")
    tail = r.get("tail_airfoil", "")
    if tail and tail != "symmetric":
        label += f" + {tail}"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    ax1.plot(alpha, L, "b-o", markersize=4)
    ax1.set_xlabel("AoA (deg)")
    ax1.set_ylabel("Lift (N)")
    ax1.set_title(f"Lift vs AoA — {label}")
    ax1.grid(True, alpha=0.3)

    ax2.plot(alpha, D, "r-o", markersize=4)
    ax2.set_xlabel("AoA (deg)")
    ax2.set_ylabel("Drag (N)")
    ax2.set_title(f"Drag vs AoA — {label}")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    safe_name = label.replace(" ", "_").replace("+", "_").replace("/", "_")
    path = os.path.join(out_dir, f"LD_vs_alpha_{safe_name}.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {path}")

def _write_results_file(filepath: str, re_val: int,
                        ranked: list[dict], show_tail: bool = False,
                        best_detail: dict = None,
                        detail_ranks: list[int] = None) -> None:
    """Write top 30 ranked results and detail for selected candidates to file."""
    from aircraft_sizing import print_results
    old_stdout = sys.stdout
    with open(filepath, "w", encoding="utf-8") as f:
        sys.stdout = f

        print(f"Airfoil Sweep Results — Re = {re_val:,}")
        print(f"{'=' * 60}\n")

        print(f"Top {min(FILE_TOP_N, len(ranked))} Candidates")
        print_table(ranked, top_n=FILE_TOP_N, show_tail=show_tail)

        if best_detail is not None:
            label = best_detail['airfoil']
            if best_detail.get('tail_airfoil') and best_detail['tail_airfoil'] != 'symmetric':
                label += f" + {best_detail['tail_airfoil']}"
            print(f"\n{'=' * 52}")
            print(f"  Full sizing output — best candidate: {label}")
            print_results(best_detail)

        if detail_ranks:
            for rank in detail_ranks:
                if 1 <= rank <= len(ranked):
                    r = ranked[rank - 1]
                    label = r['airfoil']
                    if r.get('tail_airfoil') and r['tail_airfoil'] != 'symmetric':
                        label += f" + {r['tail_airfoil']}"
                    print(f"\n{'=' * 52}")
                    print(f"  Full sizing output — #{rank}: {label}")
                    print_results(r)

        sys.stdout = old_stdout


# ---------------------------------------------------------------------------
# All inputs from shared config
# ---------------------------------------------------------------------------
from config import AIRCRAFT, AIRFOIL, SWEEP

BASE = {**AIRCRAFT, **SWEEP, "re": AIRFOIL["re_db"]}

RE = BASE["re"]

# Default symmetric tail for stage 1
DEFAULT_TAIL = dict(
    tail_cl_0         = 0.0,
    tail_cl_alpha_2d  = 6.28,
    tail_cd_0         = 0.010,
    tail_cm_ac        = 0.0,
)

# ---------------------------------------------------------------------------
# Excluded airfoils (matched by substring, case-insensitive)
# ---------------------------------------------------------------------------
EXCLUDE = SWEEP.get("exclude", [])


# ---------------------------------------------------------------------------
# Helpers
def _is_excluded(name: str) -> bool:
    """Check if an airfoil name matches any exclusion substring."""
    name_lower = name.lower()
    return any(ex.lower() in name_lower for ex in EXCLUDE)


def _load_oratex_concave_set() -> set[str]:
    """Load the set of airfoil names flagged as concave (Oratex-incompatible).

    Reads oratex_compatibility.csv produced by download_geometries.py.
    Returns a set of lowercase airfoil names that have concave lower surfaces.
    """
    csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "oratex_compatibility.csv")
    if not os.path.exists(csv_path):
        print(f"  WARNING: {csv_path} not found — run download_geometries.py first")
        print(f"  Oratex filter disabled.\n")
        return set()

    concave = set()
    with open(csv_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("concave", "").strip().upper() == "YES":
                concave.add(row["name"].strip().lower())
    return concave


# ---------------------------------------------------------------------------
def _wing_kwargs(entry: tuple) -> dict:
    """Build wing keyword args from an airfoil tuple."""
    name, cl_0, cla_2d, cl_max, cd_0, cm_ac = entry[:6]
    kw = dict(
        wing_cl_0        = cl_0,
        wing_cl_alpha_2d = cla_2d,
        wing_cl_max      = cl_max,
        wing_cd_0        = cd_0,
        wing_cm_ac       = cm_ac,
    )
    if len(entry) >= 9 and entry[6] is not None and entry[6] > 0:
        kw["wing_cd_min"]         = entry[6]
        kw["wing_alpha_minD_deg"] = entry[7]
        kw["wing_k_drag"]         = entry[8]
    return kw


def _tail_kwargs(entry: tuple) -> dict:
    """Build tail keyword args from an airfoil tuple."""
    name, cl_0, cla_2d, cl_max, cd_0, cm_ac = entry[:6]
    kw = dict(
        tail_cl_0        = cl_0,
        tail_cl_alpha_2d = cla_2d,
        tail_cd_0        = cd_0,
        tail_cm_ac       = cm_ac,
    )
    if len(entry) >= 9 and entry[6] is not None and entry[6] > 0:
        kw["tail_cd_min"]         = entry[6]
        kw["tail_alpha_minD_deg"] = entry[7]
        kw["tail_k_drag"]         = entry[8]
    return kw


# ---------------------------------------------------------------------------
# Stage 1 — Wing sweep (fixed symmetric tail)
# ---------------------------------------------------------------------------
# Keys in BASE that are sweep-only config, not passed to size_aircraft()
_SWEEP_ONLY_KEYS = {
    "re", "oratex_skin",
    "min_cl_max", "max_cd_0", "min_thickness", "max_thickness",
    "min_abs_thickness", "tail_max_camber", "tail_max_cd_0",
    "min_tail_abs_thickness", "top_wing", "skip_tail",
    "top", "detail", "detail_rank",
    "w_area", "w_drag", "w_itail", "w_stall",
    "exclude",
}


def run_wing_sweep(wing_airfoils: list, base: dict,
                   tail_kw: dict = None) -> list[dict]:
    """Sweep wing airfoils with a fixed tail configuration."""
    if tail_kw is None:
        tail_kw = DEFAULT_TAIL
    sizing_base = {k: v for k, v in base.items() if k not in _SWEEP_ONLY_KEYS}
    results = []
    for entry in wing_airfoils:
        try:
            r = size_aircraft(**sizing_base, **_wing_kwargs(entry), **tail_kw)
            r["airfoil"] = entry[0]
            r["wing_cl_max_2d"] = entry[3]
            r["wing_cd_0_2d"]   = entry[4]
            r["wing_cm_ac_2d"]  = entry[5]
            r["wing_thickness_pct"] = entry[9] if len(entry) > 9 else 0.0
            r["tail_airfoil"]   = "symmetric"
            results.append(r)
        except Exception as e:
            print(f"  [SKIP] {entry[0]}: {e}")
    return results


# ---------------------------------------------------------------------------
# Stage 2 — Tail sweep (for each top wing, sweep tail airfoils)
# ---------------------------------------------------------------------------
def run_tail_sweep(top_wing: list[tuple], tail_airfoils: list,
                   base: dict) -> list[dict]:
    """For each top wing airfoil, sweep all tail candidates."""
    sizing_base = {k: v for k, v in base.items() if k not in _SWEEP_ONLY_KEYS}
    results = []
    total = len(top_wing) * len(tail_airfoils)
    done = 0
    for w_entry in top_wing:
        w_kw = _wing_kwargs(w_entry)
        for t_entry in tail_airfoils:
            t_kw = _tail_kwargs(t_entry)
            done += 1
            try:
                r = size_aircraft(**sizing_base, **w_kw, **t_kw)
                r["airfoil"]        = w_entry[0]
                r["wing_cl_max_2d"] = w_entry[3]
                r["wing_cd_0_2d"]   = w_entry[4]
                r["wing_cm_ac_2d"]  = w_entry[5]
                r["wing_thickness_pct"] = w_entry[9] if len(w_entry) > 9 else 0.0
                r["tail_airfoil"]   = t_entry[0]
                r["tail_thickness_pct"] = t_entry[9] if len(t_entry) > 9 else 0.0
                results.append(r)
            except Exception:
                pass
        print(f"    {done}/{total} combinations evaluated...", flush=True)
    return results


# ---------------------------------------------------------------------------
# Post-sizing filter — minimum absolute wing thickness
# ---------------------------------------------------------------------------
def filter_min_abs_thickness(results: list[dict],
                             min_thick_m: float) -> list[dict]:
    """Remove candidates whose absolute wing thickness is below the minimum.

    Absolute thickness = thickness_pct/100 * chord_wing.
    """
    if min_thick_m <= 0:
        return results
    filtered = []
    for r in results:
        abs_thick = r.get("wing_thickness_pct", 0.0) / 100.0 * r["chord_wing"]
        if abs_thick >= min_thick_m:
            filtered.append(r)
    return filtered


def filter_min_tail_thickness(results: list[dict],
                              min_thick_m: float,
                              tail_chord: float) -> list[dict]:
    """Remove candidates whose absolute tail thickness is below the minimum.

    Absolute thickness = tail_thickness_pct/100 * tail_chord.
    """
    if min_thick_m <= 0:
        return results
    filtered = []
    for r in results:
        abs_thick = r.get("tail_thickness_pct", 0.0) / 100.0 * tail_chord
        if abs_thick >= min_thick_m:
            filtered.append(r)
    return filtered


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------
def score(results: list[dict], weights: dict = None) -> list[dict]:
    """
    Composite score combining wing area (weight proxy), total drag, tail
    incidence magnitude, and a stall safety penalty.

    Lower score = better.

    Normalise each metric to [0, 1] across candidates, then weight.
    Weights are read from BASE dict (w_area, w_drag, w_itail, w_stall).
    """
    if weights is None:
        weights = BASE
    w_a = weights.get("w_area",  0.40)
    w_d = weights.get("w_drag",  0.30)
    w_i = weights.get("w_itail", 0.20)
    w_m = weights.get("w_stall", 0.10)
    def norm(vals):
        lo, hi = min(vals), max(vals)
        if hi == lo:
            return [0.5] * len(vals)
        return [(v - lo) / (hi - lo) for v in vals]

    S_wing  = [r["S_wing_eff"]  for r in results]
    D_total = [r["D_total"]     for r in results]
    i_tail  = [abs(r["tail_incidence_deg"]) for r in results]
    M_pen   = [max(0.0, r["M_stall"]) for r in results]

    nS  = norm(S_wing)
    nD  = norm(D_total)
    nI  = norm(i_tail)
    nM  = norm(M_pen)

    for i, r in enumerate(results):
        r["score"] = w_a*nS[i] + w_d*nD[i] + w_i*nI[i] + w_m*nM[i]

    return sorted(results, key=lambda r: r["score"])


# ---------------------------------------------------------------------------
# Print comparison table
# ---------------------------------------------------------------------------
def print_table(results: list[dict], top_n: int = None,
                show_tail: bool = False) -> None:
    if top_n is not None:
        results = results[:top_n]

    # Auto-size name columns to fit the longest name in the displayed results
    WN = max(len("Wing Airfoil"), max((len(r["airfoil"]) for r in results), default=12))
    if show_tail:
        TN = max(len("Tail Airfoil"), max((len(r["tail_airfoil"]) for r in results), default=12))
    else:
        TN = 0

    nums = (f"{'S_wing':>7} {'Span':>5} {'Chord':>5} {'t_abs':>5}"
            f" {'i_tail':>6} {'D_tot':>6} {'L/D':>5}"
            f" {'CLmax':>6} {'M_stall':>7} {'Stall':>8} {'Score':>5}")
    units = (f"{'(m2)':>7} {'(m)':>5} {'(m)':>5} {'(cm)':>5}"
             f" {'(deg)':>6} {'(N)':>6} {'':>5}"
             f" {'':>6} {'(N.m)':>7} {'':>8} {'':>5}")

    if show_tail:
        hdr = f"  {'#':<4} {'Wing Airfoil':<{WN}} | {'Tail Airfoil':<{TN}} | {nums}"
        unit_line = f"  {'':4} {'':>{WN}} | {'':>{TN}} | {units}"
    else:
        hdr = f"  {'#':<4} {'Airfoil':<{WN}} | {nums}"
        unit_line = f"  {'':4} {'':>{WN}} | {units}"

    sep = "-" * len(hdr)

    print(f"\n{'=' * len(hdr)}")
    print("  Airfoil Sweep — Wing Sizing Comparison")
    print(f"  Ranked by composite score ({BASE['w_area']:.0%} wing area, {BASE['w_drag']:.0%} drag, {BASE['w_itail']:.0%} tail, {BASE['w_stall']:.0%} stall safety)")
    print(f"{'=' * len(hdr)}")
    print(hdr)
    print(unit_line)
    print(f"  {sep}")

    for rank, r in enumerate(results, 1):
        stall_flag = "NOSE-UP!" if r["M_stall"] > 0 else "ok"
        t_abs_cm = r.get("wing_thickness_pct", 0.0) / 100.0 * r["chord_wing"] * 100.0
        data = (f"{r['S_wing_eff']:>7.3f} {r['span_wing']:>5.2f} {r['chord_wing']:>5.3f} {t_abs_cm:>5.1f}"
                f" {r['tail_incidence_deg']:>6.2f} {r['D_total']:>6.1f} {r['LD']:>5.1f}"
                f" {r['CL_max_3d']:>6.3f} {r['M_stall']:>7.2f} {stall_flag:>8} {r['score']:>5.3f}")
        if show_tail:
            print(f"  {rank:<4} {r['airfoil']:<{WN}} | {r['tail_airfoil']:<{TN}} | {data}")
        else:
            print(f"  {rank:<4} {r['airfoil']:<{WN}} | {data}")

    print(f"  {sep}")
    print(f"\n  Best overall: {results[0]['airfoil']}"
          + (f" + {results[0]['tail_airfoil']}" if show_tail else ""))
    print(f"  Smallest wing: {min(results, key=lambda r: r['S_wing_eff'])['airfoil']}")
    print(f"  Lowest drag:   {min(results, key=lambda r: r['D_total'])['airfoil']}")
    print(f"  Best L/D:      {max(results, key=lambda r: r['LD'])['airfoil']}\n")


def print_top_detail(results: list[dict], n: int = 3) -> None:
    """Print the full sizing output for the top N candidates."""
    from aircraft_sizing import print_results
    print(f"\n{'='*52}")
    print(f"  Full sizing output - top {n} candidates")
    for r in results[:n]:
        label = r['airfoil']
        if r.get('tail_airfoil') and r['tail_airfoil'] != 'symmetric':
            label += f" + {r['tail_airfoil']}"
        print(f"\n  >>> {label} <<<")
        print_results(r)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Two-stage airfoil sweep through sizing model")

    # Reynolds number
    parser.add_argument("--re", type=int, default=RE,
                        help=f"Target Reynolds number (default: {RE})")

    # Wing filters (defaults from BASE dict)
    parser.add_argument("--min-cl-max", type=float, default=BASE["min_cl_max"],
                        help=f"Wing filter: minimum cl_max (default: {BASE['min_cl_max']})")
    parser.add_argument("--max-cd-0", type=float, default=BASE["max_cd_0"],
                        help=f"Wing filter: maximum cd_0 (default: {BASE['max_cd_0']})")
    parser.add_argument("--min-thickness", type=float, default=BASE["min_thickness"],
                        help=f"Wing filter: minimum thickness %% (default: {BASE['min_thickness']})")
    parser.add_argument("--max-thickness", type=float, default=BASE["max_thickness"],
                        help=f"Wing filter: maximum thickness %% (default: {BASE['max_thickness']})")
    parser.add_argument("--min-abs-thickness", type=float, default=BASE["min_abs_thickness"],
                        help=f"Wing filter: minimum absolute thickness [m] (default: {BASE['min_abs_thickness']})")

    # Stage control
    parser.add_argument("--top-wing", type=int, default=BASE["top_wing"],
                        help=f"Number of top wing candidates for stage 2 (default: {BASE['top_wing']})")
    parser.add_argument("--tail-max-camber", type=float, default=BASE["tail_max_camber"],
                        help=f"Tail filter: max camber %% (default: {BASE['tail_max_camber']})")
    parser.add_argument("--tail-max-cd-0", type=float, default=BASE["tail_max_cd_0"],
                        help=f"Tail filter: max cd_0 (default: {BASE['tail_max_cd_0']})")
    parser.add_argument("--min-tail-abs-thickness", type=float, default=BASE["min_tail_abs_thickness"],
                        help=f"Tail filter: minimum absolute thickness [m] (default: {BASE['min_tail_abs_thickness']})")
    parser.add_argument("--skip-tail", action="store_true", default=BASE["skip_tail"],
                        help="Skip stage 2 (wing sweep only)")

    # Output
    parser.add_argument("--top", type=int, default=BASE["top"],
                        help="Only show top N results in table")
    parser.add_argument("--detail", type=int, default=BASE["detail"],
                        help="Print full sizing output for top N candidates")
    parser.add_argument("--detail-rank", type=int, nargs="+", default=BASE["detail_rank"],
                        help="Print full sizing output for specific rank(s), e.g. --detail-rank 3 8")
    args = parser.parse_args()

    re_val = closest_re(args.re)
    db_path = db_file_for_re(re_val)
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            f"sweep_results_re{re_val}.txt")

    # -----------------------------------------------------------------------
    # Stage 1: Wing sweep
    # -----------------------------------------------------------------------
    print(f"Loading airfoil database: {db_path}")
    print(f"Reynolds number: {re_val:,}")

    # Oratex filter: load concave airfoil set if enabled in BASE
    oratex_concave = set()
    if BASE.get("oratex_skin", False):
        oratex_concave = _load_oratex_concave_set()
        if oratex_concave:
            print(f"Oratex filter: ON ({len(oratex_concave)} concave airfoils will be excluded)")
    else:
        print(f"Oratex filter: OFF")
    print()

    wing_airfoils = [af for af in load_db_as_tuples(
        db_path=db_path,
        min_cl_max=args.min_cl_max,
        max_cd_0=args.max_cd_0,
        min_thickness=args.min_thickness,
        max_thickness=args.max_thickness,
    ) if not _is_excluded(af[0])
       and af[0].strip().lower() not in oratex_concave]
    print(f"Stage 1: Wing sweep ({len(wing_airfoils)} candidates, fixed symmetric tail)")

    if not wing_airfoils:
        print("No wing airfoils match the filters.")
        sys.exit(1)

    raw_wing = run_wing_sweep(wing_airfoils, BASE)
    if args.min_abs_thickness > 0:
        before = len(raw_wing)
        raw_wing = filter_min_abs_thickness(raw_wing, args.min_abs_thickness)
        print(f"  Absolute thickness filter (>= {args.min_abs_thickness*100:.1f} cm): "
              f"{before} -> {len(raw_wing)} candidates")
    ranked_wing = score(raw_wing)

    print(f"  Stage 1 complete — top {min(args.top_wing, len(ranked_wing))} wing candidates selected\n")
    print_table(ranked_wing, top_n=args.top_wing)

    if args.skip_tail:
        _write_results_file(out_path, re_val, ranked_wing,
                            best_detail=ranked_wing[0],
                            detail_ranks=args.detail_rank)
        if args.detail > 0:
            print_top_detail(ranked_wing, args.detail)
        for rank in args.detail_rank:
            if 1 <= rank <= len(ranked_wing):
                print_top_detail([ranked_wing[rank - 1]], 1)
                _save_sweep_plot(ranked_wing[rank - 1], os.path.dirname(out_path))
        _save_sweep_plot(ranked_wing[0], os.path.dirname(out_path))
        print(f"\nResults saved to: {out_path}")
        sys.exit(0)

    # -----------------------------------------------------------------------
    # Stage 2: Tail sweep for top wing candidates
    # -----------------------------------------------------------------------
    tail_airfoils = [af for af in load_db_as_tuples(
        db_path=db_path,
        max_cd_0=args.tail_max_cd_0,
        max_camber=args.tail_max_camber,
    ) if not _is_excluded(af[0])
       and af[0].strip().lower() not in oratex_concave]
    print(f"Stage 2: Tail sweep ({len(tail_airfoils)} tail candidates "
          f"x {min(args.top_wing, len(ranked_wing))} wing candidates "
          f"= {len(tail_airfoils) * min(args.top_wing, len(ranked_wing)):,} combinations)")

    if not tail_airfoils:
        print("No tail airfoils match the filters. Try increasing --tail-max-camber.")
        sys.exit(1)

    # Extract top wing entries (need original tuples for _wing_kwargs)
    top_wing_names = {r["airfoil"] for r in ranked_wing[:args.top_wing]}
    top_wing_tuples = [af for af in wing_airfoils if af[0] in top_wing_names]

    raw_combo = run_tail_sweep(top_wing_tuples, tail_airfoils, BASE)
    if args.min_abs_thickness > 0:
        raw_combo = filter_min_abs_thickness(raw_combo, args.min_abs_thickness)
    if args.min_tail_abs_thickness > 0:
        before = len(raw_combo)
        raw_combo = filter_min_tail_thickness(
            raw_combo, args.min_tail_abs_thickness, BASE["tail_chord"])
        print(f"  Tail thickness filter (>= {args.min_tail_abs_thickness*100:.1f} cm): "
              f"{before} -> {len(raw_combo)} combinations")
    ranked_combo = score(raw_combo)

    print(f"\n  Stage 2 complete — {len(ranked_combo)} valid combinations\n")
    print_table(ranked_combo, top_n=args.top, show_tail=True)

    # Write results file (top 30 combos + best candidate detail)
    _write_results_file(out_path, re_val, ranked_combo, show_tail=True,
                        best_detail=ranked_combo[0],
                        detail_ranks=args.detail_rank)
    if args.detail > 0:
        print_top_detail(ranked_combo, args.detail)
    for rank in args.detail_rank:
        if 1 <= rank <= len(ranked_combo):
            print_top_detail([ranked_combo[rank - 1]], 1)
            _save_sweep_plot(ranked_combo[rank - 1], os.path.dirname(out_path))
    _save_sweep_plot(ranked_combo[0], os.path.dirname(out_path))
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
