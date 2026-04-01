"""
fetch_airfoil_polars.py
-----------------------
Fetch XFOIL polar data from airfoiltools.com and extract the key 2D
aerodynamic parameters needed by aircraft_sizing.py / airfoil_sweep.py.

For each airfoil the script downloads the polar CSV at the closest
available Reynolds number, then extracts:
  - cl_0        : lift coefficient at alpha = 0
  - cl_alpha_2d : lift-curve slope [per radian]
  - cl_max      : maximum lift coefficient
  - cd_0        : drag coefficient at alpha = 0
  - cm_ac       : pitching moment at alpha = 0 (approx. cm_ac for thin airfoils)

Available Reynolds numbers on airfoiltools.com:
  50 000, 100 000, 200 000, 500 000, 1 000 000

Usage:
  python fetch_airfoil_polars.py              # fetch all candidates
  python fetch_airfoil_polars.py --re 500000  # at Re = 500k
  python fetch_airfoil_polars.py --sweep      # fetch + run airfoil_sweep

Requires: requests  (pip install requests)
"""

import sys
import os
import math
import argparse
from io import StringIO

try:
    import requests
except ImportError:
    print("Please install requests:  pip install requests")
    sys.exit(1)

# -------------------------------------------------------------------------
# Airfoil slug mapping
# -------------------------------------------------------------------------
# Maps a human-readable name to the airfoiltools.com slug (the part before
# the Reynolds number in the polar key).
# To add an airfoil: find its page on airfoiltools.com and copy the slug
# from the URL, e.g.  airfoil/details?airfoil=<SLUG>
AIRFOIL_SLUGS = {
    "NACA 2412":  "naca2412-il",
    "NACA 4412":  "naca4412-il",
    "NACA 4415":  "naca4415-il",
    "NACA 4421":  "naca4421-il",
    "NACA 23012": "naca23012-il",
    "NACA 23015": "naca23015-il",
    "Clark Y":    "clarky-il",
    "S1223":      "s1223-il",
    "E423":       "e423-il",
    "FX 63-137":  "fx63137-il",
    "Eppler 387": "e387-il",
    "NACA 0012":  "naca0012-il",
    "NACA 2415":  "naca2415-il",
}

AVAILABLE_RE = [50000, 100000, 200000, 500000, 1000000]

CSV_URL = "http://airfoiltools.com/polar/csv?polar=xf-{slug}-{re}"


# -------------------------------------------------------------------------
# Fetch & parse
# -------------------------------------------------------------------------
def fetch_polar_csv(slug: str, re: int, timeout: float = 15.0) -> str:
    """Download raw polar CSV text from airfoiltools.com."""
    url = CSV_URL.format(slug=slug, re=re)
    resp = requests.get(url, headers={"User-Agent": "Mozilla/5.0"}, timeout=timeout)
    resp.raise_for_status()
    # Check that we actually got CSV, not an HTML error page
    if resp.text.strip().startswith("<!DOCTYPE") or "<html" in resp.text[:200]:
        raise ValueError(f"Got HTML instead of CSV - polar may not exist: {url}")
    return resp.text


def parse_polar(csv_text: str) -> dict:
    """
    Parse airfoiltools.com polar CSV into a dict with arrays:
      alpha, cl, cd, cdp, cm, top_xtr, bot_xtr
    Also extracts metadata from the header.
    """
    lines = csv_text.strip().splitlines()
    meta = {}
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith("Alpha,"):
            data_start = i + 1
            break
        parts = line.split(",", 1)
        if len(parts) == 2:
            meta[parts[0].strip()] = parts[1].strip()

    if data_start is None:
        raise ValueError("Could not find data header in polar CSV")

    alpha, cl, cd, cdp, cm = [], [], [], [], []
    for line in lines[data_start:]:
        parts = line.split(",")
        if len(parts) < 5:
            continue
        alpha.append(float(parts[0]))
        cl.append(float(parts[1]))
        cd.append(float(parts[2]))
        cdp.append(float(parts[3]))
        cm.append(float(parts[4]))

    return {"meta": meta, "alpha": alpha, "cl": cl, "cd": cd, "cdp": cdp, "cm": cm}


def _fit_drag_polar(alpha_deg: list, cd: list) -> tuple[float, float, float]:
    """
    Fit a quadratic drag polar:  cd = cd_min + k_drag * (alpha - alpha_minD)^2

    Anchored fit: cd_min and alpha_minD are taken directly from the data
    (the point with lowest cd in the -4..+12 deg range), then k_drag is
    fitted by least-squares around that anchor.  This guarantees cd_min
    is always a real measured value (never negative).

    Alpha is in degrees throughout (stored and used in degrees).

    Returns (cd_min, alpha_minD_deg, k_drag).
    """
    # Filter to pre-stall range for a clean fit
    a_filt, cd_filt = [], []
    for a, d in zip(alpha_deg, cd):
        if -4.0 <= a <= 12.0:
            a_filt.append(a)
            cd_filt.append(d)

    n = len(a_filt)
    if n < 4:
        return (float("nan"), float("nan"), float("nan"))

    # Anchor: find the actual minimum drag point in the data
    idx_min = min(range(n), key=lambda i: cd_filt[i])
    cd_min = cd_filt[idx_min]
    alpha_minD = a_filt[idx_min]

    # Fit k_drag by least-squares: minimize sum((cd_i - cd_min) - k*(a_i - a0)^2)^2
    # Solution: k = sum((cd_i - cd_min) * (a_i - a0)^2) / sum((a_i - a0)^4)
    num = 0.0
    den = 0.0
    for a, d in zip(a_filt, cd_filt):
        da = a - alpha_minD
        num += (d - cd_min) * da * da
        den += da ** 4

    if den < 1e-30:
        return (float("nan"), float("nan"), float("nan"))

    k_drag = num / den
    if k_drag <= 0:
        return (float("nan"), float("nan"), float("nan"))

    return (cd_min, alpha_minD, k_drag)


def extract_params(polar: dict) -> dict:
    """
    Extract aerodynamic parameters from an XFOIL polar:
      cl_0, cl_alpha_2d, cl_max, alpha_stall, cd_0, cm_ac,
      cd_min, alpha_minD_deg, k_drag

    cl_alpha_2d: linear regression of cl vs alpha in the linear range
                 (alpha from -4 to +8 degrees), converted to per-radian.
    cl_0:        cl at alpha = 0 (interpolated from the linear fit).
    cl_max:      maximum cl in the polar.
    cd_0:        cd at alpha = 0 (interpolated).
    cm_ac:       cm at alpha = 0 (interpolated, good approx for cm_ac).
    cd_min, alpha_minD_deg, k_drag:  quadratic drag polar fit
                 cd = cd_min + k_drag * (alpha_deg - alpha_minD_deg)^2
    """
    alpha = polar["alpha"]
    cl = polar["cl"]
    cd = polar["cd"]
    cm = polar["cm"]

    # --- cl_max and alpha_stall ---
    cl_max = max(cl)
    idx_max = cl.index(cl_max)
    alpha_stall = alpha[idx_max]

    # --- Linear fit for cl(alpha) in the linear range ---
    # Use points between -4 and +8 degrees
    fit_a, fit_cl = [], []
    for a, c in zip(alpha, cl):
        if -4.0 <= a <= 8.0:
            fit_a.append(math.radians(a))
            fit_cl.append(c)

    if len(fit_a) < 3:
        raise ValueError("Not enough data points in -4..+8 deg range for linear fit")

    # Simple least-squares: cl = cl_alpha * alpha_rad + cl_0
    n = len(fit_a)
    sum_x = sum(fit_a)
    sum_y = sum(fit_cl)
    sum_xx = sum(x * x for x in fit_a)
    sum_xy = sum(x * y for x, y in zip(fit_a, fit_cl))

    cl_alpha_2d = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x ** 2)
    cl_0 = (sum_y - cl_alpha_2d * sum_x) / n

    # --- cd_0 and cm_ac: interpolate at alpha = 0 ---
    cd_0 = _interp_at_zero(alpha, cd)
    cm_ac = _interp_at_zero(alpha, cm)

    # --- Drag polar fit: cd = cd_min + k_drag * (alpha - alpha_minD)^2 ---
    cd_min, alpha_minD_deg, k_drag = _fit_drag_polar(alpha, cd)

    return {
        "cl_0": round(cl_0, 4),
        "cl_alpha_2d": round(cl_alpha_2d, 3),
        "cl_max": round(cl_max, 4),
        "alpha_stall": round(alpha_stall, 2),
        "cd_0": round(cd_0, 5),
        "cm_ac": round(cm_ac, 4),
        "cd_min": round(cd_min, 6) if not math.isnan(cd_min) else float("nan"),
        "alpha_minD_deg": round(alpha_minD_deg, 4) if not math.isnan(alpha_minD_deg) else float("nan"),
        "k_drag": round(k_drag, 6) if not math.isnan(k_drag) else float("nan"),
    }


def _interp_at_zero(alpha: list, vals: list) -> float:
    """Linear interpolation of vals at alpha = 0."""
    for i in range(len(alpha) - 1):
        if alpha[i] <= 0.0 <= alpha[i + 1]:
            t = (0.0 - alpha[i]) / (alpha[i + 1] - alpha[i])
            return vals[i] + t * (vals[i + 1] - vals[i])
    # Fallback: find closest to 0
    idx = min(range(len(alpha)), key=lambda i: abs(alpha[i]))
    return vals[idx]


def closest_re(target_re: int) -> int:
    """Return the closest available Reynolds number."""
    return min(AVAILABLE_RE, key=lambda r: abs(r - target_re))


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def fetch_all(target_re: int = 1000000, airfoils: dict = None) -> list[tuple]:
    """
    Fetch polars for all airfoils and return a list of tuples matching
    the format used in airfoil_sweep.py:
      (name, cl_0, cl_alpha_2d, cl_max, cd_0, cm_ac)
    """
    if airfoils is None:
        airfoils = AIRFOIL_SLUGS

    re = closest_re(target_re)
    print(f"Fetching polars at Re = {re:,} (closest to {target_re:,})")
    print(f"Source: airfoiltools.com (XFOIL, Ncrit=9)\n")

    results = []
    for name, slug in airfoils.items():
        try:
            csv_text = fetch_polar_csv(slug, re)
            polar = parse_polar(csv_text)
            params = extract_params(polar)
            results.append((
                name,
                params["cl_0"],
                params["cl_alpha_2d"],
                params["cl_max"],
                params["cd_0"],
                params["cm_ac"],
            ))
            print(f"  OK  {name:<15} cl_0={params['cl_0']:+.3f}  "
                  f"cla={params['cl_alpha_2d']:.2f}/rad  "
                  f"cl_max={params['cl_max']:.3f}  "
                  f"cd_0={params['cd_0']:.5f}  "
                  f"cm_ac={params['cm_ac']:+.4f}")
        except Exception as e:
            print(f"  FAIL {name:<15} {e}")

    return results


def print_as_python(results: list[tuple], re: int) -> None:
    """Print the results as a copy-pasteable Python list for airfoil_sweep.py."""
    print(f"\n# --- Auto-fetched from airfoiltools.com, Re = {re:,} ---")
    print("AIRFOILS = [")
    print("    # name           cl_0   cla_2d  cl_max  cd_0    cm_ac")
    for name, cl_0, cla, cl_max, cd_0, cm_ac in results:
        print(f'    ("{name}",{cl_0:>10.2f}, {cla:>6.2f}, {cl_max:>7.2f}, '
              f'{cd_0:>7.3f}, {cm_ac:>7.3f}),')
    print("]")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch airfoil polars from airfoiltools.com")
    parser.add_argument("--re", type=int, default=1000000,
                        help="Target Reynolds number (default: 1000000)")
    parser.add_argument("--sweep", action="store_true",
                        help="After fetching, run the airfoil sweep with live data")
    args = parser.parse_args()

    re_actual = closest_re(args.re)
    results = fetch_all(args.re)

    if not results:
        print("\nNo airfoils fetched successfully.")
        sys.exit(1)

    print_as_python(results, re_actual)

    if args.sweep:
        print("\n" + "=" * 60)
        print("  Running airfoil sweep with fetched data...")
        print("=" * 60)
        sys.path.insert(0, os.path.dirname(__file__))
        from airfoil_sweep import run_sweep, score, print_table, BASE
        raw = run_sweep(results, BASE)
        ranked = score(raw)
        print_table(ranked)
