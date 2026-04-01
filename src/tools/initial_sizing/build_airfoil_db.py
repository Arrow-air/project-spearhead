"""
build_airfoil_db.py
-------------------
Scrape ALL airfoils from airfoiltools.com, fetch their XFOIL polars,
extract key aerodynamic parameters, and save to a CSV database.

The database can then be loaded by airfoil_sweep.py to sweep any/all
airfoils through the sizing model without re-fetching from the web.

Output: airfoil_db.csv with columns:
  slug, name, re, cl_0, cl_alpha_2d, cl_max, cd_0, cm_ac,
  thickness_pct, camber_pct

Usage:
  python build_airfoil_db.py                    # build full DB at Re=1M
  python build_airfoil_db.py --re 500000        # at Re=500k
  python build_airfoil_db.py --letter n         # only letter N (NACA etc.)
  python build_airfoil_db.py --resume           # resume interrupted build
  python build_airfoil_db.py --filter "naca"    # only slugs containing "naca"

Requires: requests  (pip install requests)
"""

import sys
import os
import re as regex
import csv
import time
import argparse
import math

try:
    import requests
except ImportError:
    print("Please install requests:  pip install requests")
    sys.exit(1)

# Reuse core functions from fetch_airfoil_polars
sys.path.insert(0, os.path.dirname(__file__))
from fetch_airfoil_polars import (
    fetch_polar_csv, parse_polar, extract_params, closest_re, AVAILABLE_RE,
)

SEARCH_URL = "http://airfoiltools.com/search/list?page={letter}&no={page_no}"
DETAIL_URL = "http://airfoiltools.com/airfoil/details?airfoil={slug}"
HEADERS = {"User-Agent": "Mozilla/5.0"}
DB_DIR = os.path.dirname(os.path.abspath(__file__))
DB_FILE = os.path.join(DB_DIR, "airfoil_db.csv")  # default, overridden by --re


def db_file_for_re(re_val: int) -> str:
    """Return the database filename for a given Reynolds number."""
    return os.path.join(DB_DIR, f"airfoil_db_re{re_val}.csv")

LETTERS = list("abcdefghijklmnopqrstuvwxyz")
AIRFOILS_PER_PAGE = 20

CSV_COLUMNS = [
    "slug", "name", "re", "cl_0", "cl_alpha_2d", "cl_max",
    "alpha_stall", "cd_0", "cm_ac", "cd_min", "alpha_minD_deg", "k_drag",
    "thickness_pct", "camber_pct",
]


# -------------------------------------------------------------------------
# Scrape airfoil slugs from airfoiltools.com
# -------------------------------------------------------------------------
def scrape_slugs_for_letter(letter: str) -> list[dict]:
    """
    Scrape all airfoil slugs, names, thickness, and camber for a given letter page.
    Returns list of dicts with 'slug', 'name', 'thickness', 'camber' keys.
    """
    airfoils = {}
    for page_no in range(100):  # safety limit
        url = SEARCH_URL.format(letter=letter, page_no=page_no)
        try:
            resp = requests.get(url, headers=HEADERS, timeout=15)
            resp.raise_for_status()
        except Exception as e:
            print(f"    Error fetching page {letter}/{page_no}: {e}")
            break

        # Extract airfoil slugs and names from the page
        # Name is in h3 tag: <h3>(slug) AIRFOIL NAME</h3>
        matches = regex.findall(
            r'<h3>\(([a-z0-9_.-]+)\)\s*(.+?)</h3>',
            resp.text, regex.I
        )
        new_count = 0
        for slug, name in matches:
            slug = slug.strip()
            name = name.strip()
            if slug not in airfoils:
                # Extract thickness and camber from the same page
                # Pattern: Max thickness XX.X% ... Max camber YY.Y%
                # These appear in the cell2 td after the airfoil entry
                thick_pat = regex.search(
                    rf'{regex.escape(slug)}.*?Max thickness\s+([0-9.]+)\s*%',
                    resp.text, regex.S | regex.I
                )
                camber_pat = regex.search(
                    rf'{regex.escape(slug)}.*?Max camber\s+([0-9.]+)\s*%',
                    resp.text, regex.S | regex.I
                )
                thick = float(thick_pat.group(1)) if thick_pat else float("nan")
                camber = float(camber_pat.group(1)) if camber_pat else float("nan")
                airfoils[slug] = {
                    "name": name, "thickness": thick, "camber": camber
                }
                new_count += 1

        if new_count == 0:
            break  # no new airfoils on this sub-page => done

        time.sleep(0.3)  # be polite

    return [{"slug": s, "name": d["name"],
             "thickness": d["thickness"], "camber": d["camber"]}
            for s, d in sorted(airfoils.items())]


def scrape_all_slugs(letters: list[str] = None) -> list[dict]:
    """Scrape airfoil slugs across all letter pages."""
    if letters is None:
        letters = LETTERS

    all_airfoils = []
    for letter in letters:
        print(f"  Scraping letter {letter.upper()}...", end="", flush=True)
        found = scrape_slugs_for_letter(letter)
        print(f" {len(found)} airfoils")
        all_airfoils.extend(found)
        time.sleep(0.3)

    # Deduplicate by slug
    seen = set()
    unique = []
    for af in all_airfoils:
        if af["slug"] not in seen:
            seen.add(af["slug"])
            unique.append(af)

    return unique


def scrape_thickness_camber(slug: str) -> tuple[float, float]:
    """
    Scrape thickness and camber from the airfoil detail page.
    Returns (thickness_pct, camber_pct) or (NaN, NaN) on failure.
    """
    nan = float("nan")
    try:
        url = DETAIL_URL.format(slug=slug)
        resp = requests.get(url, headers=HEADERS, timeout=10)
        resp.raise_for_status()

        thick = regex.search(r'Max thickness[^0-9]*([0-9.]+)\s*%', resp.text)
        camber = regex.search(r'Max camber[^0-9]*([0-9.]+)\s*%', resp.text)

        t = float(thick.group(1)) if thick else nan
        c = float(camber.group(1)) if camber else nan
        return t, c
    except Exception:
        return nan, nan


# -------------------------------------------------------------------------
# Database I/O
# -------------------------------------------------------------------------
def load_existing_db(db_path: str = None) -> dict[str, dict]:
    """Load existing database rows keyed by slug."""
    path = db_path or DB_FILE
    if not os.path.exists(path):
        return {}
    existing = {}
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            existing[row["slug"]] = row
    return existing


def save_db(rows: list[dict], db_path: str = None) -> None:
    """Save database rows to CSV."""
    path = db_path or DB_FILE
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in sorted(rows, key=lambda r: r["slug"]):
            writer.writerow(row)


# -------------------------------------------------------------------------
# Build database
# -------------------------------------------------------------------------
def build_db(target_re: int = 1000000, letters: list[str] = None,
             slug_filter: str = None, resume: bool = False,
             delay: float = 0.2) -> None:
    """
    Main database builder.
    1. Scrape all airfoil slugs from airfoiltools.com
    2. For each, fetch the XFOIL polar at the target Re
    3. Extract parameters and save to CSV
    """
    re_val = closest_re(target_re)
    out_path = db_file_for_re(re_val)
    print(f"Building airfoil database at Re = {re_val:,}")
    print(f"Output: {out_path}\n")

    # Step 1: Scrape slugs
    print("Step 1: Scraping airfoil list from airfoiltools.com...")
    all_airfoils = scrape_all_slugs(letters)
    print(f"  Total unique airfoils found: {len(all_airfoils)}\n")

    # Apply filter
    if slug_filter:
        all_airfoils = [af for af in all_airfoils
                        if slug_filter.lower() in af["slug"].lower()
                        or slug_filter.lower() in af["name"].lower()]
        print(f"  After filter '{slug_filter}': {len(all_airfoils)} airfoils\n")

    # Load existing DB for resume
    existing = load_existing_db(out_path) if resume else {}
    if resume and existing:
        print(f"  Resuming: {len(existing)} airfoils already in database\n")

    # Step 2: Fetch polars and extract params
    print("Step 2: Fetching polars and extracting parameters...")
    rows = list(existing.values()) if resume else []
    ok_count = len(existing) if resume else 0
    fail_count = 0
    skip_count = 0

    for i, af in enumerate(all_airfoils):
        slug = af["slug"]
        name = af["name"]

        # Skip if already in DB (resume mode)
        if resume and slug in existing:
            skip_count += 1
            continue

        thick = af.get("thickness", float("nan"))
        camber = af.get("camber", float("nan"))
        progress = f"[{i+1}/{len(all_airfoils)}]"
        try:
            csv_text = fetch_polar_csv(slug, re_val)
            polar = parse_polar(csv_text)
            params = extract_params(polar)

            _nan = float("nan")
            cd_min = params.get("cd_min", _nan)
            alpha_minD_deg = params.get("alpha_minD_deg", _nan)
            k_drag = params.get("k_drag", _nan)

            row = {
                "slug": slug,
                "name": name,
                "re": re_val,
                "cl_0": f"{params['cl_0']:.4f}",
                "cl_alpha_2d": f"{params['cl_alpha_2d']:.3f}",
                "cl_max": f"{params['cl_max']:.4f}",
                "alpha_stall": f"{params['alpha_stall']:.2f}",
                "cd_0": f"{params['cd_0']:.5f}",
                "cm_ac": f"{params['cm_ac']:.4f}",
                "cd_min": f"{cd_min:.6f}" if not math.isnan(cd_min) else "",
                "alpha_minD_deg": f"{alpha_minD_deg:.4f}" if not math.isnan(alpha_minD_deg) else "",
                "k_drag": f"{k_drag:.6f}" if not math.isnan(k_drag) else "",
                "thickness_pct": f"{thick:.2f}" if not math.isnan(thick) else "",
                "camber_pct": f"{camber:.2f}" if not math.isnan(camber) else "",
            }
            rows.append(row)
            ok_count += 1

            print(f"  {progress} OK   {slug:<25} cl_max={params['cl_max']:.3f}  "
                  f"a_stall={params['alpha_stall']:.1f}deg  "
                  f"cd_0={params['cd_0']:.5f}  cm={params['cm_ac']:+.4f}  "
                  f"t={thick:.1f}%")

        except requests.exceptions.HTTPError as e:
            if "404" in str(e):
                # No polar at this Re - common for many airfoils
                fail_count += 1
                if fail_count <= 20 or fail_count % 50 == 0:
                    print(f"  {progress} SKIP {slug:<25} no polar at Re={re_val:,}")
            else:
                fail_count += 1
                print(f"  {progress} FAIL {slug:<25} {e}")
        except Exception as e:
            fail_count += 1
            print(f"  {progress} FAIL {slug:<25} {e}")

        # Save periodically (every 50 airfoils)
        if ok_count % 50 == 0 and ok_count > 0:
            save_db(rows, out_path)
            print(f"    ... saved checkpoint ({ok_count} airfoils)")

        time.sleep(delay)  # be polite to the server

    # Final save
    save_db(rows, out_path)

    print(f"\nDone!")
    print(f"  Airfoils with polars: {ok_count}")
    print(f"  Failed / no polar:    {fail_count}")
    if skip_count:
        print(f"  Skipped (resume):     {skip_count}")
    print(f"  Database saved to:    {out_path}")
    return out_path


# -------------------------------------------------------------------------
# Load DB for use in airfoil_sweep
# -------------------------------------------------------------------------
def load_db_as_tuples(db_path: str = None,
                      min_cl_max: float = 0.0,
                      max_cd_0: float = 1.0,
                      max_thickness: float = 100.0,
                      min_thickness: float = 0.0,
                      max_camber: float = 100.0) -> list[tuple]:
    """
    Load the database CSV and return tuples in airfoil_sweep format:
      (name, cl_0, cl_alpha_2d, cl_max, cd_0, cm_ac, cd_min, alpha_minD_deg, k_drag, thickness_pct)

    Drag polar fields are None if not available in the CSV.
    Optional filters for narrowing down candidates.
    """
    if db_path is None:
        db_path = DB_FILE

    results = []
    with open(db_path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cl_max = float(row["cl_max"])
            cd_0 = float(row["cd_0"])
            thick = float(row["thickness_pct"]) if row["thickness_pct"] else 0.0

            if cl_max < min_cl_max:
                continue
            if cd_0 > max_cd_0:
                continue
            if thick > max_thickness or thick < min_thickness:
                continue
            camber = float(row["camber_pct"]) if row["camber_pct"] else 0.0
            if camber > max_camber:
                continue

            cd_min = float(row["cd_min"]) if row.get("cd_min") else None
            alpha_minD_deg = float(row["alpha_minD_deg"]) if row.get("alpha_minD_deg") else None
            k_drag = float(row["k_drag"]) if row.get("k_drag") else None

            results.append((
                row["name"],
                float(row["cl_0"]),
                float(row["cl_alpha_2d"]),
                cl_max,
                cd_0,
                float(row["cm_ac"]),
                cd_min,
                alpha_minD_deg,
                k_drag,
                thick,
            ))

    return results


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build airfoil database from airfoiltools.com")
    parser.add_argument("--re", type=int, default=1000000,
                        help="Target Reynolds number (default: 1000000)")
    parser.add_argument("--letter", type=str, default=None,
                        help="Only scrape one letter (e.g. 'n' for NACA)")
    parser.add_argument("--filter", type=str, default=None,
                        help="Only include slugs/names containing this string")
    parser.add_argument("--resume", action="store_true",
                        help="Resume interrupted build (skip existing entries)")
    parser.add_argument("--delay", type=float, default=0.2,
                        help="Delay between requests in seconds (default: 0.2)")
    parser.add_argument("--sweep", action="store_true",
                        help="After building, run airfoil sweep with the database")
    parser.add_argument("--min-cl-max", type=float, default=0.0,
                        help="Filter: minimum cl_max for sweep (default: 0)")
    parser.add_argument("--max-cd-0", type=float, default=1.0,
                        help="Filter: maximum cd_0 for sweep (default: 1.0)")
    args = parser.parse_args()

    letters = [args.letter] if args.letter else None

    out_path = build_db(
        target_re=args.re,
        letters=letters,
        slug_filter=args.filter,
        resume=args.resume,
        delay=args.delay,
    )

    if args.sweep:
        print("\n" + "=" * 60)
        print("  Running airfoil sweep with database...")
        print("=" * 60)
        airfoils = load_db_as_tuples(
            db_path=out_path,
            min_cl_max=args.min_cl_max,
            max_cd_0=args.max_cd_0,
        )
        print(f"  {len(airfoils)} airfoils pass filters\n")
        from airfoil_sweep import run_sweep, score, print_table, BASE
        raw = run_sweep(airfoils, BASE)
        ranked = score(raw)
        print_table(ranked)
