"""
config.py — Single entry point for sizing and sweep tools.

Edit the parameters below, then run: python config.py

MODE:
  1 = Run airfoil sweep (uses AIRCRAFT + SWEEP dicts)
  2 = Run sizing for specific airfoils (uses AIRCRAFT + AIRFOIL dicts)
"""
import math

# ==========================================================================
# Run mode
# ==========================================================================
MODE = 1    # 1 = Airfoil sweep, 2 = Specific Airfoil

# ==========================================================================
# Aircraft parameters (used by both sizing tool and sweep tool)
# ==========================================================================
AIRCRAFT = dict(
    # --- Atmosphere & flight ---
    altitude_m        = 1200.0,       # m
    v_cruise          = 25.0,         # m/s
    v_stall           = 15.0,         # m/s
    MTOW_kg           = 25.0,         # kg

    # --- Fuselage ---
    fuselage_width    = 0.30,         # m (diameter)
    fus_length        = 1.00,         # m (nose tip to tail)
    nose_length       = 0.30,         # m (conical nose cone)

    # --- Wing geometry ---
    wing_AR           = 8.0,
    wing_incidence    = math.radians(1.0),   # rad
    wing_x_cg         = -0.15,        # m (negative = AC forward of CG)

    # --- Tail configuration ---
    # tail_config: 1 = V-tail / inverse V-tail, 2 = conventional
    tail_config       = 1,
    tail_proj_span    = 1.10,         # m, horizontal projected span (V-tail) or full span (conventional)
    tail_chord        = 0.28,         # m
    tail_x_cg         = 1.30,         # m (positive = AC aft of CG)
    tail_dihedral     = math.radians(35.0),  # rad (config 1 only)
    downwash_factor   = 0.2,          # 0 = no downwash, 1 = full downwash

    # --- Vertical fin (config 2 only, ignored for config 1) ---
    vert_chord        = None,         # m
    vert_span         = None,         # m (height)

    # --- Drag ---
    misc_drag_factor  = 1.3,          # multiplier on clean drag (gaps, booms, roughness)
    base_drag_recovery = 0.34,        # pusher prop wake ingestion (3D Hoerner, CFD-calibrated)

    # --- Pusher propeller thrust moment ---
    thrust_N          = 30.0,         # N, cruise thrust (0 = ignore)
    z_thrust_m        = 0.15,         # m, thrust line below CG (positive = below = nose-up)

    # --- Tail incidence override ---
    # Set to a value in degrees to skip trim solver and use fixed tail incidence.
    # Set to None for normal trim solve.
    tail_incidence_override_deg = None,
)

# ==========================================================================
# Airfoil selection (Used by sizing tool. Sweep only uses Re)
# ==========================================================================
AIRFOIL = dict(
    # --- Airfoil input mode ---
    # 1 = lookup by name from database, 2 = manual coefficients
    wing_airfoil_mode = 1,
    tail_airfoil_mode = 1,

    # Mode 1: airfoil name lookup (ignored when mode = 2)
    wing_airfoil      = "CLARK Z AIRFOIL",
    tail_airfoil      = "NACA 0018",
    re_db             = 500000,

    # Mode 2: manual aero coefficients (ignored when mode = 1)
    wing_cl_0         =  0.40,
    wing_cl_alpha_2d  =  6.28,       # rad^-1
    wing_cl_max       =  1.50,
    wing_cd_0         =  0.013,
    wing_cd_min       =  0.006,
    wing_alpha_minD_deg = 3.0,
    wing_k_drag       =  0.00005,
    wing_cm_ac        = -0.101,

    tail_cl_0         =  0.0,
    tail_cl_alpha_2d  =  6.28,
    tail_cd_0         =  0.010,
    tail_cd_min       =  0.005,
    tail_alpha_minD_deg = 0.0,
    tail_k_drag       =  0.00005,
    tail_cm_ac        =  0.0,
)

# ==========================================================================
# Sweep settings (used by airfoil_sweep.py only)
# ==========================================================================
SWEEP = dict(
    # --- Skin type ---
    oratex_skin       = True,         # exclude airfoils with concave lower surfaces

    # --- Excluded airfoils ---
    exclude = [
        "GOE 780 AIRFOIL",             # Hard to manufacture shape
        "HAM-STD HS1-620 AIRFOIL",     # Unreliable data
        "USNPS4 (smoothed)",           # Not conventional for chordwise tube
        "ARA-D 13% AIRFOIL",           # Designed for Re>10M, not reliable for low Re
        "MH 82  13.31%",               # Concave upper surface
        "TsAGI R-3a (12%)",            # Concave upper surface
        "S8035 for RC aerobatic  14% thick" # No significant advantage over NACA symmetric profiles
    ],

    # --- Wing filters ---
    min_cl_max        = 0.0,          # minimum 2D cl_max
    max_cd_0          = 1.0,          # maximum 2D cd_0
    min_thickness     = 0.0,          # minimum thickness %
    max_thickness     = 100.0,        # maximum thickness %
    min_abs_thickness = 0.05,         # minimum absolute wing thickness [m]

    # --- Tail filters ---
    tail_max_camber   = 1.0,          # maximum tail camber %
    tail_max_cd_0     = 0.015,        # maximum tail cd_0
    min_tail_abs_thickness = 0.038,   # minimum absolute tail thickness [m]

    # --- Stage control ---
    top_wing          = 30,           # top wing candidates for stage 2
    skip_tail         = False,        # skip stage 2 (wing only)

    # --- Scoring weights (must sum to 1.0) ---
    w_area            = 0.40,         # wing area (weight proxy)
    w_drag            = 0.40,         # total drag
    w_itail           = 0.10,         # |tail incidence| (closer to zero = better)
    w_stall           = 0.10,         # stall penalty (nose-up moment)

    # --- Output ---
    top               = None,         # show top N in table (None = all)
    detail            = 0,            # print full sizing for top N
    detail_rank       = [],           # print full sizing for specific ranks
)


# ==========================================================================
# Runner — python config.py
# ==========================================================================
if __name__ == "__main__":
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    if MODE == 1:
        print("Running airfoil sweep...\n")
        from airfoil_sweep import main as sweep_main
        sweep_main()

    elif MODE == 2:
        print("Running sizing for specific airfoils...\n")
        from aircraft_sizing import size_aircraft, print_results
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator

        results = size_aircraft(**AIRCRAFT, **AIRFOIL)
        print_results(results)

        alpha = results["sweep_alpha_deg"]
        L = results["sweep_L"]
        D = results["sweep_D"]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        ax1.plot(alpha, L, "b-o", markersize=4, linewidth=2)
        ax1.set_xlabel("AoA (deg)", fontsize=12)
        ax1.set_ylabel("Lift (N)", fontsize=12)
        ax1.set_title("Lift vs AoA")
        ax1.xaxis.set_major_locator(MultipleLocator(1))
        ax1.grid(True, alpha=0.3)

        ax2.plot(alpha, D, "r-o", markersize=4, linewidth=2)
        ax2.set_xlabel("AoA (deg)", fontsize=12)
        ax2.set_ylabel("Drag (N)", fontsize=12)
        ax2.set_title("Drag vs AoA")
        ax2.xaxis.set_major_locator(MultipleLocator(1))
        ax2.yaxis.set_major_locator(MultipleLocator(10))
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        out_dir = os.path.dirname(os.path.abspath(__file__))
        plot_path = os.path.join(out_dir, "LD_vs_alpha.png")
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Plot saved: {plot_path}")

    else:
        print(f"Unknown MODE={MODE}. Set MODE to 1 (sweep) or 2 (sizing).")
