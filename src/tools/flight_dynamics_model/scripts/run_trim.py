"""Run a preliminary longitudinal trim and trimmed open-loop simulation."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from spearhead.force_moment import compute_force_moment_history
from spearhead.params import AircraftParams
from spearhead.plotting import plot_force_moment_breakdown, plot_open_loop
from spearhead.simulation import run_open_loop
from spearhead.trim import trim_fixed_wing_longitudinal


def main() -> None:
    params = AircraftParams()
    trim = trim_fixed_wing_longitudinal(params, target_speed=22.0, altitude=100.0)
    total = trim["force_moment_breakdown"]["total"]
    force_b = total["force_b"]
    moment_b = total["moment_b"]

    print("Spearhead preliminary longitudinal trim")
    print(f"Success: {trim['success']}")
    print(f"Message: {trim['message']}")
    print(f"alpha:    {trim['alpha_deg']: .3f} deg")
    print(f"theta:    {trim['theta_deg']: .3f} deg")
    print(f"elevator: {trim['de_deg']: .3f} deg")
    print(f"throttle: {trim['throttle']: .3f}")
    print(f"Force body [N]:  Fx={force_b[0]: .3f}, Fy={force_b[1]: .3f}, Fz={force_b[2]: .3f}")
    print(f"Moment body [N*m]: Mx={moment_b[0]: .3f}, My={moment_b[1]: .3f}, Mz={moment_b[2]: .3f}")

    sol = run_open_loop(
        t_final=5.0,
        n_points=501,
        params=params,
        initial_state=trim["x_trim"],
        control_override=trim["control_trim"],
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    plot_open_loop(sol, params, control_override=trim["control_trim"])
    history = compute_force_moment_history(
        sol.t,
        sol.y,
        params,
        control_override=trim["control_trim"],
    )
    plot_force_moment_breakdown(history)
    plt.show()


if __name__ == "__main__":
    main()

