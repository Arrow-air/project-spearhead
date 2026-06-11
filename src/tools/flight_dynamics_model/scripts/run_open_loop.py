"""Run the default open-loop Spearhead simulation."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from spearhead.params import AircraftParams
from spearhead.force_moment import compute_force_moment_history
from spearhead.plotting import plot_force_moment_breakdown, plot_open_loop
from spearhead.simulation import run_open_loop


def main() -> None:
    params = AircraftParams()
    sol = run_open_loop(params=params)
    if not sol.success:
        raise RuntimeError(sol.message)
    plot_open_loop(sol, params)
    force_moment_history = compute_force_moment_history(
        sol.t,
        sol.y,
        params,
        control_override=getattr(sol, "control_override", None),
    )
    plot_force_moment_breakdown(force_moment_history)
    plt.show()


if __name__ == "__main__":
    main()
