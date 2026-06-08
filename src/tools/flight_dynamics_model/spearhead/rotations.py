"""Rotation utilities for 3-2-1 Euler angles."""

import numpy as np


def euler_body_to_ned(phi: float, theta: float, psi: float) -> np.ndarray:
    """Return the DCM that rotates body-frame vectors into NED."""
    cphi, sphi = np.cos(phi), np.sin(phi)
    ctheta, stheta = np.cos(theta), np.sin(theta)
    cpsi, spsi = np.cos(psi), np.sin(psi)

    return np.array(
        [
            [ctheta * cpsi, sphi * stheta * cpsi - cphi * spsi, cphi * stheta * cpsi + sphi * spsi],
            [ctheta * spsi, sphi * stheta * spsi + cphi * cpsi, cphi * stheta * spsi - sphi * cpsi],
            [-stheta, sphi * ctheta, cphi * ctheta],
        ]
    )


def body_rates_to_euler_rates(phi: float, theta: float, omega_b: np.ndarray) -> np.ndarray:
    """Convert body rates [p, q, r] to 3-2-1 Euler angle rates."""
    p, q, r = omega_b
    cphi, sphi = np.cos(phi), np.sin(phi)
    ctheta = np.cos(theta)

    if abs(ctheta) < 1e-6:
        raise ValueError("Euler angle rate conversion is singular near theta = +/- 90 deg")

    ttheta = np.tan(theta)
    return np.array(
        [
            p + sphi * ttheta * q + cphi * ttheta * r,
            cphi * q - sphi * r,
            (sphi * q + cphi * r) / ctheta,
        ]
    )

