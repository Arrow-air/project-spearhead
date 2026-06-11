import numpy as np

from spearhead.rotations import euler_body_to_ned


def test_euler_body_to_ned_zero_attitude_is_identity():
    R = euler_body_to_ned(0.0, 0.0, 0.0)
    np.testing.assert_allclose(R, np.eye(3), atol=1e-12)


def test_euler_body_to_ned_is_orthonormal():
    R = euler_body_to_ned(np.deg2rad(12.0), np.deg2rad(-5.0), np.deg2rad(33.0))
    np.testing.assert_allclose(R.T @ R, np.eye(3), atol=1e-12)
    np.testing.assert_allclose(np.linalg.det(R), 1.0, atol=1e-12)

