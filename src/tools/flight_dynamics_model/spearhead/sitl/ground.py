"""Provisional four-foot compliant gear on flat ground at NED z=0."""

from dataclasses import dataclass, asdict
import numpy as np

from ..rotations import euler_body_to_ned


@dataclass(frozen=True)
class GroundContact:
    half_length_m: float = .4
    half_width_m: float = .35
    height_m: float = .3
    stiffness_n_m: float = 20000.
    damping_ns_m: float = 650.
    friction: float = .5
    tangent_damping_ns_m: float = 100.

    def __post_init__(self):
        if not all(np.isfinite(v) and v > 0 for v in asdict(self).values()):
            raise ValueError("Ground-contact constants must be positive and finite")

    @property
    def feet(self):
        return np.array([[x*self.half_length_m, y*self.half_width_m, self.height_m]
                         for x, y in ((1, 1), (1, -1), (-1, 1), (-1, -1))])

    def rest_state(self, params):
        compression = params.mass*params.g/(4*self.stiffness_n_m)
        if compression >= self.height_m/2:
            raise ValueError("Ground gear cannot support aircraft weight")
        state = np.zeros(12)
        state[2] = -self.height_m + compression
        return state

    def loads(self, state):
        rotation = euler_body_to_ned(*state[6:9])
        feet = self.feet
        points = state[:3] + feet @ rotation.T
        velocities = (state[3:6] + np.cross(state[9:12], feet)) @ rotation.T
        normals = np.where(points[:, 2] >= 0,
                           np.maximum(0., self.stiffness_n_m*points[:, 2]+self.damping_ns_m*velocities[:, 2]), 0.)
        tangent = -self.tangent_damping_ns_m*velocities[:, :2]
        tangent *= np.minimum(1., self.friction*normals/np.maximum(np.linalg.norm(tangent, axis=1), 1e-12))[:, None]
        body_forces = np.column_stack((tangent, -normals)) @ rotation
        return body_forces.sum(axis=0), np.cross(feet, body_forces).sum(axis=0), {
            "contact": bool(np.any(normals > 0)), "normal_force_n": float(normals.sum()),
            "foot_normal_force_n": normals.tolist()}
