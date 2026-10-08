"""Pinned ArduPilot JSON servo packet format (little endian, 16 or 32 PWM)."""

from dataclasses import dataclass
import json
import struct


@dataclass(frozen=True)
class ServoFrame:
    rate: int
    count: int
    pwm: tuple[int, ...]

    @classmethod
    def decode(cls, packet: bytes):
        channels = {40: (18458, 16), 72: (29569, 32)}.get(len(packet))
        if channels is None:
            raise ValueError("Servo packet must contain 40 or 72 bytes")
        magic, rate, count, *pwm = struct.unpack(f"<HHI{channels[1]}H", packet)
        if magic != channels[0] or not 1 <= rate <= 2000:
            raise ValueError("Invalid servo magic or frame rate")
        if any(value != 0 and not 800 <= value <= 2200 for value in pwm):
            raise ValueError("PWM outside supported range")
        return cls(rate, count, tuple(pwm))


def encode_sensor(sensor: dict) -> bytes:
    # allow_nan=False also checks nested arrays: firmware must never receive NaN.
    return ("\n" + json.dumps(sensor, allow_nan=False, separators=(",", ":")) + "\n").encode()
