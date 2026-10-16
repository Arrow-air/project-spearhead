"""ArduPlane 4.7.1 JSON interface to the shared Spearhead physics engine."""

FIRMWARE_TAG = "Plane-4.7.1"
FIRMWARE_COMMIT = "dbe792162d06cab66c3475fd5556bf7a120f119e"
DATASET = "tail_update_2026_10_07"


# SITL home and air-data reference: Spearhead test field near Ankara, Turkey.
# Position and MSL altitude are taken from the PT1 DataFlash GPS (19 Sep 2026 logs).
SITL_HOME = (39.858200, 32.546000, 1125.0)  # latitude deg, longitude deg, altitude m AMSL


def sitl_home_argument(heading_deg=0):
    """Return the ArduPilot --home argument for the configured test field."""
    latitude, longitude, altitude = SITL_HOME
    return f"{latitude:.6f},{longitude:.6f},{altitude:g},{heading_deg:g}"


def ardupilot_isa_density(altitude_amsl_m):
    """Troposphere density exactly as AP_Baro::get_air_density_for_alt_amsl computes it."""
    radius = 6356766.0
    geopotential = radius*altitude_amsl_m/(radius + altitude_amsl_m)
    temperature = 288.15 - 0.0065*geopotential
    return 1.225*(temperature/288.15)**(9.80665/(0.0065*287.053072) - 1.0)


def ardupilot_density_ratio(height_above_home_m):
    """AHRS air-density ratio 4.7.1 uses for motor compensation at a height above SITL_HOME.

    AP_Motors multiplies VTOL throttle and attitude outputs by 1/ratio (density
    compensation), so hover thrust derived for the FDM must include this factor.
    """
    return ardupilot_isa_density(SITL_HOME[2] + height_above_home_m)/1.225
