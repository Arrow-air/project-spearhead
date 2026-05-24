import argparse
import csv
import math
import sys
import tkinter as tk
import webbrowser
from dataclasses import dataclass
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


if getattr(sys, "frozen", False):
    APP_DIR = Path(sys.executable).resolve().parent
    MODULE_DIR = Path(getattr(sys, "_MEIPASS", APP_DIR))
else:
    APP_DIR = Path(__file__).resolve().parent
    MODULE_DIR = APP_DIR


APP_NAME = "Nondimit"
APP_DESCRIPTION = "Aerodynamic coefficient conversion tool"
APP_CREDIT = "Developed by Alperen Gundogan for Arrow Air DAO"
APP_PROJECT = "Project Spearhead"
APP_WEBSITE = "https://arrowair.com/"

ZERO_VALUE = "<zero>"

BODY_FORCE_COEFF_COLUMNS = [
    "cfx",
    "cfy",
    "cfz",
]

WIND_FORCE_COEFF_COLUMNS = [
    "cd",
    "cs",
    "cl",
]

BODY_MOMENT_COEFF_COLUMNS = [
    "cmx",
    "cmy",
    "cmz",
]

WIND_MOMENT_COEFF_COLUMNS = [
    "cr",
    "cm",
    "cn",
]

BODY_FORCE_DIM_COLUMNS = [
    "Fx",
    "Fy",
    "Fz",
]

WIND_FORCE_DIM_COLUMNS = [
    "D",
    "S",
    "L",
]

BODY_MOMENT_DIM_COLUMNS = [
    "Mx",
    "My",
    "Mz",
]

WIND_MOMENT_DIM_COLUMNS = [
    "R",
    "M",
    "N",
]

FORCE_COEFF_COLUMNS = BODY_FORCE_COEFF_COLUMNS + WIND_FORCE_COEFF_COLUMNS
MOMENT_COEFF_COLUMNS = BODY_MOMENT_COEFF_COLUMNS + WIND_MOMENT_COEFF_COLUMNS
FORCE_DIM_COLUMNS = BODY_FORCE_DIM_COLUMNS + WIND_FORCE_DIM_COLUMNS
MOMENT_DIM_COLUMNS = BODY_MOMENT_DIM_COLUMNS + WIND_MOMENT_DIM_COLUMNS

BODY_FORCE_COLUMNS = [
    "cfx",
    "cfy",
    "cfz",
    "Fx",
    "Fy",
    "Fz",
]

BODY_MOMENT_COLUMNS = [
    "cmx",
    "cmy",
    "cmz",
    "Mx",
    "My",
    "Mz",
]

WIND_FORCE_COLUMNS = [
    "cd",
    "cs",
    "cl",
    "D",
    "S",
    "L",
]

WIND_MOMENT_COLUMNS = [
    "cr",
    "cm",
    "cn",
    "R",
    "M",
    "N",
]

OUTPUT_GROUPS = [
    ("body_force", "Body Forces (Fx, Fy, Fz)", BODY_FORCE_COLUMNS),
    ("body_moment", "Body Moments (Mx, My, Mz)", BODY_MOMENT_COLUMNS),
    ("wind_force", "Wind Forces (D, S, L)", WIND_FORCE_COLUMNS),
    ("wind_moment", "Wind Moments (R, M, N)", WIND_MOMENT_COLUMNS),
]

OUTPUT_GROUP_LABELS = {key: label for key, label, _columns in OUTPUT_GROUPS}
OUTPUT_GROUP_COLUMNS = {key: columns for key, _label, columns in OUTPUT_GROUPS}
DEFAULT_OUTPUT_GROUP_ORDER = [key for key, _label, _columns in OUTPUT_GROUPS]
LEGACY_WIND_OUTPUT_COLUMNS = [
    "CD",
    "CS",
    "CL",
    "CR",
    "CM",
    "CN",
    "D_dim",
    "S_dim",
    "L_dim",
    "R_dim",
    "M_dim",
    "N_dim",
    "CFx_wind",
    "CFy_wind",
    "CFz_wind",
    "CMx_wind",
    "CMy_wind",
    "CMz_wind",
    "Fx_wind_dim",
    "Fy_wind_dim",
    "Fz_wind_dim",
    "Mx_wind_dim",
    "My_wind_dim",
    "Mz_wind_dim",
]
LEGACY_BODY_COEFF_COLUMNS = [
    "CFx",
    "CFy",
    "CFz",
    "CMx",
    "CMy",
    "CMz",
    "CFx_body",
    "CFy_body",
    "CFz_body",
    "CMx_body",
    "CMy_body",
    "CMz_body",
]
LEGACY_BODY_DIM_COLUMNS = [
    "Fx_body_dim",
    "Fy_body_dim",
    "Fz_body_dim",
    "Mx_body_dim",
    "My_body_dim",
    "Mz_body_dim",
]
GENERATED_OUTPUT_COLUMNS = set(
    FORCE_COEFF_COLUMNS
    + MOMENT_COEFF_COLUMNS
    + FORCE_DIM_COLUMNS
    + MOMENT_DIM_COLUMNS
    + LEGACY_WIND_OUTPUT_COLUMNS
    + LEGACY_BODY_COEFF_COLUMNS
    + LEGACY_BODY_DIM_COLUMNS
)
BETA_ZERO_TOLERANCE = 1e-12
ZERO_AT_BETA_COMPONENTS = [
    ("fy", "Fy"),
    ("mx", "Mx"),
    ("mz", "Mz"),
]
ZERO_AT_BETA_LABELS = {key: label for key, label in ZERO_AT_BETA_COMPONENTS}


def normalized(name):
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def detect_column(headers, exact=(), contains=()):
    clean = {normalized(header): header for header in headers}
    for item in exact:
        if item in clean:
            return clean[item]
    for header in headers:
        text = normalized(header)
        if any(item in text for item in contains):
            return header
    return ""


def read_csv_table_with_metadata(path):
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    lines = text.splitlines()
    metadata = []
    while lines and lines[0].lstrip().startswith("#"):
        metadata.append(next(csv.reader([lines[0]])))
        lines = lines[1:]
    sample = "\n".join(lines[:40])
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    rows = list(csv.DictReader(lines, dialect=dialect))
    if not rows:
        raise ValueError("input table has no rows")
    headers = list(rows[0].keys())
    if not headers:
        raise ValueError("input table has no header")
    return metadata, headers, rows


def read_csv_table(path):
    _metadata, headers, rows = read_csv_table_with_metadata(path)
    return headers, rows


def write_csv_table(path, headers, rows, metadata=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        if metadata:
            writer = csv.writer(handle)
            first = metadata[0] if isinstance(metadata, list) and metadata else None
            if isinstance(first, (list, tuple)):
                writer.writerows(metadata)
            else:
                writer.writerow(metadata)
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_float(value, name):
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} is not numeric: {value!r}") from None


def parse_metadata_row(metadata):
    parsed = {}
    if metadata and all(isinstance(item, str) for item in metadata):
        rows = [metadata]
    else:
        rows = metadata or []

    for row in rows:
        cells = [str(item).strip() for item in row]
        if len(cells) >= 3 and cells[0] == "#":
            key = cells[1]
            value = cells[2]
            if key and normalized(key) not in {"field", "parameter", "key"}:
                parsed[normalized(key)] = value.strip()
            continue
        for text in cells:
            if "=" not in text:
                continue
            key, value = text.split("=", 1)
            parsed[normalized(key)] = value.strip()
    return parsed


def metadata_value(metadata, keys, default=""):
    for key in keys:
        value = metadata.get(normalized(key))
        if value not in (None, ""):
            return value
    return default


def metadata_float(metadata, keys, name, default=None):
    value = metadata_value(metadata, keys, default=None)
    if value is None:
        if default is not None:
            return default
        raise ValueError(f"{name} is missing from export metadata")
    return as_float(value, name)


def parse_vector_text(value, name):
    text = str(value).strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if len(parts) != 3:
        raise ValueError(f"{name} must contain three comma-separated values")
    return tuple(as_float(part, f"{name} {axis}") for part, axis in zip(parts, "xyz"))


def metadata_vector(metadata, keys, name, default=None):
    value = metadata_value(metadata, keys, default=None)
    if value is None:
        if default is not None:
            return default
        raise ValueError(f"{name} is missing from export metadata")
    return parse_vector_text(value, name)


def row_value(row, column, name, default=None):
    if column == ZERO_VALUE:
        return default
    if not column:
        if default is None:
            raise ValueError(f"{name} column is not selected")
        return default
    return as_float(row.get(column, ""), name)


def decimals_from_accuracy(accuracy):
    if accuracy <= 0:
        raise ValueError("output accuracy must be positive")
    text = f"{accuracy:.12f}".rstrip("0")
    if "." not in text:
        return 0
    return len(text.split(".", 1)[1])


def format_output_value(value, accuracy):
    decimals = decimals_from_accuracy(accuracy)
    rounded = round(float(value) / accuracy) * accuracy
    return f"{rounded:.{decimals}f}"


def format_numeric_cells(row, accuracy):
    formatted = {}
    for key, value in row.items():
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            formatted[key] = value
            continue
        if not math.isfinite(numeric):
            formatted[key] = value
            continue
        formatted[key] = format_output_value(numeric, accuracy)
    return formatted


def wind_axes(alpha_deg, beta_deg):
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    ca = math.cos(alpha)
    sa = math.sin(alpha)
    cb = math.cos(beta)
    sb = math.sin(beta)

    x_wind = (-ca * cb, -sb, -sa * cb)
    y_wind = (-ca * sb, cb, -sa * sb)
    z_wind = (sa, 0.0, -ca)
    return x_wind, y_wind, z_wind


def dot(vector, axis):
    return vector[0] * axis[0] + vector[1] * axis[1] + vector[2] * axis[2]


@dataclass
class ConversionSettings:
    rho: float
    area: float
    velocity: float
    span: float
    chord: float
    output_accuracy: float
    old_center: tuple
    new_center: tuple

    @property
    def q(self):
        return 0.5 * self.rho * self.velocity * self.velocity

    @property
    def force_denominator(self):
        return self.q * self.area

    @property
    def roll_yaw_denominator(self):
        return self.force_denominator * self.span

    @property
    def pitch_denominator(self):
        return self.force_denominator * self.chord

    def metadata_row(self, mapping, output_order=None, zero_at_beta=None):
        new = ",".join(f"{value:.10g}" for value in self.new_center)
        order = ";".join(normalized_output_order(output_order))
        zero_options = ";".join(normalized_zero_at_beta_options(zero_at_beta))
        return [
            [
                "# Nondimit export",
                f"rho={self.rho:.10g}",
                f"S={self.area:.10g}",
                f"V={self.velocity:.10g}",
                f"q={self.q:.10g}",
                f"span={self.span:.10g}",
                f"chord={self.chord:.10g}",
                f"output_accuracy={self.output_accuracy:.10g}",
            ],
            [
                "# Moment center",
                f"moment_center_body_xyz=({new})",
                "body_axes=x_forward_y_right_z_down",
            ],
            [
                "# Source and options",
                f"alpha_column={mapping.get('alpha', '')}",
                f"beta_column={mapping.get('beta', '')}",
                f"output_group_order={order}",
                f"zero_at_beta0={zero_options}",
            ],
        ]


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def shifted_moment(moment, force, old_center, new_center):
    r_new_to_old = tuple(old - new for old, new in zip(old_center, new_center))
    arm_cross_force = cross(r_new_to_old, force)
    return tuple(m + delta for m, delta in zip(moment, arm_cross_force))


def append_unique(headers, names):
    out = list(headers)
    for name in names:
        if name not in out:
            out.append(name)
    return out


def normalized_output_order(output_order=None):
    order = list(output_order or DEFAULT_OUTPUT_GROUP_ORDER)
    clean = []
    for key in order:
        if key in OUTPUT_GROUP_COLUMNS and key not in clean:
            clean.append(key)
    for key in DEFAULT_OUTPUT_GROUP_ORDER:
        if key not in clean:
            clean.append(key)
    return clean


def parse_output_order_text(value):
    if not value:
        return list(DEFAULT_OUTPUT_GROUP_ORDER)
    aliases = {normalized(key): key for key in OUTPUT_GROUP_COLUMNS}
    aliases.update({normalized(label): key for key, label in OUTPUT_GROUP_LABELS.items()})
    order = []
    for raw in str(value).replace(",", ";").split(";"):
        key = aliases.get(normalized(raw))
        if key:
            order.append(key)
    return normalized_output_order(order)


def normalized_zero_at_beta_options(zero_options=None):
    if not zero_options:
        return []
    if isinstance(zero_options, dict):
        raw_options = [key for key, enabled in zero_options.items() if enabled]
    else:
        raw_options = list(zero_options)

    aliases = {normalized(key): key for key, _label in ZERO_AT_BETA_COMPONENTS}
    aliases.update({normalized(label): key for key, label in ZERO_AT_BETA_LABELS.items()})
    clean = []
    for raw in raw_options:
        key = aliases.get(normalized(raw))
        if key and key not in clean:
            clean.append(key)
    return clean


def parse_zero_at_beta_text(value):
    if not value:
        return []
    return normalized_zero_at_beta_options(str(value).replace(",", ";").split(";"))


def output_headers_for(headers, include_moments, output_order=None):
    out = [header for header in headers if header not in GENERATED_OUTPUT_COLUMNS]
    for group_key in normalized_output_order(output_order):
        if "moment" in group_key and not include_moments:
            continue
        out = append_unique(out, OUTPUT_GROUP_COLUMNS[group_key])
    return out


def validate_settings(settings):
    if settings.rho <= 0:
        raise ValueError("density must be positive")
    if settings.area <= 0:
        raise ValueError("surface area must be positive")
    if settings.velocity <= 0:
        raise ValueError("velocity must be positive")
    if settings.span <= 0 or settings.chord <= 0:
        raise ValueError("moment reference lengths must be positive")
    if settings.output_accuracy <= 0:
        raise ValueError("output accuracy must be positive")


def read_export_settings(metadata):
    parsed = parse_metadata_row(metadata)
    current_center = metadata_vector(
        parsed,
        ("moment_center_body_xyz", "new_moment_center_body_xyz", "old_moment_center_body_xyz"),
        "current moment center",
    )
    return parsed, ConversionSettings(
        rho=metadata_float(parsed, ("rho", "density"), "density"),
        area=metadata_float(parsed, ("s", "area", "surface_area"), "surface area"),
        velocity=metadata_float(parsed, ("v", "velocity"), "velocity"),
        span=metadata_float(parsed, ("span",), "span"),
        chord=metadata_float(parsed, ("chord",), "chord"),
        output_accuracy=metadata_float(parsed, ("output_accuracy",), "output accuracy", default=0.0001),
        old_center=current_center,
        new_center=current_center,
    )


def vector_from_columns(row, columns, name, row_number):
    return tuple(as_float(row.get(column, ""), f"{name} {column} at row {row_number}") for column in columns)


def body_force_from_export(row, row_number, settings):
    for dim_columns in (("Fx_body_dim", "Fy_body_dim", "Fz_body_dim"), ("Fx", "Fy", "Fz")):
        if all(column in row and row.get(column, "") != "" for column in dim_columns):
            return vector_from_columns(row, dim_columns, "body force", row_number)
    for coeff_columns in (("cfx", "cfy", "cfz"), ("CFx", "CFy", "CFz"), ("CFx_body", "CFy_body", "CFz_body")):
        if all(column in row and row.get(column, "") != "" for column in coeff_columns):
            coeff = vector_from_columns(row, coeff_columns, "body force coefficient", row_number)
            return tuple(value * settings.force_denominator for value in coeff)
    raise ValueError(f"body force columns are missing at row {row_number}")


def body_moment_from_export(row, row_number, settings):
    for dim_columns in (("Mx_body_dim", "My_body_dim", "Mz_body_dim"), ("Mx", "My", "Mz")):
        if all(column in row and row.get(column, "") != "" for column in dim_columns):
            return vector_from_columns(row, dim_columns, "body moment", row_number)
    for coeff_columns in (("cmx", "cmy", "cmz"), ("CMx", "CMy", "CMz"), ("CMx_body", "CMy_body", "CMz_body")):
        if all(column in row and row.get(column, "") != "" for column in coeff_columns):
            cmx, cmy, cmz = vector_from_columns(row, coeff_columns, "body moment coefficient", row_number)
            return (
                cmx * settings.roll_yaw_denominator,
                cmy * settings.pitch_denominator,
                cmz * settings.roll_yaw_denominator,
            )
    raise ValueError(f"body moment columns are missing at row {row_number}")


def source_component_columns(headers):
    return {
        "fy": detect_column(headers, exact=("fy", "forcey"), contains=("forcey", "yforce")),
        "mx": detect_column(headers, exact=("mx", "momentx"), contains=("momentx", "xmoment")),
        "mz": detect_column(headers, exact=("mz", "momentz"), contains=("momentz", "zmoment")),
    }


def zero_source_cell(out, column):
    if column and column != ZERO_VALUE and column in out:
        out[column] = 0.0


def apply_zero_at_beta(out, beta, f_body, m_body=None, zero_options=None, source_columns=None):
    options = normalized_zero_at_beta_options(zero_options)
    if not options or abs(beta) > BETA_ZERO_TOLERANCE:
        return f_body, m_body

    source_columns = source_columns or {}
    if "fy" in options:
        f_body = (f_body[0], 0.0, f_body[2])
        zero_source_cell(out, source_columns.get("fy"))

    if m_body is not None:
        mx, my, mz = m_body
        if "mx" in options:
            mx = 0.0
            zero_source_cell(out, source_columns.get("mx"))
        if "mz" in options:
            mz = 0.0
            zero_source_cell(out, source_columns.get("mz"))
        m_body = (mx, my, mz)

    return f_body, m_body


def write_force_and_moment_columns(out, f_body, m_body, alpha, beta, settings):
    force_den = settings.force_denominator
    roll_yaw_den = settings.roll_yaw_denominator
    pitch_den = settings.pitch_denominator
    x_axis, y_axis, z_axis = wind_axes(alpha, beta)
    f_wind = (dot(f_body, x_axis), dot(f_body, y_axis), dot(f_body, z_axis))
    m_wind = (dot(m_body, x_axis), dot(m_body, y_axis), dot(m_body, z_axis))

    out["cfx"] = f_body[0] / force_den
    out["cfy"] = f_body[1] / force_den
    out["cfz"] = f_body[2] / force_den
    out["cd"] = f_wind[0] / force_den
    out["cs"] = f_wind[1] / force_den
    out["cl"] = f_wind[2] / force_den
    out["Fx"] = f_body[0]
    out["Fy"] = f_body[1]
    out["Fz"] = f_body[2]
    out["D"] = f_wind[0]
    out["S"] = f_wind[1]
    out["L"] = f_wind[2]

    out["cmx"] = m_body[0] / roll_yaw_den
    out["cmy"] = m_body[1] / pitch_den
    out["cmz"] = m_body[2] / roll_yaw_den
    out["cr"] = m_wind[0] / roll_yaw_den
    out["cm"] = m_wind[1] / pitch_den
    out["cn"] = m_wind[2] / roll_yaw_den
    out["Mx"] = m_body[0]
    out["My"] = m_body[1]
    out["Mz"] = m_body[2]
    out["R"] = m_wind[0]
    out["M"] = m_wind[1]
    out["N"] = m_wind[2]


def convert_rows(rows, headers, mapping, settings, output_order=None, zero_at_beta=None):
    validate_settings(settings)

    required = ["alpha", "fx", "fy", "fz"]
    missing = [name for name in required if not mapping.get(name)]
    if missing:
        raise ValueError("missing required mapping: " + ", ".join(missing))

    have_moments = all(mapping.get(name) for name in ["mx", "my", "mz"])
    output_headers = output_headers_for(headers, have_moments, output_order)

    converted = []
    force_den = settings.force_denominator
    roll_yaw_den = settings.roll_yaw_denominator
    pitch_den = settings.pitch_denominator

    for index, row in enumerate(rows, start=2):
        out = dict(row)
        alpha = row_value(row, mapping.get("alpha"), f"alpha at row {index}")
        beta = row_value(row, mapping.get("beta", ZERO_VALUE), f"beta at row {index}", default=0.0)
        fx = row_value(row, mapping.get("fx"), f"Fx at row {index}")
        fy = row_value(row, mapping.get("fy"), f"Fy at row {index}")
        fz = row_value(row, mapping.get("fz"), f"Fz at row {index}")

        source_columns = {key: mapping.get(key, "") for key in ("fy", "mx", "mz")}
        f_body = (fx, fy, fz)
        f_body, _m_body = apply_zero_at_beta(out, beta, f_body, zero_options=zero_at_beta, source_columns=source_columns)
        x_axis, y_axis, z_axis = wind_axes(alpha, beta)
        f_wind = (dot(f_body, x_axis), dot(f_body, y_axis), dot(f_body, z_axis))

        out["cfx"] = f_body[0] / force_den
        out["cfy"] = f_body[1] / force_den
        out["cfz"] = f_body[2] / force_den
        out["cd"] = f_wind[0] / force_den
        out["cs"] = f_wind[1] / force_den
        out["cl"] = f_wind[2] / force_den
        out["Fx"] = f_body[0]
        out["Fy"] = f_body[1]
        out["Fz"] = f_body[2]
        out["D"] = f_wind[0]
        out["S"] = f_wind[1]
        out["L"] = f_wind[2]

        if have_moments:
            mx = row_value(row, mapping.get("mx"), f"Mx at row {index}")
            my = row_value(row, mapping.get("my"), f"My at row {index}")
            mz = row_value(row, mapping.get("mz"), f"Mz at row {index}")
            m_body = shifted_moment((mx, my, mz), f_body, settings.old_center, settings.new_center)
            _f_body, m_body = apply_zero_at_beta(
                out,
                beta,
                f_body,
                m_body,
                zero_options=zero_at_beta,
                source_columns=source_columns,
            )
            m_wind = (dot(m_body, x_axis), dot(m_body, y_axis), dot(m_body, z_axis))

            out["cmx"] = m_body[0] / roll_yaw_den
            out["cmy"] = m_body[1] / pitch_den
            out["cmz"] = m_body[2] / roll_yaw_den
            out["cr"] = m_wind[0] / roll_yaw_den
            out["cm"] = m_wind[1] / pitch_den
            out["cn"] = m_wind[2] / roll_yaw_den
            out["Mx"] = m_body[0]
            out["My"] = m_body[1]
            out["Mz"] = m_body[2]
            out["R"] = m_wind[0]
            out["M"] = m_wind[1]
            out["N"] = m_wind[2]

        converted.append(format_numeric_cells(out, settings.output_accuracy))

    return output_headers, converted


def recenter_export_rows(rows, headers, mapping, settings, output_order=None, zero_at_beta=None):
    validate_settings(settings)
    if not mapping.get("alpha"):
        raise ValueError("alpha column is not selected")

    output_headers = output_headers_for(headers, include_moments=True, output_order=output_order)
    converted = []
    source_columns = source_component_columns(headers)

    for index, row in enumerate(rows, start=3):
        out = dict(row)
        alpha = row_value(row, mapping.get("alpha"), f"alpha at row {index}")
        beta = row_value(row, mapping.get("beta", ZERO_VALUE), f"beta at row {index}", default=0.0)
        f_body = body_force_from_export(row, index, settings)
        f_body, _m_body = apply_zero_at_beta(out, beta, f_body, zero_options=zero_at_beta, source_columns=source_columns)
        m_current = body_moment_from_export(row, index, settings)
        m_body = shifted_moment(m_current, f_body, settings.old_center, settings.new_center)
        _f_body, m_body = apply_zero_at_beta(
            out,
            beta,
            f_body,
            m_body,
            zero_options=zero_at_beta,
            source_columns=source_columns,
        )

        write_force_and_moment_columns(out, f_body, m_body, alpha, beta, settings)
        converted.append(format_numeric_cells(out, settings.output_accuracy))

    return output_headers, converted


def default_mapping(headers):
    return {
        "alpha": detect_column(headers, exact=("alpha", "alphadeg", "aoa", "aoadeg"), contains=("alpha", "aoa")),
        "beta": detect_column(headers, exact=("beta", "betadeg", "sideslip", "sideslipdeg"), contains=("beta", "sideslip")) or ZERO_VALUE,
        "fx": detect_column(headers, exact=("fx", "forcex"), contains=("fx", "forcex", "xforce")),
        "fy": detect_column(headers, exact=("fy", "forcey"), contains=("fy", "forcey", "yforce")),
        "fz": detect_column(headers, exact=("fz", "forcez"), contains=("fz", "forcez", "zforce")),
        "mx": detect_column(headers, exact=("mx", "momentx"), contains=("mx", "momentx", "xmoment")),
        "my": detect_column(headers, exact=("my", "momenty"), contains=("my", "momenty", "ymoment")),
        "mz": detect_column(headers, exact=("mz", "momentz"), contains=("mz", "momentz", "zmoment")),
    }


def load_user_guide():
    for base in [APP_DIR, MODULE_DIR]:
        path = base / "README.md"
        if path.exists():
            try:
                return path.read_text(encoding="utf-8")
            except OSError:
                pass
    return (
        "# Nondimit User Guide\n\n"
        "The packaged README.md could not be loaded. Keep README.md next to "
        "Nondimit.exe, or rebuild the executable with README.md included."
    )


def add_scrolled_text(parent, text, wrap=tk.WORD):
    frame = ttk.Frame(parent)
    frame.pack(fill=tk.BOTH, expand=True)
    text_widget = tk.Text(frame, wrap=wrap, font=("Consolas", 9), height=24)
    yscroll = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=text_widget.yview)
    xscroll = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=text_widget.xview)
    text_widget.configure(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set)
    text_widget.grid(row=0, column=0, sticky="nsew")
    yscroll.grid(row=0, column=1, sticky="ns")
    xscroll.grid(row=1, column=0, sticky="ew")
    frame.rowconfigure(0, weight=1)
    frame.columnconfigure(0, weight=1)
    text_widget.insert(tk.END, text)
    text_widget.configure(state=tk.DISABLED)
    return text_widget


class Nondimit:
    def __init__(self, root):
        self.root = root
        self.root.title(APP_NAME)
        self.root.geometry("1180x760")
        self.root.minsize(980, 620)

        self.headers = []
        self.rows = []
        self.converted_headers = []
        self.converted_rows = []
        self.recenter_metadata = []
        self.recenter_headers = []
        self.recenter_rows = []
        self.recentered_headers = []
        self.recentered_rows = []

        self.input_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.recenter_input_var = tk.StringVar()
        self.recenter_output_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Idle")
        self.parameter_vars = {
            "rho": tk.StringVar(value="1.225"),
            "area": tk.StringVar(value="1.0"),
            "velocity": tk.StringVar(value="25.0"),
            "span": tk.StringVar(value="1.0"),
            "chord": tk.StringVar(value="1.0"),
            "output_accuracy": tk.StringVar(value="0.0001"),
            "old_x": tk.StringVar(value="0.0"),
            "old_y": tk.StringVar(value="0.0"),
            "old_z": tk.StringVar(value="0.0"),
            "new_x": tk.StringVar(value="0.0"),
            "new_y": tk.StringVar(value="0.0"),
            "new_z": tk.StringVar(value="0.0"),
        }
        self.recenter_parameter_vars = {
            "rho": tk.StringVar(value="1.225"),
            "area": tk.StringVar(value="1.0"),
            "velocity": tk.StringVar(value="25.0"),
            "span": tk.StringVar(value="1.0"),
            "chord": tk.StringVar(value="1.0"),
            "output_accuracy": tk.StringVar(value="0.0001"),
            "current_x": tk.StringVar(value="0.0"),
            "current_y": tk.StringVar(value="0.0"),
            "current_z": tk.StringVar(value="0.0"),
            "new_x": tk.StringVar(value="0.0"),
            "new_y": tk.StringVar(value="0.0"),
            "new_z": tk.StringVar(value="0.0"),
        }
        self.mapping_vars = {
            key: tk.StringVar()
            for key in ["alpha", "beta", "fx", "fy", "fz", "mx", "my", "mz"]
        }
        self.mapping_boxes = {}
        self.recenter_mapping_vars = {key: tk.StringVar() for key in ["alpha", "beta"]}
        self.recenter_mapping_boxes = {}
        self.output_orders = {
            "convert": list(DEFAULT_OUTPUT_GROUP_ORDER),
            "recenter": list(DEFAULT_OUTPUT_GROUP_ORDER),
        }
        self.output_order_boxes = {}
        self.zero_beta_vars = {
            key: {component: tk.BooleanVar(value=False) for component, _label in ZERO_AT_BETA_COMPONENTS}
            for key in ["convert", "recenter"]
        }
        self.preview = None
        self.recenter_preview = None

        self.build_ui()

    def build_menu(self):
        menubar = tk.Menu(self.root)
        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)
        self.root.config(menu=menubar)

    def show_about(self):
        dialog = tk.Toplevel(self.root)
        dialog.title(f"About {APP_NAME}")
        dialog.geometry("760x620")
        dialog.minsize(560, 420)
        dialog.transient(self.root)

        notebook = ttk.Notebook(dialog)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        about_tab = ttk.Frame(notebook, padding=22)
        guide_tab = ttk.Frame(notebook, padding=10)
        notebook.add(about_tab, text="About")
        notebook.add(guide_tab, text="User Guide")

        ttk.Label(about_tab, text=APP_NAME, font=("Segoe UI", 18, "bold")).pack(anchor=tk.W)
        ttk.Label(about_tab, text=APP_DESCRIPTION, font=("Segoe UI", 11)).pack(anchor=tk.W, pady=(4, 18))
        ttk.Label(about_tab, text=APP_CREDIT, font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        ttk.Label(about_tab, text=APP_PROJECT).pack(anchor=tk.W, pady=(0, 6))
        ttk.Label(about_tab, text=APP_WEBSITE).pack(anchor=tk.W, pady=(0, 18))
        ttk.Button(about_tab, text="Open Website", command=lambda: webbrowser.open(APP_WEBSITE)).pack(anchor=tk.W)

        add_scrolled_text(guide_tab, load_user_guide(), wrap=tk.NONE)

        ttk.Button(dialog, text="Close", command=dialog.destroy).pack(anchor=tk.E, padx=10, pady=(0, 10))
        dialog.focus_set()

    def build_ui(self):
        self.build_menu()
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TFrame", background="#f7f8fa")
        style.configure("Header.TFrame", background="#ffffff")
        style.configure("Header.TLabel", background="#ffffff", foreground="#111827", font=("Segoe UI", 16, "bold"))
        style.configure("TButton", font=("Segoe UI", 9))

        header = ttk.Frame(self.root, style="Header.TFrame", padding=(14, 10))
        header.pack(fill=tk.X)
        ttk.Label(header, text="Aero Coefficient Converter", style="Header.TLabel").pack(side=tk.LEFT)

        notebook = ttk.Notebook(self.root)
        notebook.pack(fill=tk.BOTH, expand=True)

        convert_tab = ttk.Frame(notebook)
        recenter_tab = ttk.Frame(notebook)
        notebook.add(convert_tab, text="Convert Raw Table")
        notebook.add(recenter_tab, text="Modify Export")

        self.build_convert_tab(convert_tab)
        self.build_recenter_tab(recenter_tab)

        status = ttk.Label(self.root, textvariable=self.status_var, anchor=tk.W, padding=(10, 5))
        status.pack(fill=tk.X)

    def build_convert_tab(self, parent):
        body = ttk.PanedWindow(parent, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(body, padding=12)
        right = ttk.Frame(body, padding=12)
        body.add(left, weight=0)
        body.add(right, weight=1)

        controls = self.scrollable_frame(left)
        self.build_controls(controls)
        self.preview = self.build_preview(right)

    def build_recenter_tab(self, parent):
        body = ttk.PanedWindow(parent, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(body, padding=12)
        right = ttk.Frame(body, padding=12)
        body.add(left, weight=0)
        body.add(right, weight=1)

        controls = self.scrollable_frame(left)
        self.build_recenter_controls(controls)
        self.recenter_preview = self.build_preview(right)

    def scrollable_frame(self, parent):
        container = ttk.Frame(parent)
        container.pack(fill=tk.BOTH, expand=True)

        canvas = tk.Canvas(container, highlightthickness=0, background="#f7f8fa", width=430)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        content = ttk.Frame(canvas)
        window = canvas.create_window((0, 0), window=content, anchor=tk.NW)

        def update_scroll_region(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))

        def update_content_width(event):
            canvas.itemconfigure(window, width=event.width)

        def on_mousewheel(event):
            canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        content.bind("<Configure>", update_scroll_region)
        canvas.bind("<Configure>", update_content_width)
        canvas.bind("<Enter>", lambda _event: canvas.bind_all("<MouseWheel>", on_mousewheel))
        canvas.bind("<Leave>", lambda _event: canvas.unbind_all("<MouseWheel>"))
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)
        return content

    def build_controls(self, parent):
        ttk.Label(parent, text="Files", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.path_row(parent, "Input", self.input_var, self.browse_input)
        self.path_row(parent, "Output", self.output_var, self.browse_output)
        ttk.Button(parent, text="Load Table", command=self.load_table).pack(fill=tk.X, pady=(6, 12))

        ttk.Separator(parent).pack(fill=tk.X, pady=8)
        ttk.Label(parent, text="Reference Values", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.number_row(parent, "Density", self.parameter_vars["rho"])
        self.number_row(parent, "Surface area", self.parameter_vars["area"])
        self.number_row(parent, "Velocity", self.parameter_vars["velocity"])
        self.number_row(parent, "Span (Mx/Mz)", self.parameter_vars["span"])
        self.number_row(parent, "Chord (My)", self.parameter_vars["chord"])
        self.number_row(parent, "Output accuracy", self.parameter_vars["output_accuracy"])

        ttk.Label(parent, text="Moment Center", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(12, 6))
        self.vector_row(
            parent,
            "Old",
            self.parameter_vars["old_x"],
            self.parameter_vars["old_y"],
            self.parameter_vars["old_z"],
        )
        self.vector_row(
            parent,
            "New",
            self.parameter_vars["new_x"],
            self.parameter_vars["new_y"],
            self.parameter_vars["new_z"],
        )

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Column Mapping", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        for label, key in [
            ("Alpha", "alpha"),
            ("Beta", "beta"),
            ("Fx", "fx"),
            ("Fy", "fy"),
            ("Fz", "fz"),
            ("Mx", "mx"),
            ("My", "my"),
            ("Mz", "mz"),
        ]:
            self.mapping_row(parent, label, key)

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Beta = 0 Cleanup", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.zero_beta_control(parent, "convert")

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Output Group Order", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.output_order_control(parent, "convert")

        ttk.Button(parent, text="Convert and Save", command=self.convert_and_save).pack(fill=tk.X, pady=(14, 6))

    def build_recenter_controls(self, parent):
        ttk.Label(parent, text="Files", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.path_row(parent, "Input", self.recenter_input_var, self.browse_recenter_input)
        self.path_row(parent, "Output", self.recenter_output_var, self.browse_recenter_output)
        ttk.Button(parent, text="Load Export", command=self.load_recenter_export).pack(fill=tk.X, pady=(6, 12))

        ttk.Separator(parent).pack(fill=tk.X, pady=8)
        ttk.Label(parent, text="Reference Values", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.number_row(parent, "Density", self.recenter_parameter_vars["rho"])
        self.number_row(parent, "Surface area", self.recenter_parameter_vars["area"])
        self.number_row(parent, "Velocity", self.recenter_parameter_vars["velocity"])
        self.number_row(parent, "Span (Mx/Mz)", self.recenter_parameter_vars["span"])
        self.number_row(parent, "Chord (My)", self.recenter_parameter_vars["chord"])
        self.number_row(parent, "Output accuracy", self.recenter_parameter_vars["output_accuracy"])

        ttk.Label(parent, text="Moment Center", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(12, 6))
        self.vector_row(
            parent,
            "Current",
            self.recenter_parameter_vars["current_x"],
            self.recenter_parameter_vars["current_y"],
            self.recenter_parameter_vars["current_z"],
            state="readonly",
            label_width=8,
        )
        self.vector_row(
            parent,
            "New",
            self.recenter_parameter_vars["new_x"],
            self.recenter_parameter_vars["new_y"],
            self.recenter_parameter_vars["new_z"],
            label_width=8,
        )

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Angle Columns", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        for label, key in [("Alpha", "alpha"), ("Beta", "beta")]:
            self.mapping_row(parent, label, key, vars_map=self.recenter_mapping_vars, boxes_map=self.recenter_mapping_boxes)

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Beta = 0 Cleanup", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.zero_beta_control(parent, "recenter")

        ttk.Separator(parent).pack(fill=tk.X, pady=12)
        ttk.Label(parent, text="Output Group Order", font=("Segoe UI", 10, "bold")).pack(anchor=tk.W, pady=(0, 6))
        self.output_order_control(parent, "recenter")

        ttk.Button(parent, text="Modify and Save", command=self.recenter_and_save).pack(fill=tk.X, pady=(14, 6))

    def path_row(self, parent, label, variable, command):
        ttk.Label(parent, text=label).pack(anchor=tk.W)
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=(2, 8))
        ttk.Entry(row, textvariable=variable, width=42).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(row, text="Browse", command=command).pack(side=tk.LEFT, padx=(6, 0))

    def number_row(self, parent, label, variable):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=3)
        ttk.Label(row, text=label, width=15).pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=variable, width=16).pack(side=tk.LEFT, fill=tk.X, expand=True)

    def vector_row(self, parent, label, x_var, y_var, z_var, state="normal", label_width=7):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=3)
        ttk.Label(row, text=label, width=label_width).pack(side=tk.LEFT)
        for axis, var in [("x", x_var), ("y", y_var), ("z", z_var)]:
            ttk.Label(row, text=axis).pack(side=tk.LEFT, padx=(4, 2))
            ttk.Entry(row, textvariable=var, width=8, state=state).pack(side=tk.LEFT)

    def mapping_row(self, parent, label, key, vars_map=None, boxes_map=None):
        vars_map = vars_map or self.mapping_vars
        boxes_map = boxes_map or self.mapping_boxes
        row = ttk.Frame(parent)
        row.pack(fill=tk.X, pady=3)
        ttk.Label(row, text=label, width=8).pack(side=tk.LEFT)
        combo = ttk.Combobox(row, textvariable=vars_map[key], state="readonly", width=28)
        combo.pack(side=tk.LEFT, fill=tk.X, expand=True)
        boxes_map[key] = combo

    def zero_beta_control(self, parent, option_key):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X)
        for component, label in ZERO_AT_BETA_COMPONENTS:
            ttk.Checkbutton(
                row,
                text=label,
                variable=self.zero_beta_vars[option_key][component],
            ).pack(side=tk.LEFT, padx=(0, 12))

    def current_zero_at_beta(self, option_key):
        return {
            component: var.get()
            for component, var in self.zero_beta_vars[option_key].items()
        }

    def set_zero_at_beta(self, option_key, options):
        enabled = set(normalized_zero_at_beta_options(options))
        for component, var in self.zero_beta_vars[option_key].items():
            var.set(component in enabled)

    def output_order_control(self, parent, order_key):
        row = ttk.Frame(parent)
        row.pack(fill=tk.X)

        listbox = tk.Listbox(row, height=4, exportselection=False)
        listbox.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.output_order_boxes[order_key] = listbox

        buttons = ttk.Frame(row)
        buttons.pack(side=tk.LEFT, padx=(6, 0), fill=tk.Y)
        ttk.Button(buttons, text="Up", command=lambda: self.move_output_group(order_key, -1)).pack(fill=tk.X)
        ttk.Button(buttons, text="Down", command=lambda: self.move_output_group(order_key, 1)).pack(fill=tk.X, pady=(4, 0))
        ttk.Button(buttons, text="Reset", command=lambda: self.reset_output_group_order(order_key)).pack(fill=tk.X, pady=(4, 0))
        self.refresh_output_order_box(order_key)

    def refresh_output_order_box(self, order_key):
        listbox = self.output_order_boxes.get(order_key)
        if not listbox:
            return
        selected = listbox.curselection()
        selected_index = selected[0] if selected else 0
        listbox.delete(0, tk.END)
        for key in self.output_orders[order_key]:
            listbox.insert(tk.END, OUTPUT_GROUP_LABELS[key])
        if self.output_orders[order_key]:
            selected_index = min(selected_index, len(self.output_orders[order_key]) - 1)
            listbox.selection_set(selected_index)

    def move_output_group(self, order_key, delta):
        listbox = self.output_order_boxes.get(order_key)
        if not listbox:
            return
        selected = listbox.curselection()
        if not selected:
            return
        index = selected[0]
        new_index = index + delta
        order = self.output_orders[order_key]
        if new_index < 0 or new_index >= len(order):
            return
        order[index], order[new_index] = order[new_index], order[index]
        self.refresh_output_order_box(order_key)
        listbox.selection_clear(0, tk.END)
        listbox.selection_set(new_index)

    def reset_output_group_order(self, order_key):
        self.output_orders[order_key] = list(DEFAULT_OUTPUT_GROUP_ORDER)
        self.refresh_output_order_box(order_key)

    def build_preview(self, parent):
        toolbar = ttk.Frame(parent)
        toolbar.pack(fill=tk.X, pady=(0, 8))
        ttk.Label(toolbar, text="Preview", font=("Segoe UI", 10, "bold")).pack(side=tk.LEFT)

        frame = ttk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)
        preview = ttk.Treeview(frame, show="headings")
        vsb = ttk.Scrollbar(frame, orient="vertical", command=preview.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=preview.xview)
        preview.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        preview.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        frame.rowconfigure(0, weight=1)
        frame.columnconfigure(0, weight=1)
        return preview

    def browse_input(self):
        path = filedialog.askopenfilename(
            title="Open dense table",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self.input_var.set(path)
            if not self.output_var.get():
                p = Path(path)
                self.output_var.set(str(p.with_name(p.stem + "_coefficients.csv")))

    def browse_output(self):
        path = filedialog.asksaveasfilename(
            title="Save coefficient table",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self.output_var.set(path)

    def browse_recenter_input(self):
        path = filedialog.askopenfilename(
            title="Open coefficient export",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self.recenter_input_var.set(path)
            if not self.recenter_output_var.get():
                p = Path(path)
                self.recenter_output_var.set(str(p.with_name(p.stem + "_modified.csv")))

    def browse_recenter_output(self):
        path = filedialog.asksaveasfilename(
            title="Save modified export",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if path:
            self.recenter_output_var.set(path)

    def load_table(self):
        path = self.input_var.get().strip()
        if not path:
            messagebox.showerror(APP_NAME, "Select an input table first.")
            return
        try:
            self.headers, self.rows = read_csv_table(path)
        except Exception as exc:
            messagebox.showerror(APP_NAME, f"Could not load table:\n{exc}")
            return

        values = [""] + self.headers
        mapping = default_mapping(self.headers)
        for key, combo in self.mapping_boxes.items():
            combo["values"] = ([ZERO_VALUE] + self.headers) if key == "beta" else values
            self.mapping_vars[key].set(mapping.get(key, ""))

        if not self.output_var.get():
            p = Path(path)
            self.output_var.set(str(p.with_name(p.stem + "_coefficients.csv")))

        self.converted_headers = []
        self.converted_rows = []
        self.update_preview(self.headers, self.rows, self.preview)
        self.status_var.set(f"Loaded {len(self.rows)} rows from {Path(path).name}")

    def settings(self):
        return ConversionSettings(
            rho=as_float(self.parameter_vars["rho"].get(), "density"),
            area=as_float(self.parameter_vars["area"].get(), "surface area"),
            velocity=as_float(self.parameter_vars["velocity"].get(), "velocity"),
            span=as_float(self.parameter_vars["span"].get(), "span"),
            chord=as_float(self.parameter_vars["chord"].get(), "chord"),
            output_accuracy=as_float(self.parameter_vars["output_accuracy"].get(), "output accuracy"),
            old_center=(
                as_float(self.parameter_vars["old_x"].get(), "old moment center x"),
                as_float(self.parameter_vars["old_y"].get(), "old moment center y"),
                as_float(self.parameter_vars["old_z"].get(), "old moment center z"),
            ),
            new_center=(
                as_float(self.parameter_vars["new_x"].get(), "new moment center x"),
                as_float(self.parameter_vars["new_y"].get(), "new moment center y"),
                as_float(self.parameter_vars["new_z"].get(), "new moment center z"),
            ),
        )

    def current_mapping(self):
        return {key: var.get() for key, var in self.mapping_vars.items()}

    def current_recenter_mapping(self):
        return {key: var.get() for key, var in self.recenter_mapping_vars.items()}

    def set_vector_vars(self, vars_map, prefix, values):
        for axis, value in zip("xyz", values):
            vars_map[f"{prefix}_{axis}"].set(f"{value:.10g}")

    def set_recenter_settings(self, settings):
        self.recenter_parameter_vars["rho"].set(f"{settings.rho:.10g}")
        self.recenter_parameter_vars["area"].set(f"{settings.area:.10g}")
        self.recenter_parameter_vars["velocity"].set(f"{settings.velocity:.10g}")
        self.recenter_parameter_vars["span"].set(f"{settings.span:.10g}")
        self.recenter_parameter_vars["chord"].set(f"{settings.chord:.10g}")
        self.recenter_parameter_vars["output_accuracy"].set(f"{settings.output_accuracy:.10g}")
        self.set_vector_vars(self.recenter_parameter_vars, "current", settings.old_center)
        self.set_vector_vars(self.recenter_parameter_vars, "new", settings.new_center)

    def recenter_settings(self):
        return ConversionSettings(
            rho=as_float(self.recenter_parameter_vars["rho"].get(), "density"),
            area=as_float(self.recenter_parameter_vars["area"].get(), "surface area"),
            velocity=as_float(self.recenter_parameter_vars["velocity"].get(), "velocity"),
            span=as_float(self.recenter_parameter_vars["span"].get(), "span"),
            chord=as_float(self.recenter_parameter_vars["chord"].get(), "chord"),
            output_accuracy=as_float(self.recenter_parameter_vars["output_accuracy"].get(), "output accuracy"),
            old_center=(
                as_float(self.recenter_parameter_vars["current_x"].get(), "current moment center x"),
                as_float(self.recenter_parameter_vars["current_y"].get(), "current moment center y"),
                as_float(self.recenter_parameter_vars["current_z"].get(), "current moment center z"),
            ),
            new_center=(
                as_float(self.recenter_parameter_vars["new_x"].get(), "new moment center x"),
                as_float(self.recenter_parameter_vars["new_y"].get(), "new moment center y"),
                as_float(self.recenter_parameter_vars["new_z"].get(), "new moment center z"),
            ),
        )

    def load_recenter_export(self):
        path = self.recenter_input_var.get().strip()
        if not path:
            messagebox.showerror(APP_NAME, "Select an export file first.")
            return
        try:
            self.recenter_metadata, self.recenter_headers, self.recenter_rows = read_csv_table_with_metadata(path)
            parsed_metadata, settings = read_export_settings(self.recenter_metadata)
        except Exception as exc:
            messagebox.showerror(APP_NAME, f"Could not load export:\n{exc}")
            return

        self.set_recenter_settings(settings)

        values = [""] + self.recenter_headers
        defaults = default_mapping(self.recenter_headers)
        metadata_alpha = metadata_value(parsed_metadata, ("alpha_column",), "")
        metadata_beta = metadata_value(parsed_metadata, ("beta_column",), "")
        alpha_column = defaults.get("alpha", "") or metadata_alpha
        beta_column = defaults.get("beta", ZERO_VALUE) or metadata_beta or ZERO_VALUE
        if alpha_column not in self.recenter_headers:
            alpha_column = metadata_alpha if metadata_alpha in self.recenter_headers else ""
        if beta_column != ZERO_VALUE and beta_column not in self.recenter_headers:
            beta_column = metadata_beta if metadata_beta in self.recenter_headers else ZERO_VALUE

        for key, combo in self.recenter_mapping_boxes.items():
            combo["values"] = ([ZERO_VALUE] + self.recenter_headers) if key == "beta" else values
        self.recenter_mapping_vars["alpha"].set(alpha_column)
        self.recenter_mapping_vars["beta"].set(beta_column or ZERO_VALUE)

        metadata_order = metadata_value(parsed_metadata, ("output_group_order",), "")
        self.output_orders["recenter"] = parse_output_order_text(metadata_order)
        self.refresh_output_order_box("recenter")
        metadata_zero = metadata_value(parsed_metadata, ("zero_at_beta0",), "")
        self.set_zero_at_beta("recenter", parse_zero_at_beta_text(metadata_zero))

        if not self.recenter_output_var.get():
            p = Path(path)
            self.recenter_output_var.set(str(p.with_name(p.stem + "_modified.csv")))

        self.recentered_headers = []
        self.recentered_rows = []
        self.update_preview(self.recenter_headers, self.recenter_rows, self.recenter_preview)
        self.status_var.set(f"Loaded {len(self.recenter_rows)} export rows from {Path(path).name}")

    def recenter_and_save(self):
        if not self.recenter_rows:
            self.load_recenter_export()
            if not self.recenter_rows:
                return
        out_path = self.recenter_output_var.get().strip()
        if not out_path:
            messagebox.showerror(APP_NAME, "Select an output path first.")
            return
        try:
            settings = self.recenter_settings()
            mapping = self.current_recenter_mapping()
            headers, recentered = recenter_export_rows(
                self.recenter_rows,
                self.recenter_headers,
                mapping,
                settings,
                output_order=self.output_orders["recenter"],
                zero_at_beta=self.current_zero_at_beta("recenter"),
            )
            write_csv_table(
                out_path,
                headers,
                recentered,
                metadata=settings.metadata_row(
                    mapping,
                    self.output_orders["recenter"],
                    self.current_zero_at_beta("recenter"),
                ),
            )
        except Exception as exc:
            messagebox.showerror(APP_NAME, f"Modify failed:\n{exc}")
            return

        self.recentered_headers = headers
        self.recentered_rows = recentered
        self.update_preview(headers, recentered, self.recenter_preview)
        self.status_var.set(f"Saved {len(recentered)} modified rows to {out_path}")

    def convert_and_save(self):
        if not self.rows:
            self.load_table()
            if not self.rows:
                return
        out_path = self.output_var.get().strip()
        if not out_path:
            messagebox.showerror(APP_NAME, "Select an output path first.")
            return
        try:
            settings = self.settings()
            mapping = self.current_mapping()
            headers, converted = convert_rows(
                self.rows,
                self.headers,
                mapping,
                settings,
                output_order=self.output_orders["convert"],
                zero_at_beta=self.current_zero_at_beta("convert"),
            )
            write_csv_table(
                out_path,
                headers,
                converted,
                metadata=settings.metadata_row(
                    mapping,
                    self.output_orders["convert"],
                    self.current_zero_at_beta("convert"),
                ),
            )
        except Exception as exc:
            messagebox.showerror(APP_NAME, f"Conversion failed:\n{exc}")
            return

        self.converted_headers = headers
        self.converted_rows = converted
        self.update_preview(headers, converted, self.preview)
        self.status_var.set(f"Saved {len(converted)} rows to {out_path}")

    def update_preview(self, headers, rows, preview):
        if preview is None:
            return
        preview.delete(*preview.get_children())
        if not headers:
            preview["columns"] = []
            return

        display_headers = list(headers[:40])
        preview["columns"] = display_headers
        for col in display_headers:
            preview.heading(col, text=col)
            preview.column(col, width=max(90, min(170, len(str(col)) * 11)), anchor=tk.W)

        for row in rows[:160]:
            values = []
            for col in display_headers:
                value = row.get(col, "")
                if isinstance(value, float):
                    value = f"{value:.8g}"
                values.append(value)
            preview.insert("", tk.END, values=values)


def self_test():
    rows = [{"Alpha(deg)": "0", "Beta(deg)": "0", "Fx": "-10", "Fy": "2", "Fz": "-30", "Mx": "1", "My": "2", "Mz": "3"}]
    headers = list(rows[0].keys())
    mapping = default_mapping(headers)
    settings = ConversionSettings(
        rho=1.0,
        area=2.0,
        velocity=10.0,
        span=4.0,
        chord=5.0,
        output_accuracy=0.0001,
        old_center=(0.0, 0.0, 0.0),
        new_center=(0.0, 0.0, 0.0),
    )
    _, converted = convert_rows(rows, headers, mapping, settings)
    row = converted[0]
    assert row["Alpha(deg)"] == "0.0000"
    assert row["Fx"] == "-10.0000"
    assert row["Fz"] == "-30.0000"
    assert row["cfx"] == "-0.1000"
    assert row["cd"] == "0.1000"
    assert row["cl"] == "0.3000"
    assert row["cmx"] == "0.0025"

    zero_options = {"fy": True, "mx": True, "mz": True}
    _, zeroed = convert_rows(rows, headers, mapping, settings, zero_at_beta=zero_options)
    zeroed_row = zeroed[0]
    assert zeroed_row["Fy"] == "0.0000"
    assert zeroed_row["Mx"] == "0.0000"
    assert zeroed_row["Mz"] == "0.0000"
    assert zeroed_row["cfy"] == "0.0000"
    assert zeroed_row["cmx"] == "0.0000"
    assert zeroed_row["cmz"] == "0.0000"
    assert zeroed_row["cr"] == "0.0000"
    assert zeroed_row["cn"] == "0.0000"

    beta_rows = [
        {"Alpha(deg)": "0", "Beta(deg)": "5", "Fx": "-10", "Fy": "2", "Fz": "-30", "Mx": "1", "My": "2", "Mz": "3"}
    ]
    _, beta_nonzero = convert_rows(beta_rows, headers, mapping, settings, zero_at_beta=zero_options)
    assert beta_nonzero[0]["Fy"] == "2.0000"
    assert beta_nonzero[0]["Mx"] == "1.0000"
    assert beta_nonzero[0]["Mz"] == "3.0000"

    shifted = ConversionSettings(
        rho=1.0,
        area=2.0,
        velocity=10.0,
        span=4.0,
        chord=5.0,
        output_accuracy=0.0001,
        old_center=(0.0, 0.0, 0.0),
        new_center=(1.0, 0.0, 0.0),
    )
    _, shifted_rows = convert_rows(rows, headers, mapping, shifted)
    assert shifted_rows[0]["My"] == "-28.0000"
    assert shifted_rows[0]["Mz"] == "1.0000"

    recentered = ConversionSettings(
        rho=1.0,
        area=2.0,
        velocity=10.0,
        span=4.0,
        chord=5.0,
        output_accuracy=0.0001,
        old_center=(1.0, 0.0, 0.0),
        new_center=(0.0, 0.0, 0.0),
    )
    shifted_headers, shifted_rows = convert_rows(rows, headers, mapping, shifted)
    _, recentered_rows = recenter_export_rows(shifted_rows, shifted_headers, mapping, recentered)
    assert recentered_rows[0]["My"] == "2.0000"
    assert recentered_rows[0]["Mz"] == "3.0000"

    renormalized = ConversionSettings(
        rho=1.0,
        area=4.0,
        velocity=10.0,
        span=4.0,
        chord=5.0,
        output_accuracy=0.0001,
        old_center=(0.0, 0.0, 0.0),
        new_center=(0.0, 0.0, 0.0),
    )
    original_headers, original_rows = convert_rows(rows, headers, mapping, settings)
    _, renormalized_rows = recenter_export_rows(original_rows, original_headers, mapping, renormalized)
    assert renormalized_rows[0]["cfx"] == "-0.0500"
    assert renormalized_rows[0]["cd"] == "0.0500"
    assert renormalized_rows[0]["cmx"] == "0.0012"

    _, zeroed_recentered = recenter_export_rows(
        shifted_rows,
        shifted_headers,
        mapping,
        recentered,
        zero_at_beta=zero_options,
    )
    assert zeroed_recentered[0]["Fy"] == "0.0000"
    assert zeroed_recentered[0]["Mx"] == "0.0000"
    assert zeroed_recentered[0]["Mz"] == "0.0000"

    order = ["wind_moment", "body_force", "body_moment", "wind_force"]
    ordered_headers, _ = convert_rows(rows, headers, mapping, settings, output_order=order)
    assert ordered_headers[2:8] == WIND_MOMENT_COLUMNS
    assert ordered_headers[8:14] == BODY_FORCE_COLUMNS

    metadata = settings.metadata_row(mapping, order, {"fy": True, "mz": True})
    parsed, loaded = read_export_settings(metadata)
    assert loaded.old_center == settings.new_center
    assert parse_output_order_text(parsed["outputgrouporder"]) == order
    assert parse_zero_at_beta_text(parsed["zeroatbeta0"]) == ["fy", "mz"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-mainloop", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    root = tk.Tk()
    Nondimit(root)
    if args.no_mainloop:
        root.update_idletasks()
        root.destroy()
    else:
        root.mainloop()


if __name__ == "__main__":
    main()
