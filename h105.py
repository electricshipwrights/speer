"""Closed-form parametric H105 low-Reynolds-number hydrofoil section.

Zero-dependency Python 3 implementation reproducing Tom Speer's H105 section
(``H105Coord.csv``) and morphing smoothly from it via intuitive physical design
parameters:

- ``thickness`` (``t/c``) and ``x_thickness`` (``x_t/c``)
- ``camber`` (``f/c``) and ``x_camber`` (``x_c/c``)
- ``le_radius`` (``r_LE/c``) and ``te_wedge_angle_deg`` (``beta_TE``)
- ``le_camber_angle_deg`` (``theta_LE``) and ``te_camber_angle_deg`` (``theta_TE``)
- ``te_thickness`` (``dz_TE/c``) and ``chord`` (``c``)

All mathematics uses only Python's standard library (``math``).
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable


# ---------------------------------------------------------------------------
# Baseline H105 constants (12 fitted Bernstein/Kulfan weights)
# ---------------------------------------------------------------------------

BT_RAW: tuple[float, ...] = (
    0.095931,
    0.101461,
    0.207744,
    0.141241,
    0.214563,
    0.105691,
)
B_NOSE_RAW: float = 0.199820

BC_RAW: tuple[float, ...] = (
    0.080249,
    -0.003853,
    0.105852,
    0.032940,
    0.154909,
)

COMB5: tuple[int, ...] = (1, 5, 10, 10, 5, 1)
COMB4: tuple[int, ...] = (1, 4, 6, 4, 1)

# Exact extrema of the baseline master curves Yt_raw(eta) and Yc_raw(eta)
ETA_T0: float = 0.379438116534871
T0: float = 0.125780798001186
ETA_C0: float = 0.645935633293794
F0: float = 0.016827654576420

# Derived H105 baseline angles and leading-edge radius
THETA_LE0_DEG: float = math.degrees(math.atan(BC_RAW[0]))
THETA_TE0_DEG: float = -math.degrees(math.atan(BC_RAW[4]))
BETA_TE0_DEG: float = 2.0 * math.degrees(math.atan(BT_RAW[5]))
R_LE0: float = 0.5 * (BT_RAW[0] ** 2) / math.hypot(1.0, BC_RAW[0])

H105_DEFAULTS: dict[str, float] = {
    "thickness": T0,
    "x_thickness": ETA_T0,
    "camber": F0,
    "x_camber": ETA_C0,
    "le_radius": R_LE0,
    "te_wedge_angle_deg": BETA_TE0_DEG,
    "le_camber_angle_deg": THETA_LE0_DEG,
    "te_camber_angle_deg": THETA_TE0_DEG,
    "te_thickness": 0.0,
}

_GOLDEN_RATIO: float = 0.5 * (math.sqrt(5.0) - 1.0)


# ---------------------------------------------------------------------------
# Baseline master functions and closed-form morphing helpers
# ---------------------------------------------------------------------------


def _yt_raw(eta: float) -> float:
    """Baseline half-thickness master curve Yt_raw(eta) on eta in [0, 1]."""
    if eta <= 0.0 or eta >= 1.0:
        return 0.0
    om = 1.0 - eta
    s = sum(
        BT_RAW[k] * COMB5[k] * (eta**k) * (om ** (5 - k)) for k in range(6)
    )
    return math.sqrt(eta) * om * s + B_NOSE_RAW * eta * (om**6)


def _dyt_raw_deta(eta: float) -> float:
    """Exact derivative dYt_raw/deta for eta in (0, 1]."""
    if eta <= 0.0:
        return math.inf
    if eta >= 1.0:
        return -BT_RAW[5]
    om = 1.0 - eta
    ds = sum(
        BT_RAW[k] * COMB5[k] * (eta**k) * (om ** (5 - k)) * (k + 0.5 - 6.5 * eta)
        for k in range(6)
    ) / math.sqrt(eta)
    dnose = B_NOSE_RAW * (om**5) * (1.0 - 7.0 * eta)
    return ds + dnose


def _yc_raw(eta: float) -> float:
    """Baseline mean camber master curve Yc_raw(eta) on eta in [0, 1]."""
    if eta <= 0.0 or eta >= 1.0:
        return 0.0
    om = 1.0 - eta
    s = sum(
        BC_RAW[k] * COMB4[k] * (eta**k) * (om ** (4 - k)) for k in range(5)
    )
    return eta * om * s


def _dyc_raw_deta(eta: float) -> float:
    """Exact derivative dYc_raw/deta on eta in [0, 1]."""
    if eta <= 0.0:
        return BC_RAW[0]
    if eta >= 1.0:
        return -BC_RAW[4]
    om = 1.0 - eta
    return sum(
        BC_RAW[k] * COMB4[k] * (eta**k) * (om ** (4 - k)) * (k + 1.0 - 6.0 * eta)
        for k in range(5)
    )


def _mobius_warp(xi: float, eta0: float, x_target: float) -> tuple[float, float]:
    """Return (eta(xi), deta/dxi) for the monotone Mobius peak-shift warp."""
    gamma = (eta0 * (1.0 - x_target)) / (x_target * (1.0 - eta0))
    denom = 1.0 + (gamma - 1.0) * xi
    eta = (gamma * xi) / denom
    deta_dxi = gamma / (denom * denom)
    return eta, deta_dxi


# ---------------------------------------------------------------------------
# Parametric H105 Section class
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class H105Section:
    """Closed-form parametric hydrofoil section morphing cleanly from H105.

    All length parameters (``thickness``, ``x_thickness``, ``camber``,
    ``x_camber``, ``le_radius``, ``te_thickness``) are dimensionless fractions
    of chord. Physical coordinates returned by geometry and export methods are
    scaled by ``chord``.
    """

    thickness: float = T0
    x_thickness: float = ETA_T0
    camber: float = F0
    x_camber: float = ETA_C0
    le_radius: float | None = None
    te_wedge_angle_deg: float | None = None
    le_camber_angle_deg: float | None = None
    te_camber_angle_deg: float | None = None
    te_thickness: float = 0.0
    chord: float = 1.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.chord) or self.chord <= 0.0:
            raise ValueError(f"chord must be positive and finite, got {self.chord}")
        if not math.isfinite(self.thickness) or self.thickness < 0.0:
            raise ValueError(
                f"thickness must be non-negative and finite, got {self.thickness}"
            )
        if not (0.05 <= self.x_thickness <= 0.95):
            raise ValueError(
                f"x_thickness must be in [0.05, 0.95], got {self.x_thickness}"
            )
        if not math.isfinite(self.camber):
            raise ValueError(f"camber must be finite, got {self.camber}")
        if not (0.05 <= self.x_camber <= 0.95):
            raise ValueError(
                f"x_camber must be in [0.05, 0.95], got {self.x_camber}"
            )
        if not math.isfinite(self.te_thickness) or self.te_thickness < 0.0:
            raise ValueError(
                f"te_thickness must be non-negative and finite, got {self.te_thickness}"
            )
        if self.le_radius is not None and (
            not math.isfinite(self.le_radius) or self.le_radius < 0.0
        ):
            raise ValueError(
                f"le_radius must be non-negative and finite, got {self.le_radius}"
            )
        if self.te_wedge_angle_deg is not None and not (
            0.0 <= self.te_wedge_angle_deg < 180.0
        ):
            raise ValueError(
                "te_wedge_angle_deg must be in [0, 180), "
                f"got {self.te_wedge_angle_deg}"
            )
        if self.le_camber_angle_deg is not None and not (
            -85.0 < self.le_camber_angle_deg < 85.0
        ):
            raise ValueError(
                "le_camber_angle_deg must be in (-85, 85), "
                f"got {self.le_camber_angle_deg}"
            )
        if self.te_camber_angle_deg is not None and not (
            -85.0 < self.te_camber_angle_deg < 85.0
        ):
            raise ValueError(
                "te_camber_angle_deg must be in (-85, 85), "
                f"got {self.te_camber_angle_deg}"
            )

        # Verify thickness remains non-negative across the interior
        for i in range(1, 100):
            yt, _ = self.half_thickness(i / 100.0)
            if yt < -1e-12:
                raise ValueError(
                    "Invalid combination of thickness / le_radius / "
                    f"te_wedge_angle_deg causes negative thickness ({yt:.6g}) "
                    f"at xi={i / 100.0:.2f}."
                )

    # ------------------------------------------------------------------
    # Resolved coefficients & effective boundary parameters
    # ------------------------------------------------------------------

    def _camber_coeffs(self) -> tuple[float, float, float, float, float]:
        """Return (gamma_c, s_c, dm_le, dm_te, tan_theta_le)."""
        gamma_c = (ETA_C0 * (1.0 - self.x_camber)) / (
            self.x_camber * (1.0 - ETA_C0)
        )
        s_c = self.camber / F0
        tan_th_le = (
            gamma_c * s_c * BC_RAW[0]
            if self.le_camber_angle_deg is None
            else math.tan(math.radians(self.le_camber_angle_deg))
        )
        tan_th_te = (
            -(s_c * BC_RAW[4]) / gamma_c
            if self.te_camber_angle_deg is None
            else math.tan(math.radians(self.te_camber_angle_deg))
        )
        dm_le = tan_th_le / gamma_c - s_c * BC_RAW[0]
        dm_te = -gamma_c * tan_th_te - s_c * BC_RAW[4]
        return gamma_c, s_c, dm_le, dm_te, tan_th_le

    def _thickness_coeffs(self) -> tuple[float, float, float, float, float]:
        """Return (gamma_t, s_t, da_le, ds_te, b0_xi)."""
        gamma_t = (ETA_T0 * (1.0 - self.x_thickness)) / (
            self.x_thickness * (1.0 - ETA_T0)
        )
        s_t = self.thickness / T0
        _, _, _, _, tan_th_le = self._camber_coeffs()
        cos_th_le = 1.0 / math.hypot(1.0, tan_th_le)

        if self.le_radius is None:
            b0_xi = math.sqrt(gamma_t) * s_t * BT_RAW[0]
            da_le = 0.0
        else:
            b0_xi = math.sqrt(2.0 * self.le_radius / cos_th_le)
            da_le = b0_xi / math.sqrt(gamma_t) - s_t * BT_RAW[0]

        if self.te_wedge_angle_deg is None:
            ds_te = 0.0
        else:
            tan_half_beta = math.tan(0.5 * math.radians(self.te_wedge_angle_deg))
            ds_te = gamma_t * (tan_half_beta + 0.5 * self.te_thickness) - (
                s_t * BT_RAW[5]
            )
        return gamma_t, s_t, da_le, ds_te, b0_xi

    @property
    def effective_le_camber_angle_deg(self) -> float:
        """Resolved leading-edge camber-line angle theta_LE [deg]."""
        _, _, _, _, tan_th_le = self._camber_coeffs()
        return math.degrees(math.atan(tan_th_le))

    @property
    def effective_te_camber_angle_deg(self) -> float:
        """Resolved trailing-edge camber-line angle theta_TE [deg]."""
        _, dyc_te = self.camber_line(1.0)
        return math.degrees(math.atan(dyc_te))

    @property
    def effective_le_radius(self) -> float:
        """Resolved dimensionless leading-edge curvature radius r_LE / c."""
        _, _, _, _, b0_xi = self._thickness_coeffs()
        _, _, _, _, tan_th_le = self._camber_coeffs()
        cos_th_le = 1.0 / math.hypot(1.0, tan_th_le)
        return 0.5 * (b0_xi * b0_xi) * cos_th_le

    @property
    def effective_te_wedge_angle_deg(self) -> float:
        """Resolved trailing-edge boat-tail angle beta_TE = 2*atan(-yt'(1)) [deg]."""
        _, dyt_te = self.half_thickness(1.0)
        return 2.0 * math.degrees(math.atan(-dyt_te))

    # ------------------------------------------------------------------
    # Camber, thickness, and contour point evaluation
    # ------------------------------------------------------------------

    def camber_line(self, xi: float) -> tuple[float, float]:
        """Return dimensionless ``(yc, dyc/dxi)`` at ``xi in [0, 1]``."""
        if xi < 0.0 or xi > 1.0:
            raise ValueError(f"xi must be in [0, 1], got {xi}")
        _, s_c, dm_le, dm_te, _ = self._camber_coeffs()
        eta_c, deta_dxi = _mobius_warp(xi, ETA_C0, self.x_camber)

        yc = s_c * _yc_raw(eta_c)
        dyc_deta = s_c * _dyc_raw_deta(eta_c)

        if eta_c <= ETA_C0:
            q = 1.0 - eta_c / ETA_C0
            yc += dm_le * eta_c * (q**3)
            dyc_deta += dm_le * (q**2) * (1.0 - 4.0 * eta_c / ETA_C0)
        else:
            sig = (eta_c - ETA_C0) / (1.0 - ETA_C0)
            yc += dm_te * (1.0 - eta_c) * (sig**3)
            dyc_deta += dm_te * (sig**2) * (3.0 - 4.0 * sig)

        return yc, dyc_deta * deta_dxi

    def half_thickness(self, xi: float) -> tuple[float, float]:
        """Return dimensionless ``(yt, dyt/dxi)`` at ``xi in [0, 1]``."""
        if xi < 0.0 or xi > 1.0:
            raise ValueError(f"xi must be in [0, 1], got {xi}")
        _, s_t, da_le, ds_te, b0_xi = self._thickness_coeffs()
        eta_t, deta_dxi = _mobius_warp(xi, ETA_T0, self.x_thickness)

        if xi == 0.0:
            if b0_xi > 0.0:
                return 0.0, math.inf
            dyt0 = deta_dxi * (s_t * B_NOSE_RAW) + 0.5 * self.te_thickness
            return 0.0, dyt0

        yt = s_t * _yt_raw(eta_t) + 0.5 * self.te_thickness * xi
        dyt_deta = s_t * _dyt_raw_deta(eta_t)

        if eta_t <= ETA_T0:
            q = 1.0 - eta_t / ETA_T0
            rt = math.sqrt(eta_t)
            yt += da_le * rt * (q**3)
            dyt_deta += da_le * (0.5 / rt) * (q**2) * (1.0 - 7.0 * eta_t / ETA_T0)
        else:
            sig = (eta_t - ETA_T0) / (1.0 - ETA_T0)
            yt += ds_te * (1.0 - eta_t) * (sig**3)
            dyt_deta += ds_te * (sig**2) * (3.0 - 4.0 * sig)

        return yt, dyt_deta * deta_dxi + 0.5 * self.te_thickness

    def _surface_point(self, xi: float, sign: float) -> tuple[float, float]:
        yc, dyc = self.camber_line(xi)
        yt, _ = self.half_thickness(xi)
        cos_th = 1.0 / math.hypot(1.0, dyc)
        sin_th = dyc * cos_th
        x = self.chord * (xi - sign * yt * sin_th)
        y = self.chord * (yc + sign * yt * cos_th)
        return x, y

    def upper_point(self, xi: float) -> tuple[float, float]:
        """Return physical ``(x, y)`` on the upper surface at ``xi in [0, 1]``."""
        return self._surface_point(xi, 1.0)

    def lower_point(self, xi: float) -> tuple[float, float]:
        """Return physical ``(x, y)`` on the lower surface at ``xi in [0, 1]``."""
        return self._surface_point(xi, -1.0)

    def point_from_angle(self, phi: float) -> tuple[float, float]:
        """Return physical ``(x, y)`` at cosine angle ``phi in [-pi, pi]``.

        ``phi >= 0`` traces the upper surface (``phi = pi`` at upper TE,
        ``phi = 0`` at the camber leading edge) and ``phi < 0`` traces the
        lower surface (``phi = -pi`` at lower TE).
        """
        xi = math.sin(0.5 * phi) ** 2
        if xi < 0.0:
            xi = 0.0
        elif xi > 1.0:
            xi = 1.0
        sign = 1.0 if phi >= 0.0 else -1.0
        return self._surface_point(xi, sign)

    def leading_edge_extrema(self) -> tuple[float, float, float]:
        """Return ``(phi_min_x, x_min, y_at_x_min)`` at the geometric nose."""
        a, b = -0.25, 0.25
        for _ in range(70):
            c = b - _GOLDEN_RATIO * (b - a)
            d = a + _GOLDEN_RATIO * (b - a)
            if self.point_from_angle(c)[0] < self.point_from_angle(d)[0]:
                b = d
            else:
                a = c
        phi_min = 0.5 * (a + b)
        x_min, y_min = self.point_from_angle(phi_min)
        return phi_min, x_min, y_min

    def _y_at_x_branch(self, x: float, upper: bool) -> float:
        phi_le, x_le, y_le = self.leading_edge_extrema()
        x_te, y_te = self.upper_point(1.0) if upper else self.lower_point(1.0)
        tol = 1e-12 * self.chord
        if x < x_le - tol or x > x_te + tol:
            raise ValueError(
                f"x={x} is outside branch domain [{x_le:.6g}, {x_te:.6g}]"
            )
        if x <= x_le + tol:
            return y_le
        if x >= x_te - tol:
            return y_te

        if upper:
            lo, hi = phi_le, math.pi
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if self.point_from_angle(mid)[0] < x:
                    lo = mid
                else:
                    hi = mid
        else:
            lo, hi = -math.pi, phi_le
            for _ in range(64):
                mid = 0.5 * (lo + hi)
                if self.point_from_angle(mid)[0] > x:
                    lo = mid
                else:
                    hi = mid
        return self.point_from_angle(0.5 * (lo + hi))[1]

    def y_upper(self, x: float) -> float:
        """Return upper-surface ordinate ``y(x)`` at physical chordwise ``x``."""
        return self._y_at_x_branch(x, upper=True)

    def y_lower(self, x: float) -> float:
        """Return lower-surface ordinate ``y(x)`` at physical chordwise ``x``."""
        return self._y_at_x_branch(x, upper=False)

    # ------------------------------------------------------------------
    # Coordinate generation, structural properties, and file export
    # ------------------------------------------------------------------

    def coordinates(
        self,
        n_points_per_side: int = 33,
        cosine_spacing: bool = True,
        include_z: bool = False,
        z: float = 0.0,
    ) -> list[tuple[float, ...]]:
        """Generate ordered foil contour points (upper TE -> LE -> lower TE).

        ``n_points_per_side`` includes both the leading edge (``xi = 0``) and
        the trailing edge (``xi = 1``), yielding ``2 * n_points_per_side - 1``
        points without duplicating the leading edge (default ``33`` -> 65 points).
        """
        if n_points_per_side < 2:
            raise ValueError(
                f"n_points_per_side must be >= 2, got {n_points_per_side}"
            )
        n = n_points_per_side
        if cosine_spacing:
            xis = [
                0.5 * (1.0 - math.cos(math.pi * i / (n - 1))) for i in range(n)
            ]
        else:
            xis = [i / (n - 1) for i in range(n)]

        pts: list[tuple[float, ...]] = []
        for xi in reversed(xis):
            x, y = self.upper_point(xi)
            pts.append((x, y, z) if include_z else (x, y))
        for xi in xis[1:]:
            x, y = self.lower_point(xi)
            pts.append((x, y, z) if include_z else (x, y))
        return pts

    def section_properties(self, n_panels: int = 400) -> dict[str, float]:
        """Compute closed-contour geometric and bending section properties.

        Uses Green's theorem over a cosine-spaced polygon of ``2 * n_panels``
        points (plus the trailing-edge closure segment when ``te_thickness > 0``).
        Returns ``area``, ``x_centroid``, ``y_centroid``, ``Ixx``, ``Iyy``,
        ``Ixy``, ``Z_upper``, ``Z_lower``, and ``Z_min`` in units of ``chord``.
        """
        if n_panels < 10:
            raise ValueError(f"n_panels must be >= 10, got {n_panels}")
        pts = self.coordinates(n_points_per_side=n_panels + 1, cosine_spacing=True)
        m = len(pts)

        area2 = 0.0
        sum_x = 0.0
        sum_y = 0.0
        ixx0 = 0.0
        iyy0 = 0.0
        ixy0 = 0.0
        y_min = math.inf
        y_max = -math.inf

        for i in range(m):
            x0, y0 = pts[i][0], pts[i][1]
            x1, y1 = pts[(i + 1) % m][0], pts[(i + 1) % m][1]
            if y0 < y_min:
                y_min = y0
            if y0 > y_max:
                y_max = y0
            cross = x0 * y1 - x1 * y0
            area2 += cross
            sum_x += (x0 + x1) * cross
            sum_y += (y0 + y1) * cross
            ixx0 += (y0 * y0 + y0 * y1 + y1 * y1) * cross
            iyy0 += (x0 * x0 + x0 * x1 + x1 * x1) * cross
            ixy0 += (
                x0 * y1 + 2.0 * x0 * y0 + 2.0 * x1 * y1 + x1 * y0
            ) * cross

        area = 0.5 * abs(area2)
        if area == 0.0:
            return {
                "area": 0.0,
                "x_centroid": 0.5 * self.chord,
                "y_centroid": 0.0,
                "Ixx": 0.0,
                "Iyy": 0.0,
                "Ixy": 0.0,
                "Z_upper": 0.0,
                "Z_lower": 0.0,
                "Z_min": 0.0,
            }

        sgn = 1.0 if area2 > 0.0 else -1.0
        xc = sum_x / (3.0 * area2)
        yc = sum_y / (3.0 * area2)
        ixx = sgn * (ixx0 / 12.0) - area * (yc * yc)
        iyy = sgn * (iyy0 / 12.0) - area * (xc * xc)
        ixy = sgn * (ixy0 / 24.0) - area * xc * yc

        dist_up = y_max - yc
        dist_lo = yc - y_min
        z_upper = ixx / dist_up if dist_up > 0.0 else 0.0
        z_lower = ixx / dist_lo if dist_lo > 0.0 else 0.0
        return {
            "area": area,
            "x_centroid": xc,
            "y_centroid": yc,
            "Ixx": ixx,
            "Iyy": iyy,
            "Ixy": ixy,
            "Z_upper": z_upper,
            "Z_lower": z_lower,
            "Z_min": min(z_upper, z_lower),
        }

    def format_csv(
        self,
        n_points_per_side: int = 33,
        cosine_spacing: bool = True,
        include_z: bool = True,
        z: float = 0.0,
        precision: int = 5,
    ) -> str:
        """Format coordinates as ``x, y[, z]`` CSV text."""
        pts = self.coordinates(
            n_points_per_side=n_points_per_side,
            cosine_spacing=cosine_spacing,
            include_z=include_z,
            z=z,
        )
        fmt = f"{{:.{precision}f}}"
        lines = [", ".join(fmt.format(v) for v in row) for row in pts]
        return "\n".join(lines) + "\n"

    def format_selig_dat(
        self,
        name: str = "H105",
        n_points_per_side: int = 65,
        cosine_spacing: bool = True,
        precision: int = 6,
    ) -> str:
        """Format coordinates as a Selig ``.dat`` airfoil file."""
        pts = self.coordinates(
            n_points_per_side=n_points_per_side,
            cosine_spacing=cosine_spacing,
            include_z=False,
        )
        fmt = f"{{: .{precision}f}}"
        lines = [name]
        for x, y in pts:
            lines.append(f"  {fmt.format(x)}  {fmt.format(y)}")
        return "\n".join(lines) + "\n"

    def to_csv(
        self,
        path: str | Path,
        n_points_per_side: int = 33,
        cosine_spacing: bool = True,
        include_z: bool = True,
        z: float = 0.0,
        precision: int = 5,
    ) -> Path:
        """Write coordinates to a CSV file."""
        out = Path(path)
        out.write_text(
            self.format_csv(
                n_points_per_side=n_points_per_side,
                cosine_spacing=cosine_spacing,
                include_z=include_z,
                z=z,
                precision=precision,
            ),
            encoding="utf-8",
        )
        return out

    def to_selig_dat(
        self,
        path: str | Path,
        name: str = "H105",
        n_points_per_side: int = 65,
        cosine_spacing: bool = True,
        precision: int = 6,
    ) -> Path:
        """Write coordinates to a Selig ``.dat`` file."""
        out = Path(path)
        out.write_text(
            self.format_selig_dat(
                name=name,
                n_points_per_side=n_points_per_side,
                cosine_spacing=cosine_spacing,
                precision=precision,
            ),
            encoding="utf-8",
        )
        return out


# ---------------------------------------------------------------------------
# Reference comparison against H105Coord.csv
# ---------------------------------------------------------------------------


def load_csv_coordinates(csv_path: str | Path) -> list[tuple[float, float]]:
    """Read ``(x, y)`` coordinates from a comma-separated coordinate file."""
    pts: list[tuple[float, float]] = []
    with Path(csv_path).open("r", encoding="utf-8") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            fields = [part.strip() for part in stripped.split(",")]
            try:
                x = float(fields[0])
                y = float(fields[1])
            except (IndexError, ValueError):
                continue
            pts.append((x, y))
    if not pts:
        raise ValueError(f"No coordinate rows parsed from {csv_path}")
    return pts


def _min_distance_to_section(
    section: H105Section, xd: float, yd: float, phi_range: tuple[float, float]
) -> float:
    lo_all, hi_all = phi_range
    n_coarse = 120
    best_phi = lo_all
    best_d2 = math.inf
    for j in range(n_coarse + 1):
        phi = lo_all + (hi_all - lo_all) * (j / n_coarse)
        xp, yp = section.point_from_angle(phi)
        d2 = (xp - xd) ** 2 + (yp - yd) ** 2
        if d2 < best_d2:
            best_d2 = d2
            best_phi = phi

    step = (hi_all - lo_all) / n_coarse
    a = max(lo_all, best_phi - 2.0 * step)
    b = min(hi_all, best_phi + 2.0 * step)
    for _ in range(55):
        c = b - _GOLDEN_RATIO * (b - a)
        d = a + _GOLDEN_RATIO * (b - a)
        xc1, yc1 = section.point_from_angle(c)
        xc2, yc2 = section.point_from_angle(d)
        d1 = (xc1 - xd) ** 2 + (yc1 - yd) ** 2
        d2 = (xc2 - xd) ** 2 + (yc2 - yd) ** 2
        if d1 < d2:
            b = d
        else:
            a = c
    xm, ym = section.point_from_angle(0.5 * (a + b))
    return math.hypot(xm - xd, ym - yd)


def compare_with_csv(
    csv_path: str | Path,
    section: H105Section | None = None,
) -> dict[str, object]:
    """Compute exact perpendicular point-to-curve error against a CSV file."""
    sec = H105Section() if section is None else section
    pts = load_csv_coordinates(csv_path)
    min_x_idx = min(range(len(pts)), key=lambda i: pts[i][0])

    errors: list[float] = []
    for idx, (xd, yd) in enumerate(pts):
        if idx < min_x_idx:
            phi_range = (0.0, math.pi)
        elif idx > min_x_idx:
            phi_range = (-math.pi, 0.0)
        else:
            phi_range = (-0.15, 0.15)
        errors.append(_min_distance_to_section(sec, xd, yd, phi_range))

    max_err = max(errors)
    max_idx = errors.index(max_err)
    rms_err = math.sqrt(sum(e * e for e in errors) / len(errors))
    return {
        "n_points": len(pts),
        "rms_distance": rms_err,
        "max_distance": max_err,
        "max_distance_index": max_idx,
        "errors": errors,
    }


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate or verify the closed-form parametric H105 hydrofoil."
    )
    p.add_argument(
        "--thickness",
        type=float,
        default=T0,
        help=f"Maximum thickness-to-chord ratio t/c (default: {T0:.6f})",
    )
    p.add_argument(
        "--x-thickness",
        type=float,
        default=ETA_T0,
        help=f"Chordwise location of max thickness (default: {ETA_T0:.6f})",
    )
    p.add_argument(
        "--camber",
        type=float,
        default=F0,
        help=f"Maximum camber-to-chord ratio f/c (default: {F0:.6f})",
    )
    p.add_argument(
        "--x-camber",
        type=float,
        default=ETA_C0,
        help=f"Chordwise location of max camber (default: {ETA_C0:.6f})",
    )
    p.add_argument(
        "--le-radius",
        type=float,
        default=None,
        help=f"Leading-edge radius r_LE/c (default scales from {R_LE0:.6f})",
    )
    p.add_argument(
        "--te-wedge-deg",
        type=float,
        default=None,
        help=f"Trailing-edge wedge angle [deg] (baseline: {BETA_TE0_DEG:.2f})",
    )
    p.add_argument(
        "--le-camber-deg",
        type=float,
        default=None,
        help=f"Leading-edge camber slope angle [deg] (baseline: {THETA_LE0_DEG:.2f})",
    )
    p.add_argument(
        "--te-camber-deg",
        type=float,
        default=None,
        help=f"Trailing-edge camber slope angle [deg] (baseline: {THETA_TE0_DEG:.2f})",
    )
    p.add_argument(
        "--te-thickness",
        type=float,
        default=0.0,
        help="Finite trailing-edge gap dz_TE/c (default: 0.0)",
    )
    p.add_argument(
        "--chord",
        type=float,
        default=1.0,
        help="Chord length scale factor (default: 1.0)",
    )
    p.add_argument(
        "--points-per-side",
        type=int,
        default=33,
        help="Cosine-spaced points per surface including LE and TE (default: 33)",
    )
    p.add_argument(
        "--format",
        choices=("csv", "dat"),
        default="csv",
        help="Output coordinate format (default: csv)",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path; prints to stdout unless --verify/--summary is set",
    )
    p.add_argument(
        "--summary",
        action="store_true",
        help="Print resolved geometric and structural section properties",
    )
    p.add_argument(
        "--verify",
        nargs="?",
        const=str(Path(__file__).with_name("H105Coord.csv")),
        default=None,
        help="Verify accuracy against reference CSV (defaults to H105Coord.csv)",
    )
    return p


def _format_summary(section: H105Section) -> str:
    props = section.section_properties()
    lines = [
        "Parametric H105 Section Summary:",
        f"  chord                  : {section.chord:.6g}",
        f"  thickness (t/c)        : {section.thickness:.6f} at xi={section.x_thickness:.6f}",
        f"  camber (f/c)           : {section.camber:.6f} at xi={section.x_camber:.6f}",
        f"  effective r_LE / c     : {section.effective_le_radius:.6f}",
        f"  effective beta_TE [deg]: {section.effective_te_wedge_angle_deg:.4f}",
        f"  effective theta_LE[deg]: {section.effective_le_camber_angle_deg:+.4f}",
        f"  effective theta_TE[deg]: {section.effective_te_camber_angle_deg:+.4f}",
        f"  TE thickness (dz_TE/c) : {section.te_thickness:.6f}",
        "Section Structural Properties (in chord units):",
        f"  area                   : {props['area']:.6g}",
        f"  centroid (x, y)        : ({props['x_centroid']:.6g}, {props['y_centroid']:.6g})",
        f"  Ixx, Iyy, Ixy          : ({props['Ixx']:.6g}, {props['Iyy']:.6g}, {props['Ixy']:.6g})",
        f"  Z_upper, Z_lower, Z_min: ({props['Z_upper']:.6g}, {props['Z_lower']:.6g}, {props['Z_min']:.6g})",
    ]
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    section = H105Section(
        thickness=args.thickness,
        x_thickness=args.x_thickness,
        camber=args.camber,
        x_camber=args.x_camber,
        le_radius=args.le_radius,
        te_wedge_angle_deg=args.te_wedge_deg,
        le_camber_angle_deg=args.le_camber_deg,
        te_camber_angle_deg=args.te_camber_deg,
        te_thickness=args.te_thickness,
        chord=args.chord,
    )

    if args.summary:
        print(_format_summary(section))

    if args.verify is not None:
        report = compare_with_csv(args.verify, section=section)
        print(
            f"Verification vs {args.verify} ({report['n_points']} points):\n"
            f"  RMS distance: {report['rms_distance']:.6f}\n"
            f"  Max distance: {report['max_distance']:.6f} "
            f"(index {report['max_distance_index']})"
        )

    if args.output is not None or (not args.summary and args.verify is None):
        if args.format == "csv":
            text = section.format_csv(n_points_per_side=args.points_per_side)
        else:
            text = section.format_selig_dat(n_points_per_side=args.points_per_side)
        if args.output is not None:
            args.output.write_text(text, encoding="utf-8")
        else:
            print(text, end="")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
