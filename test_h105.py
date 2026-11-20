"""Zero-dependency unittest suite for the parametric H105 hydrofoil module."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import h105


ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "H105Coord.csv"


def _sample_xis(n: int = 2001) -> list[float]:
    return [i / (n - 1) for i in range(n)]


def _count_interior_peaks(vals: list[float]) -> int:
    return sum(
        1
        for i in range(1, len(vals) - 1)
        if vals[i] > vals[i - 1] and vals[i] > vals[i + 1]
    )


class TestH105Section(unittest.TestCase):
    def test_baseline_matches_h105coord_csv(self) -> None:
        report = h105.compare_with_csv(CSV_PATH)
        self.assertEqual(report["n_points"], 65)
        self.assertLess(report["rms_distance"], 0.00010)
        self.assertLess(report["max_distance"], 0.00025)

    def test_thickness_and_x_thickness_exact_and_unimodal(self) -> None:
        xis = _sample_xis(4001)
        for t in (0.06, h105.T0, 0.18):
            for xt in (0.20, 0.30, h105.ETA_T0, 0.48, 0.58):
                sec = h105.H105Section(thickness=t, x_thickness=xt)
                thicks = [2.0 * sec.half_thickness(xi)[0] for xi in xis]
                peak_idx = max(range(len(thicks)), key=lambda i: thicks[i])
                self.assertAlmostEqual(thicks[peak_idx], t, places=7)
                self.assertAlmostEqual(xis[peak_idx], xt, delta=3e-4)
                self.assertEqual(_count_interior_peaks(thicks), 1)
                self.assertTrue(all(val > 0.0 for val in thicks[1:-1]))

                yt_at_xt, dyt_at_xt = sec.half_thickness(xt)
                self.assertAlmostEqual(2.0 * yt_at_xt, t, places=12)
                self.assertAlmostEqual(dyt_at_xt, 0.0, places=12)

    def test_camber_and_x_camber_exact_symmetry_and_inversion(self) -> None:
        xis = _sample_xis(4001)
        for xc in (0.25, 0.45, h105.ETA_C0, 0.75, 0.85):
            sec = h105.H105Section(x_camber=xc)
            ycs = [sec.camber_line(xi)[0] for xi in xis]
            peak_idx = max(range(len(ycs)), key=lambda i: ycs[i])
            self.assertAlmostEqual(ycs[peak_idx], h105.F0, places=7)
            self.assertAlmostEqual(xis[peak_idx], xc, delta=3e-4)
            self.assertEqual(_count_interior_peaks(ycs), 1)

            yc_at_xc, dyc_at_xc = sec.camber_line(xc)
            self.assertAlmostEqual(yc_at_xc, h105.F0, places=12)
            self.assertAlmostEqual(dyc_at_xc, 0.0, places=12)

        # Zero camber -> symmetric section
        sym = h105.H105Section(camber=0.0)
        for x in (0.01, 0.1, 0.35, 0.7, 0.95):
            self.assertAlmostEqual(sym.y_upper(x), -sym.y_lower(x), places=13)
        props = sym.section_properties()
        self.assertAlmostEqual(props["y_centroid"], 0.0, places=13)
        self.assertAlmostEqual(props["Z_upper"], props["Z_lower"], places=12)

        # Negative camber -> exact mirror reflection across x-axis
        pos = h105.H105Section()
        neg = h105.H105Section(camber=-h105.F0)
        for x in (0.02, 0.2, 0.5, 0.85, 0.98):
            self.assertAlmostEqual(neg.y_upper(x), -pos.y_lower(x), places=12)
            self.assertAlmostEqual(neg.y_lower(x), -pos.y_upper(x), places=12)

    def test_leading_edge_radius_and_g2_continuity(self) -> None:
        h = 1e-6
        for chord in (1.0, 0.35):
            for r_le in (0.0015, h105.R_LE0, 0.010):
                sec = h105.H105Section(le_radius=r_le, chord=chord)
                self.assertAlmostEqual(sec.effective_le_radius, r_le, places=12)
                target_r = r_le * chord
                for phi0 in (1e-4, -1e-4):
                    xm, ym = sec.point_from_angle(phi0 - h)
                    x0, y0 = sec.point_from_angle(phi0)
                    xp, yp = sec.point_from_angle(phi0 + h)
                    dx = (xp - xm) / (2.0 * h)
                    dy = (yp - ym) / (2.0 * h)
                    ddx = (xp - 2.0 * x0 + xm) / (h * h)
                    ddy = (yp - 2.0 * y0 + ym) / (h * h)
                    kappa = abs(dx * ddy - dy * ddx) / ((dx * dx + dy * dy) ** 1.5)
                    r_measured = 1.0 / kappa
                    rel_err = abs(r_measured - target_r) / target_r
                    self.assertLess(rel_err, 0.002)

    def test_trailing_edge_wedge_and_camber_angles(self) -> None:
        # Independent TE wedge angle preserves peak thickness & location
        for beta in (5.0, h105.BETA_TE0_DEG, 22.0):
            sec = h105.H105Section(te_wedge_angle_deg=beta)
            _, dyt1 = sec.half_thickness(1.0)
            measured_beta = 2.0 * math.degrees(math.atan(-dyt1))
            self.assertAlmostEqual(measured_beta, beta, places=9)
            self.assertAlmostEqual(sec.effective_te_wedge_angle_deg, beta, places=9)
            yt_peak, dyt_peak = sec.half_thickness(sec.x_thickness)
            self.assertAlmostEqual(2.0 * yt_peak, sec.thickness, places=12)
            self.assertAlmostEqual(dyt_peak, 0.0, places=12)

        # Independent LE and TE camber slope angles preserve peak camber & location
        for th_le in (0.0, h105.THETA_LE0_DEG, 9.0):
            sec = h105.H105Section(le_camber_angle_deg=th_le)
            _, dyc0 = sec.camber_line(0.0)
            self.assertAlmostEqual(math.degrees(math.atan(dyc0)), th_le, places=9)
            self.assertAlmostEqual(sec.effective_le_camber_angle_deg, th_le, places=9)
            yc_peak, dyc_peak = sec.camber_line(sec.x_camber)
            self.assertAlmostEqual(yc_peak, sec.camber, places=12)
            self.assertAlmostEqual(dyc_peak, 0.0, places=12)

        for th_te in (-13.0, h105.THETA_TE0_DEG, 0.0, 4.0):
            sec = h105.H105Section(te_camber_angle_deg=th_te)
            _, dyc1 = sec.camber_line(1.0)
            self.assertAlmostEqual(math.degrees(math.atan(dyc1)), th_te, places=9)
            self.assertAlmostEqual(sec.effective_te_camber_angle_deg, th_te, places=9)
            yc_peak, dyc_peak = sec.camber_line(sec.x_camber)
            self.assertAlmostEqual(yc_peak, sec.camber, places=12)
            self.assertAlmostEqual(dyc_peak, 0.0, places=12)

        # Baseline aft lower cusp: y_lower > 0 near TE, suppressed when theta_TE = 0
        base = h105.H105Section()
        no_cusp = h105.H105Section(te_camber_angle_deg=0.0)
        self.assertGreater(base.y_lower(0.975), 0.00055)
        self.assertLess(no_cusp.y_lower(0.975), 0.0)

    def test_te_thickness_and_chord_scaling(self) -> None:
        sec_unit = h105.H105Section(te_thickness=0.004, chord=1.0)
        sec_scaled = h105.H105Section(te_thickness=0.004, chord=2.5)

        xu, yu = sec_scaled.upper_point(1.0)
        xl, yl = sec_scaled.lower_point(1.0)
        self.assertAlmostEqual(math.hypot(xu - xl, yu - yl), 0.004 * 2.5, places=12)

        p1 = sec_unit.section_properties()
        p2 = sec_scaled.section_properties()
        self.assertAlmostEqual(p2["area"] / p1["area"], 2.5**2, places=10)
        self.assertAlmostEqual(p2["Z_min"] / p1["Z_min"], 2.5**3, places=10)
        self.assertAlmostEqual(p2["Ixx"] / p1["Ixx"], 2.5**4, places=10)

    def test_analytical_derivatives_smooth_across_peaks(self) -> None:
        sec = h105.H105Section(
            thickness=0.14,
            x_thickness=0.32,
            camber=0.025,
            x_camber=0.55,
            le_radius=0.006,
            te_wedge_angle_deg=15.0,
            le_camber_angle_deg=6.0,
            te_camber_angle_deg=-10.0,
            te_thickness=0.002,
        )
        h = 1e-6
        test_xis = [0.02, 0.15, 0.32, 0.45, 0.55, 0.75, 0.95]
        for xi in test_xis:
            _, dyc = sec.camber_line(xi)
            fd_yc = (sec.camber_line(xi + h)[0] - sec.camber_line(xi - h)[0]) / (
                2.0 * h
            )
            self.assertAlmostEqual(dyc, fd_yc, places=7)

            _, dyt = sec.half_thickness(xi)
            fd_yt = (
                sec.half_thickness(xi + h)[0] - sec.half_thickness(xi - h)[0]
            ) / (2.0 * h)
            self.assertAlmostEqual(dyt, fd_yt, places=7)

    def test_validation_guards(self) -> None:
        with self.assertRaises(ValueError):
            h105.H105Section(chord=0.0)
        with self.assertRaises(ValueError):
            h105.H105Section(thickness=-0.01)
        with self.assertRaises(ValueError):
            h105.H105Section(x_thickness=0.01)
        with self.assertRaises(ValueError):
            h105.H105Section(x_camber=0.99)
        with self.assertRaises(ValueError):
            h105.H105Section(le_radius=-0.001)
        with self.assertRaises(ValueError):
            h105.H105Section(te_wedge_angle_deg=180.0)
        with self.assertRaises(ValueError):
            h105.H105Section(le_camber_angle_deg=85.0)
        with self.assertRaises(ValueError):
            h105.H105Section(te_thickness=-0.001)

        sec = h105.H105Section()
        with self.assertRaises(ValueError):
            sec.y_upper(-0.01)
        with self.assertRaises(ValueError):
            sec.y_lower(1.01)

    def test_exports_and_cli(self) -> None:
        sec = h105.H105Section()
        coords = sec.coordinates(n_points_per_side=33, include_z=True)
        self.assertEqual(len(coords), 65)
        self.assertEqual(coords[0], (1.0, 0.0, 0.0))
        self.assertEqual(coords[-1], (1.0, 0.0, 0.0))

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            csv_out = sec.to_csv(tmp_path / "out.csv")
            dat_out = sec.to_selig_dat(tmp_path / "out.dat", n_points_per_side=41)
            self.assertEqual(len(csv_out.read_text(encoding="utf-8").splitlines()), 65)
            dat_lines = dat_out.read_text(encoding="utf-8").splitlines()
            self.assertEqual(dat_lines[0], "H105")
            self.assertEqual(len(dat_lines), 82)

        proc = subprocess.run(
            [
                sys.executable,
                str(ROOT / "h105.py"),
                "--verify",
                str(CSV_PATH),
                "--summary",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        self.assertIn("RMS distance:", proc.stdout)
        self.assertIn("Max distance:", proc.stdout)
        self.assertIn("Section Structural Properties", proc.stdout)

    def test_interactive_html_offline_and_js_parity(self) -> None:
        html_path = ROOT / "index.html"
        self.assertTrue(html_path.exists())
        html = html_path.read_text(encoding="utf-8")

        # Zero external network asset requests (<script src>, <link href>, @import, CSS url())
        asset_refs = re.findall(
            r"""(?:<(?:script|link|img|iframe)\b[^>]*\b(?:src|href)\s*=|@import\s+(?:url\()?|url\()\s*["']?\s*(https?://[^"'\s)>]+)""",
            html,
            flags=re.IGNORECASE,
        )
        self.assertEqual(asset_refs, [], f"Found external asset refs: {asset_refs}")
        # Any http:// in SVG XML namespace attributes is strictly the W3C SVG URI
        http_urls = re.findall(r"http://[^\s\"'<>]+", html, flags=re.IGNORECASE)
        self.assertTrue(all(u == "http://www.w3.org/2000/svg" for u in http_urls), http_urls)

        m = re.search(r"<script>\s*(.*?)\s*</script>", html, flags=re.DOTALL)
        self.assertIsNotNone(m)
        js_code = m.group(1)

        node_bin = shutil.which("node") or "/usr/bin/node"
        self.assertTrue(Path(node_bin).exists(), "Node.js required to verify index.html math engine")

        morph_opts = {
            "thickness": 0.14,
            "x_thickness": 0.31,
            "camber": 0.024,
            "x_camber": 0.52,
            "le_radius": 0.006,
            "te_wedge_angle_deg": 16.0,
            "le_camber_angle_deg": 6.0,
            "te_camber_angle_deg": -11.0,
            "te_thickness": 0.0025,
            "chord": 1.8,
        }

        probe_js = (
            js_code
            + "\n"
            + f"""
const {{ H105SectionJS, compareWithCsv, PRESETS }} = globalThis.H105Demo;
const base = new H105SectionJS();
const morph = new H105SectionJS({json.dumps(morph_opts)});
const xs = [0.02, 0.25, 0.55, 0.88, 0.97];
const phis = [];
for (let i = -30; i <= 30; i++) phis.push(i * 0.1);
const curvCheck = phis.map((p) => {{
  const c = morph.curvatureAtAngle(p);
  return [c.x, c.y, c.kappa, c.nx, c.ny];
}});
const presetCoords = PRESETS.map((pr) => new H105SectionJS(pr.params).coordinates(33));
console.log(JSON.stringify({{
  base_coords: base.coordinates(33),
  base_props: base.sectionProperties(400),
  base_csv: compareWithCsv(base),
  morph_upper: [0.0, 0.2, 0.5, 0.9, 1.0].map((xi) => morph.upperPoint(xi)),
  morph_lower: [0.0, 0.2, 0.5, 0.9, 1.0].map((xi) => morph.lowerPoint(xi)),
  morph_y_upper: xs.map((x) => morph.yUpper(x)),
  morph_y_lower: xs.map((x) => morph.yLower(x)),
  morph_props: morph.sectionProperties(400),
  curv_check: curvCheck,
  preset_coords: presetCoords
}}));
"""
        )

        with tempfile.NamedTemporaryFile("w", suffix=".js", encoding="utf-8") as tf:
            tf.write(probe_js)
            tf.flush()
            res = subprocess.run(
                [node_bin, tf.name], capture_output=True, text=True, check=True
            )

        out = json.loads(res.stdout)

        py_base = h105.H105Section()
        py_coords = py_base.coordinates(33)
        for (px, py), (jx, jy) in zip(py_coords, out["base_coords"]):
            self.assertAlmostEqual(px, jx, places=13)
            self.assertAlmostEqual(py, jy, places=13)

        py_props = py_base.section_properties(400)
        for k in ("area", "x_centroid", "y_centroid", "Ixx", "Iyy", "Ixy", "Z_min"):
            self.assertAlmostEqual(py_props[k], out["base_props"][k], places=12)

        py_rep = h105.compare_with_csv(CSV_PATH, py_base)
        self.assertAlmostEqual(py_rep["rms_distance"], out["base_csv"]["rms_distance"], places=10)
        self.assertAlmostEqual(py_rep["max_distance"], out["base_csv"]["max_distance"], places=10)

        py_morph = h105.H105Section(**morph_opts)
        for xi, j_pt in zip((0.0, 0.2, 0.5, 0.9, 1.0), out["morph_upper"]):
            px, py = py_morph.upper_point(xi)
            self.assertAlmostEqual(px, j_pt[0], places=12)
            self.assertAlmostEqual(py, j_pt[1], places=12)
        for xi, j_pt in zip((0.0, 0.2, 0.5, 0.9, 1.0), out["morph_lower"]):
            px, py = py_morph.lower_point(xi)
            self.assertAlmostEqual(px, j_pt[0], places=12)
            self.assertAlmostEqual(py, j_pt[1], places=12)

        for x, j_yu, j_yl in zip((0.02, 0.25, 0.55, 0.88, 0.97), out["morph_y_upper"], out["morph_y_lower"]):
            self.assertAlmostEqual(py_morph.y_upper(x), j_yu, places=11)
            self.assertAlmostEqual(py_morph.y_lower(x), j_yl, places=11)

        py_morph_props = py_morph.section_properties(400)
        for k in ("area", "x_centroid", "y_centroid", "Ixx", "Iyy", "Ixy", "Z_min"):
            self.assertAlmostEqual(py_morph_props[k], out["morph_props"][k], places=11)

        for row in out["curv_check"]:
            self.assertTrue(all(math.isfinite(v) for v in row))
        for coords in out["preset_coords"]:
            self.assertEqual(len(coords), 65)
            self.assertTrue(all(math.isfinite(x) and math.isfinite(y) for x, y in coords))


if __name__ == "__main__":
    unittest.main()
