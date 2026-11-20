# Parametric H105 Low-Reynolds-Number Hydrofoil Section

A zero-dependency pure-Python and closed-form mathematical parameterization of Tom Speer's **H105** low-Reynolds-number laminar hydrofoil section ([H105Coord.csv](H105Coord.csv)), designed with intuitive physical shape controls that morph smoothly and analytically from the H105 baseline.

## Interactive Browser Demo

Try the live workbench on GitHub Pages at **https://electricshipwrights.github.io/speer/**, open [index.html](index.html) directly in any browser (zero external dependencies, zero network requests, zero build step), or serve locally:

```bash
python3 -m http.server 8000
# then open http://localhost:8000/index.html
```

It includes a real-time parametric foil workbench (with 1:1 leading/trailing-edge microscope zooms, osculating nose circle, curvature comb, CSV residual comparison, and `.csv` / `.dat` / `.svg` / Python export), a 5-chapter visual mathematical walkthrough, and `@media print` styling for browser **Print -> Save as PDF**.

## Overview

The baseline H105 coordinates (`65` points from upper trailing edge `(1, 0)` around the nose to lower trailing edge `(1, 0)`) describe an aft-loaded laminar section characterized by:

- **Maximum thickness:** $t/c \approx 12.58\%$ at $x_t/c \approx 0.3794$
- **Maximum camber:** $f/c \approx 1.68\%$ at $x_c/c \approx 0.6459$
- **Leading-edge curvature radius:** $r_{LE}/c \approx 0.00459$ with camber nose angle $\theta_{LE} \approx +4.59^\circ$
- **Trailing-edge geometry:** wedge angle $\beta_{TE} \approx 12.07^\circ$ and aft-loaded camber slope $\theta_{TE} \approx -8.81^\circ$, producing H105's characteristic aft lower-surface cusp ($y_l > 0$ near $x/c \in [0.95, 1)$).

The parametric formulation combines four exact building blocks:

1. **Orthogonal thickness-camber construction** along the camber-line normal ($G^2$-continuous nose curvature without trailing-edge skew).
2. **12-constant Bernstein/Kulfan master curves** ($n_t = 5$ plus an integer-power nose slope term for half-thickness; $n_c = 4$ for mean camber) fitted to [H105Coord.csv](H105Coord.csv).
3. **Monotone Möbius domain warps** $\eta(\xi)$ mapping $[0, 1] \to [0, 1]$ to relocate the thickness and camber peaks ($x_t, x_c$) cleanly without polynomial overshoot.
4. **One-sided $C^2$ Hermite endpoint modifiers** supported strictly on $[0, \eta_0]$ (leading edge) and $[\eta_0, 1]$ (trailing edge), decoupling nose radius ($r_{LE}/c$), TE wedge angle ($\beta_{TE}$), and LE/TE camber angles ($\theta_{LE}, \theta_{TE}$) from peak thickness and peak camber.

---

## Complete Closed-Form Equations

Every quantity below can be implemented directly in Python, spreadsheets, or CAD FeatureScript.

### 1. Baseline H105 Constants

**Thickness master weights ($n_t = 5$, plus nose term):**

$$b = (0.095931,\; 0.101461,\; 0.207744,\; 0.141241,\; 0.214563,\; 0.105691), \qquad b_{\text{nose}} = 0.199820$$

**Camber master weights ($n_c = 4$):**

$$c = (0.080249,\; -0.003853,\; 0.105852,\; 0.032940,\; 0.154909)$$

**Exact baseline extrema of the master curves:**

$$\eta_{t,0} = 0.379438116534871, \qquad T_0 = 0.125780798001186$$
$$\eta_{c,0} = 0.645935633293794, \qquad F_0 = 0.016827654576420$$

**Derived baseline angles and leading-edge radius:**

$$\theta_{LE,0} = \arctan(c_0) = +4.588097^\circ, \qquad \theta_{TE,0} = -\arctan(c_4) = -8.805641^\circ$$
$$\beta_{TE,0} = 2\arctan(b_5) = 12.066500^\circ, \qquad r_{LE,0} = \frac{b_0^2}{2\sqrt{1 + c_0^2}} = 0.00458663$$

---

### 2. Baseline Master Curves on $\eta \in [0, 1]$

**Half-thickness master curve $Y_{t,\text{raw}}(\eta)$ and derivative $Y_{t,\text{raw}}'(\eta)$:**

$$Y_{t,\text{raw}}(\eta) = \sqrt{\eta}(1 - \eta)\sum_{k=0}^{5} b_k \binom{5}{k}\eta^k(1 - \eta)^{5 - k} + b_{\text{nose}}\,\eta(1 - \eta)^6$$

$$Y_{t,\text{raw}}'(\eta) = \frac{1}{\sqrt{\eta}}\sum_{k=0}^{5} b_k \binom{5}{k}\eta^k(1 - \eta)^{5 - k}\left(k + \frac{1}{2} - \frac{13}{2}\eta\right) + b_{\text{nose}}(1 - \eta)^5(1 - 7\eta) \qquad (\eta > 0)$$

**Mean camber master curve $Y_{c,\text{raw}}(\eta)$ and derivative $Y_{c,\text{raw}}'(\eta)$:**

$$Y_{c,\text{raw}}(\eta) = \eta(1 - \eta)\sum_{k=0}^{4} c_k \binom{4}{k}\eta^k(1 - \eta)^{4 - k}$$

$$Y_{c,\text{raw}}'(\eta) = \sum_{k=0}^{4} c_k \binom{4}{k}\eta^k(1 - \eta)^{4 - k}(k + 1 - 6\eta)$$

---

### 3. Chordwise Möbius Domain Warps

Given chord-line parameter $\xi \in [0, 1]$, target max-thickness location $x_t \in [0.05, 0.95]$, and target max-camber location $x_c \in [0.05, 0.95]$:

$$\gamma_t = \frac{\eta_{t,0}(1 - x_t)}{x_t(1 - \eta_{t,0})}, \qquad \eta_t(\xi) = \frac{\gamma_t\,\xi}{1 + (\gamma_t - 1)\xi}, \qquad \frac{d\eta_t}{d\xi} = \frac{\gamma_t}{\big(1 + (\gamma_t - 1)\xi\big)^2}$$

$$\gamma_c = \frac{\eta_{c,0}(1 - x_c)}{x_c(1 - \eta_{c,0})}, \qquad \eta_c(\xi) = \frac{\gamma_c\,\xi}{1 + (\gamma_c - 1)\xi}, \qquad \frac{d\eta_c}{d\xi} = \frac{\gamma_c}{\big(1 + (\gamma_c - 1)\xi\big)^2}$$

These strictly monotone mappings satisfy $\eta_t(x_t) = \eta_{t,0}$ and $\eta_c(x_c) = \eta_{c,0}$ identically.

---

### 4. Camber Line $y_c(\xi)$ with One-Sided $C^2$ Endpoint Modifiers

Let $s_c = f / F_0$ be the camber scale factor. Effective endpoint slopes are:

$$\tan\theta_{LE} = \begin{cases} \gamma_c s_c c_0, & \theta_{LE} \text{ omitted} \\ \tan(\theta_{LE,\text{rad}}), & \theta_{LE} \text{ specified} \end{cases}, \qquad \tan\theta_{TE} = \begin{cases} -\dfrac{s_c c_4}{\gamma_c}, & \theta_{TE} \text{ omitted} \\ \tan(\theta_{TE,\text{rad}}), & \theta_{TE} \text{ specified} \end{cases}$$

The $\eta_c$-space slope increments are:

$$\Delta m_{LE} = \frac{\tan\theta_{LE}}{\gamma_c} - s_c c_0, \qquad \Delta m_{TE} = -\gamma_c \tan\theta_{TE} - s_c c_4$$

Define the one-sided $C^2$ Hermite basis functions ($\sigma_c = (\eta - \eta_{c,0}) / (1 - \eta_{c,0})$):

$$\Psi_{LE}(\eta) = \begin{cases} \eta\left(1 - \dfrac{\eta}{\eta_{c,0}}\right)^3, & \eta \le \eta_{c,0} \\ 0, & \eta > \eta_{c,0} \end{cases}, \qquad \Psi_{LE}'(\eta) = \begin{cases} \left(1 - \dfrac{\eta}{\eta_{c,0}}\right)^2\left(1 - \dfrac{4\eta}{\eta_{c,0}}\right), & \eta \le \eta_{c,0} \\ 0, & \eta > \eta_{c,0} \end{cases}$$

$$\Psi_{TE}(\eta) = \begin{cases} 0, & \eta < \eta_{c,0} \\ (1 - \eta)\sigma_c^3, & \eta \ge \eta_{c,0} \end{cases}, \qquad \Psi_{TE}'(\eta) = \begin{cases} 0, & \eta < \eta_{c,0} \\ \sigma_c^2(3 - 4\sigma_c), & \eta \ge \eta_{c,0} \end{cases}$$

Then:

$$y_c(\xi) = s_c Y_{c,\text{raw}}(\eta_c) + \Delta m_{LE}\Psi_{LE}(\eta_c) + \Delta m_{TE}\Psi_{TE}(\eta_c)$$

$$\frac{dy_c}{d\xi} = \frac{d\eta_c}{d\xi}\Big[s_c Y_{c,\text{raw}}'(\eta_c) + \Delta m_{LE}\Psi_{LE}'(\eta_c) + \Delta m_{TE}\Psi_{TE}'(\eta_c)\Big]$$

---

### 5. Half-Thickness $y_t(\xi)$ with One-Sided $C^2$ Endpoint Modifiers

Let $s_t = t / T_0$ be the thickness scale factor, $\Delta z_{TE}$ be the dimensionless trailing-edge thickness gap, and $\cos\theta_{LE} = (1 + \tan^2\theta_{LE})^{-1/2}$. Because orthogonal construction preserves exact $G^2$ leading-edge curvature with $r_{LE} = \frac{1}{2} B_0^2 \cos\theta_{LE}$, the leading-edge $\sqrt{\xi}$ coefficient and $\eta_t$-space increments are:

$$B_0 = \begin{cases} \sqrt{\gamma_t}\, s_t b_0, & r_{LE} \text{ omitted} \\ \sqrt{\dfrac{2 r_{LE}}{\cos\theta_{LE}}}, & r_{LE} \text{ specified} \end{cases}, \qquad \Delta a_{LE} = \frac{B_0}{\sqrt{\gamma_t}} - s_t b_0$$

$$\Delta s_{TE} = \begin{cases} 0, & \beta_{TE} \text{ omitted} \\ \gamma_t\left(\tan\dfrac{\beta_{TE}}{2} + \dfrac{\Delta z_{TE}}{2}\right) - s_t b_5, & \beta_{TE} \text{ specified} \end{cases}$$

With one-sided $C^2$ basis functions ($\sigma_t = (\eta - \eta_{t,0}) / (1 - \eta_{t,0})$):

$$\Phi_{LE}(\eta) = \begin{cases} \sqrt{\eta}\left(1 - \dfrac{\eta}{\eta_{t,0}}\right)^3, & \eta \le \eta_{t,0} \\ 0, & \eta > \eta_{t,0} \end{cases}, \qquad \Phi_{LE}'(\eta) = \begin{cases} \dfrac{1}{2\sqrt{\eta}}\left(1 - \dfrac{\eta}{\eta_{t,0}}\right)^2\left(1 - \dfrac{7\eta}{\eta_{t,0}}\right), & 0 < \eta \le \eta_{t,0} \\ 0, & \eta > \eta_{t,0} \end{cases}$$

$$\Phi_{TE}(\eta) = \begin{cases} 0, & \eta < \eta_{t,0} \\ (1 - \eta)\sigma_t^3, & \eta \ge \eta_{t,0} \end{cases}, \qquad \Phi_{TE}'(\eta) = \begin{cases} 0, & \eta < \eta_{t,0} \\ \sigma_t^2(3 - 4\sigma_t), & \eta \ge \eta_{t,0} \end{cases}$$

the dimensionless half-thickness and derivative are:

$$y_t(\xi) = s_t Y_{t,\text{raw}}(\eta_t) + \Delta a_{LE}\Phi_{LE}(\eta_t) + \Delta s_{TE}\Phi_{TE}(\eta_t) + \xi\,\frac{\Delta z_{TE}}{2}$$

$$\frac{dy_t}{d\xi} = \frac{d\eta_t}{d\xi}\Big[s_t Y_{t,\text{raw}}'(\eta_t) + \Delta a_{LE}\Phi_{LE}'(\eta_t) + \Delta s_{TE}\Phi_{TE}'(\eta_t)\Big] + \frac{\Delta z_{TE}}{2} \qquad (\xi > 0)$$

---

### 6. Orthogonal Contour Coordinates (Scaled by Chord $C$)

Let $m(\xi) = \frac{dy_c}{d\xi}(\xi)$, $\cos\theta_c = (1 + m^2)^{-1/2}$, and $\sin\theta_c = m\cos\theta_c$. For $\xi = \sin^2(\phi/2) \in [0, 1]$ ($\phi \in [-\pi, \pi]$):

$$x_u(\xi) = C\big(\xi - y_t(\xi)\sin\theta_c\big), \qquad y_u(\xi) = C\big(y_c(\xi) + y_t(\xi)\cos\theta_c\big)$$
$$x_l(\xi) = C\big(\xi + y_t(\xi)\sin\theta_c\big), \qquad y_l(\xi) = C\big(y_c(\xi) - y_t(\xi)\cos\theta_c\big)$$

---

## Physical Design Parameters

| Symbol | Python Argument (`H105Section`) | H105 Default | Meaning & Default Behavior when `None` |
| --- | --- | --- | --- |
| $t/c$ | `thickness` | `0.125781` (`T0`) | Maximum section thickness (perpendicular to camber line) as fraction of chord |
| $x_t/c$ | `x_thickness` | `0.379438` (`ETA_T0`) | Chordwise location of maximum thickness ($\in [0.05, 0.95]$) |
| $f/c$ | `camber` | `0.016828` (`F0`) | Maximum camber-line ordinate as fraction of chord (`0` = symmetric, `< 0` = inverted) |
| $x_c/c$ | `x_camber` | `0.645936` (`ETA_C0`) | Chordwise location of maximum camber ($\in [0.05, 0.95]$) |
| $r_{LE}/c$ | `le_radius` | `None` (`0.004587`) | Leading-edge radius; scales proportionally to $\gamma_t (t/T_0)^2 \cos\theta_{LE}$ when `None` |
| $\beta_{TE}$ | `te_wedge_angle_deg` | `None` (`12.0665°`) | Trailing-edge boat-tail wedge angle $2\arctan(-y_t'(1))$ [deg]; scales with $t/(\gamma_t T_0)$ when `None` |
| $\theta_{LE}$ | `le_camber_angle_deg` | `None` (`+4.5881°`) | Leading-edge camber slope angle [deg]; scales with $\gamma_c (f/F_0)$ when `None` |
| $\theta_{TE}$ | `te_camber_angle_deg` | `None` (`-8.8056°`) | Trailing-edge camber slope angle [deg]; controls aft lower cusp; scales with $f/(\gamma_c F_0)$ when `None` |
| $\Delta z_{TE}/c$ | `te_thickness` | `0.0` | Finite trailing-edge gap normal to camber line (for machinable/composite TEs) |
| $C$ | `chord` | `1.0` | Chord length scale factor for physical coordinates and section properties |

---

## Validation Summary

### Coordinate Match vs. `H105Coord.csv` (`65` points)

Evaluated using `/usr/bin/python3 h105.py --verify H105Coord.csv --summary`:

- **RMS Euclidean distance:** `0.000086 c` (`8.6e-5 c`)
- **Maximum Euclidean distance:** `0.000223 c` (`2.23e-4 c`, at index `30`, $x/c = 0.00212$)
- **Baseline geometric properties (`chord = 1.0`):**
  - Cross-sectional area: `0.0854251`
  - Centroid: `(0.43738, 0.0122362)`
  - Area moments: $I_{xx} = 8.01669 \times 10^{-5}$, $I_{yy} = 4.41331 \times 10^{-3}$, $I_{xy} = 6.23407 \times 10^{-5}$
  - Elastic bending section modulus: $Z_{\min} = 1.25312 \times 10^{-3}$

### Low-Reynolds-Number Polar Match (NeuralFoil `xlarge`, $\alpha \in [-2^\circ, 8^\circ]$)

Comparing `H105Coord.csv` against `H105Section().coordinates(101)`:

| Reynolds Number | $\max |\Delta C_L|$ | $\max |\Delta C_D / C_D|$ | $\max |\Delta C_m|$ |
| --- | --- | --- | --- |
| `1.0e5` | `0.0072` | `0.78%` | `0.0016` |
| `2.0e5` | `0.0041` | `1.19%` | `0.0007` |
| `5.0e5` | `0.0047` | `1.28%` | `0.0009` |

---

## Usage Examples

### Python API (Zero External Dependencies)

```python
from h105 import H105Section, compare_with_csv

# 1. Exact H105 baseline
base = H105Section()
report = compare_with_csv("H105Coord.csv", base)
print(f"RMS={report['rms_distance']:.6f}, Max={report['max_distance']:.6f}")

# 2. Structural root foil: thicker (15%), slightly further forward peak, machinable 0.3% TE gap
root = H105Section(
    thickness=0.15,
    x_thickness=0.35,
    camber=0.018,
    te_thickness=0.003,
    chord=0.14,
)
props = root.section_properties()
print("Root area [m^2]:", props["area"], "Z_min [m^3]:", props["Z_min"])

# 3. High-speed symmetrical strut with custom nose radius
strut = H105Section(camber=0.0, thickness=0.10, x_thickness=0.40, le_radius=0.004)
strut.to_selig_dat("h105_strut.dat", name="H105_STRUT_10", n_points_per_side=65)

# 4. Suppress aft lower-surface cusp (flat/reflex TE) while keeping peak camber unchanged
flat_te = H105Section(te_camber_angle_deg=0.0)
print("y_lower(0.975):", flat_te.y_lower(0.975))
```

### Command-Line Interface

```bash
# Verify baseline fit and print section structural properties
python3 h105.py --verify H105Coord.csv --summary

# Export a 10% thick symmetric strut in Selig .dat format with 65 points per side
python3 h105.py --thickness 0.10 --camber 0.0 --format dat --points-per-side 65 --output strut.dat

# Export a 14% thick foil with 0.25% finite TE gap and 0.12 m chord as CSV
python3 h105.py --thickness 0.14 --te-thickness 0.0025 --chord 0.12 --output foil.csv
```

### Running the Unit Tests

```bash
python3 -m unittest discover -s . -v
```

---

## Credits & License

- **Original H105 Hydrofoil Offsets ([H105Coord.csv](H105Coord.csv)):** Designed and published by **Tom Speer**.
- **Parametric Formulation, Python Implementation & Interactive Workbench:** Released under the [MIT License](LICENSE) (`Copyright (c) 2026 Electric Shipwrights`).
