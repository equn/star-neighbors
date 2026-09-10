"""Coordinate transforms and the Milky Way potential.

Frame conventions used throughout this project
----------------------------------------------
Heliocentric Galactic Cartesian (right-handed), origin at the Sun:
    x -> Galactic center (l=0, b=0)
    y -> direction of Galactic rotation (l=90)
    z -> North Galactic Pole (b=90)

Inertial Galactic frame, origin at the Galactic center, axes *parallel* to the
above (they do not rotate). In this frame at t=0:
    Sun position P_sun = (-sqrt(R0^2 - z_sun^2), 0, +z_sun)
    Sun velocity V_sun = (U_sun, Vcirc + V_sun_pec, W_sun)

Angular momentum then points toward the South Galactic Pole, which is the
correct sense for the Milky Way (it rotates clockwise seen from the NGP).
"""

from __future__ import annotations

import numpy as np

# --- constants ---------------------------------------------------------------

# 1 AU/yr in km/s; converts proper motion (mas/yr) at parallax (mas) into km/s.
KAPPA = 4.740470446

# 1 pc/(km/s) in Myr.  (1 pc = 3.0856775814913673e13 km, 1 Myr = 3.1557e13 s)
PC_KMS_TO_MYR = 0.9777922216807891

# Gravitational constant in kpc (km/s)^2 / Msun.
G_KPC = 4.300917270036279e-6

# Solar position/motion.
R0_KPC = 8.122          # GRAVITY Collaboration 2018
Z_SUN_PC = 20.8         # Bennett & Bovy 2019
V_CIRC = 229.0          # Eilers et al. 2019, at R0
U_SUN, V_SUN_PEC, W_SUN = 11.1, 12.24, 7.25   # Schoenrich, Binney & Dehnen 2010

# ICRS -> Galactic rotation matrix (Hipparcos, ESA SP-1200 vol. 1 sect. 1.5.3).
A_G = np.array([
    [-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
    [+0.4941094278755837, -0.4448296299600112, +0.7469822444972189],
    [-0.8676661490190047, -0.1980763734312015, +0.4559837761750669],
])


def sun_phase_space() -> tuple[np.ndarray, np.ndarray]:
    """Sun's position (kpc) and velocity (km/s) in the inertial Galactic frame."""
    z_sun_kpc = Z_SUN_PC / 1000.0
    x = -np.sqrt(R0_KPC**2 - z_sun_kpc**2)
    pos = np.array([x, 0.0, z_sun_kpc])
    vel = np.array([U_SUN, V_CIRC + V_SUN_PEC, W_SUN])
    return pos, vel


# --- astrometry -> Cartesian phase space -------------------------------------

def icrs_to_galactic_cartesian(ra_deg, dec_deg, parallax_mas, pmra, pmdec, rv):
    """Convert ICRS astrometry to heliocentric Galactic Cartesian phase space.

    `pmra` must be mu_alpha* (i.e. already multiplied by cos(dec)), in mas/yr.
    `rv` is radial velocity in km/s, positive = receding.

    Returns (pos_pc, vel_kms), each shape (N, 3).
    """
    ra = np.radians(np.asarray(ra_deg, dtype=float))
    dec = np.radians(np.asarray(dec_deg, dtype=float))
    plx = np.asarray(parallax_mas, dtype=float)
    pmra = np.asarray(pmra, dtype=float)
    pmdec = np.asarray(pmdec, dtype=float)
    rv = np.asarray(rv, dtype=float)

    d_pc = 1000.0 / plx

    ca, sa = np.cos(ra), np.sin(ra)
    cd, sd = np.cos(dec), np.sin(dec)

    # ICRS unit vectors: radial, +RA, +Dec.
    r_hat = np.stack([cd * ca, cd * sa, sd], axis=-1)
    p_hat = np.stack([-sa, ca, np.zeros_like(sa)], axis=-1)
    q_hat = np.stack([-sd * ca, -sd * sa, cd], axis=-1)

    v_icrs = (
        rv[:, None] * r_hat
        + (KAPPA * pmra / plx)[:, None] * p_hat
        + (KAPPA * pmdec / plx)[:, None] * q_hat
    )
    pos_icrs = d_pc[:, None] * r_hat

    pos_gal = pos_icrs @ A_G.T
    vel_gal = v_icrs @ A_G.T
    return pos_gal, vel_gal


# --- Milky Way potential -----------------------------------------------------
# gala's MilkyWayPotential (Price-Whelan 2017): Miyamoto-Nagai disk, Hernquist
# bulge and nucleus, NFW halo.
#
# Note the deliberate inconsistency: this potential gives v_c(8.122 kpc) =
# 231.5 km/s, while the Sun's velocity above uses Eilers et al.'s measured
# 229 km/s. Adopting the measurement rather than the model's own value puts the
# Sun on a very slightly non-circular orbit. Stage 9 measures what that costs:
# 0.0014 pc at the 99th percentile, 0.0025 pc worst case among stars closing to
# within 2 pc - far below the disc-shape terms that dominate the model budget.

MN_M, MN_A, MN_B = 6.8e10, 3.0, 0.28          # Msun, kpc, kpc
BULGE_M, BULGE_C = 5.00e9, 1.0
NUC_M, NUC_C = 1.71e9, 0.07
NFW_M, NFW_RS = 5.4e11, 15.62


def acceleration(pos_kpc: np.ndarray) -> np.ndarray:
    """Galactic acceleration in (km/s)^2/kpc at positions (N,3) in kpc."""
    x, y, z = pos_kpc[..., 0], pos_kpc[..., 1], pos_kpc[..., 2]
    R2 = x * x + y * y
    r = np.sqrt(R2 + z * z)
    r = np.maximum(r, 1e-8)

    # Miyamoto-Nagai disk
    zb = np.sqrt(z * z + MN_B * MN_B)
    azb = MN_A + zb
    den = (R2 + azb * azb) ** 1.5
    gm_d = G_KPC * MN_M
    ax = -gm_d * x / den
    ay = -gm_d * y / den
    az = -gm_d * z * azb / (zb * den)

    # Hernquist spheroids (bulge + nucleus)
    for m, c in ((BULGE_M, BULGE_C), (NUC_M, NUC_C)):
        f = -G_KPC * m / (r * (r + c) ** 2)
        ax += f * x
        ay += f * y
        az += f * z

    # NFW halo
    s = r / NFW_RS
    mass_factor = np.log(1.0 + s) - s / (1.0 + s)
    f = -G_KPC * NFW_M * mass_factor / (r**3)
    ax += f * x
    ay += f * y
    az += f * z

    return np.stack([ax, ay, az], axis=-1)


def circular_velocity(R_kpc: float) -> float:
    """Circular speed in the midplane, km/s (for model validation)."""
    pos = np.array([[R_kpc, 0.0, 0.0]])
    a = acceleration(pos)[0]
    return float(np.sqrt(R_kpc * -a[0]))


# --- straight-line (naive) encounter solution --------------------------------

def linear_perihelion(pos_pc: np.ndarray, vel_kms: np.ndarray):
    """Closest approach assuming the Sun and star both move in straight lines.

    This is the model behind the reference charts. Returns (d_min_pc, t_min_Myr)
    where t is measured from today, positive = future.
    """
    v2 = np.sum(vel_kms**2, axis=-1)
    v2 = np.maximum(v2, 1e-12)
    r_dot_v = np.sum(pos_pc * vel_kms, axis=-1)
    t_pc_kms = -r_dot_v / v2                       # pc/(km/s)
    closest = pos_pc + t_pc_kms[:, None] * vel_kms
    return np.linalg.norm(closest, axis=-1), t_pc_kms * PC_KMS_TO_MYR
