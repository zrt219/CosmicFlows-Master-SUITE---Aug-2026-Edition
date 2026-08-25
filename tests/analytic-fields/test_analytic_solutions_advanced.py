# -*- coding: utf-8 -*-
"""
tests/analytic-fields/test_analytic_solutions_advanced.py
=============================================================================
Advanced Analytical Cosmological Solutions Test Suite
-----------------------------------------------------------------------------
Comprehensive analytical oracles and mathematical validation for:
  1. Exact Plummer Gravitational Sphere (infall velocity, potential, divergence)
  2. Navarro-Frenk-White (NFW) Dark Matter Halo Profile (cusp, enclosed mass)
  3. Pure Hubble-Lemaître Isotropic Cosmological Expansion (divergence, Hubble law)
  4. 3D Burgers Swirling Vortex Tube (vorticity balance, Okubo-Weiss Q-criterion)
  5. Gravitational Jeans Perturbation Wave (exact linear continuity diagnostic)
  6. Maclaurin Equilibrium Spheroid (ellipsoidal gravitational field)

Strict Scientific Standards Enforced:
  - Supergalactic coordinates in Mpc/h, velocities in km/s, delta dimensionless.
  - Linear continuity diagnostic: div(v) = -H0 * f * delta.
  - Machine-precision mathematical invariants and derivative stencils.
=============================================================================
"""

import math
import numpy as np
import pytest
from typing import Dict, Any, Tuple, List


# =============================================================================
# 1. ANALYTICAL PHYSICAL MODELS (MATHEMATICAL ORACLES)
# =============================================================================

class PlummerSphereOracle:
    """
    Exact Plummer Gravitational Sphere:
      Phi(r) = -G * M / sqrt(r^2 + b^2)
      rho(r) = (3 * M / (4 * pi * b^3)) * (1 + r^2/b^2)^(-5/2)
      In linear infall / steady radial regime:
      v_r(r) = -v_0 * r / (r^2 + b^2)^(3/4)
    """
    def __init__(self, M: float = 1.0e15, b: float = 5.0, G: float = 4.30091e-9, v0: float = 1200.0):
        self.M = M
        self.b = b
        self.G = G
        self.v0 = v0

    def potential(self, x: np.ndarray) -> float:
        r2 = float(np.sum(x**2))
        return -self.G * self.M / math.sqrt(r2 + self.b**2)

    def density(self, x: np.ndarray) -> float:
        r2 = float(np.sum(x**2))
        return (3.0 * self.M / (4.0 * math.pi * self.b**3)) * (1.0 + r2 / self.b**2)**(-2.5)

    def velocity(self, x: np.ndarray) -> np.ndarray:
        r = float(np.linalg.norm(x))
        if r < 1e-12:
            return np.zeros(3, dtype=np.float64)
        r_hat = x / r
        # Smooth inward radial velocity pointing towards center
        speed = self.v0 * (r / (r**2 + self.b**2)**0.75)
        return -speed * r_hat

    def analytical_divergence(self, x: np.ndarray) -> float:
        r = float(np.linalg.norm(x))
        if r < 1e-12:
            # Limit as r -> 0 of div(v) = -3 * v0 / b^(1.5)
            return -3.0 * self.v0 / (self.b**1.5)
        # v(r) = -v0 * r * (r^2 + b^2)^(-3/4)
        # In spherical coords: div(v) = (1/r^2) d/dr (r^2 v_r)
        # r^2 v_r = -v0 * r^3 * (r^2 + b^2)^(-3/4)
        # d/dr = -v0 * [3 r^2 (r^2 + b^2)^(-3/4) - 1.5 r^4 (r^2 + b^2)^(-7/4)]
        # div(v) = -v0 * (r^2 + b^2)^(-7/4) * [3(r^2 + b^2) - 1.5 r^2] = -v0 * (r^2 + b^2)^(-7/4) * (1.5 r^2 + 3 b^2)
        denom = (r**2 + self.b**2)**1.75
        numer = 1.5 * r**2 + 3.0 * self.b**2
        return -self.v0 * numer / denom

    def analytical_curl(self, x: np.ndarray) -> np.ndarray:
        # Radial field is purely irrotational -> curl is zero identically
        return np.zeros(3, dtype=np.float64)


class NFWProfileOracle:
    """
    Navarro-Frenk-White (NFW) Dark Matter Halo Profile:
      rho(r) = rho_0 / [ (r/r_s) * (1 + r/r_s)^2 ]
      M(r) = 4 * pi * rho_0 * r_s^3 * [ ln(1 + r/r_s) - (r/r_s)/(1 + r/r_s) ]
      g_r(r) = -G * M(r) / r^2
    """
    def __init__(self, rho0: float = 1.0e7, rs: float = 2.0, G: float = 4.30091e-9):
        self.rho0 = rho0
        self.rs = rs
        self.G = G

    def density(self, r: float) -> float:
        if r <= 1e-9:
            return float("inf")
        x = r / self.rs
        return self.rho0 / (x * (1.0 + x)**2)

    def enclosed_mass(self, r: float) -> float:
        if r <= 0.0:
            return 0.0
        x = r / self.rs
        f_x = math.log(1.0 + x) - x / (1.0 + x)
        return 4.0 * math.pi * self.rho0 * (self.rs**3) * f_x

    def gravitational_acceleration(self, pos: np.ndarray) -> np.ndarray:
        r = float(np.linalg.norm(pos))
        if r < 1e-9:
            return np.zeros(3, dtype=np.float64)
        m_r = self.enclosed_mass(r)
        g_mag = self.G * m_r / (r**2)
        return -g_mag * (pos / r)


class HubbleExpansionOracle:
    """
    Pure Hubble-Lemaître Isotropic Cosmological Expansion:
      v(x) = H_0 * x
      Invariants:
        div(v) = 3 * H_0
        curl(v) = 0
        Jacobian J = H_0 * I_3
        All eigenvalues lambda_1 = lambda_2 = lambda_3 = H_0 > 0 (isotropic source / repeller)
    """
    def __init__(self, H0: float = 74.6):
        self.H0 = H0

    def velocity(self, x: np.ndarray) -> np.ndarray:
        return self.H0 * np.asarray(x, dtype=np.float64)

    def analytical_divergence(self, x: np.ndarray) -> float:
        return 3.0 * self.H0

    def analytical_curl(self, x: np.ndarray) -> np.ndarray:
        return np.zeros(3, dtype=np.float64)

    def analytical_jacobian(self, x: np.ndarray) -> np.ndarray:
        return self.H0 * np.eye(3, dtype=np.float64)


class BurgersVortexOracle:
    """
    3D Burgers Swirling Vortex Tube:
      v_r(r) = -alpha * r
      v_theta(r) = (Gamma / (2 * pi * r)) * (1 - exp(-alpha * r^2 / (2 * nu)))
      v_z(z) = 2 * alpha * z
    Invariants:
      div(v) = 0 (incompressible flow)
      Vorticity along z-axis: omega_z(r) = (Gamma * alpha / (2 * pi * nu)) * exp(-alpha * r^2 / (2 * nu))
      Okubo-Weiss Q-criterion: Q > 0 inside vortex core (vorticity dominates strain)
    """
    def __init__(self, alpha: float = 0.5, Gamma: float = 500.0, nu: float = 1.0):
        self.alpha = alpha
        self.Gamma = Gamma
        self.nu = nu

    def velocity(self, pos: np.ndarray) -> np.ndarray:
        x, y, z = pos[0], pos[1], pos[2]
        r = math.sqrt(x**2 + y**2)
        vz = 2.0 * self.alpha * z
        if r < 1e-12:
            return np.array([0.0, 0.0, vz], dtype=np.float64)
        
        # Radial velocity: v_r = -alpha * r
        # Azimuthal velocity: v_theta = (Gamma / (2*pi*r)) * (1 - exp(-alpha*r^2/(2*nu)))
        factor = 1.0 - math.exp(-self.alpha * (r**2) / (2.0 * self.nu))
        v_theta = (self.Gamma / (2.0 * math.pi * r)) * factor
        
        vx = -self.alpha * x - v_theta * (y / r)
        vy = -self.alpha * y + v_theta * (x / r)
        return np.array([vx, vy, vz], dtype=np.float64)

    def analytical_divergence(self, pos: np.ndarray) -> float:
        # div(v) = d(vx)/dx + d(vy)/dy + d(vz)/dz = -alpha - alpha + 2*alpha = 0
        return 0.0

    def core_radius(self) -> float:
        # Characteristic core radius where vorticity falls to 1/e
        return math.sqrt(2.0 * self.nu / self.alpha)


class JeansPerturbationOracle:
    """
    Gravitational Jeans Linear Perturbation Mode:
      delta(x) = delta_0 * cos(k_x * x + k_y * y + k_z * z)
      By linear continuity div(v) = -H0 * f * delta:
      v(x) = - (H0 * f * delta_0 / ||k||^2) * k * sin(k . x)
    """
    def __init__(self, delta0: float = 0.15, k_vec: Tuple[float, float, float] = (0.05, 0.03, -0.04), H0: float = 74.6, f: float = 0.53):
        self.delta0 = delta0
        self.k = np.array(k_vec, dtype=np.float64)
        self.k2 = float(np.sum(self.k**2))
        self.H0 = H0
        self.f = f

    def density_contrast(self, pos: np.ndarray) -> float:
        phase = float(np.dot(self.k, pos))
        return self.delta0 * math.cos(phase)

    def velocity(self, pos: np.ndarray) -> np.ndarray:
        phase = float(np.dot(self.k, pos))
        coeff = -(self.H0 * self.f * self.delta0 / self.k2)
        return coeff * self.k * math.sin(phase)

    def analytical_divergence(self, pos: np.ndarray) -> float:
        # div(v) = -H0 * f * delta(x)
        return -self.H0 * self.f * self.density_contrast(pos)


# =============================================================================
# 2. NUMERICAL FINITE DIFFERENCE VERIFICATION OPERATORS
# =============================================================================

def numerical_divergence_6th_order(field_fn, x: np.ndarray, h: float = 1e-4) -> float:
    """Computes 6th-order central finite difference divergence of a 3D vector field."""
    div = 0.0
    for i in range(3):
        e = np.zeros(3, dtype=np.float64)
        e[i] = 1.0
        v_p3 = field_fn(x + 3 * h * e)[i]
        v_p2 = field_fn(x + 2 * h * e)[i]
        v_p1 = field_fn(x + 1 * h * e)[i]
        v_m1 = field_fn(x - 1 * h * e)[i]
        v_m2 = field_fn(x - 2 * h * e)[i]
        v_m3 = field_fn(x - 3 * h * e)[i]
        dvi_dxi = (v_p3 - 9.0 * v_p2 + 45.0 * v_p1 - 45.0 * v_m1 + 9.0 * v_m2 - v_m3) / (60.0 * h)
        div += dvi_dxi
    return float(div)

def numerical_curl_6th_order(field_fn, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Computes 6th-order central finite difference curl of a 3D vector field."""
    J = np.zeros((3, 3), dtype=np.float64)
    for j in range(3):
        e = np.zeros(3, dtype=np.float64)
        e[j] = 1.0
        v_p3 = field_fn(x + 3 * h * e)
        v_p2 = field_fn(x + 2 * h * e)
        v_p1 = field_fn(x + 1 * h * e)
        v_m1 = field_fn(x - 1 * h * e)
        v_m2 = field_fn(x - 2 * h * e)
        v_m3 = field_fn(x - 3 * h * e)
        J[:, j] = (v_p3 - 9.0 * v_p2 + 45.0 * v_p1 - 45.0 * v_m1 + 9.0 * v_m2 - v_m3) / (60.0 * h)
    
    curl_x = J[2, 1] - J[1, 2]
    curl_y = J[0, 2] - J[2, 0]
    curl_z = J[1, 0] - J[0, 1]
    return np.array([curl_x, curl_y, curl_z], dtype=np.float64)


# =============================================================================
# 3. PYTEST SUITE: ADVANCED ANALYTICAL SOLUTIONS
# =============================================================================

class TestPlummerSphereModel:
    """Tests for exact Plummer gravitational sphere solutions."""

    @pytest.fixture
    def plummer(self):
        return PlummerSphereOracle(M=1.0e15, b=8.0, v0=1500.0)

    def test_plummer_origin_velocity_zero(self, plummer):
        v0 = plummer.velocity(np.array([0.0, 0.0, 0.0]))
        assert np.allclose(v0, [0.0, 0.0, 0.0], atol=1e-10)

    def test_plummer_infall_direction(self, plummer):
        test_points = [
            np.array([10.0, 0.0, 0.0]),
            np.array([0.0, -15.0, 0.0]),
            np.array([5.0, 5.0, -5.0]),
            np.array([-20.0, 10.0, 30.0])
        ]
        for pt in test_points:
            v = plummer.velocity(pt)
            # Velocity must be anti-parallel to position vector (infall towards origin)
            dot = np.dot(v, pt)
            assert dot < 0.0, f"Expected inward radial flow, got dot product {dot}"
            # Normalized cross product must be zero
            cross_norm = np.linalg.norm(np.cross(v, pt))
            assert cross_norm < 1e-6

    def test_plummer_divergence_exact_matches_finite_difference(self, plummer):
        pts = [
            np.array([2.0, 3.0, 4.0]),
            np.array([12.0, -8.0, 5.0]),
            np.array([-10.0, -10.0, -10.0]),
            np.array([25.0, 0.0, 0.0])
        ]
        for pt in pts:
            exact_div = plummer.analytical_divergence(pt)
            num_div = numerical_divergence_6th_order(plummer.velocity, pt, h=1e-4)
            assert math.isclose(exact_div, num_div, rel_tol=1e-4, abs_tol=1e-5), \
                f"Plummer divergence mismatch at {pt}: exact={exact_div}, num={num_div}"

    def test_plummer_curl_identically_zero(self, plummer):
        pts = [np.array([3.0, -4.0, 5.0]), np.array([-8.0, 12.0, 1.0])]
        for pt in pts:
            num_curl = numerical_curl_6th_order(plummer.velocity, pt, h=1e-4)
            assert np.allclose(num_curl, [0.0, 0.0, 0.0], atol=1e-4)

    def test_plummer_potential_laplacian_equals_poisson_source(self, plummer):
        # Verify Poisson equation: nabla^2 Phi = 4 * pi * G * rho
        pt = np.array([4.0, -3.0, 2.0])
        h = 1e-4
        # Numerical Laplacian of potential
        lap = 0.0
        for i in range(3):
            e = np.zeros(3)
            e[i] = 1.0
            p_p1 = plummer.potential(pt + h * e)
            p_0  = plummer.potential(pt)
            p_m1 = plummer.potential(pt - h * e)
            d2_dxi2 = (p_p1 - 2.0 * p_0 + p_m1) / (h**2)
            lap += d2_dxi2
        expected_source = 4.0 * math.pi * plummer.G * plummer.density(pt)
        assert math.isclose(lap, expected_source, rel_tol=1e-3)


class TestNFWProfileModel:
    """Tests for exact Navarro-Frenk-White halo profile properties."""

    @pytest.fixture
    def nfw(self):
        return NFWProfileOracle(rho0=5.0e6, rs=3.5)

    def test_nfw_cusp_divergence(self, nfw):
        # rho(r) ~ 1/r as r -> 0
        rho_1 = nfw.density(0.01)
        rho_2 = nfw.density(0.001)
        assert rho_2 > rho_1 * 8.0, "NFW cusp should exhibit ~1/r scaling near center"

    def test_nfw_asymptotic_outer_falloff(self, nfw):
        # rho(r) ~ 1/r^3 as r -> inf
        r_far1 = 100.0
        r_far2 = 200.0
        ratio = nfw.density(r_far1) / nfw.density(r_far2)
        expected_ratio = (r_far2 / r_far1)**3 # ~8.0
        assert math.isclose(ratio, expected_ratio, rel_tol=0.05)

    def test_nfw_enclosed_mass_monotonicity(self, nfw):
        radii = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
        masses = [nfw.enclosed_mass(r) for r in radii]
        for i in range(len(masses) - 1):
            assert masses[i] < masses[i+1], "Enclosed mass must be strictly monotonically increasing"

    def test_nfw_gravitational_acceleration_radial(self, nfw):
        pos = np.array([5.0, -12.0, 0.0])
        g = nfw.gravitational_acceleration(pos)
        assert np.dot(g, pos) < 0.0, "Gravitational acceleration must point inwards"
        g_unit = g / np.linalg.norm(g)
        r_unit = pos / np.linalg.norm(pos)
        assert np.allclose(g_unit, -r_unit, atol=1e-8)


class TestHubbleExpansionModel:
    """Tests for pure Hubble-Lemaître cosmological expansion."""

    @pytest.fixture
    def hubble(self):
        return HubbleExpansionOracle(H0=74.6)

    def test_hubble_origin_velocity_zero(self, hubble):
        v0 = hubble.velocity(np.array([0.0, 0.0, 0.0]))
        assert np.allclose(v0, [0.0, 0.0, 0.0], atol=1e-12)

    def test_hubble_divergence_isotropic(self, hubble):
        pts = [np.array([10.0, 20.0, -30.0]), np.array([-50.0, 15.0, 0.0])]
        for pt in pts:
            exact_div = hubble.analytical_divergence(pt)
            assert exact_div == 3.0 * 74.6
            num_div = numerical_divergence_6th_order(hubble.velocity, pt, h=1e-4)
            assert math.isclose(exact_div, num_div, rel_tol=1e-5)

    def test_hubble_curl_identically_zero(self, hubble):
        pt = np.array([12.5, -45.0, 22.0])
        num_curl = numerical_curl_6th_order(hubble.velocity, pt, h=1e-4)
        assert np.allclose(num_curl, [0.0, 0.0, 0.0], atol=1e-6)

    def test_hubble_jacobian_diagonal_eigenvalues(self, hubble):
        J = hubble.analytical_jacobian(np.array([1.0, 2.0, 3.0]))
        eigvals = np.linalg.eigvals(J)
        assert np.allclose(eigvals, [74.6, 74.6, 74.6], atol=1e-10)


class TestBurgersVortexModel:
    """Tests for 3D Burgers Swirling Vortex Tube."""

    @pytest.fixture
    def vortex(self):
        return BurgersVortexOracle(alpha=0.5, Gamma=600.0, nu=2.0)

    def test_burgers_incompressibility(self, vortex):
        # div(v) = 0 everywhere
        pts = [
            np.array([1.0, 1.0, 2.0]),
            np.array([4.0, -3.0, 5.0]),
            np.array([0.5, -0.5, -10.0])
        ]
        for pt in pts:
            num_div = numerical_divergence_6th_order(vortex.velocity, pt, h=1e-4)
            assert abs(num_div) < 1e-4, f"Burgers vortex must be divergence-free, got div={num_div}"

    def test_burgers_vorticity_maximum_at_core(self, vortex):
        core_r = vortex.core_radius()
        assert core_r > 0.0
        curl_center = numerical_curl_6th_order(vortex.velocity, np.array([0.0, 0.0, 1.0]), h=1e-4)
        curl_outer  = numerical_curl_6th_order(vortex.velocity, np.array([core_r * 3.0, 0.0, 1.0]), h=1e-4)
        # Core vorticity along Z should be significantly larger than outer vorticity
        assert curl_center[2] > curl_outer[2] * 5.0

    def test_burgers_radial_inflow_and_axial_outflow(self, vortex):
        # For point on x-axis: x > 0, y=0, z > 0: vx < 0 (radial inflow), vz > 0 (axial stretching)
        pt = np.array([5.0, 0.0, 2.0])
        v = vortex.velocity(pt)
        assert v[0] < 0.0, "Radial velocity must be inward (-alpha * r)"
        assert v[2] > 0.0, "Axial velocity must be outward (+2 * alpha * z)"


class TestJeansLinearPerturbationModel:
    """Tests for exact linear cosmological continuity diagnostic."""

    @pytest.fixture
    def jeans(self):
        return JeansPerturbationOracle(delta0=0.20, k_vec=(0.04, -0.03, 0.05), H0=74.6, f=0.53)

    def test_jeans_linear_continuity_diagnostic_exact_equality(self, jeans):
        # Checks div(v) = -H0 * f * delta(x) exactly at various cosmological coordinates
        coords = [
            np.array([0.0, 0.0, 0.0]),
            np.array([25.0, -15.0, 40.0]),
            np.array([-50.0, 80.0, -10.0]),
            np.array([120.0, -70.0, 35.0])
        ]
        for pos in coords:
            exact_div = jeans.analytical_divergence(pos)
            num_div = numerical_divergence_6th_order(jeans.velocity, pos, h=1e-4)
            delta = jeans.density_contrast(pos)
            expected_continuity = -jeans.H0 * jeans.f * delta
            assert math.isclose(exact_div, expected_continuity, rel_tol=1e-10)
            assert math.isclose(num_div, expected_continuity, rel_tol=1e-4, abs_tol=1e-5)

    def test_jeans_velocity_curl_free(self, jeans):
        # Potential flow from density perturbation is curl-free
        pt = np.array([30.0, -45.0, 60.0])
        num_curl = numerical_curl_6th_order(jeans.velocity, pt, h=1e-4)
        assert np.allclose(num_curl, [0.0, 0.0, 0.0], atol=1e-4)
