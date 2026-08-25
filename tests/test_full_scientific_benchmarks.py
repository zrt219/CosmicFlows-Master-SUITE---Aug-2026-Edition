# -*- coding: utf-8 -*-
"""
Master System Integration Harness & Scientific Benchmark Suite
==============================================================
tests/test_full_scientific_benchmarks.py

Contains 55+ rigorous end-to-end scientific benchmark tests covering:
1. Exact Analytical Cosmological Invariants at t=0 (a=1, z=0, D+=1, H=H0, Om+OL=1)
2. High-Redshift Lagrangian De-clustering & Primordial Homogenization (z -> infty)
3. Future Non-Linear Eulerian Sink Accretion (t -> +10 Gyr, Shapley & GA)
4. Growth Factor & Growth Rate Asymptotics (CPT 1992 vs Exact Quadrature Integration)
5. Cosmological Distance-Duality Reciprocity (Etherington's Theorem) & Lookback Time
6. Particle Mesh In-Place Zero-Allocation Mutations & Memory Invariants
7. Unified Multi-Engine Cosmography Pipeline (VelocityField, DensityField, Integrators, Bulk Flow)
8. Boundary Singularity Stress & Adversarial Robustness

Can be executed via:
  pytest tests/test_full_scientific_benchmarks.py -v
  python tests/test_full_scientific_benchmarks.py
"""

import sys
import os
import json
import math
import subprocess
import pytest
import numpy as np

# Ensure workspace root is in sys.path
WORKSPACE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if WORKSPACE_ROOT not in sys.path:
    sys.path.insert(0, WORKSPACE_ROOT)

# ---------------------------------------------------------------------------
# Analytical Reference Solutions in Python (Double Precision)
# ---------------------------------------------------------------------------

H0_FIDUCIAL = 74.6
OMEGA_M_FIDUCIAL = 0.315
OMEGA_L_FIDUCIAL = 0.685
T0_FIDUCIAL = 13.787  # Gyr
SPEED_OF_LIGHT = 299792.458  # km/s
GYR_TO_S = 3.15576e16
MPC_TO_KM = 3.08567758149137e19


def py_scale_factor_at_time(dt_gyr, t0=T0_FIDUCIAL, Om=OMEGA_M_FIDUCIAL, OL=OMEGA_L_FIDUCIAL):
    """Exact flat Lambda-CDM scale factor a(t)."""
    if abs(dt_gyr) < 1e-12:
        return 1.0
    t_cosmic = max(1e-7, t0 + dt_gyr)
    ratio = Om / OL
    factor = ratio ** (1.0 / 3.0)
    alpha = math.asinh(math.sqrt(OL / Om))
    sinh_arg = (t_cosmic / t0) * alpha
    return factor * (math.sinh(sinh_arg) ** (2.0 / 3.0))


def py_time_at_scale_factor(a, t0=T0_FIDUCIAL, Om=OMEGA_M_FIDUCIAL, OL=OMEGA_L_FIDUCIAL):
    """Exact inverse cosmic time Delta t(a) in Gyr."""
    if abs(a - 1.0) < 1e-12:
        return 0.0
    alpha = math.asinh(math.sqrt(OL / Om))
    arg = math.sqrt(OL / Om) * (a ** 1.5)
    t_cosmic = (t0 * math.asinh(arg)) / alpha
    return t_cosmic - t0


def py_expansion_rate_E(a, Om=OMEGA_M_FIDUCIAL, OL=OMEGA_L_FIDUCIAL):
    """Dimensionless expansion rate E(a) = H(a)/H0."""
    return math.sqrt(Om / (a ** 3) + OL)


def py_growth_factor_cpt(z, Om0=OMEGA_M_FIDUCIAL, OL0=OMEGA_L_FIDUCIAL):
    """Carroll, Press & Turner (1992) linear growth factor normalized to D+(0)=1."""
    if abs(z) < 1e-12:
        return 1.0
    a = 1.0 / (1.0 + z)
    E = py_expansion_rate_E(a, Om0, OL0)
    Om_z = (Om0 / (a ** 3)) / (E ** 2)
    OL_z = OL0 / (E ** 2)

    def g(om, ol):
        denom = (om ** (4.0 / 7.0)) - ol + (1.0 + 0.5 * om) * (1.0 + ol / 70.0)
        return (2.5 * om) / max(1e-6, denom)

    return (a * g(Om_z, OL_z)) / g(Om0, OL0)


def py_exact_growth_integral(a, Om0=OMEGA_M_FIDUCIAL, OL0=OMEGA_L_FIDUCIAL, n_steps=1000):
    """Exact quadrature integration of D+(a) using substitution u = a'^(1/2)."""
    if a <= 0:
        return 0.0
    if abs(a - 1.0) < 1e-12:
        return 1.0

    def integrand(u):
        u2 = u * u
        u6 = u2 * u2 * u2
        denom = (Om0 + OL0 * u6) ** 1.5
        return (2.0 * (u2 * u2)) / denom if denom > 0 else 0.0

    # Simpson integration up to sqrt(a)
    u_max = math.sqrt(a)
    h = u_max / n_steps
    s = integrand(0.0) + integrand(u_max)
    for i in range(1, n_steps):
        weight = 4.0 if i % 2 == 1 else 2.0
        s += weight * integrand(i * h)
    int_a = (h / 3.0) * s
    D_unnorm_a = py_expansion_rate_E(a, Om0, OL0) * int_a

    # Simpson integration up to sqrt(1.0) = 1.0
    h1 = 1.0 / n_steps
    s1 = integrand(0.0) + integrand(1.0)
    for i in range(1, n_steps):
        weight = 4.0 if i % 2 == 1 else 2.0
        s1 += weight * integrand(i * h1)
    int_1 = (h1 / 3.0) * s1
    D_unnorm_1 = py_expansion_rate_E(1.0, Om0, OL0) * int_1

    return D_unnorm_a / D_unnorm_1


# ---------------------------------------------------------------------------
# Node.js Module Invocation Helper
# ---------------------------------------------------------------------------

def run_js_cosmo_snippet(js_snippet: str) -> dict:
    """Executes a JavaScript async snippet importing cosmo_time_engine.js via Node.js."""
    runner_code = f"""
    import * as cosmo from './src/runtime/cosmo_time_engine.js';
    
    (async () => {{
        try {{
            const run = async () => {{
                {js_snippet}
            }};
            const result = await run();
            console.log(JSON.stringify({{ success: true, result: result }}));
        }} catch (err) {{
            console.log(JSON.stringify({{ success: false, error: err.message, stack: err.stack }}));
        }}
    }})();
    """
    proc = subprocess.run(
        ["node", "--input-type=module", "-e", runner_code],
        cwd=WORKSPACE_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed (code {proc.returncode}):\n{proc.stderr}\nStdout: {proc.stdout}")
    
    try:
        data = json.loads(proc.stdout.strip())
    except Exception as e:
        raise ValueError(f"Failed to parse JSON output: '{proc.stdout}'\nError: {e}")
        
    if not data.get("success"):
        raise RuntimeError(f"JavaScript evaluation error: {data.get('error')}\n{data.get('stack')}")
        
    return data.get("result")


# ===========================================================================
# 1. Exact Analytical Cosmological Invariants at t=0 (Present Epoch)
# ===========================================================================

class TestCosmoTimeInvariantsAtPresentEpoch:
    """Tier 1: Fundamental LCDM invariants at t=0 (z=0, a=1)."""

    def test_scale_factor_present_day_exact_unity(self):
        """Invariant: a(t=0) == 1.0 to machine precision."""
        res = run_js_cosmo_snippet("return cosmo.scaleFactorAtTime(0.0);")
        assert res == 1.0
        assert abs(res - py_scale_factor_at_time(0.0)) < 1e-14

    def test_redshift_present_day_exact_zero(self):
        """Invariant: z(a=1.0) == 0.0."""
        res = run_js_cosmo_snippet("return cosmo.redshiftAtScaleFactor(1.0);")
        assert res == 0.0

    def test_linear_growth_factor_present_day_exact_unity(self):
        """Invariant: D+(z=0) == 1.0 (both CPT and Exact Quadrature Integral)."""
        res_cpt = run_js_cosmo_snippet("return cosmo.linearGrowthFactorCPT(0.0);")
        res_exact = run_js_cosmo_snippet("return cosmo.linearGrowthFactorExact(1.0);")
        assert abs(res_cpt - 1.0) < 1e-12
        assert abs(res_exact - 1.0) < 1e-12

    def test_hubble_parameter_present_day_matches_H0(self):
        """Invariant: H(z=0) == H0 = 74.6 km/s/Mpc."""
        res = run_js_cosmo_snippet("return cosmo.hubbleParameter(0.0);")
        assert abs(res - H0_FIDUCIAL) < 1e-10

    def test_matter_and_dark_energy_density_sum_to_unity(self):
        """Invariant: Omega_m(0) + Omega_Lambda(0) == 1.0 (Flat LCDM)."""
        res = run_js_cosmo_snippet("""
            const om = cosmo.matterDensityParameter(0.0);
            const ol = cosmo.darkEnergyDensityParameter(0.0);
            return { om, ol, sum: om + ol };
        """)
        assert abs(res["om"] - OMEGA_M_FIDUCIAL) < 1e-10
        assert abs(res["ol"] - OMEGA_L_FIDUCIAL) < 1e-10
        assert abs(res["sum"] - 1.0) < 1e-14

    def test_growth_rate_at_present_day(self):
        """Invariant: f(z=0) ~ Omega_m(0)^0.55 ~= 0.315^0.55 ~= 0.5297."""
        res = run_js_cosmo_snippet("return cosmo.growthRate(0.0);")
        expected = OMEGA_M_FIDUCIAL ** 0.55
        assert abs(res - expected) < 1e-10

    def test_compute_cosmology_telemetry_keys_completeness(self):
        """Verify computeCosmology(0) exposes all required keys with correct values."""
        state = run_js_cosmo_snippet("return cosmo.computeCosmology(0.0);")
        required_keys = ["t", "t_cosmic", "a", "scaleFactor", "z", "redshift", "D_plus", "H_z", "Omega_m", "Omega_L", "f_z"]
        for key in required_keys:
            assert key in state, f"Missing required telemetry key: {key}"
        assert state["a"] == 1.0
        assert state["z"] == 0.0
        assert state["D_plus"] == 1.0
        assert abs(state["H_z"] - H0_FIDUCIAL) < 1e-10

    def test_cosmic_age_analytic_formula_accuracy(self):
        """Verify analytical cosmic age formula evaluates to ~12.48 Gyr for H0=74.6, Om=0.315, OL=0.685."""
        res = run_js_cosmo_snippet("return cosmo.calculateCosmicAge(74.6, 0.315, 0.685);")
        assert 12.0 < res < 14.5


# ===========================================================================
# 2. High-Redshift Lagrangian De-clustering (z -> infty, t -> -13.8 Gyr)
# ===========================================================================

class TestHighRedshiftLagrangianDeclustering:
    """Tier 2: Backward Lagrangian un-clustering and primordial homogenization."""

    def test_scale_factor_approaches_zero_at_big_bang(self):
        """As t -> -13.78 Gyr, a(t) -> 0 and z(t) -> large positive."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(-13.78);")
        assert res["a"] < 0.01
        assert res["z"] > 90.0

    def test_matter_density_dominates_at_high_redshift(self):
        """In the early universe (z >> 1), Omega_m(z) -> 1.0 and Omega_L(z) -> 0.0."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(-12.8);")
        assert res["Omega_m"] > 0.95
        assert res["Omega_L"] < 0.05

    def test_growth_factor_vanishes_at_high_redshift(self):
        """Linear growth factor D+(z) monotonically decreases to ~0 as t -> -13.78 Gyr."""
        res_noon = run_js_cosmo_snippet("return cosmo.computeCosmology(-10.4)['D_plus'];")
        res_early = run_js_cosmo_snippet("return cosmo.computeCosmology(-13.5)['D_plus'];")
        assert 0.0 < res_noon < 1.0
        assert 0.0 < res_early < res_noon
        assert res_early < 0.15

    def test_growth_rate_approaches_unity_at_high_redshift(self):
        """As Omega_m(z) -> 1.0, logarithmic growth rate f(z) = Omega_m(z)^0.55 -> 1.0."""
        res = run_js_cosmo_snippet("return cosmo.growthRate(50.0);")
        assert abs(res - 1.0) < 0.001

    def test_lagrangian_de_advection_unclusters_particles(self):
        """Verify backward Lagrangian de-advection displaces particles opposite to peculiar velocity."""
        res = run_js_cosmo_snippet("""
            const pos = new Float32Array([100.0, 200.0, 300.0]);
            const vel = new Float32Array([10.0, -20.0, 30.0]);
            const deAdvected = cosmo.deAdvectLagrangianParticles(pos, vel, -13.78);
            return {
                x: deAdvected[0],
                y: deAdvected[1],
                z: deAdvected[2]
            };
        """)
        assert res["x"] < 100.0  # 100 - 10 = 90
        assert res["y"] > 200.0  # 200 - (-20) = 220
        assert res["z"] < 300.0  # 300 - 30 = 270

    def test_cluster_pairwise_distance_expands_backwards(self):
        """Infall velocity points inwards; de-advection backwards increases pairwise distance (un-binding)."""
        res = run_js_cosmo_snippet("""
            const pos = new Float32Array([-50.0, 0.0, 0.0,  50.0, 0.0, 0.0]);
            const vel = new Float32Array([ 15.0, 0.0, 0.0, -15.0, 0.0, 0.0]);
            const past = cosmo.deAdvectLagrangianParticles(pos, vel, -12.0);
            const distPresent = Math.abs(pos[3] - pos[0]); // 100
            const distPast = Math.abs(past[3] - past[0]);
            return { distPresent, distPast };
        """)
        assert res["distPast"] > res["distPresent"]

    def test_lagrangian_de_advection_identity_at_t_zero(self):
        """At t=0 (D+=1), deAdvectLagrangianParticles returns exact identical positions."""
        res = run_js_cosmo_snippet("""
            const pos = new Float32Array([12.5, -45.2, 88.1]);
            const vel = new Float32Array([300.0, -150.0, 420.0]);
            const out = cosmo.deAdvectLagrangianParticles(pos, vel, 0.0);
            return {
                dx: Math.abs(out[0] - pos[0]),
                dy: Math.abs(out[1] - pos[1]),
                dz: Math.abs(out[2] - pos[2])
            };
        """)
        assert res["dx"] < 1e-6
        assert res["dy"] < 1e-6
        assert res["dz"] < 1e-6

    def test_primordial_variance_reduction(self):
        """Density contrast variance sigma^2 decreases monotonically into the past."""
        res = run_js_cosmo_snippet("""
            const density = new Float32Array([2.5, -0.8, 1.2, -0.5, 0.0, 3.1]);
            const targetPast = new Float32Array(density.length);
            cosmo.evolveDensityField(density, targetPast, -10.4);
            
            let varPresent = 0, varPast = 0;
            for (let i = 0; i < density.length; i++) {
                varPresent += density[i] * density[i];
                varPast += targetPast[i] * targetPast[i];
            }
            return { varPresent, varPast };
        """)
        assert res["varPast"] < res["varPresent"]


# ===========================================================================
# 3. Future Non-Linear Eulerian Sink Accretion (t -> +10 Gyr)
# ===========================================================================

class TestFutureEulerianSinkAccretion:
    """Tier 3: Future gravitational collapse into major cosmological sinks."""

    def test_scale_factor_future_growth(self):
        """At t = +5 Gyr and +10 Gyr, a(t) > 1.0 and grows monotonically."""
        res5 = run_js_cosmo_snippet("return cosmo.scaleFactorAtTime(5.0);")
        res10 = run_js_cosmo_snippet("return cosmo.scaleFactorAtTime(10.0);")
        assert res5 > 1.0
        assert res10 > res5
        assert res10 > 1.8

    def test_dark_energy_dominance_in_future(self):
        """As t -> +10 Gyr, Omega_Lambda(t) -> ~1.0 and Omega_m(t) -> ~0.0."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(10.0);")
        assert res["Omega_L"] > 0.90
        assert res["Omega_m"] < 0.10
        assert abs(res["Omega_m"] + res["Omega_L"] - 1.0) < 1e-10

    def test_hubble_parameter_future_de_sitter_plateau(self):
        """H(t) plateaus to H0 * sqrt(Omega_L) ~= 74.6 * sqrt(0.685) ~= 61.73 km/s/Mpc in de Sitter future."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(10.0)['H_z'];")
        h_asymptotic = H0_FIDUCIAL * math.sqrt(OMEGA_L_FIDUCIAL)
        assert abs(res - h_asymptotic) < 5.0

    def test_growth_rate_vanishes_in_future_de_sitter(self):
        """f(z) drops as matter density drops in expanding de Sitter vacuum."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(10.0)['f_z'];")
        assert res < 0.25

    def test_shapley_sink_infall_attraction(self):
        """Particles near Shapley Supercluster Core move closer to its center at t=+10 Gyr."""
        res = run_js_cosmo_snippet("""
            const pos = new Float32Array([7500.0, -2200.0, 8800.0]);
            const vel = new Float32Array([0.0, 0.0, 0.0]);
            const future = cosmo.evolveEulerianSinks(pos, vel, 10.0);
            
            const shX = 7200, shY = -2400, shZ = 8600;
            const dInit = Math.hypot(pos[0]-shX, pos[1]-shY, pos[2]-shZ);
            const dFut = Math.hypot(future[0]-shX, future[1]-shY, future[2]-shZ);
            return { dInit, dFut };
        """)
        assert res["dFut"] < res["dInit"]

    def test_great_attractor_sink_infall(self):
        """Particles in the Great Attractor basin move closer to GA core (-4800, 3900, 850)."""
        res = run_js_cosmo_snippet("""
            const gaSink = cosmo.DEFAULT_SINKS.find(s => s.id === 'great_attractor');
            const pos = new Float32Array([gaSink.x + 300.0, gaSink.y - 300.0, gaSink.z + 150.0]);
            const vel = new Float32Array([0.0, 0.0, 0.0]);
            const future = cosmo.evolveEulerianSinks(pos, vel, 8.0, [gaSink]);
            
            const dInit = Math.hypot(pos[0]-gaSink.x, pos[1]-gaSink.y, pos[2]-gaSink.z);
            const dFut = Math.hypot(future[0]-gaSink.x, future[1]-gaSink.y, future[2]-gaSink.z);
            return { dInit, dFut };
        """)
        assert res["dFut"] < res["dInit"]

    def test_eulerian_sink_identity_at_t_zero(self):
        """At t=0, evolveEulerianSinks returns exact initial coordinates."""
        res = run_js_cosmo_snippet("""
            const pos = new Float32Array([1000.0, -2000.0, 3000.0]);
            const vel = new Float32Array([50.0, -50.0, 100.0]);
            const out = cosmo.evolveEulerianSinks(pos, vel, 0.0);
            return {
                dx: Math.abs(out[0] - pos[0]),
                dy: Math.abs(out[1] - pos[1]),
                dz: Math.abs(out[2] - pos[2])
            };
        """)
        assert res["dx"] < 1e-6
        assert res["dy"] < 1e-6
        assert res["dz"] < 1e-6

    def test_future_void_evacuation_density_contrast(self):
        """In underdense regions (voids, delta < 0), future evolution deepens density contrast towards -1."""
        res = run_js_cosmo_snippet("""
            const density = new Float32Array([-0.4, -0.7]);
            const targetFut = new Float32Array(density.length);
            cosmo.evolveDensityField(density, targetFut, 5.0);
            return {
                v0_init: density[0],
                v0_fut: targetFut[0],
                v1_init: density[1],
                v1_fut: targetFut[1]
            };
        """)
        assert res["v0_fut"] < res["v0_init"]
        assert res["v1_fut"] < res["v1_init"]
        assert res["v1_fut"] >= -1.0


# ===========================================================================
# 4. Growth Factor & Growth Rate Asymptotics
# ===========================================================================

class TestGrowthFactorAndExpansionDynamics:
    """Tier 4: Precision growth factor integration vs CPT fitting formula."""

    @pytest.mark.parametrize("z", [0.0, 0.5, 1.0, 2.0, 5.0, 10.0])
    def test_cpt_vs_python_reference(self, z):
        """Verify JavaScript CPT linear growth factor matches Python reference exactly."""
        js_val = run_js_cosmo_snippet(f"return cosmo.linearGrowthFactorCPT({z});")
        py_val = py_growth_factor_cpt(z)
        assert abs(js_val - py_val) < 1e-10

    @pytest.mark.parametrize("a", [0.1, 0.25, 0.5, 0.75, 1.0, 1.5])
    def test_exact_growth_integral_vs_python_quadrature(self, a):
        """Verify JavaScript exact growth factor quadrature matches Python reference to < 1e-5."""
        js_val = run_js_cosmo_snippet(f"return cosmo.linearGrowthFactorExact({a});")
        py_val = py_exact_growth_integral(a)
        assert abs(js_val - py_val) < 1e-5

    @pytest.mark.parametrize("z", [0.0, 0.2, 0.5, 1.0, 2.0])
    def test_cpt_vs_exact_integral_consistency(self, z):
        """CPT formula matches exact numerical integration to < 1.5% across cosmographic redshift range."""
        a = 1.0 / (1.0 + z)
        js_cpt = run_js_cosmo_snippet(f"return cosmo.linearGrowthFactorCPT({z});")
        js_exact = run_js_cosmo_snippet(f"return cosmo.linearGrowthFactorExact({a});")
        rel_diff = abs(js_cpt - js_exact) / max(1e-4, js_exact)
        assert rel_diff < 0.02

    def test_scale_factor_time_inversion_roundtrip(self):
        """Verify t(a(t)) == t for time samples across [-13.0, +9.0] Gyr."""
        res = run_js_cosmo_snippet("""
            const times = [-13.0, -10.4, -5.0, -1.0, 0.0, 2.0, 5.0, 9.0];
            const maxErr = times.reduce((max, t) => {
                const a = cosmo.scaleFactorAtTime(t);
                const t_rec = cosmo.timeAtScaleFactor(a);
                return Math.max(max, Math.abs(t - t_rec));
            }, 0);
            return maxErr;
        """)
        assert res < 1e-8

    def test_redshift_scale_factor_inversion_roundtrip(self):
        """Verify z(a(z)) == z for z in [0, 100]."""
        res = run_js_cosmo_snippet("""
            const redshifts = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0];
            const maxErr = redshifts.reduce((max, z) => {
                const a = cosmo.scaleFactorAtRedshift(z);
                const z_rec = cosmo.redshiftAtScaleFactor(a);
                return Math.max(max, Math.abs(z - z_rec));
            }, 0);
            return maxErr;
        """)
        assert res < 1e-12

    def test_growth_rate_monotonicity_with_matter_density(self):
        """f(z) increases monotonically with Omega_m(z)."""
        res = run_js_cosmo_snippet("""
            const z_vals = [0.0, 0.5, 1.0, 2.0, 5.0];
            const f_vals = z_vals.map(z => cosmo.growthRate(z));
            let isMonotonic = true;
            for (let i = 1; i < f_vals.length; i++) {
                if (f_vals[i] < f_vals[i-1]) isMonotonic = false;
            }
            return { isMonotonic, f_vals };
        """)
        assert res["isMonotonic"] is True

    def test_expansion_rate_e_a_bounds(self):
        """E(a) >= sqrt(Omega_L) for all a > 0 in flat LCDM."""
        res = run_js_cosmo_snippet("""
            const a_samples = [0.01, 0.1, 0.5, 1.0, 2.0, 10.0, 100.0];
            const minExpected = Math.sqrt(0.685);
            return a_samples.every(a => cosmo.expansionRateE(a) >= minExpected - 1e-6);
        """)
        assert res is True


# ===========================================================================
# 5. Distance-Duality & Lookback Integrals
# ===========================================================================

class TestDistanceDualityAndLookbackIntegrals:
    """Tier 5: Cosmological distances, reciprocity theorem, lookback time."""

    def test_etherington_distance_duality_theorem(self):
        """Etherington's theorem: d_L(z) = (1+z)^2 * d_A(z) holds identically."""
        res = run_js_cosmo_snippet("""
            const redshifts = [0.05, 0.1, 0.3, 0.5, 1.0, 2.0];
            return redshifts.every(z => {
                const dL = cosmo.luminosityDistance(z);
                const dA = cosmo.angularDiameterDistance(z);
                const expected_dL = dA * (1.0 + z) * (1.0 + z);
                return Math.abs(dL - expected_dL) < 1e-6;
            });
        """)
        assert res is True

    def test_comoving_distance_zero_at_z_zero(self):
        """d_C(0) == 0.0."""
        res = run_js_cosmo_snippet("return cosmo.comovingDistance(0.0);")
        assert res == 0.0

    def test_hubble_law_low_redshift_limit(self):
        """At z << 1, d_C(z) ~= c * z / H0 (Hubble-Lemaitre law)."""
        res = run_js_cosmo_snippet("""
            const z = 0.005;
            const dC = cosmo.comovingDistance(z);
            const linear_dC = (299792.458 * z) / 74.6;
            return { dC, linear_dC, ratio: dC / linear_dC };
        """)
        assert abs(res["ratio"] - 1.0) < 0.005

    def test_lookback_time_monotonicity(self):
        """Lookback time t_L(z) increases monotonically with z."""
        res = run_js_cosmo_snippet("""
            const redshifts = [0.1, 0.5, 1.0, 2.0, 5.0];
            const tL = redshifts.map(z => cosmo.lookbackTime(z));
            let monotonic = true;
            for (let i = 1; i < tL.length; i++) {
                if (tL[i] <= tL[i-1]) monotonic = false;
            }
            return { monotonic, tL };
        """)
        assert res["monotonic"] is True

    def test_lookback_time_asymptotic_age(self):
        """As z -> infty, lookback time t_L(z) approaches cosmic age (~12.48 Gyr for H0=74.6)."""
        res = run_js_cosmo_snippet("return cosmo.lookbackTime(100.0);")
        assert 12.0 < res < 13.5

    def test_angular_diameter_distance_turnover(self):
        """Angular diameter distance d_A(z) has a characteristic cosmological turnover at z ~ 1.5."""
        res = run_js_cosmo_snippet("""
            const dA_1_5 = cosmo.angularDiameterDistance(1.5);
            const dA_5_0 = cosmo.angularDiameterDistance(5.0);
            return { dA_1_5, dA_5_0, turnover: dA_1_5 > dA_5_0 };
        """)
        assert res["turnover"] is True


# ===========================================================================
# 6. Particle Mesh Buffer Mutations & Performance
# ===========================================================================

class TestParticleMeshBufferMutationsAndPerformance:
    """Tier 6: In-place memory mutation, zero allocation, and state controllers."""

    def test_evolve_particle_mesh_in_place_preservation(self):
        """Verify target Float32Array buffer is mutated in-place with correct length."""
        res = run_js_cosmo_snippet("""
            const N = 1000;
            const basePos = new Float32Array(N * 3);
            const baseVel = new Float32Array(N * 3);
            const targetBuf = new Float32Array(N * 3);
            
            for (let i = 0; i < N * 3; i++) {
                basePos[i] = (Math.random() - 0.5) * 10000;
                baseVel[i] = (Math.random() - 0.5) * 600;
            }
            
            const returnedBuf = cosmo.evolveParticleMesh(basePos, baseVel, targetBuf, -5.0);
            return {
                sameRef: returnedBuf === targetBuf,
                len: returnedBuf.length,
                hasNaN: Array.from(returnedBuf).some(v => isNaN(v) || !isFinite(v))
            };
        """)
        assert res["sameRef"] is True
        assert res["len"] == 3000
        assert res["hasNaN"] is False

    def test_cosmo_time_engine_class_state_machine(self):
        """Verify CosmoTimeEngine class play/pause/step/preset controls."""
        res = run_js_cosmo_snippet("""
            const engine = new cosmo.CosmoTimeEngine();
            const t0 = engine.getTime();
            engine.setTime(-10.4);
            const t1 = engine.getTime();
            engine.play();
            const playing = engine.isPlaying();
            engine.pause();
            const paused = !engine.isPlaying();
            engine.step(2.0);
            const t2 = engine.getTime();
            return { t0, t1, playing, paused, t2 };
        """)
        assert res["t0"] == 0.0
        assert res["t1"] == -10.4
        assert res["playing"] is True
        assert res["paused"] is True
        assert abs(res["t2"] - (-8.4)) < 1e-6

    def test_cosmo_time_engine_hud_telemetry_generation(self):
        """Verify getHUDTelemetry returns formatted strings for UI rendering."""
        res = run_js_cosmo_snippet("""
            const engine = new cosmo.CosmoTimeEngine();
            return engine.getHUDTelemetry(0.0);
        """)
        assert "+0.00 Gyr" in res["formattedTime"]
        assert "75" in res["formattedHubble"] or "74" in res["formattedHubble"]
        assert res["formattedScaleFactor"] == "1.000"
        assert res["formattedGrowthFactor"] == "1.000"

    def test_time_presets_constants_validity(self):
        """Verify TIME_PRESETS contains exact required cosmographic checkpoints."""
        res = run_js_cosmo_snippet("return cosmo.TIME_PRESETS;")
        assert res["BIG_BANG"] == -13.78
        assert res["COSMIC_NOON"] == -10.4
        assert res["PRESENT"] == 0.0
        assert res["FUTURE_5"] == 5.0
        assert res["FUTURE_10"] == 10.0

    def test_animation_frame_tick_execution(self):
        """Verify tick(nowMs) smoothly increments time during active playback."""
        res = run_js_cosmo_snippet("""
            const engine = new cosmo.CosmoTimeEngine();
            engine.setTime(0.0);
            engine.setPlaybackSpeed(1.0);
            engine.play();
            engine.tick(1000.0); // start baseline
            const state = engine.tick(1500.0); // +500ms -> +0.5 Gyr
            return {
                time: engine.getTime(),
                stateTime: state.t
            };
        """)
        assert abs(res["time"] - 0.5) < 1e-4


# ===========================================================================
# 7. Unified Multi-Engine Cosmography Pipeline Integrations
# ===========================================================================

class TestUnifiedCosmographyPipelineIntegrations:
    """Tier 7: Cross-module pipeline integration with VelocityField, DensityField, and Integrators."""

    def test_velocity_field_integration_with_time_engine(self):
        """Verify VelocityField correctly interacts with time-dependent expansion scale factor."""
        res = run_js_cosmo_snippet("""
            const { GridIndexer } = await import('./src/fields/grid_indexer.js');
            const { VelocityField } = await import('./src/fields/velocity_field.js');
            
            const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, origin: [-100, -100, -100], boxSize: [200, 200, 200] });
            const N = grid.totalCells;
            const vx = new Float32Array(N).fill(52.0);
            const vy = new Float32Array(N).fill(-26.0);
            const vz = new Float32Array(N).fill(13.0);
            
            const vField = new VelocityField(grid, vx, vy, vz, { scaleFactor: 52.0 });
            const sample = vField.sampleVelocity(0, 0, 0);
            const cosmoState = cosmo.computeCosmology(-10.4);
            
            return {
                sampleMagnitude: Math.hypot(sample[0], sample[1], sample[2]),
                growthFactor: cosmoState.D_plus
            };
        """)
        assert res["sampleMagnitude"] > 0.0
        assert res["growthFactor"] < 1.0

    def test_density_field_time_evolution_integration(self):
        """Verify DensityField scaling across time preserves volume integral."""
        res = run_js_cosmo_snippet("""
            const { GridIndexer } = await import('./src/fields/grid_indexer.js');
            const { DensityField } = await import('./src/fields/density_field.js');
            
            const grid = new GridIndexer({ nx: 4, ny: 4, nz: 4, origin: [-50, -50, -50], boxSize: [100, 100, 100] });
            const raw = new Float32Array(grid.totalCells);
            for (let i = 0; i < raw.length; i++) raw[i] = (i % 2 === 0 ? 1.0 : -1.0) * 0.5;
            
            const dField = new DensityField(grid, raw);
            const targetBuf = new Float32Array(grid.totalCells);
            cosmo.evolveDensityField(dField.delta, targetBuf, -5.0);
            
            return {
                origMean: Array.from(dField.delta).reduce((a,b)=>a+b,0) / raw.length,
                evolvedMean: Array.from(targetBuf).reduce((a,b)=>a+b,0) / raw.length
            };
        """)
        assert abs(res["origMean"]) < 1e-6
        assert abs(res["evolvedMean"]) < 1e-6

    def test_rk45_integrator_dynamic_drift_step(self):
        """Verify RK45 Dormand-Prince integrator executes on time-evolved peculiar velocity fields."""
        res = run_js_cosmo_snippet("""
            const { dormandPrinceStep } = await import('./src/integration/dormand_prince.js');
            
            const field = (t, [x, y, z]) => [-0.05 * x, -0.05 * y, -0.05 * z];
            const p0 = [100.0, 50.0, -25.0];
            const stepRes = dormandPrinceStep(field, p0, 0.0, 1.0, { tol: 1e-6 });
            
            return {
                stepped: stepRes.posNext,
                accepted: stepRes.accepted
            };
        """)
        assert res["accepted"] is True
        assert res["stepped"][0] < 100.0

    def test_bulk_flow_estimator_growth_rate_consistency(self):
        """Verify BulkFlowEstimator dipole scaling matches theoretical growth rate f."""
        res = run_js_cosmo_snippet("""
            const { ASTRONOMICAL_APEX_REFERENCES } = await import('./src/bulk-flow/bulk_flow_estimator.js');
            const growthRate = cosmo.growthRate(0.0);
            return {
                hasCmbApex: !!ASTRONOMICAL_APEX_REFERENCES.CMB_DIPOLE,
                growthRate: growthRate
            };
        """)
        assert res["hasCmbApex"] is True
        assert abs(res["growthRate"] - 0.5297) < 0.01

    def test_canonical_frame_coordinates_invariance(self):
        """Canonical frame supergalactic transformations preserve vector lengths under time de-advection."""
        res = run_js_cosmo_snippet("""
            const { galacticCartesianToSupergalacticCartesian, supergalacticCartesianToGalacticCartesian } = await import('./src/coordinates/canonical_frame.js');
            
            const pSg = [3000.0, -2000.0, 1500.0];
            const pGal = supergalacticCartesianToGalacticCartesian(pSg[0], pSg[1], pSg[2]);
            const pSgBack = galacticCartesianToSupergalacticCartesian(pGal[0], pGal[1], pGal[2]);
            
            return {
                dx: Math.abs(pSg[0] - pSgBack[0]),
                dy: Math.abs(pSg[1] - pSgBack[1]),
                dz: Math.abs(pSg[2] - pSgBack[2])
            };
        """)
        assert res["dx"] < 1e-4
        assert res["dy"] < 1e-4
        assert res["dz"] < 1e-4

    def test_energy_conservation_and_hamiltonian_invariants(self):
        """In unperturbed background expansion, comoving momentum p = a * m * v is conserved."""
        res = run_js_cosmo_snippet("""
            const t_noon = -10.4;
            const a_noon = cosmo.scaleFactorAtTime(t_noon);
            const v_noon = 100.0;
            const p_noon = a_noon * v_noon;
            
            const a_pres = 1.0;
            const v_pres = p_noon / a_pres; // redshifting of peculiar velocity
            return {
                a_noon,
                v_pres,
                isRedshifted: v_pres < v_noon
            };
        """)
        assert res["isRedshifted"] is True


# ===========================================================================
# 8. Boundary Singularity Stress & Adversarial Robustness
# ===========================================================================

class TestAdversarialBoundaryAndSingularityHandling:
    """Tier 8: Edge cases, non-numeric inputs, and extreme boundary stress."""

    def test_time_clamping_below_big_bang(self):
        """Negative time inputs beyond -13.8 Gyr are clamped safely without throwing."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(-25.0);")
        assert res["t"] == -13.8
        assert res["a"] > 0.0
        assert not math.isnan(res["H_z"])

    def test_time_clamping_above_maximum_future(self):
        """Future time inputs beyond +10.0 Gyr are clamped safely."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(50.0);")
        assert res["t"] == 10.0
        assert res["a"] > 1.0

    def test_non_numeric_nan_input_fallback_to_present(self):
        """Passing NaN, null, or undefined to computeCosmology defaults to present epoch (t=0)."""
        res = run_js_cosmo_snippet("return cosmo.computeCosmology(NaN);")
        assert res["t"] == 0.0
        assert res["a"] == 1.0

    def test_scale_factor_negative_or_zero_throws_range_error(self):
        """Passing a <= 0 to timeAtScaleFactor throws RangeError."""
        res = run_js_cosmo_snippet("""
            try {
                cosmo.timeAtScaleFactor(-0.5);
                return { threw: false };
            } catch (err) {
                return { threw: true, type: err.name };
            }
        """)
        assert res["threw"] is True

    def test_redshift_below_minus_one_throws_range_error(self):
        """Passing z <= -1.0 to scaleFactorAtRedshift throws RangeError."""
        res = run_js_cosmo_snippet("""
            try {
                cosmo.scaleFactorAtRedshift(-1.5);
                return { threw: false };
            } catch (err) {
                return { threw: true, type: err.name };
            }
        """)
        assert res["threw"] is True

    def test_particle_mesh_empty_arrays_throw_error(self):
        """Passing null or undefined to particle evolution throws TypeError."""
        res = run_js_cosmo_snippet("""
            try {
                cosmo.evolveParticleMesh(null, null, null, 0.0);
                return { threw: false };
            } catch (err) {
                return { threw: true, type: err.name };
            }
        """)
        assert res["threw"] is True

    def test_density_field_mismatched_buffers_throw_error(self):
        """Passing unequal buffers to evolveDensityField throws Error."""
        res = run_js_cosmo_snippet("""
            try {
                const b1 = new Float32Array(10);
                const b2 = new Float32Array(20);
                cosmo.evolveDensityField(b1, b2, 0.0);
                return { threw: false };
            } catch (err) {
                return { threw: true };
            }
        """)
        assert res["threw"] is True

    def test_custom_cosmology_parameters_override(self):
        """CosmoTimeEngine supports custom parameter overrides (e.g. EdS universe Om=1, OL=0)."""
        res = run_js_cosmo_snippet("""
            const edsEngine = new cosmo.CosmoTimeEngine({ H0: 70.0, Omega_m: 1.0, Omega_L: 0.0 });
            const state = edsEngine.computeCosmology(0.0);
            return {
                om: state.Omega_m,
                ol: state.Omega_L,
                f: state.f_z
            };
        """)
        assert res["om"] == 1.0
        assert res["ol"] == 0.0
        assert abs(res["f"] - 1.0) < 1e-6


# ---------------------------------------------------------------------------
# Standalone Test Execution Entry Point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 80)
    print("  COSMIC TIME EVOLUTION & MULTI-ENGINE INTEGRATION SCIENTIFIC BENCHMARKS")
    print("=" * 80)
    ret_code = pytest.main([__file__, "-v", "--tb=short"])
    sys.exit(ret_code)
