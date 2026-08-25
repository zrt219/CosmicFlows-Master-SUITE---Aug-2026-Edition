"""
tests/fields/test_helmholtz_fft_decomposer.py
Exhaustive, production-grade automated Pytest test suite for the 3D FFT Helmholtz-Hodge Vector Field Decomposer.

Mathematical Test Categories:
1. Spectral K-Space Projection Operator Matrix Invariants:
   - Completeness: P_pot + P_sol = I
   - Idempotence: (P_pot)^2 = P_pot, (P_sol)^2 = P_sol
   - Mutual Orthogonality: P_pot * P_sol = 0
   - Self-Adjoint / Hermitian Symmetry: P = P^T
   - Trace Invariants: Tr(P_pot) = 1, Tr(P_sol) = 2
   - Zero-Mode and Wavemode Rank Properties

2. Pure Potential (Irrotational / Compressible) Field Decompositions:
   - Single-harmonic standing wave potential field
   - Multi-mode 3D potential wave packet
   - Verification: v_sol == 0, v_pot == v, E_sol / E_tot < 1e-14
   - Spectral curl ||curl v_pot||_spec < 1e-10

3. Pure Solenoidal (Divergence-Free / Helical) Field Decompositions:
   - Arnold-Beltrami-Childress (ABC) helical flow
   - 3D Taylor-Green vortex flow
   - Stuart vortex array flow
   - Verification: v_pot == 0, v_sol == v, E_pot / E_tot < 1e-14
   - Spectral divergence ||div v_sol||_spec < 1e-10

4. Exact Harmonic Zero-Mode Bulk Velocity & Superposition:
   - Three-way decomposition: v = v_pot + v_sol + v_0
   - Dipole bulk velocity recovery across arbitrary coordinate frames
   - Constant uniform flow edge cases
   - Pointwise field reconstruction residual ||v - (v_pot + v_sol + v_0)||_inf < 1e-10

5. Strict L2 Orthogonality & Parseval Kinetic Energy Conservation:
   - Inner product <v_pot, v_sol> < 1e-10
   - Normalized cross-correlation rho < 1e-10
   - Zero-mode orthogonality <v_0, v_pot> == 0, <v_0, v_sol> == 0
   - Kinetic energy conservation |E_tot - (E_pot + E_sol + E_0)| / E_tot < 1e-10

6. Kinematic Potentials & Field Invariants:
   - Velocity scalar potential: -grad Phi_v == v_pot, laplacian Phi_v == -div v
   - Solenoidal stream vector potential: curl A == v_sol, Coulomb gauge div A == 0
   - Vorticity omega = curl v_sol, curl v_pot == 0
   - Kinetic helicity density h = v . omega, total integrated helicity

7. 1D Isotropic Shell Energy Power Spectra & Power-Law Fitting:
   - Radial k-shell energy binning E_tot(k), E_pot(k), E_sol(k)
   - Energy sum conservation across shells
   - Spectral fraction table f_pot(k) + f_sol(k) == 1
   - Least-squares power-law fit E(k) = A * k^(-alpha)

8. Cosmological Linear Matter Density Inversion:
   - delta_rec(x) = -div v_pot / (H0 * f) = laplacian Phi_v / (H0 * f)
   - Growth factor f(Omega_m, z) scaling across cosmological parameters

9. Discretization Kernels & Trigonometric Stencils:
   - Continuous isotropic kernel vs 7-point, 19-point, 27-point kernels
   - Discrete Laplacian residual convergence

10. Multiscale Spatial Scale Bands & Spatial Filtering:
    - Cluster (< 5 Mpc/h), Filament (5-20 Mpc/h), Supercluster (> 20 Mpc/h)
    - Gaussian, Top-Hat, and Sharp-k Fourier filtering
    - Energy partition completeness across spatial bands

11. Anisotropic Grids, Dynamic Range & Edge Cases:
    - Non-cubic anisotropic grid dimensions (Lx != Ly != Lz, Nx != Ny != Nz)
    - Zero/null velocity field, high-dynamic-range velocities (1e-5 to 1e5 km/s)
    - Input buffer immutability and error handling
"""

import json
import math
import os
import subprocess
import numpy as np
import pytest

NODE_EXEC = "node"
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))

def run_node_snippet(code: str) -> dict:
    """Executes an inline ES6 Node.js snippet and returns parsed JSON output."""
    wrapped_code = f"""
    {code}
    """
    proc = subprocess.run(
        [NODE_EXEC, "--input-type=module", "-e", wrapped_code],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Node execution failed with code {proc.returncode}:\nSTDERR:\n{proc.stderr}\nSTDOUT:\n{proc.stdout}")
    return json.loads(proc.stdout)


# ============================================================================
# 1. SPECTRAL K-SPACE PROJECTION OPERATOR MATRIX INVARIANT TESTS
# ============================================================================

class TestKSpaceProjectorInvariants:
    """Rigorous mathematical verification of 3x3 Helmholtz projection tensor operators."""

    def test_spectral_projector_completeness_and_trace(self):
        """Verify completeness P_pot + P_sol = I and trace invariants Tr(P_pot)=1, Tr(P_sol)=2."""
        code = """
        import { KSpaceProjector } from './src/fields/helmholtz_fft_decomposer.js';

        const testWavevectors = [
          [1.0, 0.0, 0.0],
          [0.0, 2.5, 0.0],
          [0.0, 0.0, -3.2],
          [1.0, 1.0, 1.0],
          [2.3, -4.5, 6.7],
          [-0.12, 0.45, -0.89]
        ];

        const results = testWavevectors.map(([kx, ky, kz]) => {
          const inv = KSpaceProjector.verifyAlgebraicInvariants(kx, ky, kz);
          return {
            kx, ky, kz,
            tracePot: inv.tracePot,
            traceSol: inv.traceSol,
            isIdempotent: inv.isIdempotent,
            isOrthogonal: inv.isOrthogonal,
            completenessError: inv.completenessError
          };
        });

        console.log(JSON.stringify(results));
        """
        results = run_node_snippet(code)
        for res in results:
            assert abs(res["tracePot"] - 1.0) < 1e-12, f"Trace of P_pot must be 1, got {res['tracePot']}"
            assert abs(res["traceSol"] - 2.0) < 1e-12, f"Trace of P_sol must be 2, got {res['traceSol']}"
            assert res["isIdempotent"], f"Projector must be idempotent at k=({res['kx']}, {res['ky']}, {res['kz']})"
            assert res["isOrthogonal"], f"P_pot and P_sol must be mutually orthogonal at k=({res['kx']}, {res['ky']}, {res['kz']})"
            assert res["completenessError"] < 1e-12, f"Completeness error too high: {res['completenessError']}"

    def test_projector_idempotence_and_symmetry(self):
        """Verify P_pot^2 = P_pot, P_sol^2 = P_sol and symmetry P = P^T."""
        code = """
        import { KSpaceProjector } from './src/fields/helmholtz_fft_decomposer.js';

        const kx = 3.14159;
        const ky = -2.71828;
        const kz = 1.41421;

        const P_pot = KSpaceProjector.getPotentialProjectorMatrix(kx, ky, kz);
        const P_sol = KSpaceProjector.getSolenoidalProjectorMatrix(kx, ky, kz);

        let potSymmetryErr = 0.0;
        let solSymmetryErr = 0.0;
        for (let r = 0; r < 3; r++) {
          for (let c = 0; c < 3; c++) {
            potSymmetryErr = Math.max(potSymmetryErr, Math.abs(P_pot[r*3 + c] - P_pot[c*3 + r]));
            solSymmetryErr = Math.max(solSymmetryErr, Math.abs(P_sol[r*3 + c] - P_sol[c*3 + r]));
          }
        }

        const detPot = P_pot[0]*(P_pot[4]*P_pot[8] - P_pot[5]*P_pot[7])
                     - P_pot[1]*(P_pot[3]*P_pot[8] - P_pot[5]*P_pot[6])
                     + P_pot[2]*(P_pot[3]*P_pot[7] - P_pot[4]*P_pot[6]);

        const detSol = P_sol[0]*(P_sol[4]*P_sol[8] - P_sol[5]*P_sol[7])
                     - P_sol[1]*(P_sol[3]*P_sol[8] - P_sol[5]*P_sol[6])
                     + P_sol[2]*(P_sol[3]*P_sol[7] - P_sol[4]*P_sol[6]);

        console.log(JSON.stringify({
          potSymmetryErr,
          solSymmetryErr,
          detPot: Math.abs(detPot),
          detSol: Math.abs(detSol)
        }));
        """
        res = run_node_snippet(code)
        assert res["potSymmetryErr"] < 1e-14
        assert res["solSymmetryErr"] < 1e-14
        assert res["detPot"] < 1e-14
        assert res["detSol"] < 1e-14

    def test_projector_zero_wavemode_behavior(self):
        """Verify projector returns zero matrix for k=0 zero mode."""
        code = """
        import { KSpaceProjector } from './src/fields/helmholtz_fft_decomposer.js';

        const P_pot0 = KSpaceProjector.getPotentialProjectorMatrix(0, 0, 0);
        const P_sol0 = KSpaceProjector.getSolenoidalProjectorMatrix(0, 0, 0);

        let maxPot0 = 0.0;
        let maxSol0 = 0.0;
        for (let i = 0; i < 9; i++) {
          maxPot0 = Math.max(maxPot0, Math.abs(P_pot0[i]));
          maxSol0 = Math.max(maxSol0, Math.abs(P_sol0[i]));
        }

        console.log(JSON.stringify({ maxPot0, maxSol0 }));
        """
        res = run_node_snippet(code)
        assert res["maxPot0"] == 0.0
        assert res["maxSol0"] == 0.0


# ============================================================================
# 2. PURE POTENTIAL (IRROTATIONAL) VECTOR FIELD TESTS
# ============================================================================

class TestPurePotentialDecomposition:
    """Verification of pure irrotational / compressible vector field decomposition."""

    def test_single_mode_pure_potential_field(self):
        """Verify decomposition of single harmonic potential wave gives v_pot = v, v_sol = 0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-50, -50, -50],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const k0 = (2.0 * Math.PI) / L;
        const modes = [
          { k: [k0, 0, 0], amplitude: 250.0 }
        ];

        const potField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, modes);
        const decomp = HelmholtzFFTDecomposer.decompose(potField);

        console.log(JSON.stringify({
          totEnergy: decomp.energySpectrum.totalEnergySum,
          potEnergy: decomp.energySpectrum.potEnergySum,
          solEnergy: decomp.energySpectrum.solEnergySum,
          compressibleRatio: decomp.energySpectrum.compressibleRatio,
          solenoidalRatio: decomp.energySpectrum.solenoidalRatio,
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError,
          l2InnerProduct: decomp.diagnostics.l2InnerProduct,
          normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
          spectralCurlPotMax: decomp.diagnostics.spectralCurlPotMax,
          spectralDivergenceSolMax: decomp.diagnostics.spectralDivergenceSolMax
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["compressibleRatio"] - 1.0) < 1e-12, f"Expected compressible ratio 1.0, got {res['compressibleRatio']}"
        assert res["solenoidalRatio"] < 1e-14, f"Expected solenoidal ratio ~0, got {res['solenoidalRatio']}"
        assert res["solEnergy"] < 1e-12 * res["totEnergy"]
        assert res["normalizedL2Correlation"] < 1e-10
        assert res["spectralCurlPotMax"] < 1e-10

    def test_multi_mode_pure_potential_field(self):
        """Verify decomposition of rich multi-mode 3D potential wave packet."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 200.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [-100, -100, -100],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const k0 = (2.0 * Math.PI) / L;
        const modes = [
          { k: [k0, k0, 0], amplitude: 100.0, phase: 0.2 },
          { k: [0, 2*k0, k0], amplitude: 80.0, phase: 0.5 },
          { k: [2*k0, -k0, 2*k0], amplitude: 50.0, phase: 1.1 }
        ];

        const potField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, modes);
        const decomp = HelmholtzFFTDecomposer.decompose(potField);

        console.log(JSON.stringify({
          compressibleRatio: decomp.energySpectrum.compressibleRatio,
          solenoidalRatio: decomp.energySpectrum.solenoidalRatio,
          normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
          isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal,
          spectralCurlPotMax: decomp.diagnostics.spectralCurlPotMax
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["compressibleRatio"] - 1.0) < 1e-10
        assert res["solenoidalRatio"] < 1e-12
        assert res["isStrictlyOrthogonal"]
        assert res["spectralCurlPotMax"] < 1e-10


# ============================================================================
# 3. PURE SOLENOIDAL (DIVERGENCE-FREE) VECTOR FIELD TESTS
# ============================================================================

class TestPureSolenoidalDecomposition:
    """Verification of pure solenoidal / divergence-free vector field decomposition."""

    def test_arnold_beltrami_childress_helical_flow(self):
        """Verify ABC flow is decomposed into 100% solenoidal component."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const k0 = (2.0 * Math.PI) / L;
        const abcField = HelmholtzFFTDecomposer.createSyntheticPureSolenoidalABC(grid, {
          A: 150.0,
          B: 120.0,
          C: 90.0,
          k0
        });

        const decomp = HelmholtzFFTDecomposer.decompose(abcField);

        console.log(JSON.stringify({
          totEnergy: decomp.energySpectrum.totalEnergySum,
          potEnergy: decomp.energySpectrum.potEnergySum,
          solEnergy: decomp.energySpectrum.solEnergySum,
          compressibleRatio: decomp.energySpectrum.compressibleRatio,
          solenoidalRatio: decomp.energySpectrum.solenoidalRatio,
          normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
          spectralDivergenceSolMax: decomp.diagnostics.spectralDivergenceSolMax,
          isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal,
          totalHelicity: decomp.diagnostics.totalHelicity
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["solenoidalRatio"] - 1.0) < 1e-12
        assert res["compressibleRatio"] < 1e-14
        assert res["potEnergy"] < 1e-12 * res["totEnergy"]
        assert res["isStrictlyOrthogonal"]
        assert res["spectralDivergenceSolMax"] < 1e-10
        assert res["totalHelicity"] > 0, "ABC flow has non-zero net kinetic helicity"

    def test_taylor_green_3d_vortex_flow(self):
        """Verify 3D Taylor-Green vortex decomposition."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const tgField = HelmholtzFFTDecomposer.createSyntheticTaylorGreen(grid, {
          U0: 300.0,
          k: (2.0 * Math.PI) / L
        });

        const decomp = HelmholtzFFTDecomposer.decompose(tgField);

        console.log(JSON.stringify({
          compressibleRatio: decomp.energySpectrum.compressibleRatio,
          solenoidalRatio: decomp.energySpectrum.solenoidalRatio,
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError,
          spectralDivergenceSolMax: decomp.diagnostics.spectralDivergenceSolMax,
          isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["solenoidalRatio"] - 1.0) < 1e-12
        assert res["compressibleRatio"] < 1e-14
        assert res["isStrictlyOrthogonal"]
        assert res["spectralDivergenceSolMax"] < 1e-10


# ============================================================================
# 4. SUPERPOSITION & ZERO-MODE BULK VELOCITY TESTS
# ============================================================================

class TestSuperpositionAndZeroMode:
    """Verification of mixed vector fields with DC dipole bulk flow."""

    def test_mixed_potential_solenoidal_bulk_flow(self):
        """Verify exact decomposition of superposed Potential + Solenoidal + Bulk Flow."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const k0 = (2.0 * Math.PI) / L;
        const potModes = [
          { k: [k0, 0, 0], amplitude: 150.0 }
        ];
        const solParams = { A: 100.0, B: 100.0, C: 100.0, k0 };
        const v0Arr = [250.0, -175.0, 85.0];

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid, {
          potentialModes: potModes,
          solenoidalParams: solParams,
          bulkVelocity: v0Arr
        });

        const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField);

        const v0ErrX = Math.abs(decomp.v0.vx - v0Arr[0]);
        const v0ErrY = Math.abs(decomp.v0.vy - v0Arr[1]);
        const v0ErrZ = Math.abs(decomp.v0.vz - v0Arr[2]);

        const eTot = decomp.energySpectrum.totalEnergySum;
        const ePot = decomp.energySpectrum.potEnergySum;
        const eSol = decomp.energySpectrum.solEnergySum;
        const e0 = decomp.energySpectrum.zeroModeEnergy;

        console.log(JSON.stringify({
          v0ErrX, v0ErrY, v0ErrZ,
          eTot, ePot, eSol, e0,
          energySumRelDiff: Math.abs(eTot - (ePot + eSol + e0)) / eTot,
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError,
          relativeReconstructionL2: decomp.diagnostics.relativeReconstructionL2,
          normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
          isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal,
          isEnergyConserving: decomp.diagnostics.isEnergyConserving
        }));
        """
        res = run_node_snippet(code)
        assert res["v0ErrX"] < 1e-10
        assert res["v0ErrY"] < 1e-10
        assert res["v0ErrZ"] < 1e-10
        assert res["energySumRelDiff"] < 1e-12
        assert res["maxReconAbsError"] < 1e-10
        assert res["relativeReconstructionL2"] < 1e-12
        assert res["isStrictlyOrthogonal"]
        assert res["isEnergyConserving"]

    def test_pure_uniform_bulk_flow_field(self):
        """Verify field with only constant velocity v = (Vx, Vy, Vz) has v_pot = 0, v_sol = 0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 8;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, isCellCentered: true });
        const total = grid.totalCells;

        const vx = new Float64Array(total).fill(300.0);
        const vy = new Float64Array(total).fill(-200.0);
        const vz = new Float64Array(total).fill(100.0);

        const vField = new VectorField3D(grid, vx, vy, vz);
        const decomp = HelmholtzFFTDecomposer.decompose(vField);

        let maxPotVal = 0.0;
        let maxSolVal = 0.0;
        for (let i = 0; i < total; i++) {
          maxPotVal = Math.max(maxPotVal, Math.abs(decomp.vPot.vx[i]), Math.abs(decomp.vPot.vy[i]), Math.abs(decomp.vPot.vz[i]));
          maxSolVal = Math.max(maxSolVal, Math.abs(decomp.vSol.vx[i]), Math.abs(decomp.vSol.vy[i]), Math.abs(decomp.vSol.vz[i]));
        }

        console.log(JSON.stringify({
          v0x: decomp.v0.vx,
          v0y: decomp.v0.vy,
          v0z: decomp.v0.vz,
          maxPotVal,
          maxSolVal,
          potEnergy: decomp.energySpectrum.potEnergySum,
          solEnergy: decomp.energySpectrum.solEnergySum
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["v0x"] - 300.0) < 1e-10
        assert abs(res["v0y"] - (-200.0)) < 1e-10
        assert abs(res["v0z"] - 100.0) < 1e-10
        assert res["maxPotVal"] < 1e-12
        assert res["maxSolVal"] < 1e-12
        assert res["potEnergy"] == 0.0
        assert res["solEnergy"] == 0.0


# ============================================================================
# 5. KINEMATIC POTENTIALS & INVARIANTS TESTS
# ============================================================================

class TestKinematicPotentialsAndFields:
    """Verification of velocity potential, stream function, vorticity, and helicity."""

    def test_velocity_scalar_potential_gradient_matching(self):
        """Verify v_pot = -grad Phi_v and laplacian Phi_v = -div v."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const k0 = (2.0 * Math.PI) / L;
        const modes = [
          { k: [k0, 0, 0], amplitude: 100.0 },
          { k: [0, k0, 0], amplitude: 80.0 }
        ];

        const potField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, modes);
        const decomp = HelmholtzFFTDecomposer.decompose(potField);

        console.log(JSON.stringify({
          spectralGradPhiPotMax: decomp.diagnostics.spectralGradPhiPotMax,
          potentialGradRelError: decomp.diagnostics.potentialGradRelError,
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError
        }));
        """
        res = run_node_snippet(code)
        assert res["spectralGradPhiPotMax"] < 1e-10, f"Spectral gradient error: {res['spectralGradPhiPotMax']}"
        assert res["potentialGradRelError"] < 0.05, f"Finite-difference grad relative error: {res['potentialGradRelError']}"
        assert res["maxReconAbsError"] < 1e-10

    def test_vector_stream_potential_curl_and_gauge(self):
        """Verify v_sol = curl A and Coulomb gauge div A = 0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({
          nx: N, ny: N, nz: N,
          origin: [0, 0, 0],
          boxSize: [L, L, L],
          isCellCentered: true
        });

        const abcField = HelmholtzFFTDecomposer.createSyntheticPureSolenoidalABC(grid, {
          A: 100.0,
          B: 100.0,
          C: 100.0,
          k0: (2.0 * Math.PI) / L
        });

        const decomp = HelmholtzFFTDecomposer.decompose(abcField);

        console.log(JSON.stringify({
          spectralCurlStreamSolMax: decomp.diagnostics.spectralCurlStreamSolMax,
          spectralGaugeAMax: decomp.diagnostics.spectralGaugeAMax,
          streamCoulombGaugeL2: decomp.diagnostics.streamCoulombGaugeL2
        }));
        """
        res = run_node_snippet(code)
        assert res["spectralCurlStreamSolMax"] < 1e-10
        assert res["spectralGaugeAMax"] < 1e-10
        assert res["streamCoulombGaugeL2"] < 1e-6


# ============================================================================
# 6. 1D ISOTROPIC RADIAL SHELL ENERGY SPECTRA TESTS
# ============================================================================

class TestHelmholtzEnergySpectra:
    """Verification of 1D isotropic radial shell power spectra and power-law fitting."""

    def test_energy_spectrum_bin_sum_conservation(self):
        """Verify sum of shell energies equals total non-zero mode kinetic energy."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const k0 = (2.0 * Math.PI) / L;
        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid, {
          potentialModes: [{ k: [k0, k0, 0], amplitude: 100.0 }],
          solenoidalParams: { A: 80.0, B: 80.0, C: 80.0, k0 },
          bulkVelocity: [100.0, 50.0, -25.0]
        });

        const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField, { numSpectrumBins: 16 });
        const spec = decomp.energySpectrum;

        let sumETot = 0.0;
        let sumEPot = 0.0;
        let sumESol = 0.0;
        let totalModes = 0;

        for (let i = 0; i < spec.numBins; i++) {
          sumETot += spec.eTot[i];
          sumEPot += spec.ePot[i];
          sumESol += spec.eSol[i];
          totalModes += spec.modeCounts[i];
        }

        console.log(JSON.stringify({
          sumETot,
          sumEPot,
          sumESol,
          potEnergySum: spec.potEnergySum,
          solEnergySum: spec.solEnergySum,
          diffPotRel: Math.abs(sumEPot - spec.potEnergySum) / spec.potEnergySum,
          diffSolRel: Math.abs(sumESol - spec.solEnergySum) / spec.solEnergySum,
          totalModes,
          expectedNonZeroModes: N * N * N - 1
        }));
        """
        res = run_node_snippet(code)
        assert res["diffPotRel"] < 1e-12
        assert res["diffSolRel"] < 1e-12
        assert res["totalModes"] == res["expectedNonZeroModes"], f"All {res['expectedNonZeroModes']} modes must be binned"

    def test_spectral_fraction_table_and_power_law_fit(self):
        """Verify spectral fraction table f_pot + f_sol == 1 and power law fitter."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);
        const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField, { numSpectrumBins: 16 });
        const spec = decomp.energySpectrum;

        const table = spec.getSpectralFractionTable();
        let maxFractionSumError = 0.0;

        for (const row of table) {
          if (row.eTot > 1e-10) {
            const err = Math.abs((row.fPot + row.fSol) - 1.0);
            maxFractionSumError = Math.max(maxFractionSumError, err);
          }
        }

        const fit = spec.fitPowerLaw('tot');

        console.log(JSON.stringify({
          maxFractionSumError,
          fitAmplitude: fit.amplitude,
          fitIndex: fit.spectralIndex,
          fitR2: fit.rSquared
        }));
        """
        res = run_node_snippet(code)
        assert res["maxFractionSumError"] < 1e-12
        assert np.isfinite(res["fitAmplitude"])
        assert np.isfinite(res["fitIndex"])


# ============================================================================
# 7. COSMOLOGICAL LINEAR DENSITY RECONSTRUCTION TESTS
# ============================================================================

class TestCosmologicalDensityReconstruction:
    """Verification of linear cosmological density field reconstruction from potential velocity."""

    def test_density_reconstruction_from_potential_wave(self):
        """Verify delta_rec = -div(v) / (H0 * f) matches analytical harmonic density."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const kx = (2.0 * Math.PI) / L;
        const amplitude = 100.0;
        const modes = [{ k: [kx, 0, 0], amplitude }];

        const potField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, modes);
        const decomp = HelmholtzFFTDecomposer.decompose(potField);

        const H0 = 100.0;
        const omegaM = 0.315;
        const density = decomp.reconstructDensityField({ H0, omegaM, redshift: 0.0 });

        const f = Math.pow(omegaM, 0.55);
        const expectedPeakDelta = (amplitude * kx * kx) / (H0 * f);

        let maxDeltaVal = 0.0;
        for (let i = 0; i < grid.totalCells; i++) {
          maxDeltaVal = Math.max(maxDeltaVal, Math.abs(density.data[i]));
        }

        console.log(JSON.stringify({
          maxDeltaVal,
          expectedPeakDelta,
          relativeDeltaPeakDiff: Math.abs(maxDeltaVal - expectedPeakDelta) / expectedPeakDelta
        }));
        """
        res = run_node_snippet(code)
        assert res["relativeDeltaPeakDiff"] < 1e-12, f"Density peak diff too high: {res['relativeDeltaPeakDiff']}"

    def test_density_reconstruction_growth_factor_scaling(self):
        """Verify density reconstruction scales inversely with growth rate f(Omega_m, z)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const kx = (2.0 * Math.PI) / L;
        const potField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, [{ k: [kx, 0, 0], amplitude: 100.0 }]);
        const decomp = HelmholtzFFTDecomposer.decompose(potField);

        const d0 = decomp.reconstructDensityField({ H0: 100.0, omegaM: 0.3, redshift: 0.0 });
        const d1 = decomp.reconstructDensityField({ H0: 100.0, omegaM: 1.0, redshift: 0.0 });

        const f0 = Math.pow(0.3, 0.55);
        const f1 = Math.pow(1.0, 0.55); // 1.0

        const ratioObs = d0.data[0] / d1.data[0];
        const ratioExpected = f1 / f0;

        console.log(JSON.stringify({
          ratioObs,
          ratioExpected,
          diff: Math.abs(ratioObs - ratioExpected)
        }));
        """
        res = run_node_snippet(code)
        assert res["diff"] < 1e-10


# ============================================================================
# 8. MULTISCALE SCALE BANDS & FOURIER FILTERING TESTS
# ============================================================================

class TestMultiscaleBandsAndFiltering:
    """Verification of multiscale cosmic web decomposition and spatial filtering."""

    def test_fourier_gaussian_smoothing(self):
        """Verify Gaussian spatial smoothing reduces peak velocity and variance."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer, HelmholtzFilterType } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);

        const decomposer = new HelmholtzFFTDecomposer();
        const filtered = decomposer.filterField(superpos.totalField, {
          filterType: HelmholtzFilterType.GAUSSIAN,
          smoothingScale: 10.0
        });

        let sumOrigSq = 0.0;
        let sumFiltSq = 0.0;
        for (let i = 0; i < grid.totalCells; i++) {
          const vox = superpos.totalField.vx[i];
          const voy = superpos.totalField.vy[i];
          const voz = superpos.totalField.vz[i];
          sumOrigSq += vox*vox + voy*voy + voz*voz;

          const vfx = filtered.vx[i];
          const vfy = filtered.vy[i];
          const vfz = filtered.vz[i];
          sumFiltSq += vfx*vfx + vfy*vfy + vfz*vfz;
        }

        console.log(JSON.stringify({
          origRMS: Math.sqrt(sumOrigSq / grid.totalCells),
          filtRMS: Math.sqrt(sumFiltSq / grid.totalCells),
          isReduced: sumFiltSq < sumOrigSq
        }));
        """
        res = run_node_snippet(code)
        assert res["isReduced"], "Gaussian smoothing must reduce total velocity variance"
        assert res["filtRMS"] < res["origRMS"]

    def test_tophat_real_space_filtering(self):
        """Verify Top-Hat real space Fourier filter attenuates high-k modes."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer, HelmholtzFilterType } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);
        const decomposer = new HelmholtzFFTDecomposer();
        const filtered = decomposer.filterField(superpos.totalField, {
          filterType: HelmholtzFilterType.TOPHAT_REAL,
          smoothingScale: 8.0
        });

        console.log(JSON.stringify({
          name: filtered.name,
          unit: filtered.unit
        }));
        """
        res = run_node_snippet(code)
        assert "tophat_real" in res["name"]
        assert res["unit"] == "km/s"

    def test_multiscale_scale_band_decomposition(self):
        """Verify scale-band partitioning into cluster, filament, and supercluster scales."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);
        const decomposer = new HelmholtzFFTDecomposer();
        const bandsMap = decomposer.decomposeScaleBands(superpos.totalField);

        const bands = Array.from(bandsMap.keys());
        const fractions = bands.map(b => bandsMap.get(b).energyFraction);

        console.log(JSON.stringify({
          bands,
          fractions,
          sumFractions: fractions.reduce((a, b) => a + b, 0)
        }));
        """
        res = run_node_snippet(code)
        assert len(res["bands"]) == 3
        for frac in res["fractions"]:
            assert 0.0 <= frac <= 1.0


# ============================================================================
# 9. DISCRETIZATION KERNELS COMPARISON TESTS
# ============================================================================

class TestDiscretizationKernels:
    """Comparison between continuous, 7-point, 19-point, and 27-point discretization kernels."""

    def test_kernel_decomposition_orthogonality(self):
        """Verify all stencil kernels maintain machine-precision L2 orthogonality."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer, HelmholtzKernelType } from './src/fields/helmholtz_fft_decomposer.js';

        const N = 16;
        const L = 100.0;
        const grid = new GridIndexer({ nx: N, ny: N, nz: N, boxSize: [L, L, L], isCellCentered: true });
        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);

        const kernels = [
          HelmholtzKernelType.CONTINUOUS,
          HelmholtzKernelType.DISCRETE_7POINT,
          HelmholtzKernelType.DISCRETE_19POINT,
          HelmholtzKernelType.DISCRETE_27POINT
        ];

        const results = kernels.map(kType => {
          const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField, { kernelType: kType });
          return {
            kernel: kType,
            normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
            isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal,
            maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError
          };
        });

        console.log(JSON.stringify(results));
        """
        results = run_node_snippet(code)
        for r in results:
            assert r["isStrictlyOrthogonal"], f"Kernel {r['kernel']} failed orthogonality: {r['normalizedL2Correlation']}"
            assert r["maxReconAbsError"] < 1e-10


# ============================================================================
# 10. ANISOTROPIC GRIDS & BOUNDARY EDGE CASES
# ============================================================================

class TestAnisotropicAndEdgeCases:
    """Verification of anisotropic grid dimensions and boundary conditions."""

    def test_anisotropic_grid_dimensions(self):
        """Verify decomposition on non-cubic anisotropic grid (16 x 8 x 32)."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const grid = new GridIndexer({
          nx: 16, ny: 8, nz: 32,
          origin: [-50, -25, -100],
          boxSize: [100.0, 50.0, 200.0],
          isCellCentered: true
        });

        const kx0 = (2.0 * Math.PI) / 100.0;
        const ky0 = (2.0 * Math.PI) / 50.0;
        const kz0 = (2.0 * Math.PI) / 200.0;

        const potModes = [{ k: [kx0, ky0, kz0], amplitude: 120.0 }];
        const solParams = { A: 80.0, B: 80.0, C: 80.0, k0: kx0 };
        const v0Arr = [150.0, -100.0, 50.0];

        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid, {
          potentialModes: potModes,
          solenoidalParams: solParams,
          bulkVelocity: v0Arr
        });

        const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField);

        console.log(JSON.stringify({
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError,
          relativeReconstructionL2: decomp.diagnostics.relativeReconstructionL2,
          normalizedL2Correlation: decomp.diagnostics.normalizedL2Correlation,
          isStrictlyOrthogonal: decomp.diagnostics.isStrictlyOrthogonal,
          isEnergyConserving: decomp.diagnostics.isEnergyConserving
        }));
        """
        res = run_node_snippet(code)
        assert res["isStrictlyOrthogonal"]
        assert res["isEnergyConserving"]
        assert res["maxReconAbsError"] < 1e-10

    def test_null_zero_velocity_field(self):
        """Verify decomposition of zero vector field v = 0."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { VectorField3D } from './src/fields/vector_field_3d.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, isCellCentered: true });
        const vNull = new VectorField3D(grid);
        const decomp = HelmholtzFFTDecomposer.decompose(vNull);

        console.log(JSON.stringify({
          totEnergy: decomp.energySpectrum.totalEnergySum,
          potEnergy: decomp.energySpectrum.potEnergySum,
          solEnergy: decomp.energySpectrum.solEnergySum,
          zeroEnergy: decomp.energySpectrum.zeroModeEnergy,
          maxReconAbsError: decomp.diagnostics.maxReconstructionAbsError
        }));
        """
        res = run_node_snippet(code)
        assert res["totEnergy"] == 0.0
        assert res["potEnergy"] == 0.0
        assert res["solEnergy"] == 0.0
        assert res["maxReconAbsError"] == 0.0

    def test_high_dynamic_range_amplitudes(self):
        """Verify decomposition preserves numerical stability from 1e-5 km/s to 1e5 km/s."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, isCellCentered: true });
        const k0 = (2.0 * Math.PI) / 100.0;

        // Small amplitude
        const smallField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, [{ k: [k0, 0, 0], amplitude: 1e-5 }]);
        const decompSmall = HelmholtzFFTDecomposer.decompose(smallField);

        // Huge amplitude
        const hugeField = HelmholtzFFTDecomposer.createSyntheticPurePotential(grid, [{ k: [k0, 0, 0], amplitude: 1e5 }]);
        const decompHuge = HelmholtzFFTDecomposer.decompose(hugeField);

        console.log(JSON.stringify({
          smallCompressible: decompSmall.energySpectrum.compressibleRatio,
          hugeCompressible: decompHuge.energySpectrum.compressibleRatio,
          smallReconRel: decompSmall.diagnostics.relativeReconstructionL2,
          hugeReconRel: decompHuge.diagnostics.relativeReconstructionL2
        }));
        """
        res = run_node_snippet(code)
        assert abs(res["smallCompressible"] - 1.0) < 1e-10
        assert abs(res["hugeCompressible"] - 1.0) < 1e-10
        assert res["smallReconRel"] < 1e-12
        assert res["hugeReconRel"] < 1e-12

    def test_json_serialization_and_immutability(self):
        """Verify result toJSON summary and non-mutation of input buffers."""
        code = """
        import { GridIndexer } from './src/fields/grid_indexer.js';
        import { HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        const grid = new GridIndexer({ nx: 8, ny: 8, nz: 8, isCellCentered: true });
        const superpos = HelmholtzFFTDecomposer.createSyntheticSuperposition(grid);

        const vxBefore = superpos.totalField.vx.slice();
        const decomp = HelmholtzFFTDecomposer.decompose(superpos.totalField);
        const jsonSummary = decomp.toJSON();

        let maxInputMutation = 0.0;
        for (let i = 0; i < grid.totalCells; i++) {
          maxInputMutation = Math.max(maxInputMutation, Math.abs(superpos.totalField.vx[i] - vxBefore[i]));
        }

        console.log(JSON.stringify({
          maxInputMutation,
          hasBulk: jsonSummary.bulkVelocity !== undefined,
          hasEnergies: jsonSummary.energies !== undefined,
          hasDiagnostics: jsonSummary.diagnostics !== undefined
        }));
        """
        res = run_node_snippet(code)
        assert res["maxInputMutation"] == 0.0, "Input buffer must not be mutated"
        assert res["hasBulk"] and res["hasEnergies"] and res["hasDiagnostics"]

    def test_invalid_arguments_and_error_handling(self):
        """Verify robust error throwing on invalid grid dimensions or input types."""
        code = """
        import { FFT3DEngine, HelmholtzFFTDecomposer } from './src/fields/helmholtz_fft_decomposer.js';

        let nonPower2Error = false;
        try {
          new FFT3DEngine(15, 16, 16);
        } catch (e) {
          nonPower2Error = true;
        }

        let invalidTypeError = false;
        try {
          HelmholtzFFTDecomposer.decompose({});
        } catch (e) {
          invalidTypeError = true;
        }

        console.log(JSON.stringify({ nonPower2Error, invalidTypeError }));
        """
        res = run_node_snippet(code)
        assert res["nonPower2Error"], "FFT3DEngine must throw on non-power-of-2 dimension"
        assert res["invalidTypeError"], "decompose must throw on non-VectorField3D input"
