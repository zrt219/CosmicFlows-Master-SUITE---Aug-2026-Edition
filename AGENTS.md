# Rule: Cosmological Field Dimensions and Attributions
- Never add spatial position vectors (Mpc/h) to velocity vectors (km/s) without a physical time step Delta t.
- Enforce the x52.0 scale factor on CF4 displacement fields.
- Use (SGZ, SGY, SGX) -> (SGX, SGY, SGZ) canonical stride transformations.
- Cite official publications (Courtois et al. 2023, Dupuy & Courtois 2023) for all CF4/watershed computations.

# Scientific Invariants & Validation Gates
1. **Dimensional Strictness**:
   - Supergalactic Cartesian spatial coordinates are strictly measured in Mpc/h.
   - Peculiar velocities are strictly measured in km/s.
   - Cross-dimensional operations must throw a runtime TypeError.

2. **Cosmicflows-4 Axis Order & Scaling**:
   - Official IP2I FITS binary grids use (SGZ, SGY, SGX) storage order in row-major layout.
   - Internal canonical indexing must always map to (SGX, SGY, SGZ).
   - Wiener Filter displacement fields Psi require exact x52.0 velocity scale factor enforcement: v = 52.0 * Psi.

3. **High-Order Interpolation & Differential Operators**:
   - Tricubic 64-point local Hermite splines provide C^1 continuous field interpolation and analytical spatial gradients, Hessians, and Jacobians.
   - Linear perturbation theory continuity div(v) ~= -H0 * f * delta is evaluated with code-enforced residual gates.

4. **Table A.1 Official Watershed Taxonomy**:
   - Watershed basins must map to Dupuy & Courtois (2023) Table A.1 official taxonomy:
     1: Laniakea, 2: Apus, 3: Hercules, 4: Lepus, 5: Perseus-Pisces, 6: Shapley, 7: SDSS-1a, 8: SDSS-2a.

5. **W3C PROV-O Lineage**:
   - Every scientific dataset, intermediate field, and derived calculation carries NIST SHA-256 digests and W3C PROV-O JSON-LD lineage graph metadata.
