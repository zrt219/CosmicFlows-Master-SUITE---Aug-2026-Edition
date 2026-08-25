# Rule: Cosmological Field Dimensions, CF4 Semantics, and Scientific Provenance

## 1. Physical Dimensions Are Strict
* Supergalactic Cartesian spatial positions are represented in **Mpc/h** when operating on the CosmicFlows reconstruction grids.
* Peculiar velocities are represented in **km/s**.
* Density contrast (\delta) is dimensionless.
* Never add, subtract, or otherwise combine quantities with incompatible physical dimensions without an explicit physically valid transformation.
* In particular, never perform:
  `position_Mpc_h + velocity_km_s`
  directly.
* Any propagation of position using velocity must define the physical/comoving coordinate convention and an explicit time interval or integration parameter with the required unit conversion.
* Dimensionally invalid operations must fail loudly rather than be silently coerced.

## 2. CosmicFlows-4 Grid Ordering
For the official public CF4 density and velocity grids distributed by the CosmicFlows/IP2I project:
* **Raw file ordering:** `(SGZ, SGY, SGX)`
* **Canonical ZRT representation:** `(SGX, SGY, SGZ)`
All indexing, interpolation, derivative, visualization, streamline, watershed, and topology operations must pass through one tested canonical axis-mapping layer.
No subsystem may independently guess or redefine FITS axis ordering.

## 3. CF4 Velocity Scale Factor
The official CosmicFlows project specifies that the values in the public CF4 velocity grids and their associated velocity-error products must be multiplied by:
$$52.0$$
before being interpreted as physical peculiar velocities.
Apply this scaling:
* only to validated CF4 velocity values;
* to the associated velocity errors where applicable;
* exactly once.
Do not apply ×52 to:
* SGX/SGY/SGZ coordinates;
* density contrast;
* density uncertainties;
* watershed labels;
* arbitrary displacement/potential fields merely because they use related reconstruction mathematics.
Every loaded field must record whether the official ×52 transformation has already been applied. Attempting to apply it twice must fail.

## 4. Linear-Theory Continuity Diagnostic
Where linear perturbation theory is appropriate, evaluate consistency using:
$$\nabla\cdot\mathbf v \approx -aHf\delta$$
At the present epoch this reduces to:
$$\nabla\cdot\mathbf v \approx -H_0f\delta$$
when $a=1$ and $H=H_0$.
This relation is a **linear-regime scientific validation diagnostic**, not an exact identity in strongly nonlinear regions.
Residual gates must therefore record:
* redshift/scale factor;
* $H$;
* growth rate $f$;
* spatial region;
* smoothing/resolution;
* tolerance;
* whether the tested region is expected to satisfy linear theory.

## 5. Watershed IDs Are Reconstruction-Specific
Never use one universal hard-coded basin-ID table for every CosmicFlows watershed file.
The official watershed files follow the ordering of **Dupuy & Courtois (2023), Table A.1**, ignoring entries that are blank for the corresponding reconstruction.
For the **ungrouped CF4 galaxy basin-of-attraction reconstruction**, the verified sequence begins:
1. Laniakea
2. Apus
3. Hercules
4. Lepus
5. Perseus-Pisces
6. Shapley
7. SDSS-1a
8. SDSS-2a
9. SDSS-2b
The grouped reconstruction has a different nonblank sequence.
Therefore every watershed lookup table must be keyed by at least:
`dataset + grouped/ungrouped + attraction/repulsion + publication/version`
The mapping must originate from the authoritative publication or validated dataset metadata.
Never infer basin identity from nearest famous-object coordinates.

## 6. Interpolation Claims Must Be Implementation-Validated
Tricubic/Hermite interpolation may be used when implemented and independently validated, but its properties are **ZRT numerical implementation properties**, not intrinsic CosmicFlows dataset invariants.
If a 64-point tricubic Hermite implementation is used, tests must demonstrate its claimed:
* continuity class;
* interpolation accuracy;
* analytical gradient accuracy;
* Hessian/Jacobian accuracy;
* boundary behavior;
* numerical stability.
Do not state $C^1$ continuity or derivative accuracy merely because the method is labelled “tricubic.”

## 7. Scientific Provenance
Every canonical dataset, transformed dataset, intermediate field, numerical run, derived result, and publication export must carry provenance.
Use:
* **SHA-256** digests for file/result integrity;
* W3C **PROV-O** concepts for lineage;
* JSON-LD as an interchange serialization where appropriate.
At minimum preserve:
`Entity → Activity → Derived Entity`
relationships, including: source dataset, source SHA-256, software version, Git commit, algorithm version, parameters, units, coordinate frame, transformations, generating activity, timestamps, and associated citations.
If the 2024 PROV-JSONLD serialization specification is followed specifically, identify it accurately as the W3C PROV-JSONLD Member Submission rather than claiming Recommendation status.

## 8. Required Attribution
All computations using official CF4 density/velocity products must cite the corresponding CosmicFlows publication:
- **Courtois et al. (2023), A&A 670, L15** (DOI: `10.1051/0004-6361/202245331`)
All computations using the published CF4 watershed products must additionally cite:
- **Dupuy & Courtois (2023), A&A 678, A176** (DOI: `10.1051/0004-6361/202346802`)
Derived ZRT outputs must distinguish clearly between:
* original CosmicFlows data;
* published CosmicFlows watershed products;
* ZRT-recomputed quantities;
* ZRT experimental analyses.

## 9. Failure Principle
When dataset semantics, dimensions, axis ordering, scale factors, basin mapping, provenance, or attribution cannot be verified:
**Do not guess.**
Set the affected computation to:
`BLOCKED — SCIENTIFIC SEMANTICS NOT VERIFIED`
until the ambiguity is resolved.
