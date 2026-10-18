# Moderate-spin convergence candidate (18 October 2026)

This directory retains one interior-extension candidate and its validation.
It does not change the default quasi-isotropic or trumpet examples.

The masses, spin vectors, boosts, separation, conformal decomposition, map and
exterior equations are those of `trumpet_moderate`. Only the outer radius of the
interior transition changes from 0.5 to 0.85 of the screened isolated radius;
the flat plateau remains at 0.15. Acceptance requires fresh horizon enclosure.

The legacy `.6 mass` shell overlaps this enlarged modified region. Its original
results remain in `moderate-wide-pilot-v2` on Perlmutter. The declared exterior
observer uses `max(legacy radius, 1.1 inner_max)` for the two fixed-mass near
shells. The 1.5/3 isolated-horizon-radius shells, all bulk coordinates, derivative
stencil, normalization and acceptance thresholds remain unchanged. All three
levels use exactly this profile. Near norms are not the same spatial norm as the
old narrow-transition campaign. Modified-region diagnostics remain separate.

`plan.json` pins the native images, grid sequence, solver options and gates.
`exterior-control.json` checks equation/JVP identity outside the new support.
`coarse160.json` records the coarse result and immutable original checkpoint.

Production steps, each inside an allocated compute node:

1. `pilot_perlmutter.sh`: cold 160x320x16 solve, with the legacy observer retained
   as historical evidence. The driver now uses the correct declared exterior
   profile for new runs.
2. `reobserve.py`: sampling-only observer of the saved coarse checkpoint, with
   the correct exterior profile. It does not solve again or rewrite the data.
3. `sequence_perlmutter.sh`: continue the same image/free data at192x384x16 and
   224x448x16; reuse the independent strict three-grid assessor.
4. `fine_perlmutter.sh`: continue to240x480x16 and assess192/224/240,
   directly testing the previously stalled resolution.
5. `horizons_perlmutter.sh`: requires a passing assessment; validates portable
   sampling and AthenaK initial-time horizons at orders8/12/16 plus quadrature.

The scripts use the existing isolated Perlmutter root under
`/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005`. The workflow bundle must
include `python/`, `examples/`, all top-level `validation/*.py` dependencies and
this directory. No new native library is needed for this candidate.

No broad backend matrix is repeated. Multigrid is deferred: earlier stalled
residuals were grid-scale, and the existing preconditioner solves its global
finite-difference approximation. The current continuation must establish a new
bottleneck before an additional coarse correction is warranted.
