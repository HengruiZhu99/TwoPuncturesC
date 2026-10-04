# Preserved completed physical trials

Further performance work remains stopped under the existing HUMAN OVERRIDE.
This record preserves work already completed in allocation59293910. No new
solve, finder search, performance check or acceptance upgrade was performed
during preservation. The allocation was released after its workers finished.

The separate HS99UU recipe uses the current modal basis, f width0.2 and no
inner/far attenuation. Only its 80×160×16 attempt ran. The elliptic solve
retains20 Newton/1195 Krylov steps, weighted maximum5.43978e-14 and a Newton
line-search failure. Independent physical sampling then rejected a nonpositive
solved conformal factor. The original incomplete_attempt, solve coefficients,
exception log and exit1 are retained. No completed physical row, finer grid,
charge measurement, import or horizon result is synthesized. Historical
HS99UU evidence has no acceptance transfer to this changed recipe/basis.

The Γ10 trial uses the unchanged128×256×8 physical free data and old failed
checkpoint as a diagnostic initial guess. Restart80→200 and maximum
Krylov2400→4800 are the declared solver changes. The fresh solve reaches the
unchanged1e-14 internal tolerance in2 additional Newton/574 Krylov steps,
with weighted maximum7.91335e-15. Independent near-hole H/M RMS remain
0.00128942/0.217375; bulk H/M RMS are2.14954e-6/0.00121164. Physical gates
remain failed. This retry from an existing iterate is not a standalone timing
comparison or evidence that restart alone explains the difference.

The fresh Γ10 diagnostic checkpoint passes separate-process producer/pure CPU
sampling, including metric derivatives, with zero measured differences.
AthenaK imports it at time/cycle zero with ADM/Z4c error4.59847e-16 and no
finder constructed. These checks qualify data interchange only. No measured
Γ10 component horizon, binary irreducible mass or separation calibration is
claimed.

`inventory.json` binds172 original output and workflow files by size and SHA256.
`archive-receipt.json` binds the complete49,254,760-byte archive (SHA256
`b823ac1827517fa7e0b837d4dfec04b1d9cc66f90c08638c8b5cd1404a80a946`).
The remote archive and original files remain under
`/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/physical-trials-preserved-20261003-v1`.
Local copies are verified separately in `local-verification.json`; large raw
arrays, checkpoints and archives remain ignored rather than entering Git.
The source preparation receipts and run scripts are also retained in
`../extreme_hs99uu_20261003` and `../extreme_gamma10_gmres200_20261003`.
The Γ10 initial checkpoint SHA is
`b5a6f50deff32fb97c44507e6f11586709ed138bfd7419201488677cb2179313`;
its original bytes remain in the earlier full physical archive described by
`../extreme_physics_20261003/retention.json`.
