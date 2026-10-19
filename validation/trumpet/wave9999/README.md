# Near-threshold ingoing packet on Liu Kerr, seed chi=0.9999

Fixed single-hole vacuum CTT free data, mass=1 and spin=(0,0,.9999).
The native rotating electric l=m=2 packet uses amplitude=.001, cycles=4,
r0=12, stationary cutoff2--4, propagation=+1 (ingoing), handedness=+1
(prograde). The seed threshold m Omega_H=.9859568136154433 and carrier
omega=.9760972454792889, 1% below it. A finite four-cycle packet has a
bandwidth; its entire spectrum is not restricted below threshold. This is
an initial-data test, with no evolution or amplification claim.

The inactive mass-zero focus is at z=-32; this changes the coordinate map,
not the physical free data. All levels start from zero corrections, keep
one formulation, tolerance=1e-12, modal-preconditioned GMRES, and use the
existing qualified CUDA image. Host long-double geometry and independent
physical sampling use16 threads. No backend benchmark sweep is repeated.

## Reproduction

Use the production C1 focus and polar-tau build with row power3,
compensated modal projection and batched modal transfer. The separate CPU
reader must use identical continuous-basis flags. Set OMP_NUM_THREADS=16,
OPENBLAS_NUM_THREADS=1, OMP_PROC_BIND=spread, OMP_PLACES=cores and
PYTHONPATH=python:examples:validation. Run the CUDA command on one shared
Perlmutter GPU with32 allocated CPUs.

```sh
python validation/test_wave22_gauge.py --library /absolute/cpu/libHiSpID.so
python validation/wave22_study.py --library /absolute/cuda/libHiSpID.so --output /new/study --levels 96 128 160 --focus-separation 32 --chi .9999 --frequency-ratio .99 --amplitude .001 --cycles 4 --r0 12 --multipole-radius 24.874096994494394 --extended-observers
python validation/wave22_final_check.py --library /absolute/cuda/libHiSpID.so --checkpoint /new/study/n160.hispid --study /new/study/result.json --output /new/final-check.json
python validation/portable_sampler_migration.py --checkpoint /new/study/n160.hispid --producer-library /absolute/cuda/libHiSpID.so --consumer-library /absolute/cpu/libHiSpID.so --output /new/migration.json
python /absolute/AthenaK/tst/test_suite/z4c/check_hispid_binary.py --single --executable /absolute/AthenaK/athena --checkpoint /new/study/n160.hispid --migration-proof /new/migration.json --output /new/horizon --levels 8,12,16 --geometry-threads 16 --strict-expansion --harmonic-storage factorized --guess-scale 1.00001 --domain-half-width 48 --reuse-shapes
```

The fixed54 exterior observers include radii1.1 and1.5 times the seed horizon
radius and sample the packet to its outer tail. Cartesian FD6 constraints
are raw unweighted point RMS, not the native scaled residual. Resolution
means N_r=96,128,160 and N_polar=2N_r; N_phi=16 stays fixed. Decreasing errors
here establish convergence at these observers under radial/polar refinement,
not a global norm or an independent azimuthal convergence study.

Metric multipoles at the packet midpoint use the actual seed spin in the
stationary lapse/shift. After subtracting the unperturbed Kerr multipole,
H=r delta(gamma22), the ingoing/outgoing coordinate diagnostics are
(dtH+drH)/2 and (dtH-drH)/2. They are asymptotic metric characteristics,
not Weyl energy fluxes. Coordinate horizon spin in a nonaxisymmetric
geometry is not an approximate-Killing-vector spin.

Prior moderate-spin results and incompatible-start failures are preserved
separately under ../wave22/. Large checkpoints remain on Perlmutter; the
receipt binds their hashes and the producer/reader/executable identities.

## Completed producer checks

All cold starts converge in3 Newton steps, with102/101/107 Krylov iterations.
Hamiltonian RMS:5.04918e-8,2.79218e-8,1.64984e-8.
Momentum RMS:7.28171e-7,4.07185e-7,2.60974e-7.
Observed radial/polar orders are about2.0--2.4, not exponential.
Successive sampled metric and K differences decrease. Scaling the finest
compatible correction by1.05 reconverges to the same sampled physical fields
within8.70e-12 (metric) and1.34e-12 (K). Auxiliary nodal-P differences are
2.92e-5; all nodal coefficients are not claimed identical. Halving FD6 steps
changes H RMS by0.57%, with maximum H difference5.19e-9.
The separate producer/pure-CPU reader comparison is bitwise at1602 points
for all9 fields, including metric derivatives.

The corrected high-spin gauge regression has maximum mode error1.6e-17;
the old chi=.9 diagnostic assumption would cause8.65e-8 dt-mode error.
The final coordinate-characteristic outgoing/incoming amplitude ratio is
.00269037 at the packet midpoint. It remains predominantly ingoing under
this metric diagnostic; no curvature flux or evolution was computed.

The AthenaK schedule uses lmax8/12/16 and a final quadrature increase32->48.
Each search has unchanged mass-stabilization tolerance1e-12 and expansion
RMS tolerance1e-7, with an interior-ball enclosure check. The first initial
radius is1.00001 times the Kerr seed horizon radius; later searches use
previous passed coefficients only as initial guesses. Near-extremality
slows the first FastFlow mean-radius mode. The finder samples the imported
native tensors directly; its mesh import/ADM-Z4c checks run at initial time
with zero evolution steps.

## Completed AthenaK checks

All lmax8/12/16 searches and the final ntheta48 quadrature run pass.
Final A=25.48816340684645, Mirr=0.7120891092880841,
MH=1.000000000118978, chi_coord=0.9998999997620237,
expansion RMS=9.817196e-08. The enclosure margin after the
empirical refinement buffer is.05353573249. Area changes by<5.2e-15
relative and spin vector by<1.0e-15 in the last angular checks; these do not
estimate the physical continuum error of the initial data. The first
near-extreme flow takes686s and375 updates; later shapes warm-start the
remaining searches. No tolerances were relaxed.

Perlmutter shared GPU allocation59973854 (one GPU,32 CPUs) was canceled
after completion and confirmed absent from squeue. Large checkpoints and
raw sampler artifacts are retained under
/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005/.
