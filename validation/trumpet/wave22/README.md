# Kerr metric/time-derivative free data

One fixed stationary Liu Kerr (M=1, seed chi=.9) plus the rotating electric
l=m=2 SpEC packet: amplitude=.001, omega=.5, cycles=1, r0=6, stationary
cutoff2--4. The active hole is at0 and the inactive focus is at z=-12.
The reference Kerr lapse/shift and conformal factor specify the ADM/CTT split.
Callbacks supply physical delta-gamma and fixed-coordinate delta-gamma-dot
jets. They must vanish inside/near the throat and obey asymptotically flat
falloff. The default zero correction is the qualified initial guess.

## Reproduce

Use the existing production C1 focus/polar-tau build, row power3, compensated
modal projection and batched modal transfers. Build one CUDA producer and a
pure reference CPU reader with identical continuous basis flags. Geometry
setup and reading use host long double; the coupled solve uses one A100.
Set OMP_NUM_THREADS=16 and OPENBLAS_NUM_THREADS=1 before running.

```sh
export PYTHONPATH=python:examples:validation
python validation/test_free_data.py --library /absolute/path/libHiSpID.so
python validation/wave22_study.py --library /absolute/path/cuda/libHiSpID.so --output /new/study --focus-separation 12 --levels 48 64 96
python validation/wave22_final_check.py --library /absolute/path/cuda/libHiSpID.so --checkpoint /new/study/n96.hispid --study /new/study/result.json --output /new/final-check.json
python validation/portable_sampler_migration.py --checkpoint /new/study/n96.hispid --producer-library /absolute/path/cuda/libHiSpID.so --consumer-library /absolute/path/cpu/libHiSpID.so --output /new/migration.json
```

In the AthenaK fork, build z4c/hispid with the new HiSpID header and the CPU
reader library, then use:

```sh
python tst/test_suite/z4c/check_hispid_binary.py --single --executable /absolute/path/athena --checkpoint /new/study/n96.hispid --migration-proof /new/migration.json --output /new/horizons --levels 8,12,16 --geometry-threads 16 --strict-expansion --harmonic-storage factorized --guess-scale 1.00001 --domain-half-width 24 --reuse-shapes
```

The exact successful commands are retained in horizon-evidence.

## Evidence and limits

map12-study.json is the unchanged three-resolution campaign. Finest raw
physical H/M RMS are4.004e-8/5.224e-7 at30 fixed exterior points. final-check
recovers the same sampled physical fields from a compatible nonzero start
to7.6e-14 and checks the FD6 step. migration proves bitwise reading at1602
points; horizons and horizon-evidence retain the actual AthenaK finder
searches, input bindings, area/spin measures and angular refinement checks.

The data remain preliminary: point norms are not global bounds. The physical
gdot is changed by constraint projection. Neither superradiant evolution nor
perturbed chi=.9999 data are validated here. incompatible-start.json retains
a failed constant-P start which violates homogeneous focus rows.

coarse.json and map2-study.json preserve the smaller-coordinate-scale
underresolution study. The first n160 attempt hit its16000MiB allocation
guard; the budget was explicitly raised to32768MiB. No formulation or
packet parameters changed. The later grid-scale adjustment changes only
the inactive focus, not the physical free data.

Native images, raw unknowns/checkpoints and sampler witnesses remain under
/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005 on Perlmutter. Source/native
image hashes and the exact receipt bind the selected evidence. No native
binary or large checkpoint is committed. The shared one-GPU allocation was
cancelled after the checks.
