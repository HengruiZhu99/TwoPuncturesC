# Physical radial-resolution qualification, chi=0.9999

This supplements `../wave9999/`; all previously published evidence is retained.
The original finite-observer and horizon checks did not resolve the packet
waveform. Nine physical radial rays, including probes on either side of the
inactive focus, reveal 29.5% carrier value and42.9% radial-derivative errors on
the old160x320x16 grid. Physical radius varies both prolate coordinates.

The physical free data are unchanged. A new purely numerical map uses inactive
mass-zero focus z=-16 and radial stretch1. The carrier audit on384/480/576
radial points,128 polar points and16 azimuthal points passes only at576:
value5.32e-5, derivative5.37e-4, packet-core phase1.12e-4rad. All are below1e-3.
Normalization is to the carrier peak, not to the Kerr background or pointwise
zeros; phase is checked above10% peak. These are finite-ray diagnostics, not
a global tensor/constraint bound. No extra backend performance matrix is run.

The dense solved-wave test is independent: Kerr-subtracted metric m=2 Fourier
rings at257 physical radii and three polar angles,16 azimuths. It compares
amplitude, phase, radial/time derivatives and approximate ingoing/outgoing
metric characteristics across resolutions, normalized to the wave. The
analytic free-wave oracle checks that diagnostic itself. Diagnostic checkpoints
from nonconverged native solves cannot qualify. The outer native tolerance
remains1e-12; inner inexact-Newton tolerances are separate and recorded.

Reproduce the carrier audit from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python validation/wave22_physical_radial_audit.py \
  --study validation/trumpet/wave9999/study.json \
  --output /absolute/path/to/physical-radial-audit.json
```

The matching library basis is `modal_P_C2tauC4_focus_map_v9_r1_k2`.
An old reader with `_r0.2_k2` is incompatible; sampler migration must bind
explicitly to the new producer/reader images. Horizon measurements on the
old checkpoint do not transfer to the new one.

An exploratory second derivative of the *unconstrained scalar interpolant*
fails the0.1% test, even with672/768 radial nodes. Its raw data
are retained separately. This interpolant has no Cartesian focus compatibility,
unlike the native field basis; native free-data derivatives are analytic.
Consequently this result is not a bound on native tensor second derivatives,
and the value/first-derivative carrier pass must not be described as a full
elliptic-source or physical-constraint resolution bound. Dense physical
constraints and solved-wave differences remain required.

## Current acceptance

The revised native sequence fails on its first384x128x16 grid. The final
production attempt uses restart200, maximum4000 inner iterations, fixed
inner relative tolerance.5 and the unchanged1e-12 nonlinear bound. It ends
at9.57e-11 after4516 total Krylov iterations; the last true linear relative
residual is.642. Physical M RMS worsens from3.90e-9 to2.12e-8 at the54
observers, even while the native residual decreases. The output remains
**diagnostic**. The480/576 grids have carrier-audit results only. No revised
finest horizon validation is claimed.

`run-study.sh` records the actual resumed production attempt. Its initial
checkpoint SHA is recorded in the JSON. `run-initial-map-control.sh` and
`run-inner-continuation.sh` record its earlier continuation chain; the native
C++ sources are unchanged from the pinned baseline. To start a new cold
attempt, use fresh output paths and remove `--initial-checkpoint`; this
changes the trajectory, and must be reported as such.

`run-diagnostics.sh` records residual localization plus dense profiles of the
**original saved solutions**, without re-solving those cases. It is diagnostic,
not acceptance of a revised solution. `planned-finest-checks.sh` is a
not-executed template requiring all three revised native solves to pass.
Its profiles, migration and AthenaK checks must not be listed as completed.
The current profile tool tests both m=+/-2 modes and requires decreasing
three-grid differences as well as1e-3 finest-pair wave-relative errors.
Its finest FD6 packet check includes97 radii on five rays (both poles,
mu=+/- .95,0), raw maxima below1e-8 and RMS reduction versus the free trial.

## Completed original-data results

`original-profiles.json/npz` contain three resolutions of the original saved
solutions, on the original r0.2 map. **They fail.** The finest-pair worst
wave-normalized changes are metric1.8718e-6, time derivative.1047114 and
ingoing characteristic.0523327. All four profile differences decrease, but
time derivative/ingoing errors remain far above1e-3. The packet FD6 maximum
H=6.36e-8 and M=6.86e-7 also fail1e-8. Source metadata preserve each
original preliminary checkpoint; this does not promote its acceptance.
The proposed576x128 grid passes only the carrier proxy, not these solved-data
checks. The figure explicitly distinguishes those two data sets.

The completed source used for these profiles and both analytic oracles has
SHA7f94dd8cdc2047510e5dc64869c745f64d082a968ff533b979386749e9694690.
`profile-driver-used.patch` reconstructs that source from the committed
current driver; the latter additionally checks conformal_choice,
inner_flatten and attenuation_power in its formulation inventory. These
fields are confirmed identical in `original-profile-headers.json`; no
numerical sampling or acceptance thresholds changed. The patch's resulting
SHA is verified. `profile-oracle-current.json` tests the current source.

`residual-modes.json` and `residual-mode-localization.json` localize the
revised failed solve:99.56% of squared **scaled native** residual lies in
m0--3; maxima are near the inner transition and horizon. This is not a
physical-constraint bound or proof of the cause. No new near-extreme
AthenaK horizon is claimed. Allocation59975943 was released after diagnostics.
