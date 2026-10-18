# Compactified-endpoint and interior regularity study

This directory retains successful controls and failed numerical candidates.
A small weighted collocation residual is not physical binary acceptance.
Current outcome: first-order focus compatibility improves finest bulk momentum
21.9 times and passes all absolute bounds, but the strict three-grid bulk trend
still fails. Second-order compatibility is rejected at160/192. It remains an
explicit experimental option; the first-order implementation is retained.
The fixed three-grid formulation is QI-wide (not the R0=M trumpet):
160x320x16, 192x384x16, 224x448x16, generic moderate spins/boosts.
All comparisons use identical physical free data/windows within that sequence.

- `plan.json`, `stable-plan.json`: stable background arithmetic; original-map
  prefix fails bulk trend. Its finer setup was stopped.
- `quadratic-failed/`: squared infinity distance fails its first linear solve.
- `dual-plan.json`, `dual-endpoint/`: positive-density dual endpoint map;
  160-grid physical bulk momentum fails the absolute bound. No refinement.
- `focus-plan.json`, `focus-tau/`: necessary first Cartesian compatibility at
  both flat-interior foci. All three absolute bounds pass; strict bulk trend
  fails. Finest bulk momentum improves21.9 times over the original sequence.
- `fd6-plan.json`, `focus-fd6/`: predeclared higher-accuracy independent
  observations of those *saved* fields, plus matched BY. Original FD4 verdict
  remains failed. Exact-data and manufactured controls, and exact FD4 replays,
  bind observer accuracy; there are no additional elliptic solves here.
- `focus-c2-compatibility.json`: read-only remaining second-order defects of
  the first-compatibility sequence. `focus2-plan.json` declares the next
  fixed-formulation candidate before solving.
- `focus2-plan.json`, `focus2-coarse/`, `focus2-refinement/`:
  second-order Cartesian compatibility, unchanged free data and compactification.
  The rejected row-scaling preflight and corrected controls remain separate.
  Read each actual assessment; these names do not assert acceptance.
- `by-plan.json`, `by-matched/`: matched-count vanilla-equation BY physical H
  controls. Native stopping gates fail at roundoff; retain diagnostic status.
  Counts, input charges and observation points match, but free metrics and
  horizon properties differ. This is not a performance comparison.
- `by-capacity/`: explicit checked allocation-budget extension, with all33
  baseline arrays/derivative exports bitwise preserved on a supported case.
- `liu-seed-*.json`, `liu-horizon/`: independent tensor pullbacks/vacuum
  controls and actual AthenaK isolated chi=.99 horizon. LES is a wormhole
  coordinate family with finite formal extremal coordinate radius, not a
  trumpet or exact-extremal qualification. `focus-tau/liu160/` is a separate
  moderate binary pilot, not a three-grid validation.
- `*reader*.json`: actual standalone C++ reader and sampler, not evolution.
- `*_perlmutter.sh`: exact stage commands. Build scripts require a compute
  allocation. GPU physics uses one shared A100; host geometry and observers
  use16 CPU threads. There is no broad backend benchmark matrix.

Raw checkpoints and unknown arrays are retained on Perlmutter under
`/pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005`, using the directory names
in these scripts. Large solve archives are not placed in Git. Result JSON
records their SHA256, fixed configuration, native-image identity and source
hashes. Small physical observer NPZs and failure records remain here.
No gate is relaxed by an exact-seed, import, algebra, or horizon control.

The first/second conditions share code through `HISPID_FOCUS_ORDER=1|2`.
Order1 is the default when focus tau is explicitly enabled. Current order-one
weights and retained boundary states are checked in `final-order-control/`.
The stronger condition's isolated algebra controls pass, but its actual
160/192 bulk momentum errors1.05e-6/1.78e-6 fail physical acceptance.
