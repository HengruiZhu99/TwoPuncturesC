# Retained extreme physical investigation

The human override stopped further performance work. This directory preserves
the first chi=.99 and Gamma10 physical investigations without changing any
failed acceptance flag. All seven binary rows remain failed.

`originals/` contains 105 unchanged metadata, log, input and sampler witness
files. `local-verification.json` binds their hashes and repeats the two coarse
binary sampler array comparisons locally: 2454 points each, every field and
metric derivative identical. This qualifies data interchange only.

`retention.json` binds all 151 original physical artifacts, including raw
constraint arrays and seven diagnostic checkpoints, in the full archive at
`/pscratch/sd/h/hzhu/codex-hispid-kokkos-20261002/extreme-human-override-20261003/physical-results-preserved-20261003.tar.gz`.
That archive is 718557306 bytes, SHA256
`f4c3cae2bd30ccc4f830292f1cb67465f0addd5a00e168db8783a2a663cb8ceb`.
The original files also remain on Perlmutter. The full raw archive is not
committed or claimed as a completed local backup; its compact metadata and
sampler witnesses are retained here. All 21 raw bindings from the seven
binary rows were verified before archival, and all original files were
rehashed after archival.

The spin CPU replay under `originals/spin99-coarse-reference-replay` used
the preserved original script. `binding-supplement.json` is explicitly a
post-run audit of its coefficient bits, images, basis, maps and norm, without
a repeated numerical replay or modification of the original result. The
new `../diagnose_reference_iterate.py` performs those checks before creating
a context for future use. No elliptic solve occurs in that diagnostic.

The Gamma import-only record under `originals/gamma10-coarse-import-only`
passes the roundtrip and time/cycle-zero checks with no finder constructed.
Spin component and Gamma seed horizon failures remain explicit. Exact Kerr
seed horizon controls are separate from binary acceptance.

The report compilation receipts refer to the existing open document. The
earlier successful source snapshot is preserved alongside the final compiled
snapshot. Allocation hold/resume/release receipts confirm that only the task
controller was briefly held and allocation59284155 was released after all
authorized science and retention steps finished.
