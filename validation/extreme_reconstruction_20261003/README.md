# Preserved Cartesian reconstruction diagnosis

No new performance campaign or elliptic solve was run for this record. Existing
performance tables and all seven failed extreme binary rows remain unchanged.
The artifacts here retain the completed local ARM inspection of the coarse
aligned-spin chi=0.99 iterate; they do not establish physical acceptance.

At three retained off-grid points, direct Cartesian polynomial Hessians give
Hamiltonian residuals -41.482618, -43.992324 and -164.122065. The original finite-difference values are -41.482592, -43.992332 and -164.123124. At the three
selected collocation nodes, Cartesian Hamiltonian residuals are between
3.7e-11 and 2.9e-9 in magnitude. This supports a between-node defect at these
points. It does not establish its sole cause or validate the full binary.

The unchanged six-point sampler comparison **fails** its 1e-12 criterion:
maximum scaled difference 5.02e-10 in metric derivatives. ARM long double has
53 mantissa digits. Native seed/tensor algebra is reused, and momentum outputs
are local-frame components; the current checkpoint frame is identity.
The executed wrapper checked reconstructed-field finiteness only; the separate
post-run retention receipt confirms all retained comparison arrays are finite.
GSL/CBLAS dependency paths were recorded without run-time byte hashes. No
complete run-time image qualification is claimed.

`cartesian-coarse-spin/result.json`, `Cartesian.json`, `points.txt`, `run.log`
and `build-receipt.json` are unchanged original artifacts.
`retention-arm-v1.json` separately maps their exact source, image, input and
output hashes to copies in `retained-arm-v1`. This is a post-run retention
receipt, not retroactive preflight qualification. Source snapshots are tracked;
large input/image copies remain local and ignored. The original coarse
checkpoint/solve arrays are also retained by the prior full remote physics
archive described in `../extreme_physics_20261003/retention.json`.

The previously started Cartesian Hessian control completed successfully:
normalized errors at N12/N24/N40 are 2.05e-18, 7.64e-17 and 1.10e-15.
Its log and source/image post-run hashes are retained. This is a derivative
control, with no performance or extreme-binary acceptance transfer.
