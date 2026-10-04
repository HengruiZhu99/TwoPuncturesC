# Fixed failed HS99UU iterate: scalar positivity

The retained current-basis80×160×16 HS99UU attempt failed its Newton line
search and then physical sampling rejected a nonpositive scalar. These
records reconstruct that same saved coefficient state in pure Python; no
native image is loaded and no new solve is performed.

`validation/inspect_scalar_positivity.py` uses barycentric interpolation of
mapped modal P and the analytic unboosted quasi-isotropic Kerr seed scalar.
It guards the supported zero-boost, symmetric-center, choice0/no-far-window
configuration. Two manufactured controls cover modal normalization across
the maps and the Schwarzschild quasi-isotropic limit.

Among the66 original near/bulk samples and their61 offsets (4026 points),
122 scalar values are nonpositive. The minimum is-2.0979918717274204 at a
bulk sample offset; the central sample minimum is-2.017660586629952. This
confirms that positivity failure occurs in the bulk, beyond a near-throat
FD issue. It does not establish the sole cause, nodal positivity, a new
solution, or physical acceptance.

`first-executed-source.py`, `first-result.json` and `first-run-note.json`
preserve the initial result and NumPy small-dot-product warnings unchanged.
The observer then replaced that dot product with an explicit Cartesian sum
and required finite values. The separately retained `verified-result.json`
was obtained with Python warnings promoted to errors, without warnings,
and gives the same numerical extrema/counts. Both versions are retained;
the original record has not been relabeled. Input/state/source hashes are
bound by the results. `preservation.json` is a post-run byte audit only.
