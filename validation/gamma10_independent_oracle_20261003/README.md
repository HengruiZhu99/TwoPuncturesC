# Retained independent Gamma10 seed diagnosis

These completed records are supplemental diagnostics. Original seed/binary
failure flags and all stopping/acceptance tolerances remain unchanged.

`selected-exterior-comparison.json` compares eight selected exterior points
from the existing1512-point producer witness with the independent Decimal80
boosted-Schwarzschild oracle in `validation/boosted_schwarzschild_oracle.py`.
The oracle derives the boosted isotropic four-metric and slice independently
of the native Kerr/graph jets. Analytic Cartesian metric derivatives and a
separate fourth-order Decimal derivative of K supply its constraint check.
The exact binary64 mass, velocity, center and coordinates are retained.
All eight points meet the unchanged scaled1e-12 bound: maximum metric,
extrinsic-curvature and metric-gradient differences are5.54e-15,3.33e-15
and1.37e-14. No native library was loaded or rerun for this comparison.
This is not a complete native-seed check or physical acceptance.

`producer.npz` is kept locally and ignored by Git; `witness-retention.json`
binds its original bytes together with the unchanged checkpoint/migration
proof. It is also part of the earlier retained remote exact-seed evidence.
The comparison binds the executed oracle and test sources; eight controls
were reported passing by the independent reviewer before preservation.

`fd-calibration.json` retains the separate single-point experiment: exact
Decimal fields rounded once to binary64 and fed to the unchanged Cartesian
FD observer. Shrinking the step eventually violates the1e-7 Hamiltonian
bound despite accurate input fields. This demonstrates a local rounding
floor, without attributing the complete old seed failure or waiving its gate.
`source-calibration-546b737e.py` preserves that experiment's exact earlier
source. `harness-layout-note.json` records the initial helper shape rejection
before comparison; it was not a native field failure.
The calibration did not capture runtime numerical-library hashes, so its
source binding is not a complete runtime provenance record.

`preservation.json` records post-run byte verification only. It must not be
described as a preflight check or a new numerical run.
