# Overintegrated meridional operator candidate (not physically qualified)

Motivation: the frozen 16→32 azimuthal probe leaves native physical-equivalent
bulk momentum at about 1e-12, whereas the fully off-grid observer gives about
1e-7. Stable barycentric sampling leaves the latter unchanged. This candidate
tests meridional aliasing rather than adding another performance campaign.

The opt-in `HISPID_OVERINTEGRATED_GALERKIN` applies only to trumpet solves.
At N radial/polar Chebyshev roots, interpolate the modal P unknowns to
Q=ceil(3N/2) roots in each meridional direction. Evaluate the unchanged raw
four-equation CTT residual there. Restrict it with the discrete Chebyshev L2
projection R=(N/Q) B^T, then apply the original N-grid row weights and axis tau
conditions. Apply precisely the same transfers to the analytic JVP. This is
not an assertion that 3/2 oversampling exactly dealiasses all nonpolynomial
coefficients or powers. Physical acceptance still requires the independent
observer and the predeclared three-grid bounds. The original 0.15/0.5 interior
and original moderate observer are retained for this candidate.

One manufactured tensor control tests polynomial interpolation, round-trip
projection, and removal of an aliased quadratic high mode. Errors are
2.60e-13, 1.67e-13, and 6.22e-15. The initial prototype instead projected the
weighted equations: its small solve stalled at 2400 Krylov iterations with
relative residual 0.388; retained as control/cpu.json. Moving row equilibration
after projection converges in 5 Newton/3049 Krylov steps (37.72 s CPU). The
same small control takes only tens of Krylov steps with the original discrete
operator, so this is not yet a performance improvement. The unweighted
projection JVP agrees with centered differences to 7.40e-12/5.46e-12, and QI
residual/JVP remain bitwise equal to the retained baseline. All controls are
implementation diagnostics, not binary physical acceptance.

GPU transfers use preallocated tensor scratch and uploaded rectangular
matrices; the raw quadrature context has no preconditioner or Krylov basis.
Host seed geometry remains the qualified extended-precision builder. The
memory accounting removes a nonexistent simultaneous host Krylov bank from
Kokkos execution, while retaining host outer vectors, device basis/factors,
quadrature caches, and setup/upload peaks. The 64 GiB context budget and the
free-device-memory guard remain enforced. Reference-path memory bounds are
unchanged.

The production-resolution pilot is gated by one CPU/CUDA small-field
comparison (scaled difference <1e-10) and uses one shared 80 GiB GPU, OMP16.
It starts from the old original-formulation 192 checkpoint only as an explicit
migrated initial guess. It must solve again and remeasure every physical
quantity. Subsequent 224/240 and AthenaK runs require a passing physical pilot; a small
control alone does not qualify the binary.

Perlmutter source/build: source-galerkin-20261018 / build-galerkin-20261018
under /pscratch/sd/h/hzhu/codex-hispid-trumpet-20261005.
Job 59955085 exited during CMake generation: the source bundle omitted test
sources referenced by CMake. No solve ran. Corrected bundle includes tests/;
job 59955120, output galerkin-pilot-v2, is the first numerical attempt.

The default/OFF trumpet residual and JVP are bitwise identical to the previous
compensated-projection image on the frozen small witness. QI compatibility is
also bitwise exact. The CUDA control converges in5Newton/3131Krylov,6.758solve
seconds; CPU/CUDA sampled-field scaled difference is5.25419e-13 (<1e-10).
The CUDA image is bff0cb239d8698f1be586c734671c946ec2b5bb088f7ac4d807ce2e9699e8d87.
These controls establish the implementation comparison, not physical acceptance.

Job59955120 completed14:38 and released its allocation. The192 pilot fails its first linear solve: 2400 iterations, true relative residual 0.5059823. No new accepted Newton update. Setup 211.720s, solve 444.507s, total 768.299s. The raw physical observer and failed status are retained in evidence192/. No224/240 or horizon run follows. This candidate is rejected as a convergence remedy.
