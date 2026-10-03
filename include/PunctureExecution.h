#ifndef PUNCTURE_EXECUTION_H
#define PUNCTURE_EXECUTION_H
#include "PunctureKrylov.h"
#ifdef __cplusplus
extern "C" {
#endif
/* The historical entry points always default to REFERENCE. KOKKOS denotes
 * the execution space compiled into this image (Serial, OpenMP or CUDA).
 * Runtime/geometry setup is CPU work; linear vectors/operators stay resident.
 * No implicit fallback to a different execution space or preconditioner. */
/* Kokkos API operations and initialization are serialized process-wide.
 * Do not mutate/destroy a context concurrently with a caller using it.
 * Allocation peaks include Kokkos allocations only; pinned host memory is
 * host, managed/unified spaces are counted by their logical device space. */
enum { PUNCTURE_REFERENCE=0, PUNCTURE_KOKKOS=1 };
const char *Puncture_execution_name(void);
/* Initialize once, before creating contexts. threads=0 uses Kokkos defaults.
 * Returns -1 for a non-Kokkos image or conflicting repeated initialization. */
int Puncture_execution_initialize(int threads);
int Puncture_execution_concurrency(void); /* actual host execution concurrency */
unsigned long long Puncture_execution_device_free_bytes(void); /* 0 for CPU */
/* JSON describing the actual visible CUDA device; -1 on CPU/query failure. */
int Puncture_execution_device_description(char *buffer,int capacity);
typedef struct {
  double operator_seconds,precondition_seconds,setup_seconds,linear_seconds;
  double transfer_seconds,nonlinear_seconds;
  unsigned long long resident_bytes,peak_workspace_bytes;
  unsigned long long host_to_device_bytes,device_to_host_bytes;
  unsigned long long kokkos_host_current_bytes,kokkos_host_peak_bytes;
  unsigned long long kokkos_device_current_bytes,kokkos_device_peak_bytes;
  unsigned long long nonlinear_calls,workspace_creations;
  int memory_tracking_available;
} PunctureExecutionStats;
int Puncture_execution_statistics(PunctureExecutionStats *);
void Puncture_execution_reset_statistics(void);
/* Separate from the historical TP_SolverStats ABI. Counts cover Newton,
 * including original-residual confirmation and any subsequent polishing. */
typedef struct {
  int workspace_creations,device_residual_calls,original_confirmation_calls,polishing_steps;
  double device_candidate_linf,original_confirmation_linf;
} TP_ExecutionStats;
void TP_solver_get_execution_statistics(TP_ExecutionStats *);
#ifdef PUNCTURES_KOKKOS
struct TP_Modal;
struct DERIVS;
typedef struct Puncture_BYWorkspace Puncture_BYWorkspace;
Puncture_BYWorkspace *Puncture_kokkos_BY_create(int na,int nb,int np);
void Puncture_kokkos_BY_destroy(Puncture_BYWorkspace *);
int Puncture_kokkos_BY_residual(Puncture_BYWorkspace *,struct DERIVS *v,double *F,struct DERIVS *u);
int Puncture_kokkos_BY_linear(Puncture_BYWorkspace *,const double *physical_u,
                            struct TP_Modal *,const double *rhs,double *solution,
                            const PK_Options *,PK_Result *,PK_Monitor,void *);
int Puncture_kokkos_BY(int na,int nb,int np,const double *physical_u,
                       struct TP_Modal *,const double *rhs,double *solution,
                       const PK_Options *,PK_Result *,PK_Monitor,void *);
#endif
#ifdef __cplusplus
}
#endif
#endif
