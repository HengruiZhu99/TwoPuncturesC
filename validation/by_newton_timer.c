/* Timing wrapper only: compile this instead of TP_Newton.o. */
#define Newton benchmark_Newton_impl
#include "../src/TP_Newton.c"
#undef Newton
#include <time.h>

static double elapsed_newton;
static double clock_seconds(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9*t.tv_nsec;
}
void Newton(int nvar, int n1, int n2, int n3, derivs *v,
            double tol, int itmax) {
  const double start = clock_seconds();
  benchmark_Newton_impl(nvar, n1, n2, n3, v, tol, itmax);
  elapsed_newton += clock_seconds() - start;
}
double benchmark_newton_seconds(void) { return elapsed_newton; }
