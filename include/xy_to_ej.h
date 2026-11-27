#ifndef ORCA_TRIPOLAR_XY_TO_EJ_H
#define ORCA_TRIPOLAR_XY_TO_EJ_H

#include <stddef.h>
#include "internal_errors.h"
#include "solver.h"

orca_int_error_t xy_to_ej(double x, double y, double *ej_out, int *niter_out);
orca_int_error_t xy_to_ej_core(double x, double y, solver_t solver, int maxiter, double *ej_out, int *niter_out);

void fill_aby0_cache(void);

#endif /* ORCA_TRIPOLAR_XY_TO_EJ_H */
