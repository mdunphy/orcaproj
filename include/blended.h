/* blended.h — Interface for blended evaluator/solver stitching PCHIP and Chebyshev regions. */

#ifndef ORCA_TRIPOLAR_BLENDED_H
#define ORCA_TRIPOLAR_BLENDED_H

#include "internal_errors.h"

orca_int_error_t blended_eval_ij(double i, double j, double *theta_out,
                                 double *dtheta_di, double *dtheta_dj);

orca_int_error_t blended_find_i(double j, double theta, double lb, double ub,
                                double tol, int maxiter, double *i_out, int *niter);

#endif /* ORCA_TRIPOLAR_BLENDED_H */
