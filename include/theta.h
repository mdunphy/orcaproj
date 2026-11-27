#ifndef ORCA_TRIPOLAR_THETA_H
#define ORCA_TRIPOLAR_THETA_H

#include <stddef.h>

int ij_to_theta_vector(const double *ig_in,
                       const double *ej_in,
                       size_t n,
                       double *theta_out,
                       double *dtheta_di,
                       double *dtheta_dj,
                       int *s_out);

int ij_to_theta(double ig, double ej, double *theta_out, double *dtheta_di, double *dtheta_dj);

#endif /* ORCA_TRIPOLAR_THETA_H */
