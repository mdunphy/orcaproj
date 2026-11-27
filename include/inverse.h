#ifndef ORCA_TRIPOLAR_INVERSE_H
#define ORCA_TRIPOLAR_INVERSE_H

#include <stddef.h>

int ij_to_lonlat_vector(const double *i_in,
                        const double *j_in,
                        size_t n,
                        double *lon_out,
                        double *lat_out,
                        int *s_out);

int ij_to_lonlat(double i, double j, double *lon_out, double *lat_out);

#endif /* ORCA_TRIPOLAR_INVERSE_H */
