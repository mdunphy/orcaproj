#ifndef ORCA_TRIPOLAR_FORWARD_H
#define ORCA_TRIPOLAR_FORWARD_H

#include <stddef.h>

int lonlat_to_ij_vector(const double *lon_deg,
                        const double *lat_deg,
                        size_t n,
                        double *i_out,
                        double *j_out,
                        int *s_out);

int lonlat_to_ij(double lon_deg, double lat_deg,
                 double *i_out, double *j_out);

#endif /* ORCA_TRIPOLAR_FORWARD_H */
