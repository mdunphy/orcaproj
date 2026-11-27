#ifndef ORCA_DATA_H
#define ORCA_DATA_H
#include <stdint.h>

int load_coeff_file(const char *path, uint64_t expected_count, double **out_buf, const double **out_const);
void unload_coeff_file(double **buf, const double **pub);

#endif /* ORCA_DATA_H */
