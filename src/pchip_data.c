#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "pchip_data.h"
#include "data.h"

// Read-only pointer to the coeffs
const double *pchip_coeffs = NULL;

// Local pointer to the coeffs
static double *g_coeffs = NULL;

int pchip_load_coeffs(const char *path) {
    return load_coeff_file(
        path,
        pchip_total_doubles,
        &g_coeffs,
        &pchip_coeffs
    );
}

void pchip_unload_coeffs(void) {
    unload_coeff_file(&g_coeffs, &pchip_coeffs);
}
