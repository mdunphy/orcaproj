#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <assert.h>

#include "chebcheb_data.h"
#include "data.h"

// Read-only pointer to the coeffs
const double *chebcheb_coeffs = NULL;

// Local pointer to the coeffs
static double *g_coeffs = NULL;


int chebcheb_load_coeffs(const char *path) {
    if (chebcheb_num_models > 64) {
        fprintf(stderr, "chebcheb error: Row-cache bitmask only supports up to 64 models!\n");
        abort();
    }
    return load_coeff_file(path, chebcheb_total_doubles, &g_coeffs, &chebcheb_coeffs);
}

void chebcheb_unload_coeffs(void) {
    unload_coeff_file(&g_coeffs, &chebcheb_coeffs);
}

static void extract_vertical_edge(cheb1d_edge_t *E)
{
    // Extract a 1D chebyshev series in j along i=-1 or +1 for the current tile E
    const int k = E->tile_index;
    const int x0 = E->x0;
    const chebcheb_model_t *M = &chebcheb_models[k];
    const int Px = M->Px;
    const int Py = M->Py;
    const double *C = get_C_ptr(k);

    const int nx = Px + 1;
    const int ny = Py + 1;

    E->Py   = Py;
    E->jmin = M->jmin;
    E->jmax = M->jmax;


    // Contract for f(i0, j)
    for (int q = 0; q < ny; ++q) {
        double sum = 0.0;

        if (x0 == 1) {
            for (int p = 0; p < nx; ++p)
                sum += C[p*ny + q];
        }
        else if (x0 == -1) {
            double s = 1.0;
            for (int p = 0; p < nx; ++p) {
                sum += s * C[p*ny + q];
                s = -s;
            }
        }
        E->f[q] = sum;
    }


    // Contract for df/di(x0, y)
    const double dx_di = 2.0 / (M->imax - M->imin);
    for (int q = 0; q < ny; ++q) {
        double sum = 0.0;
        double Tp_prime;

        for (int p = 0; p < nx; ++p) {
            if (x0 == 1) {
                Tp_prime = (p == 0) ? 0.0 : (double)(p*p);
            } else {
                if (p == 0) Tp_prime = 0.0;
                else Tp_prime = ((p & 1) ? +1.0 : -1.0) * (double)(p*p);
            }
            sum += C[p*ny + q] * Tp_prime;
        }
        // chain rule: df/di = df/dx * dx/di
        E->dfdi[q] = sum * dx_di;
    }

    // Contract for df/dj(x0, y)
    // df/dj as U-series coefficients
    {
        const double dy_dj = 2.0 / (M->jmax - M->jmin);
        for (int q = 1; q < ny; ++q)
            E->dfdj[q - 1] = q * E->f[q] * dy_dj;
        E->dfdj[ny-1] = 0.0;
    }
}

void chebcheb_init_vertical_edges(void)
{
    for (int n = 0; n < chebcheb_num_edge_tiles_UL; ++n) {
        extract_vertical_edge(&edge_tiles_UL[n]);
    }
    for (int n = 0; n < chebcheb_num_edge_tiles_UR; ++n) {
        extract_vertical_edge(&edge_tiles_UR[n]);
    }
}
