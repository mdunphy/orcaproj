/* chebcheb.c — C evaluator for exported chebcheb (2D Chebyshev) data */

#include <math.h>
#include <string.h>
#include <stdint.h>

#include "chebcheb.h"
#include "chebcheb_data.h"
#include "common.h"

int find_covering_models(double i, double j, int out_k[2])
{
    int n = 0;

    const chebcheb_strip_t *strips = chebcheb_strips;

    for (int s = 0; s < chebcheb_num_strips && n < 2; ++s) {

        const chebcheb_strip_t *S = &strips[s];

        if (j < S->jmin) break;
        if (j > S->jmax) continue;

        int k0 = S->k0;
        int k1 = S->k1;

        for (int k = k0; k <= k1; ++k) {

            const chebcheb_model_t *M = &chebcheb_models[k];

            double imin = M->imin;
            double imax = M->imax;

            if (i > imax) continue;
            if (i <= imin && imin != 0.0) continue;
            if (i <  imin && imin == 0.0) continue;

            out_k[n++] = k;
            break;
        }
    }

    return n;
}

static inline double to_pm1(double x, double xmin, double xmax)
{
    double d = xmax - xmin;
    if (d==0.0) return -1.0;
    return -1.0 + 2.0 * (x - xmin) / d;
}


// 1D Clenshaw for Chebyshev T basis
static inline double cheb_clenshaw_T(const double *restrict c, int n, double x)
{
    if (n <= 0) return 0.0;
    if (n == 1) return c[0];
    double b1 = 0.0, b2 = 0.0;
    const double x2 = x + x;
    for (int m = n - 1; m >= 1; --m) {
        //double b0 = x2 * b1 - b2 + c[m];
        double b0 = fma(x2, b1, c[m] - b2);
        b2 = b1; b1 = b0;
    }
    //return x * b1 - b2 + c[0];
    return fma(x, b1, c[0] - b2);
}


// Compute t[k] = T_k(y), k = 0..Ny-1
static inline void cheb_T_forward(double y, int Ny, double *restrict t)
{
    if (Ny <= 0) return;
    t[0] = 1.0;
    if (Ny == 1) return;

    t[1] = y;
    const double y2 = y + y;

    for (int k = 1; k < Ny - 1; ++k) {
        t[k + 1] = fma(y2, t[k], -t[k - 1]);
    }
}


static inline double cheb_clenshaw_U(const double *c, int n, double y)
{
    if (n <= 0) return 0.0;
    // U_n recurrence:  U_{n}(y) = 2*y*U_{n+1} - U_{n+2} + c[n]
    double u0 = 0.0;
    double u1 = 0.0;
    for (int k = n - 1; k >= 0; --k) {
        double u2 = u1;
        u1 = u0;
        u0 = 2.0*y*u1 - u2 + c[k];
    }
    return u0;
}

///* Contract a 2-D Chebyshev coefficient block along the j-dimension.
//static inline void cheb_contract_y_clenshaw(const double * restrict C,
//                                            int nx, int ny,
//                                            double y,
//                                            double * restrict ax)
//{
//    for (int p = 0; p < nx; ++p) {
//        const double *colp = C + (size_t)p * (size_t)ny;
//        ax[p] = cheb_clenshaw_T(colp, ny, y);
//    }
//}

// This is the same as cheb_contract_y_clenshaw, but faster
static inline void cheb_contract_y_forward(const double * restrict C,
                                           int nx, int ny,
                                           double y,
                                           double * restrict ax)
{
    // Compute Chebyshev basis at y
    double Ty[CHEBCHEB_MAX_AX];
    cheb_T_forward(y, ny, Ty);

    // Apply dot products, 4 at a time
    for (int p = 0; p < nx; ++p) {

        const double * restrict row = C + (size_t)p * (size_t)ny;

        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;
        double s3 = 0.0;

        int q = 0;

        for (; q + 3 < ny; q += 4) {
            s0 += row[q]     * Ty[q];
            s1 += row[q + 1] * Ty[q + 1];
            s2 += row[q + 2] * Ty[q + 2];
            s3 += row[q + 3] * Ty[q + 3];
        }

        double s = s0 + s1 + s2 + s3;

        for (; q < ny; ++q) {
            s += row[q] * Ty[q];
        }

        ax[p] = s;
    }
}

// Reset the row cache
void chebcheb_init_rowcache(row_cache_t *R, double j0)
{
    R->filled = 0;
    R->j0 = j0;
}

// Evaluate model k at coordinates (i,j)
double eval_cheb_tile(int k, double i, double j, row_cache_t *R)
{
    const uint64_t mask = (uint64_t)1 << k; // bitmask to mark cache present per tile
    if (R && (R->filled & mask)) {
        // Cache object provided AND tile k is in the cache
        row_cache_entry_t *E = &R->cache[k];
        const double x = i * E->scale_x + E->shift_x;
        return cheb_clenshaw_T(E->ax, E->nx, x);
    }
    // No cache object provided OR cache object provided but tile k not in cache
    const chebcheb_model_t *M = &chebcheb_models[k];
    const int nx = M->Px + 1;
    const int ny = M->Py + 1;
    const double imin  = M->imin;
    const double imax  = M->imax;
    const double x = to_pm1(i, imin, imax);
    const double y = to_pm1(j, M->jmin, M->jmax);
    row_cache_entry_t *E = NULL;
    if (R) {
        E = &R->cache[k];
        E->nx = nx;
        E->scale_x = 2.0 / (imax - imin);
        E->shift_x = -(imax + imin) / (imax - imin);
    }
    double ax_temp[CHEBCHEB_MAX_AX];
    double *ax = R ? E->ax : ax_temp;
    const double *C = get_C_ptr(k);
    cheb_contract_y_forward(C, nx, ny, y, ax); // Contract in y
    if (R) R->filled |= mask; // mark tile entry as present
    return cheb_clenshaw_T(ax, nx, x); // Evalate at x
}


/* Differentiate Chebyshev T-series coefficients c[0..n-1] into dc[0..n-2], where n >=1. */
static inline void cheb_diff_1d(const double *restrict c, int n, double *restrict dc)
{
    if (n <= 1) {
        dc[0] = 0.0;
        return;
    }
    const int P = n-1;
    dc[P-1] = 2.0 * P * c[P];
    if (P > 1)
        dc[P-2] = 2.0 * (P-1) * c[P-1];

    for (int k = P-3; k >= 0; --k)
        dc[k] = dc[k+2] + 2.0 * (k+1) * c[k+1];

    dc[0] *= 0.5;
}


static inline void cheb_contract_x_clenshaw(const double * restrict C,
                                            int nx, int ny,
                                            double x,
                                            double * restrict ay)
{
    for (int q = 0; q < ny; ++q) {
        double col[CHEBCHEB_MAX_AX];
        for (int p = 0; p < nx; ++p)
            col[p] = C[p*ny + q];

        ay[q] = cheb_clenshaw_T(col, nx, x);
    }
}

// Compute df/di and df/dj at (i,j) for tile k
void eval_tile_grad(int k, double i, double j,
                    double *fx, double *fy)
{
    const chebcheb_model_t *M = &chebcheb_models[k];
    const int nx = M->Px + 1;
    const int ny = M->Py + 1;

    const double *C = get_C_ptr(k);

    const double x = to_pm1(i, M->imin, M->imax);
    const double y = to_pm1(j, M->jmin, M->jmax);

    // df/di
    double ax[CHEBCHEB_MAX_AX] = {0};
    double dax[CHEBCHEB_MAX_AX];

    cheb_contract_y_forward(C, nx, ny, y, ax);
    cheb_diff_1d(ax, nx, dax);

    const double dxi = 2.0 / (M->imax - M->imin);
    *fx = cheb_clenshaw_T(dax, nx-1, x) * dxi;

    // df/dj
    double ay[CHEBCHEB_MAX_AX];
    double day[CHEBCHEB_MAX_AX];

    cheb_contract_x_clenshaw(C, nx, ny, x, ay);
    cheb_diff_1d(ay, ny, day);

    const double dyi = 2.0 / (M->jmax - M->jmin);
    *fy = cheb_clenshaw_T(day, ny-1, y) * dyi;
}

// Compute f and derivatives along the edge tiles that connect to PCHIP tiles
void cheb1d_eval_f_and_dfdi_and_dfdj(const cheb1d_edge_t *E,
                                     double j,
                                     double *f,
                                     double *dfdi,
                                     double *dfdj)
{
    const int ny = E->Py + 1;
    const double y = to_pm1(j, E->jmin, E->jmax);
    *f    = cheb_clenshaw_T(E->f,    ny, y);
    *dfdi = cheb_clenshaw_T(E->dfdi, ny, y);
    *dfdj = cheb_clenshaw_U(E->dfdj, ny, y);
}
