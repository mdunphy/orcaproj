#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

#include "xy_to_ej.h"
#include "functions.h"
#include "solver.h"
#include "errors.h"
#include "common.h"


// Goal:  given x,y we need to find ej (ellipse j, ej = jg - 290)
// Stategy: apply bisection using cached a,b and y0 values to
//          reduce the bracket, then apply a root solver to refine the root.
// SOLVER_CONFIG 1 and 2 are working well, with 2 being generally faster
// SOLVER_CONFIG 3 is faster still but not as robust
// SOLVER_CONFIG 4 is a hybrid, newton for low J and ridders for high J
// TODO: rigorous testing to converge on a single optimal config

#define SOLVER_CONFIG 2

#if SOLVER_CONFIG == 1
    #define NBISECT 10
    #define SOLVER  SOLVER_BRENT
    #define MAXITER 100

#elif SOLVER_CONFIG == 2
    #define NBISECT 10
    #define SOLVER  SOLVER_RIDDERS
    #define MAXITER 100

#elif SOLVER_CONFIG == 3
    #define NBISECT 12
    #define SOLVER  SOLVER_NEWTON
    #define MAXITER 4

#elif SOLVER_CONFIG == 4
    #define NBISECT 14
    #define SOLVER  SOLVER_HYBRID
    #define MAXITER 30
#else
    #error "Invalid SOLVER_CONFIG"
#endif

#define NPTS (((size_t)1 << NBISECT) + 1)

typedef struct {
    double ej, a2, b2, y0;
} ABY0vals;

static ABY0vals aby0_cache[NPTS];

void fill_aby0_cache(void) {
    const double ejmin=bounds.EJ_MIN;
    const double ejmax=bounds.EJ_MAX;
    const size_t nseg = (size_t)1 << NBISECT;
    const double dj = (ejmax - ejmin) / (double)nseg;
    double a,b,y0;
    for (size_t k = 0; k < NPTS; k++) {
        double ej = ejmin + (double)k * dj;
        if (k == NPTS - 1)  ej = ejmax;
        aby0(ej, &a, &b, &y0);
        aby0_cache[k].a2 = a*a;
        aby0_cache[k].b2 = b*b;
        aby0_cache[k].y0 = y0;
        aby0_cache[k].ej = ej;
    }
}

typedef struct {
    double x,y;
    int x0flag;
} FsXYCtx;

// Wrapper for fsfun
static double fsfun_wrap(double pj, void *pctx)
{
    const FsXYCtx *c = (const FsXYCtx *)pctx;
    if (c->x0flag) {
        return fsfun_fast_x0(pj, c->y);
    }
    return fsfun_fast(pj, c->x, c->y);
}

// Wrapper for fsfun based on cache
static inline double fsfun_cache(int x0flag, size_t k, double px, double py)
{
    const double a2 = aby0_cache[k].a2;
    const double y0 = aby0_cache[k].y0;
    const double b2 = aby0_cache[k].b2;
    if (x0flag) {
        return fsfun_a2y0(a2, b2, y0, py);
    }
    return fsfun_a2b2y0(a2, b2, y0, px, py);
}

// Wrapper for fsfun with derivative
static inline void fsfun_wrap_with_deriv(double j, void *pctx, double *f, double *df)
{
    const FsXYCtx *c = (const FsXYCtx *)pctx;
    if (c->x0flag)
        fsfun_fast_x0_with_deriv(j, c->y, f, df);
    else
        fsfun_fast_with_deriv(j, c->x, c->y, f, df);
}


static void bisect(int x0flag, double ftol, double x, double y, double *jlo, double *jhi, double *flo, double *fhi, size_t *ilo, size_t *ihi)
{
    size_t lo = 0;
    size_t hi = NPTS - 1;
    double f_lo = fsfun_cache(x0flag, lo, x, y);
    double f_hi = fsfun_cache(x0flag, hi, x, y);

    if (fabs(f_lo) < ftol) {
        *jlo = aby0_cache[lo].ej;
        *jhi = aby0_cache[lo+1].ej;
        *ilo = lo;
        *ihi = lo+1;
        *flo = f_lo;
        *fhi = fsfun_cache(x0flag, lo+1, x, y);
        return;
    }
    if (fabs(f_hi) < ftol) {
        *jlo = aby0_cache[hi-1].ej;
        *jhi = aby0_cache[hi].ej;
        *ilo = hi-1;
        *ihi = hi;
        *flo = fsfun_cache(x0flag, hi-1, x, y);
        *fhi = f_hi;
        return;
    }

    if (f_lo * f_hi > 0.0) {
        double xjlo=aby0_cache[lo].ej;
        double xjhi=aby0_cache[hi].ej;
        printf("BISECT NO ROOT!: x    %+.24f         y %+.24f\n",x,y);
        printf("BISECT NO ROOT!: xjlo %+.24f      xjhi %+.24f\n",xjlo,xjhi);
        printf("BISECT NO ROOT!: flo  %+.24f      fhi  %+.24f\n",f_lo,f_hi);
        *jlo = NAN;
        *jhi = NAN;
        *ilo = lo;
        *ihi = hi;
        *flo = f_lo;
        *fhi = f_hi;
        return;
    }
    for (size_t i = 0; i < NBISECT; ++i) {
        size_t mid = (lo + hi) >> 1;
        double fmid = fsfun_cache(x0flag, mid, x, y);
        if (f_lo * fmid <= 0.0) {
            hi = mid;
            f_hi = fmid;
        } else {
            lo = mid;
            f_lo = fmid;
        }
    }
    *jlo = aby0_cache[lo].ej;
    *jhi = aby0_cache[hi].ej;
    *ilo = lo;
    *ihi = hi;
    *flo = f_lo;
    *fhi = f_hi;
}

// call the core with the selected solver
orca_int_error_t xy_to_ej(double x, double y, double *ej_out, int *niter_out) {
    return xy_to_ej_core(x, y, SOLVER, MAXITER, ej_out, niter_out);
}

orca_int_error_t xy_to_ej_core(double x, double y, solver_t solver, int maxiter, double *ej_out, int *niter_out)
{
    if (!ej_out || !niter_out) return ORCA_INT_EFAULT;

    const double ftol0= 1e-15;
    const double ftol = 1e-15;
    const double xtol = 1e-15;
    double a0 = bounds.EJ_MIN;
    double b0 = bounds.EJ_MAX;
    double f1, f2;
    size_t lo=0, hi=NPTS-1;
    *niter_out = 0;

    // Select function to invert
    int x0flag = (fabs(x) < xtol) && (y < canadian_pole.y-xtol || y > asian_pole.y+xtol);
    //x0flag=0;
    FsXYCtx ctx = {.x = x, .y = y, .x0flag = x0flag};

    // Conduct the fast cache-base bisection
    bisect(x0flag, ftol, x, y, &a0, &b0, &f1, &f2, &lo, &hi);

    // Endpoint roots
    if (fabs(f1) < ftol0) {
        *ej_out = a0;
        return ORCA_INT_OK;
    }
    if (fabs(f2) < ftol0) {
        *ej_out = b0;
        return ORCA_INT_OK;
    }

    // Bracket test
    if (f1 * f2 > 0.0) {
        *ej_out = -4;
        //printf("X %.10g Y %.10g EJ %.20f |EFSFUN 1 %.20f |EFSFUN2 %.20f | a0 %.20f b0 %.20f\n",x, y, f1, f2, a0, b0);
        return ORCA_INT_ROOT_NO_BRACKET;
    }

    // Now solve the root
    double r = a0;
    int iters = 0;
    int rc=0;

    if (solver == SOLVER_BRENT) {
        rc = brent_solve(fsfun_wrap, &ctx, a0, b0, ftol, maxiter, &r, &iters, &f1, &f2);
    }
    else if (solver == SOLVER_RIDDERS) {
        rc = ridder_solve(fsfun_wrap, &ctx, a0, b0, ftol, maxiter, &r, &iters, &f1, &f2);
    }
    else if ( (solver == SOLVER_NEWTON) || (solver == SOLVER_HYBRID)) {
        // Compute the secant, fallback to midpoint if secant fails
        double denom = f2 - f1;
        double jc;
        if (fabs(denom) > 1e-14)
            jc = a0 - f1*(b0-a0)/denom;
        else
            jc = 0.5*(a0 + b0);

        // If we are not in the first or last interval, try an inverse cubic as an improvement over secant
        if ( (lo>0) && (hi<NPTS-1) ) {
            double fm = fsfun_cache(x0flag, lo-1, x, y);
            double fp = fsfun_cache(x0flag, hi+1, x, y);
            double jm = aby0_cache[lo-1].ej;
            double jp = aby0_cache[hi+1].ej;
            double jc_cubic = inverse_cubic_guess(jm,a0,b0,jp,fm,f1,f2,fp);
            if (isfinite(jc_cubic) && jc_cubic > a0 && jc_cubic < b0) jc=jc_cubic;
        }

        // This shouldn't be needed but put it anyway
        if (jc <= a0 || jc >= b0)
            jc = 0.5*(a0 + b0);

        else if (solver == SOLVER_NEWTON) {
            rc = newton_refine_n(fsfun_wrap_with_deriv,&ctx,maxiter,ftol,a0,b0,&f1,&f2,jc,&r,&iters);
        }
        else if (solver == SOLVER_HYBRID) {
            if (b0 > 205) rc = brent_solve(fsfun_wrap, &ctx, a0, b0, ftol, maxiter, &r, &iters, &f1, &f2);
            else rc = newton_refine_n(fsfun_wrap_with_deriv,&ctx,maxiter,ftol,a0,b0,&f1,&f2,jc,&r,&iters);
        }
    }
    else return ORCA_INT_EFAULT;
    *niter_out = iters;

    if (rc != 0) {
        *ej_out = r;
        return ORCA_INT_ROOT_FAIL;
    }

    *ej_out = r;
    return ORCA_INT_OK;
}
