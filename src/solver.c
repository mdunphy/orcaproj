#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "solver.h"


int brent_solve(
    fsfun_t f, void *ctx,
    double a, double b, double tol, int maxiter,
    double *root_out, int *n_iter,
    const double *ffa, const double *ffb)
{
    double fa, fb;

    if (ffa) fa = *ffa; else fa = f(a, ctx);
    if (ffb) fb = *ffb; else fb = f(b, ctx);

    *n_iter = 0;

    if (fa * fb > 0.0)
        return 2;  // not a bracketing interval

    double c = a;
    double fc = fa;
    double d = b - a;
    double e = d;

    for (int iter = 0; iter < maxiter; iter++) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        // GSL-style effective tolerance
        double tol_act = 2.0 * DBL_EPSILON * fabs(b) + tol;
        double m = 0.5 * (c - b);

        // debug
//        printf("iter %2d:  a=% .12e  b=% .12e  c=% .12e  fa=% .3e  fb=% .3e  fc=% .3e\n",
//               iter, a, b, c, fa, fb, fc);

        if (fabs(c - b) <= 2.0 * tol_act || fb == 0.0) {
            *root_out = b;
            *n_iter = iter;
            return 0;
        }

        if (fabs(e) < tol_act || fabs(fa) <= fabs(fb)) {
            e = m;
            d = e;  // bisection
        } else {
            double p, q, s;
            s = fb / fa;
            if (a == c) {
                // Secant method
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                // Inverse quadratic interpolation
                double q1 = fa / fc;
                double r  = fb / fc;
                p = s * (2.0 * m * q1 * (q1 - r) - (b - a) * (r - 1.0));
                q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.0)
                q = -q;
            else
                p = -p;

            // GSL-style safeguard for interpolation
            if ((2.0 * p < 3.0 * m * q - fabs(tol_act * q)) &&
                (p < fabs(0.5 * e * q))) {
                e = d;
                d = p / q;
            } else {
                e = m;
                d = e;
            }
        }

        a = b;
        fa = fb;

        // Avoid zero-length steps
        if (fabs(d) > tol_act)
            b += d;
        else
            b += (m > 0.0 ? tol_act : -tol_act);

        fb = f(b, ctx);

        // Maintain bracket
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            e = b - a;
            d = e;
        }

        *n_iter = iter + 1;
    }

    *root_out = b;
    return 1;  // max iterations reached
}

int ridder_solve(
    fsfun_t f, void *ctx,
    double a, double b, double tol, int maxiter,
    double *root_out, int *n_iter,
    const double *ffa, const double *ffb)
{
    double fa = ffa ? *ffa : f(a, ctx);
    double fb = ffb ? *ffb : f(b, ctx);

    *n_iter = 0;

    // Require a valid bracket
    if (fa * fb > 0.0)
        return 2;  // invalid bracket

    // how many iterations we allow Ridder to do before handing off to Brent
    int max_ridder_iter = maxiter / 2;

    double left  = a;
    double right = b;

    for (int iter = 0; iter < max_ridder_iter; iter++) {

        double m  = 0.5 * (left + right);
        double fm = f(m, ctx);

        double tol_act = 2.0 * DBL_EPSILON * fabs(m) + tol;

        // Brent-style interval test: if already very tight, let Brent polish or accept midpoint
        if (fabs(right - left) <= 2.0 * tol_act) {
            int brent_iters = 0;
            int rc = brent_solve(f, ctx, left, right, tol,
                                 maxiter - iter, root_out, &brent_iters, &fa, &fb);
            *n_iter = iter + brent_iters;
            return rc;
        }

        // Function convergence test
        if (fabs(fm) < tol) {
            *root_out = m;
            *n_iter   = iter;
            return 0;
        }

        // Ridder discriminant
        double discr = fm * fm - fa * fb;

        if (discr <= 0.0) {
            // Fallback to bisection update if discr is degenerate
            if (fa * fm < 0.0) {
                right = m; fb = fm;
            } else {
                left  = m; fa = fm;
            }
            *n_iter = iter + 1;
            continue;
        }

        double s  = sqrt(discr);
        double sg = copysign(1.0, fa - fb);
        double x  = m + (m - left) * sg * fm / s;
        double fx = f(x, ctx);

        // Update bracket
        if (fm * fx < 0.0) {
            left  = m;  fa = fm;
            right = x;  fb = fx;
        } else if (fa * fx < 0.0) {
            right = x;  fb = fx;
        } else {
            left  = x;  fa = fx;
        }

        // Brent-style end conditions
        if (fabs(right - left) <= 2.0 * tol_act || fabs(fx) < tol) {
            // Handoff to Brent
            int brent_iters = 0;
            int rc = brent_solve(f, ctx, left, right, tol,
                                 maxiter - (iter + 1), root_out, &brent_iters, &fa, &fb);
            *n_iter = (iter + 1) + brent_iters;
            return rc;
        }

        *n_iter = iter + 1;
    }

    // If we finish the Ridder loop without satisfying criteria, let Brent finish from current bracket
    {
        int brent_iters = 0;
        int rc = brent_solve(f, ctx, left, right, tol,
                             maxiter - *n_iter, root_out, &brent_iters, &fa, &fb);
        *n_iter += brent_iters;
        return rc;
    }
}

double inverse_cubic_guess(
    double j0, double j1, double j2, double j3,
    double f0, double f1, double f2, double f3)
{
    double x = 0.0;

    x += j0 * (-f1)*(-f2)*(-f3) / ((f0-f1)*(f0-f2)*(f0-f3));
    x += j1 * (-f0)*(-f2)*(-f3) / ((f1-f0)*(f1-f2)*(f1-f3));
    x += j2 * (-f0)*(-f1)*(-f3) / ((f2-f0)*(f2-f1)*(f2-f3));
    x += j3 * (-f0)*(-f1)*(-f2) / ((f3-f0)*(f3-f1)*(f3-f2));

    return x;
}

int newton_refine_n(
    fsfun_df_t func,
    void *ctx,
    int nsteps,
    double ftol,
    double a,
    double b,
    const double *ffa,
    const double *ffb,
    double jc,
    double *root,
    int *niters)
{
    double j = jc;
    *niters = 0;

    double fa, fb;
    if (ffa) fa = *ffa; else func(a, ctx, &fa, &(double){0});
    if (ffb) fb = *ffb; else func(b, ctx, &fb, &(double){0});

    // Initial evaluation at j
    double f, df;
    func(j, ctx, &f, &df);

    for (int i = 0; i < nsteps; i++) {

        if (!isfinite(f) || !isfinite(df) || fabs(df) <= DBL_EPSILON)
            break;

        // Apply Newton method
        double jn = j - f / df;
        int use_newton = isfinite(jn) && jn > a && jn < b;

        // Bisect if Newton's guess was not usable
        if (!use_newton) {
            jn = 0.5 * (a + b);
        }

        // Evaluate
        double fn, dfn;
        func(jn, ctx, &fn, use_newton ? &dfn : &(double){0});

        // Maintain bracket
        if (fa * fn <= 0.0) {
            b = jn; fb = fn;
        } else {
            a = jn; fa = fn;
        }

        j = jn;
        f = fn;
        df = dfn;

        (*niters)++;

        if (fabs(f) < ftol)
            break;
    }

    *root = j;
    return 0;
}
