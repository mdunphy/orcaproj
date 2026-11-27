#ifndef SOLVER_H
#define SOLVER_H

#include <stddef.h>

typedef enum {
    SOLVER_BRENT,
    SOLVER_RIDDERS,
    SOLVER_NEWTON,
    SOLVER_HYBRID,
} solver_t;

typedef double (*fsfun_t)(double, void *);
typedef void   (*fsfun_df_t)(double, void *, double *, double *);


int brent_solve(
    fsfun_t f, void *ctx,
    double a, double b, double tol, int maxiter,
    double *root_out, int *n_iter,
    const double *ffa, const double *ffb);

int ridder_solve(
    fsfun_t f, void *ctx,
    double a, double b, double tol, int maxiter,
    double *root_out, int *n_iter,
    const double *ffa, const double *ffb);


double inverse_cubic_guess(
    double j0, double j1, double j2, double j3,
    double f0, double f1, double f2, double f3);


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
    int *niters);

#endif /* SOLVER_H */
