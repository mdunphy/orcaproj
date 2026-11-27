#ifndef ORCA_TRIPOLAR_FUNCTIONS_H
#define ORCA_TRIPOLAR_FUNCTIONS_H

#include <math.h>

typedef struct {
    double THETA_MIN;
    double THETA_MAX;
    double IG_MIN;
    double IG_MAX;
    double JG_MIN;
    double JG_MAX;
    double I_MIN;
    double I_MAX;
    double EJ_MIN;
    double EJ_MAX;
    double IEX_L;
    double IEX_R;
    double J_TOL;
} orca_bounds_t;

typedef struct {
    double lon;
    double lat;
    double y;
} orca_pole_t;

// TODO: combine all of these into a single struct to consolidate global namespace
extern const orca_pole_t asian_pole;
extern const orca_pole_t canadian_pole;
extern const orca_bounds_t bounds;
extern const double rlampo;


double fsfun(double pj, double px, double py);
double fsfun_fast(double ej, double px, double py);
double fsfun_fast_x0(double ej, double py);

void aby0(double ej, double *a, double *b, double *y0);

void fsfun_fast_with_deriv(double ej, double px, double py, double *out_val, double *out_der);
void fsfun_fast_x0_with_deriv(double ej, double py, double *out_val, double *out_der);

double fsa(double pj);
double fsb(double pj);
double fsy0p(double pj);
double fsy0(double pj);
double fsj(double pj);
double fsej(double ej);
double phi2j(double lat);
double fslam(double pi);
double lam2i(double lam);


void fslamphixy(double lon_deg, double lat_deg, double *x_out, double *y_out);
void fsxylamphi(double px, double py, double *plam_deg, double *pphi_deg);


// Quick wrappers
static inline double fsfun_a2b2y0(double a2, double b2, double y0, double px, double py)
{
    const double dy = py - y0;
    //return a2*px*px + b2*(dy*dy) + a2*b2;
    return a2 * px * px + b2 * dy * dy - a2 * b2;

}

static inline double fsfun_a2y0(double a2, double b2, double y0, double py)
{
    const double dy = py - y0;
    return dy*dy -a2 + 0.0*b2;
}

#endif /* ORCA_TRIPOLAR_FUNCTIONS_H */
