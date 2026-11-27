#include <math.h>
#include "common.h"
#include "functions.h"
#include <stdio.h>

static const double rad       = M_PI / 180.0;        // radians per degree

const double rlampo = 253.0;  // central meridian (degrees)
static const double rdx       =   0.5;
static const double rjpnth    = 220.5;  // j index of the northern most ellipse (north fold)
static const double rjpeqt    = 250.0;
static const double rjpeq     = 291.0;
static const double rjpnthm1  = 219.5;
static const double ra1       =   1.0e-5;
static const double ra2       =   2.60269747511345;
static const double rb1       =   0.1973264041404091161706746;
static const double rb2       =   1.0;
static const double ry0       =  -0.40858412281692569; //ry0 = fsphiy(rphify0) / (fsphiy(zy0pjn) - fsa(rjpnth))

static const double fsa_rjpnth  = 0.2882633979681122893268252; // fsa(rjpnth)
static const double fsy0_rjpnth = 0.0757068362980900660952343; // fsy0(rjpnth)

const orca_pole_t asian_pole = {
    .lon = 73.0,
    .lat = 50.0,
    .y = fsy0_rjpnth + fsa_rjpnth
};

const orca_pole_t canadian_pole = {
    .lon = -107.0,
    .lat =   66.0,
    .y = fsy0_rjpnth - fsa_rjpnth
};


const orca_bounds_t bounds = {
    .THETA_MIN =     M_PI_2,
    .THETA_MAX = 3.0*M_PI_2,
    .IG_MIN    =   0.0,
    .IG_MAX    = 720.0,
    .JG_MIN    =   1.0,
    .JG_MAX    = 510.5,
    .I_MIN     =   0.0,
    .I_MAX     = 360.0,
    .EJ_MIN    =   1.0,
    .EJ_MAX    = rjpnth, // =220.5
    .J_TOL     = 1.0e-12
};

static inline double fsfa(double ej);
static inline double fsfb(double ej);
static inline double fsphia(double ej);
static inline double fsphib(double ej);
static inline double theta(double pj);
static inline double fsphiy(double pphi);

double fsfun(double ej, double px, double py) {
    double a = fsa(ej);
    double b = fsb(ej);
    double y0 = fsy0(ej);
    double a2 = a * a;
    double b2 = b * b;
    double dy = py - y0;
    return a2 * px * px + b2 * dy * dy - a2 * b2;
}

double fsfun_fast(double ej, double px, double py)
{
    // Wrapper for optimized aby0
    double a,b,y0;
    aby0(ej,&a,&b,&y0);

    const double a2 = a * a;
    const double b2 = b * b;
    const double dy = py - y0;
    //return a2 * (px*px - b2) + b2 * (dy*dy);
    return a2 * px * px + b2 * dy * dy - a2 * b2;

}

double fsfun_fast_x0(double ej, double py)  // for b=0
{
    // Wrapper for optimized aby0
    double a,b,y0;
    aby0(ej,&a,&b,&y0);
    const double dy = py - y0;
    return dy*dy -a*a;
}


void aby0(double ej, double *restrict a, double *restrict b, double *restrict y0)
{
    // Optimized calculation of A,B,Y0
    // Makes use of some trig identities
    const double ejm1 = ej - 1.0;
    const double x_fa  = ejm1 * ra2 / rjpnthm1;
    const double x3_fa = x_fa * x_fa * x_fa;
    const double denom_fa = 3.0 * (x3_fa / (ra1 * ra2) + 1.0);
    const double fa = (atan(x_fa) + x3_fa / denom_fa) * rjpnthm1 / ra2;

    // fsphia(ej) = exp(-theta_a)
    const double theta_a = rdx * rad * (fa + rjpeq - rjpeqt);
    const double aval = exp(-theta_a);
    *a=aval;

    const double t_fb = ejm1 / rjpnthm1;
    const double x_fb = t_fb * t_fb * t_fb;
    const double fsfb_val = 1.0 + rb1 * tanh(x_fb);
    const double theta_b = rdx * rad * (ejm1 + rjpeq - rjpeqt);
    const double s = exp(-theta_b);
#if 1
    // use tanh = (1-s*s)/(1+s*s)
    const double s2 = s * s;
    const double tanh_theta_b = (1.0 - s2) / (1.0 + s2);
#else
    const double tanh_theta_b = tanh(theta_b);
#endif
    const double fsy0p_b    = asin(tanh_theta_b) / rad;
    const double fsphib_val = fsfb_val * fsy0p_b;
    *b = tan(M_PI_4 - rad * 0.5 * fsphib_val);

    *y0 = (s - aval) * ry0;
}



static inline double fsphiy(double pphi) {
    // stereographic phi to radius, as in fslamphixy()
    return tan(M_PI / 4.0 - rad * pphi / 2.0);
}

static inline double fsphia(double ej) {
    // ej in [1, rjpnth]
    return fsy0p(fsfa(ej) + rjpeq);
}

static inline double fsphib(double ej) {
    // ej in [1, rjpnth]
    double pj = fsej(ej);
    return fsfb(ej) * fsy0p(pj);
}

double fsa(double ej) {
    // ej in [1, rjpnth]
    return fsphiy(fsphia(ej));
}

double fsb(double ej) {
    // ej in [1, rjpnth]
    return fsphiy(fsphib(ej));
}

double fsy0p(double pj) {
    // pj is in [1, 510.5]
    // This is inverse Mercator function
    return asin(tanh(theta(pj))) / rad;
}

static inline double fsfa(double ej) {
    // ej in [1, rjpnth]
    const double x = ra2 * (ej - 1.0) / rjpnthm1;
    /* arctan(x) + x^3 / (3 * (x^3 / (ra1*ra2) + 1)) all scaled by rjpnthm1/ra2 */
    const double x3 = x * x * x;
    const double denom = 3.0 * (x3 / (ra1 * ra2) + 1.0);
    return (atan(x) + x3 / denom) * rjpnthm1 / ra2;
}

static inline double fsfb(double ej) {
    // ej in [1, rjpnth]
    const double t = ( (ej - 1.0) / rjpnthm1 );
    const double x = t * t * t; /* ((j - 1)/rjpnthm1)^3 */
    return 1.0 + rb1 * tanh(rb2 * x);
}

double fsy0(double ej) {
    // ej in [1, rjpnth]
    double pj = fsej(ej);
    return (fsphiy(fsy0p(pj)) - fsa(ej)) * ry0;
}

static inline double theta(double pj) {
    return rdx * rad * (pj - rjpeqt);
}

double phi2j(double lat)
{
    const double phi = lat * rad;
    const double theta = atanh(sin(phi));
    const double j = rjpeqt + theta / (rdx * rad);
    return j;
}

/* Compute x0 from projection index pi */
static double fsx0p(double pi)
{
    return 90.0 - rdx * (pi - 1.0);
}

///* Derivative of fsx0p with respect to pi */
//static double fsdx0p(double pi)
//{
//    (void)pi; /* unused */
//    return -rdx;
//}

/* Compute longitude offset from projection index pi */
double fslam(double pi)
{
    return 90.0 - fsx0p(pi);
}

/* Convert longitude to i index */
double lam2i(double lam)
{
    return lam / rdx + 1.0;
}

double fsj(double pj)
{
    return pj - rjpeq +1;  //
}

double fsej(double ej)
{
    return ej -1 + rjpeq;
}

// Forward north-polar stereographic, R=0.5
void fslamphixy(double lon_deg, double lat_deg, double *x_out, double *y_out)
{
    const double phi  = lat_deg * rad;
    const double dlam = (lon_deg - rlampo) * rad;
    const double rho  = tan(M_PI_4 - 0.5 * phi);
    const double x    =  rho * sin(dlam);
    const double y    = -rho * cos(dlam);
    if (x_out) *x_out = x;
    if (y_out) *y_out = y;
}


// Inverse stereographic (x,y) to (lon,lat) in degrees
void fsxylamphi(double px, double py, double *plam_deg, double *pphi_deg)
{
    const double lon_rad = atan2(px, -py) + (rlampo * rad);
    const double t2 = px*px + py*py;
    const double t  = sqrt(t2);
    const double c = 2.0 * atan(t);
    const double lat_rad = (M_PI_2) - c;
    *pphi_deg = lat_rad / rad;
    *plam_deg = lon_rad / rad;
}


void fsfun_fast_with_deriv(double ej, double px, double py,
                           double *out_val, double *out_der)
{
    const double inv_j = 1.0 / rjpnthm1;

    const double k_theta = rdx * rad;

    const double x_fa  = ra2 * (ej - 1.0) * inv_j;
    const double x2_fa = x_fa * x_fa;
    const double x3_fa = x2_fa * x_fa;

    const double denom_fa = 3.0 * (x3_fa / (ra1 * ra2) + 1.0);

    const double fa =
        (atan(x_fa) + x3_fa / denom_fa) * rjpnthm1 / ra2;

    const double dx_fa     = ra2 * inv_j;
    const double dx3_fa    = 3.0 * x2_fa * dx_fa;
    const double ddenom_fa = 3.0 * dx3_fa / (ra1 * ra2);

    const double d_atan = 1.0 / (1.0 + x2_fa);
    const double d_frac =
        (dx3_fa * denom_fa - x3_fa * ddenom_fa) /
        (denom_fa * denom_fa);

    const double dfa =
        (d_atan + d_frac) * (rjpnthm1 / ra2) * dx_fa;

    const double theta_a  = k_theta * (fa + rjpeq - rjpeqt);
    const double dtheta_a = k_theta * dfa;

    const double aval = exp(-theta_a);
    const double da   = -aval * dtheta_a;

    const double t_fb  = (ej - 1.0) * inv_j;
    const double t2_fb = t_fb * t_fb;
    const double t3_fb = t2_fb * t_fb;

    const double u_fb    = rb2 * t3_fb;
    const double tanh_fb = tanh(u_fb);
    const double fsfb    = 1.0 + rb1 * tanh_fb;

    /* derivative of fsfb */
    const double dt_fb    = inv_j;
    const double dt3_fb   = 3.0 * t2_fb * dt_fb;
    const double sech2_fb = 1.0 - tanh_fb * tanh_fb;

    const double dfsfb = rb1 * sech2_fb * (rb2 * dt3_fb);

    // Use asin(tanh(theta)) = pi/2 - 2*atan(exp(-theta))

    const double theta_b  = k_theta * (ej - 1.0 + rjpeq - rjpeqt);
    const double dtheta_b = k_theta;

    const double em_b = exp(-theta_b);

    const double fsy0p_b =
        (M_PI_2 - 2.0 * atan(em_b)) / rad;

    const double dfsy0p_b =
        (2.0 * em_b * dtheta_b) /
        ((1.0 + em_b * em_b) * rad);


    const double fsphib  = fsfb * fsy0p_b;
    const double dfsphib = dfsfb * fsy0p_b + fsfb * dfsy0p_b;


    const double arg_b = M_PI_4 - 0.5 * rad * fsphib;
    const double bval  = tan(arg_b);

    const double darg_b = -0.5 * rad * dfsphib;
    const double db     = (1.0 + bval * bval) * darg_b;

    const double fsphiy_y0 = em_b;
    const double y0        = (fsphiy_y0 - aval) * ry0;

    const double d_fsphiy_y0 = -em_b * dtheta_b;
    const double dy0         = (d_fsphiy_y0 - da) * ry0;

    const double a2 = aval * aval;
    const double b2 = bval * bval;

    const double dy = py - y0;

    const double F = a2*(px*px - b2) + b2*dy*dy;

    const double dF =
        2.0 * aval * da * (px * px - b2) +
        2.0 * bval * db * (dy * dy - a2) -
        2.0 * b2 * dy * dy0;

    *out_val = F;
    *out_der = dF;
}


void fsfun_fast_x0_with_deriv(double ej, double py,
                              double *out_val, double *out_der)
{
    const double inv_j = 1.0 / rjpnthm1;
    const double k_theta = rdx * rad;

    const double x  = ra2 * (ej - 1.0) * inv_j;
    const double x2 = x * x;
    const double x3 = x2 * x;

    const double denom = 3.0 * (x3 / (ra1 * ra2) + 1.0);

    const double fa =
        (atan(x) + x3 / denom) * rjpnthm1 / ra2;

    const double dx   = ra2 * inv_j;
    const double dx3  = 3.0 * x2 * dx;
    const double dden = 3.0 * dx3 / (ra1 * ra2);

    const double d_atan = 1.0 / (1.0 + x2);

    const double d_frac =
        (dx3 * denom - x3 * dden) /
        (denom * denom);

    const double dfa =
        (d_atan + d_frac) * (rjpnthm1 / ra2) * dx;

    const double theta_a  = k_theta * (fa + rjpeq - rjpeqt);
    const double dtheta_a = k_theta * dfa;

    const double aval = exp(-theta_a);
    const double da   = -aval * dtheta_a;

    const double theta_b =
        k_theta * (ej - 1.0 + rjpeq - rjpeqt);

    const double fsphiy_y0 = exp(-theta_b);

    const double y0 = (fsphiy_y0 - aval) * ry0;

    const double d_fsphiy_y0 =
        -fsphiy_y0 * k_theta;

    const double dy0 =
        (d_fsphiy_y0 - da) * ry0;

    const double dy = py - y0;

    const double F =
        dy * dy - aval * aval;

    const double dF =
        -2.0 * dy * dy0
        -2.0 * aval * da;

    *out_val = F;
    *out_der = dF;
}
