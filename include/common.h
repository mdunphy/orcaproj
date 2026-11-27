#ifndef ORCA_TRIPOLAR_COMMON_H
#define ORCA_TRIPOLAR_COMMON_H

#include <math.h>
#include <assert.h>

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4  0.78539816339744830962
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

static inline double wraplon(double lon)
{
    if (!isfinite(lon))
        return lon; // unclear what to do here
    double w = fmod(lon + 180.0, 360.0);
    if (w < 0.0) w += 360.0;
    return w - 180.0;
}

#endif /* ORCA_TRIPOLAR_COMMON_H */
