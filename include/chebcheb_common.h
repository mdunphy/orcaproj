#ifndef CHEBCHEB_COMMON_H
#define CHEBCHEB_COMMON_H

// Maximum degree permitted (must be >= the max over Px, Py in the cheb tables)
#ifndef CHEBCHEB_MAX_DEG
#define CHEBCHEB_MAX_DEG 39
#endif

// Temporary buffer size, needs to be one larger than max degree
#ifndef CHEBCHEB_MAX_AX
#define CHEBCHEB_MAX_AX (CHEBCHEB_MAX_DEG + 1)
#endif

/* Maximum number of tiles; if using>64 we'll need to rework the rowcache bitmask */
#ifndef CHEBCHEB_MAX_MODELS
#define CHEBCHEB_MAX_MODELS 64
#endif

#endif /* CHEBCHEB_COMMON_H */
