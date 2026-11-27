#ifndef ORCA_TRIPOLAR_INTERNAL_ERRORS_H
#define ORCA_TRIPOLAR_INTERNAL_ERRORS_H

/* Internal error codes */

typedef enum {
    ORCA_INT_OK = 0,      /* Success                              */

    // Numerical evaluation issues
    ORCA_INT_NAN,         /* NaN encountered internally           */
    ORCA_INT_INF,         /* Inf encountered internally           */
    ORCA_INT_BAD_VALUE,   /* Non-finite or nonsensical result     */

    // Tile/model geometry issues
    ORCA_INT_NO_TILE_COVERAGE,  /* (i,j) falls in no available tile */

    // Root solver failures
    ORCA_INT_ROOT_NO_BRACKET,   /* f(a) and f(b) do not bracket      */
    ORCA_INT_ROOT_MAXITER,      /* maximum iterations reached        */
    ORCA_INT_ROOT_FAIL,         /* invalid function evaluations      */

    // Internal pointer/reference problems
    ORCA_INT_EFAULT             /* NULL pointer or corrupt context   */

} orca_int_error_t;

#endif /* ORCA_TRIPOLAR_INTERNAL_ERRORS_H */
