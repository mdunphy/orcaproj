#ifndef ORCA_ERRORS_H
#define ORCA_ERRORS_H

typedef enum {
    ORCA_OK = 0,       /* Success                                   */

    //  User/API input errors
    ORCA_NULL_PTR,     /* NULL pointer passed by caller             */
    ORCA_INVALID_ARG,  /* inconsistent or invalid argument set      */
    ORCA_INIT_FAILED,  /* orca_init() was not successful            */
    ORCA_NAN_INPUT,    /* NaN input coordinate                      */
    ORCA_INF_INPUT,    /* Inf input coordinate                      */
    ORCA_BAD_RANGE,    /* lat/lon outside supported range           */
    ORCA_DOMAIN,       /* valid numeric input, but outside map      */

    // Root solver conditions
    ORCA_ROOT_NO_BRACKET,  /* function does not change sign         */
    ORCA_ROOT_MAXITER,     /* maximum iterations exceeded           */

    // Failures of the internal mapping machinery
    ORCA_EINTERNAL          /* internal failure; not caller’s fault */

} orca_error_t;

#endif /* ORCA_ERRORS_H */
