#ifndef ORCA_ERROR_MAP_H
#define ORCA_ERROR_MAP_H

#include "errors.h"
#include "internal_errors.h"

// Map internal error codes to external error codes
// This is only used internally
static inline orca_error_t
orca_map_internal_error(orca_int_error_t e)
{
    switch (e) {

        case ORCA_INT_OK:
            return ORCA_OK;

        case ORCA_INT_ROOT_NO_BRACKET:
            return ORCA_ROOT_NO_BRACKET;

        case ORCA_INT_ROOT_MAXITER:
            return ORCA_ROOT_MAXITER;

        case ORCA_INT_NAN:
            return ORCA_NAN_INPUT;

        case ORCA_INT_INF:
            return ORCA_INF_INPUT;

        case ORCA_INT_NO_TILE_COVERAGE:
            return ORCA_DOMAIN;

        case ORCA_INT_BAD_VALUE:
        case ORCA_INT_EFAULT:
        case ORCA_INT_ROOT_FAIL:
        default:
            return ORCA_EINTERNAL;
    }
}

#endif /* ORCA_ERROR_MAP_H */
