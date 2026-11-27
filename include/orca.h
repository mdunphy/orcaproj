// Main interface to include if using this as a library
// The other files should not be included directly
// TODO: reorganize the header files into private/public
#ifndef ORCA_H
#define ORCA_H

#ifdef __cplusplus
extern "C" {
#endif

#include "orca_init.h"
#include "errors.h"
#include "forward.h"
#include "inverse.h"
#include "theta.h"

#ifdef __cplusplus
}
#endif

#endif /* ORCA_H */