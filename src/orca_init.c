#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "orca_init.h"
#include "errors.h"
#include "chebcheb_data.h"
#include "pchip_data.h"
#include "xy_to_ej.h"

static int orca_initialized = 0;

#ifndef ORCA05_DATA_PATH
#define ORCA05_DATA_PATH "data"
#endif

#define ORCA_PATH_MAX 1024
static char cheb_path_buf[ORCA_PATH_MAX];
static char pchip_path_buf[ORCA_PATH_MAX];

static const char *cheb_path = NULL;
static const char *pchip_path = NULL;

static void build_paths_from_dir(const char *dir)
{
    snprintf(cheb_path_buf,  sizeof(cheb_path_buf),  "%s/chebcheb_data.bin", dir);
    snprintf(pchip_path_buf, sizeof(pchip_path_buf), "%s/pchip_data.bin",    dir);

    cheb_path  = cheb_path_buf;
    pchip_path = pchip_path_buf;
}

void orca_set_data_dir(const char *dir)
{
    if (!dir || !dir[0])
        return;
    build_paths_from_dir(dir);
}

static void ensure_paths_initialized(void)
{
    if (cheb_path && pchip_path)
        return;

    const char *dir = getenv("ORCA05_DATA_PATH");
    if (!dir || !dir[0])
        dir = ORCA05_DATA_PATH;

    build_paths_from_dir(dir);
}


int orca_init(void)
{
    if (orca_initialized)
        return ORCA_OK;

    ensure_paths_initialized();

    if (chebcheb_load_coeffs(cheb_path) != 0) {
        fprintf(stderr, "orca_init: failed loading Cheb coefficients (%s)\n", cheb_path);
        return ORCA_INIT_FAILED;
    }

    if (pchip_load_coeffs(pchip_path) != 0) {
        fprintf(stderr, "orca_init: failed loading PCHIP coefficients (%s)\n", pchip_path);
        chebcheb_unload_coeffs();   // rollback
        return ORCA_INIT_FAILED;
    }

    chebcheb_init_vertical_edges();
    fill_aby0_cache();

    orca_initialized = 1;
    return ORCA_OK;
}
