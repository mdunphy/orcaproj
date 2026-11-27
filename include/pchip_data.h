#ifndef PCHIP_DATA_H
#define PCHIP_DATA_H
#include <stdint.h>

typedef struct {
    int ny, nx;              /* tile shape (rows, cols)    */
    double jmin, jmax;       /* vertical physical range    */
    double imin, imax;       /* horizontal physical range  */
    uint64_t off;            /* offset (in doubles)        */
    uint64_t cnt;            /* number of coefficients     */
} pchip_tile_t;

extern const uint64_t pchip_total_doubles;
extern const pchip_tile_t pchip_tiles[];
extern const int pchip_num_tiles;
extern const double *pchip_coeffs;

int  pchip_load_coeffs(const char *path);
void pchip_unload_coeffs(void);

#endif /* PCHIP_DATA_H */
