#include <stdlib.h>
#include <string.h>

#include "cfacP.h"

cfac_t *cfac_new(void)
{
    cfac_t *cfac;
    
    cfac = malloc(sizeof(cfac_t));
    if (cfac) {
        memset(cfac, 0, sizeof(cfac_t));
    }
    
    return cfac;
}

void cfac_free(cfac_t *cfac)
{
    if (cfac) {
        unsigned int k;

        for (k = 0; k <= cfac->anum; k++) {
            ArrayFree(&cfac->levels_per_ion[k]);
        }
        free(cfac->levels_per_ion);

        free(cfac);
    }
}

