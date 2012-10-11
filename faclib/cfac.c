#include <stdlib.h>
#include <string.h>

#include "cfacP.h"

cfac_t *cfac_new(void)
{
    cfac_t *cfac;
    unsigned int i;
    
    cfac = malloc(sizeof(cfac_t));
    if (!cfac) {
        return NULL;
    }
    memset(cfac, 0, sizeof(cfac_t));

    /* init config groups */
    cfac->n_groups = 0;
    cfac->cfg_groups = malloc(MAX_GROUPS*sizeof(CONFIG_GROUP));
    if (!cfac->cfg_groups) {
        cfac_free(cfac);
        return NULL;
    }
    for (i = 0; i < MAX_GROUPS; i++) {
        strcpy(cfac->cfg_groups[i].name, "_all_");
        cfac->cfg_groups[i].n_cfgs = 0;
        ArrayInit(&(cfac->cfg_groups[i].cfg_list),
            sizeof(CONFIG), CONFIGS_BLOCK, FreeConfigData, NULL);
    }

    /* init config symmetries */
    cfac->symmetry_list = malloc(MAX_SYMMETRIES*sizeof(SYMMETRY));
    if (!cfac->symmetry_list) {
        cfac_free(cfac);
        return NULL;
    }
    for (i = 0; i < MAX_SYMMETRIES; i++) {
        cfac->symmetry_list[i].n_states = 0;
        ArrayInit(&(cfac->symmetry_list[i].states),
            sizeof(STATE), STATES_BLOCK, NULL, NULL);
    }


    cfac->levels = malloc(sizeof(ARRAY));
    if (!cfac->levels) {
        cfac_free(cfac);
        return NULL;
    }
    ArrayInit(cfac->levels, sizeof(LEVEL), LEVELS_BLOCK,
        FreeLevelData, InitLevelData);

    cfac->eblevels = malloc(sizeof(ARRAY));
    if (!cfac->eblevels) {
        cfac_free(cfac);
        return NULL;
    }
    ArrayInit(cfac->eblevels, sizeof(LEVEL), LEVELS_BLOCK,
        FreeLevelData, InitLevelData);

    cfac->ecorrections = malloc(sizeof(ARRAY));
    if (!cfac->ecorrections) {
        cfac_free(cfac);
        return NULL;
    }
    ArrayInit(cfac->ecorrections, sizeof(ECORRECTION), 512, NULL, NULL);
    
    return cfac;
}


void cfac_free(cfac_t *cfac)
{
    if (cfac) {
        unsigned int i;

        for (i = 0; i < cfac->n_groups; i++) {
            ArrayFree(&(cfac->cfg_groups[i].cfg_list));
        }
        free(cfac->cfg_groups);

        for (i = 0; i < MAX_SYMMETRIES; i++) {
            if (cfac->symmetry_list[i].n_states > 0) {
                ArrayFree(&(cfac->symmetry_list[i].states));
            }
        }
        free(cfac->symmetry_list);

        for (i = 0; i <= cfac->anum; i++) {
            ArrayFree(&cfac->levels_per_ion[i]);
        }
        free(cfac->levels_per_ion);
        
        ArrayFree(cfac->levels);
        free(cfac->levels);
        ArrayFree(cfac->eblevels);
        free(cfac->eblevels);

        ArrayFree(cfac->ecorrections);
        free(cfac->ecorrections);
        
        free(cfac);
    }
}
