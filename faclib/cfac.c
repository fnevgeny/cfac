#include <stdlib.h>
#include <string.h>

#include "cfacP.h"

static int cfac_init_radial(cfac_t *cfac)
{
    int i, ndim, blocks[5];

    cfac->n_awgrid = 1;
    cfac->awgrid[0]= EPS3;
    
    cfac->optimize_control.stabilizer = OPTSTABLE;
    cfac->optimize_control.tolerance  = OPTTOL;
    cfac->optimize_control.maxiter    = OPTNITER;
    cfac->optimize_control.screened_charge = 1.0; 
    cfac->optimize_control.screened_kl     = 1;
    cfac->optimize_control.iprint          = OPTPRINT;

    SetSlaterCut(cfac, -1, -1);

    cfac->qed.se  = QEDSE;
    cfac->qed.vp  = QEDVP;
    cfac->qed.nms = QEDNMS;
    cfac->qed.sms = QEDSMS;
    cfac->qed.br  = QEDBREIT;
    
    cfac->orbitals = malloc(sizeof(ARRAY));
    if (!cfac->orbitals) {
        return 1;
    }
    if (ArrayInit(cfac->orbitals, sizeof(ORBITAL), ORBITALS_BLOCK,
      FreeOrbitalData, InitOrbitalData) < 0) {
        return 1;
    }

    cfac->slater_array    = malloc(sizeof(MULTI));
    cfac->breit_array     = malloc(sizeof(MULTI));
    cfac->residual_array  = malloc(sizeof(MULTI));
    cfac->vinti_array     = malloc(sizeof(MULTI));
    cfac->qed1e_array     = malloc(sizeof(MULTI));
    cfac->multipole_array = malloc(sizeof(MULTI));
    cfac->moments_array   = malloc(sizeof(MULTI));
    cfac->gos_array       = malloc(sizeof(MULTI));
    cfac->yk_array        = malloc(sizeof(MULTI));
    if (!cfac->slater_array    ||
        !cfac->breit_array     ||
        !cfac->residual_array  ||
        !cfac->vinti_array     ||
        !cfac->qed1e_array     || 
        !cfac->multipole_array ||
        !cfac->moments_array   ||
        !cfac->gos_array       ||   
        !cfac->yk_array) {
        return 1;
    }

    ndim = 5;
    for (i = 0; i < ndim; i++) {
        blocks[i] = MULTI_BLOCK6;
    }
    MultiInit(cfac->slater_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);

    ndim = 5;
    for (i = 0; i < ndim; i++) {
        blocks[i] = MULTI_BLOCK5;
    }
    MultiInit(cfac->breit_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);

    ndim = 2;
    for (i = 0; i < ndim; i++) {
        blocks[i] = MULTI_BLOCK2;
    }
    MultiInit(cfac->residual_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);
    MultiInit(cfac->vinti_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);
    MultiInit(cfac->qed1e_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);

    ndim = 3;
    for (i = 0; i < ndim; i++) {
        blocks[i] = MULTI_BLOCK4;
    }
    MultiInit(cfac->multipole_array,
        sizeof(double *), ndim, blocks, FreeMultipole, InitPointerData);

    ndim = 3;
    for (i = 0; i < ndim; i++) {
        blocks[i] = MULTI_BLOCK3;
    }
    MultiInit(cfac->moments_array,
        sizeof(double), ndim, blocks, NULL, InitDoubleData);
    MultiInit(cfac->gos_array,
        sizeof(double *), ndim, blocks, FreeMultipole, InitPointerData);
    MultiInit(cfac->yk_array,
        sizeof(SLATER_YK), ndim, blocks, FreeYkData, InitYkData);

    return 0;
}

static void cfac_free_radial(cfac_t *cfac)
{
    AVERAGE_CONFIG *acfg = &cfac->average_config;
    
    ArrayFree(cfac->orbitals);
    free(cfac->orbitals);
    
    MultiFree(cfac->slater_array);
    free(cfac->slater_array);
    
    MultiFree(cfac->breit_array);
    free(cfac->breit_array);
    
    MultiFree(cfac->residual_array);
    free(cfac->residual_array);
    
    MultiFree(cfac->vinti_array);
    free(cfac->vinti_array);
    
    MultiFree(cfac->qed1e_array);
    free(cfac->qed1e_array);
    
    MultiFree(cfac->multipole_array);
    free(cfac->multipole_array);
    
    MultiFree(cfac->moments_array);
    free(cfac->moments_array);
    
    MultiFree(cfac->gos_array);
    free(cfac->gos_array);
    
    MultiFree(cfac->yk_array);
    free(cfac->yk_array);
    
    if (acfg->n_shells > 0) {      
        free(acfg->n);
        free(acfg->kappa);
        free(acfg->nq);
    }
}

cfac_t *cfac_new(void)
{
    cfac_t *cfac;
    unsigned int i;
    
    cfac = malloc(sizeof(cfac_t));
    if (!cfac) {
        return NULL;
    }
    memset(cfac, 0, sizeof(cfac_t));

    cfac->confint = 0;
    
    cfac->angz_maxn = 0;
    cfac->angz_cut = ANGZCUT;
    cfac->mix_cut = MIXCUT;
    cfac->mix_cut2 = MIXCUT2;

    cfac->sym_pp = -1;
    cfac->sym_njj = 0;
    cfac->sym_jj = NULL;

    /* init config groups */
    cfac->n_groups = 0;
    cfac->cfg_groups = malloc(MAX_GROUPS*sizeof(CONFIG_GROUP));
    if (!cfac->cfg_groups) {
        cfac_free(cfac);
        return NULL;
    }
    memset(cfac->cfg_groups, 0, MAX_GROUPS*sizeof(CONFIG_GROUP));
    for (i = 0; i < MAX_GROUPS; i++) {
        CONFIG_GROUP *cfg = &cfac->cfg_groups[i];
        strcpy(cfg->name, "_all_");
        ArrayInit(&cfg->cfg_list,
            sizeof(CONFIG), CONFIGS_BLOCK, FreeConfigData, InitConfigData);
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
    
    cfac->potential = malloc(sizeof(POTENTIAL));
    if (!cfac->potential) {
        cfac_free(cfac);
        return NULL;
    }
    memset(cfac->potential, 0, sizeof(POTENTIAL));
    SetBoundary(cfac, 0, 1.0, -1.0);
    SetRadialGrid(cfac, DMAXRP, -1.0, -1.0, -1.0);

    if (cfac_init_radial(cfac)) {
        cfac_free(cfac);
        return NULL;
    }
    
    if (cfac_init_recouple(cfac)) {
        cfac_free(cfac);
        return NULL;
    }
    
    if (cfac_init_coulomb(cfac)) {
        cfac_free(cfac);
        return NULL;
    }
    
    cfac->nhams = 0;
    cfac->hams = malloc(MAX_HAMS*sizeof(SHAMILTON));
    if (!cfac->hams) {
        cfac_free(cfac);
        return NULL;
    }
    memset(cfac->hams, 0, MAX_HAMS*sizeof(SHAMILTON));


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


    cfac->angz_array = malloc(sizeof(ANGZ_DATUM)*MAX_HAMS2);
    if (!cfac->angz_array) {
        cfac_free(cfac);
        return NULL;
    }
    memset(cfac->angz_array, 0, sizeof(ANGZ_DATUM)*MAX_HAMS2);
    
    cfac->angzxz_array = malloc(sizeof(ANGZ_DATUM)*MAX_HAMS2);
    if (!cfac->angzxz_array) {
        cfac_free(cfac);
        return NULL;
    }
    memset(cfac->angzxz_array, 0, sizeof(ANGZ_DATUM)*MAX_HAMS2);
    
    cfac->transition_options.gauge = DGAUGE;
    cfac->transition_options.mode  = DMODE;
    cfac->transition_options.max_e = ERANK;
    cfac->transition_options.max_m = MRANK;
    cfac->transition_options.fr_interpolate = 1;

    return cfac;
}


void cfac_free(cfac_t *cfac)
{
    if (cfac) {
        unsigned int i;

        for (i = 0; i < cfac->n_groups; i++) {
            CONFIG_GROUP *cfg = &cfac->cfg_groups[i];
            ArrayFree(&cfg->cfg_list);
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
        
        free(cfac->potential);
        
        cfac_free_radial(cfac);
        
        ArrayFree(cfac->levels);
        free(cfac->levels);
        ArrayFree(cfac->eblevels);
        free(cfac->eblevels);

        ArrayFree(cfac->ecorrections);
        free(cfac->ecorrections);

        cfac_free_recouple(cfac);
        
        cfac_free_coulomb(cfac);
        
        FreeHamsArray(cfac);
        free(cfac->hams);
        
        if (cfac->angz_array) {
            for (i = 0; i < MAX_HAMS2; i++) {
                FreeAngZDatum(&(cfac->angz_array[i]));
            }
            free(cfac->angz_array);
        }
        
        if (cfac->angzxz_array) {
            for (i = 0; i < MAX_HAMS2; i++) {
                FreeAngZDatum(&(cfac->angzxz_array[i]));
            }
            free(cfac->angzxz_array);
        }
        
        if (cfac->angmz_array) {
            for (i = 0; i < MAX_HAMS2; i++) {
                FreeAngZDatum(&(cfac->angmz_array[i]));
            }
            free(cfac->angmz_array);
        }
        
        if (cfac->sym_jj) {
            free(cfac->sym_jj);
        }
        
        free(cfac);
    }
}
