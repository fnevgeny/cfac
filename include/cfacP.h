#ifndef __CFACP_H_
#define __CFACP_H_

#include "cfac.h"
#include "array.h"
#include "config.h"
#include "orbital.h"
#include "radial.h"
#include "structure.h"

typedef struct {
  char symbol[5];
  unsigned int anum;
  double mass;
  double rn;                  /* effective radius of the nucleus             */
} cfac_nucleus_t;

struct _cfac_t {
    cfac_nucleus_t nucleus;

    CONFIG_GROUP *cfg_groups; /* a list of configuration groups              */
    int n_groups;             /* number of configuration groups present      */

    double ef, bf, eb_angle;  /* electric, magnetic ield, and angle between  */
    double e1[3];             /* spherical components of the E field         */
    double b0, b1[3], b2[5];  /* 0th, 1st, and 2nd-order tensors of B        */

    SYMMETRY *symmetry_list;  /* list of symmetries. The i-th symmetry has
                                 j = floor(i/2) and parity = mod(i, 2).      */

    ARRAY *levels_per_ion;

    ARRAY *orbitals;          /* array of orbitals                           */
    int n_orbitals;           /* total number of orbitals                    */
    int n_continua;           /* number of continuum orbitals                */
 
    AVERAGE_CONFIG acfg;      /* average config for potential optimization   */

    POTENTIAL *potential;     /* potential                                   */


    MULTI *slater_array;
    MULTI *breit_array;
    MULTI *vinti_array;
    MULTI *qed1e_array;
    MULTI *residual_array;
    MULTI *multipole_array; 
    MULTI *moments_array;
    MULTI *gos_array;
    MULTI *yk_array;

    struct {
      double stabilizer;
      double tolerance;       /* tolerance for self-consistency              */
      int maxiter;            /* max iter. for self-consistency              */
      double screened_charge; 
      int screened_kl;
      int n_screen;
      int *screened_n;
      int iprint;             /* printing information in each iteration.     */
      int iset;
    } optimize_control;
    
    struct {
      int kl0;
      int kl1;
    } slater_cut;

    struct {
      int se;
      int vp;
      int nms;
      int sms;
      int br;
    } qed;

    int n_awgrid;
    double awgrid[MAXNTE];

    struct {
        int max_rank;         /* the maximum rank of the operators allowed.  */
        MULTI *int_shells;    /* interacting shell information array         */
    } recouple;

    struct {
        int n_h;
        int kl_h;
        int n_h_max;
        int kl_h_max;
        ARRAY *dipole_array;
    } coulomb;
    
    SHAMILTON *hams;          /* symmetry Hamiltonians                       */
    int nhams;                /* number of them in use                       */

    ARRAY *levels;            /* levels                                      */
    int n_levels;             /* number of levels                            */
    
    ARRAY *eblevels;          /* levels when calculated with fields          */
    int n_eblevels;           /* number of eblevels                          */

    ARRAY *ecorrections;      /* energy corrections                          */
    int ncorrections;         /* number of energy corrections                */

    int confint;              /* configuration interaction flag (-1,0,1,2,3) */

    int angz_maxn;            /* max PQN above which angular mix between
                                 different bound configurations is ignored   */
    double angz_cut;          /* threshold of angular mixing                 */
    
    double mix_cut;           /* threshold mixing (relative to the leading
                                 component); bases with weaker mixings are
                                 not recoupled                               */
    
    double mix_cut2;          /* _lower_ threshold for recoupling between
                                 different configurations                    */

    int sym_pp;               /* symmetry parity                             */
    int *sym_jj;              /* sorted array of user-defined 2*J symmetries */
    int sym_njj;              /* length of the above array                   */


    ANGZ_DATUM *angz_array;   /* angular coefficients                        */
    ANGZ_DATUM *angzxz_array; /* ZxZ angular coefficients                    */
    ANGZ_DATUM *angmz_array;  /* precalculated angular coefficients          */

    ANGULAR_FROZEN ang_frozen;/* angular coefficients for frozen states      */

    struct {
        int gauge;            /* gauge (Coulomb/Babushkin)                   */
        int mode;             /* mode (relativistic/non-relativistic)        */
        int max_e;            /* maximum rank of electric multipoles         */
        int max_m;            /* maximum rank of magnetic multipoles         */
        int fr_interpolate;   /* use interpolation for FR m. elements        */
    } tr_opts;
};

typedef struct {
    int cbindex[CBMULT][CBMULT+1];
    double *cb[MAXNE][MAXNTE][MAXNE][MAXNCB];
} cfac_cbcache_t;

/* config.c */
void
FreeConfigData(void *p);
void
InitConfigData(void *p, int n);

/* recouple */
int
cfac_init_recouple(cfac_t *cfac);
void
cfac_free_recouple(cfac_t *cfac);

/* coulomb.c */
int
cfac_init_coulomb(cfac_t *cfac);
void
cfac_free_coulomb(cfac_t *cfac);
void
cfac_cbcache_init(cfac_cbcache_t *cbcache);
void
cfac_cbcache_free(cfac_cbcache_t *cbcache);

/* radial.c */
void
FreeOrbitalData(void *p);
void
InitOrbitalData(void *p, int n);
void
FreeMultipole(void *p);
void
FreeYkData(void *p);
void
InitYkData(void *p, int n);

/* structure.c */
void
FreeLevelData(void *p);
void
FreeHamsArray(cfac_t *cfac);
void
FreeAngZDatum(ANGZ_DATUM *ap);
void
InitLevelData(void *p, int n);
void
cfac_hamiltonian_free(HAMILTON *h);

#endif /* __CFACP_H_ */
