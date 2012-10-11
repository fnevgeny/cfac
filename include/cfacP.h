#ifndef __CFACP_H_
#define __CFACP_H_

#include "cfac.h"
#include "array.h"
#include "config.h"
#include "structure.h"

typedef struct {
  char symbol[5];
  double atomic_number;
  double mass;
  double rn;
} cfac_nucleus_t;

struct _cfac_t {
    unsigned int anum;
    cfac_nucleus_t nucleus;

    CONFIG_GROUP *cfg_groups; /* a list of configuration groups              */
    int n_groups;             /* number of groups present                    */

    double ef, bf, eb_angle;  /* electric, magnetic ield, and angle between  */
    double e1[3];             /* spherical components of the E field         */
    double b0, b1[3], b2[5];  /* 0th, 1st, and 2nd-order tensors of B        */


    SYMMETRY *symmetry_list;  /* list of symmetries. The i-th symmetry has
                                 j = floor(i/2) and parity = mod(i, 2).      */

    ARRAY *levels_per_ion;
    
    HAMILTON *hamiltonian;    /* Hamiltonian                                 */


    ARRAY *levels;            /* levels                                      */
    int n_levels;             /* number of levels                            */
    
    ARRAY *eblevels;          /* levels when calculated with fields          */
    int n_eblevels;           /* number of eblevels                          */

    ARRAY *ecorrections;      /* energy corrections                          */
    int ncorrections;         /* number of energy corrections                */
};

/* config.c */
void
FreeConfigData(void *p);

/* structure.c */
void
FreeLevelData(void *p);
void
InitLevelData(void *p, int n);

#endif /* __CFACP_H_ */
