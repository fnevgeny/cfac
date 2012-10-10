#ifndef __CFACP_H_
#define __CFACP_H_

#include "cfac.h"
#include "array.h"
#include "config.h"

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

    SYMMETRY *symmetry_list;  /* list of symmetries. The i-th symmetry has
                                 j = floor(i/2) and parity = mod(i, 2).      */

    ARRAY *levels_per_ion;
};

/* config.c */
void
FreeConfigData(void *p);

#endif /* __CFACP_H_ */
