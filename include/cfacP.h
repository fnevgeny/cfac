#ifndef __CFACP_H_
#define __CFACP_H_

#include "cfac.h"

typedef struct {
  char symbol[5];
  double atomic_number;
  double mass;
  double rn;
} cfac_nucleus_t;

struct _cfac_t {
    cfac_nucleus_t nucleus;
};


#endif /* __CFACP_H_ */
