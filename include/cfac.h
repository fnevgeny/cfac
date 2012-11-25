#ifndef __CFAC_H_
#define __CFAC_H_

typedef struct _cfac_t cfac_t;

/* cfac.c */
cfac_t *
cfac_new(void);

void 
cfac_free(cfac_t *cfac);

/* nucleus.c */
int
cfac_set_atom(cfac_t *cfac, char *s, double z, double mass, double rn);
double
cfac_get_atomic_number(const cfac_t *cfac);
double
cfac_get_atomic_mass(const cfac_t *cfac);
double
cfac_get_atomic_rn(const cfac_t *cfac);
const char *
cfac_get_atomic_symbol(const cfac_t *cfac);
double
cfac_get_atomic_effective_z(const cfac_t *cfac, double r);

#endif /* __CFAC_H_ */
