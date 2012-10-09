#ifndef __CFAC_H_
#define __CFAC_H_

typedef struct _cfac_t cfac_t;

/* cfac.c */
cfac_t *
cfac_new(void);

void 
cfac_free(cfac_t *cfac);

/* nucleus.c */
int SetAtom(cfac_t *cfac, char *s, double z, double mass, double rn);
double GetAtomicNumber(const cfac_t *cfac);
double GetAtomicMass(const cfac_t *cfac);
double GetAtomicR(const cfac_t *cfac);
const char *GetAtomicSymbol(const cfac_t *cfac);
double GetAtomicEffectiveZ(const cfac_t *cfac, double r);

#endif /* __CFAC_H_ */
