#ifndef _TRANSITION_H_
#define _TRANSITION_H_

typedef struct _TRANSITION_ {
    int nup;
    int nlo;
    LEVEL *lup;
    LEVEL *llo;
    double e;
} TRANSITION;

int GetTransition(const cfac_t *cfac,
    int nlo, int nup, TRANSITION *tr, int *swapped);
int SetTransitionCut(cfac_t *cfac, double c);
double GetTransitionCut(const cfac_t *cfac);
void SetTransitionMode(cfac_t *cfac, int m);
void SetTransitionGauge(cfac_t *cfac, int m);
void SetTransitionMaxE(cfac_t *cfac, int m);
void SetTransitionMaxM(cfac_t *cfac, int m);
void SetTransitionOptions(cfac_t *cfac, int gauge, int mode, int max_e, int max_m);
int GetTransitionGauge(const cfac_t *cfac);
int GetTransitionMode(const cfac_t *cfac);
int TRMultipole(cfac_t *cfac, double *strength, double *energy,
		int m, int low, int up);
int OverlapLowUp(int nlow, int *low, int nup, int *up);
int SaveTransition(cfac_t *cfac, int nlow, int *low, int nup, int *up,
		   char *fn, int multipole);
int SaveTransitionEB(cfac_t *cfac, int nlow, int *low, int nup, int *up,
		     char *fn, int multipole);
int GetLowUpEB(cfac_t *cfac, int *nlow, int **low, int *nup, int **up, 
	       int nlow0, int *low0, int nup0, int *up0);
int PolarizeCoeff(char *ifn, char *ofn, int i0, int i1);

int TRMultipoleEB(cfac_t *cfac, double *strength, double *energy, int m,
                  int lower, int upper);

#endif
