#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

/*************************************************************
  Header for module "coulomb". 
  This module calculates quatities related to the H-like ions.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

void    SetHydrogenicNL(cfac_t *cfac, int n, int kl, int nm, int klm);
void    GetHydrogenicNL(const cfac_t *cfac, int *n, int *kl, int *nm, int *klm);

double  HydrogenicDipole(const cfac_t *cfac, double z, int n0, int kl0, 
			int n1, int kl1);
double HydrogenicExpectation(double z, int m, int n, int kl);
double HydrogenicSelfEnergy(double z, int n, int k);
double  CoulombPhaseShift(double z, double e, int kappa);
int CoulombMultip(char *fn, double z, double te, double e1,
		  int k, int q0, int q1, int m);
double *GetCoulombBethe(const cfac_cbcache_t *cbcache,
    int ie2, int ite, int ie1, int t, int q);
double  GetCoulombBetheAsymptotic(double te, double e1);
int     CoulombBetheTail(int n, double *w, int nkl, double *kl, double *tcb);
int     PrepCoulombBethe(cfac_cbcache_t *cbcache,
                         int ne2, int nte, int ne1, double z,
			 double *e2, double *te, double *e1,
			 int nkl, double *kl, int mode);
#endif
