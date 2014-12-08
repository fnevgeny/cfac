#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include "consts.h"

typedef struct _POTENTIAL_ {
  int mode;
  
  int flag;               /* radial grid completeness             */
  
  int maxrp;              /* used length of the [MAXRP] arrays    */
  double ratio;           /* incr. ratio of radial mesh near 0    */
  double asymp;           /* number of mesh points per oscillation
                             wavelength for high-n orbitals       */
  double rmin;            /* starting point of the radial mesh    */
  
  double hxs;
  double N;               /* effective number of electrons        */
  double lambda, a;       /* optimization parameters for Vc       */
  double ar, br;          /* parameters for the transformation    */
  int ib, nb, ib1; 
  double bqp;             /* boundary condition                   */
  
  int r_core;             /* core radius index                    */
  unsigned int nmax;      /* above this PQN, Coulomb WF for
                             r > rad[r_core] are assumed          */
  
  double rad[MAXRP];      /* radial grid; rad[0] = rmin/at.number */

  double Z[MAXRP];        /* nuclear charge distribution          */
  
  double dr_drho[MAXRP];  /* dr/d\rho                             */
  double dr_drho2[MAXRP];
  
  double Vc[MAXRP];       /* optimized potential                  */
  double dVc[MAXRP];      /* its first derivative                 */
  double dVc2[MAXRP];     /* its second derivative                */
  
  double U[MAXRP];        /* direct interaction part of potential */
  double dU[MAXRP];       /* its first derivative                 */
  double dU2[MAXRP];      /* its second derivative                */
  
  double W[MAXRP];
  
  double uehling[MAXRP];  /* the Uehling potential                */
  
  double veff[MAXRP];
} POTENTIAL;

typedef struct _ORBITAL_ {
  int n;          /* PQN; 0 => free; <0 => BasisOuter                */
  int kappa;      /* relativistic angular QN (l - j)*(2*j + 1)       */
  double energy;
  double qr_norm;
  double phase;
  double *wfun;   /* radial wave-function, wfun[0] ... wfun[ilast-1] */
  int ilast;
} ORBITAL;

int GetNMax(const POTENTIAL *pot);
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa);
int RadialSolver(const cfac_t *cfac, ORBITAL *orb);
int RadialBasis(ORBITAL *orb, POTENTIAL *pot);
int RadialRydberg(ORBITAL *orb, POTENTIAL *pot);
int RadialBound(ORBITAL *orb, POTENTIAL *pot);
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot);
int RadialFree(ORBITAL *orb, POTENTIAL *pot);
double InnerProduct(int i1, int n, 
		    double *p1, double *p2, POTENTIAL *pot);
void Differential(double *p, double *dp, int i1, int i2);
int SetOrbitalRGrid(cfac_t *cfac);
double GetRFromRho(double rho, double a, double b, double r0);
int SetPotentialZ(cfac_t *cfac, double c);
int SetPotentialUehling(cfac_t *cfac, int vp);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n, double *u);
int SetPotentialW (POTENTIAL *pot, double e, int kappa);
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot);

#endif
