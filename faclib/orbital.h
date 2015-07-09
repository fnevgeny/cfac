/*
 *   FAC - Flexible Atomic Code
 *   Copyright (C) 2001-2015 Ming Feng Gu
 *   Portions Copyright (C) 2010-2015 Evgeny Stambulchik
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 * 
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 * 
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ORBITAL_H_
#define _ORBITAL_H_

#include "consts.h"

typedef struct _POTENTIAL_ {
  int flag;               /* radial grid completeness             */
  
  unsigned int anum;      /* nuclear charge (cfac->nucleus.anum)  */
  
  int maxrp;              /* used length of the [MAXRP] arrays    */
  double rmin;            /* starting point of the radial mesh    */
  double ratio;           /* incr. ratio of radial mesh near 0    */
  double asymp;           /* number of mesh points per oscillation
                             wavelength for high-n orbitals       */
  
  double Navg;            /* total number of electrons in the
                             average config (cfac->acfg)          */

  double ar, br;          /* parameters for the transformation
                             \rho = ar*sqrt(r) + br*log(r)        */

  int ib, nb, ib1; 
  
  int r_core;             /* core radius index                    */
  unsigned int nmax;      /* above this PQN, Coulomb WF for
                             r > rad[r_core] are assumed          */
  
  double rad[MAXRP];      /* radial grid; rad[0] = rmin/at.number */

  double dr_drho[MAXRP];  /* dr/d\rho                             */
  double dr_drho2[MAXRP]; /* square root of the above             */
  
  double Vn[MAXRP];        /* nucleus potential                   */
  
  double lambda, a;       /* optimization parameters for Vc       */
  double Vc[MAXRP];       /* optimized central potential          */
  double dVc[MAXRP];      /* its first derivative d(Vc)/dr        */
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
void Differential(double *p, double *dp, int i1, int i2);
int SetOrbitalRGrid(POTENTIAL *pot);
int SetPotentialZ(cfac_t *cfac);
int SetPotentialUehling(cfac_t *cfac, int vp);
int SetPotentialVc(POTENTIAL *pot);
int SetPotentialU(POTENTIAL *pot, int n);
int SetPotentialW (POTENTIAL *pot, double e, int kappa);
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot);

#endif
