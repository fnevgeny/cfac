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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>

#include "cfacP.h"
#include "coulomb.h"
#include "interpolation.h"
#include "cf77.h"

static const int max_iteration = 512;
static const double wave_zero = 1E-10;

static int SetVEffective(int kl, POTENTIAL *pot);
static int TurningPoints(int n, double e, const POTENTIAL *pot);
static int CountNodes(double *p, int i1, int i2);
static int IntegrateRadial(double *p, double e, const POTENTIAL *pot, 
			   int i1, double p1, int i2, double p2, int q);
static double Amplitude(double *p, double e, int kl, POTENTIAL *pot, int i1);
static int Phase(double *p, POTENTIAL *pot, int i1, double p0);
static int DiracSmall(ORBITAL *orb, POTENTIAL *pot, int i2);

int GetNMax(const POTENTIAL *pot) {
  return pot->nmax;
}
 
static double EnergyH(double z, double n, int ka)
{
  double a, np;
  a = FINE_STRUCTURE_CONST*z;
  
  np = n + sqrt(ka*ka - a*a) - abs(ka);

  return (1.0/sqrt(1.0 + a*a/(np*np)) - 1.0)/FINE_STRUCTURE_CONST2;
}
    
double RadialDiracCoulomb(int npts, double *p, double *q, double *r,
			  double z, int n, int kappa) {
  int i, iordr1, iordr2, k, nr;
  double a, alfa, an1, an2, argr, b, bn, bign, bignmk, nrfac;
  double eps, fac, facn, fden, ff, fg, fk, f1, f2, gamma;
  double ovlfac, rgamm1, rgamm2, rho, rhon, twogp1, zalfa;
  double ta[MAXRP], tb[MAXRP];
  double energy, t;
  
  alfa = FINE_STRUCTURE_CONST;
  nr = n - abs(kappa);
  fk = abs(kappa);
  zalfa = z*alfa;
  gamma = sqrt(fk*fk - zalfa*zalfa);
  twogp1 = gamma*2.0 + 1.0;
  bign = sqrt(n*n - 2.0*(n - fk)*(fk - gamma));
  t = zalfa/(gamma + n - fk);
  eps = 1.0/sqrt(1.0 + t*t);

  energy = -(1.0 - eps)/FINE_STRUCTURE_CONST2;
  
  nrfac = 0.0;
  for (i = 1; i <= nr; i++) nrfac += log(i);
  
  argr = twogp1 + nr;
  rgamm1 = gsl_sf_lngamma(argr);
  argr = twogp1;
  rgamm2 = gsl_sf_lngamma(argr);

  /*
  fac = - sqrt(rgamm1) / (rgamm2*sqrt((double)nrfac)) *
    sqrt(z/(2.0*bign*bign*(bign-kappa)));
  */
  fac = -exp(0.5*rgamm1 - rgamm2 - 0.5*nrfac)*sqrt(z/(2.0*bign*bign*(bign-kappa)));

  if (kappa > 0) fac = -fac;

  fg = fac*sqrt(1.0 + eps);
  ff = fac*sqrt(1.0 - eps);
  
  if (nr == 0) {
    iordr1 = 0;
    iordr2 = 0;
  } else {
    iordr1 = nr - 1;
    iordr2 = nr;
  }

  fac = 1.0;
  facn = 1.0;
  a = -nr;
  an1 = a + 1.0;
  an2 = a;
  b = twogp1;
  bn = b;
  
  k = 0;
  
  while (1) {
    k = k + 1;
    fden = 1.0/(facn*bn);
    if (k <= iordr1) {
      ta[k] = an1 * fden;
    }

    if (k > iordr2) break;
    tb[k] = an2*fden;
    a += 1.0;
    an1 *= a + 1.0;
    an2 *= a;
    b += 1.0;
    bn *= b;
    fac += 1.0;
    facn *= fac;
  }

  p[0] = 0.0;
  q[0] = 0.0;
  fac = 2.0*z/bign;
  bignmk = bign - kappa;
  for (i = 0; i < npts; i++) {
    rho = fac * r[i];
    rhon = rho;
    k = 0;
    f1 = 1.0;
    f2 = 1.0;
    while (1) {
      k = k+1;
      if (k <= iordr1) {
	f1 += ta[k]*rhon;
      }
      if (k > iordr2) break;
      f2 += tb[k]*rhon;
      rhon *= rho;
    }
    f1 *= nr;
    f2 *= bignmk;
    ovlfac = exp(-0.5*rho)*pow(rho, gamma);
    p[i] = fg*ovlfac * (f1-f2);
    q[i] = ff*ovlfac * (f1+f2);
  }

  return energy;
}

int RadialSolver(const cfac_t *cfac, ORBITAL *orb) {
  int ierr;
  int nm, km, k;
  POTENTIAL *pot = cfac->potential;

  if (orb->n > 0) {
    if (orb->n == 1000000) {
      ierr = RadialFreeInner(orb, pot);
    } else {
      if (pot->ib > 0 && orb->n > pot->nb) {
	ierr = RadialBasis(orb, pot);
      } else {
	GetHydrogenicNL(cfac, NULL, NULL, &nm, &km);
	k = GetLFromKappa(orb->kappa);
	k /= 2;
	if (orb->n > nm || k > km) {
	  double z = cfac_get_atomic_number(cfac);
	  if (pot->Navg > 0) z -= pot->Navg - 1;
	  orb->energy = EnergyH(z, (double)(orb->n), orb->kappa);
	  orb->ilast = -1;
	  if (orb->wfun) free(orb->wfun);
	  orb->wfun = NULL;
	  return 0;
	}
	if (orb->n < pot->nmax) {
	  ierr = RadialBound(orb, pot);
	} else {
	  ierr = RadialRydberg(orb, pot);
	}
      }
    }
  } else if (orb->n == 0) {
    ierr = RadialFree(orb, pot);
  } else {
    ierr = RadialBasisOuter(orb, pot);
  }
  return ierr;
}

int LastMaximum(double *p, int i1, int i2) {
  int i;

  /* skip equal last points */
  while (i2 > i1 && p[i2] == p[i2 - 1]) {
    i2--;
  }
  
  if (p[i2] > p[i2-1]) {
    for (i = i2-1; i > i1; i--) {
      if (p[i-1] > p[i]) break;
    }
  } else if (p[i2] < p[i2-1]) {
    for (i = i2-1; i > i1; i--) {
      if (p[i-1] < p[i]) break;
    }
  } else {
    i = i2;
  }
  if (i == i1) {    
    i = i2;
  }
  return i;
} 

void Differential(double *p, double *dp, int i1, int i2) {
  double c = 1.0/24.0;
  double b;
  int i;

  b = -50.0*p[i1];
  b += 96.0*p[i1+1];
  b -= 72.0*p[i1+2];
  b += 32.0*p[i1+3];
  b -= 6.0 *p[i1+4];
  dp[i1] = b*c;

  b = -6.0*p[i1];
  b -= 20.0*p[i1+1];
  b += 36.0*p[i1+2];
  b -= 12.0*p[i1+3];
  b += 2.0 *p[i1+4];
  dp[i1+1] = b*c;

  b = -50.0*p[i2];
  b += 96.0*p[i2-1];
  b -= 72.0*p[i2-2];
  b += 32.0*p[i2-3];
  b -= 6.0 *p[i2-4];
  dp[i2] = -b*c; 

  b = -6.0*p[i2];
  b -= 20.0*p[i2-1];
  b += 36.0*p[i2-2];
  b -= 12.0*p[i2-3];
  b += 2.0 *p[i2-4];
  dp[i2-1] = -b*c;

  for (i = i1 + 2; i < i2-1; i++) {
    b = 2.0*(p[i-2] - p[i+2]);
    b += 16.0*(p[i+1] - p[i-1]);
    dp[i] = b*c;
  }
}    

static double DpDr(int kappa, int i1, double e, const POTENTIAL *pot, double b) {
  double p[MAXRP], bqp, p1, p2 = 0.0;
  int k;
  
  bqp = (b + kappa)*FINE_STRUCTURE_CONST;
  bqp /= (2.0*pot->rad[i1]);
  for (k = i1-1; k <= i1+1; k++) {
    p[k] = 0.5*FINE_STRUCTURE_CONST2*(e - (pot->U[k]+pot->Vc[k]));
    p[k] += 1.0;
    if (k == i1) {
      p2 = 2.0*p[k]*bqp/FINE_STRUCTURE_CONST;
      p2 -= kappa/pot->rad[i1];
      p2 *= pot->dr_drho[i1];
    }
    p[k] = sqrt(p[k]);
  }
  p1 = 0.5*(p[i1+1] - p[i1-1])/p[i1];
  p2 -= p1;
  p1 = 0.5*(pot->dr_drho2[i1+1] - pot->dr_drho2[i1-1]);
  p1 /= pot->dr_drho2[i1];
  bqp = p2 - p1;
  
  return bqp;
}
  
int RadialBasisOuter(ORBITAL *orb, POTENTIAL *pot) {  
  double e, de, delta, emin, emax, dr, ke;
  double *p, p1, qo, qi, bqp, bqp1;
  int n, i, kl, nr, nodes, nnodes, niter;
  int i2, i2m1, i2m2, i2p1, i2p2;
  int ib0, ib1;

  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  n = -orb->n;
  nr = n - kl - 1;
  if (kl < 0 || kl >= n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  
  ib0 = pot->ib-1;
  ib1 = pot->ib1+1;
  niter = 0;
  
  SetPotentialW(pot, 0.0, orb->kappa);
  SetVEffective(kl, pot);
  emin = 1E30;
  for (i = pot->ib; i <= pot->ib1; i++) {
    if (pot->veff[i] < emin) emin = pot->veff[i];
  }
  dr = pot->rad[pot->ib1]-pot->rad[pot->ib];
  ke = TWO_PI*(2*nr+5)/dr;
  emax = 0.5*ke*ke;
  ke = TWO_PI*(nr+1)/dr;
  e = 0.5*ke*ke;
  de = 0.5*e;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    bqp = DpDr(orb->kappa, pot->ib, e, pot, 0.0);
    bqp1 = DpDr(orb->kappa, pot->ib1, e, pot, 0.0);
    i2 = TurningPoints(orb->n, e, pot);
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
    if (nodes < 0) {
      return -1;
    }
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2);
      nodes = CountNodes(p, ib0, i2);
    }
    nnodes = IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
    if (nnodes < 0) {
      return -1;
    }
    nodes += nnodes;
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -2;
  }
  
  niter = 0;
  de = 0.25*ke*ke;
  while (de > ENEABSERR) {
    de *= 0.25;
    if (niter == max_iteration) break;
    while (nodes == nr) {
      niter++;
      if (niter == max_iteration) break;
      e += de;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      bqp = DpDr(orb->kappa, pot->ib, e, pot, 0.0);
      bqp1 = DpDr(orb->kappa, pot->ib1, e, pot, 0.0);
      i2 = TurningPoints(orb->n, e, pot);
      nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
      if (nodes < 0) {
        return -1;
      }
      if (nodes > 1) {
	i2 = LastMaximum(p, pot->ib+20, i2);
	nodes = CountNodes(p, ib0, i2);
      }
      nnodes = IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
      if (nnodes < 0) {
        return -1;
      }
      nodes += nnodes;
    }
    e -= de;
    nodes = nr;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -3;
  }
  niter = 0;
  emax = e;
  if (nr > 0) {
    de = 0.25*ke*ke;
    while (de > ENEABSERR) {
      de *= 0.25;
      if (niter == max_iteration) break;
      while (nodes == nr) {
	niter++;
	if (niter == max_iteration) break;
	e -= de;
	SetPotentialW(pot, e, orb->kappa);
	SetVEffective(kl, pot);
	bqp = DpDr(orb->kappa, pot->ib, e, pot, 0.0);
	bqp1 = DpDr(orb->kappa, pot->ib1, e, pot, 0.0);
	i2 = TurningPoints(orb->n, e, pot);
	nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2, 1.0, 2);
        if (nodes < 0) {
          return -1;
        }
	if (nodes > 1) {
	  i2 = LastMaximum(p, pot->ib+20, i2);
	  nodes = CountNodes(p, ib0, i2);
	}
	nnodes = IntegrateRadial(p, e, pot, i2, 1.0, ib1, bqp1, 1);
        if (nnodes < 0) {
          return -1;
        }
        nodes += nnodes;
      }
      e += de;
      nodes = nr;
    }  
    if (niter == max_iteration) {
      printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	     nodes, nr);
      free(p);
      return -4;
    }
    emin = e;
  }
  e = 0.5*(emin+emax);
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot); 
    bqp = DpDr(orb->kappa, pot->ib, e, pot, 0.0);
    bqp1 = DpDr(orb->kappa, pot->ib1, e, pot, 0.0);
    i2 = TurningPoints(orb->n, e, pot);
    i2p2 = i2 + 2;
    nodes = IntegrateRadial(p, e, pot, ib0, bqp, i2p2, 1.0, 2);
    if (nodes < 0) {
      return -1;
    }
    for (i = ib0; i <= i2p2; i++) {
      p[i] *= pot->dr_drho2[i];
    }
    if (nodes > 1) {
      i2 = LastMaximum(p, pot->ib+20, i2);
    }
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;    
    i2p2 = i2 + 2;
    p1 = p[i2m2]/pot->dr_drho2[i2m2];
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p[i2];
    if (IntegrateRadial(p, e, pot, i2m2, p1, ib1, bqp1, 1) < 0) {
      return -1;
    }
    for (i = i2m2; i <= ib1; i++) {
      p[i] *= pot->dr_drho2[i];
    }
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2];
    delta = qo - qi;
    /*
    printf("%3d %2d %2d %3d %15.8E %15.8E %15.8E %12.5E %12.5E %12.5E\n",
	   niter, nodes, k, i2, e, emin, emax, qo, qi, delta);
    */
    if (delta < -EPS8) {
      emax = e;
      e = 0.5*(emin + e);
    } else if (delta > EPS8) {
      emin = e;
      e = 0.5*(emax + e);
    } else {
      break;
    }
  }
  
  if (niter == max_iteration) {
    printf("Max iteration before finding solution in RadialBasisOuter %d %d\n",
	   nodes, nr);
    free(p);
    return -5;
  }
  
  if (p[pot->ib] < 0) {
    for (i = ib0; i <= ib1; i++) {
      p[i] = -p[i];
    }
  }
  orb->ilast = pot->ib1;
  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;
  if (pot->flag == -1) {
    DiracSmall(orb, pot, ib1);
  }

  return 0;
}

int RadialBasis(ORBITAL *orb, POTENTIAL *pot) {
  double z, z0, e, de, delta, emin, emax;
  double *p, p1, p2 = 0.0, qo, qi, bqp;
  int i, k, kl, nr, nodes, niter;
  int i2, i2m1, i2m2, i2p1, i2p2, i1;
  
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;

  nr = orb->n - kl - 1;

  niter = 0;

  z0 = pot->anum;
  z = z0 - (pot->Navg - 1);

  if (kl > 0) {
    SetPotentialW(pot, 0.0, orb->kappa);
    SetVEffective(kl, pot);
    emin = 1E30;
    for (i = 0; i < pot->maxrp; i++) {
      if (pot->veff[i] < emin) emin = pot->veff[i];
    }
    emin += ENEABSERR;
  } else {
    emin = EnergyH(z, orb->n, orb->kappa);
  }

  e = 1.1*fabs(emin);
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 < 0) {
      nodes = -1;
    } else {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    if (nodes > nr) break;
    e = e*2.0;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasis %d %d\n",
	   nodes, nr);
    free(p);
    return -2;
  }
  emax = e;

  if (kl == 0) {
    double emin0 = 5*EnergyH(z0, orb->n, orb->kappa);
    e = emin;      
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	nodes = -1;
      } else {
	nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
        if (nodes < 0) {
          return -1;
        }
      }
      if (nodes < nr) break;
      e = e*2.0;
      if (nr == 0 && e < emin0) break;
    }
    emin = e;
    if (niter == max_iteration) {
      printf("Max iteration before finding correct nodes in RadialBasis %d %d\n",
	     nodes, nr);
      free(p);
      return -3;
    }
  }
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    i2 = TurningPoints(orb->n, e, pot);
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
    if (nodes < 0) {
      return -1;
    }
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBasis %d %d\n",
	   nodes, nr);
    free(p);
    return -4;
  }

  de = 0.5*(emax - emin);
  while (de > ENEABSERR) {
    de *= 0.25;
    while (nodes == nr) {
      e += de;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    e -= de;
    nodes = nr;
  }
  niter = 0;
  
  bqp = orb->kappa*FINE_STRUCTURE_CONST;
  bqp /= (2.0*pot->rad[pot->ib]);    
  emax = e;
  de = emax - emin;
  while (de > ENEABSERR) {
    de *= 0.25;
    while (nodes == nr) {
      e -= de;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	nodes = -1;
	break;
      }
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    e += de;
    nodes = nr;
  }
  emin = e;
  e = 0.5*(emin+emax);
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot); 
    i2 = TurningPoints(orb->n, e, pot); 
    if (i2 == pot->ib) {
      i2m1 = i2 - 1;
      i2m2 = i2 - 2;
      i2p1 = i2 + 1;
      i2p2 = i2 + 2;      
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      for (k = i2m2-1; k <= i2p2; k++) {
	p[k] *= pot->dr_drho2[k];
	p1 = e - pot->Vc[k] - pot->U[k];
	p1 = sqrt(1.0 + FINE_STRUCTURE_CONST2*p1*0.5);
	p[k] *= p1;
      }
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	    + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
      p2 = qo/(p[i2]*pot->dr_drho[i2]) + orb->kappa/pot->rad[i2];
      p1 = e - pot->Vc[i2] - pot->U[i2];
      p1 = 1.0 + 0.5*FINE_STRUCTURE_CONST2*p1;
      p2 *= 0.5*FINE_STRUCTURE_CONST/p1;
      qo = p2 - bqp;
      if (qo > EPS8) {
	emin = e;
	e = 0.5*(emax + e);
      } else if (qo < -EPS8) {
	emax = e;
	e = 0.5*(emin + e);
      } else {
	break;
      }
    } else {
      i1 = pot->ib;
      for (k = i1-1; k <= i1+1; k++) {
	p[k] = 0.5*FINE_STRUCTURE_CONST2*(e - (pot->U[k]+pot->Vc[k]));
	p[k] += 1.0;
	if (k == i1) {
	  p2 = 2.0*p[k]*bqp/FINE_STRUCTURE_CONST;
	  p2 -= orb->kappa/pot->rad[i1];
	  p2 *= pot->dr_drho[i1];
	}
	p[k] = sqrt(p[k]);
      }
      p1 = 0.5*(p[i1+1] - p[i1-1])/p[i1];
      p2 -= p1;
      p1 = 0.5*(pot->dr_drho2[i1+1] - pot->dr_drho2[i1-1]);
      p1 /= pot->dr_drho2[i1];
      p2 -= p1;
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      for (i = 0; i <= i2p2; i++) {
	p[i] *= pot->dr_drho2[i];
      }
      if (i2 > i1-10) {
	i2 = i1 - 10;
	i2 = LastMaximum(p, 0, i2);
      }
      i2m1 = i2 - 1;
      i2m2 = i2 - 2;
      i2p1 = i2 + 1;
      i2p2 = i2 + 2;
      p1 = p[i2m2]/pot->dr_drho2[i2m2];
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	    + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
      qo /= p[i2];
      if (IntegrateRadial(p, e, pot, i2m2, p1, i1+1, p2, 1) < 0) {
        return -1;
      }
      for (i = i2m2; i <= pot->ib+1; i++) {
	p[i] *= pot->dr_drho2[i];
      }
      qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	    - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
      qi /= p[i2];
      delta = qo - qi;      
      if (delta < -EPS8) {
	emax = e;
	e = 0.5*(emin + e);
      } else if (delta > EPS8) {
	emin = e;
	e = 0.5*(emax + e);
      } else {
	break;
      }
    }
  }

  if (niter == max_iteration) {
    printf("Max iteration before finding solution in RadialBasis %d %d\n",
	   nodes, nr);
    free(p);
    return -7;
  }

  if (i2 == pot->ib) {
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    if (nodes < 0) {
      return -1;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }      
  } else {
    i2p2 = pot->ib + 1;
  }
  if (IsOdd(nodes)) {
    for (i = 0; i <= i2p2; i++) {
      p[i] = -p[i];
    }
  }
  orb->ilast = pot->ib;
  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;

  if (pot->flag == -1) {
    DiracSmall(orb, pot, i2p2);
  }

  return 0;
}

static double InnerProduct(int i1, int i2, const double *p1, const double *p2,
    const POTENTIAL *pot) {
    int k;
    double dwork[MAXRP];

    for (k = i1; k <= i2; k++) {
        dwork[k] = p1[k]*p2[k]*pot->dr_drho[k];
    }
    
    return Simpson(dwork, i1, i2);
}

int RadialBound(ORBITAL *orb, POTENTIAL *pot) {
  double z, z0, e, de, ep, delta, emin, emax;
  double *p, norm2, fact, p1, p2, qo, qi;
  int i, kl, nr, nodes = 0, niter;
  int i2, i2m1, i2m2, i2p1, i2p2;
  
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;

  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;

  nr = orb->n - kl - 1;

  z0 = pot->anum;
  z = z0 - (pot->Navg - 1);

  if (kl > 0) {
    SetPotentialW(pot, 0.0, orb->kappa);
    SetVEffective(kl, pot);
    emin = 1E30;
    for (i = 0; i < pot->maxrp; i++) {
      if (pot->veff[i] < emin) emin = pot->veff[i];
    }
    emin += ENEABSERR;
  } else {
    emin = EnergyH(z, orb->n, orb->kappa);
  }
  
  e = 0.5*emin;
  niter = 0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 > 0) {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      if (nodes > nr) break;
    }
    e *= 0.5;
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBound %d %d %d\n",
	   nodes, nr, i2);
    free(p);
    return -2;
  }
  emax = e;

  niter = 0;
  if (kl == 0) {
    double emin0 = 1.2*EnergyH(z0, orb->n, orb->kappa);
    e = emin;   
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	nodes = -1;
      } else {
	nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
        if (nodes < 0) {
          return -1;
        }
      }
      if (nodes < nr) break;
      e = e*2.0;
      if (nr == 0 && e < emin0) break;
    }
    emin = e;
    if (niter == max_iteration) {
      printf("Max iteration before finding correct nodes in RadialBound %d %d\n",
	     nodes, nr);
      free(p);
      return -3;
    }
  }
  
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    i2 = TurningPoints(orb->n, e, pot);
    if (i2 < 0) {
      nodes = -1;
    } else {
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    if (nodes == nr) break;
    else if (nodes > nr) {
      emax = e;
      e = 0.5*(emin + e);
    } else {
      emin = e;
      e = 0.5*(emax + e);
    }
  }
  if (niter == max_iteration) {
    printf("Max iteration before finding correct nodes in RadialBound %d %d\n",
	   nodes, nr);
    free(p);
    return -4;
  }
  
  de = 0.5*(emax - emin);
  while (de > ENEABSERR) {
    de *= 0.25;
    while (nodes == nr) {
      e += de;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    e -= de;
    nodes = nr;
  }

  niter = 0;
  de = 0.0;
  while (niter < max_iteration) {
    niter++;
    SetPotentialW(pot, e, orb->kappa);
    SetVEffective(kl, pot);
    i2 = TurningPoints(orb->n, e, pot);
    i2p2 = i2 + 2;
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    if (nodes < 0) {
      return -1;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, 0, i2);
    i2m1 = i2 - 1;
    i2m2 = i2 - 2;
    i2p1 = i2 + 1;
    i2p2 = i2 + 2;
    p1 = p[i2m2]/pot->dr_drho2[i2m2];
    p2 = p[i2];    
    qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m1]
	  + 40.0*p[i2] + 60.0*p[i2p1] - 6.0*p[i2p2])/120.0;
    qo /= p2*pot->dr_drho[i2];
    norm2 = InnerProduct(0, i2, p, p, pot);
    if (IntegrateRadial(p, e, pot, i2m2, p1, pot->maxrp-1, 0.0, 0) < 0) {
      return -1;
    }
    for (i = i2m2; i < pot->maxrp; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    qi = (6.0*p[i2m2] - 60.0*p[i2m1] - 40.0*p[i2] + 120.0*p[i2p1]
	  - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
    qi /= p[i2]*pot->dr_drho[i2];
    fact = p2/p[i2];
    norm2 += fact*fact*InnerProduct(i2, pot->maxrp-1, p, p, pot);

    delta = 0.5*p2*p2*(qo - qi)/norm2;
    e = e + delta;
    ep = fabs(e)*ENERELERR;
    if (ep > ENEABSERR) ep = ENEABSERR;
    if (fabs(delta) < ep) break;
  }
  if (niter == max_iteration) {
    printf("Max iteration reached in RadialBound\n");
    free(p);
    return -5;
  }
  ep = sqrt(norm2);
  fact = 1.0/ep;
  if (IsOdd(nodes)) {
    fact = -fact;
  }
 
  ep *= wave_zero;     
  for (i = pot->maxrp-1; i >= 0; i--) {
    if (fabs(p[i]) > ep) break;
  }
  if (IsEven(i)) i++;
  orb->ilast = i;       
  for (i = 0; i < pot->maxrp; i++) {    
    p[i] *= fact;
  }

  nodes = 0;
  for (i = orb->ilast; i >= 0; i--) {
    if (fabs(p[i]) > EPS5) break;
  }
  p1 = p[i];
  for (; i > 0; i--) {
    p2 = fabs(p[i]);
    if (p2 > wave_zero) {
      if ((p1 > 0 && p[i] < 0) ||
	  (p1 < 0 && p[i] > 0)) {
	nodes++;
	p1 = p[i];
      }
    }
  }
  if (nodes != nr) {
    printf("RadialBound: No. nodes changed in iteration %d %d\n", nodes, nr);
    free(p);
    return -6;
  }

  orb->energy = e;
  orb->wfun = p;
  orb->qr_norm = 1.0;

  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1);
  }
 
  return 0;
}

int RadialRydberg(ORBITAL *orb, POTENTIAL *pot) {
#define ME 20
  double z, e, e0;
  int i, kl, niter, ierr;
  double x0, pp, qq, ppi, qqi;
  int i2, i2p, i2m, i2p2, i2m2, nodes, nr;
  double qo, qi, norm2, delta, *p;
  double ep, p1, p2, fact;
  double en[ME], dq[ME], dn, zero=0.0;
  int j, np, nme, one=1;

  z = pot->anum - (pot->Navg - 1);
  kl = orb->kappa;
  kl = (kl < 0)? (-kl-1):kl;
  if (kl < 0 || kl >= orb->n) {
    printf("Invalid orbital angular momentum, L=%d, %d %d\n", 
	   kl, orb->n, orb->kappa);
    return -1;
  }
  e = EnergyH(z, orb->n, orb->kappa);
  p = malloc(sizeof(double)*2*pot->maxrp);
  if (!p) return -1;
  nr = orb->n - kl - 1;

  SetPotentialW(pot, e, orb->kappa);
  SetVEffective(kl, pot);
  i2 = TurningPoints(orb->n, e, pot);
  if (i2 < pot->maxrp - 1.5*pot->asymp) {
    niter = 0;
    dn = orb->n;
    dn = z*z/(dn*dn*dn);
    while (1) {
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      if (nodes != nr) {
	e += 0.5*(nr-nodes)*dn;
      } else {
	break;
      }
    }
    dn *= 0.01;
    while (nodes == nr) {
      e += dn;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
    }
    e -= 1.25*dn;
    while (niter < max_iteration) {
      niter++;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      i2 = TurningPoints(orb->n, e, pot);
      if (i2 < 0) {
	printf("The orbital angular momentum = %d too high\n", kl);
	return -2;
      }
      i2p2 = i2 + 2;
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      for (i = 0; i <= i2p2; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      i2 = LastMaximum(p, 0, i2);
      i2p = i2 + 1;
      i2m = i2 - 1;
      i2p2 = i2 + 2;
      i2m2 = i2 - 2;
      p1 = p[i2m2]/pot->dr_drho2[i2m2];
      p2 = p[i2];
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	    + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
      qo /= p2*pot->dr_drho[i2];
      norm2 = InnerProduct(0, i2, p, p, pot);
      
      if (IntegrateRadial(p, e, pot, i2m2, p1, pot->maxrp-1, 0.0, 0) < 0) {
        return -1;
      }
      for (i = i2m2; i < pot->maxrp; i++) {
	p[i] = p[i] * pot->dr_drho2[i];
      }
      qi = (6.0*p[i2m2] - 60.0*p[i2m] - 40.0*p[i2] + 120.0*p[i2p]
	    - 30.0*p[i2p2] + 4.0*p[i2p2+1])/120.0;
      qi /= p[i2]*pot->dr_drho[i2];
      fact = p2/p[i2];
      norm2 += fact*fact*InnerProduct(i2, pot->maxrp-1, p, p, pot);
      
      delta = 0.5*p2*p2*(qo - qi)/norm2;
      e = e + delta;
      ep = fabs(e)*ENERELERR;
      if (ep > ENEABSERR) ep = ENEABSERR;
      if (fabs(delta) < ep) {
	fact = 1.0/sqrt(norm2);
	if (IsOdd(nodes)) fact = -fact;
	for (i = 0; i < pot->maxrp; i++) {
	  p[i] *= fact;
	}
	break;
      }
    }
    if (niter == max_iteration) {
      printf("Max iteration reached in RadialRydberg\n");
      free(p);
      return -3;
    }    
    pp = 0.0;
  } else {
    if (pot->Navg <= 1) j = 3;
    else j = 1;
    i2p2 = pot->maxrp - j*pot->asymp - 5;

    if (i2p2 >= pot->maxrp) {
      free(p);
      printf("Consider enlarging radial grid\n");
      return -5;
    }
    
    for (i2 = pot->r_core; i2 < i2p2; i2++) {
      if (pot->veff[i2+1] > pot->veff[i2]) break;
    }
    i2 += j*pot->asymp;
    i2p2 = i2 + 2;

    if (i2p2 >= pot->maxrp) {
      free(p);
      printf("Consider enlarging radial grid\n");
      return -5;
    }

    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    if (nodes < 0) {
      return -1;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, pot->r_core, i2);
    i2p = i2 + 1;
    i2m = i2 - 1;
    i2p2 = i2 + 2;
    i2m2 = i2 - 2;

    dn = 1.2/(ME-1.0);
    for (j = 0; j < ME; j++) {
      double n_eff = orb->n - 0.95 + dn*j;
      e = EnergyH(z, n_eff, orb->kappa);
      en[j] = e;
      SetPotentialW(pot, e, orb->kappa);
      SetVEffective(kl, pot);
      nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
      if (nodes < 0) {
        return -1;
      }
      p2 = p[i2];
      qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	    + 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
      qo /= p2;
      x0 = pot->rad[i2];
      DCOUL(z, e, orb->kappa, x0, &pp, &qq, &ppi, &qqi, &ierr);
      if (ierr) {
        printf("DCOUL failed with ierr = %d\n", ierr);
        return -4;
      }
      p1 = (qq/pp)*2.0*x0/FINE_STRUCTURE_CONST - orb->kappa;
      qi = DpDr(orb->kappa, i2, e, pot, p1);
      delta = qo-qi;
      dq[j] = delta;
    }
    for (j = 0; j < ME - 1; j++) {
      if (dq[j] > 0 && dq[j+1] < dq[j]) break;
    }
    i = j;
    for (; j < ME - 1; j++) {
      if (dq[j+1] >= dq[j]) break;
    }
    nme = j - i + 1;
    for (np = i; np <= j; np++) {
      dq[np] = -dq[np];
    }
    UVIP3P(nme, &(dq[i]), &(en[i]), one, &zero, &e);
    i2p2 = pot->maxrp-1;
    nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
    if (nodes < 0) {
      return -1;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] = p[i] * pot->dr_drho2[i];
    }
    i2 = LastMaximum(p, pot->r_core, i2);
    x0 = pot->rad[i2];
    DCOUL(z, e, orb->kappa, x0, &pp, &qq, &ppi, &qqi, &ierr);
    if (ierr) {
      printf("DCOUL failed with ierr = %d\n", ierr);
      return -4;
    }
    norm2 = pp;
    fact = norm2/p[i2];
    if (IsOdd(nodes)) {
      fact = -fact;
    }
    for (i = 0; i <= i2p2; i++) {
      p[i] *= fact;
    }
  }

  i = pot->maxrp-1;
  orb->ilast = i;
  orb->energy = e;
  orb->wfun = p;
  
  e0 = InnerProduct(0, pot->maxrp-1, p, p, pot);
  orb->qr_norm = 1.0/e0;

  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1);
    if (pp != 0) {
      fact = fabs(pp/p[i2]);
      for (i = 0; i < pot->maxrp-1; i++) {
	orb->wfun[i] *= fact;
	orb->wfun[i+pot->maxrp] *= fact;
      }
    }
  }

  return 0;
#undef ME
}
  
int RadialFreeInner(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl;
  int i2;
  double *p, e;

  e = orb->energy;
  kl = orb->kappa;
  if (orb->kappa == 0) {
    printf("Kappa == 0 in RadialFreeInner\n");
    return -1;
  }
  if (pot->ib1 <= 0 && pot->ib <= 0) {
    printf("Boundary not set\n");
    return -2;
  }
  SetPotentialW(pot, e, kl);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  SetVEffective(kl, pot);

  i2 = pot->ib1;
  i2 = Max(i2,pot->ib) + 2;
  if (IntegrateRadial(p, e, pot, 0, 0.0, i2, 1.0, 0) < 0) {
    return -1;
  }
  for (i = i2; i >= 0; i--) {
    p[i] *= pot->dr_drho2[i];
  }
  orb->ilast = i2;
  orb->wfun = p;
  orb->phase = 0.0;
  
  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, i2);
  }

  return 0;
}

/* note that the free states are normalized to have asymptotic 
   amplitude of 1/sqrt(k), */
int RadialFree(ORBITAL *orb, POTENTIAL *pot) {
  int i, kl, nodes;
  int i2, i2p, i2m, i2p2, i2m2;
  double *p, po, qo, e, po1;
  double dfact, da, cs, si, phase0;

  e = orb->energy;
  if (e < 0.0) { 
    printf("Energy < 0 in RadialFree\n");
    return -1;
  }
  kl = orb->kappa;
  if (orb->kappa == 0) {
    printf("Kappa == 0 in RadialFree\n");
    return -1;
  }
  SetPotentialW(pot, e, kl);
  kl = (kl < 0)? (-kl-1):kl;  
  p = malloc(2*pot->maxrp*sizeof(double));
  if (!p) return -1;
  SetVEffective(kl, pot);
  
  i2 = TurningPoints(0, e, pot);
  i2m = i2 - 1;
  i2p = i2 + 1;
  i2m2 = i2 - 2;
  i2p2 = i2 + 2;

  nodes = IntegrateRadial(p, e, pot, 0, 0.0, i2p2, 1.0, 0);
  if (nodes < 0) {
    return -1;
  }

  for (i = i2p2; i >= 0; i--) {
    p[i] *= pot->dr_drho2[i];
  }
  dfact = 1.0 / pot->dr_drho[i2];
  qo = (-4.0*p[i2m2-1] + 30.0*p[i2m2] - 120.0*p[i2m]
	+ 40.0*p[i2] + 60.0*p[i2p] - 6.0*p[i2p2])/120.0;
  qo *= dfact;
  po = p[i2];
  po1 = p[i2p];

  da = Amplitude(p, e, orb->kappa, pot, i2);
  if (!da) {
    return -1;
  }
  cs = p[i2] * qo - da*po;
  si = po / p[i2];
  dfact = (si*si + cs*cs);
  dfact = 1.0/sqrt(dfact);
  phase0 = atan2(si, cs);
  if (phase0 < 0) phase0 += TWO_PI;
    
  if (IsOdd(nodes)) {
    phase0 = (phase0 < M_PI)?(phase0 + M_PI):(phase0-M_PI);
    dfact = -dfact;
  }
  Phase(p, pot, i2, phase0);
  
  p[i2] = po;
  p[i2p] = po1;
  for (i = 0; i < i2p2; i++) {
    p[i] *= dfact;
  }
    
  orb->ilast = i2p;
  orb->wfun = p;
  orb->phase = 0.0;

  orb->qr_norm = 1.0;
  
  if (pot->flag == -1) {
    DiracSmall(orb, pot, -1);
  }

  return 0;
}

/*
** The storage of amplitude and phase in the wfun array is arranged:
** large component: large[i]=amplitude_i, large[i+1]=phase_i,
** small component: small[i]=a_cos_i, small[i+1]=a_sin_i,
** so: P[i] = large[i]*sin(large[i+1])
** and Q[i] = small[i]*cos(large[i+1])+small[i+1]*sin(large[i+1])
*/
int DiracSmall(ORBITAL *orb, POTENTIAL *pot, int i2) {
  int i, i0, i1, kappa;
  double xi, e, *p, a, b;
  double _dwork[MAXRP];
  double _dwork1[MAXRP];
  double _dwork2[MAXRP];

  if (orb->n >= 0) i0 = 0;
  else i0 = pot->ib-1;
  e = orb->energy;
  kappa = orb->kappa;
  p = orb->wfun;
  if (i2 < 0) i2 = orb->ilast;
  i1 = orb->ilast+1;

  for (i = i0; i <= i2; i++) {
    xi = e - pot->Vc[i] - pot->U[i];
    xi = xi*FINE_STRUCTURE_CONST2*0.5;
    _dwork[i] = 1.0 + xi;
    _dwork1[i] = sqrt(_dwork[i])*p[i];
    _dwork2[i] = 1.0/(24.0*pot->dr_drho[i]);
  }
  for (i = i0; i < i1; i++) {
    if (fabs(_dwork1[i]) < 1E-16) {
      p[i+pot->maxrp] = 0.0;
      continue;
    } 
    if (i == i0) {
      b = -50.0*_dwork1[i0];
      b += 96.0*_dwork1[i0+1];
      b -= 72.0*_dwork1[i0+2];
      b += 32.0*_dwork1[i0+3];
      b -= 6.0 *_dwork1[i0+4];
    } else if (i == i0+1) {
      b = -6.0*_dwork1[i0];
      b -= 20.0*_dwork1[i0+1];
      b += 36.0*_dwork1[i0+2];
      b -= 12.0*_dwork1[i0+3];
      b += 2.0 *_dwork1[i0+4];
    } else if (i == i2) {
      b = -50.0*_dwork1[i];
      b += 96.0*_dwork1[i-1];
      b -= 72.0*_dwork1[i-2];
      b += 32.0*_dwork1[i-3];
      b -= 6.0 *_dwork1[i-4];
      b = -b; 
    } else if (i == i2-1) {
      b = -6.0*_dwork1[i+1];
      b -= 20.0*_dwork1[i];
      b += 36.0*_dwork1[i-1];
      b -= 12.0*_dwork1[i-2];
      b += 2.0 *_dwork1[i-3];
      b = -b;
    } else {
      b = 2.0*(_dwork1[i-2] - _dwork1[i+2]);
      b += 16.0*(_dwork1[i+1] - _dwork1[i-1]);
    }
    b *= _dwork2[i];
    b += _dwork1[i]*kappa/pot->rad[i];
    b /= (2.0*_dwork[i]);
    p[i] = _dwork1[i];
    p[i+pot->maxrp] = b*FINE_STRUCTURE_CONST;
  }
  p[pot->maxrp] = p[pot->maxrp+1];

  if (orb->n != 0) {
    for (i = i1; i < pot->maxrp; i++) {
      p[i] = 0;
      p[i+pot->maxrp] = 0.0;
    }
    a = InnerProduct(i0, i1-1, p+pot->maxrp, p+pot->maxrp, pot);
    b = InnerProduct(i0, i1-1, p, p, pot);
    a *= orb->qr_norm;
    b *= orb->qr_norm;
    a = sqrt(a+b);
    orb->qr_norm = a/sqrt(b);
    a = 1.0/a;
    for (i = i0; i < i1; i++) {
      p[i] *= a;
      p[i+pot->maxrp] *= a;
    }
    for (i = 0; i < i0; i++) {
      p[i] = 0.0;
      p[i+pot->maxrp] = 0.0;
    }
    for (i = i1+pot->maxrp; i < 2*pot->maxrp; i++) {
      p[i] = 0.0;
    }
    if (i0 > 0) {
      p[i0-1] = 0.0;
      p[i0+pot->maxrp-1] = 0.0;
    }
    return 0;
  }
  
  for (i = i1; i < pot->maxrp; i += 2) {
    xi = e - pot->Vc[i] - pot->U[i];
    xi = xi*FINE_STRUCTURE_CONST2*0.5;
    _dwork[i] = 1.0 + xi;
    _dwork1[i] = sqrt(_dwork[i])*p[i];
    _dwork2[i] = 0.25/(pot->dr_drho[i]);
  }

  for (i = i1; i < pot->maxrp; i += 2) {
    if (i == i1) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i+2] - _dwork1[i+4];
    } else if (i == pot->maxrp-2) {
      b = -3.0*_dwork1[i] + 4.0*_dwork1[i-2] - _dwork1[i-4];
      b = -b;
    } else {
      b = _dwork1[i+2] - _dwork1[i-2];
    }
    b *= _dwork2[i];
    b += _dwork1[i]*kappa/pot->rad[i];
    a = _dwork1[i]/(p[i]*p[i]);
    p[i] = _dwork1[i];
    p[i+pot->maxrp] = FINE_STRUCTURE_CONST*a/(2.0*_dwork[i]); 
    p[i+1+pot->maxrp] = FINE_STRUCTURE_CONST*b/(2.0*_dwork[i]);
  }

  b = FINE_STRUCTURE_CONST2*orb->energy;
  orb->qr_norm = sqrt((1.0 + b)/(1.0 + 0.5*b));
  return 0;
}

typedef struct {
    double t0;
    double w0;
    double s;
    double e;
} ode_amp_params;


static int ode_func(double t, const double y[], double ydot[], void *params)
{
  ode_amp_params *p = params;

  double w = 2.0*(p->e - (p->w0 + (t - p->t0)*p->s)/t);
  
  ydot[0] = y[1];
  ydot[1] = 1.0/(y[0]*y[0]*y[0]) - y[0]*w;
  
  return GSL_SUCCESS;
}

double Amplitude(double *p, double e, int ka, POTENTIAL *pot, int i0) {
  int i, n, kl1;
  double a, b, xi, r2 = 0.0, r3 = 0.0;
  double z, dk, r0, r_o[MAXRP], v_o[MAXRP], w, v1;
  
  double y[2], rtol = EPS4, atol = EPS12, h0 = -0.01;

  gsl_odeiv2_system ode_sys;
  gsl_odeiv2_driver *ode_drv;
  int gsl_status;
  ode_amp_params params;
  
  ode_sys.function  = ode_func;
  ode_sys.jacobian  = NULL;
  ode_sys.dimension = 2;
  ode_sys.params    = &params;

  ode_drv = gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rk8pd,
                                          h0, atol, rtol);
  
  n = pot->maxrp - 1;
  z = pot->anum - (pot->Navg - 1);
  kl1 = ka*(ka+1);

  dk = sqrt(2.0*e*(1.0+0.5*FINE_STRUCTURE_CONST2*e));
  dk = EPS5*e*dk;

  /* calculate (r & v) grids for the outer region */
  r_o[0] = pot->rad[n];
  v_o[0] = pot->veff[n];
  for (i = 1; i < pot->maxrp; i++) {
    r_o[i] = r_o[i-1]*1.05;
    r2 = r_o[i]*r_o[i];
    r3 = r2*r_o[i];
    a = -z/r_o[i];
    b = e - a;
    xi = sqrt(1.0 + 0.5*FINE_STRUCTURE_CONST2*b);
    xi = xi*xi;
    b = b*b;
    v1 = z/r2;
    w = (-2.0*v1/r_o[i] + 0.75*FINE_STRUCTURE_CONST2*v1*v1/xi - 2*ka*v1/r_o[i]);
    w /= 4.0*xi;
    v_o[i] = a + 0.5*kl1/r2 - 0.5*FINE_STRUCTURE_CONST2*(b - w);
    
    a = z/r2;
    b = kl1/r3;
    if (a < dk && b < dk) {
      break;
    }
  }

  a = FINE_STRUCTURE_CONST2*e;
  b = FINE_STRUCTURE_CONST*z;
  w = 2.0*(e - v_o[i]);

  /* initial values */
  r0   = r_o[i];
  y[0] = pow(w, -0.25);
  y[1] = 0.5*(y[0]/w)*(z*(1+a)/r2-(kl1-b*b)/r3);

  /* integration over the outer region */
  i--;
  for (; i >= 0; i--) {
    params.t0 = r0;
    params.w0 = v_o[i+1]*r0;
    params.s  = (v_o[i]*r_o[i] - params.w0)/(r_o[i] - r0);
    params.e  = e;

    gsl_status = gsl_odeiv2_driver_apply(ode_drv, &r0, r_o[i], y);
    if (gsl_status != GSL_SUCCESS) {
       printf ("ORI ODE error %d\n", gsl_status);
       return 0;
    }
  }
  p[n] = y[0];
  
  /* integration over the inner region */
  for (i = n-1; i >= i0; i--) {
    double r = pot->rad[i];
    params.t0 = r0;
    params.w0 = pot->veff[i+1]*r0;
    params.s  = (pot->veff[i]*r - params.w0)/(r - r0);
    params.e  = e;

    gsl_status = gsl_odeiv2_driver_apply(ode_drv, &r0, r, y);
    if (gsl_status != GSL_SUCCESS) {
       printf ("IRI ODE error %d\n", gsl_status);
       return 0;
    }
    p[i] = y[0];
  }

  gsl_odeiv2_driver_free(ode_drv);

  return y[1];
}    

int Phase(double *p, POTENTIAL *pot, int i1, double phase0) {
  int i;
  double fact;
  double _dwork[MAXRP];
  
  fact = 1.0/3.0;

  for (i = i1; i < pot->maxrp; i++) {
    _dwork[i] = 1.0/p[i];
    _dwork[i] *= _dwork[i];
    _dwork[i] *= pot->dr_drho[i];
  }

  i = i1+1;
  p[i] = phase0;
  for (i = i1+3; i < pot->maxrp; i += 2) {
    p[i] = p[i-2] + (_dwork[i-3] + 4.0*_dwork[i-2] + _dwork[i-1])*fact;
  }
  return 0;
}

static int SetVEffective(int kl, POTENTIAL *pot) {
    int i;
    double kl1 = 0.5*kl*(kl+1);

    for (i = 0; i < pot->maxrp; i++) {
        double r = pot->rad[i];
        double r2 = r*r;
        
        pot->veff[i] = pot->Vc[i] + pot->U[i] + pot->W[i] + kl1/r2;
    }

    return 0;
}

static int TurningPoints(int n, double e, const POTENTIAL *pot) {
  int i, i2;

  if (n == 0) {
    for (i = 10; i < pot->maxrp-5; i++) {
      double r, a, x, de = e - pot->veff[i];
      if (de <= 0) continue;
      r = pot->rad[i];
      a = 40/(pot->ar/sqrt(r) + 2*pot->br/r);
      x = TWO_PI/sqrt(2*de);
      if (x < a) break;
    }
    i2 = i-2;
    if (IsOdd(i2)) (i2)--;
  } else if (n > 0) {
    for (i = pot->maxrp-1; i > 10; i--) {
      if (e - pot->veff[i] > wave_zero) break;
    }
    if (i <= 10) return -2;
    i2 = pot->maxrp-5;
    i2 = Min(i2, i);
    if (pot->ib && n > pot->nb) {
      if (i2 > pot->ib) {
	i2 = pot->ib;
      }
    }
  } else {
    for (i = pot->ib; i < pot->ib1; i++) {
      if (e - pot->veff[i] > wave_zero) break;
    }
    i2 = i+20;
    for (i = pot->ib1-20; i > i2; i--) {
      if (e - pot->veff[i] > wave_zero) break;
    }
    i2 = i;
  }
  
  return i2;
}

static int CountNodes(double *p, int i1, int i2) {
  int n, i;
  double a, p0;

  n = 0;
  i = i2;
  p0 = p[i];
  for (; i >= i1; i--) {
    a = fabs(p[i]);
    if (a > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	n++;
	p0 = p[i];
      }
    }
  }

  return n;
}

/* Coefficient at the RHS of Eq.(16) in gu:2008a */
static double get_dirac_rhs(const POTENTIAL *pot, double e, unsigned int i)
{
    double a = pot->ar, b = pot->br;
    double r = pot->rad[i];
    double x, y, z, t;

    t = a*sqrt(r) + 2*b;
    y = 1/(t*t);

    z = 3.0/4.0*a*a*r + 5*a*b*sqrt(r) + 4*b*b;

    x = 8*(pot->veff[i] - e)*r*r + z*y;

    return x*y;
}

/* Solve the radial DE using Numerov's method from radial index i1 to i2;
   store result in p[i1]..p[i2]. Return node count.
   q = {0,1,2} determines the boundary conditions (value/derivative) for WF
   matching: both, left, or right */
static int IntegrateRadial(double *p, double e, const POTENTIAL *pot,
			   int i1, double p1, int i2, double p2, int q) {
#define NN  4 /* stride of ABAND */
  double x, p0;
  int i, info, m, k, nodes;
  double ABAND[NN*MAXRP];
  gsl_vector_view bv;
  gsl_vector *xv, *dv, *supv, *subv;
  
  m = i2 - i1 - 1;
  if (m < 0) return 0;
  
  if (i2 >= pot->maxrp || i2 < 0) {
    printf("Radial index is outside of array bounds in IntegrateRadial(): %d\n",
        i2);
    return -1;
  }
  
  if (m == 0) {
    if (q == 0) {
      p[i1] = p1;
      p[i2] = p2;
    } else if (q == 1) {
      p[i1] = p1;
      p[i2] = (1 + p2)*p1;
    } else if (q == 2) {
      p[i1] = (1 - p1)*p2;
      p[i2] = p2;
    }
    return 0;
  }

  for (i = i1 + 1; i < i2; i++) {
    p[i] = 0.0;
  }

  for (i = i1 + 1, k = 2; i < i2; i++, k += NN) {
    x = get_dirac_rhs(pot, e, i)/12;

    ABAND[k-1] = 1 - x;
    ABAND[k]   = -2*(1 + 5*x);
    ABAND[k+1] = 1 - x;
  }

  x = get_dirac_rhs(pot, e, i2)/12;
  if (q == 1) {
    k -= NN;
    ABAND[k] += 2*p2*2;
    k -= NN;
    ABAND[k+1] += 2;
  } else {
    p[i2] = p2;
    p[i2-1] += -(1 - x)*p2;
  }
    
  x = get_dirac_rhs(pot, e, i1)/12;
  if (q == 2) {
    k = 2;
    ABAND[k] -= 2*(1 - x)*p1;
    k += NN;
    ABAND[k-1] += (1 - x);
  } else {
    p[i1] = p1;
    p[i1+1] += -(1 - x)*p1;
  }

  bv   = gsl_vector_view_array(p+i1+1, m);
  xv   = gsl_vector_alloc(m);
  dv   = gsl_vector_alloc(m);
  supv = gsl_vector_alloc(m - 1);
  subv = gsl_vector_alloc(m - 1);
  
  for (i = 0; i < m; i++) {
    k = NN*(i + 1) - 2;
    gsl_vector_set(dv, i, ABAND[k]);
    if (i < m - 1) {
      k = NN*(i + 1) - 1;
      gsl_vector_set(subv, i, ABAND[k]);
      k = NN*(i + 2) - 3;
      gsl_vector_set(supv, i, ABAND[k]);
    }
  }
  
  info = gsl_linalg_solve_tridiag(dv, supv, subv, &bv.vector, xv);

  gsl_vector_free(dv);
  gsl_vector_free(supv);
  gsl_vector_free(subv);

  gsl_vector_memcpy(&bv.vector, xv);
  
  gsl_vector_free(xv);

  if (info) {
    printf("Error in integrating the radial equation: %d\n", info);
    return -1;
  }

  if (q == 1) {
    p[i2] = 2*p2*p[i2-1] + p[i2-2];
  } else
  if (q == 2) {
    p[i1] = -2*p1*p[i1+1] + p[i1+2];
  }

  /* now count nodes */
  nodes = 0;
  p0 = p[i2];
  for (i = i2; i >= i1; i--) {
    if (fabs(p[i]) > wave_zero) {
      if ((p0 > 0 && p[i] < 0) ||
	  (p0 < 0 && p[i] > 0)) {
	nodes++;
	p0 = p[i];
      }
    }
  }
    
  return nodes;
}
#undef NN

/* solving "a*sqrt(r) + b*log(r) = rho" for r using r0 as an initial guess */
static double GetRFromRho(double rho, double a, double b, double r0) {
    double e, d1;
    int i;

    e = 1.0;
    i = 0;
    while (fabs(e) > 1E-10) {
        if (i > 100) {
            printf("Newton iteration failed to converge in GetRFromRho\n");
            return -1.0;
        }

        d1 = a*sqrt(r0);
        e = (d1 + b*log(r0) - rho)/(0.5*d1 + b);
        r0 *= (1.0 - e);

        i++;
    }

    return r0;
}

int SetOrbitalRGrid(const cfac_t *cfac, POTENTIAL *pot) {
  double ratio = cfac->rratio, asymp = cfac->rasymp;
  int maxrp = cfac->maxrp;
  int i;  
  double rmin, z, d1, d2;
  double a = 0.0, b, c = 0.0, rmax = 0.0;
  double rho, drho;

  if (pot->maxrp != maxrp) {
    pot->rad      = realloc(pot->rad,      maxrp*sizeof(double));
    pot->dr_drho  = realloc(pot->dr_drho,  maxrp*sizeof(double));
    pot->dr_drho2 = realloc(pot->dr_drho2, maxrp*sizeof(double));
    pot->Vn       = realloc(pot->Vn,       maxrp*sizeof(double));
    pot->Vc       = realloc(pot->Vc,       maxrp*sizeof(double));
    pot->dVc      = realloc(pot->dVc,      maxrp*sizeof(double));
    pot->dVc2     = realloc(pot->dVc2,     maxrp*sizeof(double));
    pot->U        = realloc(pot->U,        maxrp*sizeof(double));
    pot->dU       = realloc(pot->dU,       maxrp*sizeof(double));
    pot->dU2      = realloc(pot->dU2,      maxrp*sizeof(double));
    pot->W        = realloc(pot->W,        maxrp*sizeof(double));
    pot->uehling  = realloc(pot->uehling,  maxrp*sizeof(double));
    pot->veff     = realloc(pot->veff,     maxrp*sizeof(double));
  }
  
  if (!pot->rad      ||
      !pot->dr_drho  ||
      !pot->dr_drho2 ||
      !pot->Vn       ||
      !pot->Vc       ||
      !pot->dVc      ||
      !pot->dVc2     ||
      !pot->U        ||
      !pot->dU       ||
      !pot->dU2      ||
      !pot->W        ||
      !pot->uehling  ||
      !pot->veff) {
    printf("Memory allocation failed in SetOrbitalRGrid()\n");
    return -1;
  }

  pot->maxrp = maxrp;
  
  pot->anum = cfac_get_atomic_number(cfac);
  pot->asymp = asymp;
  
  rmin = cfac->rmin/pot->anum;
  
  z = pot->anum;
  if (pot->Navg > 0) {
    z -= pot->Navg - 1;
  }
  
  if (pot->flag == 0) pot->flag = -1; 

  if (asymp > 0 && ratio > 0) {
    a = asymp*sqrt(2*z)/M_PI;
    c = 1/log(ratio);
    d2 = pot->maxrp - 10 + a*sqrt(rmin) + c*log(rmin);
    rmax = d2*d2/(a*a);
    d1 = 1.0;
    while (d1 > EPS3) {
      double r1 = (d2 - c*log(rmax))/a;
      r1 *= r1;
      d1 = fabs(r1/rmax-1.0);
      rmax = r1;
    }
  } else if (ratio > 0) {
    rmax = -asymp;
    c = 1/log(ratio);
    a = (pot->maxrp - 15 - c*log(rmax/rmin))/(sqrt(rmax) - sqrt(rmin));
  } else if (asymp > 0) {
    rmax = -ratio;
    a = asymp*sqrt(2*z)/M_PI;
    c = (pot->maxrp - 15 + a*(sqrt(rmin) - sqrt(rmax)))/log(rmax/rmin);
  }     
  pot->nmax = sqrt(rmax*z)/2.0;

  d1 = log(rmax/rmin);
  d2 = sqrt(rmax) - sqrt(rmin);
  b = (pot->maxrp - 1 - a*d2)/d1;
  if (b < c) {
    printf("Not enough radial mesh points, ");
    printf("enlarge to at least %d\n", (int) (1 + a*d2 + c*d1));
    return -1;
  }

  drho = (b*d1 + a*d2)/(pot->maxrp - 1);
  pot->rad[0] = rmin;

  rho = a*sqrt(rmin) + b*log(rmin);
  for (i = 1; i < pot->maxrp; i++) {
    rho += drho;
    pot->rad[i] = GetRFromRho(rho, a, b, pot->rad[i-1]);
    if (pot->rad[i] < 0) {
      return -1;
    }
  }

  pot->ar = a;
  pot->br = b;

  for (i = 0; i < pot->maxrp; i++) {
    double r = pot->rad[i];
    pot->dr_drho[i] = 2*r/(a*sqrt(r) + 2*b);
    pot->dr_drho2[i] = sqrt(pot->dr_drho[i]);
  }

  return 0;
}

int SetPotentialZ(cfac_t *cfac) {
    int i;
    POTENTIAL *pot = cfac->potential;

    for (i = 0; i < pot->maxrp; i++) {
        pot->Vn[i] = cfac_get_nucleus_potential(cfac, pot->rad[i]);
    }
    
    return 0;
}

/* Set central potential */
int SetPotentialVc(POTENTIAL *pot) {
    int i;
    double ncore = pot->Navg - 1;
    
    /* nucleus */
    for (i = 0; i < pot->maxrp; i++) {
        pot->Vc[i] = pot->Vn[i];
    }
    
    /* screening core electrons */
    if (ncore > 0 && (pot->a > 0 || pot->lambda > 0)) {
        for (i = 0; i < pot->maxrp; i++) {
            double v;
            double r = pot->rad[i];
        
            v = ncore/r*(1 - exp(-pot->lambda*r)/(1 + r*pot->a));
            pot->Vc[i] += v;
        }
    }

    Differential(pot->Vc, pot->dVc, 0, pot->maxrp-1);
    for (i = 0; i < pot->maxrp; i++) {
        pot->dVc[i] /= pot->dr_drho[i];
    }
    Differential(pot->dVc, pot->dVc2, 0, pot->maxrp-1);
    for (i = 0; i < pot->maxrp; i++) {
        pot->dVc2[i] /= pot->dr_drho[i];
    }

    return 0;
}

int SetPotentialU(POTENTIAL *pot, int n) {
  int i;
  
  if (n < 0) {
    for (i = 0; i < pot->maxrp; i++) { 
      pot->U[i] = 0.0;
      pot->dU[i] = 0.0;
      pot->dU2[i] = 0.0;
    }
    return 0;
  }

  Differential(pot->U, pot->dU, 0, pot->maxrp-1);
  for (i = 0; i < pot->maxrp; i++) {
    pot->dU[i] /= pot->dr_drho[i];
  }
  Differential(pot->dU, pot->dU2, 0, pot->maxrp-1);
  for (i = 0; i < pot->maxrp; i++) {
    pot->dU2[i] /= pot->dr_drho[i];
  }
  return 0;
}

int SetPotentialW(POTENTIAL *pot, double e, int kappa) {
    int i;

    for (i = 0; i < pot->maxrp; i++) {
        double v, dv, dv2, xi, r, t, x, y;
        
        r = pot->rad[i];
        
        v   = pot->Vc[i]   + pot->U[i];
        dv  = pot->dVc[i]  + pot->dU[i];
        dv2 = pot->dVc2[i] + pot->dU2[i];

        xi = e - v;
        t = 1 + 0.5*xi*FINE_STRUCTURE_CONST2;
        x = 0.75*dv*dv*FINE_STRUCTURE_CONST2/t;
        y = -2*kappa*dv/r;
        
        pot->W[i] = -0.5*FINE_STRUCTURE_CONST2*(xi*xi - (x + y + dv2)/(4*t));
    }

    return 0;
}

static double Poly(double x, double n, double c[]) {
  int i;
  double t, r;

  t = 1.0;
  r = 0.0;
  for (i = 0; i < n; i++) {
    r += c[i]*t;
    t *= x;
  }
  
  return r;
}
  
static double UehlingK0(double x) {
  double a[8] = {+0.88357293375,
		 -0.28259817381,
		 -0.58904879578,
		 +0.12500133434,
		 -0.032729913852,
		 +0.0082888574511,
		 -0.0000103277658,
		 +0.0000636436689};
  double b[3] = {-319.999594328,
		 +2.53900995981,
		 1.0};
  double c[3] = {-319.999594333,
		 +2.53901020662,
		 0.0};
  double d[5] = {5.018065179,
		 71.51891262,
		 211.6209929,
		 31.40327478,
		 -1.0};
  double e[5] = {2.669207401,
		 51.72549669,
		 296.9809720,
		 536.4324164,
		 153.5335924};

  double r, x2;

  if (x <= 0.0) {
    return a[0];
  }
  
  if (x <= 1.0) {
    r = Poly(x, 8, a);
    if (x) {
      x2 = x*x;
      r += x*log(x)*Poly(x2, 3, b)/Poly(x2, 3, c);
    }
  } else {
    r = exp(-x)/pow(x, 1.5);
    x2 = 1.0/x;
    r *= Poly(x2, 5, d)/Poly(x2, 5, e);
  }

  return r;
}

static double UehlingK1(double x) {
  double a[8] = {-0.71740181754,
		 1.1780972274,
		 -0.37499963087,
		 0.13089675530,
		 -0.038258286439,
		 -0.0000242972873,
		 -0.0003592014867,
		 -0.0000171700907};
  double b[3] = {-64.0514843293,
		 0.711722714285,
		 1.0};
  double c[3] = {64.0514843287,
		 -0.711722686403,
		 0.0008042207748};
  double d[5] = {217.2386409,
		 1643.364528,
		 2122.244512,
		 1.0};
  double e[5] = {115.5589983,
		 1292.191441,
		 3831.198012,
		 2904.410075,
		 0.0};

  double r, x2;

  if (x <= 0.0) {
    return 0.0;
  }
  
  if (x <= 1.0) {
    r = Poly(x, 8, a);
    x2 = x*x;
    r += log(x)*Poly(x2, 3, b)/Poly(x2, 3, c);
  } else {
    r = exp(-x)/pow(x, 1.5);
    x2 = 1.0/x;
    r *= Poly(x2, 5, d)/Poly(x2, 5, e);
  }

  return r;
}

double UehlingL0(double x) {
  double f[6] = {1.990159,
		 -2.397605,
		 1.046471,
		 -0.367066,
		 0.063740,
		 -0.037058};
  double g[3] = {0.751198,
		 0.128889,
		 0.020886};
  double h[2] = {-0.444444,
		 -0.003472};
  
  double logx, x2, r;
  
  if (x <= 0.0) {
    return f[0];
  }
  
  if (x <= 2.0) {
    logx = log(x);
    r = Poly(x, 6, f);
    if (x) {
      x2 = x*x;
      r += x*Poly(x2, 3, g)*logx;
      x2 = x2*x2;
      r += x*Poly(x2, 2, h)*logx*logx;
    }
  } else {
    r = 0.0;
  }
  
  return r;
}

double UehlingL1(double x) {
  double f[6] = {1.646407,
		 -2.092942,
		 0.962310,
		 -0.254960,
		 0.164404,
		 0.0};
  double g[3] = {0.137691,
		 -0.416667,
		 -0.097486};
  double h[2] = {0.444444,
		 0.017361};
  
  double logx, x2, r;
  
  if (x <= 0.0) {
    return 0.0;
  }
  
  if (x <= 2.0) {
    logx = log(x);
    r = Poly(x, 6, f);
    x2 = x*x;
    r += Poly(x2, 3, g)*logx;
    x2 = x2*x2;
    r += Poly(x2, 2, h)*logx*logx;
  } else {
    r = 0.0;
  }
  
  return r;
}

int SetPotentialUehling(cfac_t *cfac, int vp) {
  int i, j, i1;
  double a, b, r0, r, rm, rp, rn;
  double v, d, c, z;
  POTENTIAL *pot = cfac->potential;
  double _dwork[MAXRP];

  if (vp <= 0) return 0;
  
  z = cfac_get_atomic_number(cfac);

  r0 = 3.86159E-11/RBOHR;
  a = -2.0*z*FINE_STRUCTURE_CONST/(3.0*M_PI);
  b = -z*FINE_STRUCTURE_CONST2/(M_PI*M_PI);
  
  for (i = 0; i < pot->maxrp-1; i++) {
    r = pot->rad[i]*2.0/r0;
    pot->uehling[i] = a*UehlingK1(r);
    if (vp > 1) {
      pot->uehling[i] += b*UehlingL1(r);
    }
    pot->uehling[i] /= pot->rad[i];
  }
  
  rn = cfac_get_atomic_rn(cfac);

  if (rn) {
    a = -2.0*r0*FINE_STRUCTURE_CONST/3.0;
    b = -r0*FINE_STRUCTURE_CONST2/M_PI;
    v = 4.0*M_PI*rn*rn*rn/3.0;
    d = z/v;
    for (i = 0; i < pot->maxrp-1; i++) {
      if (pot->rad[i] > rn) break;
    }
    i1 = i;
    for (i = 0; i < pot->maxrp-1; i++) {
      r = pot->rad[i];
      for (j = 0; j < i1; j++) {
	rp = r + pot->rad[j];
	rp = 2.0*rp/r0;
	rm = fabs(r - pot->rad[j]);
	rm = 2.0*rm/r0;
	_dwork[j] = a*(UehlingK0(rm) - UehlingK0(rp));
	if (vp > 1) {
	  _dwork[j] += b*(UehlingL0(rm) - UehlingL0(rp));
	}
	_dwork[j] *= pot->rad[j]*pot->dr_drho[j]*d;
      }
      j = i1 - 1;
      v = Simpson(_dwork, 0, j) / pot->rad[i];
      rp = r + rn;
      rp = 2.0*rp/r0;
      rm = fabs(r - rn);
      rm = 2.0*rm/r0;
      _dwork[i1] = a*(UehlingK0(rm) - UehlingK0(rp));
      if (vp > 1) {
	_dwork[i1] += b*(UehlingL0(rm) - UehlingL0(rp));
      }
      _dwork[i1] *= rn*d;
      _dwork[j] /= pot->dr_drho[j];
      c = 0.5*(_dwork[j] + _dwork[i1])*(rn - pot->rad[j]);
      v += c/pot->rad[i];
      if (fabs(1.0 - v/pot->uehling[i]) < EPS3) break;
      pot->uehling[i] = v;
    }
  }
  return 0;
}
