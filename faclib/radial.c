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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_multimin.h>

#include "cfacP.h"
#include "coulomb.h"
#include "recouple.h"
#include "angular.h"
#include "grid.h"
#include "interpolation.h"

#define Large(orb) ((orb)->wfun)
#define Small(orb) ((orb)->wfun + potential->maxrp)

void InitOrbitalData(void *p, int n) {
  ORBITAL *d;
  int i;
  
  d = p;
  for (i = 0; i < n; i++) {
    d[i].wfun = NULL;
    d[i].phase = 0.0;
    d[i].ilast = -1;
  }
}

void InitYkData(void *p, int n) {
  SLATER_YK *d;
  int i;

  d = p;
  for (i = 0; i < n; i++) {
    d[i].npts = -1;
    d[i].yk = NULL;
  }
}

static int FreeSimpleArray(MULTI *ma) {
  MultiFreeData(ma);
  return 0;
}

void FreeMultipole(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

void FreeYkData(void *p) {
  SLATER_YK *dp;
  
  dp = p;
  if (dp->npts > 0) free(dp->yk);
}

int FreeMultipoleArray(cfac_t *cfac) {
  MultiFreeData(cfac->multipole_array);
  return 0;
}

int FreeMomentsArray(cfac_t *cfac) {
  MultiFreeData(cfac->moments_array);
  return 0;
}

int FreeGOSArray(cfac_t *cfac) {
  MultiFreeData(cfac->gos_array);
  return 0;
}

int FreeYkArray(cfac_t *cfac) {
  MultiFreeData(cfac->yk_array);
  return 0;
}
  
void SetSlaterCut(cfac_t *cfac, int k0, int k1) {
  if (k0 > 0) {
    cfac->slater_cut.kl0 = 2*k0;
  } else {
    cfac->slater_cut.kl0 = 1000000;
  } 
  if (k1 > 0) {
    cfac->slater_cut.kl1 = 2*k1;
  } else {
    cfac->slater_cut.kl1 = 1000000;
  }
}

void SetSE(cfac_t *cfac, int n) {
  cfac->qed.se = n;
}

void SetVP(cfac_t *cfac, int n) {
  cfac->qed.vp = n;
}

void SetBreit(cfac_t *cfac, int n) {
  cfac->qed.br = n;
}

void SetMS(cfac_t *cfac, int nms, int sms) {
  cfac->qed.nms = nms;
  cfac->qed.sms = sms;
}

int SetAWGrid(cfac_t *cfac, int n, double awmin, double awmax) {
  if (awmin < 1E-3) {
    awmin = 1E-3;
    awmax = awmax + 1E-3;
  }
  cfac->n_awgrid = SetTEGrid(cfac->awgrid, NULL, n, awmin, awmax);

  return 0;
}

void SetOptimizeMaxIter(cfac_t *cfac, int m) {
  cfac->optimize_control.maxiter = m;
}

void SetOptimizeStabilizer(cfac_t *cfac, double m) {
  if (m < 0) {
    cfac->optimize_control.iset = 0;
  } else {
    cfac->optimize_control.iset = 1;
    cfac->optimize_control.stabilizer = m;
  }
}

void SetOptimizeTolerance(cfac_t *cfac, double c) {
  cfac->optimize_control.tolerance = c;
}

void SetOptimizePrint(cfac_t *cfac, int m) {
  cfac->optimize_control.iprint = m;
}

void SetOptimizeControl(cfac_t *cfac, double tolerance, double stabilizer, 
			int maxiter, int iprint) {
  cfac->optimize_control.maxiter = maxiter;
  cfac->optimize_control.stabilizer = stabilizer;
  cfac->optimize_control.tolerance = tolerance;
  cfac->optimize_control.iprint = iprint;  
  cfac->optimize_control.iset = 1;
}

void SetScreening(cfac_t *cfac, int n_screen, int *screened_n, 
		  double screened_charge, int kl) {
  cfac->optimize_control.screened_n = screened_n;
  cfac->optimize_control.screened_charge = screened_charge;
  cfac->optimize_control.n_screen = n_screen;
  cfac->optimize_control.screened_kl = kl;
}

int SetRadialGrid(cfac_t *cfac,
    int maxrp, double ratio, double asymp, double rmin) {
  POTENTIAL *potential = cfac->potential;
  if (maxrp > MAXRP) {
    printf("MAXRP must be <= %d\n", MAXRP);
    printf("to enlarge the limit, change MAXRP in consts.h\n");
    return -1;
  }
  if (maxrp < 0) maxrp = DMAXRP;
  potential->maxrp = maxrp;
  if (asymp < 0 && ratio < 0) {
    asymp = GRIDASYMP;
    ratio = GRIDRATIO;
  }
  if (rmin <= 0) rmin = GRIDRMIN;
  potential->rmin = rmin;
  if (ratio == 0) potential->ratio = GRIDRATIO;
  else potential->ratio = ratio;
  if (asymp == 0) potential->asymp = GRIDASYMP;
  else potential->asymp = asymp;
  potential->flag = 0;
  return 0;
}

static void AdjustScreeningParams(POTENTIAL *potential, const double *u) {
  int i;
  double c;
  
  c = 0.5*u[potential->maxrp-1];
  for (i = 0; i < potential->maxrp; i++) {
    if (u[i] > c) break;
  }
  potential->lambda = log(2.0)/potential->rad[i];
}

/* w - electron density distribution of the average config */
static int PotentialHX(const cfac_t *cfac, double *u, double *w) {
  int i, j, k1, jmax, m, jm;
  ORBITAL *orb1;
  double large, small, a, b, c, d0, d1, d;
  POTENTIAL *potential = cfac->potential;
  const AVERAGE_CONFIG *acfg = &cfac->acfg;
  double yk[MAXRP];
  
  if (potential->Navg < 1+EPS3) return -1;

  for (m = 0; m < potential->maxrp; m++) {
    u[m] = 0.0;
    w[m] = 0.0;
  }
  
  jmax = -1;
  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(cfac, acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(cfac, k1);
    if (orb1->wfun == NULL) continue;
    for (m = 0; m <= orb1->ilast; m++) {
      large = Large(orb1)[m];
      small = Small(orb1)[m];
      w[m] += acfg->nq[i]*(large*large + small*small);
    }
    if (jmax < orb1->ilast) jmax = orb1->ilast;    
  }
  if (jmax < 0) return jmax;

  for (i = 0; i < acfg->n_shells; i++) {
    k1 = OrbitalExists(cfac, acfg->n[i], acfg->kappa[i], 0.0);
    if (k1 < 0) continue;
    orb1 = GetOrbital(cfac, k1);
    if (orb1->wfun == NULL) continue;
    GetYk(cfac, 0, yk, orb1, orb1, k1, k1, 1);
    for (m = 0; m <= jmax; m++) {    
      u[m] += acfg->nq[i]*yk[m];
      if (w[m]) {
	large = Large(orb1)[m];
	small = Small(orb1)[m];  
	b = large*large + small*small;
	a = yk[m]*acfg->nq[i]*b/w[m];
	u[m] -= a;
      }
    }
  }

  c = potential->Navg - 1;
  for (jm = jmax; jm >= 10; jm--) {
    if (fabs(w[jm]) > EPS6 && c > u[jm] && u[jm] > u[jm-1]) {
      break;
    }
  }
  d0 = log(potential->rad[jm-1]);
  d1 = log(potential->rad[jm]);
  a = log(c - u[jm-1]);
  b = log(c - u[jm]);
  d = (b-a)/(d1-d0);    
  for (j = jm+1; j < potential->maxrp; j++) {
    u[j] = d*(log(potential->rad[j]/potential->rad[jm])) + b;
    u[j] = c - exp(u[j]);
  }

  for (m = jmax; m > 50; m--) {
    if (fabs(u[m]-c) > EPS6) break;
  }
  potential->r_core = m+1;

  return jmax;
}

static double SetPotential(cfac_t *cfac, int iter, double *v) {
  int jmax, i, j, k;
  double *u, a, b, c, r;
  POTENTIAL *potential = cfac->potential;
  AVERAGE_CONFIG *acfg = &cfac->acfg;
  double w[MAXRP];

  u = potential->U;

  jmax = PotentialHX(cfac, u, w);

  if (jmax > 0) {
    if (iter < 3) {
      r = 1.0;
      for (j = 0; j < potential->maxrp; j++) {
	v[j] = u[j];
      }
    } else {
      r = 0.0;
      k = 0;
      a = cfac->optimize_control.stabilizer;
      b = 1.0 - a;
      for (j = 0; j < potential->maxrp; j++) {
	if (u[j] + 1.0 != 1.0) {
	  r += fabs(1.0 - v[j]/u[j]);
	  k++;
	}
	u[j] = b*v[j] + a*u[j];
	v[j] = u[j];
      }
      r /= k;
    }
    
    AdjustScreeningParams(potential, u);
    SetPotentialVc(potential);
    
    for (j = 0; j < potential->maxrp; j++) {
      u[j] = u[j]/potential->rad[j] + potential->Vn[j] - potential->Vc[j];
    }
    SetPotentialU(potential, 0);
  } else {
    if (potential->Navg < 1.0+EPS3) {
      SetPotentialVc(potential);
      SetPotentialU(potential, -1);
      return 0.0;
    }
    r = cfac_get_atomic_number(cfac);
    b = (1.0 - 1.0/potential->Navg);
    for (i = 0; i < acfg->n_shells; i++) {
      a = acfg->nq[i];
      c = acfg->n[i];
      c = r/(c*c);
      for (j = 0; j < potential->maxrp; j++) {
	u[j] += a*b*(1.0 - exp(-c*potential->rad[j]));
      }
    }
    
    AdjustScreeningParams(potential, u);
    SetPotentialVc(potential);
    
    for (j = 0; j < potential->maxrp; j++) {
      u[j] = u[j]/potential->rad[j] + potential->Vn[j] - potential->Vc[j];
    }
    SetPotentialU(potential, 0);
    
    r = 1.0;
  }
  
  return r;
}

int GetPotential(const cfac_t *cfac, char *fn) {
  const AVERAGE_CONFIG *acfg = &(cfac->acfg);
  FILE *f;
  int i;
  double u[MAXRP], w[MAXRP];
  POTENTIAL *potential = cfac->potential;

  f = fopen(fn, "w");
  if (!f) return -1;
  
  fprintf(f, "#   Navg = %10.3E\n", potential->Navg);
  fprintf(f, "# Lambda = %10.3E\n", potential->lambda);
  fprintf(f, "#      A = %10.3E\n", potential->a);
  fprintf(f, "#     ar = %10.3E\n", potential->ar);
  fprintf(f, "#     br = %10.3E\n", potential->br);
  fprintf(f, "#     rb = %10.3E\n", potential->rad[potential->ib]);
  fprintf(f, "#    rb1 = %10.3E\n", potential->rad[potential->ib1]);
  fprintf(f, "# r_core = %10.3E\n", potential->rad[potential->r_core]);
  fprintf(f, "#     nb = %d\n",     potential->nb);

  PotentialHX(cfac, u, w);

  fprintf(f, "# Mean configuration:\n");
  for (i = 0; i < acfg->n_shells; i++) {
    fprintf(f, "# %2d %2d\t%10.3E\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
  }
  fprintf(f, "\n\n");
  for (i = 0; i < potential->maxrp; i++) {
    fprintf(f, "%5d %14.8E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
	    i, potential->rad[i], potential->Vn[i], 
	    potential->Vc[i], potential->U[i], u[i], w[i]);
  }

  fclose(f);  
  
  return 0;
}

double GetResidualZ(const cfac_t *cfac) {
  double z;
  POTENTIAL *potential = cfac->potential;
  
  z = cfac_get_atomic_number(cfac);
  if (potential->Navg > 0) z -= potential->Navg - 1;
  return z;
}

double GetRMax(cfac_t *cfac) {
  return cfac->potential->rad[cfac->potential->maxrp-10];
}

int SetAverageConfig(cfac_t *cfac, int nshells, int *n, int *kappa, double *nq) {
  int i;
  if (nshells <= 0) return -1;
  cfac->acfg.kappa = realloc(cfac->acfg.kappa, 
					 sizeof(int)*nshells);
  cfac->acfg.nq = realloc(cfac->acfg.nq, 
					 sizeof(double)*nshells);
  cfac->acfg.n = realloc(cfac->acfg.n, 
				     sizeof(int)*nshells);
  for (i = 0; i < nshells; i++) {
    cfac->acfg.n[i] = n[i];
    cfac->acfg.kappa[i] = kappa[i];
    cfac->acfg.nq[i] = nq[i];
  }
  cfac->acfg.n_shells = nshells;
  cfac->acfg.n_cfgs = 1;
  return 0;
}

static int OptimizeLoop(cfac_t *cfac) {
  double tol, a, b;
  ORBITAL orb_old, *orb;
  int i, k, iter, no_old;
  AVERAGE_CONFIG *acfg = &cfac->acfg;
  double vbuf[MAXRP];
  
  no_old = 0;
  iter = 0;
  tol = 1.0; 
  while (tol > cfac->optimize_control.tolerance) {
    if (iter > cfac->optimize_control.maxiter) break;
    a = SetPotential(cfac, iter, vbuf);
    FreeYkArray(cfac);
    tol = 0.0;
    for (i = 0; i < acfg->n_shells; i++) {
      k = OrbitalExists(cfac, acfg->n[i], acfg->kappa[i], 0.0);
      if (k < 0) {
	orb_old.energy = 0.0;
	orb = GetNewOrbital(cfac);
	orb->kappa = acfg->kappa[i];
	orb->n = acfg->n[i];
	orb->energy = 1.0;
	no_old = 1;	
      } else {
	orb = GetOrbital(cfac, k);
	if (orb->wfun == NULL) {
	  orb_old.energy = 0.0;
	  orb->energy = 1.0;
	  orb->kappa = acfg->kappa[i];
	  orb->n = acfg->n[i];
	  no_old = 1;	
	} else {
	  orb_old.energy = orb->energy; 
	  if (orb->wfun) free(orb->wfun);
	  no_old = 0;
	}
      }

      if (SolveDirac(cfac, orb) < 0) {
	return -1;
      }
      
      if (no_old) { 
	tol = 1.0;
	continue;
      } 
      b = fabs(1.0 - orb_old.energy/orb->energy);
      if (tol < b) tol = b;
    }
    if (cfac->optimize_control.iprint) {
      printf("%4d %13.5E %13.5E\n", iter, tol, a);
    }
    if (tol < a) tol = a;
    iter++;
  }

  return iter;
}

int OptimizeRadial(cfac_t *cfac, int ng, int *kg, double *weight) {
  AVERAGE_CONFIG *acfg = &(cfac->acfg);
  double z;
  int iter, i;
  POTENTIAL *potential = cfac->potential;
  
  if (ng > 0) {
    if (ng > 1) {
      printf("\nWarning: more than 1 configuration groups");
      printf(" are used in OptimizeRadial.\n");
      printf("It is usually best to use the lowest lying configuration group.\n");
      printf("Make sure that you know what you are doing.\n\n");
    }
    if (acfg->n_shells > 0) {      
      acfg->n_cfgs = 0;
      acfg->n_shells = 0;
      free(acfg->n);
      free(acfg->kappa);
      free(acfg->nq);
      acfg->n = NULL;
      acfg->nq = NULL;
      acfg->kappa = NULL;
    }
    MakeAverageConfig(cfac, ng, kg, weight); 
  } else {
    if (acfg->n_shells <= 0) {
      printf("No average configuration exist.\n");
      printf("Specify with AvgConfig, ");
      printf("or give config groups to OptimizeRadial.\n");
      return -1;
    }
  }
  
  potential->Navg = 0.0;
  for (i = 0; i < acfg->n_shells; i++) {
    if (cfac->optimize_control.iprint) {
      printf("%d %d %f\n", acfg->n[i], acfg->kappa[i], acfg->nq[i]);
    }
    potential->Navg += acfg->nq[i];
  }

  /* setup the radial grid if not yet */
  if (potential->flag == 0) {
    SetOrbitalRGrid(potential);
  }

  SetPotentialZ(cfac);
  z = cfac_get_atomic_number(cfac);
  if (potential->Navg > 0.0) {
    z -= potential->Navg - 1;
  }
  potential->a = 0.0;
  potential->lambda = 0.5*z;
  if (potential->Navg > 1) {
    potential->r_core = potential->maxrp - 5;
  } else {
    potential->r_core = 50;
  }

  if (cfac->optimize_control.iset == 0) {
    cfac->optimize_control.stabilizer = 0.25 + 0.75*(z/cfac_get_atomic_number(cfac));
  }

  iter = OptimizeLoop(cfac);

  if (potential->uehling[0] == 0.0) {
    SetPotentialUehling(cfac, cfac->qed.vp);
  }

  return iter;
}      


static double EnergyFunc(const gsl_vector *v, void *params) {
  double lambda, a, avg;
  cfac_t *cfac = (cfac_t *) params;
  POTENTIAL *potential = cfac->potential;
  
  lambda = gsl_vector_get(v, 0);
  a      = gsl_vector_get(v, 1);

  if (a < -EPS10) return 0.0;
  if (lambda <= 0.0) return 0.0;

  potential->lambda = lambda;
  potential->a = a;
  SetPotentialVc(potential);
  ReinitRadial(cfac, 1);
  ClearOrbitalTable(cfac, 0);
  avg = AverageEnergyAvgConfig(cfac);

  /* printf("x[0]=%g, x[1]=%g\n", lambda, a); */

  return avg;
}


double subplx(cfac_t *cfac, double (* my_f) (const gsl_vector *v, void *params),
    int n, double xtol, int maxfun, const gsl_vector *ss, gsl_vector *x,
    int *iter)
{
    const gsl_multimin_fminimizer_type *T = 
      gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function minex_func;

    int status;
    double size;

    *iter = 0;
    
    /* Initialize method and iterate */
    minex_func.n = n;
    minex_func.f = my_f;
    minex_func.params = cfac;

    s = gsl_multimin_fminimizer_alloc(T, n);
    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

    do {
        (*iter)++;
        status = gsl_multimin_fminimizer_iterate(s);

        if (status) 
          break;

        size = gsl_multimin_fminimizer_size(s);
        
        status = gsl_multimin_test_size(size, xtol);
    } while (status == GSL_CONTINUE && *iter < maxfun);
    
    gsl_vector_memcpy(x, s->x);

    gsl_multimin_fminimizer_free(s);

    return status;
}


int RefineRadial(cfac_t *cfac, int maxfun, int msglvl) {
  int n, status, nfe = 0;
  double xtol;
  double f0, f;
  gsl_vector *x, *ss;
  POTENTIAL *potential = cfac->potential;
  
  if (maxfun <= 0) maxfun = 250;
  xtol = EPS3;
  n = 2;

  x = gsl_vector_alloc(n);
  gsl_vector_set(x, 0, potential->lambda);
  gsl_vector_set(x, 1, potential->a);

  ss = gsl_vector_alloc(n);
  gsl_vector_set_all(ss, 0.01);

  f0 = EnergyFunc(x, cfac);
  if (msglvl > 0) {
    printf("%10.3E %10.3E %15.8E\n", potential->lambda, potential->a, f0);
  }

  status = subplx(cfac, EnergyFunc, n, xtol, maxfun, ss, x, &nfe);

  f = EnergyFunc(x, cfac);
  if (msglvl > 0) {
    printf("%10.3E %10.3E %15.8E %d\n", potential->lambda, potential->a, f, nfe);
  }
  if (status != GSL_SUCCESS) {
    if (f > f0) {
      printf("Error in RefineRadial: %d\n", status);
      return status;
    } else if (msglvl > 0) {
      printf("Warning in RefineRadial: %d\n", status);
    }
  }
  
  return 0;
}

int SolveDirac(const cfac_t *cfac, ORBITAL *orb) {
  int err;
  POTENTIAL *potential = cfac->potential;

  err = 0;  
  potential->flag = -1;
  err = RadialSolver(cfac, orb);
  if (err) { 
    printf("Error ocuured in RadialSolver, %d\n", err);
    printf("%d %d %10.3E\n", orb->n, orb->kappa, orb->energy);
    exit(1);
  }

  return err;
}

int WaveFuncTable(cfac_t *cfac, char *s, int n, int kappa, double e) {
  int i, k;
  FILE *f;
  ORBITAL *orb;
  double z, a, ke, y;
  POTENTIAL *potential = cfac->potential;

  e /= HARTREE_EV;
  k = OrbitalIndex(cfac, n, kappa, e);
  if (k < 0) return -1;
  f = fopen(s, "w");
  if (!f) return -1;
  
  orb = GetOrbitalSolved(cfac, k);
  
  fprintf(f, "#      n = %2d\n", n);
  fprintf(f, "#  kappa = %2d\n", kappa);
  fprintf(f, "# energy = %15.8E\n", orb->energy*HARTREE_EV);
  fprintf(f, "#     vc = %15.8E\n", MeanPotential(cfac, k, k)*HARTREE_EV);
  if (n != 0) {
    fprintf(f, "\n\n");
    if (n < 0) k = potential->ib;
    else k = 0;
    for (i = k; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n", 
	      i, potential->rad[i], 
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
  } else {
    a = GetPhaseShift(cfac, k);
    while(a < 0) a += TWO_PI;
    a -= (int)(a/TWO_PI);
    fprintf(f, "#  phase = %15.8E\n", a);
    fprintf(f, "\n\n");
    z = GetResidualZ(cfac);
    e = orb->energy;
    a = FINE_STRUCTURE_CONST2 * e;
    ke = sqrt(2.0*e*(1.0+0.5*a));
    y = (1.0+a)*z/ke; 
    for (i = 0; i <= orb->ilast; i++) {
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n", 
	      i, potential->rad[i],
	      (potential->Vc[i])*potential->rad[i],
	      potential->U[i] * potential->rad[i],
	      Large(orb)[i], Small(orb)[i]); 
    }
    for (; i < potential->maxrp; i += 2) {
      a = ke * potential->rad[i];
      a = a + y*log(2.0*a);
      a = Large(orb)[i+1] - a;
      a = a - ((int)(a/(TWO_PI)))*TWO_PI;
      if (a < 0) a += TWO_PI;
      fprintf(f, "%-4d %14.8E %13.6E %13.6E %13.6E %13.6E\n",
	      i, potential->rad[i],
	      Large(orb)[i], Large(orb)[i+1], 
	      Small(orb)[i], a);
    }
  }

  fclose(f);

  return 0;
}

static double PhaseRDependent(double x, double eta, double b) {
  double tau, tau2, y, y2, t, a1, a2, sb;
  
  y = 1.0/x;
  y2 = y*y;
  tau2 = 1.0 + 2.0*eta*y - b*y2;
  tau = x*sqrt(tau2);
  
  t = eta*log(x+tau+eta) + tau - eta;
  if (b > 0.0) {
    sb = sqrt(b);
    a1 = b - eta*x;
    a2 = tau*sb;
    tau2 = atan2(a1, a2);
    tau2 -= atan2(-eta, sb);
    t += sb*tau2;
  } else if (b < 0.0) {
    b = -b;
    sb = sqrt(b);
    a1 = 2.0*(b+eta*x)/(sb*b*x);
    a2 = 2.0*tau/(b*x);
    tau2 = log(a1+a2);
    tau2 -= log(2.0*(eta/sb + 1.0)/b);
    t -= sb*tau2;
  }
  
  return t;
}

double GetPhaseShift(cfac_t *cfac, int k) {
  ORBITAL *orb;
  double phase1, r, y, z, ke, e, a, b1;
  int i;
  POTENTIAL *potential = cfac->potential;

  orb = GetOrbitalSolved(cfac, k);
  if (orb->n > 0) return 0.0;

  if (orb->phase) return orb->phase;

  z = GetResidualZ(cfac);
  e = orb->energy;
  a = FINE_STRUCTURE_CONST2 * e;
  ke = sqrt(2.0*e*(1.0 + 0.5*a));
  y = (1.0 + a)*z/ke;

  i = potential->maxrp - 1;  
  phase1 = orb->wfun[i];
  r = potential->rad[i-1];  
  b1 = orb->kappa;
  b1 = b1*(b1+1.0) - FINE_STRUCTURE_CONST2*z*z;
 
  a = ke * r;
  b1 = PhaseRDependent(a, y, b1);
  phase1 = phase1 - b1;
  
  orb->phase = phase1;

  return phase1;  
}

int GetNumBounds(const cfac_t *cfac) {
  return cfac->n_orbitals - cfac->n_continua;
}

int GetNumOrbitals(const cfac_t *cfac) {
  return cfac->n_orbitals;
}

int GetNumContinua(const cfac_t *cfac) {
  return cfac->n_continua;
}

int OrbitalIndex(cfac_t *cfac, int n, int kappa, double energy) {
  int i, j;
  ORBITAL *orb;
  int resolve_dirac;

  resolve_dirac = 0;
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (n == 0) {
      if (orb->n == 0 &&
	  orb->kappa == kappa && 
	  orb->energy > 0.0 &&
	  fabs(orb->energy - energy) < EPS10) {
	if (orb->wfun == NULL) {
	  if (RestoreOrbital(i) == 0) return i;
	  else {
	    resolve_dirac = 1;
	    break;
	  }
	}
	return i;
      }
    } else if (orb->n == n && orb->kappa == kappa) {
      if (orb->wfun == NULL) {
	if (RestoreOrbital(i) == 0) return i;
	else {
	  resolve_dirac = 1;
	  break;
	}
      }
      return i;
    }
  }
  
  if (!resolve_dirac) {
    orb = GetNewOrbital(cfac);
  } 

  orb->n = n;
  orb->kappa = kappa;
  orb->energy = energy;
  j = SolveDirac(cfac, orb);
  if (j < 0) {
    printf("Error occured in solving Dirac eq. err = %d\n", j);
    exit(1);
  }
  
  if (n == 0 && !resolve_dirac) {
    cfac->n_continua++;
  }
  return i;
}

int OrbitalExists(const cfac_t *cfac, int n, int kappa, double energy) {
  int i;
  ORBITAL *orb;
  
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (n == 0) {
      if (orb->kappa == kappa &&
	  fabs(orb->energy - energy) < EPS10) 
	return i;
    } else if (orb->n == n && orb->kappa == kappa) {
      return i;
    }
  }
  return -1;
}

int AddOrbital(cfac_t *cfac, ORBITAL *orb) {

  if (orb == NULL) return -1;

  orb = ArrayAppend(cfac->orbitals, orb);
  if (!orb) {
    printf("Not enough memory for cfac->orbitals array\n");
    exit(1);
  }

  if (orb->n == 0) {
    cfac->n_continua++;
  }
  cfac->n_orbitals++;
  return cfac->n_orbitals - 1;
}

ORBITAL *GetOrbital(const cfac_t *cfac, int k) {
  return ArrayGet(cfac->orbitals, k);
}

ORBITAL *GetOrbitalSolved(const cfac_t *cfac, int k) {
  ORBITAL *orb;
  int i;
  
  orb = ArrayGet(cfac->orbitals, k);
  if (orb->wfun == NULL) {
    i = SolveDirac(cfac, orb);
    if (i < 0) {
      printf("Error occured in solving Dirac eq. err = %d\n", i);
      exit(1);
    }
  }
  return orb;
}

ORBITAL *GetNewOrbital(cfac_t *cfac) {
  ORBITAL *orb;

  orb = ArrayAppend(cfac->orbitals, NULL);
  if (!orb) {
    printf("Not enough memory for cfac->orbitals array\n");
    exit(1);
  }

  cfac->n_orbitals++;
  return orb;
}

void FreeOrbitalData(void *p) {
  ORBITAL *orb;

  orb = p;
  if (orb->wfun) free(orb->wfun);
  orb->wfun = NULL;
  orb->phase = 0.0;
  orb->ilast = -1;
}

int ClearOrbitalTable(cfac_t *cfac, int m) {
  ORBITAL *orb;
  int i;

  if (m == 0) {
    cfac->n_orbitals = 0;
    cfac->n_continua = 0;
    ArrayFree(cfac->orbitals);
  } else {
    for (i = cfac->n_orbitals-1; i >= 0; i--) {
      orb = GetOrbital(cfac, i);
      if (orb->n == 0) {
	cfac->n_continua--;
      }
      if (orb->n > 0) {
	cfac->n_orbitals = i+1;
	ArrayTrim(cfac->orbitals, i+1);
	break;
      }
    }
  }
  return 0;
}

int SaveOrbital(int i) {
  return 0;
}

int RestoreOrbital(int i) {
  return -1;
}

int FreeOrbital(cfac_t *cfac, int i) {
  ORBITAL *orb;
  orb = GetOrbital(cfac, i);
  FreeOrbitalData((void *)orb);
  return 0;
}

int SaveAllContinua(cfac_t *cfac, int mode) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (orb->n == 0 && orb->wfun != NULL) {
      if (SaveOrbital(i) < 0) return -1;
      if (mode) {
	FreeOrbital(cfac, i);
      }
    }
  }
  return 0;
}

int SaveContinua(cfac_t *cfac, double e, int mode) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (orb->n == 0 && 
	orb->wfun != NULL &&
	fabs(orb->energy - e) < EPS3) {
      if (SaveOrbital(i) < 0) return -1;
      if (mode) FreeOrbital(cfac, i);
    }
  }
  return 0;
}

int FreeAllContinua(cfac_t *cfac) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (orb->n == 0 && orb->wfun != NULL) {
      FreeOrbital(cfac, i);
    }
  }
  return 0;
}

int FreeContinua(cfac_t *cfac, double e) {
  int i;
  ORBITAL *orb;
  for (i = 0; i < cfac->n_orbitals; i++) {
    orb = GetOrbital(cfac, i);
    if (orb->n == 0 && 
	orb->wfun != NULL &&
	fabs(orb->energy - e) < EPS3) {
      FreeOrbital(cfac, i);
    }
  }
  return 0;
}

int ConfigEnergy(cfac_t *cfac, int m, int mr, int ng, int *kg) {
  CONFIG_GROUP *g;
  CONFIG *cfg;
  int k, i;

  if (m == 0) {
    if (ng == 0) {
      ng = GetNumGroups(cfac);
      for (k = 0; k < ng; k++) {
	OptimizeRadial(cfac, 1, &k, NULL);
	if (mr > 0) RefineRadial(cfac, mr, 0);
	g = GetGroup(cfac, k);
	for (i = 0; i < g->n_cfgs; i++) {
	  cfg = ArrayGet(&(g->cfg_list), i);
	  cfg->energy = AverageEnergyConfig(cfac, cfg);
	}
	ReinitRadial(cfac, 1);
	ClearOrbitalTable(cfac, 0);
      }
    } else {
      OptimizeRadial(cfac, ng, kg, NULL);
      if (mr) RefineRadial(cfac, mr, 0);
      for (k = 0; k < ng; k++) {
	g = GetGroup(cfac, kg[k]);
	for (i = 0; i < g->n_cfgs; i++) {
	  cfg = ArrayGet(&(g->cfg_list), i);
	  if (cfg->energy == 0) {
	    cfg->energy = AverageEnergyConfig(cfac, cfg);
	  }
	}
      }
      ReinitRadial(cfac, 1);
      ClearOrbitalTable(cfac, 0);
    }
  } else {
    ng = GetNumGroups(cfac);
    for (k = 0; k < ng; k++) {
      g = GetGroup(cfac, k);
      for (i = 0; i < g->n_cfgs; i++) {
	cfg = ArrayGet(&(g->cfg_list), i);
	if (cfg->energy != 0) {
	  cfg->delta = cfg->energy - AverageEnergyConfig(cfac, cfg);
	}
      }
    }
  }
  return 0;
}

/* calculate the total configuration average energy of a group. */
double TotalEnergyGroup(cfac_t *cfac, int kg) {
  CONFIG_GROUP *g;
  ARRAY *c;
  CONFIG *cfg;
  int t;
  double total_energy;

  g = GetGroup(cfac, kg);
  c = &(g->cfg_list);
  
  total_energy = 0.0;
  for (t = 0; t < g->n_cfgs; t++) {
    cfg = ArrayGet(c, t);
    total_energy += AverageEnergyConfig(cfac, cfg);
  }
  return total_energy;
}

double ZerothEnergyConfig(cfac_t *cfac, CONFIG *cfg) {
  int i, n, nq, kappa, k;
  double r, e;

  r = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    nq = (cfg->shells[i]).nq;
    kappa = (cfg->shells[i]).kappa;
    k = OrbitalIndex(cfac, n, kappa, 0.0);
    e = GetOrbital(cfac, k)->energy;
    r += nq * e;
  }
  return r;
}

double ZerothResidualConfig(cfac_t *cfac, CONFIG *cfg) {
  int i, n, nq, kappa, k;
  double r, e;

  r = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    nq = (cfg->shells[i]).nq;
    kappa = (cfg->shells[i]).kappa;
    k = OrbitalIndex(cfac, n, kappa, 0.0);
    ResidualPotential(cfac, &e, k, k);
    r += nq * e;
  }
  return r;
}
  
/* calculate the average energy of a configuration */
double AverageEnergyConfig(cfac_t *cfac, CONFIG *cfg) {
  int i, j, n, kappa, nq, np, kappap, nqp;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r, e, am;
  
  am = AMU * cfac_get_atomic_mass(cfac);
 
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = (cfg->shells[i]).n;
    kappa = (cfg->shells[i]).kappa;
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = (cfg->shells[i]).nq;
    k = OrbitalIndex(cfac, n, kappa, 0.0);
    
    if (nq > 1) {
      t = 0.0;

      for (kk = 1; kk <= j2; kk += 1) {
	y = 0;
	if (IsEven(kk)) {
	  Slater(cfac, &y, k, k, k, k, kk, 0);
	}
	if (cfac->qed.br < 0 || n <= cfac->qed.br) {
	  y += Breit(cfac, k, k, k, k, kk, kl, kl, kl, kl);
	}
	if (y) {
	  q = W3j(j2, 2*kk, j2, -1, 0, 1);
	  t += y * q * q ;
	}
      }
      Slater(cfac, &y, k, k, k, k, 0, 0);
      if (cfac->qed.br < 0 || (n > 0 && n <= cfac->qed.br)) {
	y += Breit(cfac, k, k, k, k, 0, kl, kl, kl, kl);
      }
      b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);

    } else {
      b = 0.0;
    }

    t = 0.0;
    for (j = 0; j < i; j++) {
      int maxn;
      double bi;
      
      np = (cfg->shells[j]).n;
      kappap = (cfg->shells[j]).kappa;
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = (cfg->shells[j]).nq;
      kp = OrbitalIndex(cfac, np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      maxn = Max(n, np);
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 2) {
	int kk2 = kk/2;
	y = 0.0;
	if (IsEven((kl+klp+kk)/2)) {
	  Slater(cfac, &y, k, kp, kp, k, kk2, 0);
	  if (kk == 2 && cfac->qed.sms) {
	    double v = Vinti(cfac, k, kp)*Vinti(cfac, kp, k);
	    y -= v/am;
	  }
        }
        if (cfac->qed.br < 0 || maxn <= cfac->qed.br) {
          y += Breit(cfac, k, kp, kp, k, kk2, kl, klp, klp, kl);
        }
	if (y) {
	  q = W3j(j2, kk, j2p, -1, 0, 1);
	  a += y * q * q;
 	}
      }
      y = 0.0;
      Slater(cfac, &y, k, kp, k, kp, 0, 0);
      bi = 0.0;
      if (cfac->qed.br < 0 || maxn <= cfac->qed.br) {
	bi = Breit(cfac, k, kp, k, kp, 0, kl, klp, kl, klp);
	y += bi;
      }      
      t += nqp * (y - a);
    }

    ResidualPotential(cfac, &y, k, k);
    e = GetOrbital(cfac, k)->energy;
    e += QED1E(cfac, k, k);
    r = nq * (b + t + e + y);
    x += r;
  }

  return x;
}

/* calculate the average energy of the average configuration */
double AverageEnergyAvgConfig(cfac_t *cfac) {
  int i, j, n, kappa, np, kappap;
  int k, kp, kk, kl, klp, kkmin, kkmax, j2, j2p;
  double x, y, t, q, a, b, r, nq, nqp, r0, r1;
  AVERAGE_CONFIG *cfg = &cfac->acfg;
 
  r0 = 0.0;
  r1 = 0.0;
  x = 0.0;
  for (i = 0; i < cfg->n_shells; i++) {
    n = cfg->n[i];
    kappa = cfg->kappa[i];
    kl = GetLFromKappa(kappa);
    j2 = GetJFromKappa(kappa);
    nq = cfg->nq[i];
    k = OrbitalIndex(cfac, n, kappa, 0.0);
    
    t = 0.0;
    for (kk = 2; kk <= j2; kk += 2) {
      Slater(cfac, &y, k, k, k, k, kk, 0);
      q = W3j(j2, 2*kk, j2, -1, 0, 1);
      t += y * q * q ;
    }
    Slater(cfac, &y, k, k, k, k, 0, 0);
    b = ((nq-1.0)/2.0) * (y - (1.0 + 1.0/j2)*t);
    
    t = 0.0;
    for (j = 0; j < i; j++) {
      np = cfg->n[j];
      kappap = cfg->kappa[j];
      klp = GetLFromKappa(kappap);
      j2p = GetJFromKappa(kappap);
      nqp = cfg->nq[j];
      kp = OrbitalIndex(cfac, np, kappap, 0.0);

      kkmin = abs(j2 - j2p);
      kkmax = (j2 + j2p);
      if (IsOdd((kkmin + kl + klp)/2)) kkmin += 2;
      a = 0.0;
      for (kk = kkmin; kk <= kkmax; kk += 4) {
	Slater(cfac, &y, k, kp, kp, k, kk/2, 0);
	q = W3j(j2, kk, j2p, -1, 0, 1);
	a += y * q * q;
      }
      Slater(cfac, &y, k, kp, k, kp, 0, 0);

      t += nqp * (y - a);
    }

    ResidualPotential(cfac, &y, k, k);
    a = GetOrbital(cfac, k)->energy;
    r = nq * (b + t + a + y);
    r0 += nq*y;
    r1 += nq*(b+t);
    x += r;
  }

  /*printf("%12.5E %12.5E %15.8E\n", r0, r1, x);*/
  return x;
}

/* calculate the expectation value of the residual potential:
   -Z/r - v0(r), where v0(r) is central potential used to solve 
   dirac equations. the orbital index must be valid, i.e., their 
   radial equations must have been solved. */
int ResidualPotential(const cfac_t *cfac, double *s, int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, z, *p1, *p2, *q1, *q2;
  POTENTIAL *potential = cfac->potential;
  double dwork[MAXRP];

  orb1 = GetOrbitalSolved(cfac, k0);
  orb2 = GetOrbitalSolved(cfac, k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    *s = 0.0;
    return 0;
  }

  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  p = (double *) MultiSet(cfac->residual_array, index, NULL);
  if (p && *p) {
    *s = *p;
    return 0;
  } 

  *s = 0.0;
 
  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      dwork[i] = potential->Vn[i] - z;
      dwork[i] *= potential->dr_drho[i];
      dwork[i] *= p1[i]*p2[i] + q1[i]*q2[i];
    }
    *s = Simpson(dwork, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->U[i];
      z += potential->Vc[i];
      dwork[i] = potential->Vn[i] - z;
    }
    IntegrateS(potential, dwork, orb1, orb2, INT_P1P2pQ1Q2, s, -1);
  }
  *p = *s;
  return 0;
}

double MeanPotential(cfac_t *cfac, int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  double z, *p1, *p2, *q1, *q2;
  POTENTIAL *potential = cfac->potential;
  double dwork[MAXRP];

  orb1 = GetOrbitalSolved(cfac, k0);
  orb2 = GetOrbitalSolved(cfac, k1);
  if (!orb1 || !orb2) return -1;
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (orb1->n < 0 || orb2->n < 0) {
    p1 = Large(orb1);
    p2 = Large(orb2);
    q1 = Small(orb1);
    q2 = Small(orb2);
    for (i = potential->ib; i <= potential->ib1; i++) {
      z = potential->U[i] + potential->Vc[i];
      dwork[i] = z*potential->dr_drho[i]*(p1[i]*p2[i] + q1[i]*q2[i]);
    }
    z = Simpson(dwork, potential->ib, potential->ib1);
  } else {
    for (i = 0; i < potential->maxrp; i++) {
      z = potential->U[i] + potential->Vc[i];
      dwork[i] = z;
    }
    IntegrateS(potential, dwork, orb1, orb2, INT_P1P2pQ1Q2, &z, -1);
  }

  return z;
}

/* \int R(k1)R(k2)r^m */
double RadialMoments(const cfac_t *cfac, int m, int k1, int k2) {
  int index[3];
  int npts, i0, i;
  ORBITAL *orb1, *orb2;
  double *q, r, *p1, *p2, *q1, *q2;
  int n1, n2;
  int kl1, kl2;
  int nh, klh;
  POTENTIAL *potential = cfac->potential;
  double yk[MAXRP];
  
  orb1 = GetOrbitalSolved(cfac, k1);
  orb2 = GetOrbitalSolved(cfac, k2);
  n1 = orb1->n;
  n2 = orb2->n;
  kl1 = GetLFromKappa(orb1->kappa)/2;
  kl2 = GetLFromKappa(orb2->kappa)/2;

  GetHydrogenicNL(cfac, &nh, &klh, NULL, NULL);

  if (n1 > 0 && n2 > 0 && potential->ib <= 0) {
    if ((n1 > nh && n2 > nh) || 
	(kl1 > klh && kl2 > klh) ||
	orb1->wfun == NULL || 
	orb2->wfun == NULL) {
      if (n1 == n2 && kl1 == kl2) {
	double z = GetResidualZ(cfac);
	r = HydrogenicExpectation(z, m, n1, kl1);
	if (r) {
	  return r;
	}
      } else if (m == 1) {
	r = HydrogenicDipole(cfac, n1, kl1, n2, kl2);
	return r;
      }
    }
  }

  if (potential->ib <= 0 && n1 == n2 && m > 0 && n1 > GetNMax(potential)) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }

  if (m >= 0) {
    index[0] = 2*m;
  } else {
    index[0] = -2*m-1;
  }
    
  if (k1 < k2) {
    index[1] = k1;
    index[2] = k2;
  } else {
    index[1] = k2;
    index[2] = k1;
  }
  
  q = (double *) MultiSet(cfac->moments_array, index, NULL);
 
  if (*q) {
    return *q;
  } 

  if (n1 < 0 || n2 < 0) {
    i0 = potential->ib;
    npts = potential->ib1;
    p1 = Large(orb1);
    q1 = Small(orb1);
    p2 = Large(orb2);
    q2 = Small(orb2);
    for (i = i0; i <= npts; i++) {
      r = p1[i]*p2[i] + q1[i]*q2[i];
      r *= potential->dr_drho[i];
      yk[i] = pow(potential->rad[i], m)*r;
    }
    r = Simpson(yk, i0, npts);
  } else {    
    npts = potential->maxrp-1;
    if (n1 != 0) npts = Min(npts, orb1->ilast);
    if (n2 != 0) npts = Min(npts, orb2->ilast);

    for (i = 0; i <= npts; i++) {
      yk[i] = pow(potential->rad[i], m);
    }
    r = 0.0;
    IntegrateS(potential, yk, orb1, orb2, INT_P1P2pQ1Q2, &r, m);
  }
  
  *q = r;
  return r;
}

double MultipoleRadialNR(cfac_t *cfac, int m, int k1, int k2, int gauge) {
  int i, p, t;
  ORBITAL *orb1, *orb2;
  double r;
  int kappa1, kappa2;

  orb1 = GetOrbital(cfac, k1);
  orb2 = GetOrbital(cfac, k2);
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;

  r = 0.0;
  if (m == 1) {
    /* use the relativistic version is just simpler */
    printf("should call MultipoleRadialFR instead\n");
  } else if (m > 1) {
    t = kappa1 + kappa2;
    p = m - t;
    if (p && t) {
      r = RadialMoments(cfac, m-1, k1, k2);
      r *= p*t;
      r /= sqrt(m*(m+1.0));
      r *= -0.5 * FINE_STRUCTURE_CONST;
      for (i = 2*m - 1; i > 0; i -= 2) r /= i;
      r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
    }
  } else if (m < 0) {
    m = -m;
    if (gauge != G_BABUSHKIN) {
      printf("the velocity form is not implemented yet\n");
      abort();
    }

    r = RadialMoments(cfac, m, k1, k2);
    r *= sqrt((m+1.0)/m);
    for (i = 2*m - 1; i > 1; i -= 2) r /= i;
    r *= ReducedCL(GetJFromKappa(kappa1), 2*m, GetJFromKappa(kappa2));
  }

  return r;
}

int MultipoleRadialFRGrid(cfac_t *cfac,
    double **p0, int m, int k1, int k2, int gauge) {
  double q, ip, ipm, im, imm;
  int kappa1, kappa2;
  int am, t;
  int index[4];
  ORBITAL *orb1, *orb2;
  double x, a, r, rp, ef, **p1;
  int n, i, j, npts;
  double rcl;
  POTENTIAL *potential = cfac->potential;
  double dwork1[MAXRP], dwork2[MAXRP];

  if (m == 0) return 0;
  
  if (m >= 0) {
    index[0] = 2*m;
    am = m;
  } else {
    index[0] = -2*m-1;
    am = -m;
  }
  index[1] = k1;
  index[2] = k2;

  p1 = MultiSet(cfac->multipole_array, index, NULL);
  if (*p1) {
    *p0 = *p1;
    return cfac->n_awgrid;
  }

  orb1 = GetOrbitalSolved(cfac, k1);
  orb2 = GetOrbitalSolved(cfac, k2);
  
  kappa1 = orb1->kappa;
  kappa2 = orb2->kappa;
  rcl = ReducedCL(GetJFromKappa(kappa1), abs(2*m), 
		  GetJFromKappa(kappa2));

  ef = Max(orb1->energy, orb2->energy);
  if (ef > 0.0) {
    ef *= FINE_STRUCTURE_CONST;
  } else {
    ef = 0.0;
  }

  *p1 = malloc(sizeof(double)*cfac->n_awgrid);
  
  npts = potential->maxrp-1;
  if (orb1->n > 0) npts = Min(npts, orb1->ilast);
  if (orb2->n > 0) npts = Min(npts, orb2->ilast);
  r = 0.0;

  for (i = 0; i < cfac->n_awgrid; i++) {
    r = 0.0;
    a = cfac->awgrid[i];
    (*p1)[i] = 0.0;
    if (ef > 0.0) a += ef;
    if (m > 0) {
      t = kappa1 + kappa2;
      if (t) {
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = m;
	  dwork1[j] = gsl_sf_bessel_jl(n, x);
	}
	IntegrateS(potential, dwork1, orb1, orb2, INT_P1Q2pQ1P2, &r, 0);
	r *= t;
	r *= (2*m + 1.0)/sqrt(m*(m+1.0));
	r /= pow(a, m);
	(*p1)[i] = r*rcl;
      }
    } else {
      if (gauge == G_COULOMB) {
	t = kappa1 - kappa2;
	q = sqrt(am/(am+1.0));
	for (j = 0; j <= npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  dwork1[j] = gsl_sf_bessel_jl(n, x);
	  n = am-1;
	  dwork2[j] = gsl_sf_bessel_jl(n, x);
	}
	r = 0.0;
	rp = 0.0;
	if (t) {
	  IntegrateS(potential, dwork1, orb1, orb2, INT_P1Q2pQ1P2, &ip, 0);
	  IntegrateS(potential, dwork2, orb1, orb2, INT_P1Q2pQ1P2, &ipm, 0);
	  r = t*ip*q - t*ipm/q;
	}
	if (k1 != k2) {
	  IntegrateS(potential, dwork1, orb1, orb2, INT_P1Q2mQ1P2, &im, 0);
	  IntegrateS(potential, dwork2, orb1, orb2, INT_P1Q2mQ1P2, &imm, 0);
	  rp = (am + 1.0)*im*q + am*imm/q;
	}
	r += rp;
	if (am > 1) r /= pow(a, am-1);
	(*p1)[i] = r*rcl;
      } else if (gauge == G_BABUSHKIN) {
	t = kappa1 - kappa2;
	for (j = 0; j < npts; j++) {
	  x = a*potential->rad[j];
	  n = am+1;
	  dwork1[j] = gsl_sf_bessel_jl(n, x);
	  n = am;
	  dwork2[j] = gsl_sf_bessel_jl(n, x);
	}
	if (t) {
	  IntegrateS(potential, dwork1, orb1, orb2, INT_P1Q2pQ1P2, &ip, 0);
	  r = t*ip;
	}
	if (k1 != k2) {
	  IntegrateS(potential, dwork1, orb1, orb2, INT_P1Q2mQ1P2, &im, 0);
	} else {
	  im = 0.0;
	}
	IntegrateS(potential, dwork2, orb1, orb2, INT_P1P2pQ1Q2, &imm, 0);
	rp = (am + 1.0) * (imm + im);
	q = (2*am + 1.0)/sqrt(am*(am+1.0));
	q /= pow(a, am);
	r *= q;
	rp *= q;
	(*p1)[i] = (r+rp)*rcl;
      }
    }
  }

  *p0 = *p1;
  return cfac->n_awgrid;
}

double MultipoleRadialFR(cfac_t *cfac,
    double aw, int m, int k1, int k2, int gauge) {
  int n;
  ORBITAL *orb1, *orb2;
  double *y, ef, r;

  orb1 = GetOrbitalSolved(cfac, k1);
  orb2 = GetOrbitalSolved(cfac, k2);
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    if (m == -1) {
      return MultipoleRadialNR(cfac, m, k1, k2, gauge);
    } else {
      return 0.0;
    }
  }
  
  n = MultipoleRadialFRGrid(cfac, &y, m, k1, k2, gauge);
  if (n == 0) return 0.0;
  
  ef = Max(orb1->energy, orb2->energy);
  if (ef > 0.0) {
    ef *= FINE_STRUCTURE_CONST;
  } else {
    ef = 0.0;
  }
  if (cfac->n_awgrid > 1) {
    if (ef > 0) aw += ef;
  }
  
  r = InterpolateMultipole(aw, n, cfac->awgrid, y);
  if (gauge == G_COULOMB && m < 0) r /= aw;

  return r;
}

double *GeneralizedMoments(cfac_t *cfac, int k1, int k2, int m) {
  ORBITAL *orb1, *orb2;
  int n1, i, nk;
  double x, r, r0;
  double *p1, *p2, *q1, *q2;
  int index[3], t;
  double **p, k, *kg;
  double amin, amax, kmin, kmax;
  double _phase[MAXRP];
  double _dphase[MAXRP];
  POTENTIAL *potential = cfac->potential;
  double yk[MAXRP];
  
  index[0] = m;
  if (k1 > k2) {
    index[1] = k2;
    index[2] = k1;
    orb1 = GetOrbitalSolved(cfac, k2);
    orb2 = GetOrbitalSolved(cfac, k1);
  } else {
    index[1] = k1;
    index[2] = k2;
    orb1 = GetOrbitalSolved(cfac, k1);
    orb2 = GetOrbitalSolved(cfac, k2);
  }

  p = (double **) MultiSet(cfac->gos_array, index, NULL);
  if (*p) {
    return *p;
  }

  nk = NGOSK;
  *p = malloc(sizeof(double)*nk*2);
  kg = *p + nk;

  if (orb1->wfun == NULL || orb2->wfun == NULL || 
      (orb1->n <= 0 && orb2->n <= 0)) {
    for (t = 0; t < nk*2; t++) {
      (*p)[t] = 0.0;
    }
    return *p;
  }
  
  p1 = Large(orb1);
  p2 = Large(orb2);
  q1 = Small(orb1);
  q2 = Small(orb2);

  amin = sqrt(2.0*fabs(orb1->energy));
  amax = sqrt(2.0*fabs(orb2->energy));
  if (amin < amax) {
    kmin = amin;
    kmax = amax;
  } else {
    kmin = amax;
    kmax = amin;
  }
  kmin = log(0.05*kmin);
  kmax = log(50.0*kmax);
  r = (kmax - kmin)/(nk-1.0);
  kg[0] = kmin;
  for (i = 1; i < nk; i++) {
    kg[i] = kg[i-1] + r;
  }
  
  if (orb1->n > 0 && orb2->n > 0) {
    n1 = Min(orb1->ilast, orb2->ilast);
  
    for (i = 0; i <= n1; i++) {
      _phase[i] = (p1[i]*p2[i] + q1[i]*q2[i])*potential->dr_drho[i];
    }
    
    if (m == 0) {
      if (k1 == k2) r0 = 1.0;
      else if (orb1->n != orb2->n) r0 = 0.0;
      else {
	if (orb1->kappa + orb2->kappa != -1) r0 = 0.0;
	else {
	  r0 = Simpson(_phase, 0, n1);
	}
      }
    } else {
      r0 = 0.0;
    }
    
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	_dphase[i] = gsl_sf_bessel_jl(m, x);
	_dphase[i] *= _phase[i];
      }
      r = Simpson(_dphase, 0, n1);
      
      (*p)[t] = (r - r0)/k;
    }
  } else {
    if (orb1->n > 0) n1 = orb1->ilast;
    else n1 = orb2->ilast;
    for (t = 0; t < nk; t++) {
      k = exp(kg[t]);      
      for (i = 0; i <= n1; i++) {
	x = k * potential->rad[i];
	yk[i] = gsl_sf_bessel_jl(m, x);
      }
      IntegrateS(potential, yk, orb1, orb2, INT_P1P2pQ1Q2, &r, 0);
      (*p)[t] = r/k;
    }
  }
  return *p;
}

void PrintGeneralizedMoments(cfac_t *cfac, char *fn, int m, int n0, int k0, 
			     int n1, int k1, double e1) {
  FILE *f;
  int i0, i1, i;
  double *g, *x;

  i0 = OrbitalIndex(cfac, n0, k0, 0.0);
  e1 /= HARTREE_EV;
  i1 = OrbitalIndex(cfac, n1, k1, e1);
  f = fopen(fn, "w");
  if (f == NULL) {
    printf("cannot open file %s\n", fn);
    return;
  }
  g = GeneralizedMoments(cfac, i0, i1, m);
  x = g + NGOSK;
  for (i = 0; i < NGOSK; i++) {
    fprintf(f, "%15.8E %15.8E %15.8E\n", x[i], exp(x[i]), g[i]);
  }
  fclose(f);
}
  
double InterpolateMultipole(double aw, int n, double *x, double *y) {
  double r;
  int nd;

  if (n == 1) {
    r = y[0];
  } else {
    nd = 1;
    UVIP3P(n, x, y, nd, &aw, &r);
  }

  return r;
}

int SlaterTotal(cfac_t *cfac,
    double *sd, double *se, int *j, int *ks, int k, int mode) {
  int t, kk, tt, maxn;
  int tmin, tmax;
  double e, a, d, a1, a2, am;
  int kl0, kl1, kl2, kl3;
  int k0, k1, k2, k3;
  int js[4];
  ORBITAL *orb0, *orb1, *orb2, *orb3;

  k0 = ks[0];
  k1 = ks[1];
  k2 = ks[2];
  k3 = ks[3];
  kk = k/2;

  maxn = 0;
  orb0 = GetOrbitalSolved(cfac, k0);
  orb1 = GetOrbitalSolved(cfac, k1);
  orb2 = GetOrbitalSolved(cfac, k2);
  orb3 = GetOrbitalSolved(cfac, k3);

  if (orb0->n <= 0) {
    maxn = -1;
  } else if (orb0->n > maxn) {
    maxn = orb0->n;
    if (orb1->n <= 0) {
      maxn = -1;
    } else if (orb1->n > maxn) {
      maxn = orb1->n;
      if (orb2->n <= 0) {
	maxn = -1;
      } else if (orb2->n > maxn) {
	maxn = orb2->n;
	if (orb3->n <= 0) {
	  maxn = -1;
	} else if (orb3->n > maxn) {
	  maxn = orb3->n;
	}
      }
    }
  }

  if (orb0->wfun == NULL || orb1->wfun == NULL ||
      orb2->wfun == NULL || orb3->wfun == NULL) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  kl0 = GetLFromKappa(orb0->kappa);
  kl1 = GetLFromKappa(orb1->kappa);
  kl2 = GetLFromKappa(orb2->kappa);
  kl3 = GetLFromKappa(orb3->kappa);

  if (orb1->n < 0 || orb3->n < 0) {
    mode = 2;
  } else {
    if (kl1 > cfac->slater_cut.kl0 && kl3 > cfac->slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }
    if (kl0 > cfac->slater_cut.kl0 && kl2 > cfac->slater_cut.kl0) {
      if (se) {
	*se = 0.0;
	se = NULL;
      }
    }  
    if (kl1 > cfac->slater_cut.kl1 && kl3 > cfac->slater_cut.kl1) {
      mode = 2;
    }
    if (kl0 > cfac->slater_cut.kl1 && kl2 > cfac->slater_cut.kl1) {
      mode = 2;
    }
  }
  if (cfac->qed.br == 0 && IsOdd((kl0+kl1+kl2+kl3)/2)) {
    if (sd) *sd = 0.0;
    if (se) *se = 0.0;
    return 0;
  }

  if (j) {
    memcpy(js, j, sizeof(int)*4);
  } else {
    js[0] = 0;
    js[1] = 0;
    js[2] = 0;
    js[3] = 0;
  }

  if (js[0] <= 0) js[0] = GetJFromKappa(orb0->kappa);
  if (js[1] <= 0) js[1] = GetJFromKappa(orb1->kappa);
  if (js[2] <= 0) js[2] = GetJFromKappa(orb2->kappa);
  if (js[3] <= 0) js[3] = GetJFromKappa(orb3->kappa);  

  am = AMU * cfac_get_atomic_mass(cfac);
  if (sd) {
    d = 0.0;
    if (Triangle(js[0], js[2], k) && Triangle(js[1], js[3], k)) {
      if (IsEven((kl0+kl2)/2+kk) && IsEven((kl1+kl3)/2+kk)) {	
	Slater(cfac, &d, k0, k1, k2, k3, kk, mode);
	if (kk == 1 && cfac->qed.sms && maxn > 0) {
	  a1 = Vinti(cfac, k0, k2);
	  a2 = Vinti(cfac, k1, k3);
	  d -= a1 * a2 / am;
	}
      }

      if (cfac->qed.br < 0 || (maxn > 0 && maxn <= cfac->qed.br)) {
        d += Breit(cfac, k0, k1, k2, k3, kk, kl0, kl1, kl2, kl3);
      }
      if (d) {
	a1 = ReducedCL(js[0], k, js[2]);
	a2 = ReducedCL(js[1], k, js[3]);
	d *= a1*a2;
	if (k0 == k1 && k2 == k3) d *= 0.5;
      }
    }
    *sd = d;
  }
  
  if (!se) goto EXIT;

  if (abs(mode) == 2) {
    *se = 0.0;
    goto EXIT;
  }
  *se = 0.0;
  if (k0 == k1 && (orb0->n > 0 || orb1->n > 0)) goto EXIT;
  if (k2 == k3 && (orb2->n > 0 || orb3->n > 0)) goto EXIT;
  tmin = abs(js[0] - js[3]);
  tt = abs(js[1] - js[2]);
  tmin = Max(tt, tmin);
  tmax = js[0] + js[3];
  tt = js[1] + js[2];
  tmax = Min(tt, tmax);
  tmax = Min(tmax, GetMaxRank(cfac));
  if (IsOdd(tmin)) tmin++;
  
  for (t = tmin; t <= tmax; t += 2) {
    a = W6j(js[0], js[2], k, js[1], js[3], t);
    if (fabs(a) > EPS30) {
      e = 0.0;
      if (IsEven((kl0+kl3+t)/2) && IsEven((kl1+kl2+t)/2)) {
	Slater(cfac, &e, k0, k1, k3, k2, t/2, mode);
	if (t == 2 && cfac->qed.sms && maxn > 0) {
	  e -= Vinti(cfac, k0, k3) * Vinti(cfac, k1, k2) / am;
	}
      }
      if (cfac->qed.br < 0 || (maxn > 0 && maxn <= cfac->qed.br)) {
	e += Breit(cfac, k0, k1, k3, k2, t/2, kl0, kl1, kl3, kl2);
      }
      if (e) {
	e *= ReducedCL(js[0], t, js[3]); 
	e *= ReducedCL(js[1], t, js[2]);
	e *= a * (k + 1.0);
	if (IsOdd(t/2 + kk)) e = -e;
	*se += e;
      }
    }
  }

 EXIT:
  return 0;
}

double SelfEnergyRatio(POTENTIAL *potential, ORBITAL *orb) {
  int i, npts;
  double p[MAXRP], q[MAXRP], z;
  double *large, *small;
  double a, b;
  
  if (orb->wfun == NULL) return 1.0;

  for (npts = 0; npts < potential->maxrp; npts++) {
    if (potential->uehling[npts] > -EPS4) break;
  }
  
  z = potential->anum;
  RadialDiracCoulomb(npts, p, q, potential->rad, z, 
			 orb->n, orb->kappa);
  large = Large(orb);
  small = Small(orb);  
  for (i = 0; i < npts; i++) {
    p[i] = (p[i]*p[i] + q[i]*q[i])*potential->dr_drho[i];
    p[i] *= potential->uehling[i];
    q[i] = (large[i]*large[i] + small[i]*small[i])*potential->dr_drho[i];
    q[i] *= potential->uehling[i];
  }
  a = Simpson(p, 0, npts-1);
  b = Simpson(q, 0, npts-1);

  return b/a;
}

double QED1E(cfac_t *cfac, int k0, int k1) {
  POTENTIAL *potential = cfac->potential;
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, r, a;
  double dwork[MAXRP];

  if (cfac->qed.nms == 0 && cfac->qed.vp == 0) {
    if (cfac->qed.se == 0 || k0 != k1) {
      return 0.0;
    }
  }

  orb1 = GetOrbitalSolved(cfac, k0);
  orb2 = GetOrbitalSolved(cfac, k1);

  if (orb1->n <= 0 || orb2->n <= 0) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }
  
  if (k0 > k1) {
    index[0] = k1;
    index[1] = k0;
  } else {
    index[0] = k0;
    index[1] = k1;
  }
  
  p = (double *) MultiSet(cfac->qed1e_array, index, NULL);
  if (p && *p) {
    return *p;
  }

  r = 0.0;
  
  if (cfac->qed.nms > 0) {
    for (i = 0; i < potential->maxrp; i++) {
      dwork[i] = potential->U[i] + potential->Vc[i];
    }
    a = 0.0;
    IntegrateS(potential, dwork, orb1, orb2, INT_P1P2pQ1Q2, &a, -1);
    a = -a;
    if (k0 == k1) a += orb1->energy;
    a /= (AMU * cfac_get_atomic_mass(cfac));
    r += a;
    for (i = 0; i < potential->maxrp; i++) {
      dwork[i] = (orb1->energy - (potential->U[i] + potential->Vc[i]))*
                 (orb2->energy - (potential->U[i] + potential->Vc[i]));
    }
    IntegrateS(potential, dwork, orb1, orb2, INT_P1P2pQ1Q2, &a, -1);
    a *= FINE_STRUCTURE_CONST2/(2.0 * AMU * cfac_get_atomic_mass(cfac));
    r += a;
  }

  if (cfac->qed.vp > 0) {
    a = 0.0;
    IntegrateS(potential, potential->uehling, orb1, orb2, INT_P1P2pQ1Q2, &a, 0);
    r += a;
  }

  if (k0 == k1 && (cfac->qed.se < 0 || orb1->n <= cfac->qed.se)) {
    if (potential->ib <= 0 || orb1->n <= potential->nb) {
      a = HydrogenicSelfEnergy(cfac_get_atomic_number(cfac), 
			       orb1->n, orb1->kappa);
      if (a) {
	a *= SelfEnergyRatio(potential, orb1);
	r += a;
      }
    }
  }
  *p = r;
  return r;
}
  
double Vinti(cfac_t *cfac, int k0, int k1) {
  int i;
  ORBITAL *orb1, *orb2;
  int index[2];
  double *p, *large0, *small0, *large1, *small1;
  int ka0, ka1;
  double a, b, r;
  POTENTIAL *potential = cfac->potential;
  double dwork[MAXRP];

  orb1 = GetOrbitalSolved(cfac, k0);
  orb2 = GetOrbitalSolved(cfac, k1);

  if (orb1->n <= 0 || orb2->n <= 0) {
    return 0.0;
  }
  if (orb1->wfun == NULL || orb2->wfun == NULL) {
    return 0.0;
  }
  
  index[0] = k0;
  index[1] = k1;
  
  p = (double *) MultiSet(cfac->vinti_array, index, NULL);
  if (p && *p) {
    return *p;
  }

  ka0 = orb1->kappa;
  ka1 = orb2->kappa;
  large0 = Large(orb1);
  large1 = Large(orb2);
  small0 = Small(orb1);
  small1 = Small(orb2);
  a = 0.5*(ka0*(ka0+1.0) - ka1*(ka1+1.0));
  b = 0.5*(-ka0*(-ka0+1.0) + ka1*(-ka1+1.0));
  r = 0.0;

  Differential(large1, dwork, 0, potential->maxrp-1);
  for (i = 0; i < potential->maxrp; i++) {
    dwork[i] = large0[i]*dwork[i] - a*large0[i]*large1[i]/potential->rad[i];
    dwork[i] *= potential->dr_drho[i];
  }
  r += Simpson(dwork, 0, potential->maxrp-1);
  
  Differential(small1, dwork, 0, potential->maxrp-1);
  for (i = 0; i < potential->maxrp; i++) {
    dwork[i] = small0[i]*dwork[i] - b*small0[i]*small1[i]/potential->rad[i];
    dwork[i] *= potential->dr_drho[i];
  }
  r += Simpson(dwork, 0, potential->maxrp-1);
  
  *p = r;

  return r;
}

double BreitC(cfac_t *cfac, int n, int m, int k, int k0, int k1, int k2, int k3) {
  int ka0, ka1, ka2, ka3, kb, kp;
  double r, b, c;
  
  ka0 = GetOrbital(cfac, k0)->kappa;
  ka1 = GetOrbital(cfac, k1)->kappa;
  ka2 = GetOrbital(cfac, k2)->kappa;
  ka3 = GetOrbital(cfac, k3)->kappa;
  if (k == m) {
    r = -(ka0 + ka2)*(ka1 + ka3);
    if (r) r /= (m*(m+1.0));
  } else if (k == (m + 1)) {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m + 2.0)/(2.0*(2.0*m + 1.0));
    c = -(m - 1.0)/((2.0*m+1.0)*(2.0*m+2.0));
    switch (n) {
    case 0:
      r = (k + kb)*(b + c*kp);
      break;
    case 1:
      r = (k + kp)*(b + c*kb);
      break;
    case 2:
      r = (k - kb)*(b - c*kp);
      break;
    case 3:
      r = (k - kp)*(b - c*kb);
      break;
    case 4:
      r = -(k + kb)*(b - c*kp);
      break;
    case 5:
      r = -(k - kp)*(b + c*kb);
      break;
    case 6:
      r = -(k - kb)*(b + c*kp);
      break;
    case 7:
      r = -(k + kp)*(b - c*kb);
      break;
    default:
      r = 0;
    }
  } else {
    kb = ka2 - ka0;
    kp = ka3 - ka1;
    b = (m - 1.0)/(2.0*(2.0*m + 1.0));
    c = (m + 2.0)/(2.0*m*(2.0*m + 1.0));
    switch (n) {
    case 0:
      r = (b + c*kb)*(kp - k - 1.0);
      break;
    case 1:
      r = (b + c*kp)*(kb - k - 1.0);
      break;
    case 2:
      r = (b - c*kb)*(-kp - k - 1.0);
      break;
    case 3:
      r = (b - c*kp)*(-kb - k - 1.0);
      break;
    case 4:
      r = -(b + c*kb)*(-kp - k - 1.0);
      break;
    case 5:
      r = -(b - c*kp)*(kb - k - 1.0);
      break;
    case 6:
      r = -(b - c*kb)*(kp - k - 1.0);
      break;
    case 7:
      r = -(b + c*kp)*(-kb - k - 1.0);
      break;
    default:
      r = 0;
    }
  }

  return r;
}

double BreitS(cfac_t *cfac, int k0, int k1, int k2, int k3, int k) {
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  int index[5], i;
  double *p, r;
  double dwork1[MAXRP], dwork2[MAXRP];
  POTENTIAL *potential = cfac->potential;
  
  index[0] = k0;
  index[1] = k1;
  index[2] = k2;
  index[3] = k3;
  index[4] = k;

  p = (double *) MultiSet(cfac->breit_array, index, NULL);
  if (p && *p) {
    r = *p;
  } else {
    orb0 = GetOrbitalSolved(cfac, k0);
    orb1 = GetOrbitalSolved(cfac, k1);
    orb2 = GetOrbitalSolved(cfac, k2);
    orb3 = GetOrbitalSolved(cfac, k3);
    if (!orb0 || !orb1 || !orb2 || !orb3) return 0.0;
    
    for (i = 0; i < potential->maxrp; i++) {
      dwork1[i] = pow(potential->rad[i], k);
    }
    
    IntegrateF(potential, dwork1, orb0, orb1, INT_P1Q2, dwork2, 0);
    
    for (i = 0; i < potential->maxrp; i++) {
      dwork2[i] /= dwork1[i]*potential->rad[i];
    }

    IntegrateS(potential, dwork2, orb2, orb3, INT_P1Q2, &r, 0);
    *p = r;
  }

  return r;
}

double BreitI(cfac_t *cfac, int n, int k0, int k1, int k2, int k3, int m) {
  double r;

  switch (n) {
  case 0:
    r = BreitS(cfac, k0, k2, k1, k3, m);
    break;
  case 1:
    r = BreitS(cfac, k1, k3, k0, k2, m);
    break;
  case 2:
    r = BreitS(cfac, k2, k0, k3, k1, m);
    break;
  case 3:
    r = BreitS(cfac, k3, k1, k2, k0, m);
    break;
  case 4:
    r = BreitS(cfac, k0, k2, k3, k1, m);
    break;
  case 5:
    r = BreitS(cfac, k3, k1, k0, k2, m);
    break;
  case 6:
    r = BreitS(cfac, k2, k0, k1, k3, m);
    break;
  case 7:
    r = BreitS(cfac, k1, k3, k2, k0, m);
    break;
  default:
    r = 0.0;
  }

  return r;
}

double Breit(cfac_t *cfac, int k0, int k1, int k2, int k3, int k,
	     int kl0, int kl1, int kl2, int kl3) {
  int m, m0, m1, n;
  double a, c, r;
  
  m0 = k - 1;
  if (m0 < 0) m0 = 0;
  m1 = k + 1;
  r = 0.0;
  for (m = m0; m <= m1; m++) {
    if (IsEven((kl0+kl2)/2 + m) || IsEven((kl1+kl3)/2 + m)) continue;
    for (n = 0; n < 8; n++) {
      c = BreitC(cfac, n, m, k, k0, k1, k2, k3);
      a = BreitI(cfac, n, k0, k1, k2, k3, m);
      r += a*c;
    }
  }

  return r;
}

/* calculate the slater integral of rank k */
int Slater(const cfac_t *cfac,
    double *s, int k0, int k1, int k2, int k3, int k, int mode) {
  int index[5];
  double *p;
  int ilast, i, npts, m;
  ORBITAL *orb0, *orb1, *orb2, *orb3;
  double norm;
  POTENTIAL *potential = cfac->potential;
  double yk[MAXRP];

  index[0] = k0;
  index[1] = k1;
  index[2] = k2;
  index[3] = k3;
  index[4] = k;  

  if (abs(mode) < 2) {
    SortSlaterKey(index);
    p = (double *) MultiSet(cfac->slater_array, index, NULL);
  } else {
    p = NULL;
  }
  if (p && *p) {
    *s = *p;
  } else {
    orb0 = GetOrbitalSolved(cfac, k0);
    orb1 = GetOrbitalSolved(cfac, k1);
    orb2 = GetOrbitalSolved(cfac, k2);
    orb3 = GetOrbitalSolved(cfac, k3);
    *s = 0.0;
    if (!orb0 || !orb1 || !orb2 || !orb3) return -1;  

    npts = potential->maxrp;
    switch (mode) {
    case 0: /* fall through to case 1 */
    case 1: /* full relativistic with distorted free cfac->orbitals */
      GetYk(cfac, k, yk, orb0, orb2, k0, k2, 1); 
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	yk[i] = (yk[i]/potential->rad[i]);
      }
      IntegrateS(potential, yk, orb1, orb3, INT_P1P2pQ1Q2, s, 0);
      break;
    
    case -1: /* quasi relativistic with distorted free cfac->orbitals */
      GetYk(cfac, k, yk, orb0, orb2, k0, k2, 2);
      if (orb1->n > 0) ilast = orb1->ilast;
      else ilast = npts-1;
      if (orb3->n > 0) ilast = Min(ilast, orb3->ilast);
      for (i = 0; i <= ilast; i++) {
	yk[i] /= potential->rad[i];
      }
      IntegrateS(potential, yk, orb1, orb3, INT_P1P2, s, 0);

      norm  = orb0->qr_norm;
      norm *= orb1->qr_norm;
      norm *= orb2->qr_norm;
      norm *= orb3->qr_norm;

      *s *= norm;
      break;

    case 2: /* separable coulomb interaction, orb0, orb2 is inner part */
      m = k;
      *s = RadialMoments(cfac, m, k0, k2);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(cfac, m, k1, k3);
      }
      break;

    case -2: /* separable coulomb interaction, orb1, orb3 is inner part  */
      m = k;
      *s = RadialMoments(cfac, m, k1, k3);
      if (*s != 0.0) {
	m = -m-1;
	*s *= RadialMoments(cfac, m, k0, k2);      
      }
      break;
      
    default:
      break;
    }      

    if (p) *p = *s;
  }

  return 0;
}


/* reorder the orbital index appears in the slater integral, so that it is
   in a form: a <= b <= d, a <= c, and if (a == b), c <= d. */ 
void SortSlaterKey(int *kd) {
  int i;

  if (kd[0] > kd[2]) {
    i = kd[0];
    kd[0] = kd[2];
    kd[2] = i;
  }

  if (kd[1] > kd[3]) {
    i = kd[1];
    kd[1] = kd[3];
    kd[3] = i;
  }
  
  if (kd[0] > kd[1]) {
    i = kd[0];
    kd[0] = kd[1];
    kd[1] = i;
    i = kd[2];
    kd[2] = kd[3];
    kd[3] = i;
  } else if (kd[0] == kd[1]) {
    if (kd[2] > kd[3]) {
      i = kd[2];
      kd[2] = kd[3];
      kd[3] = i;
    }
  }
}

/*
** this is a better version of Yk than GetYk0.
** note that on exit, rk contains r^k, which is used in GetYk
*/      
static int GetYk1(POTENTIAL *potential,
    int k, double *yk, double *rk,
    const ORBITAL *orb1, const ORBITAL *orb2, int type) {
  int i, ilast;
  double r0, a;
  double dwork1[MAXRP];
  double dwork2[MAXRP];
  
  ilast = Min(orb1->ilast, orb2->ilast);
  r0 = sqrt(potential->rad[0]*potential->rad[ilast]);  
  for (i = 0; i < potential->maxrp; i++) {
    dwork1[i] = pow(potential->rad[i]/r0, k);
  }
  IntegrateF(potential, dwork1, orb1, orb2, type, dwork2, 0);
  a = pow(r0, k);
  for (i = 0; i < potential->maxrp; i++) {
    yk[i] = dwork2[i]/dwork1[i];
    rk[i] = dwork1[i]*a;
  }  
  for (i = 0; i < potential->maxrp; i++) {
    dwork1[i] = (r0/potential->rad[i])/dwork1[i];
  }
  IntegrateF(potential, dwork1, orb1, orb2, type, dwork2, -1);
  for (i = 0; i < potential->maxrp; i++) {
    yk[i] += dwork2[i]/dwork1[i];
  }
      
  return 0;
}
      
int GetYk(const cfac_t *cfac, int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, 
	  int k1, int k2, RadIntType type) {
  int i, i0, i1, n;
  double a, b, a2, b2, max, max1;
  int index[3];
  SLATER_YK *syk;
  POTENTIAL *potential = cfac->potential;
  double dwork[MAXRP];

  if (k1 <= k2) {
    index[0] = k1;
    index[1] = k2;
  } else {
    index[0] = k2;
    index[1] = k1;
  }
  index[2] = k;

  syk = (SLATER_YK *) MultiSet(cfac->yk_array, index, NULL);
  if (syk->npts < 0) {
    GetYk1(potential, k, yk, dwork, orb1, orb2, type);
    max = 0;
    for (i = 0; i < potential->maxrp; i++) {
      dwork[i] *= yk[i];
      a = fabs(dwork[i]); 
      if (a > max) max = a;
    }
    max1 = max*EPS5;
    max = max*EPS4;
    a = dwork[potential->maxrp-1];
    for (i = potential->maxrp-2; i >= 0; i--) {
      if (fabs(dwork[i] - a) > max1) {
	break;
      }
    }
    i1 = i;
    for (i = i1; i >= 0; i--) {      
      b = fabs(a - dwork[i]);
      dwork[i] = log(b);
      if (b > max) {
	break;
      }
    }
    i0 = i;
    if (i0 == i1) {
      i0--;
      b = fabs(a - dwork[i0]);
      dwork[i0] = log(b);
    }
    syk->coeff[0] = a;    
    syk->npts = i0+1;
    syk->yk = malloc(sizeof(float)*(syk->npts));
    for (i = 0; i < syk->npts ; i++) {
      syk->yk[i] = yk[i];
    }
    n = i1 - i0 + 1;
    a = 0.0;
    b = 0.0;
    a2 = 0.0;
    b2 = 0.0;
    for (i = i0; i <= i1; i++) {      
      max = (potential->rad[i]-potential->rad[i0]);
      a += max;
      b += dwork[i];
      a2 += max*max;
      b2 += dwork[i]*max;
    }
    syk->coeff[1] = (a*b - n*b2)/(a*a - n*a2);       
    if (syk->coeff[1] >= 0) {
      i1 = i0 + (i1-i0)*0.3;
      if (i1 == i0) i1 = i0 + 1;
      for (i = i0; i <= i1; i++) {      
	max = (potential->rad[i]-potential->rad[i0]);
	a += max;
	b += dwork[i];
	a2 += max*max;
	b2 += dwork[i]*max;
      }
      syk->coeff[1] = (a*b - n*b2)/(a*a - n*a2);  
    }
    if (syk->coeff[1] >= 0) {
      syk->coeff[1] = -10.0/(potential->rad[i1]-potential->rad[i0]);
    } 
  } else {
    for (i = syk->npts-1; i < potential->maxrp; i++) {
      dwork[i] = pow(potential->rad[i], k);
    }
    for (i = 0; i < syk->npts; i++) {
      yk[i] = syk->yk[i];
    }
    i0 = syk->npts-1;
    a = syk->yk[i0]*dwork[i0];
    for (i = syk->npts; i < potential->maxrp; i++) {
      b = potential->rad[i] - potential->rad[i0];
      b = syk->coeff[1]*b;
      if (b < -20) {
	yk[i] = syk->coeff[0];
      } else {
	yk[i] = (a - syk->coeff[0])*exp(b);
	yk[i] += syk->coeff[0];
      }
      yk[i] /= dwork[i];
    }    
  }
      
  return 0;
}  

/*
 * integrate a function given by f with two cfac->orbitals.
 * type indicates the type of integral
 * id indicates whether integrate inward (-1) or outward (0)
 * if last_only is set, only the end point is returned in x,
 * otherwise, the whole function is returned (in x[])
 */
static int _Integrate(POTENTIAL *potential,
    const double *f, const ORBITAL *orb1, const ORBITAL *orb2,
	              RadIntType type, double *x, int id, int last_only)
{
  int i1, i2, ilast;
  int mode;
  double *r, _dwork[MAXRP], ext;

  if (last_only) {
    r = _dwork;
  } else {
    r = x;
  }
  
  /* zero-fill r */
  memset(r, 0, potential->maxrp*sizeof(double));

  /* first, the overlapping region */
  ilast = Min(orb1->ilast, orb2->ilast);
  if (id >= 0) {
    mode = 0;
  } else {
    mode = -1;
  }
  IntegrateSubRegion(potential, 0, ilast, f, orb1, orb2, type, r, mode, last_only);
  
  if (orb1->ilast == ilast && orb1->n == 0) {
    /* orb2 extends beyond orb1 AND orb1 belongs to a free electron */
    i1 = ilast + 1;
    i2 = orb2->ilast;
    if (i2 > i1) {
      if (type == INT_P1Q2) {
	type = INT_Q1P2;
	mode = 1;
        IntegrateSubRegion(potential, i1, i2, f, orb2, orb1, type, r, mode, last_only);
      } else {
	mode = 2;
	IntegrateSubRegion(potential, i1, i2, f, orb1, orb2, type, r, mode, last_only);
      }
      i2--;
    }
    
    /* if orb2 is a free orbital, continue till the last point of potential */
    if (orb2->n == 0) {
      i1 = orb2->ilast + 1;
      i2 = potential->maxrp - 1;
      mode = 3;
      IntegrateSubRegion(potential, i1, i2, f, orb1, orb2, type, r, mode, last_only);
      i2--;
    }
  } else
  if (orb2->ilast == ilast && orb2->n == 0) {
    /* orb1 extends beyond orb2 AND orb2 belongs to a free electron */
    i1 = ilast + 1;
    i2 = orb1->ilast;
    if (i2 > i1) {
      /* TODO: why no INT_P1Q2/INT_Q1P2 here??? */
      mode = 1;
      IntegrateSubRegion(potential, i1, i2, f, orb1, orb2, type, r, mode, last_only);
      i2--;
    }
    /* if orb1 is a free orbital, continue till the last point of potential */
    if (orb1->n == 0) {
      i1 = orb1->ilast + 1;
      i2 = potential->maxrp - 1;
      mode = 3;
      IntegrateSubRegion(potential, i1, i2, f, orb1, orb2, type, r, mode, last_only);
      i2--;
    }
  } else {
    i2 = ilast;
  }
  
  ext = r[i2];
  
  /* finally, extend till the last point of the potential if t < 0
     (and if not there already) */
  if (last_only) {
    *x = ext;
  } else {
    int i;
    
    if (id >= 0) {
      /* constant, no longer varying */
      for (i = i2+1; i < potential->maxrp; i++) {
	r[i] = ext;
      }
    } else {
      if (i2 > ilast) {
	/* the outermost part is zero (inward integration!) */
        /* NB: it should be zero-filled already, but some of IntegrateSubRegion
           above corrupts at least r[i2+1]; type=-2, mode=1 certainly */
        for (i = i2+1; i < potential->maxrp; i++) {
          r[i] = 0.0;
	}
	
        for (i = ilast + 1; i <= i2; i++) {
	  r[i] = ext - r[i];
	}
	for (i = 0; i <= ilast; i++) {
	  /* in the overlap region the inward integration is taken care
             by IntegrateSubRegion(potential, ), hence, + r[i] contrary to above */
          r[i] = r[ilast+1] + r[i];
	}
      }
    }
  }
  
  return 0;
}

/* Integrate and return a single (last) point */
int IntegrateS(POTENTIAL *potential, const double *f, const ORBITAL *orb1, const ORBITAL *orb2, 
	      RadIntType type, double *r, int id)
{
    return _Integrate(potential, f, orb1, orb2, type, r, id, 1);
}

/* Integrate and return the whole anti-derivative array */
int IntegrateF(POTENTIAL *potential, const double *f, const ORBITAL *orb1, const ORBITAL *orb2, 
	      RadIntType type, double x[], int id)
{
    return _Integrate(potential, f, orb1, orb2, type, x, id, 0);
}


static void AddEvenPoints(double *r, double *r1, int i0, int i1, int last_only) {
  int i;

  if (last_only) {
    r[i1-1] += r1[i1-1];
  } else {
    for (i = i0; i < i1; i++) {
      r[i] += r1[i];
    }
  }
}

/*
 * Radial integral (orb1|f|orb2) over subregion
 * between points with indices i0 and i1.
 * t - type (as in Integrate(), plus 7 - Q1*P2)
 * r - the integral array
 * mode of integration:
 * mode = -1: same as 0 but inward integration
 * mode =  0: overlap region
 * mode =  1:
 * mode =  2: same as 1 but orb1 and orb2 swapped AND ???
 * mode =  3:
 */
int IntegrateSubRegion(POTENTIAL *potential, int i0, int i1, 
		       const double *f,
                       const ORBITAL *orb1, const ORBITAL *orb2,
		       RadIntType type, double *r,
                       int mode, int last_only) {
  int i, j, ip, i2, id;
  double *large1, *large2, *small1, *small2;
  double *x, *y, *r1, *x1, *x2, *y1, *y2;
  double a, b, e1, e2, a2, r0 = 0.0;
  double _dwork1[MAXRP];
  double _dwork2[MAXRP];
  double _dwork3[MAXRP];
  double _dwork4[MAXRP];
  double _dwork5[MAXRP];
  double _dwork6[MAXRP];
  double _dwork7[MAXRP];
  double _phase[MAXRP];
  double _dphase[MAXRP];
  double _dphasep[MAXRP];

  if (i1 <= i0) return 0;

  x  = _dwork1;
  y  = _dwork2;
  x1 = _dwork3;
  x2 = _dwork4;
  y1 = _dwork5;
  y2 = _dwork6;
  r1 = _dwork7;
  i2 = i1;

  if (mode == -1) {
      mode = 0;
      id = -1;
  } else {
      id = 0;
  }

  switch (mode) {
  case 0:  /* mode =  0 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    switch (type) {
    case INT_P1P2pQ1Q2:
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] += small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] += (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] += (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] += (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case INT_P1P2:
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = large2[i]*a*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = large2[i]*a*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case INT_Q1Q2:
      for (i = i0; i <= i1; i++) {
	x[i] = small1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = (small1[i]*b + small1[ip]*a)*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = (small1[i]*b+small1[ip]*a)*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case INT_P1Q2pQ1P2:
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] += small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] += (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] += (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] += large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] += large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case INT_P1Q2mQ1P2:
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] -= small1[i] * large2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      } 
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] -= (small1[i]*b + small1[ip]*a)*large2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] -= (small1[i]*b+small1[ip]*a)*large2[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] -= large2[i]*a*small1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] -= large2[i]*a*(small1[i]*e2+small1[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    case INT_P1Q2:
      for (i = i0; i <= i1; i++) {
	x[i] = large1[i] * small2[i];
	x[i] *= f[i]*potential->dr_drho[i];
      }
      if (i1 == orb1->ilast && orb1->n == 0 && i < potential->maxrp) {
	if (i <= orb2->ilast) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  x[i] = large1[i]*a*small2[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb2->n == 0) {
	  ip = i+1;
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  e1 = sin(large2[ip]);
	  e2 = cos(large2[ip]);
	  x[i] = large1[i]*a*(small2[i]*e2+small2[ip]*e1);
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      } else if (i1 == orb2->ilast && orb2->n == 0 && i < potential->maxrp) {
	if (i <= orb1->ilast) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  x[i] = (small2[i]*b + small2[ip]*a)*large1[i];
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	} else if (orb1->n == 0) {
	  ip = i+1;
	  a = sin(large2[ip]);
	  b = cos(large2[ip]);
	  e1 = sin(large1[ip]);
	  e2 = cos(large1[ip]);
	  x[i] = (small2[i]*b+small2[ip]*a)*large1[i]*e1;
	  x[i] *= f[i]*potential->dr_drho[i];
	  i2 = i;
	}
      }
      break;
    default: /* error */
      return -1;
    }      
    NewtonCotes(r+i0, x+i0, i2-i0, last_only, id);
    break;

  case 1: /* mode = 1 */
    if (type == INT_P1Q2) { /* type INT_P1Q2 needs special treatments */
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    } else if (type == INT_Q1P2) {
      large1 = Large(orb1);
      large2 = Large(orb2);
      small1 = Small(orb1);
      small2 = Small(orb2);
      e1 = orb1->energy;
      e2 = orb2->energy;
      a2 = 0.5*FINE_STRUCTURE_CONST2;
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * large2[i];
	x[j] *= f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;	
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] = b * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, NULL, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    }
  case 2: /* mode = 2 */
    /* type INT_P1Q2/INT_Q1P2 is treated in mode = 1 */
    if (mode == 2) {
      const ORBITAL *tmp = orb1;
      orb1 = orb2;
      orb2 = tmp;
    }
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    switch (type) {
    case INT_P1P2pQ1Q2:
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] += small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case INT_P1P2:
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * large2[i];
	x[j] *= f[i]; 
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  a = large1[i]*a;
	  x[j] = a * large2[i];
	  x[j] *= f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, NULL, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case INT_Q1Q2:
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = small1[i] * small2[ip];
	x[j] *= f[i];
	y[j] = small1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  x[j] += b * small2[ip];
	  x[j] *= f[i];
	  y[j] = b * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;  
    case INT_P1Q2pQ1P2:
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] += small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) {
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] += b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    case INT_P1Q2mQ1P2:
      j = 0;
      for (i = i0; i <= i1; i+= 2) {
	ip = i+1;
	x[j] = large1[i] * small2[ip];
	x[j] -= small1[i] * large2[i];
	x[j] *= f[i];
	y[j] = large1[i] * small2[i] * f[i];
	_phase[j] = large2[ip];
	_dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	j++;
      }
      if (i < potential->maxrp) {
	ip = i+1;
	if (i > orb1->ilast && orb1->n == 0) { 
	  a = sin(large1[ip]);
	  b = cos(large1[ip]);
	  b = small1[i]*b + small1[ip]*a;
	  a = large1[i]*a;
	  x[j] = a * small2[ip];
	  x[j] -= b * large2[i];
	  x[j] *= f[i];
	  y[j] = a * small2[i] * f[i];
	  _phase[j] = large2[ip];
	  _dphase[j] = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	    /(large2[i]*large2[i]);
	  j++;
	  i2 = i;
	}
      }
      if (mode == 2) {
	r0 = r[i0];
	r[i0] = 0.0;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      if (mode == 2) {
	if (last_only) {
	  if (IsOdd(i2)) r[i2-1] = r0 - r[i2-1];
	  else r[i2] = r0 - r[i2];
	} else {
	  for (i = i0; i <= i2; i++) {
	    r[i] = r0 - r[i];
	  }
	}
      }
      if (IsOdd(i2)) r[i2] = r[i2-1];
      break;
    default:
      return -1;
    }
    break;

  case 3: /* mode = 3 */
    large1 = Large(orb1);
    large2 = Large(orb2);
    small1 = Small(orb1);
    small2 = Small(orb2);
    e1 = orb1->energy;
    e2 = orb2->energy;
    a2 = 0.5*FINE_STRUCTURE_CONST2;
    r1[i0] = 0.0;
    switch (type) {
    case INT_P1P2pQ1Q2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;	
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];	
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y2[j] += large1[i] * large2[i];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, NULL, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;	
    case INT_P1P2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y2[j] = -large1[i] * large2[i];
	y2[j] *= 0.5*f[i];
	y[j] = y2[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      } 
      IntegrateSinCos(potential, j, NULL, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, NULL, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;
    case INT_Q1Q2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = small1[i] * small2[ip];
	x2[j] = small1[ip] * small2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y1[j] = small1[i] * small2[i];
	y2[j] = small1[ip] * small2[ip];
	y[j] = y1[j] - y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = -x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = y1[j] + y2[j];
	y[j] *= 0.5*f[i];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;
    case INT_P1Q2pQ1P2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] -= small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];	
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;
    case INT_P1Q2mQ1P2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x2[j] = small1[i] * large2[i];
	x[j] = x1[j] - x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -large1[i] * small2[ip];
	y[j] += small1[ip] * large2[i];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j] + x2[j];
	x[j] *= 0.5*f[i];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;
    case INT_P1Q2:
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x1[j] = large1[i] * small2[i];
	x1[j] *= 0.5*f[i];
	x[j] = x1[j];
	y[j] = -large1[i] * small2[ip];
	y[j] *= 0.5*f[i];
	y2[j] = y[j];
	_phase[j] = large1[ip] + large2[ip];
	a = (1.0 + a2*(e1-potential->U[i]-potential->Vc[i]))
	  /(large1[i]*large1[i]);
	b = (1.0 + a2*(e2-potential->U[i]-potential->Vc[i]))
	  /(large2[i]*large2[i]);
	_dphase[j] = a + b;
	_dphasep[j] = a - b;
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphase, i0, r, last_only);
      j = 0;
      for (i = i0; i <= i1; i += 2) {
	ip = i+1;
	x[j] = x1[j];
	y[j] = -y2[j];
	_phase[j] = large1[ip] - large2[ip];
	j++;
      }
      IntegrateSinCos(potential, j, x, y, _phase, _dphasep, i0, r1, last_only);
      AddEvenPoints(r, r1, i0, i1, last_only);
      break;
    default:
      return -1;
    }
  default:
    return -1; 
  }

  return 0;
}

int IntegrateSinCos(POTENTIAL *potential, int j, double *x, double *y, 
		    double *phase, double *dphase, 
		    int i0, double *r, int last_only) {
  int i, k, m, n, q, s, i1;
  double si0, si1 = 0.0, cs0, cs1 = 0.0;
  double is0 = 0.0, is1 = 0.0, is2 = 0.0, is3 = 0.0;
  double ic0 = 0.0, ic1 = 0.0, ic2 = 0.0, ic3 = 0.0;
  double his0 = 0.0, his1 = 0.0, his2 = 0.0, his3 = 0.0;
  double hic0 = 0.0, hic1 = 0.0, hic2 = 0.0, hic3 = 0.0;
  double a0, a1, a2, a3;
  double d, p, dr;
  double *z, *u, *w, *u1, *w1;
  double xa1[MAXRP], xa2[MAXRP], xa3[MAXRP];
  double ya1[MAXRP], ya2[MAXRP], ya3[MAXRP];
  double _dwork1[MAXRP], _dwork2[MAXRP], _dwork3[MAXRP];

  w = _dwork1;
  z = _dwork2;
  u = _dwork3;
  
  if (phase[j-1] < 0) {
    for (i = 0; i < j; i++) {
      phase[i] = -phase[i];
      dphase[i] = -dphase[i];
      if (x) x[i] = -x[i];
      /* no need to invert y since cosine is an even function */
    }
  }
  
  for (i = 1, k = i0+2; i < j; i++, k += 2) {
    double h = phase[i] - phase[i-1];
    z[k] = 0.0;
    if (x != NULL) z[k] += x[i]*sin(phase[i]);
    if (y != NULL) z[k] += y[i]*cos(phase[i]);
    z[k] *= potential->dr_drho[k];
    if (i < 4) {
      if (h > 0.8) break;
    } else {
      if (h > 0.4) break;
    }
  }
  
  if (i > 1) {
    z[i0] = 0.0;
    if (x != NULL) z[i0] += x[0]*sin(phase[0]);
    if (y != NULL) z[i0] += y[0]*cos(phase[0]);
    z[i0] *= potential->dr_drho[i0];
    if (i == j) {
      i1 = i;
    } else {
      i1 = i + 1;
    }
    u1 = u + i1;
    w1 = w + i1;
    for (m = 0, n = i0; m < i; m++, n += 2) {
      u[m] = potential->rad[n];
      w[m] = potential->rad[n+1];
      z[n+1] = 0.0;
    }
    if (i1 > i) {
      u[i] = potential->rad[n];
    }
    UVIP3P(i1, u, phase, i, w, w1);
    if (x) {
      UVIP3P(i1, u, x, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*sin(w1[m]);
      }
    }
    if (y) {
      UVIP3P(i1, u, y, i, w, u1);
      for (m = 0, n = i0+1; m < i; m++, n += 2) {
        z[n] += u1[m]*cos(w1[m]);
      }
    }
    for (m = 0, n = i0+1; m < i; m++, n += 2) {
      z[n] *= potential->dr_drho[n];
    }
    NewtonCotes(r+i0, z+i0, k-2-i0, last_only, 0);
  }

  q = i-1;
  m = j-q;
  if (m < 2) {
    return 0;
  }

  for (n = 1; n <= 5; n++) {
    i1 = q-n;
    if (i1 < 0 || phase[i1] >= phase[i1+1]) {
      i1++;
      break;
    }
  }
  if (!last_only) {
    for (n = i1, s = i0; n < j; n++, s += 2) {
      u[n] = potential->rad[s];
      z[n] = potential->rad[s+1];
    }
    UVIP3P(j-i1, u+i1, phase+i1, m, z+q, w+q);
  }

  /* interpolate (calculate spline coefficients of) x & y on the phase axis */
  if (x != NULL) {
    for (n = i1; n < j; n++) {
      x[n] /= dphase[n];
    }
    UVIP3C(j-i1, phase+i1, x+i1, xa1, xa2, xa3);
  }
  if (y != NULL) {
    for (n = i1; n < j; n++) {
      y[n] /= dphase[n];
    }
    UVIP3C(j-i1, phase+i1, y+i1, ya1, ya2, ya3);
  }

  /* analytic piecewise integration of cubic spline with sin or cos */
  si0 = sin(phase[i-1]);
  cs0 = cos(phase[i-1]);
  for (; i < j; i++, k += 2) {
    double h;
    if (!last_only) {
      dr = w[i-1] - phase[i-1];
      si1 = sin(w[i-1]);
      cs1 = cos(w[i-1]);
      his0 = -(cs1 - cs0);
      hic0 = si1 - si0;
      p = dr;
      his1 = -p * cs1 + hic0;
      hic1 = p * si1 - his0;
      p *= dr; 
      his2 = -p * cs1 + 2.0*hic1;
      hic2 = p * si1 - 2.0*his1;
      p *= dr;
      his3 = -p * cs1 + 3.0*hic2;
      hic3 = p * si1 - 3.0*his2;
      r[k-1] = r[k-2];
    }
    if (x != NULL || y != NULL) {
        d = phase[i] - phase[i-1];
        si1 = sin(phase[i]);
        cs1 = cos(phase[i]);
        is0 = -(cs1 - cs0);
        ic0 = si1 - si0;
        p = d;
        is1 = -p * cs1 + ic0;
        ic1 = p * si1 - is0;
        p *= d;
        is2 = -p * cs1 + 2.0*ic1;
        ic2 = p * si1 - 2.0*is1; 
        p *= d;
        is3 = -p * cs1 + 3.0*ic2;
        ic3 = p * si1 - 3.0*is2;
    }
    r[k] = r[k-2];
    if (x != NULL) {
      a0 =   x[i-1];
      a1 = xa1[i-1-i1]; 
      a2 = xa2[i-1-i1];
      a3 = xa3[i-1-i1];
      h = a0*is0 + a1*is1 + a2*is2 + a3*is3;
      r[k] += h;
      if (!last_only) {
        h = a0*his0 + a1*his1 + a2*his2 + a3*his3;
        r[k-1] += h;
      }
    }
    if (y != NULL) {
      a0 =   y[i-1];
      a1 = ya1[i-1-i1];
      a2 = ya2[i-1-i1];
      a3 = ya3[i-1-i1];
      h = a0*ic0 + a1*ic1 + a2*ic2 + a3*ic3;
      r[k] += h;
      if (!last_only) {
        h = a0*hic0 + a1*hic1 + a2*hic2 + a3*hic3;
        r[k-1] += h;
      }
    }
    si0 = si1;
    cs0 = cs1;
  }

  return 0;
}

int InitRadial(cfac_t *cfac) {
  return 0;
}

int ReinitRadial(cfac_t *cfac, int m) {
  if (m < 0) return 0;
  SetSlaterCut(cfac, -1, -1);
  ClearOrbitalTable(cfac, m);
  FreeSimpleArray(cfac->slater_array);
  FreeSimpleArray(cfac->breit_array);
  FreeSimpleArray(cfac->residual_array);
  FreeSimpleArray(cfac->qed1e_array);
  FreeSimpleArray(cfac->vinti_array);
  FreeMultipoleArray(cfac);
  FreeMomentsArray(cfac);
  FreeYkArray(cfac);
  if (m < 2) {
    FreeGOSArray(cfac);
    if (m == 0) {
      if (cfac->optimize_control.n_screen > 0) {
	free(cfac->optimize_control.screened_n);
	cfac->optimize_control.n_screen = 0;
      }
      cfac->potential->flag = 0;
      
      SetRadialGrid(cfac, DMAXRP, -1.0, -1.0, -1.0);
      cfac->potential->uehling[0] = 0.0;
      
      cfac->n_awgrid = 1;
      cfac->awgrid[0] = EPS3;
    }
  }
  return 0;
}
