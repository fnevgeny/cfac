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
#include <time.h>
#include <math.h>

#include "consts.h"
#include "cfacP.h"
#include "radial.h"
#include "angular.h"
#include "dbase.h"
#include "transition.h"

void SetTransitionMode(cfac_t *cfac, int m) {
  cfac->tr_opts.mode = m;
}

void SetTransitionGauge(cfac_t *cfac, int m) {
  cfac->tr_opts.gauge = m;
}

void SetTransitionMaxE(cfac_t *cfac, int m) {
  cfac->tr_opts.max_e = m;
}

void SetTransitionMaxM(cfac_t *cfac, int m) {
  cfac->tr_opts.max_m = m;
}

void SetTransitionOptions(cfac_t *cfac, int gauge, int mode, 
			  int max_e, int max_m) {
  cfac->tr_opts.gauge = gauge;
  cfac->tr_opts.mode = mode;
  cfac->tr_opts.max_e = max_e;
  cfac->tr_opts.max_m = max_m;
}

int GetTransitionGauge(const cfac_t *cfac) {
  return cfac->tr_opts.gauge;
}

int GetTransitionMode(const cfac_t *cfac) {
  return cfac->tr_opts.mode;
}

static int _TRMultipole(cfac_t *cfac, double *rme, double *energy,
		int m, int lower, int upper) {
  int m2;
  int p1, p2, j1, j2;
  LEVEL *lev1, *lev2;
  double s, r, aw;
  int nz, i;
  ANGULAR_ZMIX *ang;
  
  *rme = 0.0;

  lev1 = GetLevel(cfac, lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(cfac, upper);
  if (lev2 == NULL) return -1;
  
  if (GetNumElectrons(cfac, lower) != GetNumElectrons(cfac, upper)) {
    return -1;
  }
  
  *energy = lev2->energy - lev1->energy;
  if (*energy < 0.0) {
    return -1;
  }
    
  aw = FINE_STRUCTURE_CONST*fabs(*energy);

  DecodePJ(lev1->pj, &p1, &j1);
  DecodePJ(lev2->pj, &p2, &j2);
  if (j1 == 0 && j2 == 0) return 0;

  m2 = 2*abs(m);    
  if (!Triangle(j1, j2, m2)) return 0;
  if (m > 0 && IsEven(p1+p2+m)) return 0;
  if (m < 0 && IsOdd(p1+p2-m)) return 0;    

  s = 0.0;

  nz = AngularZMix(cfac, &ang, lower, upper, m2, m2);
  if (nz <= 0) {
    return -1;
  }

  for (i = 0; i < nz; i++) {
    if (ang[i].k != m2) continue;
    if (cfac->tr_opts.mode == M_NR && m != 1) {
      r = MultipoleRadialNR(cfac, m, ang[i].k0, ang[i].k1, 
			    cfac->tr_opts.gauge);
    } else {
      r = MultipoleRadialFR(cfac, aw, m, ang[i].k0, ang[i].k1,
			    cfac->tr_opts.gauge);
    }
    s += r * ang[i].coeff;
  }
  if (nz > 0) {
    free(ang);	
  }  

  *rme = s;
  
  return 0;
}

typedef struct {
  int m;
  int valid;
  double energy;
  double rme;
} TRANS_T;

typedef struct {
  unsigned int dim;
  TRANS_T *transitions;
} TRM_CACHE_T;

static TRM_CACHE_T *trm_cache = NULL;

static TRM_CACHE_T *TRMultipole_cache_new(unsigned int dim)
{
  TRM_CACHE_T *cache;
  
  cache = calloc(1, sizeof(TRM_CACHE_T));
  if (!cache) {
    return NULL;
  }
  
  cache->transitions = calloc(dim*dim, sizeof(TRANS_T));
  if (!cache->transitions) {
    free(cache);
    return NULL;
  }
  
  cache->dim = dim;
  
  return cache;
}

static void TRMultipole_cache_free(TRM_CACHE_T *cache)
{
  if (!cache) {
    return;
  }
  
  if (cache->transitions) {
    free(cache->transitions);
  }
  
  free(cache);
}

/* If energy is not NULL, it is assigned trans. energy; */
int TRMultipole(cfac_t *cfac, double *rme, double *energy,
		int m, int lower, int upper) {
  double dE = 0.0;
  
  TRANS_T *trans = NULL;
  
  int res;
  
  if (trm_cache &&
      lower >= 0 && lower < trm_cache->dim &&
      upper >= 0 && upper < trm_cache->dim) {
    trans = &trm_cache->transitions[trm_cache->dim*upper + lower];
    if (trans->m == m && trans->valid) {
      if (energy) {
        *energy = trans->energy;
      }
      *rme = trans->rme;
      
      return 0;
    }
  }
  
  res = _TRMultipole(cfac, rme, &dE, m, lower, upper);
  if (energy) {
    *energy = dE;
  }
  
  if (trans) {
    trans->m      = m;
    trans->energy = dE;
    trans->rme    = *rme;
    trans->valid  = 1;
  }
  
  return res;
}  

int TRMultipoleEB(cfac_t *cfac, double *strength, double *energy, int m, int lower, int upper) {
  LEVEL *lev1, *lev2;
  int i1, m2, q;

  lev1 = GetEBLevel(cfac, lower);
  if (lev1 == NULL) return -1;
  lev2 = GetEBLevel(cfac, upper);
  if (lev2 == NULL) return -1;
  
  *energy = lev2->energy - lev1->energy;
  if (*energy <= 0.0) return -1;
  
  m2 = 2*abs(m);
  
  for (q = 0; q <= m2; q++) strength[q] = 0.0;

  for (i1 = 0; i1 < lev1->n_basis; i1++) {
    LEVEL *plev1;
    int ilev1, mlev1, j1, p1;
    int i2;

    if (lev1->mixing[i1] == 0) continue;

    DecodeBasisEB(lev1->basis[i1], &ilev1, &mlev1);
    plev1 = GetLevel(cfac, ilev1);
    DecodePJ(plev1->pj, &p1, &j1);
    
    for (i2 = 0; i2 < lev2->n_basis; i2++) {
      LEVEL *plev2;
      int ilev2, mlev2, j2, p2;
      double r, a, c;

      if (lev2->mixing[i2] == 0) continue;

      c = lev1->mixing[i1]*lev2->mixing[i2];
      DecodeBasisEB(lev2->basis[i2], &ilev2, &mlev2);
      plev2 = GetLevel(cfac, ilev2);
      DecodePJ(plev2->pj, &p2, &j2);
      
      if (TRMultipole(cfac, &r, NULL, m, ilev1, ilev2) != 0) {
        continue;
      }
      
      a = W3j(j1, m2, j2, -mlev1, mlev1-mlev2, mlev2);
      if (IsOdd((j1-mlev1)/2)) a = -a;
     
      q = (mlev1-mlev2)/2+abs(m);
      strength[q] += c*r*a;
      /*
      printf("%d %d %d %d %2d %2d %2d %10.3E %10.3E %10.3E %10.3E %10.3E\n",
	     lower, upper, ilev1, ilev2, mlev1, mlev2, q-1,c, a, r, c*a*r, 
	     strength[q]);
      */
    }
  }

  return 0;
}

int SaveTransitionEB0(cfac_t *cfac, int nlow, int *low, int nup, int *up, 
		      char *fn, int m) {
  int k, i, j, nq;
  double emin, emax, e0, s[101], et;
  F_HEADER fhdr;
  TRF_HEADER tr_hdr;
  TRF_RECORD r;
  LEVEL *lev1, *lev2;
  FILE *f;
  
  if (nlow <= 0 || nup <= 0) return -1;
  if (m == 1 || cfac->tr_opts.mode == M_FR) {
    k = 0;
    emin = 1E10;
    emax = 1E-10;
    for (i = 0; i < nlow; i++) {
      lev1 = GetEBLevel(cfac, low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetEBLevel(cfac, up[j]);
	e0 = lev2->energy - lev1->energy;
	if (e0 > 0) k++;
	if (e0 < emin && e0 > 0) emin = e0;
	if (e0 > emax) emax = e0;
      }
    }
      
    if (k == 0) {
      return 0;
    }
    
    emin *= FINE_STRUCTURE_CONST;
    emax *= FINE_STRUCTURE_CONST;
    e0 = 2.0*(emax-emin)/(emin+emax);
    
    FreeMultipoleArray(cfac);
    if (e0 < EPS3) {
      SetAWGrid(cfac, 1, emin, emax);
    } else if (e0 < 1.0) {
      SetAWGrid(cfac, 2, emin, emax);
    } else {
      SetAWGrid(cfac, 3, emin, emax);
    }
  }
  fhdr.type = DB_TRF;
  strcpy(fhdr.symbol, cfac_get_atomic_symbol(cfac));
  fhdr.atom = cfac_get_atomic_number(cfac);  
  lev1 = GetEBLevel(cfac, low[0]);
  DecodeBasisEB(lev1->pb, &i, &j);  
  tr_hdr.nele = GetNumElectrons(cfac, i);
  tr_hdr.multipole = m;
  tr_hdr.gauge = GetTransitionGauge(cfac);
  if (m == 1) { /* always FR for M1 transitions */
    tr_hdr.mode = M_FR;
  } else {
    tr_hdr.mode = GetTransitionMode(cfac);
  }
  nq = 2*abs(m) + 1;
  r.strength = (float *) malloc(sizeof(float)*nq);
  GetFields(cfac, &tr_hdr.bfield, &tr_hdr.efield, &tr_hdr.fangle);
    
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);

  for (j = 0; j < nup; j++) {
    for (i = 0; i < nlow; i++) {
      k = TRMultipoleEB(cfac, s, &et, m, low[i], up[j]);
      if (k != 0) continue;
      e0 = 0.0;
      for (k = 0; k < nq; k++) {
	r.strength[k] = s[k];
	if (s[k]) e0 = s[k];
      }
      if (e0 == 0.0) continue;
      r.lower = low[i];
      r.upper = up[j];
      WriteTRFRecord(f, &r);
    }
  }
  
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);
  free(r.strength);

  return 0;
}

double OscillatorStrength(int m, double e, double s, double *ga) {
  double aw, x;

  aw = FINE_STRUCTURE_CONST*e;
  if (m != 0) {
    int m2 = 2*abs(m);
    x = s*s*e*pow(aw, m2 - 2)/(m2 + 1);
  } else {
    x = s;
  }
  if (ga) {
    *ga = x*2.0*pow(aw,2)*FINE_STRUCTURE_CONST;
  }  

  return x;
}  

static int CompareInt(const void *a1, const void *a2) {
  int *i1, *i2;
  
  i1 = (int *) a1;
  i2 = (int *) a2;
  return (*i1 - *i2);
}

int OverlapLowUp(int nlow, int *low, int nup, int *up) {
  int i, j, n;
  int *lowinup, *upinlow, *icom;

  lowinup = (int *) malloc(sizeof(int)*nlow);
  upinlow = (int *) malloc(sizeof(int)*nup);

  for (i = 0; i < nlow; i++) {
    lowinup[i] = -1;
  }
  for (i = 0; i < nup; i++) {
    upinlow[i] = -1;
  }
  qsort(low, nlow, sizeof(int), CompareInt);
  if (up != low) {
    qsort(up, nup, sizeof(int), CompareInt);
  }
  for (i = 0; i < nlow; i++) {
    lowinup[i] = IBisect(low[i], nup, up);    
    if (lowinup[i] >= 0) {
      upinlow[lowinup[i]] = i;
    }
  }
  icom = (int *) malloc(sizeof(int)*nlow);
  n = 0;
  for (i = 0; i < nlow; i++) {
    if (lowinup[i] >= 0) icom[n++] = low[i];
  }
  j = 0;
  for (i = 0; i < nlow; i++) {
    if (lowinup[i] < 0) {
      low[j++] = low[i];
    }
  }
  for (i = 0; i < n; i++) {
    low[j++] = icom[i];
  }
  j = 0;
  for (i = 0; i < nup; i++) {
    if (upinlow[i] < 0) {
      up[j++] = up[i];
    }
  }
  for (i = 0; i < n; i++) {
    up[j++] = icom[i];
  }
  
  free(lowinup);
  free(upinlow);
  free(icom);

  return n;
}

/* 
 * Find complements and intersection of unsorted sets I and F:
 * K = I\F; L = F\I; M = I*F.
 * Return bool flag whether an allocation occured (and hence resultant
 * sets need to be free'd.)
 */
int cfac_overlap_if(unsigned ni, unsigned *i, unsigned nf, unsigned *f,
    unsigned *nk, unsigned **k, unsigned *nl, unsigned **l,
    unsigned *nm, unsigned **m)
{
    int j;
    int i_min, i_max, f_min, f_max, sn_min, sn_max;
    unsigned n_k = 0, n_l = 0, n_m = 0;
    char *s;
    int allocated = 0;
    
    *nm = *nk = *nl = 0;
    
    if (ni == 0 || nf == 0) {
        return allocated;
    }
    
    if (ni == nf && i == f) {
        /* degenerate case I = F */
        
        *nm = ni;
        *m  = i;
        
        return allocated;
    }
    
    i_min = i_max = i[0];
    f_min = f_max = f[0];
    
    for (j = 1; j < ni; j++) {
        if (i[j] < i_min) {
            i_min = i[j];
        }
        if (i[j] > i_max) {
            i_max = i[j];
        }
    }
    
    for (j = 1; j < nf; j++) {
        if (f[j] < f_min) {
            f_min = f[j];
        }
        if (f[j] > f_max) {
            f_max = f[j];
        }
    }
    
    if (i_max < f_min || i_min > f_max) {
        /* no overlap */
        
        *nk = ni;
        *k  = i;
        *nl = nf;
        *l  = f;
        
        return allocated;
    }
    
    sn_min = (i_min < f_min) ? i_min:f_min;
    sn_max = (i_max > f_max) ? i_max:f_max;
    
    s = calloc(sn_max - sn_min + 1, 1);
    
    for (j = 0; j < ni; j++) {
        s[i[j] - sn_min] |= 1;
    }
    for (j = 0; j < nf; j++) {
        s[f[j] - sn_min] |= 2;
    }
    
    for (j = sn_min; j <= sn_max; j++) {
        switch (s[j - sn_min]) {
        case 1:
            n_k++;
            break;
        case 2:
            n_l++;
            break;
        case 3:
            n_m++;
            break;
        }
    }
    
    *nk = n_k;
    *nl = n_l;
    *nm = n_m;

    *k = malloc(sizeof(int)*(*nk));
    *l = malloc(sizeof(int)*(*nl));
    *m = malloc(sizeof(int)*(*nm));
    
    allocated = 1;

    n_k = n_m = n_l = 0;
    for (j = sn_min; j <= sn_max; j++) {
        switch (s[j - sn_min]) {
        case 1:
            (*k)[n_k++] = j;
            break;
        case 2:
            (*l)[n_l++] = j;
            break;
        case 3:
            (*m)[n_m++] = j;
            break;
        }
    }

    free(s);
    
    return allocated;
}



typedef struct {
  TR_RECORD r;
  TR_EXTRA rx;
  int ks[2];
} TR_DATUM;


static int CompareNRConfig(const void *p1, const void *p2) {
  CONFIG *c1, *c2;

  c1 = (CONFIG *) p1;
  c2 = (CONFIG *) p2;
  if (c1->nnrs > c2->nnrs) return 1;
  else if (c1->nnrs < c2->nnrs) return -1;
  else {
    return memcmp(c1->nrs, c2->nrs, sizeof(int)*c1->nnrs);
  }
}


static int CompareTRDatum(const void *p1, const void *p2) {
  TR_DATUM *r1, *r2;

  r1 = (TR_DATUM *) p1;
  r2 = (TR_DATUM *) p2;
  if ((r1->r).upper < (r2->r).upper) return -1;
  else if ((r1->r).upper > (r2->r).upper) return 1;
  else {
    if ((r1->r).lower < (r2->r).lower) return -1;
    else if ((r1->r).lower > (r2->r).lower) return 1;
    else return 0;
  }
}

static double FKB(cfac_t *cfac, int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(cfac, ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(cfac, kb)->kappa, &jb, &ib);

  if (!Triangle(ia, k, ia) || !Triangle(ib, k, ib)) return 0.0;
  a = W3j(ja, k, ja, 1, 0, -1)*W3j(jb, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0;
  Slater(cfac, &b, ka, kb, ka, kb, k/2, 0);
  
  b *= a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;

  return b;
}

static double GKB(cfac_t *cfac, int ka, int kb, int k) {
  int ja, jb, ia, ib;
  double a, b;
  
  GetJLFromKappa(GetOrbital(cfac, ka)->kappa, &ja, &ia);
  GetJLFromKappa(GetOrbital(cfac, kb)->kappa, &jb, &ib);
  
  if (IsOdd((ia+k+ib)/2) || !Triangle(ia, k, ib)) return 0.0;
  a = W3j(ja, k, jb, 1, 0, -1);
  if (fabs(a) < EPS30) return 0.0;
  Slater(cfac, &b, ka, kb, kb, ka, k/2, 0);
  
  b *= a*a*(ja+1.0)*(jb+1.0);
  if (IsEven((ja+jb)/2)) b = -b;
  return b;
}

static double ConfigEnergyVarianceParts0(cfac_t *cfac,
    SHELL *bra, int ia, int ib, int m2, int p) {
  int ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb;
  double a, b, c, d, e;

  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(cfac, bra[ia].n, bra[ia].kappa, 0);
  if (p > 0) {
    jb = GetJFromKappa(bra[ib].kappa);
    kb = OrbitalIndex(cfac, bra[ib].n, bra[ib].kappa, 0);
  }
  e = 0.0;
  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = -a + W6j(ja, ja, k, ja, ja, kp);
	if (k == kp) b += 1.0/(k+1.0);
	b *= a*FKB(cfac, ka, ka, k)*FKB(cfac, ka, ka, kp);
	e += b;
      }
    }
    break;
  case 1:
    k0 = 4;
    k1 = 2*ja;
    kp0 = 4;
    kp1 = 2*jb;
    kp1 = Min(kp1, k1);
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 4) {
	b = a - W6j(ja, ja, kp, ja, ja, k);
	if (k == kp) b -= 1.0/(k+1.0);
	b *= W6j(ja, ja, kp, jb, jb, m2);
	b /= 0.5*ja;
	b *= FKB(cfac, ka, ka, k)*FKB(cfac, ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb)/2+1)) e = -e;
    break;
  case 2:
    k0 = 4;
    k1 = 2*ja;
    kp0 = abs(ja-jb);
    kp1 = ja + jb;
    a = 1.0/(ja*(ja+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, ja);
	b = -b*b;
	b += W6j(jb, jb, k, ja, ja, m2)*W6j(jb, jb, k, ja, ja, kp);
	if (m2 == kp) {
	  b -= (1.0/(m2+1.0)-1.0/(jb+1.0))*a;
	} else {
	  b += a/(jb+1.0);
	}
	b /= 0.5*ja;
	b *= FKB(cfac, ka, ka, k)*GKB(cfac, ka, kb, kp);
	e += b;
      }
    }
    if (IsOdd((ja+jb)/2+1)) e = -e;
    break;
  case 3:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      for (kp = k0; kp <= k1; kp += 4) {
	b = 0.0;
	if (k == kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, ja, k, ja, m2, jb, kp, jb, jb);
	b -= W6j(ja, jb, m2, jb, ja, k)*W6j(ja, jb, m2, jb, ja, kp)/ja;
	b /= ja;
	b *= FKB(cfac, ka, kb, k)*FKB(cfac, ka, kb, kp);
	e += b;
      }
    }
    break;
  case 4:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0.0;
	if (k == kp) b += 1.0/((k+1.0)*(jb+1.0));
	b -= W9j(ja, jb, k, jb, m2, ja, kp, ja, jb);
	c = -1.0/(jb+1.0);
	d = c;
	if (k == m2) {
	  c += 1.0/(m2+1.0);
	}
	if (kp == m2) {
	  d += 1.0/(m2+1.0);
	}
	b -= c*d/ja;
	b /= ja;
	b *= GKB(cfac, ka, kb, k)*GKB(cfac, ka, kb, kp);
	e += b;
      }
    }
    break;
  case 5:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k1, k);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = -W6j(jb, jb, k, ja, ja, kp)/(jb+1.0);
	c = W6j(k, kp, m2, ja, jb, jb)*W6j(k, kp, m2, jb, ja, ja);
	if (IsOdd((ja+jb)/2)) b += c;
	else b -= c;
	c = -1.0/(jb+1.0);
	if (kp == m2) c += 1.0/(m2+1.0);
	c *= W6j(ja, jb, m2, jb, ja, k);
	b += c/ja;
	b /= 0.5*ja;
	b *= FKB(cfac, ka, kb, k)*GKB(cfac, ka, kb, kp);
	e += b;
      }
    }
    break;
  case 6:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      b = a/(k+1.0);
      c = FKB(cfac, ka, kb, k);
      b *= c*c;
      e += b;
    }
    break;
  case 7:
    k0 = abs(ja-jb);
    k1 = ja+jb;
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = k0; kp <= k1; kp += 2) {
	b = 0;
	if (k == kp) {
	  b += 1.0/(k+1.0);
	}
	b -= a;
	c = GKB(cfac, ka, kb, k);
	d = GKB(cfac, ka, kb, kp);
	e += a*b*c*d;
      }
    }
    break;
  case 8:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    kp0 = abs(ja-jb);
    kp1 = ja+jb;
    kp0 = k0;
    kp1 = k1;
    a = 1.0/((ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, ja, kp, ja, jb, k);
	if (fabs(b) < EPS30) continue;
	b *= 2.0*a;
	if (IsOdd(kp/2)) b = -b;
	c = FKB(cfac, ka, kb, k);
	d = GKB(cfac, ka, kb, kp);
	e += b*c*d;
      }
    }
    break;
  }

  return e;
}

static double ConfigEnergyVarianceParts1(cfac_t *cfac, SHELL *bra, int i,
					 int ia, int ib, int m2, int p) {
  int js, ja, jb, k, kp, k0, k1, kp0, kp1, ka, kb, ks;
  double a, b, e;
  
  js = GetJFromKappa(bra[i].kappa);
  ks = OrbitalIndex(cfac, bra[i].n, bra[i].kappa, 0);
  ja = GetJFromKappa(bra[ia].kappa);
  ka = OrbitalIndex(cfac, bra[ia].n, bra[ia].kappa, 0);
  jb = GetJFromKappa(bra[ib].kappa);
  kb = OrbitalIndex(cfac, bra[ib].n, bra[ib].kappa, 0);
  e = 0.0;

  switch (p) {
  case 0:
    k0 = 4;
    k1 = 2*ja;
    k = 2*jb;
    k1 = Min(k, k1);
    k = 2*js;
    k1 = Min(k, k1);
    for (k = k0; k <= k1; k += 4) {
      b = W6j(ja, ja, k, jb, jb, m2);
      if (fabs(b) < EPS30) continue;
      b *= 2.0/((k+1.0)*(js+1.0));
      b *= FKB(cfac, ks, ka, k)*FKB(cfac, ks, kb, k);
      if (IsOdd((ja+jb)/2)) b = -b;
      e += b;
    }
    break;
  case 1:
    k0 = abs(js-ja);
    k1 = js+ja;
    kp0 = abs(js-jb);
    kp1 = js+jb;
    a = 1.0/((js+1.0)*(ja+1.0)*(jb+1.0));
    for (k = k0; k <= k1; k += 2) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(k, kp, m2, jb, ja, js);
	b = -b*b + a;
	b /= (js+1.0);
	b *= 2.0*GKB(cfac, ks, ka, k)*GKB(cfac, ks, kb, kp);
	if (IsOdd((ja+jb)/2+1)) b = -b;
	e += b;
      }
    }
    break;
  case 2:
    k0 = 4;
    k1 = 2*js;
    k = 2*ja;
    k1 = Min(k, k1);
    kp0 = abs(js-jb);
    kp1 = js+jb;
    for (k = k0; k <= k1; k += 4) {
      for (kp = kp0; kp <= kp1; kp += 2) {
	b = W6j(jb, jb, k, ja, ja, m2);
	b *= W6j(jb, jb, k, js, js, kp);
	if (fabs(b) < EPS30) continue;
	b /= (js+1.0);
	if (IsOdd((ja+jb+kp)/2)) b = -b;
	b *= 2.0*FKB(cfac, ks, ka, k)*GKB(cfac, ks, kb, kp);
	e += b;
      }
    }
    break;
  }

  return e;
}

double ConfigEnergyVariance(cfac_t *cfac,
    int ns, SHELL *bra, int ia, int ib, int m2) {
  int i, js, p;
  double e, a, b, c;
  
  e = 0.0;
  for (i = 0; i < ns; i++) {
    js = GetJFromKappa(bra[i].kappa);
    a = bra[i].nq;
    b = js+1.0 - bra[i].nq;
    if (i == ia) {
      a -= 1.0;
    }
    if (i == ib) {
      b -= 1.0;
    }
    if (a == 0.0 || b == 0.0) continue;
    a = a*b;
    b = 0.0;
    if (i == ia) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(cfac, bra, ia, ib, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else if (i == ib) {
      for (p = 0; p < 6; p++) {
	c = ConfigEnergyVarianceParts0(cfac, bra, ib, ia, m2, p);
	b += c;
      }
      b /= js-1.0;
    } else {
      for (p = 6; p < 9; p++) {
	c = ConfigEnergyVarianceParts0(cfac, bra, i, ia, m2, p);
	b += c;
	c = ConfigEnergyVarianceParts0(cfac, bra, i, ib, m2, p);
	b += c;
      }
      c = ConfigEnergyVarianceParts1(cfac, bra, i, ia, ib, m2, 0);
      b += c;
      c = ConfigEnergyVarianceParts1(cfac, bra, i, ia, ib, m2, 1);
      b += c;
      c = ConfigEnergyVarianceParts1(cfac, bra, i, ia, ib, m2, 2);
      b += c;
      c = ConfigEnergyVarianceParts1(cfac, bra, i, ib, ia, m2, 2);
      b += c;
      b /= js;
    }

    e += a*b;
  }
    
  if (e < 0.0) e = 0.0;
  return e;
}

double ConfigEnergyShift(cfac_t *cfac,
    int ns, SHELL *bra, int ia, int ib, int m2) {
  double qa, qb, a, b, c, e;
  int ja, jb, k, kmin, kmax;
  int k0, k1;

  qa = bra[ia].nq;
  qb = bra[ib].nq;
  ja = GetJFromKappa(bra[ia].kappa);
  jb = GetJFromKappa(bra[ib].kappa);
  if (qa == 1 && qb == 0) e = 0.0;
  else {
    e = (qa-1.0)/ja - qb/jb;
    if (e != 0.0) {
      kmin = 4;
      kmax = 2*ja;
      k = 2*jb;
      kmax = Min(kmax, k);
      a = 0.0;
      k0 = OrbitalIndex(cfac, bra[ia].n, bra[ia].kappa, 0);
      k1 = OrbitalIndex(cfac, bra[ib].n, bra[ib].kappa, 0);
      for (k = kmin; k <= kmax; k += 4) {
	b = W6j(k, ja, ja, m2, jb, jb);
	if (fabs(b) > EPS30) {
	  a -= b*FKB(cfac, k0, k1, k);
	}
      }
      kmin = abs(ja-jb);
      kmax = ja + jb;
      c = 1.0/((ja+1.0)*(jb+1.0));
      for (k = kmin; k <= kmax; k += 2) {
	if (k == m2) {
	  b = 1.0/(m2+1.0)-c;
	} else {
	  b = -c;
	}
	a += b*GKB(cfac, k0, k1, k);
      }
      if (IsEven((ja+jb)/2)) a = -a;
      e *= a;
    }
  }

  return e;
}

static double ConfigEnergyShiftCI(cfac_t *cfac, int nrs0, int nrs1) {
  int n0, k0, q0, n1, k1, q1;
  int j0, j1, s0, s1, kmax, k;
  double a, g, pk, w, sd, q;

  UnpackNRShell(&nrs0, &n0, &k0, &q0);
  UnpackNRShell(&nrs1, &n1, &k1, &q1);
  q = (q0-1.0)/(2.0*k0+1.0) - q1/(2.0*k1+1.0);
  if (q == 0.0) return q;
  g = 0.0;
  for (j0 = k0-1; j0 <= k0+1; j0 += 2) {
    if (j0 < 0) continue;
    s0 = OrbitalIndex(cfac, n0, GetKappaFromJL(j0, k0), 0);
    for (j1 = k1-1; j1 <= k1+1; j1 += 2) {
      if (j1 < 0) continue;
      s1 = OrbitalIndex(cfac, n1, GetKappaFromJL(j1, k1), 0);
      w = W6j(j0, k0, 1, k1, j1, 2);
      w = w*w*(j0+1.0)*(j1+1.0)*0.5;
      pk = ((j0+1.0)*(j1+1.0))/((k0+1.0)*(k1+1.0)*4.0) - w;
      Slater(cfac, &sd, s0, s1, s0, s1, 0, 0);
      g += sd*pk;
      kmax = 2*j0;
      k = 2*j1;
      kmax = Min(k, kmax);
      for (k = 4; k <= kmax; k += 4) {	
	pk = W3j(k0, k, k0, 0, 0, 0);
	if (fabs(pk) > 0) {
	  pk *= W3j(k1, k, k1, 0, 0, 0);
	  if (fabs(pk) > 0) {
	    pk *= 0.25;
	    a = W6j(j0, 2, j1, k1, 1, k0);
	    if (fabs(a) > 0) {
	      a = a*a;
	      a *= W6j(j0, k, j0, k0, 1, k0);
	      if (fabs(a) > 0) {
		a *= W6j(j1, k, j1, k1, 1, k1);
		if (fabs(a) > 0) {
		  a *= W6j(j0, k, j0, j1, 2, j1);
		  a *= 2.0*(j0+1.0)*(j1+1.0)*(k0+1.0)*(k1+1.0);
		}
	      }
	    }
	    a += W6j(k0, k, k0, k1, 2, k1);
	    pk *= a;
	  }
	  if (fabs(pk) > 0) {
	    Slater(cfac, &sd, s0, s1, s0, s1, k/2, 0);
	    g += pk*sd*(j0+1.0)*(j1+1.0);
	  }
	}
      }
      pk = W3j(k0, 2, k1, 0, 0, 0);
      if (fabs(pk) > 0) {
	pk = pk*pk;
	a = W6j(j0, 2, j1, k1, 1, k0);
	a *= a;
	pk *= a;
	pk *= (k0+1.0)*(k1+1.0)/3.0;
	pk *= (1.0 - 0.5*(j0+1.0)*(j1+1.0)*a);
	if (fabs(pk) > 0) {
	  Slater(cfac, &sd, s0, s1, s1, s0, 1, 0);
	  g += pk*sd*(j0+1.0)*(j1+1.0);
	}
      }
    }
  }

  return q*g;
}

static int TRMultipoleUTA(cfac_t *cfac, double *strength, TR_EXTRA *rx,
		   int m, int lower, int upper, int *ks) {
  int m2, ns, k0, k1, q1, q2;
  int p1, p2, j1, j2, ia, ib;
  LEVEL *lev1, *lev2;
  double r, aw, e0;
  INTERACT_DATUM *idatum;
  
  *ks = 0;
  
  *strength = 0.0;
  lev1 = GetLevel(cfac, lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(cfac, upper);
  if (lev2 == NULL) return -1;

  p1 = lev1->pj;
  p2 = lev2->pj;
  if (m > 0 && IsEven(p1+p2+m)) return -1;
  if (m < 0 && IsOdd(p1+p2+m)) return -1;
  
  idatum = NULL;
  ns = GetInteract(cfac, &idatum, NULL, NULL, lev1->iham, lev2->iham,
		   lev1->pb, lev2->pb, 0, 0, 0);
  if (ns <= 0) return -1;
  if (idatum->s[0].index < 0 || idatum->s[3].index >= 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  if (idatum->s[0].nq_bra > idatum->s[0].nq_ket) {
    ia = ns-1-idatum->s[0].index;
    ib = ns-1-idatum->s[1].index;
    j1 = idatum->s[0].j;
    j2 = idatum->s[1].j;
    q1 = idatum->s[0].nq_bra;
    q2 = idatum->s[1].nq_bra;
    k0 = OrbitalIndex(cfac, idatum->s[0].n, idatum->s[0].kappa, 0.0);
    k1 = OrbitalIndex(cfac, idatum->s[1].n, idatum->s[1].kappa, 0.0);
    PackNRShell(ks, idatum->s[0].n, idatum->s[0].kl, q1);
    if (idatum->s[0].kappa > 0) ks[0] |= 0x01000000;
    PackNRShell(ks+1, idatum->s[1].n, idatum->s[1].kl, q2);
    if (idatum->s[1].kappa > 0) ks[1] |= 0x01000000;
  } else {
    ia = ns-1-idatum->s[1].index;
    ib = ns-1-idatum->s[0].index;
    j1 = idatum->s[1].j;
    j2 = idatum->s[0].j;
    q1 = idatum->s[1].nq_bra;
    q2 = idatum->s[0].nq_bra;
    k0 = OrbitalIndex(cfac, idatum->s[1].n, idatum->s[1].kappa, 0.0);
    k1 = OrbitalIndex(cfac, idatum->s[0].n, idatum->s[0].kappa, 0.0);
    PackNRShell(ks, idatum->s[1].n, idatum->s[1].kl, q1);
    if (idatum->s[1].kappa > 0)  ks[0] |= 0x01000000;
    PackNRShell(ks+1, idatum->s[0].n, idatum->s[0].kl, q2);
    if (idatum->s[0].kappa > 0)  ks[1] |= 0x01000000;
  }

  m2 = 2*abs(m);
  if (!Triangle(j1, j2, m2)) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  e0 = lev2->energy - lev1->energy;
  rx->de = 0.0;
  if (m < 0) {
    rx->de += ConfigEnergyShift(cfac, ns, idatum->bra, ia, ib, m2);
    rx->sdev = sqrt(ConfigEnergyVariance(cfac, ns, idatum->bra, ia, ib, m2));
  } else {
    rx->sdev = 0.0;
  }
  aw = FINE_STRUCTURE_CONST*(e0 + rx->de);
  if (aw < 0.0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }
  
  if (cfac->tr_opts.mode == M_NR && m != 1) {
    r = MultipoleRadialNR(cfac, m, k0, k1, cfac->tr_opts.gauge);
  } else {
    r = MultipoleRadialFR(cfac, aw, m, k0, k1, cfac->tr_opts.gauge);
  }

  *strength = sqrt((lev1->ilev+1.0)*q1*(j2+1.0-q2)/((j1+1.0)*(j2+1.0)))*r;
  
  free(idatum->bra);
  free(idatum);
  return 0;
}


/* save radiative transitions; low & up are assumed NOT overlapping */
static int crac_save_rtrans0(cfac_t *cfac,
    unsigned nlow, const unsigned *low, unsigned nup, const unsigned *up,
    int mpole, int mode,
    cfac_tr_sink_t sink, void *udata) {
  unsigned i, j, ntr;
  double emin, emax;

  if (!nlow || !nup) return 0;

  emin = 0.0;
  emax = 0.0;
  ntr = 0;
  for (i = 0; i < nlow; i++) {
    LEVEL *llev = GetLevel(cfac, low[i]);
    if (!llev) {
      return -1;
    }
    for (j = 0; j < nup; j++) {
      double dE;
      LEVEL *ulev = GetLevel(cfac, up[j]);
      if (!ulev) {
        return -1;
      }
      dE = ulev->energy - llev->energy;
      
      if (dE < 0) {
        continue;
      }
      if (!emin || dE < emin) {
        emin = dE;
      }
      if (dE > emax) {
        emax = dE;
      }
      ntr++;
    }
  }
  
  if (!ntr) {
    return 0;
  }

  if (mode == M_FR && cfac->tr_opts.fr_interpolate) {
      double e0;
      emin *= FINE_STRUCTURE_CONST;
      emax *= FINE_STRUCTURE_CONST;
      e0 = 2.0*(emax-emin)/(emin+emax);

      FreeMultipoleArray(cfac);
      if (e0 < EPS3) {
        SetAWGrid(cfac, 1, emin, emax);
      } else if (e0 < 1.0) {
        SetAWGrid(cfac, 2, emin, emax);
      } else {
        SetAWGrid(cfac, 3, emin, emax);
      }
  }
  
  if (cfac->uta) {
    TR_DATUM *rd;
    double gf;
    double e0;
    int ic0, ic1, nic0, nic1, *nc0, *nc1, j0, j1, ntr;
    int imin, imax, jmin, jmax, nrs0, nrs1, ir, ir0;
    double ep, em, wp, wm, w0, de, cp, cm;
    CONFIG *c0, *c1;
    LEVEL *lev1, *lev2;
    int k;
    const int mj = 0xFF000000, mn = 0xFFFFFF;

    // qsort(low, nlow, sizeof(int), CompareNRLevel);
    nc0 = malloc(sizeof(int)*nlow);
    ic0 = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(cfac, low[i]);
      c1 = GetConfigFromGroup(cfac, lev1->iham, lev1->pb);
      if (i > 0 && CompareNRConfig(c1, c0)) {
	nc0[ic0++] = i;
      }
      c0 = c1;
    }
    nc0[ic0] = nlow;
    nic0 = ic0+1;
    
    
    if (up != low) {
      // qsort(up, nup, sizeof(int), CompareNRLevel);
      nc1 = malloc(sizeof(int)*nup);
      ic1 = 0;
      for (i = 0; i < nup; i++) {
	lev1 = GetLevel(cfac, up[i]);
	c1 = GetConfigFromGroup(cfac, lev1->iham, lev1->pb);
	if (i > 0 && CompareNRConfig(c1, c0) != 0) {
	  nc1[ic1++] = i;
	}
	c0 = c1;
      }
      nc1[ic1] = nup;
      nic1 = ic1+1;
    } else {
      nc1 = nc0;
      nic1 = nic0;
    }
    
    
    imin = 0;
    for (ic0 = 0; ic0 < nic0; ic0++) {
      imax = nc0[ic0];
      jmin = 0;
      lev1 = GetLevel(cfac, low[imin]);
      c0 = GetConfigFromGroup(cfac, lev1->iham, lev1->pb);
      for (ic1 = 0; ic1 < nic1; ic1++) {
	jmax = nc1[ic1];
	lev2 = GetLevel(cfac, up[jmin]);
	c1 = GetConfigFromGroup(cfac, lev2->iham, lev2->pb);
	ir = 0;
	ntr = (jmax-jmin)*(imax-imin);
	rd = malloc(sizeof(TR_DATUM)*ntr);
	ep = 0.0;
	em = 0.0;
	e0 = 0.0;
	wp = 0.0;
	wm = 0.0;
	w0 = 0.0;
	ir0 = -1;
	for (i = imin; i < imax; i++) {
	  for (j = jmin; j < jmax; j++) {
	    k = TRMultipoleUTA(cfac,
                &gf, &(rd[ir].rx), mpole, low[i], up[j], rd[ir].ks);
	    if (k != 0) {
	      rd[ir].r.lower = -1;
	      rd[ir].r.upper = -1;
	      ir++;
	      continue;
	    }
	    ir0 = ir;
	    rd[ir].r.lower = low[i];
	    rd[ir].r.upper = up[j];
	    rd[ir].r.rme = gf;
	    rd[ir].rx.sci = 1.0;
	    if (mpole == -1) {
              double energy = GetLevel(cfac, up[j])->energy -
                              GetLevel(cfac, low[i])->energy + rd[ir].rx.de;

              gf = OscillatorStrength(mpole, energy,
				      rd[ir].r.rme, NULL);
	      j0 = rd[ir].ks[0]&mj;
	      j1 = rd[ir].ks[1]&mj;
	      if (j0==0 && j1==0) {
		wp += gf;
		ep += gf*energy;
	      } else if (j0 && j1) {
		wm += gf;
		em += gf*energy;
	      }
	      e0 += gf*energy;
	      w0 += gf;
              
              printf("%g %g %g -- %g %g %g\n", e0, ep, em, w0, wp, wm);
	    }
	    ir++;
	  }
	}
	if (wm > 0.0 && ir0 >= 0) {
	  nrs0 = 0;
	  nrs1 = 0;
	  for (i = 0; i < c0->nnrs; i++) {
	    if (c0->nrs[i]>>8 == (rd[ir0].ks[0]&mn)>>8) nrs0 = c0->nrs[i];
	    else if (c0->nrs[i]>>8 == (rd[ir0].ks[1]&mn)>>8) nrs1 = c1->nrs[i];
	  }
	  if (nrs0 == 0) {
	    nrs0 = rd[ir0].ks[0] & mn;
	  }
	  if (nrs1 == 0) {
	    nrs1 = rd[ir0].ks[1] & mn;
	  }
	  ep /= wp;
	  em /= wm;
	  e0 /= w0;
	  de = ConfigEnergyShiftCI(cfac, nrs0, nrs1);
	  cm = 1.0 + de/(e0 - ep);
	  cp = 1.0 + de/(e0 - em);
	  if (cm < EPS3) {
	    cp = (wp+wm)/wp;
	    cm = 0.0;
	  } else if (cp < EPS3) {
	    cm = (wp+wm)/wm;
	    cp = 0.0;
	  }
	  ir = 0;
	  for (i = imin; i < imax; i++) {
	    for (j = jmin; j < jmax; j++) {
	      if (rd[ir].r.lower < 0) {
		ir++;
		continue;
	      }
	      j0 = rd[ir].ks[0]&mj;
	      j1 = rd[ir].ks[1]&mj;
	      if (j0==0 && j1==0) {
		rd[ir].rx.sci = cp;
	      } else if (j0 && j1) {
		rd[ir].rx.sci = cm;
	      }
	      ir++;
	    }
	  }
	}
	qsort(rd, ntr, sizeof(TR_DATUM), CompareTRDatum);
	ir = 0;
	for (ir = 0; ir < ntr; ir++) {
          cfac_rtrans_data_t rtdata;
	  if (rd[ir].r.lower < 0) {
	    continue;
	  }
	  // WriteTRRecord(f, &(rd[ir].r), &(rd[ir].rx));
          rtdata.fi = rd[ir].r.lower;
          rtdata.ii = rd[ir].r.upper;
          rtdata.rme = rd[ir].r.rme;
          
          rtdata.uta_de = rd[ir].rx.de;
          rtdata.uta_sd = rd[ir].rx.sdev;
          rtdata.uta_ci = rd[ir].rx.sci;
          
          if (sink(cfac, &rtdata, udata) != 0) {
            return -1;
          }
	}
	free(rd);
	jmin = jmax;
      }
      imin = imax;
    }
    free(nc0);
    if (up != low) free(nc1);
    
    return 0;
  }

  for (j = 0; j < nup; j++) {
    LEVEL *ulev = GetLevel(cfac, up[j]);
    for (i = 0; i < nlow; i++) {
      double rme;
      cfac_rtrans_data_t rtdata;

      if (mode == M_FR && !cfac->tr_opts.fr_interpolate) {
        LEVEL *llev = GetLevel(cfac, low[i]);
        double dE = ulev->energy - llev->energy;
        if (dE < 0) {
          continue;
        }
        
        FreeMultipoleArray(cfac);
        SetAWGrid(cfac, 1, dE*FINE_STRUCTURE_CONST, dE*FINE_STRUCTURE_CONST);
      }
      
      if (TRMultipole(cfac, &rme, NULL, mpole, low[i], up[j]) != 0) {
        continue;
      }
      if (fabs(rme) < EPS30) continue;
      
      rtdata.fi = low[i];
      rtdata.ii = up[j];
      rtdata.rme = rme;
      if (sink(cfac, &rtdata, udata) != 0) {
        return -1;
      }
    }
  }

  return 0;
}

int crac_calculate_rtrans(cfac_t *cfac,
    unsigned nlow, unsigned *low, unsigned nup, unsigned *up,
    int mpole, int mode,
    cfac_tr_sink_t sink, void *udata) {
    
    unsigned int nk, nl, nm;
    unsigned int *k, *l, *m;
    int allocated, res;

    if (nlow <= 0 || nup <= 0) {
        return -1;
    }

    allocated = cfac_overlap_if(nlow, low, nup, up, &nk, &k, &nl, &l, &nm, &m);

    /* K <-> F */
    res = crac_save_rtrans0(cfac, nk, k, nup, up, mpole, mode, sink, udata);
    if (res != 0) {
        return -1;
    }
    res = crac_save_rtrans0(cfac, nup, up, nk, k, mpole, mode, sink, udata);
    if (res != 0) {
        return -1;
    }
    /* M <-> M */
    res = crac_save_rtrans0(cfac, nm, m, nm, m, mpole, mode, sink, udata);
    if (res != 0) {
        return -1;
    }
    /* M <-> L */
    res = crac_save_rtrans0(cfac, nm, m, nl, l, mpole, mode, sink, udata);
    if (res != 0) {
        return -1;
    }
    res = crac_save_rtrans0(cfac, nl, l, nm, m, mpole, mode, sink, udata);
    if (res != 0) {
        return -1;
    }

    if (allocated) {
        free(k);
        free(l);
        free(m);
    }
    
    return 0;
}

int GetLowUpEB(const cfac_t *cfac, int *nlow, int **low, int *nup, int **up, 
	       int nlow0, const int *low0, int nup0, const int *up0) {
  int i, n;
 
  n = GetNumEBLevels(cfac);
  if (n == 0) return -1;

  *low = malloc(sizeof(int)*n);
  *up = malloc(sizeof(int)*n);
  *nlow = 0;
  *nup = 0;
  
  for (i = 0; i < n; i++) {
    int j, ilev, mlev;
    LEVEL *lev = GetEBLevel(cfac, i);
    DecodeBasisEB(lev->pb, &ilev, &mlev);
    for (j = 0; j < nlow0; j++) {
      if (low0[j] == ilev) {
	(*low)[(*nlow)++] = i;	
	break;
      }      
    }
    for (j = 0; j < nup0; j++) {
      if (up0[j] == ilev) {
	(*up)[(*nup)++] = i;	
	break;
      }
    }
  }

  return 0;
}

int SaveTransitionEB(cfac_t *cfac, int nlow0, int *low0, int nup0, int *up0,
		     char *fn, int m) {
  int n, nlow, *low, nup, *up, nc;

  n = GetLowUpEB(cfac, &nlow, &low, &nup, &up, nlow0, low0, nup0, up0);
  if (n == -1) return 0;
  
  trm_cache = TRMultipole_cache_new(cfac_get_num_levels(cfac));

  nc = OverlapLowUp(nlow, low, nup, up);
  SaveTransitionEB0(cfac, nc, low+nlow-nc, nc, up+nup-nc, fn, m);
  SaveTransitionEB0(cfac, nc, low+nlow-nc, nup-nc, up, fn, m);
  SaveTransitionEB0(cfac, nup-nc, up, nc, low+nlow-nc, fn, m);
  SaveTransitionEB0(cfac, nlow-nc, low, nup, up, fn, m);
  SaveTransitionEB0(cfac, nup, up, nlow-nc, low, fn, m);

  free(low);
  free(up);
  ReinitRadial(cfac, 1);
  
  TRMultipole_cache_free(trm_cache);
  trm_cache = NULL;

  return 0;
}

int PolarizeCoeff(char *ifn, char *ofn, int i0, int i1) {
  FILE *f1, *f2 = NULL;
  int n, i, t, tp, k, q, s, sp, m, mp, m2;
  double a, c, e;
  F_HEADER fh;
  TRF_HEADER h;
  TRF_RECORD r;
  EN_SRECORD *mem_en_table;
  int swp, mem_en_table_size;
  
  mem_en_table = GetMemENFTable(&mem_en_table_size);
  if (mem_en_table == NULL) {
    printf("Energy table has not been built in memory.\n");
    return -1;
  }

  f1 = fopen(ifn, "r");
  if (f1 == NULL) {
    printf("cannot open file %s\n", ifn);
    return -1;
  }  

  n = ReadFHeader(f1, &fh, &swp);
  if (n == 0) {
    printf("File %s is not in FAC Binary format\n", ifn);
    goto DONE;
  }
  
  if (fh.type != DB_TRF || fh.nblocks == 0) {
    printf("File %s is not of DB_TRF type\n", ifn);
    goto DONE;
  }
  
  if (strcmp(ofn, "-") == 0) {
    f2 = stdout;
  } else {
    f2 = fopen(ofn, "w");
  }
  if (f2 == NULL) {
    printf("cannot open file %s\n", ofn);
    goto DONE;
  }
    
  while (1) {
    n = ReadTRFHeader(f1, &h, swp);
    if (n == 0) break;
    m2 = 2*abs(h.multipole);
    for (i = 0; i < h.ntransitions; i++) {
      n = ReadTRFRecord(f1, &r, swp, &h);
      if ((r.lower == i0 || i0 < 0) && (r.upper == i1 || i1 < 0)) {
	e = mem_en_table[r.upper].energy - mem_en_table[r.lower].energy;
	e = FINE_STRUCTURE_CONST*e;
	e = e*e*e;
	e *= RATE_AU/(4.0*M_PI);
	for (t = -1; t <= 1; t += 2) {
	  for (tp = -1; tp <= 1; tp += 2) {
	    for (k = 0; k <= m2; k++) {
	      for (q = -k; q <= k; q++) {
		c = 0.0;
		for (s = 0; s <= m2; s++) {
		  m = s - m2/2;
		  mp = m + q;
		  sp = mp + m2/2;
		  if (sp < 0 || sp > m2) continue;
		  a = r.strength[s]*r.strength[sp];
		  a *= W3j(m2, m2, 2*k, 2*m, -2*mp, 2*q);
		  a *= W3j(m2, m2, 2*k, 2*t, -2*tp, 2*(tp-t));
		  if (IsOdd(abs(mp-tp))) a = -a;
		  c += a;
		}
		c *= 2.0*k + 1.0;
		c *= e;
		fprintf(f2, "%4d %4d %2d %2d %2d %2d %2d %15.8E\n",
			r.lower, r.upper, t, tp, k, q, tp-t, c);
	      }
	    }
	  }
	}
      }
      free(r.strength);
    }
  }

 DONE:
  fclose(f1);
  
  if (f2) {
    if (f2 != stdout) {
      fclose(f2);
    } else {
      fflush(f2);
    }
  }

  return 0;
}
  
int GetTransition(const cfac_t *cfac,
    int nlo, int nup, TRANSITION *tr, int *swapped)
{
    if (!tr) {
        return -1;
    }
    
    tr->llo = GetLevel(cfac, nlo);
    tr->lup = GetLevel(cfac, nup);
    if (!tr->llo || !tr->lup) {
        return -1;
    }
    
    tr->e = tr->lup->energy - tr->llo->energy;
    if (tr->e < 0) {
        LEVEL *lbuf;
        tr->e = -tr->e;
        lbuf = tr->llo;
        tr->llo = tr->lup;
        tr->lup = lbuf;
        
        tr->nup = nlo;
        tr->nlo = nup;
        
        *swapped = 1;
    } else {
        tr->nup = nup;
        tr->nlo = nlo;
        
        *swapped = 0;
    }
    
    return 0;
}
