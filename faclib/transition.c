#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "global.h"
#include "consts.h"
#include "cfacP.h"
#include "radial.h"
#include "angular.h"
#include "structure.h"
#include "dbase.h"
#include "transition.h"

/* the following options controll the computation methods of OS.
   gauge = 1 coulomb gauge (velocity form)
         = 2 babushkin gauge (length form)

   mode = 0 use relativistic expression for radial integrals.
        = 1 use non-relativistic approximation.
   
   max_e, the maximum rank of electric multipoles.
   max_m, the maximum rank of magnetic multipoles.
*/
static struct {
  int gauge;
  int mode;
  int max_e;
  int max_m;
  double eps0;
  double eps;
} transition_option = {DGAUGE, DMODE, ERANK, MRANK, TRCUT0, TRCUT};

typedef struct {
  TR_RECORD r;
  TR_EXTRA rx;
  int ks[2];
} TR_DATUM;

int SetTransitionCut(double c0, double c) {
  if (c0 >= 0) {
    transition_option.eps0 = c0;
  } else {
    transition_option.eps0 = TRCUT0;
  }
  if (c >= 0) {
    transition_option.eps = c;
  } else {
    transition_option.eps = TRCUT;
  }
  return 0;
}

double GetTransitionCut(void) {
  return transition_option.eps;
}

void SetTransitionMode(int m) {
  transition_option.mode = m;
}

void SetTransitionGauge(int m) {
  transition_option.gauge = m;
}

void SetTransitionMaxE(int m) {
  transition_option.max_e = m;
}

void SetTransitionMaxM(int m) {
  transition_option.max_m = m;
}

void SetTransitionOptions(int gauge, int mode, 
			  int max_e, int max_m) {
  transition_option.gauge = gauge;
  transition_option.mode = mode;
  transition_option.max_e = max_e;
  transition_option.max_m = max_m;
}

int GetTransitionGauge(void) {
  return transition_option.gauge;
}

int GetTransitionMode(void) {
  return transition_option.mode;
}

/* If strict flag is set, calculate only if energy > 0 */
static int _TRMultipole(double *strength, double *energy,
		int m, int lower, int upper, int strict) {
  int m2;
  int p1, p2, j1, j2;
  LEVEL *lev1, *lev2;
  double s, r, aw;
  int nz, i;
  ANGULAR_ZMIX *ang;
  
  *strength = 0.0;

  lev1 = GetLevel(cfac, lower);
  if (lev1 == NULL) return -1;
  lev2 = GetLevel(cfac, upper);
  if (lev2 == NULL) return -1;
  
  *energy = lev2->energy - lev1->energy;
  
  if (GetNumElectrons(cfac, lower) != GetNumElectrons(cfac, upper)) {
    return -1;
  }
  
  if (strict && *energy <= 0.0) {
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
    if (transition_option.mode == M_NR && m != 1) {
      r = MultipoleRadialNR(m, ang[i].k0, ang[i].k1, 
			    transition_option.gauge);
    } else {
      r = MultipoleRadialFR(aw, m, ang[i].k0, ang[i].k1,
			    transition_option.gauge);
    }
    s += r * ang[i].coeff;
  }
  if (nz > 0) {
    free(ang);	
  }  

  *strength = s;
  
  return 0;
}

typedef struct {
  int m;
  int valid;
  double energy;
  double strength;
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

/* If energy is not NULL, it is assigned trans. energy;
   if strict flag is set, calculate only if *energy > 0 */
int TRMultipole(double *strength, double *energy,
		int m, int lower, int upper, int strict) {
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
      *strength = trans->strength;
      
      return 0;
    }
  }
  
  res = _TRMultipole(strength, &dE, m, lower, upper, strict);
  if (energy) {
    *energy = dE;
  }
  
  if (trans) {
    trans->m        = m;
    trans->energy   = dE;
    trans->strength = *strength;
    trans->valid    = 1;
  }
  
  return res;
}  

int TRMultipoleEB(double *strength, double *energy, int m, int lower, int upper) {
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
      
      if (TRMultipole(&r, NULL, m, ilev1, ilev2, 0) != 0) {
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

int SaveTransitionEB0(int nlow, int *low, int nup, int *up, 
		      char *fn, int m) {
  int k, i, j, nq;
  double emin, emax, e0, s[101], et;
  F_HEADER fhdr;
  TRF_HEADER tr_hdr;
  TRF_RECORD r;
  LEVEL *lev1, *lev2;
  FILE *f;
  
  if (nlow <= 0 || nup <= 0) return -1;
  if (m == 1 || transition_option.mode == M_FR) {
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
    
    FreeMultipoleArray();
    if (e0 < EPS3) {
      SetAWGrid(1, emin, emax);
    } else if (e0 < 1.0) {
      SetAWGrid(2, emin, emax);
    } else {
      SetAWGrid(3, emin, emax);
    }
  }
  fhdr.type = DB_TRF;
  strcpy(fhdr.symbol, GetAtomicSymbol(cfac));
  fhdr.atom = GetAtomicNumber(cfac);  
  lev1 = GetEBLevel(cfac, low[0]);
  DecodeBasisEB(lev1->pb, &i, &j);  
  tr_hdr.nele = GetNumElectrons(cfac, i);
  tr_hdr.multipole = m;
  tr_hdr.gauge = GetTransitionGauge();
  if (m == 1) { /* always FR for M1 transitions */
    tr_hdr.mode = M_FR;
  } else {
    tr_hdr.mode = GetTransitionMode();
  }
  nq = 2*abs(m) + 1;
  r.strength = (float *) malloc(sizeof(float)*nq);
  GetFields(cfac, &tr_hdr.bfield, &tr_hdr.efield, &tr_hdr.fangle);
    
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);

  for (j = 0; j < nup; j++) {
    for (i = 0; i < nlow; i++) {
      k = TRMultipoleEB(s, &et, m, low[i], up[j]);
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
      
int SaveTransition0(int nlow, int *low, int nup, int *up, 
		    char *fn, int m) {
  int i, j, k, jup;
  FILE *f;
  LEVEL *lev1, *lev2;
  TR_RECORD r;
  TR_HEADER tr_hdr;
  F_HEADER fhdr;
  double *s, *et, *a, trd, gf;
  double e0, emin, emax;

#ifdef PERFORM_STATISTICS
  STRUCT_TIMING structt;
  ANGULAR_TIMING angt;
  RECOUPLE_TIMING recouplet;
  RAD_TIMING radt;
#endif
  
  if (nlow <= 0 || nup <= 0) return -1;
  k = 0;
  emin = 1E10;
  emax = 1E-10;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(cfac, low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(cfac, up[j]);
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
    
  FreeMultipoleArray();
  if (e0 < EPS3) {
    SetAWGrid(1, emin, emax);
  } else if (e0 < 1.0) {
    SetAWGrid(2, emin, emax);
  } else {
    SetAWGrid(3, emin, emax);
  }
  
  fhdr.type = DB_TR;
  strcpy(fhdr.symbol, GetAtomicSymbol(cfac));
  fhdr.atom = GetAtomicNumber(cfac);
  tr_hdr.nele = GetNumElectrons(cfac, low[0]);
  tr_hdr.multipole = m;
  tr_hdr.gauge = GetTransitionGauge();
  if (m == 1) { /* always FR for M1 transitions */
    tr_hdr.mode = M_FR;
  } else {
    tr_hdr.mode = GetTransitionMode();
  }
  f = OpenFile(fn, &fhdr);
  InitFile(f, &fhdr, &tr_hdr);
    
  a = malloc(sizeof(double)*nlow);
  s = malloc(sizeof(double)*nlow);
  et = malloc(sizeof(double)*nlow);
  for (j = 0; j < nup; j++) {
    jup = LevelTotalJ(cfac, up[j]);
    trd = 0.0;
    for (i = 0; i < nlow; i++) {
      a[i] = 0.0;
      k = TRMultipole(s+i, et+i, m, low[i], up[j], 1);
      if (k != 0) continue;
      gf = OscillatorStrength(m, et[i], s[i], &(a[i]));
      a[i] /= jup+1.0;
      trd += a[i];
    } 
    if (trd < 1E-30) continue;
    r.upper = up[j];
    for (i = 0; i < nlow; i++) {
      if (a[i] < (transition_option.eps * trd)) continue;
      r.lower = low[i];
      r.strength = s[i];
      WriteTRRecord(f, &r, NULL);
    }
  }

  free(a);
  free(s);
  free(et);

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

#ifdef PERFORM_STATISTICS
  GetStructTiming(&structt);
  fprintf(perform_log, "AngZMix: %6.1E, AngZFB: %6.1E, AngZxZFB: %6.1E, SetH: %6.1E DiagH: %6.1E\n",
	  ((double) (structt.angz_mix))/CLOCKS_PER_SEC,
	  ((double) (structt.angz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxz_fb))/CLOCKS_PER_SEC,
	  ((double) (structt.set_ham))/CLOCKS_PER_SEC,
	  ((double) (structt.diag_ham))/CLOCKS_PER_SEC);
  fprintf(perform_log, "AngZS: %ld %ld %6.1E, %6.1E, AngZFBS: %6.1E, AngZxZFBS: %6.1E, AddZ: %6.1E, AddZxZ: %6.1E\n",
	  structt.n_angz_states, structt.n_angz_states_load,
	  ((double) (structt.angz_states))/CLOCKS_PER_SEC, 
	  ((double) (structt.angz_states_load))/CLOCKS_PER_SEC,
	  ((double) (structt.angzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.angzxzfb_states))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angz))/CLOCKS_PER_SEC,
	  ((double) (structt.add_angzxz))/CLOCKS_PER_SEC);

  GetAngularTiming(&angt);
  fprintf(perform_log, "W3J: %6.1E, W6J: %6.1E, W9J: %6.1E\n", 
	  ((double)angt.w3j)/CLOCKS_PER_SEC, 
	  ((double)angt.w6j)/CLOCKS_PER_SEC, 
	  ((double)angt.w9j)/CLOCKS_PER_SEC);
  GetRecoupleTiming(&recouplet);
  fprintf(perform_log, "AngZ: %6.1E, AngZxZ: %6.1E, Interact: %6.1E\n",
	  ((double)recouplet.angz)/CLOCKS_PER_SEC,
	  ((double)recouplet.angzxz)/CLOCKS_PER_SEC,
	  ((double)recouplet.interact)/CLOCKS_PER_SEC);
  GetRadTiming(&radt);
  fprintf(perform_log, "Dirac: %d, %6.1E, 1E: %6.1E, Slater: %6.1E, 2E: %6.1E\n", 
	  GetNumContinua(),
	  ((double)radt.dirac)/CLOCKS_PER_SEC, 
	  ((double)radt.radial_1e)/CLOCKS_PER_SEC,
	  ((double)radt.radial_slater)/CLOCKS_PER_SEC,
	  ((double)radt.radial_2e)/CLOCKS_PER_SEC);
  fprintf(perform_log, "\n");
  fflush(perform_log);
#endif /* PERFORM_STATISTICS */

  return 0;
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
  
int SaveTransition(int nlow, int *low, int nup, int *up,
		   char *fn, int m) {
  int n, *alev = NULL, i, nc;
  
  n = 0;
  if (nlow == 0 || nup == 0) {
    n = GetNumLevels(cfac);
    if (n <= 0) return -1;
    alev = malloc(sizeof(int)*n);
    if (!alev) return -1;
    
    for (i = 0; i < n; i++) alev[i] = i;

    if (nlow == 0) {
      nlow = n; 
      low = alev;
    }
    if (nup == 0) {
      nup = n;
      up = alev;
    }
  }
  if (nlow <= 0 || nup <= 0) return -1;

  nc = OverlapLowUp(nlow, low, nup, up);
  SaveTransition0(nc, low+nlow-nc, nc, up+nup-nc, fn, m);
  SaveTransition0(nc, low+nlow-nc, nup-nc, up, fn, m);
  SaveTransition0(nup-nc, up, nc, low+nlow-nc, fn, m);
  SaveTransition0(nlow-nc, low, nup, up, fn, m);
  SaveTransition0(nup, up, nlow-nc, low, fn, m);

  if (n > 0) free(alev);
  ReinitRadial(1);

  return 0;
}
  
int GetLowUpEB(int *nlow, int **low, int *nup, int **up, 
	       int nlow0, int *low0, int nup0, int *up0) {  
  int i, j, ilev, mlev, n;
  LEVEL *lev;
 
  n = GetNumEBLevels(cfac);
  if (n == 0) return -1;

  *low = malloc(sizeof(int)*n);
  *up = malloc(sizeof(int)*n);
  *nlow = 0;
  *nup = 0;
  for (i = 0; i < n; i++) {
    lev = GetEBLevel(cfac, i);
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

int SaveTransitionEB(int nlow0, int *low0, int nup0, int *up0,
		     char *fn, int m) {
  int n, nlow, *low, nup, *up, nc;

  n = GetLowUpEB(&nlow, &low, &nup, &up, nlow0, low0, nup0, up0);
  if (n == -1) return 0;
  
  trm_cache = TRMultipole_cache_new(GetNumLevels(cfac));

  nc = OverlapLowUp(nlow, low, nup, up);
  SaveTransitionEB0(nc, low+nlow-nc, nc, up+nup-nc, fn, m);
  SaveTransitionEB0(nc, low+nlow-nc, nup-nc, up, fn, m);
  SaveTransitionEB0(nup-nc, up, nc, low+nlow-nc, fn, m);
  SaveTransitionEB0(nlow-nc, low, nup, up, fn, m);
  SaveTransitionEB0(nup, up, nlow-nc, low, fn, m);

  free(low);
  free(up);
  ReinitRadial(1);
  
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
  
