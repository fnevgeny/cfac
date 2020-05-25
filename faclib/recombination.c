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

#include "cfacP.h"
#include "consts.h"
#include "parser.h"
#include "coulomb.h"
#include "structure.h"
#include "radial.h"
#include "angular.h"
#include "grid.h"
#include "interpolation.h"
#include "transition.h"
#include "dbase.h"
#include "recombination.h"
#include "coulrad.h"

static int qk_mode;
static double qk_fit_tolerance;

static int egrid_type = -1;
static int usr_egrid_type = -1;

static int n_egrid = 0;
static double egrid[MAXNE];
static double log_egrid[MAXNE];
static double xegrid[MAXNE];
static double log_xegrid[MAXNE];
static double egrid_min;
static double egrid_max;
static int egrid_limits_type = 0;

static int n_usr = 0;
static double usr_egrid[MAXNUSR];
static double log_usr[MAXNUSR];

static int n_tegrid = 0;
static double tegrid[MAXNTE];
static double log_te[MAXNTE];

static MULTI *pk_array;
static MULTI *qk_array;

#define MAXAIM 1024
#define NPARAMS 3
static ARRAY *hyd_qk_array;

static struct {
  int n_spec;
  int n_frozen;
  int max_kl;
  int kl_interp;
  int nkl0;
  int nkl;
  int pw_limits[2];
  int kl[MAXNKL+1];
  int kappa0[(MAXNKL+1)*2];
} pw_scratch = {RECNSPEC, RECNFROZEN, 
		RECLMAX, RECLMAX,
		0, 0, {0, RECLMAX}, {0.0}, {0.0}};

double ai_cut = AICUT;

static REC_COMPLEX rec_complex[MAX_COMPLEX];
int n_complex = 0;

static void FreeRecPkData(void *p) {
  double *dp;
  dp = *((double **) p);
  free(dp);
  *((double **) p) = NULL;
}

int SetAICut(double c) {
  ai_cut = c;
  return 0;
}

int SetRRTEGrid(int n, double emin, double emax) {
  n_tegrid = SetTEGrid(tegrid, log_te, n, emin, emax);
  return n_tegrid;
}

int SetRRTEGridDetail(int n, double *x) {
  n_tegrid = SetTEGridDetail(tegrid, log_te, n, x);
  return n_tegrid;
}

int SetUsrPEGridType(int type) {
  if (type >= 0) usr_egrid_type = type;
  return 0;
}

int SetPEGridLimits(double min, double max, int type) {
  if (min <= 0) egrid_min = 0.05;
  else egrid_min = min;
  if (max <= 0) egrid_max = 8.0;
  else egrid_max = max;
  egrid_limits_type = type;
  return 0;
}

int SetPEGridDetail(int n, double *xg) {
  n_egrid = SetEGridDetail(egrid, log_egrid, n, xg);
  return n_egrid;
}

int SetPEGrid(int n, double emin, double emax, double eth) {
  n_egrid = SetEGrid(egrid, log_egrid, n, emin, emax, eth);
  return n_egrid;
}

int SetUsrPEGridDetail(int n, double *xg) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGridDetail(usr_egrid, log_usr, n, xg);
  return n_usr;
}
  						      
int SetUsrPEGrid(int n, double emin, double emax, double eth) {
  if (n > MAXNUSR) {
    printf("Max # of grid points reached \n");
    return -1;
  }
  n_usr = SetEGrid(usr_egrid, log_usr, n, emin, emax, eth);
  return n_usr;
}

int AddRecPW(int n, int step) {
  int i, i2, kl, kl2;
  for (i = pw_scratch.nkl0; i < n+pw_scratch.nkl0; i++) {
    i2 = i*2;
    kl = pw_scratch.kl[i-1] + step;
    kl2 = kl*2;
    pw_scratch.kl[i] = kl;
    pw_scratch.kappa0[i2] = GetKappaFromJL(kl2-1, kl2);
    pw_scratch.kappa0[i2+1] = GetKappaFromJL(kl2+1, kl2);
    if (kl > pw_scratch.max_kl) break;
  }
  pw_scratch.nkl0 = i;
  return 0;
}

int SetRecSpectator(int n_frozen, int n_spec) {
  if (n_frozen > 0) pw_scratch.n_frozen = n_frozen;
  if (n_spec > 0) pw_scratch.n_spec = n_spec;
  return 0;
}

int SetRecQkMode(int m, double tol) {
  if (m == QK_DEFAULT) qk_mode = QK_FIT;
  else qk_mode = m;
  if (tol > 0.0) qk_fit_tolerance = tol;
  return 0;
}

int SetRecPWOptions(int kl_interp, int max_kl) {
  int k, j, m;

  pw_scratch.kl_interp = kl_interp;
  pw_scratch.max_kl = max_kl;
  pw_scratch.kl[0] = 0;
  pw_scratch.nkl0 = 1;
  pw_scratch.kappa0[0] = 0;
  pw_scratch.kappa0[1] = -1;

  AddRecPW(pw_scratch.kl_interp, 1);
  k = 2;
  j = 2;
  m = pw_scratch.kl[pw_scratch.nkl0-1];
  while (m+k <= max_kl) {
    AddRecPW(j, k);
    m = pw_scratch.kl[pw_scratch.nkl0-1];
    k *= 2;
  }
  pw_scratch.nkl = pw_scratch.nkl0;
  return 0;
}

int SetRecPWLimits(int min, int max) {
  pw_scratch.pw_limits[0] = min;
  pw_scratch.pw_limits[1] = max;
  return 0;
}

int ConstructRecGroupName(char *rgn, char *gn, int n) {
  sprintf(rgn, "+%d__%s", n, gn);
  return 0;
}

int IsRecombinedGroup(cfac_t *cfac, int i) {
  char *s;
  int n;
  s = GetGroup(cfac, i)->name;
  if (s[0] != '+') return 0;
  n = strtol(s, NULL, 10);
  return n;
}

int RecStates(cfac_t *cfac, int n, int k, int *kg, char *fn) {
  int i, j, m, nsym, nlevels, ncfgs, kg0, t;
  ARRAY *clist;
  CONFIG *rcfg, *c;
  SHELL ns;
  CONFIG_GROUP *g;
  char *gn, rgn[GROUP_NAME_LEN + 20];
  int nm;

  nm = 0;
  for (i = 0; i < k; i++) {
    g = GetGroup(cfac, kg[i]);
    clist = &(g->cfg_list);
    for (t = 0; t < g->n_cfgs; t++) {
      c = (CONFIG *) ArrayGet(clist, t);
      if (c->shells[0].n > nm) nm = c->shells[0].n;
    }
  }
  if (n < nm) return 0;

  nm++;

  if (pw_scratch.n_spec > nm) nm = pw_scratch.n_spec;

  if (n >= nm) {
    i = RecStatesFrozen(cfac, n, k, kg, fn);
    return i;
  }

  ns.n = n;
  ns.nq = 1;
  ncfgs = 0;
  for (i = 0; i < k; i++) {
    kg0 = kg[i];
    g = GetGroup(cfac, kg0);
    gn = g->name;
    ConstructRecGroupName(rgn, gn, n);
    if ((kg[i] = GroupExists(cfac, rgn)) >= 0) continue;
    kg[i] = AddGroup(cfac, rgn);
    if (kg[i] < 0) {
      printf("Can not add more Groups\n");
      return -2;
    }
    ArrayAppend(rec_complex[n_complex].rg, kg+i);
    clist = &(g->cfg_list);
    for (t = 0; t < g->n_cfgs; t++) {
      c = (CONFIG *) ArrayGet(clist, t);
      for (j = 1; j < 2*pw_scratch.nkl0; j++) {
	if (pw_scratch.kl[j/2] >= n) break;
	ns.kappa = pw_scratch.kappa0[j];
	m = CompareShell(c->shells, &ns);
	if (m > 0) continue;
	if (m == 0) {
	  if (ShellClosed(c->shells)) continue;
	  else {
	    rcfg = malloc(sizeof(CONFIG));
	    rcfg->n_shells = c->n_shells;
	    rcfg->shells = malloc(sizeof(SHELL)*rcfg->n_shells);
	    memcpy(rcfg->shells, c->shells, sizeof(SHELL)*rcfg->n_shells);
	    rcfg->shells[0].nq += 1;
	  }
	} else {
	  rcfg = malloc(sizeof(CONFIG));
	  rcfg->n_shells = c->n_shells + 1;
	  rcfg->shells = malloc(sizeof(SHELL)*rcfg->n_shells);
	  memcpy(rcfg->shells+1, c->shells, sizeof(SHELL)*c->n_shells);
	  memcpy(rcfg->shells, &ns, sizeof(SHELL));
	}
	
	if (Couple(rcfg) < 0) return -3;
	if (AddConfigToList(cfac, kg[i], rcfg) < 0) return -4;
	ncfgs++;
      }
    }
  }
  if (ncfgs == 0) return -1;
  nlevels = cfac_get_num_levels(cfac);
  nsym = MAX_SYMMETRIES;
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    HAMILTON *h = ConstructHamilton(cfac, i, k, kg, 0, NULL);
    if (!h) continue;
    j = DiagonalizeHamilton(cfac, h);
    if (!j) {
      if (AddToLevels(cfac, h, 0, NULL) != 0) {
        return -2;
      }
    }
    cfac_hamiltonian_free(h);
    if (j < 0) return -1;
  }
  rec_complex[n_complex].s1 = cfac_get_num_levels(cfac)-1;
  n_complex++;
  SortLevels(cfac, nlevels, -1, 0);
  SaveLevels(cfac, fn, nlevels, -1);

  return 0;
}

int RecStatesFrozen(cfac_t *cfac, int n, int k, int *kg, char *fn) {
  int i, j, m, nlevels, nsym, nstates, ko;
  int kl2, j1, j2, p1, p, jmin, jmax, tj;
  LEVEL *lev;
  STATE *s;
  SYMMETRY *sym;
  int i0, i1, t, nt;

  nstates = 0;

  i0 = 0;
  if (n_complex == 0) {
    nt = 1;
    i1 = cfac_get_num_levels(cfac);
    rec_complex[0].s0 = i1;
  } else {
    nt = n_complex;
  }
  for (t = 0; t < nt; t++) {
    i1 = rec_complex[t].s0;
    for (i = i0; i < i1; i++) {
      lev = GetLevel(cfac, i);
      m = lev->pb;
      sym = GetSymmetry(cfac, lev->pj);
      s = GetSymmetryState(sym, m);
      if (!InGroups(s->kgroup, k, kg)) {
	continue;
      }
      j1 = lev->pj;
      DecodePJ(j1, &p1, &j1);
      
      for (j = 1; j < 2*pw_scratch.nkl0; j++) {
	kl2 = pw_scratch.kl[j/2];
	p = p1 + kl2;
	if (kl2 >= n) break;
	ko = OrbitalIndex(cfac, n, pw_scratch.kappa0[j], 0.0);
	j2 = GetJFromKappa(pw_scratch.kappa0[j]);
	jmin = abs(j2 - j1);
	jmax = j2 + j1;
	for (tj = jmin; tj <= jmax; tj += 2) {
	  if (AddStateToSymmetry(cfac, -(i+1), ko, tj, p, tj) != 0) {
            return -1;
          }
	  nstates++;
	}
      }
    }  
    i0 = rec_complex[t].s1+1;
  }

  if (nstates == 0) return -1;
  nsym = MAX_SYMMETRIES;
  nlevels = cfac_get_num_levels(cfac);
  rec_complex[n_complex].n = n;
  rec_complex[n_complex].s0 = nlevels;
  for (i = 0; i < nsym; i++) {
    HAMILTON *h;
    if (n >= pw_scratch.n_frozen) {
      i0 = 0;
      for (t = 0; t < nt; t++) {
	i1 = rec_complex[t].s0;
	for (j = i0; j < i1; j++) {
	  lev = GetLevel(cfac, j);
	  m = lev->pb;
	  sym = GetSymmetry(cfac, lev->pj);
	  s = GetSymmetryState(sym, m);
	  if (!InGroups(s->kgroup, k, kg)) continue;
	  h = ConstructHamiltonFrozen(cfac, i, j, NULL, n, 0, NULL);
	  if (!h) continue;
	  if (DiagonalizeHamilton(cfac, h) < 0) return -2;
	  if (AddToLevels(cfac, h, 0, NULL) != 0) {
            return -2;
          }
	}
	i0 = rec_complex[t].s1+1;
      }
    } else {
      h = ConstructHamiltonFrozen(cfac, i, k, kg, n, 0, NULL);
      if (!h) continue;
      if (DiagonalizeHamilton(cfac, h) == 0) {
        if (AddToLevels(cfac, h, 0, NULL) != 0) {
          return -2;
        }
        cfac_hamiltonian_free(h);
      } else {
        cfac_hamiltonian_free(h);
        return -2;
      }
    }
  }
  rec_complex[n_complex].s1 = cfac_get_num_levels(cfac)-1;
  n_complex++;
  SortLevels(cfac, nlevels, -1, 0);
  SaveLevels(cfac, fn, nlevels, -1);

  return 0;
}

void RRRadialQkHydrogenicParams(int np, double *p, 
				double z0, int n, int kl) {
#define NNE 12
  double **qk;
  double *t;
  double d, eth;
  double pnc[NNE], e[NNE];
  static double x[NNE], logx[NNE];
  double fvec[NNE], fjac[NNE*NPARAMS];
  double tol;
  int i, j, ne;
  static int iopt = 2;

  qk = ArraySet(hyd_qk_array, n, NULL);
  if (*qk == NULL) {
    gsl_coulomb_fb *rfb;
    tol = 1E-4;
    *qk = malloc(sizeof(double)*n*np);
    ne = NNE;
    if (iopt == 2) {
      x[0] = 1.01;
      logx[0] = log(x[0]);
      logx[ne-1] = log(80.0);
      d = (logx[ne-1] - logx[0])/(ne-1.0);
      for (i = 1; i < ne; i++) {
	logx[i] = logx[i-1] + d;
	x[i] = exp(logx[i]);
      }
      i = 200;
      iopt = 0;
    }
    
    eth = 0.5/(n*n);
    for (i = 0; i < ne; i++) {
      e[i] = (x[i]-1.0)*eth;
    }
    
    rfb = gsl_coulomb_fb_alloc(n, ne, e);
    
    t = *qk;    
    for (i = 0; i < n; i++) {
      for (j = 0; j < ne; j++) {
        pnc[j] = gsl_coulomb_fb_get_dfdE(rfb, i, -1, j)/x[j];
      }
      t[0] = pnc[0];
      t[1] = 3.0*(i+1);
      t[2] = 1.0;
      NLSQFit(np, t, tol, fvec, fjac, 
		     ne, x, logx, pnc, pnc, RRRadialQkFromFit, &i);
      t += np;
    }
    
    gsl_coulomb_fb_free(rfb);
  }

  /* FIXME: reduced mass */
  t = (*qk) + kl*np;
  p[0] = t[0]/(z0*z0);
  for (i = 1; i < np; i++) {
    p[i] = t[i];
  }
#undef NNE
}  

void RRRadialQkFromFit(int np, double *p, int n, double *x, double *logx, 
		       double *y, double *dy, int ndy, void *extra) {
  int kl, i, k;
  double t, s;

  kl = *((int *) extra);
  for (i = 0; i < n; i++) {
    t = (1 + p[2])/(sqrt(x[i]) + p[2]);
    if (t <= 0.0) {
      s = 0.0;
    } else {
      s = pow(x[i], - 4.5 - kl)*pow(t*sqrt(x[i]), p[1]);
    }
    
    if (ndy <= 0) {
      y[i] = p[0]*s;
    } else {
      k = i;
      dy[k] = s;

      k += ndy;
      dy[k] = p[0]*s*(logx[i]/2 + log(t));
      if (np == 3) {
	k += ndy;
	dy[k] = p[0]*s*p[1]*(sqrt(x[i]) - 1)/((1 + p[2])*(sqrt(x[i]) + p[2]));
      }
    }
  }
}

int RRRadialMultipoleTable(cfac_t *cfac, double *qr, int k0, int k1, int m) {
  int index[3], k, nqk;
  double **p, *qk;
  int kappaf, jf, klf, kf;
  int ite, ie, i;
  double aw, e, pref;
  int mode, gauge;

  klf = k1;
  if (IsOdd(klf)) {
    klf++;
    jf = klf-1;
  } else {
    jf = klf+1;
  }
  kappaf = GetKappaFromJL(jf, klf);

  k = 2*abs(m);
  index[0] = k0;
  index[1] = k1;
  if (m >= 0) {
    index[2] = k;
  } else {
    index[2] = k-1;
  }

  nqk = n_tegrid*n_egrid;
  p = (double **) MultiSet(qk_array, index, NULL);
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    return 0;
  }
  
  gauge = GetTransitionGauge(cfac);
  mode = GetTransitionMode(cfac);

  *p = (double *) malloc(sizeof(double)*nqk);
  
  qk = *p;
  /* the factor 2 comes from the conitinuum norm */
  pref = sqrt(2.0);  
  for (ite = 0; ite < n_tegrid; ite++) {
    aw = FINE_STRUCTURE_CONST * tegrid[ite];        
    for (ie = 0; ie < n_egrid; ie++) {
      e = egrid[ie];
      kf = OrbitalIndex(cfac, 0, kappaf, e);
      if (mode == M_NR && m != 1) {
	qk[ie] = MultipoleRadialNR(cfac, m, k0, kf, gauge);
      } else {
	qk[ie] = MultipoleRadialFR(cfac, aw, m, k0, kf, gauge);
      }
      qk[ie] *= pref;
    }
    qk += n_egrid;
  }
  
  for (i = 0; i < nqk; i++) {
    qr[i] = (*p)[i];
  }
  return 0;
}
    
int RRRadialQkTable(cfac_t *cfac, double *qr, int k0, int k1, int m) {
  int index[3], k, nqk;
  double **p, *qk, tq[MAXNE];
  double r0, r1, tq0[MAXNE];
  ORBITAL *orb;
  int kappa0, jb0, klb02, klb0;
  int kappaf, jf, klf, kf;
  int jfmin, jfmax;
  int ite, ie, i;
  double eb, aw, e, pref;
  int mode, gauge;

  orb = GetOrbital(cfac, k0);
  kappa0 = orb->kappa;
  GetJLFromKappa(kappa0, &jb0, &klb02);
  klb0 = klb02/2;
  eb = -orb->energy;

  for (ie = 0; ie < n_egrid; ie++) {
    xegrid[ie] = 1.0 + egrid[ie]/eb;
    log_xegrid[ie] = log(xegrid[ie]);
  }

  if (m == -1) {
    int nh, klh;
    GetHydrogenicNL(cfac, &nh, &klh, NULL, NULL);
    if (klb0 > klh || orb->n > nh) {
      double hparams[NPARAMS];
      double z;
      
      if (k0 != k1) {
        return -1;
      }
      
      z = GetResidualZ(cfac);
      RRRadialQkHydrogenicParams(NPARAMS, hparams, z, orb->n, klb0);
      RRRadialQkFromFit(NPARAMS, hparams, n_egrid, 
		        xegrid, log_xegrid, tq0, NULL, 0, &klb0);
      for (ie = 0; ie < n_egrid; ie++) {
	qr[ie] = tq0[ie]/eb;
      }
      k = n_egrid;
      for (ite = 1; ite < n_tegrid; ite++) {
	for (ie = 0; ie < n_egrid; ie++) {
	  qr[k] = qr[ie];
	  k++;
	}
      }
      return 0;
    }
  }
  
  k = 2*abs(m);
  index[0] = k0;
  index[1] = k1;
  if (m >= 0) {
    index[2] = k;
  } else {
    index[2] = k-1;
  }

  nqk = n_tegrid*n_egrid;
  p = MultiSet(qk_array, index, NULL);
  if (*p) {
    for (i = 0; i < nqk; i++) {
      qr[i] = (*p)[i];
    }
    return 0;
  }

  gauge = GetTransitionGauge(cfac);
  mode = GetTransitionMode(cfac);

  *p = malloc(sizeof(double)*nqk);
  
  qk = *p;
  /* the factor 2 comes from the conitinuum norm */
  pref = 2.0/((k+1.0)*(jb0+1.0));
  
  for (ite = 0; ite < n_tegrid; ite++) {
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] = 0.0;
    }
    aw = FINE_STRUCTURE_CONST * tegrid[ite];
    jfmin = jb0 - k;
    jfmax = jb0 + k;
    for (jf = jfmin; jf <= jfmax; jf += 2) {
      for (klf = jf-1; klf <= jf+1; klf += 2) {
	if (jf <= 0 ||
	    klf < 0 ||
	    (m < 0 && IsOdd((klb02+klf+k)/2)) ||
	    (m > 0 && IsEven((klb02+klf+k)/2))) {
	  continue;
	}
	kappaf = GetKappaFromJL(jf, klf);
	for (ie = 0; ie < n_egrid; ie++) {
	  e = egrid[ie];
	  kf = OrbitalIndex(cfac, 0, kappaf, e);
	  if (mode == M_NR && m != 1) {
	    r0 = MultipoleRadialNR(cfac, m, k0, kf, gauge);
	    if (k1 == k0) {
	      r1 = r0;
	    } else {
	      r1 = MultipoleRadialNR(cfac, m, k1, kf, gauge);
	    }
	  } else {
	    r0 = MultipoleRadialFR(cfac, aw, m, k0, kf, gauge);
	    if (k1 == k0) {
	      r1 = r0;
	    } else {
	      r1 = MultipoleRadialFR(cfac, aw, m, k1, kf, gauge);
	    }
	  }	  
	  tq[ie] += r0*r1;
	}  
      }
    }

    for (ie = 0; ie < n_egrid; ie++) {
      qk[ie] = tq[ie]*pref;
    }
    qk += n_egrid;
  }
  
  for (i = 0; i < nqk; i++) {
    qr[i] = (*p)[i];
  }
  return 0;
}
  
int RRRadialMultipole(cfac_t *cfac, double *rqc, double te, int k0, int k1, int m) {
  int i, j, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXNE];

  i = RRRadialMultipoleTable(cfac, rqe, k0, k1, m);
  if (i < 0) return -1;
  
  if (n_tegrid == 1) {
    for (i = 0; i < n_egrid; i++) {
      rqc[i] = rqe[i];
    }
  } else {
    nd = 1;
    x0 = te;
    for (i = 0; i < n_egrid; i++) {
      j = i;
      for (k = 0; k < n_tegrid; k++) {
	rq[k] = rqe[j];
	j += n_egrid;
      }
      UVIP3P(n_tegrid, tegrid, rq, nd, &x0, &rqc[i]);
    }
  }
  return 0;
}
  
int RRRadialQk(cfac_t *cfac, double *rqc, double te, int k0, int k1, int m) {
  int i, j, nd, k;
  double rq[MAXNTE];
  double x0, rqe[MAXNTE*MAXNE];

  i = RRRadialQkTable(cfac, rqe, k0, k1, m);
  if (i < 0) return -1;

  if (n_tegrid == 1) {
    for (i = 0; i < n_egrid; i++) {
      rqc[i] = rqe[i];
    }
  } else {
    nd = 1;
    x0 = te;
    for (i = 0; i < n_egrid; i++) {
      j = i;
      for (k = 0; k < n_tegrid; k++) {
	rq[k] = rqe[j];
	j += n_egrid;
      }
      UVIP3P(n_tegrid, tegrid, rq, nd, &x0, &rqc[i]);
    }
  }
  return 0;
}

int BoundFreeMultipole(cfac_t *cfac, FILE *fp, int rec, int f, int m) {
  LEVEL *lev1, *lev2;
  ANGULAR_ZFB *ang;
  ORBITAL *orb;
  int i, nz, ie, k, kb, j1, j2, n, p1, p2;
  int jmin, jmax, jt, jfmin, jfmax, jf, klf, kf, jb, klb;
  double rq[MAXNE], rqu[MAXNE], eb, a;

  lev1 = GetLevel(cfac, rec);
  lev2 = GetLevel(cfac, f);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, &p1, &j1);
  DecodePJ(j2, &p2, &j2);

  eb = (lev2->energy - lev1->energy);
  if (eb <= 0.0) return -1;
  
  eb = (lev2->energy - lev1->energy);
  if (eb <= 0.0) return -1;

  nz = AngularZFreeBound(cfac, &ang, f, rec);
  if (nz <= 0) return -1;

  k = abs(m)*2;
  jmax = j1 + k;
  jmin = abs(j1 - k);
  for (jt = jmin; jt <= jmax; jt += 2) {
    jfmin = abs(jt - j2);
    jfmax = jt + j2;
    for (jf = jfmin; jf <= jfmax; jf += 2) {
      for (klf = jf-1; klf <= jf+1; klf += 2) {    
	kf = klf;
	if (jf < klf) kf--;	
	for (ie = 0; ie < n_egrid; ie++) {
	  rqu[ie] = 0.0;
	}
	if (jf <= 0 ||
	    klf < 0 ||
	    (m < 0 && IsOdd(p1+p2+klf/2+k/2)) ||
	    (m > 0 && IsEven(p1+p2+klf/2+k/2))) {
	  continue;
	}
	for (i = 0; i < nz; i++) {
	  kb = ang[i].kb;
	  orb = GetOrbital(cfac, kb);
	  GetJLFromKappa(orb->kappa, &jb, &klb);
	  if ((m < 0 && IsOdd((klb+klf+k)/2)) ||
	      (m > 0 && IsEven((klb+klf+k)/2))) {
	    continue;
	  }
	  n = RRRadialMultipole(cfac, rq, eb, kb, kf, m);
	  if (n < 0) continue;
	  a = ang[i].coeff*sqrt(jt+1.0)*W6j(j2, jf, jt, k, j1, jb);
	  if (IsOdd((jt+j1-k)/2)) a = -a;
	  if (fabs(a) > EPS16) {
	    for (ie = 0; ie < n_egrid; ie++) {
	      rqu[ie] += a * rq[ie];
	    }
	  }
	}
	for (ie = 0; ie < n_egrid; ie++) {
	  fprintf(fp, "%5d %2d %5d %2d  %2d  %3d %3d %3d %12.5E %12.5E %12.5E\n",
		  rec, j1, f, j2, m, jt, klf, jf, eb, egrid[ie], rqu[ie]);
	}
      }
    }
  }

  fprintf(fp, "\n");
  free(ang);

  return 0;
}

int BoundFreeOS(cfac_t *cfac, double *rqu, double *rqc, double *eb, 
		int rec, int f, int m, int iuta) {
  LEVEL *lev1, *lev2;
  ORBITAL *orb = NULL;
  double rq[MAXNE], tq[MAXNE];
  double a, b, d, eb0 = 0.0, z;
  int nkl = 0, nq = 0, k;
  int ie, c;
  int kb = 0, jb, klb;
  
  ANGULAR_ZFB *ang;
  int i, j, kbp, jbp, nz = 0;
  double amax;

  INTERACT_DATUM *idatum;
  int j1 = 0, ns, q1 = 0;

  lev1 = GetLevel(cfac, rec);
  lev2 = GetLevel(cfac, f);

  *eb = (lev2->energy - lev1->energy);
  if (*eb <= 0.0) return -1;

  if (iuta) {
    idatum = NULL;
    ns = GetInteract(cfac, &idatum, NULL, NULL, lev2->uta_cfg_g, lev1->uta_cfg_g,
		     lev2->uta_g_cfg, lev1->uta_g_cfg, 0, 0, 1);  
    if (ns <= 0) return -1;
    if (idatum->s[1].index < 0 || idatum->s[3].index >= 0) {
      free(idatum->bra);
      free(idatum);
      return -1;
    }

    j1 = idatum->s[1].j;
    q1 = idatum->s[1].nq_ket;
    kb = OrbitalIndex(cfac, idatum->s[1].n, idatum->s[1].kappa, 0.0);
    orb = GetOrbital(cfac, kb);
    eb0 = -(orb->energy);
    GetJLFromKappa(orb->kappa, &jb, &klb);
    klb /= 2;
  } else {
    nz = AngularZFreeBound(cfac, &ang, f, rec);
    if (nz <= 0) return -1;
  }

  c = 2*abs(m) - 2;

  for (ie = 0; ie < n_egrid; ie++) {
    tq[ie] = 0.0;
  }
  
  if (iuta) {
    k = RRRadialQk(cfac, rq, *eb, kb, kb, m);
    if (k < 0) {
      free(idatum->bra);
      free(idatum);
      return -1;
    };

    nq = orb->n;
    nkl = klb;
    for (ie = 0; ie < n_egrid; ie++) {
      tq[ie] += (jb+1.0)*rq[ie];
    }
  } else {
    amax = 0.0;
    for (i = 0; i < nz; i++) {
      kb = ang[i].kb;
      orb = GetOrbital(cfac, kb);
      jbp = orb->kappa;
      GetJLFromKappa(jbp, &jb, &klb);
      klb /= 2;
      for (j = 0; j <= i; j++) {
        kbp = ang[j].kb;
        jbp = GetOrbital(cfac, kbp)->kappa;
        jbp = GetJFromKappa(jbp);
        if (jbp != jb) continue;
        k = RRRadialQk(cfac, rq, *eb, kb, kbp, m);
        if (k < 0) continue;
        a = ang[i].coeff*ang[j].coeff;
        if (j != i) {
	  a *= 2;
        } else {
	  if (a > amax) {
	    nkl = klb;
	    nq = orb->n; 
	    amax = a;
	    eb0 = -(orb->energy);
	  }
        }
        for (ie = 0; ie < n_egrid; ie++) {
	  tq[ie] += a*rq[ie];
        }
      }
    }
  }

  if (qk_mode == QK_FIT) {
    z = GetResidualZ(cfac);
    RRRadialQkHydrogenicParams(NPARAMS, rqc, z, nq, nkl);
    for (ie = 0; ie < n_egrid; ie++) {
      xegrid[ie] = 1.0 + egrid[ie]/eb0;
      log_xegrid[ie] = log(xegrid[ie]);
    }

    for (ie = n_egrid-2; ie > 2; ie--) {
      a = log(tq[ie+1]/tq[ie]);
      b = xegrid[ie+1]/xegrid[ie];
      d = (sqrt(xegrid[ie]) + rqc[2])/(sqrt(xegrid[ie+1]) + rqc[2]);
      b = log(b);
      d = log(d);
      z = (a + (4.5+nkl)*b)/(0.5*b+d);
      if (a < 0 && z > 0) {
	rqc[1] = z;
	break;
      }
    }
    RRRadialQkFromFit(NPARAMS, rqc, n_egrid, xegrid, log_xegrid, 
		      rq, NULL, 0, &nkl);
    ie++;
    a = eb0*tq[ie]/rq[ie];
    rqc[0] *= a;
    rqc[3] = eb0;
    for (ie++; ie < n_egrid; ie++) {
      tq[ie] = a*(rq[ie]/eb0);
    }
    for (ie = 0; ie < n_egrid; ie++) {
      a = (*eb) + egrid[ie];
      rqu[ie] = tq[ie]*a;
    }
  } else {
    for (ie = 0; ie < n_egrid; ie++) { 
      a = *eb + egrid[ie];
      tq[ie] *= a;
      if (c) {
	a *= FINE_STRUCTURE_CONST;
	tq[ie] *= pow(a, c);
      }
    }
    if (qk_mode == QK_INTERPOLATE) {
      for (ie = 0; ie < n_egrid; ie++) {
	tq[ie] = log(tq[ie]);
      }
      UVIP3P(n_egrid, log_egrid, tq, n_usr, log_usr, rqu);
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = exp(rqu[ie]);
      }
    } else {
      for (ie = 0; ie < n_usr; ie++) {
	rqu[ie] = tq[ie];
      }
    }
  }      
 
  if (iuta) {
    d = lev1->uta_g*(q1/(j1+1.0));
    for (ie = 0; ie < n_usr; ie++) {
      rqu[ie] *= d;
    }
    rqc[0] *= d;

    free(idatum->bra);
    free(idatum);
  } else {
    free(ang);
  }

  return nkl;
}

int AutoionizeRate(cfac_t *cfac, double *rate, double *e, int rec, int f, int msub) {  
  LEVEL *lev1, *lev2;
  ANGULAR_ZxZMIX *ang;
  ANGULAR_ZFB *zfb;
  STATE *st;
  int k, nz, nzfb, ik, i, j1, j2, ij, kappaf, ip;
  int jf, k0, k1, kb, njf, nkappaf, klf, jmin, jmax;
  double *p, r, s, log_e, a;
  double *ai_pk, ai_pk0[MAXNE];
  int nt, m1, m2, m;
  int kappafp, jfp, klfp, dkl;

  *rate = 0.0;
  lev1 = GetLevel(cfac, rec);
  lev2 = GetLevel(cfac, f);

  if (GetLevNumElectrons(cfac, lev1) != GetLevNumElectrons(cfac, lev2) + 1) {
    return -1;
  }
  
  log_e = log(*e);

  i = lev1->pb;

  j1 = lev1->pj;
  j2 = lev2->pj;

  st = GetSymmetryState(GetSymmetry(cfac, j1), i);
  if (st->kgroup < 0) {
    k = GetOrbital(cfac, st->kcfg)->kappa;
  } else {
    k = (GetConfig(cfac, st)->shells[0]).kappa;
  }
  k = GetLFromKappa(k);
  k = k/2;
  if (k < pw_scratch.pw_limits[0] || k > pw_scratch.pw_limits[1]) return -1;

  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);

  jmin = abs(j1-j2);
  jmax = j1+j2;
  njf = (jmax - jmin)/2 + 1;
  nkappaf = njf*2;
  p = malloc(sizeof(double)*nkappaf);
  for (ip = 0; ip < nkappaf; ip++) p[ip] = 0.0;

  nz = AngularZxZFreeBound(cfac, &ang, f, rec);
  nt = 1;
  if (nz > 0) {
    for (i = 0; i < nz; i++) {
      jf = ang[i].k0;
      kb = ang[i].k1;
      k0 = ang[i].k2;
      k1 = ang[i].k3;
      kappafp = GetOrbital(cfac, kb)->kappa;
      klfp = GetLFromKappa(kappafp);
      kappafp = GetOrbital(cfac, k0)->kappa;
      klfp += GetLFromKappa(kappafp);
      kappafp = GetOrbital(cfac, k1)->kappa;
      klfp += GetLFromKappa(kappafp);      
      ij = (jf - jmin);
      for (ik = -1; ik <= 1; ik += 2) {
	klf = jf + ik;  	
	if (IsOdd((klfp+klf)/2)) continue;
	kappaf = GetKappaFromJL(jf, klf);
	AIRadialPk(cfac, &ai_pk, k0, k1, kb, kappaf, ang[i].k);
	if (n_egrid > 1) {
	  UVIP3P(n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	} else {
	  s = ai_pk[0];
	}
	ip = (ik == -1)? ij:(ij+1);
	p[ip] += s*ang[i].coeff;
      }
    }    
    free(ang);
  }
  
  nzfb = AngularZFreeBound(cfac, &zfb, f, rec);
  if (nzfb > 0) {
    for (i = 0; i < nzfb; i++) {
      kb = zfb[i].kb;
      kappaf = GetOrbital(cfac, kb)->kappa;
      GetJLFromKappa(kappaf, &jf, &klf);
      ij = jf - jmin;
      ik = klf - jf;
      AIRadial1E(cfac, ai_pk0, kb, kappaf);
      if (n_egrid > 1) {
	UVIP3P(n_egrid, log_egrid, ai_pk0, nt, &log_e, &s);
      } else {
	s = ai_pk0[0];
      }
      ip = (ik == -1)?ij:(ij+1);
      if (j2 > j1) {
	if (IsEven(ij/2)) s = -s;
      } else {
	if (IsOdd(ij/2)) s = -s;
      }
      p[ip] += s*zfb[i].coeff;
    }
    free(zfb);
  }
  if (nz <= 0 && nzfb <= 0) return -1;

  if (!msub) {
    r = 0.0;
    for (i = 0; i < nkappaf; i++) {
      r += p[i]*p[i];
    }
    /* the prefactor 4.0 includes the factor 2 from the continuum norm,
       otherwize, it should have been 2.0 */
    *rate = 4.0*r/(j1+1.0);
    free(p);    
    return 0;
  } else {
    k = 0;
    for (m1 = -j1; m1 <= 0; m1 += 2) {
      for (m2 = -j2; m2 <= j2; m2 += 2) {
	m = m1-m2;
	rate[k] = 0;
	rate[k+1] = 0;
	for (i = 0; i < nkappaf; i++) {
	  if (IsOdd(i)) {
	    ij = i-1;
	    klf = 1;
	  } else {
	    ij = i;
	    klf = -1;
	  }
	  jf = ij+jmin;
	  klf += jf;
	  kappaf = GetKappaFromJL(jf, klf);
	  s = W3j(j2, jf, j1, m2, m, -m1);
	  rate[k] += s*s*p[i]*p[i];
	  if (m != 1 && m != -1) continue;
	  for (ip = 0; ip < nkappaf; ip++) {
	    if (IsOdd(ip)) {
	      ij = ip-1;
	      klfp = 1;
	    } else {
	      ij = ip;
	      klfp = -1;
	    }
	    jfp = ij + jmin;
	    klfp += jfp;
	    if (ip == i) {
	      r = W3j(klf, 1, jf, 0, m, -m);
	      r = r*r*s*s*p[i]*p[i];
	      r *= (klf+1.0)*(jf+1.0);
	      rate[k+1] += r;
	    } else {
	      kappafp = GetKappaFromJL(jfp, klfp);
	      for (ik = 0; ik < n_egrid; ik++) {
		k0 = OrbitalIndex(cfac, 0, kappaf, egrid[ik]);
		k1 = OrbitalIndex(cfac, 0, kappafp, egrid[ik]);
		ai_pk0[ik] = GetPhaseShift(cfac, k0);
		ai_pk0[ik] -= GetPhaseShift(cfac, k1);
	      }	      
	      if (n_egrid > 1) {
		UVIP3P(n_egrid, log_egrid, ai_pk0, nt, &log_e, &a);
	      } else {
		a = ai_pk0[0];
	      }
	      r = W3j(klf, 1, jf, 0, m, -m);
	      r *= W3j(klfp, 1, jfp, 0, m, -m);
	      r *= W3j(j2, jfp, j1, m2, m, -m1);
	      r = r*s*p[i]*p[ip];
	      r *= sqrt((klf+1.0)*(klfp+1.0)*(jf+1.0)*(jfp+1.0));
	      r *= cos(a);
	      dkl = (jf + jfp)/2 + 1;
	      if (IsOdd(dkl)) r = -r;
	      rate[k+1] += r;
	    }
	  }
	}
	rate[k] *= 4.0;
	rate[k+1] *= 2.0*M_PI*M_PI/(*e);
	k += 2;
      }
    }
    free(p);    
    return k;
  }
}

int AutoionizeRateUTA(cfac_t *cfac, double *rate, double *e, int rec, int f) {
  INTERACT_DATUM *idatum;
  LEVEL *lev1, *lev2;
  int j0, j1, jb, ns, q0, q1, qb;
  int k0, k1, kb, kmin, kmax, jmin, jmax;
  int jf, ik, klf, kappaf, k, nt, j, jm;
  double a, b, r, s, log_e, *ai_pk;
  
  *rate = 0.0;
  lev1 = GetLevel(cfac, rec);
  lev2 = GetLevel(cfac, f);

  if (GetLevNumElectrons(cfac, lev1) != GetLevNumElectrons(cfac, lev2) + 1) {
    return -1;
  }
  
  log_e = log(*e);
  
  idatum = NULL;
  ns = GetInteract(cfac, &idatum, NULL, NULL, lev2->uta_cfg_g, lev1->uta_cfg_g,
		   lev2->uta_g_cfg, lev1->uta_g_cfg, 0, 0, 1);  
  if (ns <= 0) return -1;
  if (idatum->s[3].index < 0) {
    free(idatum->bra);
    free(idatum);
    return -1;
  }

  kb = OrbitalIndex(cfac, idatum->s[1].n, idatum->s[1].kappa, 0.0);
  k0 = OrbitalIndex(cfac, idatum->s[2].n, idatum->s[2].kappa, 0.0);
  k1 = OrbitalIndex(cfac, idatum->s[3].n, idatum->s[3].kappa, 0.0);
  j0 = idatum->s[2].j;
  j1 = idatum->s[3].j;
  jb = idatum->s[1].j;
  q0 = idatum->s[2].nq_ket;
  q1 = idatum->s[3].nq_ket;
  qb = idatum->s[1].nq_ket;

  if (idatum->s[1].index != idatum->s[3].index) {
    kmin = abs(j0-j1);
    kmax = j0 + j1;
    jmin = 1;
    jmax = j1+j0+jb;
    nt = 1;
    r = 0.0;
    for (jf = jmin; jf <= jmax; jf += 2) {
      for (ik = -1; ik <= 1; ik += 2) {
	klf = jf + ik;
	kappaf = GetKappaFromJL(jf, klf);
	for (k = kmin; k <= kmax; k += 2) {
	  if (!Triangle(j0, j1, k) || !Triangle(jb, jf, k)) continue;
	  AIRadialPk(cfac, &ai_pk, k0, k1, kb, kappaf, k);
	  if (n_egrid > 1) {
	    UVIP3P(n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	  } else {
	    s = ai_pk[0];
	  }
	  s = s*s/(k + 1.0);
	  r += s;
	}
      }
    }  
    r *= 4.0*(q1/(j1+1.0))*(qb/(jb+1.0))*((j0+1.0-q0)/(j0+1.0));
  } else {
    jm = 2*j1;
    r = 0.0;
    nt = 1;
    for (j = 0; j <= jm; j += 4) {
      jmin = abs(j-j0);
      jmax = j+j0;
      for (jf = jmin; jf <= jmax; jf += 2) {
	for (ik = -1; ik <= 1; ik += 2) {
	  klf = jf + ik;
	  kappaf = GetKappaFromJL(jf, klf);
	  kmin = abs(j0-j1);
	  kmax = j0 + j1;
	  a = 0.0;
	  for (k = kmin; k <= kmax; k += 2) {
	    if (!Triangle(jb, jf, k)) continue;
	    b = W6j(j, jf, j0, k, j1, j1);
	    if (fabs(b) < EPS30) continue;
	    AIRadialPk(cfac, &ai_pk, k0, k1, kb, kappaf, k);
	    if (n_egrid > 1) {
	      UVIP3P(n_egrid, log_egrid, ai_pk, nt, &log_e, &s);
	    } else {
	      s = ai_pk[0];
	    } 
	    a += b*s;
	  }
	  r += a*a*2.0*(j+1.0);
	}
      }
    }
    r *= 4.0*(q1/(j1+1.0))*((qb-1.0)/jb)*((j0+1.0-q0)/(j0+1.0));
  }

  *rate = r;
  
  free(idatum->bra);
  free(idatum);
  return 0;
}

int AIRadial1E(cfac_t *cfac, double *ai_pk, int kb, int kappaf) {
  int kf;
  int i;

  for (i = 0; i < n_egrid; i++) {
    kf = OrbitalIndex(cfac, 0, kappaf, egrid[i]);
    ResidualPotential(cfac, ai_pk+i, kf, kb);
  }
  return 0;
}  

int AIRadialPk(cfac_t *cfac, double **ai_pk, int k0, int k1, int kb, int kappaf, int k) {
  int i, kf;
  int ks[4];
  double e, sd, se;
  double **p;
  int index[5];

  if (kappaf > 0) {
    index[0] = 2*kappaf-1;
  } else {
    index[0] = -2*kappaf-2;
  }
  index[1] = kb;
  index[2] = k0;
  index[3] = k1;
  index[4] = k/2;

  p = (double **) MultiSet(pk_array, index, NULL);
  if (*p) {
    *ai_pk = *p;
    return 0;
  } 
 
  (*p) = (double *) malloc(sizeof(double)*n_egrid);
  *ai_pk = *p;
  for (i = 0; i < n_egrid; i++) {
    e = egrid[i];
    kf = OrbitalIndex(cfac, 0, kappaf, e);
    ks[0] = k0;
    ks[1] = kf;
    ks[2] = k1;
    ks[3] = kb;
    if (SlaterTotal(cfac, &sd, &se, NULL, ks, k, 0) != 0) {
      return -1;
    }
    (*ai_pk)[i] = sd+se;
  }

  return 0;
}

int PrepRREGrids(double e, double emax0) { 
  double rmin, rmax;
  double emin, emax;
  int j;

  if (egrid_limits_type == 0) {
    rmin = egrid_min;
    rmax = egrid_max;
  } else {
    rmin = egrid_min/e;
    rmax = egrid_max/e;
  }
  emin = rmin*e;
  emax = rmax*e;
  if (emax < emax0) {
    emax = 50.0*e;
    if (emax > emax0) emax = emax0;
  }
  egrid_type = 1;
  if (usr_egrid_type < 0) usr_egrid_type = 1;

  if (n_egrid == 0) {
    n_egrid = 6;
  }
  if (egrid[0] < 0.0) {
    SetPEGrid(n_egrid, emin, emax, e);
  }
  if (n_usr <= 0) {
    SetUsrPEGridDetail(n_egrid, egrid);
    usr_egrid_type = 1;
  } else if (usr_egrid[0] < 0) {
    SetUsrPEGrid(n_usr, emin, emax, e);
    usr_egrid_type = 1;
  }

  if (qk_mode == QK_INTERPOLATE) {
    for (j = 0; j < n_egrid; j++) {
      log_egrid[j] = egrid[j];
      if (egrid_type == 1) log_egrid[j] += e;
      log_egrid[j] = log(log_egrid[j]);
    }
    for (j = 0; j < n_usr; j++) {
      log_usr[j] = usr_egrid[j];
      if (usr_egrid_type == 1) log_usr[j] += e;
      log_usr[j] = log(log_usr[j]);
    }
  }
  return 0;
}

int SaveRRMultipole(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn, int m) {
  int i, j, k;
  FILE *f;
  LEVEL *lev1, *lev2;
  double e, emin, emax, emax0;
  double awmin, awmax;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set;
  double c, e0, e1;

  f = fopen(fn, "w");
  if (f == NULL) return -1;
  fprintf(f, "# This file contains 11 columns\n");
  fprintf(f, "#  1, Bound state index\n");
  fprintf(f, "#  2, 2*J_Bound\n");
  fprintf(f, "#  3, Ionized state index\n");
  fprintf(f, "#  4, 2*J_Ionized\n");
  fprintf(f, "#  5, Multipole type, -1=E1, 1=M1, -2=E2, 2=M2, ...\n");
  fprintf(f, "#  6, 2*J_total, coupled angular momentum of J_Ionized and J\n");
  fprintf(f, "#  7, 2*L, L is the orbital angular momentum of the photo-electron\n");
  fprintf(f, "#  8, 2*J, J is the total angular momentum of the photo-electron\n");
  fprintf(f, "#  9, E_th, Ionization threshold in Hartree\n");
  fprintf(f, "# 10, E_e, photo-electron energy in Hartree\n");
  fprintf(f, "# 11, Multipole matrix element in atomic unit\n");
  fprintf(f, "\n\n");
  
  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(cfac, low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(cfac, up[j]);
      e = lev2->energy - lev1->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (k == 0) {
    return 0;
  }
  
  if (tegrid[0] < 0) {
    te_set = 0;
  } else {
    te_set = 1;
  }
  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }
  if (usr_egrid[0] < 0) {
    usr_set = 0;
  } else {
    usr_set = 1;
  }
  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;
  n_usr0 = n_usr;

  if (egrid_limits_type == 0) {
    emax0 = 0.5*(emin + emax)*egrid_max;
  } else {
    emax0 = egrid_max;
  }
  ArrayInit(&subte, sizeof(double), 128, NULL, NULL);
  ArrayAppend(&subte, &emin);
  c = TE_MAX_MIN;
  if (!e_set || !te_set) {
    e = c*emin;
    while (e < emax) {
      ArrayAppend(&subte, &e);
      e *= c;
    }
  }
  ArrayAppend(&subte, &emax);

    
  e0 = emin;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(cfac, low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(cfac, up[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }
    emin = e0;
    emax = e1;  
    if (m == 1 || GetTransitionMode(cfac) == M_FR) {
      e = (emax - emin)/(0.5*(emin+emax));
      if (!te_set) {
	if (e < 0.1) {
	  SetRRTEGrid(1, emin, emax);
	} else if (e < 0.5) {
	  SetRRTEGrid(2, emin, emax);
	} else {
	  if (k == 2) n_tegrid = 2;
	  else if (n_tegrid0 == 0) n_tegrid = 3;
	  SetRRTEGrid(n_tegrid, emin, emax);
	}
      }
      FreeMultipoleArray(cfac);
      awmin = emin * FINE_STRUCTURE_CONST;
      awmax = emax * FINE_STRUCTURE_CONST;
      if (e < 0.3) {
	SetAWGrid(cfac, 1, awmin, awmax);
      } else if (e < 1.0) {
	SetAWGrid(cfac, 2, awmin, awmax);
      } else {
	SetAWGrid(cfac, 3, awmin, awmax);
      }
    } else {
      SetRRTEGrid(1, emin, emax);
    }
    
    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    e = 0.5*(emin + emax);
    PrepRREGrids(e, emax0);
    
    for (i = 0; i < nup; i++) {
      lev1 = GetLevel(cfac, up[i]);
      for (j = 0; j < nlow; j++) {
	lev2 = GetLevel(cfac, low[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	BoundFreeMultipole(cfac, f, low[j], up[i], m);
      }
    }
    
    ReinitRadial(cfac, 1);
    FreeRecQk();
    FreeRecPk();
    
    e0 = e1;
  }
  
  fclose(f);
  ArrayFree(&subte);
  ReinitRecombination(1);

  return 0;
}
    
int SaveRecRR(cfac_t *cfac, int nlow, int *low, int nup, int *up, 
	      char *fn, int m) {
  int i, j, k, ie, ip;
  FILE *f;
  double rqu[MAXNUSR], qc[NPARAMS+1];
  double eb;
  LEVEL *lev1, *lev2;
  RR_RECORD r;
  RR_HEADER rr_hdr;
  F_HEADER fhdr;
  double e, emin, emax, emax0;
  double awmin, awmax;
  int nq, nqk;
  ARRAY subte;
  int isub, n_tegrid0, n_egrid0, n_usr0;
  int te_set, e_set, usr_set;
  double c, e0, e1;

  if (m != -1 && GetTransitionGauge(cfac) != G_BABUSHKIN && qk_mode == QK_FIT) {
    printf("QK_FIT mode is only available to LENGTH form of E1 transitions\n");
    printf("Changing QK_FIT to QK_INTERPOLATE.\n");
    SetRecQkMode(QK_INTERPOLATE, -1.0);
  }

  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(cfac, low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(cfac, up[j]);
      e = lev2->energy - lev1->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (k == 0) {
    return 0;
  }
  
  if (tegrid[0] < 0) {
    te_set = 0;
  } else {
    te_set = 1;
  }
  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }
  if (usr_egrid[0] < 0) {
    usr_set = 0;
  } else {
    usr_set = 1;
  }
  n_tegrid0 = n_tegrid;
  n_egrid0 = n_egrid;
  n_usr0 = n_usr;

  if (egrid_limits_type == 0) {
    emax0 = 0.5*(emin + emax)*egrid_max;
  } else {
    emax0 = egrid_max;
  }
  ArrayInit(&subte, sizeof(double), 128, NULL, NULL);
  ArrayAppend(&subte, &emin);
  c = TE_MAX_MIN;
  if (!e_set || !te_set) {
    e = c*emin;
    while (e < emax) {
      ArrayAppend(&subte, &e);
      e *= c;
    }
  }
  ArrayAppend(&subte, &emax);

  if (qk_mode == QK_FIT) {
    nqk = NPARAMS+1;
    r.params = (float *) malloc(sizeof(float)*nqk);
  } else {
    nqk = 0;
  }

  fhdr.type = DB_RR;
  strcpy(fhdr.symbol, cfac_get_atomic_symbol(cfac));
  fhdr.atom = cfac_get_atomic_number(cfac);
  rr_hdr.nele = GetNumElectrons(cfac, low[0]);
  rr_hdr.qk_mode = qk_mode;
  rr_hdr.nparams = nqk;
  rr_hdr.multipole = m;
  f = OpenFile(fn, &fhdr);
  if (!f) {
    return -1;
  }
  
  e0 = emin*0.999;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(cfac, low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(cfac, up[j]);
	e = lev2->energy - lev1->energy;
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }
    emin = e0;
    emax = e1;  
    if (m == 1 || GetTransitionMode(cfac) == M_FR) {
      e = (emax - emin)/(0.5*(emin+emax));
      if (!te_set) {
	if (e < EPS3) {
	  SetRRTEGrid(1, emin, emax);
	} else if (e < 0.5) {
	  SetRRTEGrid(2, emin, emax);
	} else {
	  if (k == 2) n_tegrid = 2;
	  else if (n_tegrid0 == 0) n_tegrid = 3;
	  SetRRTEGrid(n_tegrid, emin, emax);
	}
      }
      FreeMultipoleArray(cfac);
      awmin = emin * FINE_STRUCTURE_CONST;
      awmax = emax * FINE_STRUCTURE_CONST;
      if (e < 0.3) {
	SetAWGrid(cfac, 1, awmin, awmax);
      } else if (e < 1.0) {
	SetAWGrid(cfac, 2, awmin, awmax);
      } else {
	SetAWGrid(cfac, 3, awmin, awmax);
      }
    } else {
      SetRRTEGrid(1, emin, emax);
    }
    
    n_egrid = n_egrid0;
    n_usr = n_usr0;
    if (!usr_set) usr_egrid[0] = -1.0;
    if (!e_set) egrid[0] = -1.0;
    e = 0.5*(emin + emax);
    PrepRREGrids(e, emax0);
    
    if (qk_mode == QK_FIT && n_egrid <= NPARAMS) {
      printf("n_egrid must > %d to use QK_FIT mode\n", NPARAMS);
      return -1;
    }
    rr_hdr.n_tegrid = n_tegrid;
    rr_hdr.tegrid = tegrid;
    rr_hdr.n_egrid = n_egrid;
    rr_hdr.egrid = egrid;
    rr_hdr.n_usr = n_usr;
    rr_hdr.usr_egrid = usr_egrid;
    rr_hdr.egrid_type = egrid_type;
    rr_hdr.usr_egrid_type = usr_egrid_type;
    
    r.strength = (float *) malloc(sizeof(float)*n_usr);
    
    InitFile(f, &fhdr, &rr_hdr);
    
    for (i = 0; i < nup; i++) {
      lev1 = GetLevel(cfac, up[i]);
      for (j = 0; j < nlow; j++) {
	int iuta;
	
	lev2 = GetLevel(cfac, low[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
        
        if (lev1->uta || lev2->uta) {
	  iuta = 1;
	} else {
	  iuta = 0;
	}
	nq = BoundFreeOS(cfac, rqu, qc, &eb, low[j], up[i], m, iuta);
	if (nq < 0) continue;
	
        r.b = low[j];
	r.f = up[i];
	r.kl = nq;
	
	if (qk_mode == QK_FIT) {
	  for (ip = 0; ip < nqk; ip++) {
	    r.params[ip] = (float) qc[ip];
	  }
	}
	
	for (ie = 0; ie < n_usr; ie++) {
	  r.strength[ie] = (float) rqu[ie];
	}
	WriteRRRecord(f, &r);
      }
    }      

    DeinitFile(f, &fhdr);
    
    free(r.strength);
    ReinitRadial(cfac, 1);
    FreeRecQk();
    FreeRecPk();
    
    e0 = e1;
  }

  if (qk_mode == QK_FIT) {
    free(r.params);
  }
      
  ReinitRecombination(1);

  ArrayFree(&subte);
  CloseFile(f, &fhdr);

  return 0;
}
      
int SaveAI(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn, 
	   int msub) {
  int i, j, k, t;
  LEVEL *lev1, *lev2;
  AI_RECORD r;
  AIM_RECORD r1;
  AI_HEADER ai_hdr;
  AIM_HEADER ai_hdr1;
  F_HEADER fhdr;
  double emin, emax;
  double e, s, s1[MAXAIM];
  float rt[MAXAIM];
  FILE *f;
  ARRAY subte;
  double c, e0, e1, b;
  int isub, n_egrid0;
  int e_set;

  int nc;

  if (nup <= 0 || nlow <= 0) return -1;

  nc = OverlapLowUp(nlow, low, nup, up);

  emin = 1E10;
  emax = 1E-10;
  k = 0;
  for (i = 0; i < nlow; i++) {
    lev1 = GetLevel(cfac, low[i]);
    for (j = 0; j < nup; j++) {
      lev2 = GetLevel(cfac, up[j]);
      
      if (GetLevNumElectrons(cfac, lev1) != GetLevNumElectrons(cfac, lev2) + 1) {
        continue;
      }
      
      e = lev1->energy - lev2->energy;
      if (e > 0) k++;
      if (e < emin && e > 0) emin = e;
      if (e > emax) emax = e;
    }
  }
  if (k == 0) {
    return 0;
  }

  if (egrid[0] < 0) {
    e_set = 0;
  } else {
    e_set = 1;
  }

  n_egrid0 = n_egrid;

  ArrayInit(&subte, sizeof(double), 128, NULL, NULL);
  ArrayAppend(&subte, &emin);
  c = TE_MAX_MIN;
  if (!e_set) {
    b = c*emin;
    while (b < emax) {
      ArrayAppend(&subte, &b);
      b *= c;
    }
  }
  ArrayAppend(&subte, &emax);
  
  if (!msub) {
    fhdr.type = DB_AI;
  } else {
    fhdr.type = DB_AIM;
  }
  strcpy(fhdr.symbol, cfac_get_atomic_symbol(cfac));
  fhdr.atom = cfac_get_atomic_number(cfac);
  if (!msub) {
    ai_hdr.nele = GetNumElectrons(cfac, low[0]);
    ai_hdr.emin = 0.0;
  } else {
    ai_hdr1.nele = GetNumElectrons(cfac, low[0]);
    ai_hdr1.emin = 0.0;
  }
  f = OpenFile(fn, &fhdr);
  if (!f) {
    return -1;
  }

  e0 = emin*0.999;
  for (isub = 1; isub < subte.dim; isub++) {
    e1 = *((double *) ArrayGet(&subte, isub));
    if (isub == subte.dim-1) e1 = e1*1.001;
    emin = e1;
    emax = e0;
    k = 0;
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(cfac, low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(cfac, up[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	if (e < emin) emin = e;
	if (e > emax) emax = e;
	k++;
      }
    }
    if (k == 0) {
      e0 = e1;
      continue;
    }
    if (!e_set) {
      double a = (emax-emin)/(0.5*(emax+emin));
      if (a < EPS3) {
	a = 0.5*(emin+emax);
	SetPEGrid(1, a, a, 0.0);
      } else if (a < 0.4) {
	SetPEGrid(2, emin, emax, 0.0);
      } else if (a < 1.0) {
	if (k == 2) n_egrid = 2;
	else if (n_egrid0 == 0)	n_egrid = 3;
	SetPEGrid(n_egrid, emin, emax, 0.0);
      } else {
	if (k == 2) n_egrid = 2;
	else if (n_egrid0 == 0) n_egrid = 4;
	SetPEGrid(n_egrid, emin, emax, 0.0);
      }
    }      
    
    if (!msub) {
      ai_hdr.n_egrid = n_egrid;
      ai_hdr.egrid = egrid;
      InitFile(f, &fhdr, &ai_hdr);
    } else {
      ai_hdr1.n_egrid = n_egrid;
      ai_hdr1.egrid = egrid;
      InitFile(f, &fhdr, &ai_hdr1);
    }
    
    for (i = 0; i < nlow; i++) {
      lev1 = GetLevel(cfac, low[i]);
      for (j = 0; j < nup; j++) {
	lev2 = GetLevel(cfac, up[j]);
	e = lev1->energy - lev2->energy;
	if (e < e0 || e >= e1) continue;
	if (!msub) {
	  /* FIXME: generalize AutoionizeRateUTA to treat detailed mode */
          if (lev1->uta || lev2->uta) {
	    k = AutoionizeRateUTA(cfac, &s, &e, low[i], up[j]);
          } else {
            k = AutoionizeRate(cfac, &s, &e, low[i], up[j], msub);
          }
	  if (k < 0) continue;
	  if (s < ai_cut) continue;
	  r.b = low[i];
	  r.f = up[j];
	  r.rate = s;
	  WriteAIRecord(f, &r);
	} else {
	  k = AutoionizeRate(cfac, s1, &e, low[i], up[j], msub);
	  if (k < 0) continue;
	  r1.rate = rt;
	  s = 0;
	  for (t = 0; t < k; t++) {
	    r1.rate[t] = s1[t];
	    s += s1[t];
	  }
	  if (s < ai_cut) continue;
	  r1.b = low[i];
	  r1.f = up[j];
	  r1.nsub = k;
	  WriteAIMRecord(f, &r1);
	}
      }
    }

    DeinitFile(f, &fhdr);

    ReinitRadial(cfac, 1);
    FreeRecQk();
    FreeRecPk();

    e0 = e1;
  }

  ReinitRecombination(1);
  
  ArrayFree(&subte);
  CloseFile(f, &fhdr);

  return 0;
}

int FreeRecPk(void) {
  MultiFreeData(pk_array);
  return 0;
}

int FreeRecQk(void) {
  if (qk_array->array == NULL) return 0;
  MultiFreeData(qk_array);
  return 0;
}

int InitRecombination(void) {
  int blocks[5];
  int ndim;
  int i;
  
  ndim = 5;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK5;
  pk_array = (MULTI *) malloc(sizeof(MULTI));
  MultiInit(pk_array, sizeof(double *), ndim, blocks, 
			   FreeRecPkData, InitPointerData);
  
  ndim = 3;
  for (i = 0; i < ndim; i++) blocks[i] = MULTI_BLOCK3;
  qk_array = (MULTI *) malloc(sizeof(MULTI));
  blocks[0] = 10;
  blocks[1] = 10;
  blocks[2] = 4;
  MultiInit(qk_array, sizeof(double *), ndim, blocks,
    FreeRecPkData, InitPointerData);  
  
  hyd_qk_array = (ARRAY *) malloc(sizeof(ARRAY));
  ArrayInit(hyd_qk_array, sizeof(double *), 32, NULL, InitPointerData);
  
  n_complex = 0;
  for (i = 0; i < MAX_COMPLEX; i++) {
    rec_complex[i].rg = (ARRAY *) malloc(sizeof(ARRAY));
    ArrayInit(rec_complex[i].rg, sizeof(int), 64, NULL, NULL);
  }
  n_egrid = 0;
  egrid[0] = -1.0;
  SetPEGridLimits(-1, -1, 0);
  n_usr = 0;
  usr_egrid[0] = -1.0;
  
  n_tegrid = 0;
  tegrid[0] = -1.0;

  SetRecQkMode(QK_DEFAULT, 0.1);
  SetRecPWOptions(RECLMAX, RECLMAX);
  return 0;
}

int ReinitRecombination(int m) {
  int i;

  if (m < 0) return 0;

  FreeRecQk();
  FreeRecPk();

  n_egrid = 0;
  egrid[0] = -1.0;
  SetPEGridLimits(-1, -1, 0);
  n_usr = 0;
  usr_egrid[0] = -1.0;
  
  n_tegrid = 0;
  tegrid[0] = -1.0;

  SetRecQkMode(QK_DEFAULT, 0.1);
  SetRecPWOptions(RECLMAX, RECLMAX);

  if (m > 0) return 0;

  for (i = 0; i < n_complex; i++) {
    ArrayFree(rec_complex[i].rg);
  }
  n_complex = 0;
  
  return 0;
}
