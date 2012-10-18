#include <string.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include "cfacP.h"
#include "radial.h"
#include "angular.h"
#include "dbase.h"
#include "structure.h"

static ARRAY *cfac_get_ion_levels(cfac_t *cfac, unsigned int nele)
{
    if (nele > cfac->anum) {
        return NULL;
    } else {
        return &cfac->levels_per_ion[nele];
    }
}

static int cfac_get_ion_nlevels(cfac_t *cfac, unsigned int nele)
{
    ARRAY *larray = cfac_get_ion_levels(cfac, nele);
    if (larray) {
        return larray->dim;
    } else {
        return 0;
    }
}

void GetFields(const cfac_t *cfac, double *b, double *e, double *a) {
  *b = cfac->bf;
  *e = cfac->ef;
  *a = cfac->eb_angle;
}

/* 
** if the angle a > 0, the B is Z-axis, BxE is Y-axis. and E is in X-Z plane.
** if the angle a < 0, the E is Z-axis, ExB is Y-axis, and B is in X-Z plane.
** angle a is always measured from E->B, 
*/
void SetFields(cfac_t *cfac, double b, double e, double a, int m) {
  int i, q, i1, i2, q1, q2;
  double w, mass;

  cfac->ef = e;
  cfac->bf = b;
  cfac->eb_angle = a;
  
  a = fabs(a);
  a *= M_PI/180.0;
  
  b *= MBOHR/HARTREE_EV;
  e *= RBOHR/HARTREE_EV;
  
  if (cfac->eb_angle >= 0) {
    cfac->b1[0] = cfac->b1[2] = 0.0;
    cfac->b1[1] = b;
    cfac->e1[0] = e*sin(a)/sqrt(2);
    cfac->e1[1] = e*cos(a);
    cfac->e1[2] = -cfac->e1[0];
  } else {
    cfac->e1[0] = cfac->e1[2] = 0.0;
    cfac->e1[1] = e;
    cfac->b1[0] = -b*sin(a)/sqrt(2);
    cfac->b1[1] = b*cos(a);
    cfac->b1[2] = -cfac->b1[0];
  }

  mass = GetAtomicMass(cfac);
  mass = 1.0 + 5.45683e-4/mass;
  for (i = 0; i < 3; i++) {
    cfac->b1[i] *= mass;
  }

  cfac->b0 = 0.0;
  for (i = 0; i < 5; i++) {
    cfac->b2[i] = 0.0;
  }

  if (m == 0) {
    for (i1 = 0; i1 < 3; i1++) {
      q1 = 2*(i1-1);
      if (cfac->b1[i1] == 0) continue;
      for (i2 = 0; i2 < 3; i2++) {
	q2 = 2*(i2-1);
	if (cfac->b1[i2] == 0) continue;
	w = W3j(2, 2, 0, q1, q2, 0);
	if (w) {
	  cfac->b0 += w*cfac->b1[i1]*cfac->b1[i2];
	}
	for (i = 0; i < 5; i++) {
	  q = 2*(i-2);
	  w = W3j(2, 2, 4, q1, q2, q);
	  if (w) {
	    cfac->b2[i] += w*cfac->b1[i1]*cfac->b1[i2];
	  }
	}
      }
    }
    cfac->b0 *= sqrt(3)*W6j(2, 2, 0, 2, 2, 2);
    for (i = 0; i < 5; i++) {
      cfac->b2[i] *= -sqrt(30)*W6j(2, 2, 4, 2, 2, 2);
    }
    cfac->b0 *= mass;
    for (i = 0; i < 5; i++) {
      cfac->b2[i] *= mass;
    }
  }

  /*
  printf("EB: %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E\n", 
	 e1[0], e1[1], e1[2], b1[0], b1[1], b1[2]);
  */
}

void InitLevelData(void *p, int n) {
  LEVEL *lev;
  int k;

  lev = (LEVEL *) p;
  for (k = 0; k < n; k++, lev++) {
    lev->n_basis = 0;
  }
}

int SetCILevel(cfac_t *cfac, int m) {
  cfac->confint = m;
  return 0;
}

int SetAngZCut(cfac_t *cfac, double cut) {
  if (cut >= 0) cfac->angz_cut = cut;
  else cfac->angz_cut = ANGZCUT;
  return 0;
}

int SetMixCut(cfac_t *cfac, double cut, double cut2) {
  if (cut >= 0) cfac->mix_cut = cut;
  else cfac->mix_cut = MIXCUT;
  if (cut2 >= 0) cfac->mix_cut2 = cut2;
  else cfac->mix_cut2 = MIXCUT2;
  return 0;
}

int SetAngZOptions(cfac_t *cfac, int n, double mix, double cut) {
  cfac->angz_maxn = n;
  cfac->mix_cut = mix;
  cfac->angz_cut = cut;
  return 0;
}

static void FlagClosed(cfac_t *cfac, SHAMILTON *hs) {
  int i, j, k, m1, m2;
  CONFIG *c;
  unsigned char t[MBCLOSE];

  for (i = 0; i < hs->nbasis; i++) {
    c = GetConfig(cfac, hs->basis[i]);
    for (k = 0; k < MBCLOSE; k++) {
      t[k] = 0;
    }
    for (j = 0; j < c->n_shells; j++) {
      k = ShellToInt(c->shells[j].n, c->shells[j].kappa);
      m1 = k/8;      
      if (m1 >= MBCLOSE) continue;
      m2 = k%8;
      if (ShellClosed(c->shells+j)) {
	t[m1] |= (1 << m2);
      }      
    }
    for (k = 0; k < MBCLOSE; k++) {
      if (i == 0) hs->closed[k] = t[k];
      else hs->closed[k] &= t[k];
    }
  }
}

static int IsClosedShell(const cfac_t *cfac, int ih, int k) {
  int i, j;
  
  if (ih >= MAX_HAMS) {
    abort();
  }

  i = k/8;
  if (i >= MBCLOSE) return 0;
  j = k%8;
  return (cfac->hams[ih].closed[i] & (1 << j));
}

static int ConstructHamiltonDiagonal(cfac_t *cfac,
  int isym, int k, int *kg, int m) {
  HAMILTON *h = cfac->hamiltonian;
  int i, j, t;
  SHAMILTON *hs;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r;

  DecodePJ(isym, &i, &j);

  if (k <= 0) return -1;
  sym = GetSymmetry(cfac, isym);
  if (sym == NULL) return -1;
  st = &(sym->states);
  j = 0;
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) j++;
  }
  if (j == 0) return -1;

  h->pj = isym;

  h->dim = j;
  h->n_basis = j;
  h->hsize = j;

  if (h->basis == NULL) {
    h->n_basis0 = h->n_basis;
    h->basis = (int *) malloc(sizeof(int)*(h->n_basis));
  } else if (h->n_basis > h->n_basis0) {
    h->n_basis0 = h->n_basis;
    free(h->basis);
    h->basis = (int *) malloc(sizeof(int)*h->n_basis);
  }
  if (!(h->basis)) goto ERROR;

  if (h->hamilton == NULL) {
    h->hsize0 = h->hsize;
    h->hamilton = (double *) malloc(sizeof(double)*h->hsize);
  } else if (h->hsize > h->hsize0) {
    h->hsize0 = h->hsize;
    free(h->hamilton);
    h->hamilton = (double *) malloc(sizeof(double)*h->hsize);
  }
  if (!(h->hamilton)) goto ERROR;

  j = 0;  
  for (t = 0; t < sym->n_states; t++) {
    s = (STATE *) ArrayGet(st, t);
    if (InGroups(s->kgroup, k, kg)) {
      h->basis[j] = t;
      j++;
    }
  }

  for (j = 0; j < h->dim; j++) {
    s = ArrayGet(st, h->basis[j]);
    if (m == 0) {
      r = ZerothEnergyConfig(cfac, GetConfig(cfac, s));
    } else {
      r = HamiltonElement(cfac, isym, h->basis[j], h->basis[j]);
    }
    h->hamilton[j] = r;
  }

  if (m > 0) {
    if (cfac->nhams >= MAX_HAMS) {
      printf("Number of hamiltons exceeded the maximum %d\n", MAX_HAMS);
      exit(1);
    }
    hs = &cfac->hams[cfac->nhams];
    cfac->nhams++;
    hs->pj = h->pj;
    hs->nlevs = h->dim;
    hs->nbasis = h->n_basis;
    hs->basis = malloc(sizeof(STATE *)*hs->nbasis);
    for (t = 0; t < h->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), h->basis[t]);
      hs->basis[t] = s;
    }
    FlagClosed(cfac, hs);
  }

  return 0;

 ERROR:
  printf("ConstructHamiltonDiagonal Error\n");
  return -1;
}

int CodeBasisEB(int s, int m) {
  int k;
  
  k = s + MAXLEVEB*abs(m);
  if (m < 0) k = -k;

  return k;
}

void DecodeBasisEB(int k, int *s, int *m) {
  *m = abs(k)/MAXLEVEB;
  *s = abs(k)%MAXLEVEB;
  if (k < 0) *m = -(*m);
}

int ConstructHamiltonEB(cfac_t *cfac, int n, int *ilev) {
  HAMILTON *h = cfac->hamiltonian;
  int i, j, p, k, t, m;
  double r;
  LEVEL *lev;
  int n_basis;

  h->pj = -1;
  ClearAngularFrozen(cfac);
  AngularFrozen(cfac, n, ilev, 0, NULL);

  n_basis = 0;
  for (i = 0; i < n; i++) {
    lev = GetLevel(cfac, ilev[i]);    
    DecodePJ(lev->pj, &p, &j);
    n_basis += j+1;
  }

  if (AllocHamMem(h, n_basis, n_basis) == -1) goto ERROR;
  k = 0;
  for (i = 0; i < n; i++) {
    lev = GetLevel(cfac, ilev[i]);
    DecodePJ(lev->pj, &p, &j);
    for (m = -j; m <= j; m += 2) {
      h->basis[k] = CodeBasisEB(ilev[i], m);
      k++;
    }
  }

  for (j = 0; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      r = HamiltonElementEB(cfac, i, j);
      /*
      printf("HAM: %8d %8d %15.8E\n", h->basis[i], h->basis[j], r);
      */
      h->hamilton[i+t] = r;
    }
  }

  return 0;

 ERROR:
  printf("ConstructHamiltonEB Error\n");
  return -1;
}

int ConstructHamilton(cfac_t *cfac,
    int isym, int k0, int k, int *kg, int kp, int *kgp, int md) {
  HAMILTON *h = cfac->hamiltonian;
  int i, j, j0, t, jp = 0, m1, m2, m3;
  SHAMILTON *hs;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r;

  /* 
  ** the return code -2, and -3, here distinguses it from the case where 
  ** no basis exists later.
  **/
  DecodePJ(isym, &i, &j);

  sym = GetSymmetry(cfac, isym);
  if (sym == NULL) return -1;
  h->pj = isym;

  m1 = md/100;
  t = md%100;
  m2 = t/10;
  m3 = t%10;
  
  if (m1) {
    if (k <= 0) return -1;
    if (cfac->confint == -1) {
      return ConstructHamiltonDiagonal(cfac, isym, k, kg, 1);
    }
    st = &(sym->states);
    j = 0;
    j0 = 0;
    for (t = 0; t < sym->n_states; t++) {
      s = (STATE *) ArrayGet(st, t);
      if (InGroups(s->kgroup, k0, kg)) j0++;
      if (InGroups(s->kgroup, k, kg)) j++;
    }
    if (j0 == 0) return -1;

    jp = 0;
    if (kp > 0) {
      for (t = 0; t < sym->n_states; t++) {
	s = (STATE *) ArrayGet(st, t);
	if (InGroups(s->kgroup, kp, kgp)) jp++;
      }
    }    

    if (AllocHamMem(h, j, jp+j) == -1) goto ERROR;
    
    j = 0;  
    for (t = 0; t < sym->n_states; t++) {
      s = (STATE *) ArrayGet(st, t);
      if (InGroups(s->kgroup, k, kg)) {
	h->basis[j] = t;
	j++;
      }
    }
    if (jp > 0) {  
      for (t = 0; t < sym->n_states; t++) {
	s = (STATE *) ArrayGet(st, t);
	if (kp > 0 && InGroups(s->kgroup, kp, kgp)) {
	  h->basis[j] = t;
	  j++;
	}
      }
    }
  }
  if (m2) {
    for (j = 0; j < h->dim; j++) {
      t = j*(j+1)/2;
      for (i = 0; i <= j; i++) {
	r = HamiltonElement(cfac, isym, h->basis[i], h->basis[j]);
	h->hamilton[i+t] = r;
      }
    } 
    
    if (jp > 0) {
      t = ((h->dim+1)*(h->dim))/2;
      for (i = 0; i < h->dim; i++) {
	for (j = h->dim; j < h->n_basis; j++) {
	  r = HamiltonElement(cfac, isym, h->basis[i], h->basis[j]);
	  h->hamilton[t++] = r;
	}
	ReinitRecouple(0);
	ReinitRadial(cfac, 1);
      }
      for (j = h->dim; j < h->n_basis; j++) {
	r = HamiltonElement(cfac, isym, h->basis[j], h->basis[j]);
	h->hamilton[t++] = r;
      }
      ReinitRecouple(0);
      ReinitRadial(cfac, 1);
    }
  }
  if (m3) {
    if (cfac->nhams >= MAX_HAMS) {
      printf("Number of hamiltons exceeded the maximum %d\n", MAX_HAMS);
      exit(1);
    }
    hs = &cfac->hams[cfac->nhams];
    cfac->nhams++;
    hs->pj = h->pj;
    hs->nlevs = h->dim;
    hs->nbasis = h->n_basis;
    hs->basis = malloc(sizeof(STATE *)*hs->nbasis);
    for (t = 0; t < h->n_basis; t++) {
      s = (STATE *) ArrayGet(&(sym->states), h->basis[t]);
      hs->basis[t] = s;
    }
    FlagClosed(cfac, hs);
  }

  return 0;

 ERROR:
  printf("ConstructHamilton Error\n");
  return -1;
}

int ValidBasis(cfac_t *cfac, STATE *s, int k, int *kg, int n) {
  int t, m, kb;
  LEVEL *lev;
  STATE *sp;
  SYMMETRY *sym;
  
  t = s->kgroup;
  if (t >= 0) return 0;
  
  if (n > 0) {
    kb = s->kcfg;
    if (kb < 0) return 0;
    kb = GetOrbital(cfac, kb)->n;
    if (kb != n) return 0;
  } else if (n == 0) {
    kb = s->kcfg;
    if (kb >= 0) {
      kb = GetOrbital(cfac, kb)->n;
      if (kb <= 0) return 0;
    }
  } else {
    kb = s->kcfg;
    if (kb < 0) return 0;
    kb = GetOrbital(cfac, kb)->n;
    if (kb > 0) return 0;
  }

  t = -t-1;
  if (kg) {
    lev = GetLevel(cfac, t);
    m = lev->pb;
    sym = GetSymmetry(cfac, lev->pj);
    sp = (STATE *) ArrayGet(&(sym->states), m);
    t = sp->kgroup;
    return InGroups(t, k, kg);
  } else {
    if (t == k) return 1;
    else return 0;
  }
}

int ConstructHamiltonFrozen(cfac_t *cfac,
    int isym, int k, int *kg, int n, int nc, int *kc) {
  HAMILTON *h = cfac->hamiltonian;
  int i, j, t, ncs;
  LEVEL *lev;
  ARRAY *st;
  STATE *s;
  SYMMETRY *sym;
  double r, delta;

  DecodePJ(isym, &i, &j);

  j = 0;
  ncs = 0;
  sym = GetSymmetry(cfac, isym);
  st = &(sym->states);
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if (ValidBasis(cfac, s, k, kg, n)) {
      j++;
    } else if (nc > 0) {
      if (ValidBasis(cfac, s, nc, kc, 0)) {
	j++;
	ncs++;
      }
    }
  }
  
  if (j == ncs) return -1;

  h->pj = isym;

  if (AllocHamMem(h, j, j) == -1) goto ERROR;
      
  j = 0;
  if (ncs > 0) {
    for (t = 0; t < sym->n_states; t++) { 
      s = (STATE *) ArrayGet(st, t);
      if (ValidBasis(cfac, s, nc, kc, 0)) {
	h->basis[j] = t;
	j++;
      }
    }
  }
  for (t = 0; t < sym->n_states; t++) { 
    s = (STATE *) ArrayGet(st, t);
    if (ValidBasis(cfac, s, k, kg, n)) {
      h->basis[j] = t;
      j++;
    }
  }

  for (j = ncs; j < h->dim; j++) {
    t = j*(j+1)/2;
    for (i = ncs; i <= j; i++) {
      r = HamiltonElementFrozen(cfac, isym, h->basis[i], h->basis[j]);
      h->hamilton[i+t] = r;
    }
    for (i = 0; i < j; i++) {
      delta = fabs(h->hamilton[i+t]/h->hamilton[j+t]);
      if (delta < EPS16) h->hamilton[i+t] = 0.0;
    }
  }
  
  for (j = 0; j < ncs; j++) {
    t = j*(j+1)/2;
    for (i = 0; i < j; i++) {
      h->hamilton[i+t] = 0.0;
    }
    s = (STATE *) ArrayGet(st, h->basis[j]);
    lev = GetLevel(cfac, -(s->kgroup+1));
    h->hamilton[j+t] = lev->energy;
  }
  for (i = 0; i < ncs; i++) {
    for (j = ncs; j < h->dim; j++) {
      t = j*(j+1)/2 + i;
      r = HamiltonElementFB(cfac, isym, h->basis[j], h->basis[i]);
      h->hamilton[t] = r;
    }
  }

  return 0;

 ERROR:
  return -1;
}

void AngularFrozen(cfac_t *cfac, int nts, int *ts, int ncs, int *cs) {
  int i, j, kz;
  
  cfac->ang_frozen.nts = nts;
  if (nts > 0) {
    cfac->ang_frozen.ts = malloc(sizeof(int)*nts);
    memcpy(cfac->ang_frozen.ts, ts, sizeof(int)*nts);
  }
  cfac->ang_frozen.ncs = ncs;  
  if (ncs > 0) {
    cfac->ang_frozen.cs = malloc(sizeof(int)*ncs);
    memcpy(cfac->ang_frozen.cs, cs, sizeof(int)*ncs);
  }

  cfac->ang_frozen.nz = malloc(sizeof(int)*nts*nts);
  cfac->ang_frozen.z = malloc(sizeof(ANGULAR_ZMIX *)*nts*nts);
  for (i = 0; i < nts; i++) {
    for (j = 0; j < nts; j++) {
      kz = j*nts + i;
      cfac->ang_frozen.nz[kz] = AngularZMix(cfac, &(cfac->ang_frozen.z[kz]), 
				      ts[i], ts[j], -1, -1);
    }
  }
  if (ncs > 0) {
    cfac->ang_frozen.nzfb = malloc(sizeof(int)*nts*ncs);
    cfac->ang_frozen.zfb = malloc(sizeof(ANGULAR_ZFB *)*nts*ncs);
    cfac->ang_frozen.nzxzfb = malloc(sizeof(int)*nts*ncs);
    cfac->ang_frozen.zxzfb = malloc(sizeof(ANGULAR_ZxZMIX *)*nts*ncs);
    for (i = 0; i < nts; i++) {
      for (j = 0; j < ncs; j++) {
	kz = j*nts + i;
	cfac->ang_frozen.nzfb[kz] = AngularZFreeBound(cfac,
            &(cfac->ang_frozen.zfb[kz]), ts[i], cs[j]);
	cfac->ang_frozen.nzxzfb[kz] = AngularZxZFreeBound(cfac,
            &(cfac->ang_frozen.zxzfb[kz]), ts[i], cs[j]);
      }
    }
  }
}

double HamiltonElementEB(const cfac_t *cfac, int i0, int j0) {
  HAMILTON *h = cfac->hamiltonian;
  int ib, jb;
  int si, sj, mi, mj, pi, pj, ji, jj, ti, tj, kz, nz;
  int i, m, q, q2, jorb0, korb0, jorb1, korb1;
  double r, a, b, c;
  ANGULAR_ZMIX *ang;
  LEVEL *levi, *levj;
  ORBITAL *orb0, *orb1;
  
  ib = h->basis[i0];
  jb = h->basis[j0];

  DecodeBasisEB(ib, &si, &mi);
  DecodeBasisEB(jb, &sj, &mj);
  levi = GetLevel(cfac, si);
  levj = GetLevel(cfac, sj);
  DecodePJ(levi->pj, &pi, &ji);
  DecodePJ(levj->pj, &pj, &jj);
  if (ib == jb) {
    r = levi->energy;
  } else {
    r = 0.0;
  }
  ti = IBisect(si, cfac->ang_frozen.nts, cfac->ang_frozen.ts);
  tj = IBisect(sj, cfac->ang_frozen.nts, cfac->ang_frozen.ts);
  kz = tj*cfac->ang_frozen.nts + ti;
  ang = cfac->ang_frozen.z[kz];
  nz = cfac->ang_frozen.nz[kz];
  for (i = 0; i < nz; i++) {
    if (cfac->b1[0] || cfac->b1[1] || cfac->b1[2]) {
      if (ang[i].k == 2) {
	orb0 = GetOrbital(cfac, ang[i].k0);
	orb1 = GetOrbital(cfac, ang[i].k1);
	GetJLFromKappa(orb0->kappa, &jorb0, &korb0);
	GetJLFromKappa(orb1->kappa, &jorb1, &korb1);
        if (orb0->n == orb1->n && korb0 == korb1) {
	  for (m = 0; m < 3; m++) {
	    if (cfac->b1[m] == 0) continue;
	    
            q = m-1;
	    q2 = 2*q;	
	    
            a = W3j(ji, 2, jj, -mi, -q2, mj);
	    if (a == 0.0) continue;
	    
            a *= cfac->b1[m]*ang[i].coeff;
	    if (IsOdd(abs(ji-mi+q2)/2)) a = -a;
	    
            /* (j|J|j') = sqrt(j*(j + 1)*(2j + 1))*delta(j,j') */
            if (jorb0 == jorb1) {
              b = sqrt(jorb0*(jorb0+2.0)*(jorb0+1.0))/2;
            } else {
              b = 0.0;
            }
            
            /* (s|S|s') = sqrt(3/2)*delta(s,s') */
            c = 1.0023192*sqrt((jorb0+1.0)*(jorb1+1.0))*
                          W6j(korb0, 1, jorb0, 2, jorb1, 1)*sqrt(1.5);
	    if (IsEven((korb0+jorb0+1)/2)) c = -c;
	    
            r += a*(b + c);
	  }
        }
      }      
    }
    if (cfac->e1[0] || cfac->e1[1] || cfac->e1[2]) {
      if (ang[i].k == 2) {
	orb0 = GetOrbital(cfac, ang[i].k0);
	orb1 = GetOrbital(cfac, ang[i].k1);
	GetJLFromKappa(orb0->kappa, &jorb0, &korb0);
	GetJLFromKappa(orb1->kappa, &jorb1, &korb1);
	if (IsOdd((korb0+korb1)/2)) {
	  for (m = 0; m < 3; m++) {
	    double rvme;
            if (cfac->e1[m] == 0) continue;
	    q = m-1;
	    q2 = 2*q;
	    a = W3j(ji, 2, jj, -mi, -q2, mj);
	    if (a == 0.0) continue;
	    a *= ang[i].coeff;
	    if (IsOdd(abs(ji-mi+q2)/2)) a = -a;
	    b = ReducedCL(jorb0, 2, jorb1);
	    c = RadialMoments(cfac, 1, ang[i].k0, ang[i].k1);
	    rvme = a*b*c;
            r += cfac->e1[m]*rvme;
	  }
	}    
      }
    }

    if (cfac->b0 || cfac->b2[0] || cfac->b2[1] || cfac->b2[2] ||
        cfac->b2[3] || cfac->b2[4]) {
      if (ang[i].k == 0) {
	a = W3j(ji, 0, jj, -mi, 0, mj);
	if (a) {
	  orb0 = GetOrbital(cfac, ang[i].k0);
	  orb1 = GetOrbital(cfac, ang[i].k1);
	  GetJLFromKappa(orb0->kappa, &jorb0, &korb0);
	  GetJLFromKappa(orb1->kappa, &jorb1, &korb1);
	  if (IsEven((korb0+korb1)/2)) {
	    a *= cfac->b0*ang[i].coeff;
	    b = ReducedCL(jorb0, 0, jorb1);
	    c = RadialMoments(cfac, 2, ang[i].k0, ang[i].k1);
	    if (IsOdd((ji-mi)/2)) a = -a;
	    r += a*b*c;
	  }
	}
      }
      if (ang[i].k == 4) {	
	orb0 = GetOrbital(cfac, ang[i].k0);
	orb1 = GetOrbital(cfac, ang[i].k1);
	GetJLFromKappa(orb0->kappa, &jorb0, &korb0);
	GetJLFromKappa(orb1->kappa, &jorb1, &korb1);
	if (IsEven((korb0+korb1)/2)) {
	  for (m = 0; m < 5; m++) {
	    if (cfac->b2[m] == 0) continue;
	    q = m-2;
	    q2 = q*2;
	    a = W3j(ji, 4, jj, -mi, -q2, mj);
	    if (a == 0.0) continue;
	    a *= cfac->b2[m]*ang[i].coeff;
	    b = ReducedCL(jorb0, 4, jorb1);
	    c = RadialMoments(cfac, 2, ang[i].k0, ang[i].k1);
	    if (IsOdd((ji-mi)/2)) a = -a;
	    r += a*b*c;
	  }
	}
      }
    }
  }

  return r;
}

double HamiltonElementFB(cfac_t *cfac, int isym, int isi, int isj) {
  STATE *si, *sj;
  double r, sd, se, a;
  int i, ti, tj, ji, jj, j0, k0, k1, nz1, nz2;
  int kz = 0, ks[4];
  ORBITAL *orb0, *orb1;
  ANGULAR_ZFB *a1;
  ANGULAR_ZxZMIX *a2;
  SYMMETRY *sym;
  LEVEL *lev1, *lev2;

  sym = GetSymmetry(cfac, isym);
  si = (STATE *) ArrayGet(&(sym->states), isi);
  sj = (STATE *) ArrayGet(&(sym->states), isj);
  r = 0.0;
  ti = si->kgroup;
  k0 = si->kcfg;
  orb0 = GetOrbital(cfac, k0);
  j0 = GetJFromKappa(orb0->kappa);
  tj = sj->kgroup;
  ti = -ti-1;
  tj = -tj-1;
  lev1 = GetLevel(cfac, ti);
  lev2 = GetLevel(cfac, tj);
  ji = lev1->pj;
  jj = lev2->pj;
  DecodePJ(ji, NULL, &ji);
  DecodePJ(jj, NULL, &jj);

  if (cfac->ang_frozen.nts > 0) {
    ti = IBisect(ti, cfac->ang_frozen.nts, cfac->ang_frozen.ts);
    tj = IBisect(tj, cfac->ang_frozen.ncs, cfac->ang_frozen.cs);
    kz = tj*cfac->ang_frozen.nts + ti;
    nz1 = cfac->ang_frozen.nzfb[kz];
    a1 = cfac->ang_frozen.zfb[kz];
  } else {
    nz1 = AngularZFreeBound(cfac, &a1, ti, tj);
  }
  for (i = 0; i < nz1; i++) {
    k1 = a1[i].kb;
    orb1 = GetOrbital(cfac, k1);
    if (orb0->kappa == orb1->kappa) {
      ResidualPotential(cfac, &a, k0, k1);
      a += QED1E(cfac, k0, k1);
      a *= a1[i].coeff;
      if (IsOdd((ji-jj+j0)/2)) a = -a;
      r += a;
    }
  }
  if (nz1 > 0 && cfac->ang_frozen.nts == 0) free(a1);
  if (cfac->ang_frozen.nts > 0) {
    nz2 = cfac->ang_frozen.nzxzfb[kz];
    a2 = cfac->ang_frozen.zxzfb[kz];
  } else {
    nz2 = AngularZxZFreeBound(cfac, &a2, ti, tj);
  }
  for (i = 0; i < nz2; i++) {
    if (j0 == a2[i].k0) {
      ks[0] = k0;
      ks[1] = a2[i].k2;
      ks[2] = a2[i].k1;
      ks[3] = a2[i].k3;
      SlaterTotal(cfac, &sd, &se, NULL, ks, a2[i].k, 0);
      a = (sd + se)*a2[i].coeff;
      r += a;
    }
  }
  if (nz2 > 0 && cfac->ang_frozen.nts == 0) free(a2);

  r /= sqrt(jj + 1.0);
  
  return r;
}

double HamiltonElementFrozen(cfac_t *cfac, int isym, int isi, int isj) {
  STATE *si, *sj;
  double r, r0, sd, se, a;
  int i, ti, tj, ji1, ji2, jj1, jj2, ki2, kj2, j, nz;
  int kz, ks[4];
  ORBITAL *orbi, *orbj;
  ANGULAR_ZMIX *ang;
  SYMMETRY *sym;
  LEVEL *lev1, *lev2;

  sym = GetSymmetry(cfac, isym);
  si = (STATE *) ArrayGet(&(sym->states), isi);
  sj = (STATE *) ArrayGet(&(sym->states), isj);
  r = 0.0;
  ti = si->kgroup;
  tj = sj->kgroup;
  ti = -ti-1;
  tj = -tj-1;
  orbi = GetOrbital(cfac, si->kcfg);
  orbj = GetOrbital(cfac, sj->kcfg);
  lev1 = GetLevel(cfac, ti);
  lev2 = GetLevel(cfac, tj);
  ji1 = lev1->pj;
  jj1 = lev2->pj;
  DecodePJ(ji1, NULL, &ji1);
  DecodePJ(jj1, NULL, &jj1);
  GetJLFromKappa(orbi->kappa, &ji2, &ki2);
  GetJLFromKappa(orbj->kappa, &jj2, &kj2);
  j = si->kstate;
  if (si->kgroup == sj->kgroup) { 
    if (ji2 == jj2 && ki2 == kj2) {
      ResidualPotential(cfac, &a, si->kcfg, sj->kcfg);
      r += a;
      r0 = QED1E(cfac, si->kcfg, sj->kcfg);
      r += r0;
    } 
    if (si->kcfg == sj->kcfg) {
      r += lev1->energy;
      r += orbi->energy;
    }
  }
 
  ks[1] = si->kcfg;
  ks[3] = sj->kcfg;

  if (cfac->ang_frozen.nts > 0) {
    ti = IBisect(ti, cfac->ang_frozen.nts, cfac->ang_frozen.ts);
    tj = IBisect(tj, cfac->ang_frozen.nts, cfac->ang_frozen.ts);
    kz = tj*cfac->ang_frozen.nts + ti;
    nz = cfac->ang_frozen.nz[kz];
    ang = cfac->ang_frozen.z[kz];
  } else {
    nz = AngularZMix(cfac, &ang, ti, tj, -1, -1);
  }
  a = 0.0;
  for (i = 0; i < nz; i++) {
    if (fabs(ang[i].coeff) < EPS30) continue;
    r0 = W6j(ji1, ji2, j, jj2, jj1, ang[i].k);
    if (fabs(r0) < EPS30) continue;
    ks[0] = ang[i].k0;
    ks[2] = ang[i].k1;
    SlaterTotal(cfac, &sd, &se, NULL, ks, ang[i].k, 0);
    r0 *= ang[i].coeff*(sd+se);
    a += r0;
  }

  if (IsOdd((ji2 + jj1 + j)/2)) a = -a;
  r += a;

  if (nz > 0 && cfac->ang_frozen.nts == 0) {
    free(ang);
  } 

  return r;
} 

double HamiltonElement(cfac_t *cfac, int isym, int isi, int isj) {
  double r1, r2;
  
  HamiltonElement1E2E(cfac, isym, isi, isj, &r1, &r2);
  
  return r1 + r2;
}

void HamiltonElement1E2E(cfac_t *cfac,
  int isym, int isi, int isj, double *x1, double *x2) { 
  CONFIG *ci, *cj;
  int ki, kj;
  SYMMETRY *sym;
  STATE *si, *sj;
  SHELL_STATE *sbra, *sket;
  SHELL *bra;
  int n_shells, i, j;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];
  double x, r;
  int phase;

  *x1 = 0.0;
  *x2 = 0.0;
  sym = GetSymmetry(cfac, isym);
  si = (STATE *) ArrayGet(&(sym->states), isi);
  sj = (STATE *) ArrayGet(&(sym->states), isj);
  
  ci = GetConfig(cfac, si);
  if (ci->n_shells == 0) return;
  cj = GetConfig(cfac, sj);
  if (cj->n_shells == 0) return;
  
  switch (cfac->confint) {
  case 1:
    if (ci != cj) return;
  case 2:
    if (ci->nnrs != cj->nnrs) return;
    else {
      if (memcmp(ci->nrs, cj->nrs, sizeof(int)*ci->nnrs)) return;
    }
  case 3:
    if (si->kgroup != sj->kgroup) return;
  }
    
  ki = si->kstate;
  kj = sj->kstate;

  idatum = NULL;
  n_shells = GetInteract(cfac, &idatum, &sbra, &sket, 
			 si->kgroup, sj->kgroup,
			 si->kcfg, sj->kcfg,
			 ki, kj, 0);
  if (n_shells <= 0) return;
  memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
  bra = idatum->bra;
  phase = idatum->phase;

  if (s[0].index >= 0 && s[3].index >= 0) {
    r = Hamilton2E(cfac, n_shells, sbra, sket, s);
    *x2 += r;
  } else if( s[0].index >= 0) {
    r = Hamilton1E(cfac, n_shells, sbra, sket, s);
    *x1 += r;
    for (i = 0; i < n_shells; i++) {
      s[2].index = n_shells - i - 1;
      s[3].index = s[2].index;
      s[2].n = bra[i].n;
      s[3].n = s[2].n;
      s[2].kappa = bra[i].kappa;
      s[3].kappa = s[2].kappa;
      s[2].j = GetJ(bra+i);
      s[3].j = s[2].j;
      s[2].kl = GetL(bra+i);
      s[3].kl = s[2].kl;
      s[2].nq_bra = GetNq(bra+i);
      if (s[2].index == s[0].index) {
	s[2].nq_ket = s[2].nq_bra - 1;
      } else if (s[2].index == s[1].index) {
	s[2].nq_ket = s[2].nq_bra + 1;
      } else {
	s[2].nq_ket = s[2].nq_bra;
      }
      if (s[2].nq_bra <= 0 || s[2].nq_ket <= 0 ||
	  s[2].nq_bra > s[2].j+1 || s[2].nq_ket > s[2].j+1) {
	continue;
      }
      s[3].nq_bra = s[2].nq_bra;
      s[3].nq_ket = s[2].nq_ket;
      r = Hamilton2E(cfac, n_shells, sbra, sket, s);
      *x2 += r;
    }
  } else {
    for (i = 0; i < n_shells; i++) {
      s[0].index = n_shells - i - 1;
      s[1].index = s[0].index;
      s[0].n = bra[i].n;
      s[1].n = s[0].n;
      s[0].kappa = bra[i].kappa;
      s[1].kappa = s[0].kappa;
      s[0].j = GetJ(bra+i);
      s[1].j = s[0].j;
      s[0].kl = GetL(bra+i);
      s[1].kl = s[0].kl;
      s[0].nq_bra = GetNq(bra+i);
      s[0].nq_ket = s[0].nq_bra;
      s[1].nq_bra = s[0].nq_bra;
      s[1].nq_ket = s[1].nq_bra;      
      r = Hamilton1E(cfac, n_shells, sbra, sket, s);
      *x1 += r;      
      for (j = 0; j <= i; j++) {
	s[2].nq_bra = GetNq(bra+j);
	if (j == i && s[2].nq_bra < 2) continue;
	s[2].nq_ket = s[2].nq_bra;
	s[3].nq_bra = s[2].nq_bra;
	s[3].nq_ket = s[3].nq_bra;
	s[2].index = n_shells - j - 1;
	s[3].index = s[2].index;
	s[2].n = bra[j].n;
	s[3].n = s[2].n;
	s[2].kappa = bra[j].kappa;
	s[3].kappa = s[2].kappa;
	s[2].j = GetJ(bra+j);
	s[3].j = s[2].j;
	s[2].kl = GetL(bra+j);
	s[3].kl = s[2].kl;
	r = Hamilton2E(cfac, n_shells, sbra, sket, s);
	*x2 += r;
      }
    }
  }
  /* the prefactor in the Wigner-Eckart theorem should be included in 
     the matrix element. for a scalar operator, this is [J]^{-1/2} */
  x = sqrt(sbra[0].totalJ + 1.0);
  if (IsOdd(phase)) x = -x;
  *x1 /= x;
  *x2 /= x;

  if (isi == isj) {
    *x1 += ci->delta;
  }

  free(sbra);
  free(sket);
}

int SlaterCoeff(cfac_t *cfac, char *fn, int nlevs, int *ilevs, 
		int na, SHELL *sa, int nb, SHELL *sb) {
  FILE *f;
  int m, i, j, i0, i1, k0, k1, q0, q1;
  int na2, nb2, nab2, n_shells, vnl;
  double a, *coeff;
  CONFIG *c0, *c1;
  STATE *s0, *s1;
  SYMMETRY *sym;
  LEVEL *lev;
  SHELL_STATE *sbra, *sket;
  SHELL *bra;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];  
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  char nc[LEVEL_NAME_LEN];
  
  f = fopen(fn, "w");
  if (f == NULL) return -1;
  
  na2 = 2*na;
  nb2 = 2*nb;
  nab2 = na2*nb2;
  coeff = (double *) malloc(sizeof(double)*nab2*4);

  for (m = 0;  m < nlevs; m++) {
    for (i = 0; i < nab2*4; i++) {
      coeff[i] = 0.0;
    }
    lev = GetLevel(cfac, ilevs[m]);
    sym = GetSymmetry(cfac, lev->pj);
    for (i0 = 0; i0 < lev->n_basis; i0++) {
      s0 = (STATE *) ArrayGet(&(sym->states), lev->basis[i0]);
      c0 = GetConfig(cfac, s0);
      k0 = s0->kstate;
      for (i1 = 0; i1 < lev->n_basis; i1++) {
	a = lev->mixing[i0] * lev->mixing[i1];
	if (fabs(a) < cfac->angz_cut) continue;
	s1 = (STATE *) ArrayGet(&(sym->states), lev->basis[i1]);
	c1 = GetConfig(cfac, s1);
	k1 = s1->kstate;
	idatum = NULL;
	n_shells = GetInteract(cfac, &idatum, &sbra, &sket, s0->kgroup, s1->kgroup, 
			       s0->kcfg, s1->kcfg, k0, k1, 0);
	if (n_shells <= 0) continue;
	memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
	bra = idatum->bra;
	if (IsOdd(idatum->phase)) a = -a;
	if (s[0].index >= 0 && s[3].index >= 0) {
	  AddSlaterCoeff(cfac, coeff, a, n_shells, sbra, sket, s, na, sa, nb, sb);
	} else if (s[0].index >= 0) {
	  q0 = ShellIndex(s[0].n, s[0].kappa, na, sa);
	  if (q0 >= 0) {
	    for (j = 0; j < nb; j++) {
	      i = ShellIndex(sb[j].n, sb[j].kappa, n_shells, bra);
	      if (i >= 0) {
		s[2].index = n_shells - i - 1;
		if (s[2].index == s[0].index) continue;
		s[3].index = s[2].index;
		s[2].n = bra[i].n;
		s[3].n = s[2].n;
		s[2].kappa = bra[i].kappa;
		s[3].kappa = s[2].kappa;
		s[2].j = GetJ(bra+i);
		s[3].j = s[2].j;
		s[2].kl = GetL(bra+i);
		s[3].kl = s[2].kl;
		s[2].nq_bra = GetNq(bra+i);
		if (s[2].index == s[0].index) {
		  s[2].nq_ket = s[2].nq_bra - 1;
		} else if (s[2].index == s[1].index) {
		  s[2].nq_ket = s[2].nq_bra + 1;
		} else {
		  s[2].nq_ket = s[2].nq_bra;
		}
		if (s[2].nq_bra < 0 || s[2].nq_ket < 0 ||
		    s[2].nq_bra > s[2].j+1 || s[2].nq_ket > s[2].j+1) {
		  continue;
		}
		s[3].nq_bra = s[2].nq_bra;
		s[3].nq_ket = s[2].nq_ket;
		AddSlaterCoeff(cfac, coeff, a, n_shells, sbra, sket, s, na, sa, nb, sb);
	      } else if (s[0].j == s[1].j) {
		s[2].n = sb[j].n;
		s[2].kappa = sb[j].kappa;
		s[2].nq_bra = 0;
		s[2].nq_ket = 0;
		GetJLFromKappa(s[2].kappa, &(s[2].j), &(s[2].kl));
		memcpy(s+3, s+2, sizeof(INTERACT_SHELL));
		AddSlaterCoeff(cfac, coeff, a, n_shells, sbra, sket, s, na, sa, nb, sb);
	      }
	    }
	  }
	} else {
	  for (j = 0; j < na; j++) {
	    i = ShellIndex(sa[j].n, sa[j].kappa, n_shells, bra);
	    if (i < 0 || bra[i].nq == 0) continue;
	    s[0].index = n_shells - i - 1;
	    s[1].index = s[0].index;
	    s[0].n = bra[i].n;
	    s[1].n = s[0].n;
	    s[0].kappa = bra[i].kappa;
	    s[1].kappa = s[0].kappa;
	    s[0].j = GetJ(bra+i);
	    s[1].j = s[0].j;
	    s[0].kl = GetL(bra+i);
	    s[1].kl = s[0].kl;
	    s[0].nq_bra = GetNq(bra+i);
	    s[0].nq_ket = s[0].nq_bra;
	    s[1].nq_bra = s[0].nq_bra;
	    s[1].nq_ket = s[1].nq_bra; 	    
	    for (q0 = 0; q0 < nb; q0++) {
	      q1 = ShellIndex(sb[q0].n, sb[q0].kappa, n_shells, bra);
	      if (q1 == i) continue;
	      if (q1 >= 0) {
		s[2].nq_bra = GetNq(bra+q1);
		s[2].nq_ket = s[2].nq_bra;
		s[3].nq_bra = s[2].nq_bra;
		s[3].nq_ket = s[3].nq_bra;
		s[2].index = n_shells - q1 - 1;
		s[3].index = s[2].index;
		s[2].n = bra[q1].n;
		s[3].n = s[2].n;
		s[2].kappa = bra[q1].kappa;
		s[3].kappa = s[2].kappa;
		s[2].j = GetJ(bra+q1);
		s[3].j = s[2].j;
		s[2].kl = GetL(bra+q1);
		s[3].kl = s[2].kl;
		s[2].nq_bra = GetNq(bra+q1);
		s[3].nq_bra = s[2].nq_bra;
		s[2].nq_ket = s[2].nq_bra;
		s[3].nq_ket = s[2].nq_bra;
		AddSlaterCoeff(cfac, coeff, a, n_shells, sbra, sket, s, na, sa, nb, sb);
	      } else {		
		s[2].n = sb[q0].n;
		s[2].kappa = sb[q0].kappa;
		s[2].nq_bra = 0;
		s[2].nq_ket = 0;
		GetJLFromKappa(s[2].kappa, &(s[2].j), &(s[2].kl));
		memcpy(s+3, s+2, sizeof(INTERACT_SHELL));
		AddSlaterCoeff(cfac, coeff, a, n_shells, sbra, sket, s, na, sa, nb, sb);
	      }
	    }
	  }
	}
	free(sbra);
	free(sket);	
      }
    }
    
    DecodePJ(lev->pj, &i, &j);
    a = sqrt(j+1.0);

    s0 = (STATE *) ArrayGet(&(sym->states), lev->pb);
    ConstructLevelName(cfac, name, sname, nc, &vnl, s0);
    fprintf(f, "# %6d %1d %3d   %-s\n",
	    ilevs[m], i, j, name);
    for (i0 = 0; i0 < na; i0++) {
      k0 = GetLFromKappa(sa[i0].kappa)/2;
      for (i1 = 0; i1 < nb; i1++) {
	k1 = GetLFromKappa(sb[i1].kappa)/2;
	for (q0 = 0; q0 < 2; q0++) {
	  for (q1 = 0; q1 < 2; q1++) {
	    i = (i0*2+q0)*nb2 + i1*2+q1;
	    for (j = 0; j < 4; j++) {
	      if (coeff[i+j*nab2] != 0) break;
	    }
	    if (j < 4) {
	      fprintf(f, "  %6d %2d %2d %2d %d %2d %2d %2d %d",
		      ilevs[m], sa[i0].n, sa[i0].kappa, k0, q0,
		      sb[i1].n, sb[i1].kappa, k1, q1);
	      for (j = 0; j < 4; j++) {
		fprintf(f, " %12.5E", coeff[i+j*nab2]/a);
	      }
	      fprintf(f, "\n");
	    }
	  }
	}
      }
    }
  }

  fclose(f);
  free(coeff);
  
  return 0;
}
	  
void AddSlaterCoeff(const cfac_t *cfac, double *coeff, double a, int n_shells, 
		    SHELL_STATE *sbra, SHELL_STATE *sket, 
		    INTERACT_SHELL *s, int na, SHELL *sa, int nb, SHELL *sb) {
  int t, nk, *kk, j, i0, i1, i, nk0, k, *kk0;
  int k0, k1, k2, k3, js[4], na2, nb2, nab2;
  double e, z0, *y, *ang;
  int kmax = GetMaxRank(cfac);

  if (s[0].kl != s[1].kl) return;
  if (s[2].kl != s[3].kl) return;
  if (s[0].n != s[1].n) return;
  if (s[2].n != s[3].n) return;

  k0 = ShellIndex(s[0].n, s[0].kappa, na, sa);
  k1 = ShellIndex(s[2].n, s[2].kappa, nb, sb);
  k2 = ShellIndex(s[1].n, s[1].kappa, na, sa);
  k3 = ShellIndex(s[3].n, s[3].kappa, nb, sb);
  if (k0 < 0 || k1 < 0 || k2 < 0 || k3 < 0) return;
  
  if (s[0].n == s[2].n && s[0].kl == s[2].kl) return;
  if (s[1].n == s[3].n && s[1].kl == s[3].kl) return;
  
  js[0] = s[0].j;
  js[1] = s[2].j;
  js[2] = s[1].j;
  js[3] = s[3].j;

  i0 = k0*2;
  if (k2 != k0) i0++;
  i1 = k1*2;
  if (k3 != k1) i1++;  
  na2 = 2*na;
  nb2 = 2*nb;
  nab2 = na2*nb2;
  j = i0*nb2 + i1;

  if (s[2].nq_bra > 0) {
    nk = AngularZxZ0(&ang, &kk, 0, n_shells, sbra, sket, s, kmax);
    for (i = 0; i < nk; i++) {    
      if (fabs(ang[i]) < EPS30) continue;    
      for (t = 2; t <= 4; t += 2) {
	e = W6j(js[0], js[2], kk[i], js[1], js[3], t);
	if (fabs(e) < EPS30) continue;
	e *= ReducedCL(js[0], t, js[3]); 
	e *= ReducedCL(js[1], t, js[2]);
	e *= (kk[i] + 1.0);
	if (IsOdd((t+kk[i])/2)) e = -e;
	if (fabs(e) < EPS30) continue;
	e *= a*ang[i];
	k = (t/2-1)*nab2;
	coeff[j + k] += e;
      }
    }

    if (nk > 0) {
      free(ang);
      free(kk);
    }
  }

  if (k1 == k3) {
    z0 = 0.0;
    nk0 = 1;
    k = 0; 
    kk0 = &k;
    y = &z0;
    nk0 = AngularZ(&y, &kk0, nk0, n_shells, sbra, sket, s, s+1, kmax);
    if (nk0 > 0) {
      z0 /= sqrt(s[0].j + 1.0);
      if (IsOdd((s[0].j - s[2].j)/2)) z0 = -z0;
      for (t = 2; t <= 4; t += 2) {
	e = ReducedCL(js[0], t, js[3]);
	e *= ReducedCL(js[1], t, js[2]);
	k = (t/2+1)*nab2;
	coeff[j+k] += a*z0*e;
      }
    }
  }
}

double Hamilton1E(cfac_t *cfac, int n_shells, SHELL_STATE *sbra, SHELL_STATE *sket,
		  INTERACT_SHELL *s) {
  int nk0, k;
  int *k0;
  double *x, z0, r0;
  int k1, k2;
  int kmax = GetMaxRank(cfac);

  if (s[0].j != s[1].j ||
      s[0].kl != s[1].kl) return 0.0;
  nk0 = 1;
  k = 0;
  k0 = &k;
  x = &z0;
  nk0 = AngularZ(&x, &k0, nk0, n_shells, sbra, sket, s, s+1, kmax);
  if (fabs(z0) < EPS30) return 0.0;
  k1 = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
  k2 = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
  ResidualPotential(cfac, &r0, k1, k2);
  if (k1 == k2) r0 += (GetOrbital(cfac, k1))->energy;
  r0 += QED1E(cfac, k1, k2);
  z0 *= sqrt(s[0].j + 1.0);

  r0 *= z0;
  return r0;
}

double Hamilton2E(cfac_t *cfac,
    int n_shells, SHELL_STATE *sbra, SHELL_STATE *sket, INTERACT_SHELL *s) {
  int nk0, nk, *kk, k, *kk0, i;
  double *ang;
  double se, sd, x;
  double z0, *y;
  int ks[4], js[4];
  int kmax = GetMaxRank(cfac);

  js[0] = 0;
  js[1] = 0;
  js[2] = 0;
  js[3] = 0;
  
  ks[0] = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
  ks[1] = OrbitalIndex(cfac, s[2].n, s[2].kappa, 0.0);
  ks[2] = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
  ks[3] = OrbitalIndex(cfac, s[3].n, s[3].kappa, 0.0);

  z0 = 0.0;
  nk0 = 0;

  if (ks[1] == ks[2]) {
    nk0 = 1;
    k = 0;
    kk0 = &k;
    y = &z0;
    nk0 = AngularZ(&y, &kk0, nk0, n_shells, sbra, sket, s, s+3, kmax);
    if (nk0 > 0) {
      z0 /= sqrt(s[0].j + 1.0);
      if (IsOdd((s[0].j - s[2].j)/2)) z0 = -z0;
    }
  }

  x = 0.0;    
  nk = AngularZxZ0(&ang, &kk, 0, n_shells, sbra, sket, s, kmax);
  for (i = 0; i < nk; i++) {
    sd = 0;
    se = 0;
    if (fabs(ang[i]) > EPS30) {
      SlaterTotal(cfac, &sd, &se, js, ks, kk[i], 0);
      x += ang[i] * (sd+se);
    } else if (nk0 > 0) {
      SlaterTotal(cfac, &sd, NULL, js, ks, kk[i], 0);
    }
    if (nk0 > 0) x -= z0 * sd;
  }

  if (nk > 0) {
    free(ang);
    free(kk);
  }

  return x;
}

int TestHamilton(cfac_t *cfac) {
  SYMMETRY *sym;
  int i, j, t;
  double r1;

  for (t = 0; t < MAX_SYMMETRIES; t++) {
    sym = GetSymmetry(cfac, t);
    for (i = 0; i < sym->n_states; i++) {
      for (j = 0; j < sym->n_states; j++) {
	r1 = HamiltonElement(cfac, t, i, j);
	printf("HAM: %3d %d %d %10.3E\n", t, i, j, r1);
      }
    }
  }
  return 0;
}

int DiagonalizeHamilton(cfac_t *cfac) {
  HAMILTON *h = cfac->hamiltonian;
  gsl_matrix *am, *evec;
  gsl_vector_view vv;
  gsl_eigen_symmv_workspace *wsp;
  double *w;
  double *z;
  double *mixing = NULL;
  int n, m;
  int info;
  int i, j;

  n = h->dim;
  m = h->n_basis;
  
  if (m > n) {
    printf("m > n in DiagonalizeHamilton(), %d %d\n", m, n);
    abort();
  }
  
  if (cfac->confint == -1) {
    /* no configuration interaction at all */
    mixing = h->mixing+n;
    for (i = 0; i < n; i++) {
      h->mixing[i] = h->hamilton[i];
      for (j = 0; j < n; j++) {
	if (i == j) *mixing = 1.0;
	else *mixing = 0.0;
	mixing++;
      }
    }
    return 0;
  }

  mixing = h->mixing;
  w = mixing;
  z = mixing + n;

  wsp = gsl_eigen_symmv_alloc(n);

  am   = gsl_matrix_alloc(n, n);
  evec = gsl_matrix_alloc(n, n);

  for (j = 0; j < h->dim; j++) {
    int t = j*(j+1)/2;
    for (i = 0; i <= j; i++) {
      gsl_matrix_set(am, j, i, h->hamilton[i + t]);
    }
  }

  vv = gsl_vector_view_array(w, n);

  info = gsl_eigen_symmv(am, &vv.vector, evec, wsp);
  gsl_eigen_symmv_free(wsp);
  gsl_matrix_free(am);

  gsl_eigen_symmv_sort(&vv.vector, evec, GSL_EIGEN_SORT_VAL_ASC);

  for (j = 0; j < h->dim; j++) {
    for (i = 0; i < h->dim; i++) {
      z[j*h->dim + i] = gsl_matrix_get(evec, i, j);
    }
  }
  
  gsl_matrix_free(evec);

  if (info) {
    return -1;
  } else {
    return 0;
  }
}

int AddToLevels(cfac_t *cfac, int ng, int *kg) {
  HAMILTON *h = cfac->hamiltonian;
  int i, d, j, k, t, m;
  LEVEL lev;
  SYMMETRY *sym;
  STATE *s, *s1;
  double *mix, a;
  
  if (h->basis == NULL ||
      h->mixing == NULL) return -1;
  d = h->dim;
  mix = h->mixing + d;

  if (h->pj < 0) {
    j = cfac->n_eblevels;
    for (i = 0; i < d; i++) {
      k = GetPrincipleBasis(mix, d, NULL);
      lev.energy = h->mixing[i];
      lev.pj = h->pj;
      lev.iham = -1;
      lev.ilev = j;
      lev.pb = h->basis[k];
      lev.ibasis = (short *) malloc(sizeof(short)*h->n_basis);
      lev.basis = (int *) malloc(sizeof(int)*h->n_basis);
      lev.mixing = (double *) malloc(sizeof(double)*h->n_basis);
      lev.n_basis = h->n_basis;
      for (t = 0; t < h->n_basis; t++) {
	lev.ibasis[t] = t;
	lev.basis[t] = h->basis[t];
	lev.mixing[t] = mix[t];
      }
      m = h->n_basis;
      SortMixing(0, m, &lev, NULL);
      GetPrincipleBasis(lev.mixing, m, lev.kpb);      
    
      if (ArrayAppend(cfac->eblevels, &lev) == NULL) {
	printf("Not enough memory for levels array\n");
	exit(1);
      }
      j++;
      mix += h->n_basis;
    }

    cfac->n_eblevels = j;
    if (i < d-1) return -2;
    return 0;
  }

  j = cfac->n_levels;
  sym = GetSymmetry(cfac, h->pj);  
  for (i = 0; i < d; i++) {
    k = GetPrincipleBasis(mix, d, NULL);
    s = (STATE *) ArrayGet(&(sym->states), h->basis[k]);
    if (ng > 0) {      
      if (!InGroups(s->kgroup, ng, kg)) {
	m = 0;
	if (cfac->mix_cut2 < 1.0) {
	  a = fabs(cfac->mix_cut2*mix[k]);
	  for (t = 0; t < h->n_basis; t++) {
	    if (fabs(mix[t]) >= a && t != k) {
	      s1 = (STATE *) ArrayGet(&(sym->states), h->basis[t]);
	      if (InGroups(s1->kgroup, ng, kg)) {
		m = 1;
		break;
	      }
	    }
	  }
	}
	if (m == 0) {
	  mix += h->n_basis;
	  continue;
	}
      }
    }
    lev.energy = h->mixing[i];
    lev.pj = h->pj;
    lev.iham = cfac->nhams-1;
    lev.ilev = i;
    lev.pb = h->basis[k];
    lev.ibasis = (short *) malloc(sizeof(short)*h->n_basis);
    lev.basis = (int *) malloc(sizeof(int)*h->n_basis);
    lev.mixing = (double *) malloc(sizeof(double)*h->n_basis);
    a = fabs(cfac->mix_cut * mix[k]);
    for (t = 0, m = 0; t < h->n_basis; t++) {
      if (fabs(mix[t]) < a) continue;
      lev.ibasis[m] = t;
      lev.basis[m] = h->basis[t];
      lev.mixing[m] = mix[t];
      m++;
    }
    lev.n_basis = m;
    if (m < t) {
      lev.ibasis = (short *) realloc(lev.ibasis, sizeof(short)*m);
      lev.basis = (int *) realloc(lev.basis, sizeof(int)*m);
      lev.mixing = (double *) realloc(lev.mixing, sizeof(double)*m);
    }
    SortMixing(0, m, &lev, sym);
    GetPrincipleBasis(lev.mixing, m, lev.kpb);

    if (s->kgroup < 0) {
      lev.ibase = -(s->kgroup + 1);
      lev.iham = -1;
    }

    if (ArrayAppend(cfac->levels, &lev) == NULL) {
      printf("Not enough memory for levels array\n");
      exit(1);
    }
    j++;
    mix += h->n_basis;
  }

  cfac->n_levels = j;
  if (i < d-1) return -2;

  return 0;
}

void CutMixing(cfac_t *cfac, int nlev, int *ilev, int n, int *kg, double c) {
  int i, m, t;
  SYMMETRY *sym;
  STATE *s;
  LEVEL *lev;

  for (i = 0; i < nlev; i++) {
    lev = GetLevel(cfac, ilev[i]);
    m = 0;
    sym = GetSymmetry(cfac, lev->pj);
    for (t = 0; t < lev->n_basis; t++) {
      if (fabs(lev->mixing[t]) < c) continue;
      s = (STATE *) ArrayGet(&(sym->states), lev->basis[t]);
      if (n > 0 && !InGroups(s->kgroup, n, kg)) continue;
      lev->ibasis[m] = lev->ibasis[t];
      lev->basis[m] = lev->basis[t];
      lev->mixing[m] = lev->mixing[t];
      m++;
    }
    if (m < lev->n_basis) {
      lev->n_basis = m;
      lev->ibasis = (short *) realloc(lev->ibasis, sizeof(short)*m);
      lev->basis = (int *) realloc(lev->basis, sizeof(int)*m);
      lev->mixing = (double *) realloc(lev->mixing, sizeof(double)*m);      
      SortMixing(0, m, lev, sym);
      GetPrincipleBasis(lev->mixing, m, lev->kpb);
    }
  }
}
  
static int CompareBasis(double m1, double m2, SYMMETRY *sym) { 
  if (fabs(m1) > fabs(m2)) return 1;
  else if (fabs(m1) < fabs(m2)) return -1;
  else return 0;
}  

int SortMixing(int start, int n, LEVEL *lev, SYMMETRY *sym) {
  short *ibasis;
  int *basis;
  double *mix;
  int i, j, i0, j0, t;
  int *b1, *b2, *bp;
  short *s1, *s2, *sp;
  double *m1, *m2, *mp, tmp;

  ibasis = lev->ibasis;
  basis = lev->basis;
  mix = lev->mixing;
  while (1 < n) {
    i = start;
    j = start + n - 1;
    m1 = mix + i;
    m2 = mix + j;
    b1 = basis + i;
    b2 = basis + j;
    s1 = ibasis + i;
    s2 = ibasis + j;
    mp = m2;    
    bp = b2;
    sp = s2;
    
    while (i < j) {
      while (i < j) {
	if (CompareBasis(*m1, *mp, sym) < 0) break;
	i++;
	m1 = mix + i;
	b1 = basis + i;
	s1 = ibasis + i;
      }
      while (i < j) {
	if (CompareBasis(*mp, *m2, sym) < 0) break;
	j--;
	m2 = mix + j;
	b2 = basis + j;
	s2 = ibasis + j;	
      }
      if (i < j) {
	tmp = *m1;
	*m1 = *m2;
	*m2 = tmp;
	t = *b1;
	*b1 = *b2;
	*b2 = t;
	t = *s1;
	*s1 = *s2;
	*s2 = t;
	i++;
	m1 = mix + i;
	b1 = basis + i;
	s1 = ibasis + i;
      }
    }
    if (CompareBasis(*m1, *mp, sym)) {
      tmp = *m1;
      *m1 = *mp;
      *mp = tmp;
      t = *b1;
      *b1 = *bp;
      *bp = t;
      t = *s1;
      *s1 = *sp;
      *sp = t;
    }

    i0 = i - start;
    j0 = n - i0 - 1;
    if (j0 < i0) {
      if (1 < j0) {
	SortMixing(i+1, j0, lev, sym);
      }
      n = i0;
    } else {
      if (1 < i0) {
	SortMixing(start, i0, lev, sym);
      }
      start = i+1;
      n = j0;
    }
  }
  return 0;
}
  
int AddECorrection(cfac_t *cfac, int iref, int ilev, double e, int nmin) {
  ECORRECTION c;

  c.iref = iref;
  c.ilev = ilev;
  c.e = e;
  c.nmin = nmin;
  ArrayAppend(cfac->ecorrections, &c);
  cfac->ncorrections += 1;

  return 0;
}

LEVEL *GetEBLevel(const cfac_t *cfac, int k) {
  return (LEVEL *) ArrayGet(cfac->eblevels, k);
}

LEVEL *GetLevel(const cfac_t *cfac, int k) {
  return (LEVEL *) ArrayGet(cfac->levels, k);
}

int LevelTotalJ(cfac_t *cfac, int k) {
  int i;
  i = GetLevel(cfac, k)->pj;
  DecodePJ(i, NULL, &i);
  return i;
}

int GetNumEBLevels(const cfac_t *cfac) {
  return cfac->n_eblevels;
}

int GetNumLevels(const cfac_t *cfac) {
  return cfac->n_levels;
}

int GetPrincipleBasis(double *mix, int d, int *kpb) {
  int i, k = 0, t, q, iskpb;
  double c;
  double fm;

  if (kpb) {
    for (t = 0; t < NPRINCIPLE; t++) {
      if (d <= t) {
	kpb[t] = kpb[t-1];
	continue;
      } 
      c = 0.0;
      for (i = 0; i < d; i++) {
	iskpb = 0;
	for (q = 0; q < t; q++) {
	  if (i == kpb[q]) {
	    iskpb = 1;
	    break;
	  }
	}
	if (iskpb) continue;
	fm = fabs(mix[i]);
	if (fm > c) {
	  c = fm;
	  kpb[t] = i;
	}
      }
    }
    k = kpb[0];
  } else {
    c = 0.0;
    for (i = 0; i < d; i++) {
      fm = fabs(mix[i]);
      if (fm > c) {
	c = fm;
	k = i;
      }
    }
  }

  return k;
}

int CompareLevels(cfac_t *cfac, LEVEL *lev1, LEVEL *lev2) {
  STATE *s1, *s2;
  SYMMETRY *sym1, *sym2;
  ORBITAL *orb;
  int i1, i2;
  int p1, p2, j1, j2;

  if (lev1->pj < 0 || lev2->pj < 0) {
    if (lev1->energy > lev2->energy) return 1;
    else if (lev1->energy < lev2->energy) return -1;
    return 0;
  }

  i1 = lev1->pb;
  i2 = lev2->pb;
  sym1 = GetSymmetry(cfac, lev1->pj);
  sym2 = GetSymmetry(cfac, lev2->pj);
  s1 = (STATE *) ArrayGet(&(sym1->states), i1);
  s2 = (STATE *) ArrayGet(&(sym2->states), i2);
  if (s1->kgroup < 0 && s2->kgroup < 0) {
    orb = GetOrbital(cfac, s1->kcfg);
    GetJLFromKappa(orb->kappa, &p1, &j1);
    orb = GetOrbital(cfac, s2->kcfg);
    GetJLFromKappa(orb->kappa, &p2, &j2);
    i1 = p1 - p2;
    if (i1) return i1;
    i1 = j1 - j2;
    if (i1) return i1;
    i1 = (s2->kgroup - s1->kgroup);
    if (i1) return i1;
    DecodePJ(lev1->pj, &p1, &j1);
    DecodePJ(lev2->pj, &p2, &j2);
    i1 = p1 - p2;
    if (i1) return i1;
    return (j1 - j2);
  } else {
    if (lev1->energy > lev2->energy) return 1;
    else if (lev1->energy < lev2->energy) return -1;
    else return 0;
  }
}

int SortLevels(cfac_t *cfac, int start, int n, int m) {
  int i, j, i0, j0;
  LEVEL tmp, *lev1, *lev2, *levp;

  if (m == 0) {
    if (n < 0) n = cfac->n_levels-start;
  } else {
    if (n < 0) n = cfac->n_eblevels-start;
  }
  while (1 < n) {
    i = start;
    j = start + n - 1;
    if (m == 0) {
      lev1 = GetLevel(cfac, i);
      lev2 = GetLevel(cfac, j);
    } else {
      lev1 = GetEBLevel(cfac, i);
      lev2 = GetEBLevel(cfac, j);
    }
    levp = lev2;
    
    while (i < j) {
      while (i < j) {
	if (CompareLevels(cfac, lev1, levp) > 0) break;
	i++;
	if (m == 0) {
	  lev1 = GetLevel(cfac, i);
	} else {
	  lev1 = GetEBLevel(cfac, i);
	}
      }
      while (i < j) {
	if (CompareLevels(cfac, levp, lev2) > 0) break;
	j--;
	if (m == 0) {
	  lev2 = GetLevel(cfac, j);
	} else {
	  lev2 = GetEBLevel(cfac, j);
	}
      }
      if (i < j) {
	memcpy(&tmp, lev1, sizeof(LEVEL));
	memcpy(lev1, lev2, sizeof(LEVEL));
	memcpy(lev2, &tmp, sizeof(LEVEL));
	i++;
	if (m == 0) {
	  lev1 = GetLevel(cfac, i);
	} else {
	  lev1 = GetEBLevel(cfac, i);
	}
      }
    }
    if (lev1 != levp) {
      memcpy(&tmp, lev1, sizeof(LEVEL));
      memcpy(lev1, levp, sizeof(LEVEL));
      memcpy(levp, &tmp, sizeof(LEVEL));
    }

    i0 = i - start;
    j0 = n - i0 - 1;
    if (j0 < i0) {
      if (1 < j0) {
	SortLevels(cfac, i+1, j0, m);
      }
      n = i0;
    } else {
      if (1 < i0) {
	SortLevels(cfac, start, i0, m);
      }
      start = i+1;
      n = j0;
    }
  }
  return 0;
}

int GetLevNumElectrons(cfac_t *cfac, const LEVEL *lev) {
  SYMMETRY *sym;
  STATE *s;
  CONFIG_GROUP *g;
  int nele;
  
  sym = GetSymmetry(cfac, lev->pj);
  s = (STATE *) ArrayGet(&(sym->states), lev->basis[0]);
  if (s->kgroup >= 0) {
    g = GetGroup(cfac, s->kgroup);
    nele = g->n_electrons;
  } else {
    nele = 1+GetNumElectrons(cfac, -(s->kgroup)-1);
  }    

  return nele;
}

int GetNumElectrons(cfac_t *cfac, int k)
{
    LEVEL *lev = GetLevel(cfac, k);
    
    if (!lev) {
        return -1;
    } else {
        return GetLevNumElectrons(cfac, lev);
    }
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


int SaveEBLevels(cfac_t *cfac, char *fn, int m, int n) {
  int n0, k, i, ilev, mlev, nele;
  FILE *f;
  LEVEL *lev;
  F_HEADER fhdr;
  ENF_HEADER enf_hdr;
  ENF_RECORD r;

  n0 = m;
  if (n < 0) n = cfac->n_eblevels - m;
  fhdr.type = DB_ENF;
  strcpy(fhdr.symbol, GetAtomicSymbol(cfac));
  fhdr.atom = GetAtomicNumber(cfac);
  f = OpenFile(fn, &fhdr);

  lev = GetEBLevel(cfac, n0);
  DecodeBasisEB(lev->pb, &ilev, &mlev);
  nele = GetNumElectrons(cfac, ilev);
  enf_hdr.nele = nele;
  enf_hdr.efield = cfac->ef;
  enf_hdr.bfield = cfac->bf;
  enf_hdr.fangle = cfac->eb_angle;
  InitFile(f, &fhdr, &enf_hdr);
  for (k = 0; k < n; k++) {
    i = m + k;
    lev = GetEBLevel(cfac, i);
    r.ilev = i;
    r.energy = lev->energy;
    r.pbasis = lev->pb;
    WriteENFRecord(f, &r);
  }
  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);
  
  return 0;
}  
  
int SaveLevels(cfac_t *cfac, char *fn, int m, int n) {
  STATE *s, *s1;
  SYMMETRY *sym, *sym1;
  CONFIG *cfg, *cfg1;
  SHELL_STATE *csf, *csf1;
  LEVEL *lev, *lev1;
  EN_RECORD r;
  EN_HEADER en_hdr;
  F_HEADER fhdr;
  ECORRECTION *ec;
  LEVEL_ION *gion, gion1;
  ORBITAL *orb;
  double e0, md, md1, a;
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  char nc[LEVEL_NAME_LEN];
  FILE *f;
  int i, k, p, j0;
  int nele, nele0, vnl, ib, dn, ik;
  int si, ms, mst, t, q, nk, n0;

  f = NULL;
  nele0 = -1;
  n0 = m;
  if (n < 0) n = cfac->n_levels - m;
  fhdr.type = DB_EN;
  strcpy(fhdr.symbol, GetAtomicSymbol(cfac));
  fhdr.atom = GetAtomicNumber(cfac);
  f = OpenFile(fn, &fhdr);

  for (k = 0; k < n; k++) {
    i = m + k;
    lev = GetLevel(cfac, i);
    si = lev->pb;
    sym = GetSymmetry(cfac, lev->pj);
    s = (STATE *) ArrayGet(&(sym->states), si);
    if (cfac->ncorrections > 0) {
      for (p = 0; p < cfac->ecorrections->dim; p++) {
	ec = (ECORRECTION *) ArrayGet(cfac->ecorrections, p);
	if (ec->ilev == i) {
	  if (ec->ilev == ec->iref) {
	    e0 = lev->energy;
	  } else {
	    e0 = GetLevel(cfac, ec->iref)->energy;
	  }
	  ec->e = e0 + ec->e - lev->energy;
	  lev->energy += ec->e;
	  ec->s = s;
	  ec->ilev = -(ec->ilev+1);
	  cfac->ncorrections -= 1;
	  break;
	}
      }
    }

    if (s->kgroup > 0) {
      cfg = GetConfig(cfac, s);
      nk = cfg->n_electrons-1;
      if (cfac_get_ion_nlevels(cfac, nk) == 0 || cfg->shells[0].nq > 1) {
	lev->ibase = -1;
      } else {
	csf = cfg->csfs + s->kstate;
	md = 1E30;
	lev->ibase = -1;
	dn = cfg->shells[0].n - cfg->shells[1].n;
	a = 0.0;
	if (dn < MAXDN) {
	  a = 0.0;
	  for (t = 0; t < lev->n_basis; t++) {
	    s1 = ArrayGet(&(sym->states), lev->basis[t]);
	    cfg1 = GetConfig(cfac, s1);
	    if (cfg1->shells[0].n == cfg->shells[0].n &&
		cfg1->shells[0].nq == 1) {
	      a += (lev->mixing[t])*(lev->mixing[t]);
	    }
	  }
	  a = 1.0/a;
	}
	for (ib = 0; ib < NPRINCIPLE; ib++) {
	  for (t = 0; t < cfac_get_ion_nlevels(cfac, nk); t++) {
	    gion = (LEVEL_ION *) ArrayGet(cfac->levels_per_ion+nk, t);
	    for (q = gion->imin; q <= gion->imax; q++) {
	      lev1 = GetLevel(cfac, q);
	      sym1 = GetSymmetry(cfac, lev1->pj);
	      s1 = ArrayGet(&(sym1->states), lev1->basis[lev1->kpb[ib]]);
	      cfg1 = GetConfig(cfac, s1);
	      csf1 = cfg1->csfs + s1->kstate;
	      mst = cfg1->n_shells*sizeof(SHELL_STATE);
	      ms = cfg1->n_shells*sizeof(SHELL);
	      if (cfg->n_shells == cfg1->n_shells+1 &&
		  memcmp(cfg->shells+1, cfg1->shells, ms) == 0 &&
		  memcmp(csf+1, csf1, mst) == 0) {
		if (dn < MAXDN) {
		  md1 = fabs(fabs(a*lev->mixing[lev->kpb[0]]) - 
			     fabs(lev1->mixing[lev1->kpb[ib]]));
		  if (md1 < md) {
		    md = md1;
		    lev->ibase = q;
		  }
		} else {
		  ik = OrbitalIndex(cfac, cfg->shells[0].n, cfg->shells[0].kappa, 0.0);
		  orb = GetOrbital(cfac, ik);
		  a = lev->energy - orb->energy;
		  for (p = 0; p < cfac->ecorrections->dim; p++) {
		    ec = (ECORRECTION *) ArrayGet(cfac->ecorrections, p);
		    if (-(q+1) == ec->ilev) {
		      a += ec->e;
		      break;
		    }
		  }
		  md1 = fabs(lev1->energy - a);
		  if (md1 < md) {
		    md = md1;
		    lev->ibase = q;
		  }
		}
	      }
	    }
	  }
	  if (lev->ibase >= 0) {
	    break;
	  }
	}
      }

      if (lev->ibase >= 0) {
	for (p = 0; p < cfac->ecorrections->dim; p++) {
	  ec = (ECORRECTION *) ArrayGet(cfac->ecorrections, p);
	  if (-(i+1) == ec->ilev) break;
	  if (-(lev->ibase + 1) == ec->ilev && cfg->shells[0].n >= ec->nmin) {
	    lev->energy += ec->e;
	    break;
	  }
	}
      } 
    } else {
      lev->ibase = -(s->kgroup + 1);
    }
 
    DecodePJ(lev->pj, &p, &j0);
    r.ilev = i;
    r.ibase = lev->ibase;
    r.p = p;
    r.j = j0;
    r.energy = lev->energy;

    nele = ConstructLevelName(cfac, name, sname, nc, &vnl, s);
    strncpy(r.name, name, LNAME);
    strncpy(r.sname, sname, LSNAME);
    strncpy(r.ncomplex, nc, LNCOMPLEX);
    r.name[LNAME-1] = '\0';
    r.sname[LSNAME-1] = '\0';
    r.ncomplex[LNCOMPLEX-1] = '\0';
    if (r.p == 0) {
      r.p = vnl;
    } else {
      r.p = -vnl;
    }
    if (nele != nele0) {
      if (nele0 >= 0) {
	DeinitFile(f, &fhdr);
	q = 0;
	nk = nele0;
	t = cfac_get_ion_nlevels(cfac, nk);
	if (t > 0) {
	  gion = (LEVEL_ION *) ArrayGet(cfac->levels_per_ion+nk, t-1);
	  if (gion->imax+1 == n0) {
	    gion->imax = cfac->n_levels-1;
	    q = 1;
	  }
	}
	if (q == 0) {
	  gion1.imin = n0;
	  gion1.imax = i-1;
	  ArrayAppend(cfac->levels_per_ion+nk, &gion1);
	}
      }
      n0 = i;
      nele0 = nele;
      en_hdr.nele = nele;
      InitFile(f, &fhdr, &en_hdr);
    }
    WriteENRecord(f, &r);
  }

  DeinitFile(f, &fhdr);
  CloseFile(f, &fhdr);

  q = 0;
  nk = nele0;
  if (nk >= 0) {
    t = cfac_get_ion_nlevels(cfac, nk);
    if (t > 0) {
      gion = (LEVEL_ION *) ArrayGet(cfac->levels_per_ion+nk, t-1);
      if (gion->imax+1 == n0) {
	gion->imax = cfac->n_levels-1;
	q = 1;
      }
    }
    if (q == 0 && cfac->n_levels > n0) {
      gion1.imin = n0;
      gion1.imax = cfac->n_levels-1;
      ArrayAppend(cfac->levels_per_ion+nk, &gion1);
    }
  }
  
  return 0;
}

int ConstructLevelName(cfac_t *cfac, char *name, char *sname, char *nc, 
		       int *vnl, STATE *basis) {
  int n, nq, kl, j;
  int nele, i, len;
  char symbol[20];
  char jsym;
  char ashell[16];
  CONFIG *c;
  SHELL_STATE *s = NULL;
  ORBITAL *orb;
  LEVEL *lev;
  SYMMETRY *sym;
  int si;
  int n0, kl0, nq0;

  symbol[0] = '\0';
  if (basis->kgroup < 0) {
    i = basis->kgroup;
    i = -(i + 1);
    if (basis->kcfg < 0) {
      lev = GetLevel(cfac, i);
      si = lev->pb;
      sym = GetSymmetry(cfac, lev->pj);
      basis = ArrayGet(&(sym->states), si);
      nele = ConstructLevelName(cfac, name, sname, nc, vnl, basis);
      return nele;
    } else {
      orb = GetOrbital(cfac, basis->kcfg);
      GetJLFromKappa(orb->kappa, &j, &kl);
      if (vnl) {
	*vnl = (kl/2) + 100*(orb->n);
      }
      if (name) {
	if (j < kl) jsym = '-';
	else jsym = '+';
	
	kl /= 2;
	SpecSymbol(symbol, kl);
	sprintf(name, "%5d + %d%s%c1(%d)%d ", 
		i, orb->n, symbol, jsym, j, basis->kstate);
      }
      lev = GetLevel(cfac, i);
      si = lev->pb;
      sym = GetSymmetry(cfac, lev->pj);
      basis = (STATE *) ArrayGet(&(sym->states), si);
      if (sname || nc) {
	nele = ConstructLevelName(cfac, NULL, sname, nc, NULL, basis);
	if (nc) {
	  if (nele == 0) {
	    nc[0] = '\0';
	  }
	  sprintf(ashell, "%1d*1", orb->n);
	  strcat(nc, ashell);
	}
      } else {
	nele = ConstructLevelName(cfac, NULL, NULL, NULL, NULL, basis);
      }
      return nele+1;
    }
  }

  c = GetConfig(cfac, basis);
  nele = c->n_electrons;
  if (!name && !sname && !nc) return nele;

  if (c->n_csfs > 0) {
    s = c->csfs + basis->kstate;
  }
  len = 0;
  if (name) name[0] = '\0';
  if (sname) sname[0] = '\0';
  if (nc) nc[0] = '\0';
  n0 = 0;
  kl0= -1;
  nq0 = 0;
  for (i = c->n_shells-1; i >= 0; i--) {
    UnpackShell(c->shells+i, &n, &kl, &j, &nq);
    if (j < kl) jsym = '-';
    else jsym = '+';
    kl = kl/2;
    if (name) {
      if (((nq < j+1) && nq > 0) || (i == 0 && name[0] == '\0')) {
	SpecSymbol(symbol, kl);
	if (c->n_csfs > 0) {
	  sprintf(ashell, "%1d%s%c%1d(%1d)%1d ", 
		  n, symbol, jsym, nq, s[i].shellJ, s[i].totalJ); 
	} else {
	  sprintf(ashell, "%1d%s%c%1d ", n, symbol, jsym, nq);
	}
	len += strlen(ashell);
	if (len >= LEVEL_NAME_LEN) return -1;
	strcat(name, ashell);
      }
    }
    if (sname) {
      if (n == n0 && kl == kl0) {
	nq0 += nq;
      } else {
	if (nq0 > 0 && nq0 < 2*(2*kl0+1)) {
	  SpecSymbol(symbol, kl0);
	  sprintf(ashell, "%1d%s%1d ", n0, symbol, nq0);
	  strcat(sname, ashell);
	}
	n0 = n;
	kl0 = kl;
	nq0 = nq;
      }
    }
  }
  
  if (sname && n0 > 0) {
    if ((nq0 > 0 && nq0 < 2*(2*kl0+1)) || sname[0] == '\0') {
      SpecSymbol(symbol, kl0);
      sprintf(ashell, "%1d%s%1d ", n0, symbol, nq0);
      strcat(sname, ashell);
    }
  }

  if (nc) {
    n0 = 0;
    nq0 = 0;
    for (i = c->n_shells-1; i >= 0; i--) {
      UnpackShell(c->shells+i, &n, &kl, &j, &nq);
      if (n == n0) {
	nq0 += nq;
      } else {
	if (nq0 > 0) {
	  sprintf(ashell, "%1d*%1d ", n0, nq0);
	  strcat(nc, ashell);
	}
	n0 = n;
	nq0 = nq;
      }
    }
    if (n0 > 0 && (nq0 > 0 || nc[0] == '\0')) {
      sprintf(ashell, "%1d*%1d ", n0, nq0);
      strcat(nc, ashell);
    }
  }

  if (vnl) {
    UnpackShell(c->shells, &n, &kl, &j, &nq);
    *vnl = (kl/2) + 100*n;
  }
      
  return nele;
}
    
int GetBasisTable(cfac_t *cfac, char *fn, int m) {
  FILE *f;
  int i, p, j, k, si, nsym;
  char nc[LEVEL_NAME_LEN];
  char name[LEVEL_NAME_LEN];
  char sname[LEVEL_NAME_LEN];
  ARRAY *st;
  STATE *s;
  LEVEL *lev;
  SYMMETRY *sym;

  f = fopen(fn, "w");
  if (!f) return -1;

  if (m > 0) {
    for (i = 0; i < cfac->n_eblevels; i++) {
      lev = GetEBLevel(cfac, i);
      for (k = 0; k < lev->n_basis; k++) {
	DecodeBasisEB(lev->basis[k], &p, &j);
	fprintf(f, "%6d %6d %3d %15.8E\n", 
		i, p, j, lev->mixing[k]);
      }
    }    
  } else {
    nsym = MAX_SYMMETRIES;
    fprintf(f, "============Basis Table===========================\n");
    for (i = 0; i < nsym; i++) {
      sym = GetSymmetry(cfac, i);
      DecodePJ(i, &p, &j);
      st = &(sym->states);
      if (sym->n_states <= 0) continue;
      for (k = 0; k < sym->n_states; k++) {
	s = (STATE *) ArrayGet(st, k);
	ConstructLevelName(cfac, name, sname, nc, NULL, s);
	fprintf(f, "%6d   %2d %2d   %5d %3d %5d %5d   %-20s %-20s %-20s\n",
		i, p, j, k, s->kgroup, s->kcfg, s->kstate, nc, sname, name);
      }
      fprintf(f, "\n");
    }
    
    fprintf(f, "============Mixing Coefficients===================\n");
    for (i = 0; i < cfac->n_levels; i++) {
      lev = GetLevel(cfac, i);
      sym = GetSymmetry(cfac, lev->pj);
      DecodePJ(lev->pj, &p, &j);
      for (k = 0; k < lev->n_basis; k++) {
	si = lev->basis[k];
	s = (STATE *) ArrayGet(&(sym->states), si);
	fprintf(f, "%6d   %2d %2d   %5d %3d %5d %5d   %15.8E\n", 
		i, p, j, si, s->kgroup, s->kcfg, s->kstate, lev->mixing[k]);
      }
      fprintf(f, "\n");
    }
  }
  fclose(f);
  return 0;
}

void StructureEB(cfac_t *cfac, char *fn, int n, int *ilev) {
  int k;
  
  ConstructHamiltonEB(cfac, n, ilev);

  DiagonalizeHamilton(cfac);
  
  k = cfac->n_eblevels;
  AddToLevels(cfac, 0, NULL);
  SortLevels(cfac, k, -1, 1);
  SaveEBLevels(cfac, fn, k, -1);
}

int AngularZMixStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2) {
  int kg1, kg2, kc1, kc2;
  int ns, n, p, q, nz, iz, iz1, iz2;
  int ns1, ns2, *pnz;
  int ks1, ks2, i1, i2, i2m;
  int n_shells, *k, nkk;
  double *r;
  int phase, im;
  int orb0, orb1;
  CONFIG *c1, *c2;
  STATE *s1, *s2;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  ANGULAR_ZMIX **a, *ang;
  int kmax = GetMaxRank(cfac);
  
  iz = ih1*MAX_HAMS + ih2;
  *ad = &(cfac->angz_array[iz]);
  ns = (*ad)->ns;
  if (ns < 0) {
    return -1;
  }
  
  if (ns > 0) {
    return ns;
  }

  ns1 = cfac->hams[ih1].nbasis;
  ns2 = cfac->hams[ih2].nbasis;
  (*ad)->ns = ns1*ns2;
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  iz = 0;
  iz1 = 0;
  iz2 = 0;
  a = (ANGULAR_ZMIX **) (*ad)->angz;
  pnz = (*ad)->nz;

  for (i1 = 0; i1 < ns1; i1++) {
    s1 = cfac->hams[ih1].basis[i1];
    kg1 = s1->kgroup;
    kc1 = s1->kcfg;
    ks1 = s1->kstate;    
    c1 = GetConfigFromGroup(cfac, kg1, kc1);
    if (ih1 == ih2) {
      i2m = i1;
    } else {
      i2m = 0;
    }
    for (i2 = i2m; i2 < ns2; i2++) {
      s2 = cfac->hams[ih2].basis[i2];
      kg2 = s2->kgroup;
      kc2 = s2->kcfg;
      ks2 = s2->kstate;          
      c2 = GetConfigFromGroup(cfac, kg2, kc2);
      if (ih1 == ih2) {
	iz1 = i1*ns2 + i2;
	iz2 = i2*ns1 + i1;
      }
      if (abs(c1->n_shells - c2->n_shells) > 1) {
	if (ih1 != ih2) {
	  a[iz] = NULL;
	  pnz[iz] = 0;
	  iz++;
	} else {
	  a[iz1] = NULL;
	  pnz[iz1] = 0;
	  a[iz2] = NULL;
	  pnz[iz2] = 0;
	}
	continue;
      }

      idatum = NULL; 
      n_shells = GetInteract(cfac, &idatum, &sbra, &sket, 
			     kg1, kg2, kc1, kc2, 
			     ks1, ks2, 0);
      n = 0;
      ang = NULL;
      if (n_shells <= 0) goto OUT;
      memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
      phase = idatum->phase;
      bra = idatum->bra;
      if (s[2].index >= 0 && s[0].index >= 0) {
	free(sbra);
	free(sket);
	goto OUT;
      }

      nz = ANGZ_BLOCK;
      ang = malloc(sizeof(ANGULAR_ZMIX)*nz);
      if (!ang) {
	printf("failed allocating memory for ang %d %d\n", nz, ns);
	exit(1);
      }
      if (s[0].index >= 0) {
	nkk = AngularZ(&r, &k, 0, n_shells, sbra, sket, s, s+1, kmax);
	if (nkk > 0) {
	  orb0 = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
	  orb1 = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
	  for (p = 0; p < nkk; p++) {
	    if (IsOdd(phase)) r[p] = -r[p];
	    im = AddToAngularZMix(&n, &nz, &ang, k[p], orb0, orb1, r[p]);
	  }
	  free(r);
	  free(k);
	}
      } else {
	for (q = 0; q < n_shells; q++) {	    
	  p = ShellToInt(bra[q].n, bra[q].kappa);
	  if (ih1 != ih2 || i1 != i2) {
	    if (IsClosedShell(cfac, ih1, p) &&
                IsClosedShell(cfac, ih2, p)) continue;
	  }
	  s[0].index = n_shells - q - 1;      
	  s[1].index = s[0].index;
	  s[0].n = bra[q].n;
	  s[1].n = s[0].n;
	  s[0].kappa = bra[q].kappa;
	  s[1].kappa = s[0].kappa;
	  s[0].j = GetJ(bra+q);
	  s[1].j = s[0].j;
	  s[0].kl = GetL(bra+q);
	  s[1].kl = s[0].kl;
	  s[0].nq_bra = GetNq(bra+q);
	  s[0].nq_ket = s[0].nq_bra;
	  s[1].nq_bra = s[0].nq_bra;
	  s[1].nq_ket = s[1].nq_bra;
	  nkk = AngularZ(&r, &k, 0, n_shells, sbra, sket, s, s+1, kmax);
	  if (nkk > 0) {
	    orb0 = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
	    orb1 = orb0;
	    for (p = 0; p < nkk; p++) {
	      if (fabs(r[p]) < EPS30) continue;
	      if (IsOdd(phase)) r[p] = -r[p];
	      im = AddToAngularZMix(&n, &nz, &ang, k[p], orb0, orb1, r[p]);
	    }	    
	    free(r);
	    free(k);
	  }
	}
      }
      PackAngularZMix(&n, &ang, nz);
      free(sbra);
      free(sket);
    OUT:
      if (ih1 != ih2) {
	a[iz] = ang;
	pnz[iz] = n;
	iz++;
      } else {
	a[iz1] = ang;
	pnz[iz1] = n;
	if (iz2 != iz1) {
	  if (n > 0) {
	    ang = malloc(sizeof(ANGULAR_ZMIX)*n);
	    memcpy(ang, a[iz1], sizeof(ANGULAR_ZMIX)*n);
	    AngZSwapBraKet(cfac, n, ang, 0);
	  }
	  a[iz2] = ang;
	  pnz[iz2] = n;
	}
      }
    }
  }

  return (*ad)->ns;
}

int AngZSwapBraKet(cfac_t *cfac, int nz, ANGULAR_ZMIX *ang, int p) {
  int i;
  int k0, k1;
  for (i = 0; i < nz; i++) {
    k0 = ang[i].k0;
    k1 = ang[i].k1;
    ang[i].k0 = k1;
    ang[i].k1 = k0;
    k0 = GetOrbital(cfac, k0)->kappa;
    k1 = GetOrbital(cfac, k1)->kappa;
    k0 = GetJFromKappa(k0);
    k1 = GetJFromKappa(k1);
    if (IsOdd(abs(k1-k0+p)/2)) ang[i].coeff = - ang[i].coeff;
  }
  return 0;
}
    
int AngularZFreeBoundStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2) {
  int kg1, kg2;
  int kc1, kc2;
  int n_shells, i1, i2;
  int phase;
  int j1, j2, kb;
  int *k, k0, nkk, kmax;
  int jf, jp, tf;
  double *r, r0;
  int ns1, ns2, *pnz, iz;
  int ns, ks1, ks2, n;
  STATE *s1, *s2;
  INTERACT_DATUM *idatum;
  INTERACT_SHELL s[4];
  SHELL_STATE *sbra, *sket;
  CONFIG *c1, *c2;
  ANGULAR_ZFB *ang, **a;
  
  iz = ih1*MAX_HAMS + ih2;
  *ad = &(cfac->angz_array[iz]);
  ns = (*ad)->ns;

  if (ns < 0) {
    return -1;
  }
  if (ns > 0) {
    return ns;
  }

  ns1 = cfac->hams[ih1].nbasis;
  ns2 = cfac->hams[ih2].nbasis;
  (*ad)->ns = ns1 * ns2;
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  
  kmax = GetMaxRank(cfac);

  iz = 0;
  a = (ANGULAR_ZFB **) (*ad)->angz;
  pnz = (*ad)->nz;
    
  for (i1 = 0; i1 < ns1; i1++) {
    s1 = cfac->hams[ih1].basis[i1];
    kg1 = s1->kgroup;
    kc1 = s1->kcfg;
    ks1 = s1->kstate;    
    c1 = GetConfigFromGroup(cfac, kg1, kc1);
    j1 = (c1->csfs[ks1]).totalJ;
    for (i2 = 0; i2 < ns2; i2++) {
      s2 = cfac->hams[ih2].basis[i2];
      kg2 = s2->kgroup;
      kc2 = s2->kcfg;
      ks2 = s2->kstate;          
      c2 = GetConfigFromGroup(cfac, kg2, kc2);
          
      if (abs(c1->n_shells+1 - c2->n_shells) > 1) {
	a[iz] = NULL;
	pnz[iz] = 0;
	iz++;
	continue;
      }      
      
      j2 = (c2->csfs[ks2]).totalJ;
      idatum = NULL;
      n_shells = GetInteract(cfac, &idatum, &sbra, &sket,
			     kg1, kg2, kc1, kc2, 
			     ks1, ks2, 1);
      n = 0;
      ang = NULL;
      if (n_shells <= 0) goto END;
      memcpy(s, idatum->s, sizeof(INTERACT_SHELL)*4);
      phase = idatum->phase;
      if (s[0].index < 0 || (s[0].index >= 0 && s[2].index >= 0)) {
	free(sbra);
	free(sket);
	goto END;
      }
      
      tf = 0;
      for (k0 = 0; k0 <= kmax; k0 += 2) {
	for (jf = abs(k0-s[1].j); jf <= k0+s[1].j; jf += 2) {
	  for (jp = abs(jf-j1); jp <= jf+j1; jp += 2) {
	    if (Triangle(jp, j2, k0) && Triangle(j1, j2, s[1].j)) {
	      tf = 1;
	      goto TF;
	    }
	  }
	}
      }	  
    TF:	  
      if (tf == 1) {
	s[0].j = jf;
	s[0].kl = jf + 1;
	s[0].kappa = GetKappaFromJL(s[0].j, s[0].kl);
	sbra[0].shellJ = jf;
	sbra[0].totalJ = jp; 
	k = &k0;
	r = &r0;
	nkk = AngularZ(&r, &k, 1, n_shells, sbra, sket, s, s+1, kmax);
	if (fabs(*r) < EPS30) goto END;
	if (IsOdd(phase+(jp+j2-k0)/2)) *r = -(*r);
	*r /= sqrt(jp+1.0)*W6j(j1, jf, jp, k0, j2, s[1].j);
	kb = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
	ang = (ANGULAR_ZFB *) malloc(sizeof(ANGULAR_ZFB));
	ang->kb = kb;
	ang->coeff = *r;
	n = 1;
      }
      free(sbra);
      free(sket);
    END:
      a[iz] = ang;
      pnz[iz] = n;
      iz++;
    }
  }
  
  return (*ad)->ns;
}

int AngularZxZFreeBoundStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2) {
  int kg1, kg2, kc1, kc2;
  int n_shells, i1, i2;
  int phase;
  int j1, j2, i, n, nz, p;
  int jmin, jmax, jf;
  int ns1, ns2, *pnz, iz;
  int ns, ks1, ks2;  
  INTERACT_SHELL s[4];
  SHELL *bra;
  SHELL_STATE *sbra, *sket;
  INTERACT_DATUM *idatum;
  CONFIG *c1, *c2;
  STATE *s1, *s2;
  ANGULAR_ZxZMIX **a, *ang;

  iz = ih1 * MAX_HAMS + ih2;
  *ad = &(cfac->angzxz_array[iz]);
  ns = (*ad)->ns;

  if (ns < 0) { 
    return -1;
  }
  
  if (ns > 0) {
    return ns;
  }
  
  ns1 = cfac->hams[ih1].nbasis;
  ns2 = cfac->hams[ih2].nbasis;
  (*ad)->ns = ns1*ns2;
  ns = (*ad)->ns;
  (*ad)->angz = malloc(sizeof(ANGULAR_ZxZMIX *)*ns);
  (*ad)->nz = (int *) malloc(sizeof(int)*ns);
  
  iz = 0;
  a = (ANGULAR_ZxZMIX **) (*ad)->angz;
  pnz = (*ad)->nz;  


  for (i1 = 0; i1 < ns1; i1++) {
    s1 = cfac->hams[ih1].basis[i1];
    kg1 = s1->kgroup;
    kc1 = s1->kcfg;
    ks1 = s1->kstate;    
    c1 = GetConfigFromGroup(cfac, kg1, kc1);
    j1 = (c1->csfs[ks1]).totalJ;
    for (i2 = 0; i2 < ns2; i2++) {
      s2 = cfac->hams[ih2].basis[i2];
      kg2 = s2->kgroup;
      kc2 = s2->kcfg;
      ks2 = s2->kstate;          
      c2 = GetConfigFromGroup(cfac, kg2, kc2);
      
      if (abs(c1->n_shells+1 - c2->n_shells) > 2) {
	a[iz] = NULL;
	pnz[iz] = 0;
	iz++;
	continue;
      }
   
      j2 = (c2->csfs[ks2]).totalJ;
      idatum = NULL;
      n_shells = GetInteract(cfac, &idatum, &sbra, &sket,
			     kg1, kg2, kc1, kc2, 
			     ks1, ks2, 1);
      n = 0;
      ang = NULL;
      if (n_shells <= 0) goto END;
      
      nz = ANGZxZ_BLOCK;
      ang = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
      phase = idatum->phase;
      bra = idatum->bra;
      sbra[0].totalJ = j2;
      jmin = abs(j2 - j1);
      jmax = j1 + j2;
      for (jf = jmin; jf <= jmax; jf += 2) {
	memcpy(s, idatum->s+2, sizeof(INTERACT_SHELL)*2);
	memcpy(s+2, idatum->s, sizeof(INTERACT_SHELL)*2);
	s[2].j = jf;
	s[2].kl = jf+1;
	s[2].kappa = GetKappaFromJL(s[2].j, s[2].kl);
	sbra[0].shellJ = s[2].j;
	
	if (s[0].index >= 0) {
	  AddToAngularZxZ(cfac, &n, &nz, &ang, n_shells, phase, 
			  sbra, sket, s, 1);
	} else {
	  for (i = 0; i < n_shells; i++) {
	    p = ShellToInt(bra[i].n, bra[i].kappa);
	    s[0].index = n_shells - i - 1;
	    if (s[0].index == s[2].index) continue;
	    if (s[0].index == s[3].index && s[3].nq_ket < 2) continue;
	    s[1].index = s[0].index;
	    s[0].n = bra[i].n;
	    s[1].n = s[0].n;
	    s[0].kappa = bra[i].kappa;
	    s[1].kappa = s[0].kappa;
	    s[0].j = GetJ(bra+i);
	    s[1].j = s[0].j;
	    s[0].kl = GetL(bra+i);
	    s[1].kl = s[0].kl;
	    s[0].nq_bra = GetNq(bra+i);
	    if (s[0].index == s[2].index) {
	      s[0].nq_ket = s[0].nq_bra - 1;
	    } else if (s[0].index == s[3].index) {
	      s[0].nq_ket = s[0].nq_bra + 1;
	    } else {
	      s[0].nq_ket = s[0].nq_bra;
	    }
	    s[1].nq_bra = s[0].nq_bra;
	    s[1].nq_ket = s[0].nq_ket;
	    AddToAngularZxZ(cfac, &n, &nz, &ang, n_shells, phase, 
			    sbra, sket, s, 1);
	  }
	}
      }
      
      PackAngularZxZMix(&n, &ang, nz);
      free(sbra);
      free(sket);
      
    END:	  
      a[iz] = ang;
      pnz[iz] = n;
      iz++;
    }
  }

  return (*ad)->ns;
}

int PrepAngular(cfac_t *cfac, int n1, int *is1, int n2, int *is2) {
  int i1, i2, ih1, ih2, ns1, ns2, ne1, ne2;
  int iz, is, i, nz, ns = 0;
  SYMMETRY *sym1, *sym2;
  STATE *s1, *s2;
  LEVEL *lev1, *lev2;
  ANGZ_DATUM *ad;
  
  if (cfac->angmz_array == NULL) {
    cfac->angmz_array = malloc(sizeof(ANGZ_DATUM)*MAX_HAMS2);
    if (!cfac->angmz_array) {
      printf("cannot allocate memory for cfac->angmz_array %d\n", MAX_HAMS2);
    }
    for (i = 0; i < MAX_HAMS2; i++) {
      cfac->angmz_array[i].ns = 0;
      cfac->angmz_array[i].mk = NULL;
    }
  }

  if (n2 == 0) {
    n2 = n1;
    is2 = is1;
  }
  
  for (i1 = 0; i1 < n1; i1++) {
    lev1 = GetLevel(cfac, is1[i1]);
    ih1 = lev1->iham;
    sym1 = GetSymmetry(cfac, lev1->pj);
    s1 = ArrayGet(&(sym1->states), lev1->pb);
    if (s1->kgroup < 0) continue;
    ne1 = GetGroup(cfac, s1->kgroup)->n_electrons;
    ns1 = cfac->hams[ih1].nlevs;
    for (i2 = 0; i2 < n2; i2++) {
      lev2 = GetLevel(cfac, is2[i2]);
      ih2 = lev2->iham;
      sym2 = GetSymmetry(cfac, lev2->pj);
      s2 = ArrayGet(&(sym2->states), lev2->pb);
      if (s2->kgroup < 0) continue;
      ne2 = GetGroup(cfac, s2->kgroup)->n_electrons;
      if (abs(ne2-ne1) > 1) continue;
      ns2 = cfac->hams[ih2].nlevs;
      ns = ns1*ns2;
      if (ne1 == ne2) {
	if (ih1 > ih2) {
	  iz = ih2 * MAX_HAMS + ih1;
	  is = lev2->ilev * cfac->hams[ih1].nlevs + lev1->ilev;
	} else {
	  iz = ih1 * MAX_HAMS + ih2;
	  is = lev1->ilev * cfac->hams[ih2].nlevs + lev2->ilev;
	}
      } else {
	if (ne1 > ne2) {
	  iz = ih2 * MAX_HAMS + ih1;
	  is = lev2->ilev * cfac->hams[ih1].nlevs + lev1->ilev;
	} else {
	  iz = ih1 * MAX_HAMS + ih2;
	  is = lev1->ilev * cfac->hams[ih2].nlevs + lev2->ilev;
	}
      }
      ad = &(cfac->angmz_array[iz]);
      if (ad->ns == 0) {
	if (ne1 == ne2) {
	  ad->angz = malloc(sizeof(ANGULAR_ZMIX *)*ns);
	} else {
	  ad->angz = malloc(sizeof(ANGULAR_ZFB *)*ns);
	}
	ad->nz = malloc(sizeof(int)*ns);
	for (i = 0; i < ns; i++) (ad->nz)[i] = 0;
	ad->ns = ns;
      }
      nz = (ad->nz)[is];
      if (nz != 0) continue;
      if (ne1 == ne2) {
	if (ih1 > ih2) {
	  nz = AngularZMix(cfac, (ANGULAR_ZMIX **)(&((ad->angz)[is])), 
			   is2[i2], is1[i1], -1, -1);
	} else {
	  nz = AngularZMix(cfac, (ANGULAR_ZMIX **)(&((ad->angz)[is])), 
			   is1[i1], is2[i2], -1, -1);
	}
      } else {
	if (ne1 > ne2) {
	  nz = AngularZFreeBound(cfac, (ANGULAR_ZFB **)(&((ad->angz)[is])), 
				 is2[i2], is1[i1]);
	} else {
	  nz = AngularZFreeBound(cfac, (ANGULAR_ZFB **)(&((ad->angz)[is])), 
				 is1[i1], is2[i2]);
	}	
      }
      if (nz == 0) nz = -1;
      (ad->nz)[is] = nz;
    }
  }

  return ns;
}

int AngularZFreeBound(cfac_t *cfac, ANGULAR_ZFB **ang, int lower, int upper) {
  int i, j, m; 
  int nz, n;
  double r0;
  int ih1, ih2;
  int ns, isz0, isz;
  STATE *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2;
  ANGZ_DATUM *ad;
  ANGULAR_ZFB *ang_sub;
  double mix1, mix2, sqrt_j2;
  int kg, jf, kb, ia, j1, j2;
  
  lev1 = GetLevel(cfac, lower);
  lev2 = GetLevel(cfac, upper);

  if (cfac->angmz_array) {
    ih1 = lev1->iham;
    ih2 = lev2->iham;
    if (ih1 >= 0 && ih2 >= 0) {
      isz = ih1 * MAX_HAMS + ih2;
      ad = &(cfac->angmz_array[isz]);
      nz = 0;
      if (ad->ns > 0) {
	isz0 = lev1->ilev * cfac->hams[ih2].nlevs + lev2->ilev;
	nz = (ad->nz)[isz0];
	if (nz > 0) {
	  *ang = malloc(sizeof(ANGULAR_ZFB)*nz);
	  memcpy(*ang, (ad->angz)[isz0], sizeof(ANGULAR_ZFB)*nz);
	}
      }
      if (nz != 0) return nz;
    }
  }

  sym1 = GetSymmetry(cfac, lev1->pj);
  sym2 = GetSymmetry(cfac, lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  if (sup->kgroup < 0) {
    sqrt_j2 = sqrt(j2 + 1.0);
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZFB)*nz);
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < cfac->angz_cut) {
	break;
      }
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
      kg = sup->kgroup;
      kg = -kg-1;
      if (kg == lower) {
	kb = sup->kcfg;
	jf = GetOrbital(cfac, kb)->kappa;
	jf = GetJFromKappa(jf);
	r0 = mix2*sqrt_j2;
	if (fabs(r0) < cfac->angz_cut) continue;
	if (IsEven((j2+jf-j1)/2)) r0 = -r0;
	ia = AddToAngularZFB(&n, &nz, ang, kb, r0);
      }
    }    
  } else {
    n = 0;
    nz = ANGZ_BLOCK;
    (*ang) = malloc(sizeof(ANGULAR_ZFB)*nz);
    ns = AngularZFreeBoundStates(cfac, &ad, lev1->iham, lev2->iham);
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < cfac->angz_cut) continue;
      ih1 = lev1->ibasis[i];
      isz0 = ih1 * cfac->hams[lev2->iham].nbasis;
      for (j = 0; j < lev2->n_basis; j++) {
	mix2 = lev2->mixing[j];
	if (fabs(mix2) < cfac->angz_cut) continue;
	r0 = mix1*mix2;
	if (fabs(r0) < cfac->angz_cut) continue;
	ih2 = lev2->ibasis[j];
	isz = isz0 + ih2;
	m = (ad->nz)[isz];
	if (m == 1) {
	  ang_sub = (ad->angz)[isz];
	  kb = ang_sub->kb;
	  r0 *= ang_sub->coeff;
	  if (fabs(r0) < cfac->angz_cut) continue;
	  ia = AddToAngularZFB(&n, &nz, ang, kb, r0);
	}
      }
    }
  }
  
  PackAngularZFB(&n, ang, nz);

  return n;
}

int GetBaseJ(cfac_t *cfac, STATE *s) {
  int ih;
  if (s->kgroup >= 0) return -1;
  ih = s->kgroup;
  ih = -ih-1;
  ih = GetLevel(cfac, ih)->pj;
  DecodePJ(ih, NULL, &ih);
  return ih;
}

int AngularZMix(cfac_t *cfac,
  ANGULAR_ZMIX **ang, int lower, int upper, int mink, int maxk) {
  int i, j, j1, j2, jb1, jb2;
  int kg1, kg2, kc1, kc2;
  int ih1, ih2, isz0, isz;
  int jlow, jup, kb1, kb2;
  int nz, n, ns, im;
  double r0;
  int ik, kmin, kmax, m, nmax;
  int nz_sub, nfb;
  STATE *slow, *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2, *lev;
  ANGZ_DATUM *ad;
  ANGULAR_ZMIX *ang_sub;
  ANGULAR_ZFB *afb;
  double mix1, mix2, sqrt_j12, a;
  int ignore_ryd = 0;

  lev1 = GetLevel(cfac, lower);
  lev2 = GetLevel(cfac, upper);
  if (cfac->angmz_array) {
    ih1 = lev1->iham;
    ih2 = lev2->iham;
    if (ih1 >= 0 && ih2 >= 0) {
      if (ih1 > ih2) {
	isz = ih2 * MAX_HAMS + ih1;
      } else {
	isz = ih1 * MAX_HAMS + ih2;
      }
      ad = &(cfac->angmz_array[isz]);
      nz = 0;
      if (ad->ns > 0) {
	if (ih1 > ih2) {
	  isz0 = lev2->ilev * cfac->hams[ih1].nlevs + lev1->ilev;
	} else {
	  isz0 = lev1->ilev * cfac->hams[ih2].nlevs + lev2->ilev;
	}
	nz = (ad->nz)[isz0];
	if (nz > 0) {
	  *ang = malloc(sizeof(ANGULAR_ZMIX)*nz);
	  memcpy(*ang, (ad->angz)[isz0], sizeof(ANGULAR_ZMIX)*nz);
	}
      }  
      if (nz != 0) {
	if (nz > 0 && ih1 > ih2) {
	  sym1 = GetSymmetry(cfac, lev1->pj);
	  sym2 = GetSymmetry(cfac, lev2->pj);
	  j1 = lev1->pj;
	  j2 = lev2->pj;
	  DecodePJ(j1, NULL, &j1);
	  DecodePJ(j2, NULL, &j2);
	  AngZSwapBraKet(cfac, nz, *ang, j1-j2);
	}
	return nz;
      }
    }
  }

  sym1 = GetSymmetry(cfac, lev1->pj);
  sym2 = GetSymmetry(cfac, lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  sqrt_j12 = sqrt((j1+1.0)*(j2+1.0));
  kmin = abs(j1-j2);
  kmax = j1+j2;
  if (IsOdd(kmin)) kmin++;
  if (mink >= 0) kmin = Max(kmin, mink);
  if (maxk >= 0) kmax = Min(kmax, maxk);
  if (kmax < kmin) {
    return 0;
  }

  slow = (STATE *) ArrayGet(&(sym1->states), lev1->pb);
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  kg1 = slow->kgroup;
  kg2 = sup->kgroup;
  kc1 = slow->kcfg;
  kc2 = sup->kcfg;

  if (kg1 < 0) {
    kb1 = slow->kcfg;
    n = GetOrbital(cfac, kb1)->n;
    if (cfac->angz_maxn > 0 && cfac->angz_maxn < n) ignore_ryd = 1;
    nmax = GetNMax(cfac->potential);
    if (n >= nmax && kg2 < 0) {
      kb2 = sup->kcfg;
      m = GetOrbital(cfac, kb2)->n;
      if (m == n) {
	ignore_ryd = 1;
      }
    }
  } else if (kg2 < 0) {
    kb2 = sup->kcfg;
    n = GetOrbital(cfac, kb2)->n;
    if (cfac->angz_maxn > 0 && cfac->angz_maxn < n) ignore_ryd = 1;
  }
  
  n = 0;
  nz = ANGZ_BLOCK;
  (*ang) = malloc(sizeof(ANGULAR_ZMIX)*nz);  
  if (kg1 < 0 && kg2 < 0) {      
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < cfac->angz_cut) continue;
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
      jlow = GetBaseJ(cfac, slow);
      kg1 = slow->kgroup;
      kg1 = -kg1-1;
      kb1 = slow->kcfg;
      jb1 = GetOrbital(cfac, kb1)->kappa;
      jb1 = GetJFromKappa(jb1);
      for (j = 0; j < lev2->n_basis; j++) {
	mix2 = lev2->mixing[j];
	if (fabs(mix2) < cfac->angz_cut) continue;
	a = mix1*mix2;
	if (fabs(a) < cfac->angz_cut) continue;
	sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
	jup = GetBaseJ(cfac, sup);
	kg2 = sup->kgroup;
	kg2 = -kg2-1;
	kb2 = sup->kcfg;
	jb2 = GetOrbital(cfac, kb2)->kappa;
	jb2 = GetJFromKappa(jb2);

	if (ignore_ryd && kb1 != kb2) {
	  continue;
	}

	if (kg1 == kg2) {
	  for (ik = kmin; ik <= kmax; ik += 2) {
	    r0 = W6j(jlow, jb1, j1, ik, j2, jb2);
	    if (fabs(r0) < EPS30) continue;
	    r0 *= a*sqrt_j12;
	    if (fabs(r0) < cfac->angz_cut) continue;
	    if (IsEven((j1+jb2-jlow-ik)/2+j2)) r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ik, kb1, kb2, r0);
	  }
	}
	if (kb1 == kb2){	  
	  nz_sub = AngularZMix(cfac, &ang_sub, kg1, kg2, kmin, kmax);
	  if (nz_sub <= 0) {
	    continue;
	  }
	  for (m = 0; m < nz_sub; m++) {
	    r0 = W6j(jlow, jb1, j1, j2, ang_sub[m].k, jup);
	    if (fabs(r0) < EPS30) continue;
	    r0 *= a*ang_sub[m].coeff;
	    if (fabs(r0) < cfac->angz_cut) continue;
	    r0 *= sqrt_j12;
	    if (IsOdd((jlow+jb1+j2+ang_sub[m].k)/2)) r0 = -r0;
	    im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				  ang_sub[m].k0, ang_sub[m].k1, r0);
	  }
	  if (nz_sub > 0) free(ang_sub);
	}
      }
    }
    PackAngularZMix(&n, ang, nz);
  } else if (kg1 < 0 && !ignore_ryd) {        
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < cfac->angz_cut) continue;
      slow = (STATE *) ArrayGet(&(sym1->states), lev1->basis[i]);
      kb1 = slow->kcfg;
      jb1 = GetOrbital(cfac, kb1)->kappa;
      jb1 = GetJFromKappa(jb1);
      jlow = GetBaseJ(cfac, slow);
      kg1 = slow->kgroup;
      kg1 = -kg1-1;
      nfb = AngularZFreeBound(cfac, &afb, kg1, upper);
      for (m = 0; m < nfb; m++) {
	jb2 = GetOrbital(cfac, afb[m].kb)->kappa;
	jb2 = GetJFromKappa(jb2);
	for (ik = kmin; ik <= kmax; ik += 2) {
	  r0 = W6j(jlow, jb1, j1, ik, j2, jb2);
	  if (fabs(r0) < EPS30) continue;
	  r0 *= mix1*afb[m].coeff*sqrt(j1+1.0);
	  if (fabs(r0) < cfac->angz_cut) continue;
	  if (IsOdd((j1+j2-ik)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, kb1, afb[m].kb, r0);
	}
      }
      if (nfb > 0) free(afb);
    }
    PackAngularZMix(&n, ang, nz);
  } else if (kg2 < 0 && !ignore_ryd) {
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < cfac->angz_cut) continue;
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
      kb2 = sup->kcfg;
      jb2 = GetOrbital(cfac, kb2)->kappa;
      jb2 = GetJFromKappa(jb2);
      jup = GetBaseJ(cfac, sup);
      kg2 = sup->kgroup;
      kg2 = -kg2-1;
      nfb = AngularZFreeBound(cfac, &afb, kg2, lower);
      for (m = 0; m < nfb; m++) {
	jb1 = GetOrbital(cfac, afb[m].kb)->kappa;
	jb1 = GetJFromKappa(jb1);
	for (ik = kmin; ik <= kmax; ik += 2) {
	  r0 = W6j(jup, jb2, j2, ik, j1, jb1);
	  if (fabs(r0) < EPS30) continue;
	  r0 *= mix2*afb[m].coeff*sqrt(j2+1.0);
	  if (fabs(r0) < cfac->angz_cut) continue;
	  if (IsOdd((2*j1-ik+jb1-jb2)/2)) r0 = -r0;
	  im = AddToAngularZMix(&n, &nz, ang, ik, afb[m].kb, kb2, r0);
	}
      }
      if (nfb > 0) free(afb);
    }
    PackAngularZMix(&n, ang, nz);
  } else {
    if (lev1->iham > lev2->iham) {
      lev = lev1;
      lev1 = lev2;
      lev2 = lev;
    } else {
      lev = NULL;
    }
    ns = AngularZMixStates(cfac, &ad, lev1->iham, lev2->iham);
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < cfac->angz_cut) continue;
      ih1 = lev1->ibasis[i];
      isz0 = ih1 * cfac->hams[lev2->iham].nbasis;
      for (j = 0; j < lev2->n_basis; j++) {
	mix2 = lev2->mixing[j];
	if (fabs(mix2) < cfac->angz_cut) continue;
	a = mix1*mix2;
	if (fabs(a) < cfac->angz_cut) continue;
	ih2 = lev2->ibasis[j];
	isz = isz0 + ih2;
	nz_sub = (ad->nz)[isz];
	if (nz_sub > 0) {
	  ang_sub = (ad->angz)[isz];
	  for (m = 0; m < nz_sub; m++) {
	    if (ang_sub[m].k > kmax || ang_sub[m].k < kmin) continue;
	    r0 = ang_sub[m].coeff*a;
	    if (fabs(r0) < cfac->angz_cut) continue;
	    im = AddToAngularZMix(&n, &nz, ang, ang_sub[m].k, 
				  ang_sub[m].k0, ang_sub[m].k1, r0);
	  }
	}
      }
    }
    PackAngularZMix(&n, ang, nz);

    if (lev) {
      AngZSwapBraKet(cfac, n, *ang, j1-j2);
    }

    /*
    for (i = 0; i < n; i++) {
      printf("%2d %3d %2d %2d %2d %10.3E\n", lower, upper, 
	     (*ang)[i].k, (*ang)[i].k0, (*ang)[i].k1, (*ang)[i].coeff);
    }
    */
  }

  return n;
}

int AngularZxZFreeBound(cfac_t *cfac,
  ANGULAR_ZxZMIX **ang, int lower, int upper) {
  int i, j, j1, j2;
  int nz, n;
  int kg, jf, ih1, ih2, isz0;
  int ns, isz, jmin, jmax;
  double mix1, mix2;
  STATE *sup;
  SYMMETRY *sym1, *sym2;
  LEVEL *lev1, *lev2;
  ANGZ_DATUM *ad;
  ANGULAR_ZxZMIX *ang_sub;
  ANGULAR_ZMIX *ang_z;
  int kb, jb, jup, orb0, orb1;
  double r, r0, sqrt_j2;
  int nz_sub;
  int im, m;

  lev1 = GetLevel(cfac, lower);
  lev2 = GetLevel(cfac, upper);
  sym1 = GetSymmetry(cfac, lev1->pj);
  sym2 = GetSymmetry(cfac, lev2->pj);
  j1 = lev1->pj;
  j2 = lev2->pj;
  DecodePJ(j1, NULL, &j1);
  DecodePJ(j2, NULL, &j2);
  sqrt_j2 = sqrt(j2+1.0);

  nz = ANGZxZ_BLOCK;
  n = 0;
  (*ang) = malloc(sizeof(ANGULAR_ZxZMIX)*nz);
  if (!(*ang)) return -1;
    
  sup = (STATE *) ArrayGet(&(sym2->states), lev2->pb);
  if (sup->kgroup < 0) {  
    for (j = 0; j < lev2->n_basis; j++) {
      mix2 = lev2->mixing[j];
      if (fabs(mix2) < cfac->angz_cut) continue;
      sup = (STATE *) ArrayGet(&(sym2->states), lev2->basis[j]);
      kb = sup->kcfg;
      jb = GetOrbital(cfac, kb)->kappa;
      jb = GetJFromKappa(jb);
      jup = GetBaseJ(cfac, sup);
      kg = sup->kgroup;
      kg = -kg-1;
      nz_sub = AngularZMix(cfac, &ang_z, lower, kg, -1, -1);
      if (nz_sub <= 0) {
	continue;
      }
      for (i = 0; i < nz_sub; i++) {
	r = ang_z[i].coeff*mix2*sqrt_j2;
	orb0 = ang_z[i].k0;
	orb1 = ang_z[i].k1;
	jmin = abs(j1-j2);
	jmin = Max(jmin, abs(jb-ang_z[i].k));
	jmax = j1+j2;
	jmax = Min(jmax, jb+ang_z[i].k);
	for (jf = jmin; jf <= jmax; jf += 2) {
	  r0 = W6j(j1, jf, j2, jb, jup, ang_z[i].k);
	  if (fabs(r0) < EPS30) continue;
	  if (IsOdd((jup+jf+j2)/2)) r0 = -r0;
	  r0 *= r;
	  if (fabs(r0) < cfac->angz_cut) continue;
	  im = AddToAngularZxZMix(&n, &nz, ang, ang_z[i].k, 
				  jf, kb, orb0, orb1, r0);
	}
      }
      free(ang_z);
    }
  } else {
    ns = AngularZxZFreeBoundStates(cfac, &ad, lev1->iham, lev2->iham);
    for (i = 0; i < lev1->n_basis; i++) {
      mix1 = lev1->mixing[i];
      if (fabs(mix1) < cfac->angz_cut) continue;      
      ih1 = lev1->ibasis[i];
      isz0 = ih1 * cfac->hams[lev2->iham].nbasis;
      for (j = 0; j < lev2->n_basis; j++) {
	mix2 = lev2->mixing[j];
	if (fabs(mix2) < cfac->angz_cut) continue;
	r = mix1*mix2;
	if (fabs(r) < cfac->angz_cut) continue;	
	ih2 = lev2->ibasis[j];
	isz = isz0 + ih2;
	nz_sub = (ad->nz)[isz];
	if (nz_sub > 0) {
	  ang_sub = (ad->angz)[isz];
	  for (m = 0; m < nz_sub; m++) {
	    r0 = ang_sub[m].coeff*r;
	    if (fabs(r0) < cfac->angz_cut) continue;
	    im = AddToAngularZxZMix(&n, &nz, ang, 
				    ang_sub[m].k, ang_sub[m].k0,
				    ang_sub[m].k1, ang_sub[m].k2,
				    ang_sub[m].k3, r0);
	  }
	}
      }
    }
  }

  PackAngularZxZMix(&n, ang, nz);

  return n;
}
  
int CompareAngularZxZMix(const void *c1, const void *c2) {
  ANGULAR_ZxZMIX *a1, *a2;

  a1 = (ANGULAR_ZxZMIX *) c1;
  a2 = (ANGULAR_ZxZMIX *) c2;
  
  if (a1->k > a2->k) return 1;
  else if (a1->k < a2->k) return -1;
  else {
    if (a1->k0 > a2->k0) return 1;
    else if (a1->k0 < a2->k0) return -1;
    else {
      if (a1->k1 > a2->k1) return 1;
      else if (a1->k1 < a2->k1) return -1;
      else {
	if (a1->k2 > a2->k2) return 1;
	else if (a1->k2 < a2->k2) return -1;
	else {
	  if (a1->k3 > a2->k3) return 1;
	  else if (a1->k3 < a2->k3) return -1;
	  else return 0;
	}
      }
    }
  }
}
  
int CompareAngularZMix(const void *c1, const void *c2) {
  ANGULAR_ZMIX *a1, *a2;

  a1 = (ANGULAR_ZMIX *) c1;
  a2 = (ANGULAR_ZMIX *) c2;
  
  if (a1->k > a2->k) return 1;
  else if (a1->k < a2->k) return -1;
  else {
    if (a1->k0 > a2->k0) return 1;
    else if (a1->k0 < a2->k0) return -1;
    else {
      if (a1->k1 > a2->k1) return 1;
      else if (a1->k1 < a2->k1) return -1;
      else return 0;
    }
  }
}
  
int CompareAngularZFB(const void *c1, const void *c2) {
  ANGULAR_ZFB *a1, *a2;

  a1 = (ANGULAR_ZFB *) c1;
  a2 = (ANGULAR_ZFB *) c2;
  
  if (a1->kb > a2->kb) return 1;
  else if (a1->kb < a2->kb) return -1;
  else return 0;
}

int PackAngularZxZMix(int *n, ANGULAR_ZxZMIX **ang, int nz) {
  int j, m;
  ANGULAR_ZxZMIX *p1, *p2;
  
  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZxZMIX), CompareAngularZxZMix);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZxZMix(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      memcpy(p1, p2, sizeof(ANGULAR_ZxZMIX));
    }
    j++;
    p2++;
  }
  
 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZxZMIX));
      *n = m;
    }
  }

  return 0;
}

int PackAngularZMix(int *n, ANGULAR_ZMIX **ang, int nz) {
  int j, m;
  ANGULAR_ZMIX *p1, *p2;

  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZMIX), CompareAngularZMix);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZMix(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      if (p1 != p2) {
        memcpy(p1, p2, sizeof(ANGULAR_ZMIX));
      }
    }
    j++;
    p2++;
  }
  
 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZMIX));
      *n = m;
    }
  }

  return 0;
}

int PackAngularZFB(int *n, ANGULAR_ZFB **ang, int nz) {
  int j, m;
  ANGULAR_ZFB *p1, *p2;
  
  m = *n;
  if (*n <= 1) goto OUT;
  if (*n > 2) {
    qsort((void *)(*ang), *n, sizeof(ANGULAR_ZFB), CompareAngularZFB);
  }
  m = 1;
  p1 = (*ang);
  j = 1;
  p2 = p1 + 1;
  while (j < *n) {
    if (CompareAngularZFB(p1, p2) == 0) {
      p1->coeff += p2->coeff;
    } else {
      p1++;
      m++;
      memcpy(p1, p2, sizeof(ANGULAR_ZFB));
    }
    j++;
    p2++;
  }

 OUT:
  if (*n <= 0) {
    if (nz > 0) free(*ang);
  } else {
    if (m < nz) {
      (*ang) = realloc((*ang), m*sizeof(ANGULAR_ZFB));
      *n = m;
    }
  }
  return 0;
}			    

int AddToAngularZxZ(cfac_t *cfac, int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		    int n_shells, int phase, SHELL_STATE *sbra, 
		    SHELL_STATE *sket, INTERACT_SHELL *s, int m) {
  int nkk, *k, p, im;
  double *r;
  int orb0, orb1;
  int kk0, kk1;
  int kmax = GetMaxRank(cfac);
  
  nkk = AngularZxZ0(&r, &k, 0, n_shells, sbra, sket, s, kmax);
  if (nkk > 0) {    
    if (m == 0) {
      orb0 = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
      orb1 = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
      kk0 = OrbitalIndex(cfac, s[2].n, s[2].kappa, 0.0);
      kk1 = OrbitalIndex(cfac, s[3].n, s[3].kappa, 0.0);
    } else {
      orb0 = s[2].j;      
      orb1 = OrbitalIndex(cfac, s[3].n, s[3].kappa, 0.0);
      kk0 = OrbitalIndex(cfac, s[0].n, s[0].kappa, 0.0);
      kk1 = OrbitalIndex(cfac, s[1].n, s[1].kappa, 0.0);
    }
    for (p = 0; p < nkk; p++) {
      if (fabs(r[p]) < EPS30) continue;
      if (IsOdd(phase)) r[p] = -r[p];
      im = AddToAngularZxZMix(n, nz, ang, k[p], 
			      orb0, orb1, kk0, kk1, r[p]);
    }
    free(r);
    free(k);
  }
  
  return 0;
}

int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double r) {
  int im;

  im = *n;
  (*n)++;
  if (*n > *nz) {
    *nz += ANGZxZ_BLOCK;
    *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZxZMIX));
    if (!(*ang)) {
      printf("Can't enlarge AngularZxZMix array\n");
      return -1;
    }
  }
  (*ang)[im].k = k;
  (*ang)[im].k0 = k0;
  (*ang)[im].k1 = k1;
  (*ang)[im].k2 = k2;
  (*ang)[im].k3 = k3;
  (*ang)[im].coeff = r;
  
  return 0;
}

int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff) {
  int im;
  
  im = *n;
  (*n)++;
  if (*n > *nz) {
    *nz += ANGZ_BLOCK;
    *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZMIX));
    if (!(*ang)) {
      printf("Can't enlarge AngularZMix array\n");
      return -1;
    }
  }
  (*ang)[im].k = k;
  (*ang)[im].k0 = k0;
  (*ang)[im].k1 = k1;
  (*ang)[im].coeff = coeff;  

  return 0;
}

int AddToAngularZFB(int *n, int *nz, ANGULAR_ZFB **ang,
		    int kb, double coeff) {
  int im;

  im = *n;
  (*n)++;
  if (*n > *nz) {
    *nz += ANGZ_BLOCK;
    *ang = realloc((*ang), (*nz)*sizeof(ANGULAR_ZFB));
    if (!(*ang)) {
      printf("Cannot enlarge AngularZFB array\n");
      return -1;
    }
  }
  (*ang)[im].kb = kb;
  (*ang)[im].coeff = coeff;

  return 0;
}

void FreeAngZDatum(ANGZ_DATUM *ap) {
  int i;

  for (i = 0; i < ap->ns; i++) {
    if (ap->nz[i] > 0) free(ap->angz[i]);
  }
  if (ap->ns > 0) {
    free(ap->angz);
    free(ap->nz);
    if (ap->mk) {
      for (i = 0; i < ap->ns*2; i++) {
	free(ap->mk[i]);
      }
      free(ap->mk);
    }
  }
  ap->ns = 0;
}

void FreeLevelData(void *p) {
  LEVEL *lev;
  lev = (LEVEL *) p;
  if (lev->n_basis > 0) {
    free(lev->basis);
    free(lev->ibasis);
    free(lev->mixing);
    lev->n_basis = 0;
  }
}

void FreeHamsArray(cfac_t *cfac) {
  int i;

  for (i = 0; i < cfac->nhams; i++) {
    if (cfac->hams[i].nbasis > 0) {
      free(cfac->hams[i].basis);
    }
    cfac->hams[i].nbasis = 0;
  }
  cfac->nhams = 0;
}

void ClearAngularFrozen(cfac_t *cfac) {
  int i, n;

  if (cfac->ang_frozen.nts > 0) {
    free(cfac->ang_frozen.ts);
    n = cfac->ang_frozen.nts*cfac->ang_frozen.nts;
    for (i = 0; i < n; i++) {
      if (cfac->ang_frozen.nz[i] > 0) free(cfac->ang_frozen.z[i]);
    }
    free(cfac->ang_frozen.nz);
    free(cfac->ang_frozen.z);
    cfac->ang_frozen.ts = NULL;
    cfac->ang_frozen.nz = NULL;
    cfac->ang_frozen.z = NULL;
  }
  if (cfac->ang_frozen.ncs > 0) {
    free(cfac->ang_frozen.cs);
    n = cfac->ang_frozen.nts*cfac->ang_frozen.ncs;
    for (i = 0; i < n; i++) {
      if (cfac->ang_frozen.nzfb[i] > 0) free(cfac->ang_frozen.zfb[i]);
      if (cfac->ang_frozen.nzxzfb[i] > 0) free(cfac->ang_frozen.zxzfb[i]);
    }
    free(cfac->ang_frozen.nzfb);
    free(cfac->ang_frozen.nzxzfb);
    free(cfac->ang_frozen.zfb);
    free(cfac->ang_frozen.zxzfb);
    cfac->ang_frozen.cs = NULL;
    cfac->ang_frozen.nzfb = NULL;
    cfac->ang_frozen.nzxzfb = NULL;
    cfac->ang_frozen.zfb = NULL;
    cfac->ang_frozen.zxzfb = NULL;
  }
  cfac->ang_frozen.nts = 0;
  cfac->ang_frozen.ncs = 0;
}

/* (Re)allocate Hamiltonian */
int AllocHamMem(HAMILTON *h, int hdim, int nbasis) {
  int jp, t;

  jp = nbasis - hdim;
  if (jp < 0) return -1;

  h->dim     = hdim;
  h->n_basis = nbasis;
  
  /* size of the triangular matrix, including diagonal elements */
  t = hdim*(hdim+1)/2;
  
  /* matrix partitioned as: H1[t] + B[hdim*jp] + H2[jp] */
  h->hsize = t + hdim*jp + jp;
  
  /* start reallocations as needed */
  if (h->n_basis > h->n_basis0) {
    h->basis = realloc(h->basis, sizeof(int)*(h->n_basis));
    h->n_basis0 = h->n_basis;
  }
  if (!h->basis) return -1;
  
  if (h->hsize > h->hsize0) {
    h->hsize0 = h->hsize;
    h->hamilton = realloc(h->hamilton, sizeof(double)*h->hsize);
  }
  if (!h->hamilton) return -1;

  h->msize = h->dim * h->n_basis + h->dim;  
  if (h->msize > h->msize0) {
    h->msize0 = h->msize;
    h->mixing = realloc(h->mixing, sizeof(double)*h->msize);
  }
  if (!h->mixing) return -1;

  return 0;
}
