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
#include <stdio.h>
#include <math.h>

#include "cfacP.h"
#include "grid.h"

static int AddPW(const cfac_t *cfac, int *nkl0, double *kl, double *logkl,
          int maxkl, int n, int step) {
  int i;
  for (i = *nkl0; i < n+(*nkl0); i++) {
    if (i >= MAXNKL) {
      cfac_errmsg(cfac, "Maximum partial wave grid points reached: ");
      cfac_errmsg(cfac, "%d > %d in constructing grid\n",  i, MAXNKL);
      return CFAC_FAILURE;
    }
    kl[i] = kl[i-1] + step;
    logkl[i] = log(kl[i]);
    if ((int)(kl[i]) > maxkl) break;
  }
  (*nkl0) = i;
  return CFAC_SUCCESS;
}

int SetPWGrid(const cfac_t *cfac, int *nkl0, double *kl, double *logkl,
              int maxkl, int *ns, int *n, int *step) {
  int i, m, k, j;

  if ((*ns) > 0) {
    for (i = 0; i < (*ns); i++) {
      if (AddPW(cfac, nkl0, kl, logkl, maxkl, n[i], step[i]) != CFAC_SUCCESS) {
        return CFAC_FAILURE;
      }
    }
    k = step[(*ns)-1]*2;
    j = 2;
  } else {
    (*ns) = -(*ns);
    if ((*ns) == 0) (*ns) = 8;
    if (AddPW(cfac, nkl0, kl, logkl, maxkl, (*ns), 1) != CFAC_SUCCESS) {
      return CFAC_FAILURE;
    }
    k = 2;
    j = 2;
  }

  m = kl[(*nkl0)-1];
  while (m+k <= maxkl) {
    if (AddPW(cfac, nkl0, kl, logkl, maxkl, j, k) != CFAC_SUCCESS) {
      return CFAC_FAILURE;
    }
    m = kl[(*nkl0)-1];
    if (k < 8) k *= 2;
    else k += 4;
  }
  kl[(*nkl0)] = maxkl+1;
  return CFAC_SUCCESS;
}

int SetTEGridDetail(const cfac_t *cfac,
  double *te, double *logte, int n, double *x) {
  int i;

  if (n > MAXNTE) {
    cfac_errmsg(cfac, "Max # of grid points reached \n");
    return -1;
  }

  for (i = 0; i < n; i++) {
    te[i] = x[i];
    logte[i] = log(te[i]);
  }
  return n;
}

int SetTEGrid(const cfac_t *cfac,
  double *te, double *logte, int n, double emin, double emax) {
  int i;
  double del;

  if (n < 1) {
    te[0] = -1.0;
    return 0;
  }

  if (emin < 0.0) {
    te[0] = emin;
    return n;
  }

  if (n > MAXNTE) {
    cfac_errmsg(cfac, "Max # of grid points reached \n");
    return -1;
  }

  if (n == 1) {
    te[0] = (emin + emax)/2;
    if (logte) logte[0] = log(te[0]);
    return n;
  }

  if (n == 2) {
    te[0] = emin;
    te[1] = emax;
    if (logte) {
      logte[0] = log(emin);
      logte[1] = log(emax);
    }
    return n;
  }

  if (emax < emin) {
    cfac_errmsg(cfac, "emin must > 0 and emax < emin in SetTEGrid\n");
    return -1;
  }

  del = emax - emin;
  del /= n-1.0;
  te[0] = emin;
  if (logte) logte[0] = log(emin);
  for (i = 1; i < n; i++) {
    te[i] = te[i-1] + del;
    if (logte) logte[i] = log(te[i]);
  }

  return n;
}

int SetEGridDetail(double *e, double *log_e, int n, double *xg) {
  int i;

  for (i = 0; i < n; i++) {
    e[i] = xg[i];
    log_e[i] = log(e[i]);
  }

  return n;
}

static double EFromX(const cfac_t *cfac, double x, double b) {
  double x0, e, a, d;
  int i;

  x0 = 2.0*log(1.0/FINE_STRUCTURE_CONST2);
  if (x > x0) e = x/b;
  else e = exp(x);

  for (i = 0; i < 100; i++) {
    d = log(e) + b*e - x;
    if (fabs(d/x) < EPS5) return e;
    a = 1.0/e + b;
    e = e - d/a;
  }
  cfac_errmsg(cfac, "Newton iteration failed to converge in EFromX, %10.3E %10.3E %10.3E\n",
         x, e, d);
  return e;
}

int SetEGrid(const cfac_t *cfac, double *e, double *log_e,
             int n, double emin, double emax, double eth) {
  double del, et, b;
  int i;

  if (n < 1) {
    e[0] = -1.0;
    return 0;
  }
  if (emin < 0.0) {
    e[0] = emin;
    return n;
  }

  if (emax < emin) {
    cfac_errmsg(cfac, "emin must > 0 and emax < emin in SetEGrid\n");
    return -1;
  }

  if (eth >= 0.0) {
    et = eth;
    emin += et;
    emax += et;

    b = 0.2/FINE_STRUCTURE_CONST2;
    b = log(b)/b;
    e[0] = emin;
    log_e[0] = log(emin) + b*emin;
    if (n == 1) goto DONE1;
    e[n-1] = emax;
    log_e[n-1] = log(emax) + b*emax;
    if (n == 2) goto DONE1;
    del = (log_e[n-1] - log_e[0])/(n-1.0);
    for (i = 1; i < n-1; i++) {
      log_e[i] = log_e[i-1]+del;
      e[i] = EFromX(cfac, log_e[i], b);
    }
  DONE1:
    for (i = 0; i < n; i++) {
      e[i] -= eth;
      log_e[i] = log(e[i]);
    }
  } else if (eth < 0.0) {
    et = -eth;
    e[0] = emin/(emin+et);
    if (n == 1) goto DONE2;
    e[n-1] = emax/(emax+et);
    if (n == 2) goto DONE2;
    del = (e[n-1]-e[0])/(n-1.0);
    for (i = 1; i < n-1; i++) {
      e[i] = e[i-1] + del;
    }

  DONE2:
    for (i = 0; i < n; i++) {
      e[i] = e[i]*et/(1.0-e[i]);
      log_e[i] = log(e[i]);
    }
  }

  return n;
}

int SetLinearGrid(double *x, int n, double xmin, double xmax) {
  int i;
  double d;

  if (n == 0) return 0;
  if (n == 1) {
    x[0] = 0.5*(xmin + xmax);
    return 1;
  }
  d = (xmax - xmin)/(n-1);
  x[0] = xmin;
  for (i = 1; i < n; i++) {
    x[i] = x[i-1] + d;
  }

  return n;
}
