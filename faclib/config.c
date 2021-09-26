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

/*************************************************************
  Implementation of module "config".

  This module generates electron configuations and
  carries out the angular momentum coupling.

  Author: M. F. Gu, mfgu@space.mit.edu
**************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "cfacP.h"
#include "parser.h"

/*
** VARIABLE:    spec_symbols
** TYPE:        string.
** PURPOSE:     a string orbital angular momentum symbols.
** NOTE:        the last char "*" is not part of the symbol,
**              rather, it represents any of the previous symbols.
*/
static char spec_symbols[MAX_SPEC_SYMBOLS+2] = "spdfghiklmnoqrtuvwxyz*";

static int AddConfigToSymmetry(cfac_t *cfac, int kg, int kc, CONFIG *cfg);

/*
** FUNCTION:    DistributeElectronsShell
** PURPOSE:     distribute nq electrons among the specified shells
**              to construct all possible configurations.
** INPUT:       {CONFIG **cfg},
**              pointer to a pointer of CONFIG, which holds the
**              resulting configrations on output.
**              {int ns},
**              number of shells.
**              {SHELL *shells}
**              an array of shells.
**              {int nq},
**              number of electrons to be distributed.
**              {int *maxq},
**              maxq[i] is the maximum number of electrons allowed
**              in all the shells except the the i-th shell.
**              this is to determine the minimum number of electrons
**              allowed in the i-th shell, nq-maxq[i].
** RETURN:      {int},
**              number of configurations.
** SIDE EFFECT:
** NOTE:        This is a static function only used in module "config".
*/
static int DistributeElectronsShell(CONFIG **cfg, int ns, SHELL *shell,
                                int nq, int *maxq) {
  CONFIG **cfg1, **cfg2;
  int *ncfg2;
  int qmin, qmax, j, q, t, k, ncfg;

  if (nq < 0) return 0;

  if (nq == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    (*cfg)->n_shells = 0;
    return 1;
  }

  if (nq == 1){
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    for (t = 0; t < ns; t++) {
      (*cfg)[t].n_shells = 1;
      (*cfg)[t].shells = (SHELL *) malloc(sizeof(SHELL));
      (*cfg)[t].shells[0].n = shell[t].n;
      (*cfg)[t].shells[0].kappa = shell[t].kappa;
      (*cfg)[t].shells[0].nq = 1;
    }
    return ns;
  }

  if (ns == 1) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    (*cfg)->shells[0].n = shell[0].n;
    (*cfg)->shells[0].kappa = shell[0].kappa;
    (*cfg)->shells[0].nq = nq;
    return 1;
  }

  j = GetJFromKappa(shell[0].kappa);
  qmax = j+1;
  qmax = Min(qmax, nq);
  qmin = nq - maxq[0];
  qmin = Max(qmin, 0);
  ncfg = 0;
  t = qmax-qmin+1;
  cfg1 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  cfg2 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  ncfg2 = (int *) malloc(sizeof(int)*t);
  t = 0;
  for (q = qmin; q <= qmax; q++) {
    DistributeElectronsShell(cfg1+t, 1, shell, q, NULL);
    ncfg2[t] = DistributeElectronsShell(cfg2+t, ns-1, shell+1, nq-q, maxq+1);
    ncfg += ncfg2[t];
    t++;
  }

  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  t = 0;
  k = 0;
  for (q = qmin; q <= qmax; q++) {
    for (j = 0; j < ncfg2[t]; j++) {
      (*cfg)[k].n_shells = cfg1[t]->n_shells + cfg2[t][j].n_shells;
      (*cfg)[k].shells = (SHELL *) malloc(sizeof(SHELL)*(*cfg)[k].n_shells);
      if (cfg1[t]->n_shells > 0) {
        memcpy((*cfg)[k].shells, cfg1[t]->shells, sizeof(SHELL));
      }
      if (cfg2[t][j].n_shells > 0) {
        memcpy((*cfg)[k].shells+cfg1[t]->n_shells, cfg2[t][j].shells,
               sizeof(SHELL)*(cfg2[t][j].n_shells));
      }
      if (cfg2[t][j].n_shells > 0) free(cfg2[t][j].shells);
      k++;
    }
    if (cfg1[t]->n_shells > 0) free(cfg1[t]->shells);
    free(cfg1[t]);
    free(cfg2[t]);
    t++;
  }
  free(cfg1);
  free(cfg2);
  free(ncfg2);

  return ncfg;
}

static int DistributeElectronsShellNR(CONFIG **cfg, int ns, SHELL *shell,
                                      int nq, int *maxq) {
  CONFIG **cfg1, **cfg2;
  int *ncfg2;
  int qmin, qmax, j, q, t, k, ncfg;

  if (nq < 0) return 0;

  if (nq == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    (*cfg)->n_shells = 0;
    return 1;
  }

  if (nq == 1){
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ns);
    for (t = 0; t < ns; t++) {
      (*cfg)[t].n_shells = 1;
      (*cfg)[t].shells = (SHELL *) malloc(sizeof(SHELL));
      (*cfg)[t].shells[0].n = shell[t].n;
      (*cfg)[t].shells[0].kappa = shell[t].kappa;
      (*cfg)[t].shells[0].nq = 1;
    }
    return ns;
  }

  if (ns == 1) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    (*cfg)->shells[0].n = shell[0].n;
    (*cfg)->shells[0].kappa = shell[0].kappa;
    (*cfg)->shells[0].nq = nq;
    return 1;
  }

  j = shell[0].kappa;
  qmax = 2.0*(j+1);
  qmax = Min(qmax, nq);
  qmin = nq - maxq[0];
  qmin = Max(qmin, 0);
  ncfg = 0;
  t = qmax-qmin+1;
  cfg1 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  cfg2 = (CONFIG **) malloc(sizeof(CONFIG *)*t);
  ncfg2 = (int *) malloc(sizeof(int)*t);
  t = 0;
  for (q = qmin; q <= qmax; q++) {
    DistributeElectronsShellNR(cfg1+t, 1, shell, q, NULL);
    ncfg2[t] = DistributeElectronsShellNR(cfg2+t, ns-1, shell+1, nq-q, maxq+1);
    ncfg += ncfg2[t];
    t++;
  }

  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  t = 0;
  k = 0;
  for (q = qmin; q <= qmax; q++) {
    for (j = 0; j < ncfg2[t]; j++) {
      (*cfg)[k].n_shells = cfg1[t]->n_shells + cfg2[t][j].n_shells;
      (*cfg)[k].shells = (SHELL *) malloc(sizeof(SHELL)*(*cfg)[k].n_shells);
      if (cfg1[t]->n_shells > 0) {
        memcpy((*cfg)[k].shells, cfg1[t]->shells, sizeof(SHELL));
      }
      if (cfg2[t][j].n_shells > 0) {
        memcpy((*cfg)[k].shells+cfg1[t]->n_shells, cfg2[t][j].shells,
               sizeof(SHELL)*(cfg2[t][j].n_shells));
      }
      if (cfg2[t][j].n_shells > 0) free(cfg2[t][j].shells);
      k++;
    }
    if (cfg1[t]->n_shells > 0) free(cfg1[t]->shells);
    free(cfg1[t]);
    free(cfg2[t]);
    t++;
  }
  free(cfg1);
  free(cfg2);
  free(ncfg2);

  return ncfg;
}

int ShellsFromString(const char *scfg, double *dnq, SHELL **shell) {
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int i, j, t, k, kl2, ns;
  char *s;

  SetParserQuote("[", "]");
  SetParserBreak(spec_symbols);
  SetParserWhite("");
  SetParserEscape('\0');

  next = 0;
  r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  if (r) {
    return -1;
  }
  if (quotepos == 0) {
    nn = StrSplit(token, ',');
    if (nn > 16) {
      printf("number of n's in a single shell must be <= 16\n");
      return -1;
    }
    s = token;
    for (i = 0; i < nn; i++) {
      while (*s == ' ' || *s == '\t') s++;
      n[i] = atoi(s);
      while (*s) s++;
      s++;
    }
  } else {
    nn = 1;
    n[0] = atoi(token);
  }
  if (brkpos < 0) {
    r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
    if (r) {
      return -1;
    }
  }
  if (brkpos >= 0) {
    kl[0] = brkpos;
    if (brkpos == MAX_SPEC_SYMBOLS) {
      if (n[nn-1] >= 512) {
        printf("not all L-terms are allowed for n >= %d\n", 512);
        return -1;
      }
      nkl = n[nn-1];
      for (i = 0; i < nkl; i++) {
        kl[i] = i;
      }
    } else {
      nkl = 1;
      kl[0] = brkpos;
    }
    nkappa = 0;
    for (i = 0; i < nkl; i++) {
      kl2 = 2*kl[i];
      if (kl2 > 0) {
        kappa[nkappa++] = GetKappaFromJL(kl2-1, kl2);
      }
      kappa[nkappa++] = GetKappaFromJL(kl2+1, kl2);
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq < 0) {
      printf("negative occupation number, use brackets to indicate j-1/2 shell\n");
      return -1;
    }
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else if (quotepos == 0) {
    nkl = StrSplit(token, ',');
    if (nkl > 512) {
      printf("number of L-terms must < 512\n");
      return -1;
    }
    s = token;
    nkappa = 0;
    for (k = 0; k < nkl; k++) {
      while (*s == ' ' || *s == '\t') s++;
      GetJLFromSymbol(s, &j, &kl[k]);
      kl2 = 2*kl[k];
      if (j != 1 && kl2 > 0) {
        kappa[nkappa++] = GetKappaFromJL(kl2-1, kl2);
      }
      if (j != -1) {
        kappa[nkappa++] = GetKappaFromJL(kl2+1, kl2);
      }
      while (*s) s++;
      s++;
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else {
    return -1;
  }

  *shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = GetLFromKappa(kappa[k]);
      if (kl2/2 >= n[i]) continue;
      (*shell)[t].n = n[i];
      (*shell)[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;

  if (ns == 0) {
    free(*shell);
    return -1;
  }

  return ns;
}

int ShellsFromStringNR(const char *scfg, double *dnq, SHELL **shell) {
  char token[128];
  int r, brkpos, quotepos, next;
  int nn, nkl, nkappa;
  int n[16];
  int kl[512];
  int kappa[1024];
  int i, j, t, k, kl2, ns;
  char *s;

  SetParserQuote("[", "]");
  SetParserBreak(spec_symbols);
  SetParserWhite("");
  SetParserEscape('\0');

  next = 0;
  r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
  if (r) {
    return -1;
  }
  if (quotepos == 0) {
    nn = StrSplit(token, ',');
    if (nn > 16) {
      printf("number of n's in a single shell must be <= 16\n");
      return -1;
    }
    s = token;
    for (i = 0; i < nn; i++) {
      while (*s == ' ' || *s == '\t') s++;
      n[i] = atoi(s);
      while (*s) s++;
      s++;
    }
  } else {
    nn = 1;
    n[0] = atoi(token);
  }
  if (brkpos < 0) {
    r = Parse(token, 512, scfg, &next, &brkpos, &quotepos);
    if (r) {
      return -1;
    }
  }
  if (brkpos >= 0) {
    kl[0] = brkpos;
    if (brkpos == MAX_SPEC_SYMBOLS) {
      if (n[nn-1] >= 512) {
        printf("not all L-terms are allowed for n >= %d\n", 512);
        return -1;
      }
      nkl = n[nn-1];
      for (i = 0; i < nkl; i++) {
        kl[i] = i;
      }
    } else {
      nkl = 1;
      kl[0] = brkpos;
    }
    nkappa = 0;
    for (i = 0; i < nkl; i++) {
      kl2 = 2*kl[i];
      kappa[nkappa++] = kl2;
    }
    if (scfg[next] == '+' || scfg[next] == '-') next++;
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else if (quotepos == 0) {
    nkl = StrSplit(token, ',');
    if (nkl > 512) {
      printf("number of L-terms must < 512\n");
      return -1;
    }
    s = token;
    nkappa = 0;
    for (k = 0; k < nkl; k++) {
      while (*s == ' ' || *s == '\t') s++;
      GetJLFromSymbol(s, &j, &kl[k]);
      kl2 = 2*kl[k];
      kappa[nkappa++] = kl2;
      while (*s) s++;
      s++;
    }
    *dnq = atof(&(scfg[next]));
    if (*dnq == 0 && !(isdigit(scfg[next]))) *dnq = 1;
  } else {
    return -1;
  }

  *shell = (SHELL *) malloc(sizeof(SHELL)*nn*nkappa);
  t = 0;
  for (i = nn-1; i >= 0; i--) {
    for (k = nkappa-1; k >= 0; k--) {
      kl2 = kappa[k];
      if (kl2/2 >= n[i]) continue;
      (*shell)[t].n = n[i];
      (*shell)[t].kappa = kappa[k];
      t++;
    }
  }
  ns = t;

  if (ns == 0) {
    free(*shell);
    return -1;
  }

  return ns;
}

int GetRestriction(const char *iscfg, SHELL_RESTRICTION **sr, int m) {
  int nc, i;
  double dnq;
  char *s, *p;

  char * const scfg = strdup(iscfg);
  p = scfg;
  nc = StrSplit(p, ';');
  nc--;

  if (nc > 0) {
    *sr = malloc(sizeof(SHELL_RESTRICTION)*nc);
  } else {
    *sr = NULL;
  }

  for (i = 0; i < nc; i++) {
    while (*p) p++;
    p++;
    s = p;
    (*sr)[i].nq = -1;
    while (*s) {
      switch (*s) {
      case '=':
        (*sr)[i].op = 0;
        *s = '\0';
        (*sr)[i].nq = atoi(s+1);
        break;
      case '<':
        (*sr)[i].op = -1;
        *s = '\0';
        (*sr)[i].nq = atoi(s+1);
        break;
      case '>':
        (*sr)[i].op = 1;
        *s = '\0';
        (*sr)[i].nq = atoi(s+1);
        break;
      default:
        break;
      }
      s++;
    }
    if (m == 0) {
      (*sr)[i].ns = ShellsFromString(p, &dnq, &((*sr)[i].shells));
    } else {
      (*sr)[i].ns = ShellsFromStringNR(p, &dnq, &((*sr)[i].shells));
    }
    if ((*sr)[i].nq < 0) {
      (*sr)[i].op = 0;
      (*sr)[i].nq = dnq;
    }
    p = s;
  }

  free(scfg);
  return nc;
}

int ApplyRestriction(int ncfg, CONFIG *cfg, int nc, SHELL_RESTRICTION *sr) {
  int i, j, k, t, nq, c;

  for (k = 0; k < nc; k++) {
    for (i = 0; i < ncfg; i++) {
      nq = 0;
      for (j = 0; j < cfg[i].n_shells; j++) {
        for (t = 0; t < sr[k].ns; t++) {
          if (cfg[i].shells[j].n == sr[k].shells[t].n &&
              cfg[i].shells[j].kappa == sr[k].shells[t].kappa) {
            nq += cfg[i].shells[j].nq;
          }
        }
      }
      switch (sr[k].op) {
      case -1:
        c = (nq < sr[k].nq);
        break;
      case 0:
        c = (nq == sr[k].nq);
        break;
      case 1:
        c = (nq > sr[k].nq);
        break;
      default:
        c = 0;
        break;
      }
      if (c == 0) {
        if (cfg[i].n_shells > 0) {
          cfg[i].n_shells = 0;
          free(cfg[i].shells);
          cfg[i].shells = NULL;
        }
      }
    }
  }

  i = 0;
  j = 0;
  while (i < ncfg) {
    if (cfg[i].n_shells == 0) {
      i++;
    } else {
      cfg[j].n_shells = cfg[i].n_shells;
      cfg[j].shells = cfg[i].shells;
      j++;
      i++;
    }
  }

  return j;
}

/*
** FUNCTION:    DistributeElectrons
** PURPOSE:     Decode a single string shell, distribute electrons
**              among all physical shells if the configurations
**              are not the average configurations, otherwise, the
**              average configurations is constructed with all
**              shells present, and the number of electrons returned.
** INPUT:       {CONFIG **cfg},
**              pointer to a pointer to CONFIG, which holds the
**              resulting configurations or average configurations.
**              {double *nq},
**              pointer to double, which will hold the total number of
**              electrons for the average configuration. it should be
**              passed in as NULL if the actual configurations to be
**              constructed.
**              {char *},
**              a single string shell.
** RETURN:      {int},
**              if actual configurations to be constructed (nq == NULL),
**              return the number of configurations constructed.
**              if average configuration is to be determined,
**              return the number of shells in the average configuration.
** SIDE EFFECT:
** NOTE:
*/
int DistributeElectrons(CONFIG **cfg, double *nq, const char *scfg) {
  SHELL *shell;
  int ncfg, *maxq, ns, nc, i, j, inq;
  double dnq;
  SHELL_RESTRICTION *sr;

  nc = GetRestriction(scfg, &sr, 0);

  ns = ShellsFromString(scfg, &dnq, &shell);
  if (ns <= 0) {
    if (nc > 0) {
      for (i = 0; i < nc; i++) {
        free(sr[i].shells);
      }
      free(sr);
    }
    return ns;
  }

  if (nq) {
    *nq = dnq;
    *cfg = malloc(sizeof(CONFIG));
    (*cfg)->n_shells = ns;
    (*cfg)->shells = malloc(sizeof(SHELL)*ns);
    memcpy((*cfg)->shells, shell, sizeof(SHELL)*ns);
    free(shell);
    return ns;
  }

  maxq = malloc(sizeof(int)*ns);
  maxq[ns-1] = 0;
  for (i = ns-2; i >= 0; i--) {
    j = GetJFromKappa(shell[i+1].kappa);
    maxq[i] = maxq[i+1] + j+1;
  }

  inq = (int) dnq;
  ncfg = DistributeElectronsShell(cfg, ns, shell, inq, maxq);

  free(shell);
  free(maxq);

  if (nc > 0) {
    ncfg = ApplyRestriction(ncfg, *cfg, nc, sr);
    for (i = 0; i < nc; i++) {
      free(sr[i].shells);
    }
    free(sr);
  }

  return ncfg;
}

int DistributeElectronsNR(CONFIG **cfg, const char *scfg) {
  SHELL *shell;
  int ncfg, *maxq, ns, nc, i, inq;
  double dnq;
  SHELL_RESTRICTION *sr;

  nc = GetRestriction(scfg, &sr, 1);

  ns = ShellsFromStringNR(scfg, &dnq, &shell);

  maxq = (int *) malloc(sizeof(int)*ns);
  maxq[ns-1] = 0;
  for (i = ns-2; i >= 0; i--) {
    maxq[i] = maxq[i+1] + 2*(shell[i+1].kappa + 1);
  }

  inq = (int) dnq;
  ncfg = DistributeElectronsShellNR(cfg, ns, shell, inq, maxq);

  free(shell);
  free(maxq);

  if (nc > 0) {
    ncfg = ApplyRestriction(ncfg, *cfg, nc, sr);
    for (i = 0; i < nc; i++) {
      free(sr[i].shells);
    }
    free(sr);
  }

  return ncfg;
}

/*
** FUNCTION:    GetConfigOrAverageFromString
** PURPOSE:     decode the string representation of configurations,
**              construct all possible configurations or average
**              configuration.
** INPUT:       {CONFIG **cfg}
**              holds the resuting configurations or average configuration.
**              {double **nq},
**              return fractional occupation numbers of each shell in the
**              average configuration, if it is not NULL on input.
**              {char *scfg},
**              a string representation of configurations.
** RETURN:      {int},
**              if nq == NULL, return the number of configurations
**              constructed.
**              if (nq != NULL, return the number of shells in the
**              average configuration.
** SIDE EFFECT:
** NOTE:
*/
int GetConfigOrAverageFromString(const cfac_t *cfac,
    CONFIG **cfg, double **nq, const char *iscfg) {
  CONFIG **dcfg, **p1;
  double *dnq, *p2, a, b;
  char *s;
  int *isp, ncfg, *dnc;
  int size, size_old, tmp;
  int i, t, j, k, ns;

  char * const scfg = strdup(iscfg);
  StrTrim(scfg, '\0');
  ns = QuotedStrSplit(cfac, scfg, ' ', '[', ']');
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->uta = 0;
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    free(scfg);
    return 1;
  }

  dcfg = (CONFIG **) malloc(sizeof(CONFIG *)*ns);
  dnc = (int *) malloc(sizeof(int)*ns);
  if (nq) {
    dnq = (double *) malloc(sizeof(double)*ns);
    p2 = dnq;
  } else {
    dnq = NULL;
    p2 = NULL;
  }

  s = scfg;
  isp = (int *) malloc(sizeof(int)*ns);
  isp[0] = 0;
  for (i = 1; i < ns; i++) {
    isp[i] = isp[i-1];
    while (s[isp[i]]) isp[i]++;
    isp[i]++;
  }
  p1 = dcfg;
  t = 0;
  for (i = 0; i < ns; i++) {
    s = scfg + isp[i];
    while (*s == ' ' || *s == '\t') s++;
    dnc[t] = DistributeElectrons(p1, p2, s);
    if (dnc[t] <= 0) {
      free(dcfg);
      free(dnc);
      free(isp);
      free(scfg);
      return -1;
    }
    if (dnc[t] > 1 || (*p1)->n_shells > 0) {
      p1++;
      t++;
      if (p2) p2++;
    } else {
      free(*p1);
    }
  }
  free(isp);

  ns = t;
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    free(dcfg);
    free(dnc);
    free(scfg);
    return 1;
  }

  if (!nq) {
    ncfg = dnc[0];
    for (i = 1; i < ns; i++) ncfg *= dnc[i];
    *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
    tmp = ncfg;
    p1 = dcfg + ns - 1;
    for (i = ns-1; i >= 0; i--) {
      tmp /= dnc[i];
      t = 0;
      while (t < ncfg) {
        for (j = 0; j < dnc[i]; j++) {
          for (k = 0; k < tmp; k++) {
            if (i == ns-1) {
              (*cfg)[t].n_shells = (*p1)[j].n_shells;
              size = sizeof(SHELL)*(*p1)[j].n_shells;
              (*cfg)[t].shells = (SHELL *) malloc(size);
              memcpy((*cfg)[t].shells, (*p1)[j].shells, size);
            } else {
              size_old = sizeof(SHELL)*(*cfg)[t].n_shells;
              size = sizeof(SHELL)*(*p1)[j].n_shells;
              (*cfg)[t].shells = (SHELL *) realloc((*cfg)[t].shells,
                                                   size_old+size);
              memcpy((*cfg)[t].shells+(*cfg)[t].n_shells,
                     (*p1)[j].shells, size);
              (*cfg)[t].n_shells += (*p1)[j].n_shells;
            }
            t++;
          }
        }
      }
      p1--;
    }
  } else {
    ncfg = dnc[0];
    for (i = 1; i < ns; i++) ncfg += dnc[i];
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)[0].n_shells = ncfg;
    (*cfg)[0].shells = (SHELL *) malloc(sizeof(SHELL)*ncfg);
    *nq = (double *) malloc(sizeof(double)*ncfg);
    p1 = dcfg + ns-1;
    t = 0;
    for (i = ns-1; i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < dnc[i]; j++) {
        a += GetJFromKappa((*p1)->shells[j].kappa) + 1.0;
      }
      for (j = 0; j < dnc[i]; j++) {
        b = GetJFromKappa((*p1)->shells[j].kappa) + 1.0;
        (*nq)[t] = dnq[i]*b/a;
        memcpy((*cfg)->shells+t, (*p1)->shells+j, sizeof(SHELL));
        t++;
      }
      p1--;
    }
  }

  for (i = 0; i < ns; i++) {
    if (!nq) {
      for (j = 0; j < dnc[i]; j++) {
        if ((dcfg[i][j]).n_shells > 0) {
          free((dcfg[i][j]).shells);
        }
      }
    } else {
      if (dcfg[i][0].n_shells > 0) {
        free((dcfg[i][0]).shells);
      }
    }
    free(dcfg[i]);
  }
  free(dcfg);
  free(dnc);
  free(scfg);

  return ncfg;
}

int GetConfigFromStringNR(const cfac_t *cfac, CONFIG **cfg, const char *iscfg) {
  CONFIG **dcfg, **p1;
  char *s;
  int *isp, ncfg, *dnc;
  int size, size_old, tmp;
  int i, t, j, k, ns;

  char * const scfg = strdup(iscfg);
  StrTrim(scfg, '\0');
  ns = QuotedStrSplit(cfac, scfg, ' ', '[', ']');
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    (*cfg)->shells->kappa = 0;
    return 1;
  }

  dcfg = (CONFIG **) malloc(sizeof(CONFIG *)*ns);
  dnc = (int *) malloc(sizeof(int)*ns);

  s = scfg;
  isp = (int *) malloc(sizeof(int)*ns);
  isp[0] = 0;
  for (i = 1; i < ns; i++) {
    isp[i] = isp[i-1];
    while (s[isp[i]]) isp[i]++;
    isp[i]++;
  }
  p1 = dcfg;
  t = 0;
  for (i = 0; i < ns; i++) {
    s = scfg + isp[i];
    while (*s == ' ' || *s == '\t') s++;
    dnc[t] = DistributeElectronsNR(p1, s);
    if (dnc[t] <= 0) {
      free(dcfg);
      free(dnc);
      free(isp);
      free(scfg);
      return -1;
    }
    if (dnc[t] > 1 || (*p1)->n_shells > 0) {
      p1++;
      t++;
    } else {
      free(*p1);
    }
  }
  free(isp);

  ns = t;
  if (ns == 0) {
    *cfg = (CONFIG *) malloc(sizeof(CONFIG));
    (*cfg)->n_shells = 1;
    (*cfg)->shells = (SHELL *) malloc(sizeof(SHELL));
    PackShell((*cfg)->shells, 1, 0, 1, 0);
    (*cfg)->shells->kappa = 0;
    free(dcfg);
    free(dnc);
    free(scfg);
    return 1;
  }

  ncfg = dnc[0];
  for (i = 1; i < ns; i++) ncfg *= dnc[i];
  *cfg = (CONFIG *) malloc(sizeof(CONFIG)*ncfg);
  tmp = ncfg;
  p1 = dcfg + ns - 1;
  for (i = ns-1; i >= 0; i--) {
    tmp /= dnc[i];
    t = 0;
    while (t < ncfg) {
      for (j = 0; j < dnc[i]; j++) {
        for (k = 0; k < tmp; k++) {
          if (i == ns-1) {
            (*cfg)[t].n_shells = (*p1)[j].n_shells;
            size = sizeof(SHELL)*(*p1)[j].n_shells;
            (*cfg)[t].shells = (SHELL *) malloc(size);
            memcpy((*cfg)[t].shells, (*p1)[j].shells, size);
          } else {
            size_old = sizeof(SHELL)*(*cfg)[t].n_shells;
            size = sizeof(SHELL)*(*p1)[j].n_shells;
            (*cfg)[t].shells = (SHELL *) realloc((*cfg)[t].shells,
                                                 size_old+size);
            memcpy((*cfg)[t].shells+(*cfg)[t].n_shells,
                   (*p1)[j].shells, size);
            (*cfg)[t].n_shells += (*p1)[j].n_shells;
          }
          t++;
        }
      }
    }
    p1--;
  }

  for (i = 0; i < ns; i++) {
    for (j = 0; j < dnc[i]; j++) {
      if ((dcfg[i][j]).n_shells > 0) {
        free((dcfg[i][j]).shells);
      }
    }
    free(dcfg[i]);
  }
  free(dcfg);
  free(dnc);
  free(scfg);

  return ncfg;
}

/*
** FUNCTION:    GetConfigFromString
** PURPOSE:     construct all possible cofigurations from string.
** INPUT:       {CONFIG **cfg},
**              holds the resulting configurations.
**              {char *scfg},
**              string representation of the configuraitons.
** RETURN:      {int},
**              number of the resulting configurations.
** SIDE EFFECT:
** NOTE:
*/
int GetConfigFromString(const cfac_t *cfac, CONFIG **cfg, const char *scfg) {
  return GetConfigOrAverageFromString(cfac, cfg, NULL, scfg);
}

/*
** FUNCTION:    GetAverageConfigFromString
** PURPOSE:     construct the average configuration from a string.
** INPUT:       {int **n, **kappa, double **nq},
**              a list of principle quantum numbers, angular
**              quantum numbers, and the fractional occupation
**              numbers of the resulting average configuration.
**              {char *scfg},
**              string representation of the average configuration.
** RETURN:      {int},
**              number shells in the average configuration.
** SIDE EFFECT:
** NOTE:
*/
int GetAverageConfigFromString(const cfac_t *cfac, int **n, int **kappa,
                               double **nq, const char *scfg) {
  CONFIG *cfg;
  int i, ns;

  ns = GetConfigOrAverageFromString(cfac, &cfg, nq, scfg);
  if (ns <= 0) return ns;

  *n = (int *) malloc(sizeof(int)*ns);
  *kappa = (int *) malloc(sizeof(int)*ns);

  for (i = 0; i < ns; i++) {
    (*n)[i] = cfg->shells[i].n;
    (*kappa)[i] = cfg->shells[i].kappa;
  }

  free(cfg->shells);
  free(cfg);

  return ns;
}

/*
** FUNCTION:    GetJLFromSymbol
** PURPOSE:     decode the spectroscopic symbol for a shell
**              and retrieve the j and L values.
** INPUT:       {char *s},
**              the spectroscopic symbol.
**              {int *j},
**              holds the total angular momentum of the shell,
**              either +1 or -1, indicates whether it's *kl+1 or
**              *kl-1.
**              {int *kl},
**              holds the orbital angular momentum of the shell.
** RETURN:      {int},
**               0: success.
**              -1: the symobel unrecoginized.
** SIDE EFFECT:
** NOTE:        if the "+/-" sign is not present in the symbol,
**              the j-value returned is 0, which indicates either
**              value can be taken.
*/
int GetJLFromSymbol(char *s, int *j, int *kl) {
  int i;
  char s0[17], *p;

  strncpy(s0, s, 16);
  p = s0;
  while (*p) p++;
  p--;
  if ((*p) == '+') {
    if (j) *j = 1;
    *p = '\0';
  } else if ((*p) == '-') {
    if (j) *j = -1;
    *p = '\0';
  } else {
    if (j) *j = 0;
  }

  if (kl) {
    if (isdigit(s0[0])) *kl = atoi(s0);
    else {
      for (i = 0; i < MAX_SPEC_SYMBOLS; i++) {
        if (spec_symbols[i] == s0[0]) {
          *kl = i;
          return 0;
        }
      }
      return -1;
    }
  }
  return 0;
}

/*
** FUNCTION:    SpecSymbol
** PURPOSE:     construct the spectroscopic symbol for the
**              non-relativistic shell.
** INPUT:       {char *s},
**              string holding the result.
**              {int kl},
**              orbital angular momentum of the shell.
** RETURN:      {int},
**              always 0.
** SIDE EFFECT:
** NOTE:        if kl >= MAX_SPEC_SYMBOLS, then the symbol is
**              returned as "[kl]".
*/
int SpecSymbol(char *s, int kl) {
  if (kl < MAX_SPEC_SYMBOLS) {
    s[0] = spec_symbols[kl];
    s[1] = '\0';
  } else {
    sprintf(s, "[%d]", kl);
  }
  return 0;
}

static int CompareShellInvert(const void *ts1, const void *ts2) {
  return -CompareShell(ts1, ts2);
}

/*
** FUNCTION:    Couple
** PURPOSE:     recursively construct all possible states for a Config.
** INPUT:       {CONFIG *cfg},
**              pointer to the config. to be coupled.
** RETURN:      {int},
**               0: success.
**              <0: error.
** SIDE EFFECT:
** NOTE:
*/
int Couple(CONFIG *cfg) {
  CONFIG outmost, inner;
  int errcode;

  if (cfg->n_shells == 0) {
    errcode = -1;
    goto ERROR;
  }

  if (cfg == NULL) {
    errcode = -2;
    goto ERROR;
  }

  /* make sure that the shells are sorted in inverse order */
  qsort(cfg->shells, cfg->n_shells, sizeof(SHELL), CompareShellInvert);

  if (cfg->uta) {
    int i;
    cfg->csfs = NULL;
    cfg->n_csfs = 0;
    cfg->n_electrons = 0;
    for (i = 0; i < cfg->n_shells; i++) {
      cfg->n_electrons += cfg->shells[i].nq;
    }
    return 0;
  }

  if (cfg->n_shells == 1) {
    if (GetSingleShell(cfg) < 0) {
      errcode = -3;
      goto ERROR;
    }
  } else {
    outmost.uta = cfg->uta;
    outmost.n_shells = 1;
    outmost.shells = cfg->shells;
    inner.uta = cfg->uta;
    inner.n_shells = cfg->n_shells - 1;
    inner.shells = cfg->shells + 1;
    if (Couple(&outmost) < 0) {
      errcode = -4;
      goto ERROR;
    }
    if (Couple(&inner) < 0) {
      errcode = -5;
      free(outmost.csfs);
      goto ERROR;
    }

    if (CoupleOutmost(cfg, &outmost, &inner) < 0) {
      errcode = -5;
      free(outmost.csfs);
      free(inner.csfs);
      goto ERROR;
    }

    free(outmost.csfs);
    free(inner.csfs);
  }

  return 0;

 ERROR:
  printf("****Error in Couple****\n");
  return errcode;
}

/*
** FUNCTION:    CoupleOutmost
** PURPOSE:     constructs all possible states by coupling
**              the outmost shell to the inner shells.
** INPUT:       {CONFIG *cfg},
**              pointer to the resulting configuaration.
**              {CONFIG *outmost, *inner},
**              pointer to the configurations of the outmost shell
**              and the inner shells.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT:
** NOTE:        both outmost shell and inner shells must
**              have been coupled.
*/
int CoupleOutmost(CONFIG *cfg, CONFIG *outmost, CONFIG *inner) {
  int i, j, k;
  int bytes_csf, bytes_csf_inner, bytes_csf_outmost;
  int j2_min, j2_max;
  int j2_inner, j2_outmost;
  SHELL_STATE *csf, *csf_outmost, *csf_inner;

  if (outmost->n_shells != 1) goto ERROR;
  if (cfg->n_shells != 1 + inner->n_shells) goto ERROR;

  if (inner->n_shells == 0) {
    cfg->n_csfs = outmost->n_csfs;
    cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    memcpy(cfg->csfs, outmost->csfs, cfg->n_csfs * sizeof(SHELL_STATE));
    return 0;
  }

  bytes_csf_outmost = sizeof(SHELL_STATE);
  bytes_csf_inner = inner->n_shells * sizeof(SHELL_STATE);
  bytes_csf = bytes_csf_inner + bytes_csf_outmost;

  /*************************************************************************
  First calculte the total # of possible states, and allocate memory for
  cfg->csfs.
  *************************************************************************/
  csf_outmost = outmost->csfs;
  cfg->n_csfs = 0;
  for (i = 0; i < outmost->n_csfs; i++) {
    j2_outmost = csf_outmost->totalJ;
    csf_inner = inner->csfs;
    for (j = 0; j < inner->n_csfs; j++) {
      j2_inner = csf_inner->totalJ;
      j2_min = abs(j2_outmost - j2_inner);
      j2_max = j2_outmost + j2_inner;
      cfg->n_csfs += (j2_max - j2_min) / 2 + 1;
      csf_inner += inner->n_shells;
    }
    csf_outmost ++;
  }

  cfg->csfs = malloc(cfg->n_csfs * bytes_csf);
  if (cfg->csfs == NULL) goto ERROR;
  csf = cfg->csfs;

  /** Fill in the cfg->csfs array. **/
  csf_outmost = outmost->csfs;
  for (i = 0; i < outmost->n_csfs; i++) {
    j2_outmost = csf_outmost->totalJ;
    csf_inner = inner->csfs;
    for (j = 0; j < inner->n_csfs; j++) {
      j2_inner = csf_inner->totalJ;
      j2_min = abs(j2_outmost - j2_inner);
      j2_max = j2_outmost + j2_inner;
      for (k = j2_min; k <= j2_max; k += 2) {
        memcpy(csf, csf_outmost, bytes_csf_outmost);
        csf->totalJ = k;
        memcpy(csf + 1, csf_inner, bytes_csf_inner);
        csf += cfg->n_shells;
      }
      csf_inner += inner->n_shells;
    }
    csf_outmost++;
  }

  cfg->n_electrons = outmost->n_electrons + inner->n_electrons;

  return 0;

 ERROR:
  printf("****Error in CoupleOutmost****\n");
  return -1;
}

/*
** FUNCTION:    GetSingleShell
** PURPOSE:     construct all possible states for a single shell.
** INPUT:       {CONFIG *cfg},
**              pointer to the resulting configuration for the
**              single shell.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT:
** NOTE:        for j > 9/2, no more than 2 electrons are allowed.
**              The data were taken from "Nuclear Shell Theory" by
**              AMOS de-SHALIT and IGAL TALMI.
*/
int GetSingleShell(CONFIG *cfg) {
  int j2, max_q;
  int occupation;
  SHELL_STATE *csf;
  int i;

  if (cfg->n_shells != 1) goto ERROR;

  j2 = GetJ(cfg->shells);
  if (!(IsOdd(j2)) || j2 < 0) goto ERROR;

  max_q = j2 + 1;
  occupation = GetNq(cfg->shells);
  cfg->n_electrons = occupation;
  if ((2 * occupation) > max_q) {
    occupation = max_q - occupation;
  }
  if (occupation < 0) goto ERROR;

  switch(occupation) {
  case 0: /** 0 occupation or closed shell **/
    cfg->n_csfs = 1;
    cfg->csfs = malloc(sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    PackShellState(cfg->csfs, 0, 0, 0, 0);
    break;

  case 1: /** 1 electron **/
    cfg->n_csfs = 1;
    cfg->csfs = malloc(sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    PackShellState(cfg->csfs, j2, j2, 1, 0);
    break;

  case 2: /** 2 equivelent electrons **/
    cfg->n_csfs = (j2 + 1) / 2;
    cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
    if (cfg->csfs == NULL) goto ERROR;
    csf = cfg->csfs;
    PackShellState(csf++, 0, 0, 0, 0);
    for (i = 2; i < j2; i += 2) {
      PackShellState(csf++, i*2, i*2, 2, 0);
    }
    break;

  case 3: /** 3 equivelant electrons **/
    switch(j2) {
    case 5: /** j = 5/2 **/
      cfg->n_csfs = 3;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 5, 5, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      break;
    case 7: /** j = 7/2 **/
      cfg->n_csfs = 6;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 7, 7, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      break;
    case 9: /** j = 9/2 **/
      cfg->n_csfs = 10;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 9, 9, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 7, 7, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 13, 13, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      PackShellState(csf++, 17, 17, 3, 0);
      PackShellState(csf++, 21, 21, 3, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  case 4:
    switch(j2) {
    case 7: /** j = 7/2 **/
      cfg->n_csfs = 8;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 0, 0, 0, 0);
      PackShellState(csf++, 4, 4, 2, 0);
      PackShellState(csf++, 8, 8, 2, 0);
      PackShellState(csf++, 12, 12, 2, 0);
      PackShellState(csf++, 4, 4, 4, 0);
      PackShellState(csf++, 8, 8, 4, 0);
      PackShellState(csf++, 10, 10, 4, 0);
      PackShellState(csf++, 16, 16, 4, 0);
      break;
    case 9: /** j = 9/2 **/
      cfg->n_csfs = 18;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 0, 0, 0, 0);
      PackShellState(csf++, 4, 4, 2, 0);
      PackShellState(csf++, 8, 8, 2, 0);
      PackShellState(csf++, 12, 12, 2, 0);
      PackShellState(csf++, 16, 16, 2, 0);
      PackShellState(csf++, 0, 0, 4, 0);
      PackShellState(csf++, 4, 4, 4, 0);
      PackShellState(csf++, 6, 6, 4, 0);
      PackShellState(csf++, 8, 8, 4, 0);
      PackShellState(csf++, 8, 8, 4, 1);
      PackShellState(csf++, 10, 10, 4, 0);
      PackShellState(csf++, 12, 12, 4, 0);
      PackShellState(csf++, 12, 12, 4, 1);
      PackShellState(csf++, 14, 14, 4, 0);
      PackShellState(csf++, 16, 16, 4, 0);
      PackShellState(csf++, 18, 18, 4, 0);
      PackShellState(csf++, 20, 20, 4, 0);
      PackShellState(csf++, 24, 24, 4, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  case 5:
    switch(j2) {
    case 9: /** only j = 9/2 is allowed **/
      cfg->n_csfs = 20;
      cfg->csfs = malloc(cfg->n_csfs * sizeof(SHELL_STATE));
      if (cfg->csfs == NULL) goto ERROR;
      csf = cfg->csfs;
      PackShellState(csf++, 9, 9, 1, 0);
      PackShellState(csf++, 3, 3, 3, 0);
      PackShellState(csf++, 5, 5, 3, 0);
      PackShellState(csf++, 7, 7, 3, 0);
      PackShellState(csf++, 9, 9, 3, 0);
      PackShellState(csf++, 11, 11, 3, 0);
      PackShellState(csf++, 13, 13, 3, 0);
      PackShellState(csf++, 15, 15, 3, 0);
      PackShellState(csf++, 17, 17, 3, 0);
      PackShellState(csf++, 21, 21, 3, 0);
      PackShellState(csf++, 1, 1, 5, 0);
      PackShellState(csf++, 5, 5, 5, 0);
      PackShellState(csf++, 7, 7, 5, 0);
      PackShellState(csf++, 9, 9, 5, 0);
      PackShellState(csf++, 11, 11, 5, 0);
      PackShellState(csf++, 13, 13, 5, 0);
      PackShellState(csf++, 15, 15, 5, 0);
      PackShellState(csf++, 17, 17, 5, 0);
      PackShellState(csf++, 19, 19, 5, 0);
      PackShellState(csf++, 25, 25, 5, 0);
      break;
    default:
      goto ERROR;
    }
    break;

  default:
    goto ERROR;

  }

  return 0;

 ERROR:
  printf("****Error in GetSingleShell****\n");
  return -1;
}

/*
** FUNCTION:    PackShell, UnpackShell
** PURPOSE:     pack and unpack the fields of SHELL.
** INPUT:
** RETURN:
** SIDE EFFECT:
** NOTE:
*/
void UnpackShell(SHELL *s, int *n, int *kl, int *j, int *nq) {
  *n = s->n;
  *nq = s->nq;
  *j = 2*abs(s->kappa) - 1;
  *kl = (s->kappa < 0)? (*j - 1):(*j + 1);
}

void PackShell(SHELL *s, int n, int kl, int j, int nq){
  s->n = n;
  s->nq = nq;
  s->kappa = (kl - j)*(j + 1)/2;
}

void UnpackNRShell(int *s, int *n, int *kl, int *nq) {
  *nq = (*s)&0xFF;
  *kl = 2*(((*s)>>8)&0xFF);
  *n = ((*s)>>16)&0xFF;
}

void PackNRShell(int *s, int n, int kl, int nq) {
  *s = (n<<16) | ((kl/2)<<8) | nq;
}

/*
** FUNCTION:    GetNq, GetJ, GetL
** PURPOSE:     retrieve nq, j, and L of a shell.
** INPUT:
** RETURN:
** SIDE EFFECT:
** NOTE:
*/
int GetNq(SHELL *s){
  return s->nq;
}

int GetJ(SHELL *s){
  return 2*abs(s->kappa) - 1;
}

int GetL(SHELL *s){
  int j;
  j = 2*abs(s->kappa) - 1;
  return (s->kappa < 0)? (j - 1):(j + 1);
}

/*
** FUNCTION:    ShellClosed
** PURPOSE:     determine if the shell is a closed one.
** INPUT:
** RETURN:
** SIDE EFFECT:
** NOTE:
*/
int ShellClosed(SHELL *s) {
  int j;
  j = GetJ(s);
  if (s->nq < j+1) return 0;
  return 1;
}

/*
** FUNCTION:    GetLFromKappa, GetJFromKappa,
**              GetKappaFromJL, GetJLFromKappa
** PURPOSE:     convert between kappa and JL values.
** INPUT:
** RETURN:
** SIDE EFFECT:
** NOTE:
*/
int GetLFromKappa(int kappa) {
  int j;
  j = 2*abs(kappa) - 1;
  return (kappa < 0)? (j - 1):(j + 1);
}

int GetJFromKappa(int kappa) {
  return 2*abs(kappa) - 1;
}

int GetKappaFromJL(int j, int kl) {
  if (j <= 0 || kl < 0) return 0;
  return (kl-j)*(j+1)/2;
}

void GetJLFromKappa(int kappa, int *j, int *kl) {
  *j = 2*abs(kappa) - 1;
  *kl = (kappa < 0)? (*j - 1):(*j + 1);
}

/*
** FUNCTION:    PackShellState
** PURPOSE:     pack fields of SHELL_STATE to the structure.
** INPUT:
** RETURN:
** SIDE EFFECT:
** NOTE:
*/
void PackShellState(SHELL_STATE *s, int J, int j, int nu, int Nr){
  s->totalJ = J;
  s->shellJ = j;
  s->nu = nu;
  s->Nr = Nr;
}

/*
** FUNCTION:    GroupIndex
** PURPOSE:     find the index of the group by its name.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              the group index.
** SIDE EFFECT: if the group does not exist, a new one is created.
** NOTE:
*/
int GroupIndex(cfac_t *cfac, const char *name) {
    int gid = cfac_get_config_gid(cfac, name);

    if (gid < 0) {
        gid = AddGroup(cfac, name);
    }
    return gid;
}

/*
** FUNCTION:    GroupExists
** PURPOSE:     determine if a group exists.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              >=0: the group index, if it exists.
**               <0: the group does not exist.
** SIDE EFFECT:
** NOTE:
*/
int GroupExists(const cfac_t *cfac, const char *name) {
  int i;

  for (i = cfac->n_groups - 1; i >= 0; i--) {
    if (strncmp(name, cfac->cfg_groups[i].name, GROUP_NAME_LEN) == 0)
      break;
  }
  return i;
}

/*
** FUNCTION:    AddGroup
** PURPOSE:     add a group to the group array.
** INPUT:       {char *name},
**              the group name.
** RETURN:      {int},
**              the index of the added group
** SIDE EFFECT:
** NOTE:
*/
int AddGroup(cfac_t *cfac, const char *name) {
  if (name == NULL) return -1;
  if (cfac->n_groups == MAX_GROUPS) {
    printf("Max # groups reached\n");
    return -1;
  }
  strncpy(cfac->cfg_groups[cfac->n_groups].name, name, GROUP_NAME_LEN);
  cfac->n_groups++;
  return cfac->n_groups-1;
}

/*
** FUNCTION:    GetGroup
** PURPOSE:     retrieve the pointer to the group by its index.
** INPUT:       {int k},
**              the index of the group.
** RETURN:      {CONFIG_GROUP *},
**              the pointer to the group.
** SIDE EFFECT:
** NOTE:
*/
CONFIG_GROUP *GetGroup(const cfac_t *cfac, int k) {
  if (k < 0 || k >= cfac->n_groups) return NULL;
  return cfac->cfg_groups+k;
}

/*
** FUNCTION:    GetNumGroups
** PURPOSE:     retrieve the number of groups in the array.
** INPUT:
** RETURN:      {int},
**              the number of groups.
** SIDE EFFECT:
** NOTE:
*/
int GetNumGroups(const cfac_t *cfac) {
  return cfac->n_groups;
}

/*
** FUNCTION:    GetConfig
** PURPOSE:     return a pointer to CONFIG, which a state belongs to.
** INPUT:       {STATE *s},
**              pointer to a state.
** RETURN:      {CONFIG *},
**              the configuration the state belongs to.
** SIDE EFFECT:
** NOTE:
*/
CONFIG *GetConfig(const cfac_t *cfac, const STATE *s) {
  CONFIG *c;
  int i, j;

  i = s->kgroup;
  j = s->kcfg;
  c = ArrayGet(&(cfac->cfg_groups[i].cfg_list), j);
  return c;
}

static int configs_are_equal(const CONFIG *cfg1, const CONFIG *cfg2) {
    if (cfg1->n_csfs != cfg2->n_csfs) {
        return 0;
    }
    if (cfg1->n_shells != cfg2->n_shells) {
        return 0;
    }

    if (memcmp(cfg1->csfs, cfg2->csfs, cfg1->n_csfs*sizeof(SHELL_STATE)) != 0) {
        return 0;
    }
    if (memcmp(cfg1->shells, cfg2->shells, cfg1->n_shells*sizeof(SHELL)) != 0) {
        return 0;
    }

    return 1;
}

CONFIG *GetConfigFromGroup(const cfac_t *cfac, int kg, int kc) {
  if (kg < 0 || kg >= cfac->n_groups) {
    return NULL;
  } else {
    return (CONFIG *) ArrayGet(&(cfac->cfg_groups[kg].cfg_list), kc);
  }
}

/*
** FUNCTION:    AddConfigToList
** PURPOSE:     add a configuration to the specified group,
**              and add all states to the symmetry list.
** INPUT:       {int k},
**              the group index where the config is add to.
**              {CONFIG *cfg},
**              pointer to CONFIG to be added.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT:
** NOTE:
*/
int AddConfigToList(cfac_t *cfac, int k, CONFIG *cfg) {
  ARRAY *clist;
  int n0, kl0, nq0, m, i, n, kl, j, nq, ig;
  CONFIG_GROUP *cfgr;

  if (k < 0 || k >= cfac->n_groups) return -1;

  cfgr = &cfac->cfg_groups[k];

  if (cfgr->n_cfgs == 0) {
    cfgr->n_electrons = cfg->n_electrons;
  } else if (cfgr->n_electrons != cfg->n_electrons) {
    printf("Error: AddConfigToList, Configurations in a group ");
    printf("must have the same number of electrons\n");
    return -1;
  }
  clist = &(cfgr->cfg_list);

  cfg->energy = 0.0;
  cfg->delta = 0.0;

  n0 = 0;
  kl0 = -1;
  nq0 = 0;
  m = 0;
  cfg->nrs = malloc(sizeof(int)*cfg->n_shells);
  for (i = 0; i < cfg->n_shells; i++) {
    UnpackShell(cfg->shells+i, &n, &kl, &j, &nq);
    if (n == n0 && kl == kl0) {
      nq0 += nq;
    } else {
      if (nq0 > 0) {
        PackNRShell(cfg->nrs+m, n0, kl0, nq0);
        m++;
      }
      n0 = n;
      kl0 = kl;
      nq0 = nq;
    }
  }
  if (nq0 > 0) {
    PackNRShell(cfg->nrs+m, n0, kl0, nq0);
    m++;
  }

  cfg->nnrs = m;
  if (m < cfg->n_shells) {
    cfg->nrs = realloc(cfg->nrs, sizeof(int)*m);
  }

  /* Make sure that all configs are unique */
  for (ig = 0; ig < cfac->n_groups; ig++) {
    CONFIG_GROUP *tcfgr = &cfac->cfg_groups[ig];
    ARRAY *tclist;

    if (tcfgr->n_cfgs == 0 || tcfgr->n_electrons != cfg->n_electrons) {
      continue;
    }

    tclist = &(tcfgr->cfg_list);
    for (i = 0; i < tclist->dim; i++) {
      CONFIG *tcfg = ArrayGet(tclist, i);
      if (configs_are_equal(cfg, tcfg)) {
        printf("Error: overlapping configurations detected\n");
        return -1;
      }
    }
  }

  if (ArrayAppend(clist, cfg) == NULL) return -1;
  if (cfg->n_csfs > 0) {
    if (AddConfigToSymmetry(cfac, k, cfgr->n_cfgs, cfg) != 0) {
      return -1;
    }
  }
  cfgr->n_cfgs++;

  return 0;
}

/*
** FUNCTION:    AddStateToSymmetry
** PURPOSE:     add a state to the symmetry list.
** INPUT:       {int kg, kc, kstate},
**              the group index, configuration index, and state index
**              of the state.
**              {int parity, j},
**              parity and total angular momentum of the state.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT:
** NOTE:
*/
int AddStateToSymmetry(cfac_t *cfac, int kg, int kc, int kstate, int parity, int j) {
  int k;
  STATE s;
  ARRAY *st;

  k = IsEven(parity)? 2*j : (2*j+1);
  if (k >= MAX_SYMMETRIES) {
    printf("Maximum symmetry reached: %d %d\n", MAX_SYMMETRIES, k);
    return -1;
  }

  s.kgroup = kg;
  s.kcfg = kc;
  s.kstate = kstate;
  st = &(cfac->symmetry_list[k].states);
  if (ArrayAppend(st, &s) == NULL) return -1;
  cfac->symmetry_list[k].n_states++;
  return 0;
}

int ConfigParity(CONFIG *cfg) {
  int parity, i;

  parity = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    parity += (cfg->shells[i].nq)*GetL(&(cfg->shells[i]));
  }
  parity /= 2;
  parity = IsOdd(parity);

  return parity;
}

int PackSymState(int s, int k) {
  return s*100000 + k;
}

void UnpackSymState(int st, int *s, int *k) {
  if (s) *s = st/100000;
  if (k) *k = st%100000;
}

/*
** FUNCTION:    AddConfigToSymmetry
** PURPOSE:     add all states of a configuration to the symmetry list.
** INPUT:       {int kg, kc},
**              the group index and configuration index of the config.
**              a pointer of the configuration to be added.
** RETURN:      {int},
**               0: success.
**              -1: error.
** SIDE EFFECT:
** NOTE:
*/
static int AddConfigToSymmetry(cfac_t *cfac, int kg, int kc, CONFIG *cfg) {
  int parity;
  int i, j, k, m;
  STATE s;
  ARRAY *st;

  parity = 0;
  for (i = 0; i < cfg->n_shells; i++) {
    parity += (cfg->shells[i].nq)*GetL(&(cfg->shells[i]));
  }
  parity /= 2;
  for (m = 0; m < cfg->n_csfs; m++) {
    i = m*cfg->n_shells;
    j = (cfg->csfs)[i].totalJ;
    k = IsEven(parity)? 2*j : (2*j+1);
    if (k >= MAX_SYMMETRIES) {
      printf("Maximum symmetry reached: %d %d\n", MAX_SYMMETRIES, k);
      return -1;
    }
    s.kgroup = kg;
    s.kcfg = kc;
    s.kstate = i;
    st = &(cfac->symmetry_list[k].states);
    if (ArrayAppend(st, &s) == NULL) return -1;
    cfac->symmetry_list[k].n_states++;
  }
  return 0;
}

/*
** FUNCTION:    DecodePJ
** PURPOSE:     get the parity and J value from the symmetry index.
** INPUT:       {int i},
**              the symmetry index.
**              {int *p, int *j},
**              pointer holding the parity and J values.
** RETURN:
** SIDE EFFECT:
** NOTE:        if p or j is not required, pass in a NULL pointer.
*/
void DecodePJ(int i, int *p, int *j) {
  if (p) *p = IsOdd(i);
  if (j) *j = i/2;
}

/*
** FUNCTION:    GetSymmetry
** PURPOSE:     return a pointer to the symmetry by its index.
** INPUT:       {int k},
**              the symmetry index.
** RETURN:      {SYMMETRY *},
**              pointer to the returned symmetry.
** SIDE EFFECT:
** NOTE:
*/
SYMMETRY *GetSymmetry(const cfac_t *cfac, int k) {
  if (k < 0 || k >= MAX_SYMMETRIES) return NULL;
  return cfac->symmetry_list+k;
}

STATE *GetSymmetryState(SYMMETRY *sym, int isym)
{
  return ArrayGet(&(sym->states), isym);
}

int ShellIndex(int n, int kappa, int ns, SHELL *s) {
  int i;

  for (i = 0; i < ns; i++) {
    if (s[i].n == n && s[i].kappa == kappa) return i;
  }
  return -1;
}

int ShellToInt(int n, int kappa) {
  int k;
  k = (n-1)*(n-1) + 2*abs(kappa)- ((kappa>0)?1:2);
  return k;
}

void IntToShell(int i, int *n, int *kappa) {
  int k;

  *n = ((int) sqrt(i)) + 1 ;
  k = i - ((*n)-1)*((*n)-1) + 1;
  if (IsOdd(k)) *kappa = -(k+1)/2;
  else *kappa = k/2;
}

int ConstructConfigName(char *s, int n, CONFIG *c) {
  int i, j, k, m;
  char a[16], b[32], x;

  m = 0;
  s[0] = '\0';
  for (i = c->n_shells-1; i >= 0; i--) {
    GetJLFromKappa(c->shells[i].kappa, &j, &k);
    SpecSymbol(a, k/2);
    if (j > k) x = '+';
    else x = '-';
    if (i > 0) {
      sprintf(b, "%d%s%c%d ", c->shells[i].n, a, x, c->shells[i].nq);
    } else {
      sprintf(b, "%d%s%c%d", c->shells[i].n, a, x, c->shells[i].nq);
    }
    m += strlen(b);
    if (m >= n) return -1;
    strcat(s, b);
  }
  return m;
}

void ListConfig(const cfac_t *cfac, const char *fn, int n, int *kg) {
  int i, m, j;
  CONFIG *c;
  CONFIG_GROUP *g;
  char a[2048];
  FILE *f;

  if (fn == NULL || strcmp(fn, "-") == 0) f = stdout;
  else f = fopen(fn, "w");

  m = 0;
  for (i = 0; i < n; i++) {
    g = GetGroup(cfac, kg[i]);
    for (j = 0; j < g->n_cfgs; j++) {
      c = GetConfigFromGroup(cfac, kg[i], j);
      ConstructConfigName(a, 2048, c);
      fprintf(f, "%10s %6d   %s\n", g->name, m, a);
      m++;
    }
  }

  if (f != stdout) fclose(f);
}

/*
** FUNCTION:    MakeAverageConfig
** PURPOSE:     determine the average configuration based on given
**              groups, the weight given to each group, and possible
**              screening orbitals.
** INPUT:       {int ng},
**              number of groups which determines the average config.
**              {int *kg},
**              ng elements array of groups indexes.
**              {double *weight},
**              weight given for each group.
** RETURN:      {int},
**              >=0: success, the number of shells in the average config.
**               -1: error.
** SIDE EFFECT:
** NOTE:
*/
int MakeAverageConfig(cfac_t *cfac, int ng, int *kg, double *weight) {
  AVERAGE_CONFIG *acfg = &cfac->acfg;
  int n_screen = cfac->optimize_control.n_screen;
  int *screened_n = cfac->optimize_control.screened_n;
  double screened_charge = cfac->optimize_control.screened_charge;
  int screened_kl = cfac->optimize_control.screened_kl;

  double *tnq;
  int i, j, k, kmax, n, kappa, t;
  int weight_allocated = 0;

  if (ng <= 0) return -1;

  kmax = 0;
  for (i = 0; i < ng; i++) {
    ARRAY *c = &(cfac->cfg_groups[kg[i]].cfg_list);
    for (t = 0; t < cfac->cfg_groups[kg[i]].n_cfgs; t++) {
      CONFIG *cfg = ArrayGet(c, t);
      for (j = 0; j < cfg->n_shells; j++) {
        n = cfg->shells[j].n;
        kappa = cfg->shells[j].kappa;
        k = ShellToInt(n, kappa);
        if (k > kmax) {
          kmax = k;
        }
      }
    }
  }

  tnq = calloc(kmax + 1, sizeof(double));
  if (!tnq) return -1;

  if (weight == NULL) {
    weight = malloc(sizeof(double)*ng);
    if (!weight) {
      free(tnq);
      return -1;
    }
  }

  for (i = 0; i < ng; i++) {
    weight[i] = 1.0/ng;
  }
  weight_allocated = 1;

  for (i = 0; i < ng; i++) {
    ARRAY *c = &(cfac->cfg_groups[kg[i]].cfg_list);
    double a = 1.0/cfac->cfg_groups[kg[i]].n_cfgs;
    for (t = 0; t < cfac->cfg_groups[kg[i]].n_cfgs; t++) {
      CONFIG *cfg = ArrayGet(c, t);
      for (j = 0; j < cfg->n_shells; j++) {
        n = cfg->shells[j].n;
        kappa = cfg->shells[j].kappa;
        k = ShellToInt(n, kappa);
        tnq[k] += cfg->shells[j].nq*weight[i]*a;
      }
    }
    acfg->n_cfgs += cfac->cfg_groups[kg[i]].n_cfgs;
  }

  for (i = 0, j = 0; i <= kmax; i++) {
    if (tnq[i] > EPS10) {
      j++;
    }
  }

  if (!j) {
    if (weight_allocated) {
      free(weight);
    }
    free(tnq);
    return 0;
  }

  acfg->n_shells = j;
  acfg->n = malloc(sizeof(int)*j);
  acfg->kappa = malloc(sizeof(int)*j);
  acfg->nq = malloc(sizeof(double)*j);

  if (!acfg->n || !acfg->kappa || !acfg->nq) {
    if (acfg->n) free(acfg->n);
    if (acfg->nq) free(acfg->nq);
    if (acfg->kappa) free(acfg->kappa);
    if (weight_allocated) free(weight);
    free(tnq);
    return -1;
  }

  for (i = 0, j = 0; i <= kmax; i++) {
    if (tnq[i] > EPS10) {
      IntToShell(i, &n, &kappa);
      acfg->n[j] = n;
      acfg->kappa[j] = kappa;
      acfg->nq[j] = tnq[i];
      j++;
    }
  }

  /* add in configs for screened_charge */
  if (n_screen > 0) {
    screened_charge /= (double) n_screen;
    for (i = 0; i < n_screen; i++) {
      if (screened_kl < 0) {
        t = 0;
        kappa = -1;
      } else if (screened_kl == 0) {
        t = screened_n[i];
        kappa = GetKappaFromJL(t+1, t);
      } else {
        t = screened_n[i]*2-2;
        kappa = GetKappaFromJL(t+1, t);
      }
      for (j = 0; j < acfg->n_shells; j++) {
        k = GetLFromKappa(acfg->kappa[j]);
        if (acfg->n[j] < screened_n[i]) continue;
        if (acfg->n[j] > screened_n[i]) break;
        if (k > t) break;
        if (acfg->kappa[j] == kappa) break;
      }
      if (j < acfg->n_shells &&
          acfg->n[j] == screened_n[i] &&
          acfg->kappa[j] == kappa) {
        acfg->nq[j] += screened_charge;
      } else {
        acfg->n_shells += 1;
        acfg->n = realloc(acfg->n, sizeof(int)*acfg->n_shells);
        acfg->kappa = realloc(acfg->kappa, sizeof(int)*acfg->n_shells);
        acfg->nq = realloc(acfg->nq, sizeof(double)*acfg->n_shells);
        for (k = acfg->n_shells-1; k > j; k--) {
          acfg->n[k] = acfg->n[k-1];
          acfg->kappa[k] = acfg->kappa[k-1];
          acfg->nq[k] = acfg->nq[k-1];
        }
        acfg->n[j] = screened_n[i];
        acfg->kappa[j] = kappa;
        acfg->nq[j] = screened_charge;
      }
    }
  }

  if (weight_allocated) {
    free(weight);
  }

  free(tnq);

  return j;
}

int IBisect(int b, int n, int *a) {
  int i, i0, i1;

  if (n == 0) return -1;
  i0 = 0;
  i1 = n - 1;
  while (i1 - i0 > 1) {
    i = (i0 + i1)/2;
    if (b == a[i]) return i;
    else if (b < a[i]) i1 = i;
    else i0 = i;
  }

  if (b == a[i0]) return i0;
  else if (b == a[i1]) return i1;
  else return -1;
}

/*
** FUNCTION:    InGroups
** PURPOSE:     determine if a group is within a list of groups.
** INPUT:       {int kg},
**              the index of the group to be tested.
**              {int ng, *kgroup}
**              the number and index of the groups in the list.
** RETURN:      {int},
**              0: group kg is not in the list.
**              1: group kg is in the list.
** SIDE EFFECT:
** NOTE:
*/
int InGroups(int kg, int ng, const int *kgroup) {
  int i;
  if (ng < 0) {
    printf("negative ng passed to InGroups %d\n", ng);
    abort();
  }

  for (i = 0; i < ng; i++) {
    if (kg == kgroup[i]) return 1;
  }

  return 0;
}

/*
** FUNCTION:    CompareShell
** PURPOSE:     determine which of the two shells is the inner one.
** INPUT:       {const void *s1, *s2},
**              two shells in comparison.
** RETURN:      {int},
**              -1: s1 is inside s2.
**               0: s1 and s2 are the same.
**              +1: s1 is outside s2.
** SIDE EFFECT:
** NOTE:
*/
int CompareShell(const void *ts1, const void *ts2) {
  SHELL *s1, *s2;
  int ak1, ak2;

  s1 = (SHELL *) ts1;
  s2 = (SHELL *) ts2;
  if (s1->n > s2->n) return 1;
  else if (s1->n < s2->n) return -1;
  else {
    if (s1->kappa == s2->kappa) return 0;
    else {
      ak1 = abs(s1->kappa);
      ak2 = abs(s2->kappa);
      if (ak1 > ak2) return 1;
      else if (ak1 < ak2) return -1;
      else {
        if (s1->kappa < s2->kappa) return -1;
        else return 1;
      }
    }
  }
}

void InitConfigData(void *p, int n) {
  CONFIG *c;
  int k;

  c = (CONFIG *) p;
  for (k = 0; k < n; k++, c++) {
    c->n_shells = 0;
    c->n_csfs = 0;
    c->nnrs = 0;
  }
}

void FreeConfigData(void *p) {
  CONFIG *c;

  c = (CONFIG *) p;
  if (c->n_shells > 0) {
    free(c->shells);
    c->n_shells = 0;
  }
  if (c->n_csfs > 0) {
    free(c->csfs);
    c->n_csfs = 0;
  }
  if (c->nnrs > 0) {
    free(c->nrs);
    c->nnrs = 0;
  }
}

int cfac_get_config_gid(const cfac_t *cfac, const char *cname) {
    int i;

    for (i = cfac->n_groups - 1; i >= 0; i--) {
        if (!strncmp(cname, cfac->cfg_groups[i].name, GROUP_NAME_LEN)) {
            return i;
        }
    }

    return -1;
}

int cfac_add_config(cfac_t *cfac,
    const char *gname, const char *cfg_str, int uta)
{
    CONFIG *cfgs;
    int j, gidx, ncfgs;
    char *tmpstr;

    if (!gname || !cfg_str) {
        return -1;
    }

    gidx = GroupIndex(cfac, gname);
    if (gidx < 0) {
        return -1;
    }

    tmpstr = malloc((strlen(cfg_str) + 1)*sizeof(char));
    strcpy(tmpstr, cfg_str);

    ncfgs = GetConfigFromString(cfac, &cfgs, tmpstr);
    free(tmpstr);
    for (j = 0; j < ncfgs; j++) {
        CONFIG *cfg = &cfgs[j];
        cfg->uta = uta;

        if (Couple(cfg) < 0) {
            return -1;
        }
        if (AddConfigToList(cfac, gidx, cfg) < 0) {
            return -1;
        }
    }

    if (ncfgs > 0) {
        free(cfgs);
    }

    return gidx;
}
