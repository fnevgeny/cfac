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
#include <math.h>

#include <gsl/gsl_ieee_utils.h>

#include "cfacP.h"
#include "radial.h"
#include "angular.h"
#include "coulomb.h"
#include "grid.h"
#include "structure.h"
#include "transition.h"
#include "excitation.h"
#include "ionization.h"
#include "recombination.h"
#include "dbase.h"
#include "stoken.h"

static cfac_t *cfac = NULL;

static void cfac_verinfo(void) {
  printf("cFAC-%d.%d.%d\n",
    CFAC_VERSION, CFAC_SUBVERSION, CFAC_SUBSUBVERSION);
  printf("Based on the Flexible Atomic Code (FAC) by Ming Feng Gu\n");
  printf("Maintained by Evgeny Stambulchik\n");
}

static int PPrint(int argc, char *argv[], int argt[], ARRAY *variables) {
  int i;
  for (i = 0; i < argc; i++) {
    switch (argt[i]) {
    case NUMBER:
      printf("%s", argv[i]);
      break;
    case STRING:
      printf("%s", argv[i]);
      break;
    case LIST:
      printf("[%s]", argv[i]);
      break;
    case TUPLE:
      printf("(%s)", argv[i]);
      break;
    case KEYWORD:
      printf("%s = ", argv[i]);
    }
    if (i != argc-1 && argt[i] != KEYWORD) {
      printf(", ");
    }
  }
  if (argc > 0) printf("\n");
  fflush(stdout);
  return 0;
}

static int IntFromList(char *argv, int tp, ARRAY *variables, int **k) {
  int i;
  char *v[MAXNARGS];
  int t[MAXNARGS], n;

  if (tp == LIST) {
    n = DecodeArgs(argv, v, t, variables);
    if (n > 0) {
      *k = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
        (*k)[i] = atoi(v[i]);
        free(v[i]);
      }
    }
  } else {
    n = 1;
    *k = malloc(sizeof(int));
    (*k)[0] = atoi(argv);
  }

  return n;
}

static int DecodeGroupArgs(int **kg, int n, char *argv[], int argt[],
                           ARRAY *variables) {
  char *s;
  int i, k, ng;
  char *v[MAXNARGS];
  int t[MAXNARGS], nv;

  ng = n;
  nv = 0;
  if (ng > 0) {
    if (argt[0] == LIST || argt[0] == TUPLE) {
      if (ng > 1) {
        cfac_errmsg(cfac, "there should be only one list or tuple\n");
        return -1;
      }
      ng = DecodeArgs(argv[0], v, t, variables);
      nv = ng;
    } else {
      for (i = 0; i < ng; i++) {
        v[i] = argv[i];
        t[i] = argt[i];
      }
    }
    (*kg) = malloc(sizeof(int)*ng);
    for (i = 0; i < ng; i++) {
      if (t[i] != STRING) {
        cfac_errmsg(cfac, "argument must be a group name\n");
        free((*kg));
        return -1;
      }
      s = v[i];
      k = GroupExists(cfac, s);

      if (k < 0) {
        free((*kg));
        cfac_errmsg(cfac, "group does not exist\n");
        return -1;
      }

      (*kg)[i] = k;
    }
  } else {
    ng = GetNumGroups(cfac);
    (*kg) = malloc(sizeof(int)*ng);
    for (i = 0; i < ng; i++) (*kg)[i] = i;
  }

  for (i = 0; i < nv; i++) free(v[i]);

  return ng;
}

static int SelectNeleGroups(int nele, int **kg) {
    int i, j, ng = GetNumGroups(cfac);

    (*kg) = malloc(sizeof(int)*ng);
    if (!*kg) {
        return 0;
    }

    j = 0;
    for (i = 0; i < ng; i++) {
        CONFIG_GROUP *cg = GetGroup(cfac, i);
        if (cg->n_electrons == nele) {
            (*kg)[j++] = i;
        }
    }

    (*kg) = realloc(*kg, sizeof(int)*j);

    return j;
}

static int SelectLevels(int **t, char *argv, int argt, ARRAY *variables) {
  int n, ng, *kg, i, j, k, im, m, m0;
  int nrg, *krg, nrec;
  int ig, nlevels;
  LEVEL *lev;
  SYMMETRY *sym;
  STATE *s;
  char rgn[GROUP_NAME_LEN];
  char *v[MAXNARGS], *v1[MAXNARGS];
  int at[MAXNARGS], at1[MAXNARGS], nv, nv1, rv;

  if (argt != LIST  && argt != TUPLE) return -1;
  nv = 0;
  nv1 = 0;
  rv = 0;

  n = DecodeArgs(argv, v, at, variables);
  nv = n;
  if (n > 0) {
    if (at[0] == STRING) {
      ng = DecodeGroupArgs(&kg, n, v, at, variables);
      if (ng <= 0) {
        rv = -1;
        goto END;
      }
      nlevels = cfac_get_num_levels(cfac);
      (*t) = malloc(sizeof(int)*nlevels);
      k = 0;
      for (j = 0; j < nlevels; j++) {
        lev = GetLevel(cfac, j);

        if (lev->uta) {
          ig = lev->uta_cfg_g;
        } else {
          im = lev->pb;
          sym = GetSymmetry(cfac, lev->pj);
          s = GetSymmetryState(sym, im);
          ig = s->kgroup;
        }

        if (InGroups(ig, ng, kg)) {
          (*t)[k] = j;
          k++;
        }
      }
      free(kg);
      (*t) = realloc(*t, k*sizeof(int));
      rv = k;
      goto END;
    } else if (at[0] == LIST) {
      if (n != 2) {
        cfac_errmsg(cfac, "recombined states specification unrecognized\n");
        rv = -1;
        goto END;
      }
      ng = DecodeGroupArgs(&kg, 1, v, at, variables);
      if (ng <= 0) {
        rv = -1;
        goto END;
      }
      if (at[1] == LIST) {
        m0 = 0;
        n = DecodeArgs(v[1], v1, at1, variables);
        nv1 = n;
      } else if (at[1] == NUMBER) {
        m0 = 1;
        v1[1] = v[1];
        at1[1] = at[1];
      } else {
        cfac_errmsg(cfac, "Level specification unrecognized\n");
        rv = -1;
        goto END;
      }

      nrg = ng;
      krg = malloc(sizeof(int)*nrg);
      nlevels = cfac_get_num_levels(cfac);
      (*t) = malloc(sizeof(int)*nlevels);
      if (!(*t)) {
        rv = -1;
        goto END;
      }
      k = 0;
      for (m = m0; m < n; m++) {
        if (at1[m] != NUMBER) {
          rv = -1;
          goto END;
        }
        nrec = atoi(v1[m]);
        for (i = 0; i < nrg; i++) {
          ConstructRecGroupName(rgn, GetGroup(cfac, kg[i])->name, nrec);
          krg[i] = GroupExists(cfac, rgn);
        }
        for (j = 0; j < nlevels; j++) {
          lev = GetLevel(cfac, j);
          im = lev->pb;
          sym = GetSymmetry(cfac, lev->pj);
          s = GetSymmetryState(sym, im);
          ig = s->kgroup;
          if (ig < 0) {
            if (!ValidBasis(cfac, s, ng, kg, nrec)) continue;
            (*t)[k] = j;
            k++;
          } else {
            if (InGroups(ig, nrg, krg)) {
              (*t)[k] = j;
              k++;
            }
          }
        }
      }

      free(krg);
      free(kg);
      (*t) = realloc(*t, k*sizeof(int));
      rv = k;
      goto END;
    } else {
      (*t) = malloc(sizeof(int)*n);
      for (i = 0; i < n; i++) {
        if (at[i] != NUMBER) return -1;
        (*t)[i] = atoi(v[i]);
      }
      rv = n;
      goto END;
    }
  }

 END:
  for (i = 0; i < nv; i++) free(v[i]);
  for (i = 0; i < nv1; i++) free(v1[i]);
  return rv;
}

static int SelectNeleLevels(cfac_t *cfac, int nele, int **levels)
{
    int i, j, n;

    n = cfac_get_num_levels(cfac);
    *levels = malloc(n*sizeof(int));
    if (!(*levels)) {
        return -1;
    }

    j = 0;
    for (i = 0; i < n; i++) {
        if (nele < 0 || GetNumElectrons(cfac, i) == nele) {
            (*levels)[j++] = i;
        }
    }

    return j;
}

static int PSetErrorOutput(int argc, char *argv[], int argt[], ARRAY *variables)
{
    char *fname;
    FILE *fp;

    if (argc != 1 || argt[0] != STRING) return -1;

    fname = argv[0];

    if (!strcmp(fname, "stdout")) {
        fp = stdout;
    } else
    if (!strcmp(fname, "stderr")) {
        fp = stderr;
    } else {
        fp = fopen(fname, "w");
        if (!fp) {
            cfac_errmsg(cfac, "Cannot open file '%s' for writing\n", fname);
            return -1;
        }
    }

    cfac_set_err_fp(cfac, fp);
    return 0;
}

static int PAvgConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  int ns, *n, *kappa;
  double *nq;

  if (argc != 1 || argt[0] != STRING) return -1;

  ns = GetAverageConfigFromString(cfac, &n, &kappa, &nq, argv[0]);
  if (ns <= 0) return -1;

  if (SetAverageConfig(cfac, ns, n, kappa, nq) < 0) return -1;

  free(n);
  free(kappa);
  free(nq);
  return 0;
}

static int PCheckEndian(int argc, char *argv[], int argt[], ARRAY *variables) {
  FILE *f;
  F_HEADER fh;
  int i, swp;

  if (argc == 0) {
    i = CheckEndian(NULL);
  } else {
    f = fopen(argv[0], "rb");
    if (f == NULL) {
      cfac_errmsg(cfac, "Cannot open file '%s'\n", argv[0]);
      return -1;
    }
    ReadFHeader(f, &fh, &swp);
    i = CheckEndian(&fh);
    fclose(f);
  }

  printf("Endian: %d\n", i);

  return 0;
}

static char _closed_shells[MCHSHELL] = "";
static int PClosed(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  int i, j, kl, n, nq, ncfg;
  char *p;
  char s[8], st[16];
  int ns, k;

  if (argc == 0) _closed_shells[0] = '\0';
  for (i = 0; i < argc; i++) {
    if (argt[i] != STRING) return -1;
    ns = StrSplit(argv[i], ' ');
    p = argv[i];
    for (k = 0; k < ns; k++) {
      while (*p == ' ') p++;
      ncfg = GetConfigFromStringNR(cfac, &cfg, p);
      for (j = ncfg-1; j >= 0; j--) {
        if (cfg[j].n_shells != 1) return -1;
        n = (cfg[j].shells)[0].n;
        kl = (cfg[j].shells)[0].kappa;
        nq = 2*(kl + 1);
        kl = kl/2;
        SpecSymbol(s, kl);
        sprintf(st, "%d%s%d ", n, s, nq);
        strcat(_closed_shells, st);
        free(cfg[j].shells);
      }
      if (ncfg > 0) free(cfg);
      while (*p) p++;
      p++;
    }
  }
  return 0;
}

static int PGetConfigNR(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  int i, j, t, ncfg;
  char scfg[MCHSHELL], s[16];

  for (i = 0; i < argc; i++) {
    if (argt[i] != STRING) return -1;
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, argv[i], MCHSHELL - 1);
    ncfg = GetConfigFromStringNR(cfac, &cfg, scfg);
    for (j = 0; j < ncfg; j++) {
      scfg[0] = '\0';
      for (t = cfg[j].n_shells-1; t >= 0; t--) {
        sprintf(s, "%d", (cfg[j].shells)[t].n);
        strcat(scfg, s);
        SpecSymbol(s, (cfg[j].shells)[t].kappa/2);
        strcat(scfg, s);
        if (t == 0) {
          sprintf(s, "%d", (cfg[j].shells)[t].nq);
        } else {
          sprintf(s, "%d ", (cfg[j].shells)[t].nq);
        }
        strcat(scfg, s);
      }
      printf("%s\n", scfg);
    }
    if (ncfg > 0) free(cfg);
  }

  return 0;
}

static int PConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  CONFIG *cfg;
  static char gname[GROUP_NAME_LEN + 1] = "_all_";
  int i, j, k, t, ncfg;
  char scfg[MCHSHELL];
  int iuta = cfac->uta;

  k = -2;
  for (i = 0; i < argc; i++) {
    if (argt[i] == KEYWORD) {
      if (strcmp(argv[i], "group") != 0) {
        cfac_errmsg(cfac, "The keyword must be group=gname\n");
        return -1;
      }
      if (i > argc-2) return -1;
      if (argt[i+1] != STRING) return -1;
      k = i;
    }
  }

  i = 0;

  if (k >= 0) {
    strncpy(gname, argv[k+1], GROUP_NAME_LEN);
  } else {
    if (argc == 0) {
      return -1;
    }

    if (argt[i] == NUMBER && argc >= 3) {
      iuta = atoi(argv[i]) ? 1:0;
      i++;
    }

    if (argt[i] != STRING) {
      return -1;
    }
    strncpy(gname, argv[i], GROUP_NAME_LEN);

    i++;
  }

  for (; i < argc; i++) {
    if (i == k || i == k+1) continue;
    if (argt[i] != STRING) return -1;
    strncpy(scfg, _closed_shells, MCHSHELL);
    strncat(scfg, argv[i], MCHSHELL - 1);
    ncfg = GetConfigFromString(cfac, &cfg, scfg);
    for (j = 0; j < ncfg; j++) {
      cfg[j].uta = iuta;
      if (Couple(&cfg[j]) < 0) return -1;
      t = GroupIndex(cfac, gname);
      if (t < 0) return -1;
      if (AddConfigToList(cfac, t, &cfg[j]) < 0) return -1;
    }
    if (ncfg > 0) free(cfg);
  }

  return 0;
}

static int PListConfig(int argc, char *argv[], int argt[], ARRAY *variables) {
  int k, ng, *kg;
  char *s;

  s = NULL;
  ng = 0;
  if (argc > 0) {
    s = argv[0];
    if (argc > 1) {
      ng = DecodeGroupArgs(&kg, 1, argv+1, argt+1, variables);
    }
  }
  if (ng <= 0) {
    ng = GetNumGroups(cfac);
    kg = malloc(sizeof(int)*ng);
    for (k = 0; k < ng; k++) {
      kg[k] = k;
    }
  }

  ListConfig(cfac, s, ng, kg);

  if (ng > 0) free(kg);

  return 0;
}

static int PConfigEnergy(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int m, mr, i;
  int ng, *kg;

  if (argc == 0) return -1;
  m = atoi(argv[0]);

  if (argc == 1 || m != 0) {
    if (ConfigEnergy(cfac, m, 0, 0, NULL) < 0) {
      return -1;
    }
  } else {
    mr = atoi(argv[1]);
    if (argc == 2) {
      if (ConfigEnergy(cfac, m, mr, 0, NULL) < 0) {
        return -1;
      }
    } else {
      for (i = 1; i < argc; i++) {
        ng = DecodeGroupArgs(&kg, 1, argv+i, argt+i, variables);
        if (ng < 0) return -1;
        if (ConfigEnergy(cfac, m, mr, ng, kg) < 0) {
          if (ng > 0) free(kg);
          return -1;
        }
        if (ng > 0) free(kg);
      }
    }
  }

  return 0;
}

static int PAITable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow = 0, *low = NULL, nup = 0, *up = NULL;
  char *fname;
  int nele;

  if (argc < 1 || argt[0] != STRING) {
    return -1;
  }
  fname = argv[0];

  switch (argc) {
  case 3:
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    nup  = SelectLevels(&up, argv[2], argt[2], variables);
    break;
  case 2:
    nele = atoi(argv[1]);
    if (nele < 1) {
      return -1;
    }
    nlow = SelectNeleLevels(cfac, nele, &low);
    nup = SelectNeleLevels(cfac, nele - 1, &up);
    break;
  default:
    return -1;
    break;
  }

  if (nlow == 0 || nup == 0) {
    cfac_errmsg(cfac, "Empty set of initial or final levels in AITable(), skipping\n");
    return 0;
  } else
  if (nlow < 0 || nup < 0) {
    return -1;
  }

  SaveAI(cfac, nlow, low, nup, up, fname, 0);

  if (nlow > 0) {
    free(low);
  }
  if (nup > 0) {
    free(up);
  }

  return 0;
}

static int PAITableMSub(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlow, *low, nup, *up;

  if (argc != 3) return -1;
  if (argt[0] != STRING) return -1;

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  SaveAI(cfac, nlow, low, nup, up, argv[0], 1);
  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PBasisTable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int m;

  if (argc == 0) return -1;
  if (argc > 2 || argt[0] != STRING) return -1;
  if (argc == 2) m = atoi(argv[1]);
  else m = 0;

  GetBasisTable(cfac, argv[0], m);

  return 0;
}

static int PCETable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  if (argc < 2 || argc > 3) {
    return -1;
  }
  if (argt[0] != STRING ) {
    return -1;
  }

  if (argc == 2) {
    if (argt[1] == STRING) {
      nlow = SelectLevels(&low, argv[1], argt[1], variables);
    } else {
      int nele = atoi(argv[1]);
      nlow = SelectNeleLevels(cfac, nele, &low);
      nup = nlow;
      up = low;
    }
    if (nlow <= 0) return -1;
    SaveExcitation(cfac, nlow, low, nlow, low, 0, argv[0]);
    free(low);
  } else if (argc == 3) {
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    SaveExcitation(cfac, nlow, low, nup, up, 0, argv[0]);
    free(low);
    free(up);
  } else {
    return -1;
  }

  return 0;
}

static int PCETableMSub(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  if (argc == 1) {
    if (argt[0] != STRING ) return -1;
    SaveExcitation(cfac, nlow, low, nup, up, 1, argv[0]);
  } else if (argc == 2) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    SaveExcitation(cfac, nlow, low, nlow, low, 1, argv[0]);
    free(low);
  } else if (argc == 3) {
    if (argt[0] != STRING) return -1;
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    if (nlow <= 0) return -1;
    nup = SelectLevels(&up, argv[2], argt[2], variables);
    if (nup <= 0) return -1;
    SaveExcitation(cfac, nlow, low, nup, up, 1, argv[0]);
    free(low);
    free(up);
  } else {
    return -1;
  }

  return 0;
}

static int PCITable(int argc, char *argv[], int argt[], ARRAY *variables) {
  int nlow, *low, nup, *up;

  if (argc < 2 || argc > 3) return -1;
  if (argt[0] != STRING) return -1;

  if (argc == 2) {
    if (argt[1] == NUMBER) {
      int nele = atoi(argv[1]);
      if (nele < 1) {
        return -1;
      }
      nlow = SelectNeleLevels(cfac, nele, &low);
      nup = SelectNeleLevels(cfac, nele - 1, &up);
    } else {
      return -1;
    }
  } else {
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    nup = SelectLevels(&up, argv[2], argt[2], variables);
  }

  if (nlow == 0 || nup == 0) {
    cfac_errmsg(cfac, "Empty set of initial or final levels in CITable(), skipping\n");
    return 0;
  } else
  if (nlow < 0 || nup < 0) {
    return -1;
  }

  if (SaveIonization(cfac, nlow, low, nup, up, argv[0]) < 0) return -1;

  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PCITableMSub(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlow, *low, nup, *up;

  if (argc != 3) return -1;
  if (argt[0] != STRING) return -1;
  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (SaveIonizationMSub(cfac, nlow, low, nup, up, argv[0]) < 0) return -1;

  if (nlow > 0) free(low);
  if (nup > 0) free(up);

  return 0;
}

static int PCorrectEnergy(int argc, char *argv[], int argt[],
                          ARRAY *variables) {
  if (argc == 2) {
    FILE *f;
    int nele;
    char *refname = NULL;
    if (argt[0] != STRING || argt[1] != NUMBER) {
      return -1;
    }

    nele = atoi(argv[1]);

    f = fopen(argv[0], "r");
    if (!f) {
      cfac_errmsg(cfac, "Cannot open file '%s'\n", argv[0]);
      return -1;
    }

    while (1) {
      double e;
      char strbuf[1024], *q, *name;
      if (fgets(strbuf, 1024, f) == NULL) {
        break;
      }

      if (strbuf[0] == '#') {
        continue;
      }

      if ((q = strchr(strbuf, ',')) == NULL) {
        cfac_errmsg(cfac, "Failed parsing line '%s'\n", strbuf);
        return -1;
      }

      *q = '\0';
      name = strbuf;

      int nf = sscanf(q + 1, "%lf", &e);
      if (nf != 1) {
        cfac_errmsg(cfac, "Failed parsing line '%s'\n", strbuf);
        return -1;
      }

      if (refname == NULL) {
        refname = strdup(name);
      }

      e /= HARTREE_EV;
      AddECorrection(cfac, nele, name, refname, e);
    }
    if (refname) {
      free(refname);
    }
    fclose(f);
  } else {
    return -1;
  }

  return 0;
}

static int PExit(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  exit(0);
}

static int PGetPotential(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  if (argc != 1 || argt[0] != STRING) return -1;
  GetPotential(cfac, argv[0]);
  return 0;
}

static int PInfo(int argc, char *argv[], int argt[], ARRAY *variables) {
  if (argc != 0) return -1;
  cfac_verinfo();
  return 0;
}

static int PMemENTable(int argc, char *argv[], int argt[],
                       ARRAY *variables) {

  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  MemENTable(cfac, argv[0]);

  return 0;
}

static int POptimizeRadial(int argc, char *argv[], int argt[],
                           ARRAY *variables) {
  int ng, i, k;
  int *kg;
  double z;
  double *weight;
  char *vw[MAXNARGS];
  int iw[MAXNARGS], ni;

  ni = 0;

  ng = argc;
  if (ng == 0) {
    ng = 0;
    kg = NULL;
    weight = NULL;
    goto END;
  }

  if (argt[0] == STRING) {
    weight = NULL;
    ng = DecodeGroupArgs(&kg, argc, argv, argt, variables);
    if (ng < 0) return -1;
  } else {
    ng = DecodeGroupArgs(&kg, 1, argv, argt, variables);
    if (ng < 0) return -1;

    if (argc == 1) {
      weight = NULL;
    } else {
      if (argt[1] != LIST && argt[1] != TUPLE) return -1;
      k = DecodeArgs(argv[1], vw, iw, variables);
      ni = k;
      if (k < 0 || k > ng) {
        cfac_errmsg(cfac, "weights must be a sequence\n");
        return -1;
      }
      weight = malloc(sizeof(double)*ng);
      z = 0.0;
      for (i = 0; i < k; i++) {
        if (iw[i] != NUMBER) {
          return -1;
        }
        weight[i] = atof(vw[i]);
        z += weight[i];
      }
      for (i = k; i < ng; i++) {
        if (z >= 1.0) {
          weight[i] = weight[k-1];
        } else {
          weight[i] = (1.0-z)/(ng-k);
        }
      }
    }
  }

 END:
  if (OptimizeRadial(cfac, ng, kg, weight) < 0) {
    if (kg) free(kg);
    if (weight) free(weight);
    return -1;
  }
  if (weight) free(weight);
  if (kg) free(kg);

  for (i = 0; i < ni; i++) free(vw[i]);

  return 0;
}

static int PRefineRadial(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  int maxfun, msglvl;

  maxfun = 0;
  msglvl = 0;
  if (argc > 0) {
    maxfun = atoi(argv[0]);
    if (argc > 1) {
      msglvl = atoi(argv[1]);
    }
  }

  return RefineRadial(cfac, maxfun, msglvl);
}

static int PPause(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  char s[10];

  while (1) {
    printf("Type go to continue: ");
    if (scanf("%s", s) == 1 && strcmp(s, "go") == 0) break;
  }

  return 0;
}

static int PPrintTable(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  int v;

  if (argc != 2 && argc != 3) return -1;
  if (argt[0] != STRING || argt[1] != STRING) return -1;

  v = 1;
  if (argc == 3) {
    if (argt[2] != NUMBER) return -1;
    v = atoi(argv[2]);
  }

  PrintTable(cfac, argv[0], argv[1], v);
  return 0;
}

static long unsigned int sid = 0;
static sqlite3 *db = NULL;
static char *cmdline = NULL;

static int PStoreInit(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  int reset = 0;
  if (argc != 1 && argc != 2) return -1;

  if (sid) {
    cfac_errmsg(cfac, "Store has already been initialized\n");
    return -1;
  }

  if (argc == 2) {
    if (argt[1] != NUMBER) {
      return -1;
    } else {
      reset = atoi(argv[1]);
    }
  }

  return StoreInit(cfac, argv[0], reset, &db, &sid);
}

static int PStoreTable(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  if (argc != 1) return -1;
  if (argt[0] != STRING) return -1;

  if (!sid) {
    cfac_errmsg(cfac, "Store has not been initialized yet\n");
    return -1;
  }

  return StoreTable(cfac, db, sid, argv[0]);
}

static int PStoreClose(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  if (argc != 0) return -1;

  if (!sid) {
    cfac_errmsg(cfac, "Store has not been initialized yet\n");
    return -1;
  }

  StoreClose(cfac, db, sid, cmdline);

  sid = 0;
  return 0;
}

static int PRecStates(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int ng, *kg, n;

  if (argc != 3 || argt[0] != STRING ||
      (argt[1] != LIST && argt[1] != TUPLE) ||
      argt[2] != NUMBER)
    return -1;
  ng = DecodeGroupArgs(&kg, 1, &(argv[1]), &(argt[1]), variables);
  if (ng <= 0) return -1;
  n = atoi(argv[2]);
  if (RecStates(cfac, n, ng, kg, argv[0]) < 0) {
    cfac_errmsg(cfac, "RecStates Error\n");
    free(kg);
    return -1;
  }

  free(kg);

  return 0;
}

static int PRRMultipole(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlow, *low, nup, *up, m;

  m = -1;
  if (argc != 3 && argc != 4) return -1;
  if (argt[0] != STRING) return -1;

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;
  if (argc == 4) {
    if (argt[3] != NUMBER) return -1;
    m = atoi(argv[3]);
  }
  SaveRRMultipole(cfac, nlow, low, nup, up, argv[0], m);

  free(low);
  free(up);

  return 0;
}

static int PRRTable(int argc, char *argv[], int argt[],
                    ARRAY *variables) {
  int nlow, *low, nup, *up, m;

  m = -1;

  if (argc < 2 || argc > 4 || argt[0] != STRING) {
    return -1;
  }

  if (argt[1] == NUMBER) {
    int nele = atoi(argv[1]);
    if (nele < 1) {
      return -1;
    }

    if (argc > 3) {
      return -1;
    }

    if (argc == 3) {
      if (argt[2] == NUMBER) {
        m = atoi(argv[2]);
      } else {
        return -1;
      }
    }

    nlow = SelectNeleLevels(cfac, nele, &low);
    nup = SelectNeleLevels(cfac, nele - 1, &up);
  } else {
    if (argc < 3) {
      return -1;
    }

    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    nup = SelectLevels(&up, argv[2], argt[2], variables);

    if (argc == 4) {
      if (argt[3] != NUMBER) return -1;
      m = atoi(argv[3]);
    }
  }

  if (nlow == 0 || nup == 0) {
    cfac_errmsg(cfac, "Empty set of initial or final levels in RRTable(), skipping\n");
    return 0;
  } else
  if (nlow < 0 || nup < 0) {
    return -1;
  }

  SaveRecRR(cfac, nlow, low, nup, up, argv[0], m);

  if (nlow) {
    free(low);
  }
  if (nup) {
    free(up);
  }

  return 0;
}

static int PSetAICut(int argc, char *argv[], int argt[],
                     ARRAY *variables) {
  double c;

  if (argc != 1 || argt[0] != NUMBER) return -1;

  c = atof(argv[0]);
  SetAICut(c);

  return 0;
}

static int PSetAngZOptions(int argc, char *argv[], int argt[],
                           ARRAY *variables) {
  int n;
  double c, mc;

  n = atoi(argv[0]);
  c = EPS3;
  mc = EPS3;
  if (argc > 1) {
    mc = atof(argv[1]);
    if (argc > 2) {
      c = atof(argv[2]);
    }
  }
  SetAngZOptions(cfac, n, mc, c);

  return 0;
}

static int PSetCILevel(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  int i;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  i = atoi(argv[0]);
  SetCILevel(cfac, i);

  return 0;
}

static int PSetAngZCut(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  double c;

  if (argc != 1 || argt[0] != NUMBER) return -1;
  c = atof(argv[0]);
  SetAngZCut(cfac, c);

  return 0;
}

static int PSetMixCut(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  double c, c2;

  if (argc < 1 || argc > 2) return -1;
  if (argt[0] != NUMBER) return -1;
  c = atof(argv[0]);
  c2 = -1.0;
  if (argc > 1) {
    if (argt[1] != NUMBER) return -1;
    c2 = atof(argv[1]);
  }
  SetMixCut(cfac, c, c2);

  return 0;
}

static int PSetAtom(int argc, char *argv[], int argt[],
                    ARRAY *variables) {
  double mass, rn;
  int z;

  mass = 0.0;
  z = 0;
  rn = -1.0;

  if (argc < 1 || argt[0] != STRING || argc > 4) return -1;
  if (argc > 1) {
    z = atoi(argv[1]);
    if (argc > 2) {
      mass = atof(argv[2]);
      if (argc > 3) {
        rn = atof(argv[3]);
      }
    }
  }

  if (cfac_set_atom(cfac, argv[0], z, mass, rn) < 0) return -1;

  return 0;
}

static int PSetAvgConfig(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int ns, i, j;
  int *n, *kappa;
  double *nq;
  char *vc[MAXNARGS], *vs[MAXNARGS];
  int ic[MAXNARGS], is[MAXNARGS];

  if (argc != 1 || argt[0] != LIST) {
    return -1;
  }

  ns = DecodeArgs(argv[0], vc, ic, variables);

  n = malloc(sizeof(int)*ns);
  kappa = malloc(sizeof(int)*ns);
  nq = malloc(sizeof(double)*ns);

  for (i = 0; i < ns; i++) {
    if (DecodeArgs(vc[i], vs, is, variables) != 4) {
      return -1;
    }
    n[i] = atoi(vs[0]);
    kappa[i] = atoi(vs[1]);
    if (atoi(vs[2]) > 0) kappa[i] = -(kappa[i]+1);
    nq[i] = atof(vs[3]);
    for (j = 0; j < 4; j++) free(vs[j]);
  }

  for (j = 0; j < ns; j++) free(vc[j]);

  if (SetAverageConfig(cfac, ns, n, kappa, nq) < 0) return -1;
  free(n);
  free(kappa);
  free(nq);

  return 0;
}

static int PSetCEGrid(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCEEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        free(vg[i]);
        xg[i] /= HARTREE_EV;
      }
      err = SetCEEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCEEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetAngleGrid(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int n, ng, i, err, m;
  double xg[MAXNTHETA+MAXNPHI];
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;

  if (n == 2) {
    m = atoi(argv[0]);
    if (argt[1] == NUMBER) {
      ng = atoi(argv[1]);
      if (m == 0) {
        emin = 0.0;
        emax = M_PI;
      } else {
        emin = 0.0;
        emax = TWO_PI;
      }
      err = SetAngleGrid(m, ng, emin, emax);
    } else if (argt[1] == LIST || argt[1] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        free(vg[i]);
        xg[i] /= HARTREE_EV;
      }
      err = SetAngleGridDetail(cfac, m, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 4) {
    m = atoi(argv[0]);
    ng = atoi(argv[1]);
    emin = atof(argv[2]);
    emax = atof(argv[3]);
    emin *= M_PI/180.0;
    emax *= M_PI/180.0;
    err = SetAngleGrid(m, ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetCEGridLimits(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;

  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
        type = atoi(argv[2]);
      }
    }
  }

  SetCEEGridLimits(emin, emax, type);

  return 0;
}

static int PSetCEGridType(int argc, char *argv[], int argt[],
                          ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetCEEGridType(atoi(argv[0]));
  return 0;
}

static int PSetTEGrid(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCETEGrid(cfac, ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        if (ig[i] != NUMBER) return -1;
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetCETEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetCETEGrid(cfac, ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetCEPWOptions(int argc, char *argv[], int argt[],
                           ARRAY *variables) {
  int qr, max, kl_cb;

  qr = EXCLQR;
  max = EXCLMAX;
  kl_cb = EXCLCB;

  if (argc < 1 || argc > 3) return -1;

  max = atoi(argv[0]);
  if (argc > 1) {
    qr = atoi(argv[1]);
    if (argc > 2) {
      kl_cb = atoi(argv[2]);
    }
  }

  return SetCEPWOptions(cfac, qr, max, kl_cb);
}

static int PSetCEPWGridType(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetCEPWGridType(atoi(argv[0]));

  return 0;
}

static int PSetCEPWGrid(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int ns, i;
  int n;
  int *m, *step;
  char *v1[MAXNARGS], *v2[MAXNARGS];
  int t1[MAXNARGS], t2[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] != NUMBER) return -1;
    ns = atoi(argv[0]);
    SetCEPWGrid(cfac, -ns, NULL, NULL);
  } else {
    if (argt[0] != LIST || argt[1] != LIST) return -1;
    ns = DecodeArgs(argv[0], v1, t1, variables);
    if (ns <= 0) return -1;
    if (ns != DecodeArgs(argv[1], v2, t2, variables)) return -1;
    m = malloc(ns*sizeof(int));
    step = malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      if (t1[i] != NUMBER || t2[i] != NUMBER) return -1;
      m[i] = atoi(v1[i]);
      step[i] = atoi(v2[i]);
      free(v1[i]);
      free(v2[i]);
    }
    SetCEPWGrid(cfac, ns, m, step);
    free(m);
    free(step);
  }

  return 0;
}

static int PSetCEBorn(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  double eb, x, x1, x0;

  if (argc < 1 || argc > 4) return -1;
  if (argt[0] != NUMBER) return -1;
  x0 = XBORN0;
  x1 = XBORN1;
  x = XBORN;
  if (argc > 1) {
    if (argt[1] != NUMBER) return -1;
    x = atof(argv[1]);
    if (argc > 2) {
      if (argt[2] != NUMBER) return -1;
      x1 = atof(argv[2]);
      if (argc > 3) {
        if (argt[3] != NUMBER) return -1;
        x0 = atof(argv[3]);
      }
    }
  }

  eb = atof(argv[0]);
  SetCEBorn(eb, x, x1, x0);

  return 0;
}

static int PSetCIBorn(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int x;

  if (argc != 1) return -1;
  if (argt[0] != NUMBER) return -1;

  x = atoi(argv[0]);
  SetCIBorn(x);

  return 0;
}

static int PSetCIEGrid(int argc, char *argv[], int argt[],
                       ARRAY *variables) {  int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetCIEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetCIEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetCIEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetCIEGridLimits(int argc, char *argv[], int argt[],
                             ARRAY *variables) {
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;

  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
        type = atoi(argv[2]);
      }
    }
  }

  SetCIEGridLimits(emin, emax, type);

  return 0;
}

static int PSetIEGrid(int argc, char *argv[], int argt[],
                      ARRAY *variables) {  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetIEGrid(cfac, ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        if (ig[i] != NUMBER) return -1;
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetIEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetIEGrid(cfac, ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetCIPWGrid(int argc, char *argv[], int argt[],
                        ARRAY *variables) {  int ns, i;
  int n;
  int *m, *step;
  char *v1[MAXNARGS], *v2[MAXNARGS];
  int t1[MAXNARGS], t2[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] != NUMBER) return -1;
    ns = atoi(argv[0]);
    SetCIPWGrid(cfac, -ns, NULL, NULL);
  } else {
    if (argt[0] != LIST || argt[1] != LIST) return -1;
    ns = DecodeArgs(argv[0], v1, t1, variables);
    if (ns <= 0) return -1;
    if (ns != DecodeArgs(argv[1], v2, t2, variables)) return -1;
    m = malloc(ns*sizeof(int));
    step = malloc(ns*sizeof(int));
    for (i = 0; i < ns; i++) {
      if (t1[i] != NUMBER || t2[i] != NUMBER) return -1;
      m[i] = atoi(v1[i]);
      step[i] = atoi(v2[i]);
      free(v1[i]);
      free(v2[i]);
    }
    SetCIPWGrid(cfac, ns, m, step);
    free(m);
    free(step);
  }

  return 0;
}

static int PSetCIQkMode(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int m;
  double tol;

  m = QK_DEFAULT;
  tol = -1;

  if (argc > 2) return -1;
  if (argc > 0) {
    if (argt[0] == STRING) {
      if (strcasecmp(argv[0], "cb") == 0) m = QK_CB;
      else if (strcasecmp(argv[0], "bed") == 0) m = QK_BED;
      else if (strcasecmp(argv[0],"dw") == 0) m = QK_DW;
      else return -1;
    } else if (argt[0] == NUMBER) {
      m = atoi(argv[0]);
      if (m < QK_CB) return -1;
    }
    if (argc > 1) {
      if (argt[1] != NUMBER) return -1;
      tol = atof(argv[1]);
    }
  }

  SetCIQkMode(m, tol);

  return 0;
}

static int PSetMaxRank(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }
  SetMaxRank(cfac, 2*atoi(argv[0]));

  return 0;
}

static int PSetOptimizeControl(int argc, char *argv[], int argt[],
                               ARRAY *variables) {
  int maxiter, iprint;
  double tol, s;

  iprint = 0;

  if (argc != 3 && argc != 4) return -1;
  tol = atof(argv[0]);
  s = atof(argv[1]);
  maxiter = atoi(argv[2]);
  if (argc == 4) iprint = atoi(argv[3]);

  SetOptimizeControl(cfac, tol, s, maxiter, iprint);

  return 0;
}

static int PSetHydrogenicNL(int argc, char *argv[], int argt[],
                            ARRAY *variables){
  int n, k, nm, km;

  n = -1;
  k = -1;
  nm = -1;
  km = -1;

  if (argc > 0) {
    n = atoi(argv[0]);
    if (argc > 1) {
      k = atoi(argv[1]);
      if (argc > 2) {
        nm = atoi(argv[2]);
        if (argc > 3) {
          km = atoi(argv[3]);
        }
      }
    }
  }
  SetHydrogenicNL(cfac, n, k, nm, km);

  return 0;
}

static int PSetPEGrid(int argc, char *argv[], int argt[],
                      ARRAY *variables) { int n, ng, i, err;
  double xg[MAXNE];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetPEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetPEGridDetail(ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetPEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetPEGridLimits(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  double emin, emax;
  int type;

  emin = -1;
  emax = -1;
  type = 0;

  if (argc > 0) {
    emin = atof(argv[0]);
    if (argc > 1) {
      emax = atof(argv[1]);
      if (argc > 2) {
        type = atoi(argv[2]);
      }
    }
  }

  SetPEGridLimits(emin, emax, type);

  return 0;
}

static int PSetRecPWOptions(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  int kl_interp, max_kl;

  max_kl = -1;
  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    kl_interp = atoi(argv[0]);
  } else if (argc == 2) {
    if (argt[0] != NUMBER) return -1;
    kl_interp = atoi(argv[0]);
    if (argt[1] != NUMBER) return -1;
    max_kl = atoi(argv[1]);
  } else {
    return -1;
  }

  if (max_kl < 0) max_kl = kl_interp;

  SetRecPWOptions(kl_interp, max_kl);

  return 0;
}

static int PSetRecQkMode(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int m;
  double tol;

  m = QK_DEFAULT;
  tol = -1;

  if (argc > 2) return -1;
  if (argc > 0) {
    if (argt[0] == STRING) {
      if (strcasecmp(argv[0], "exact") == 0) m = QK_EXACT;
      else if (strcasecmp(argv[0], "interpolate") == 0) m = QK_INTERPOLATE;
      else if (strcasecmp(argv[0], "fit") == 0) m = QK_FIT;
      else return -1;
    } else if (argt[0] == NUMBER) {
      m = atoi(argv[0]);
      if (m >= QK_CB) return -1;
    }
    if (argc > 1) {
      if (argt[1] != NUMBER) return -1;
      tol = atof(argv[1]);
    }
  }

  SetRecQkMode(m, tol);

  return 0;
}

static int PSetRadialGrid(int argc, char *argv[], int argt[],
                          ARRAY *variables) {
  double rmin = -1.0, ratio = -1.0, asym = -1.0;
  int maxrp;

  switch (argc) {
  case 1:
    break;
  case 4:
    ratio = atof(argv[1]);
    asym = atof(argv[2]);
    rmin = atof(argv[3]);
    break;
  default:
    return -1;
    break;
  }

  maxrp = atoi(argv[0]);

  return SetRadialGrid(cfac, maxrp, ratio, asym, rmin);
}

static int PSetRecPWLimits(int argc, char *argv[], int argt[],
                           ARRAY *variables) {
  int m1, m2;

  if (argc == 2) {
    if (argt[0] != NUMBER) return -1;
    m1 = atoi(argv[0]);
    if (argt[1] != NUMBER) return -1;
    m2 = atoi(argv[1]);
  } else {
    return -1;
  }

  SetRecPWLimits(m1, m2);

  return 0;
}

static int PSetRecSpectator(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  int n_spec, n_frozen;

  n_spec = 0;
  n_frozen = 0;

  if (argc < 1 || argc > 2) return -1;
  n_spec = atoi(argv[0]);
  if (argc > 1) {
    n_frozen = atoi(argv[1]);
  }

  if (n_frozen == 0) n_frozen = n_spec;

  SetRecSpectator(n_frozen, n_spec);

  return 0;
}

static int PSetRRTEGrid(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int i, n;
  double xg[MAXNTE];
  int ng, err;
  double emin, emax;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  n = argc;
  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetRRTEGrid(cfac, ng, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        if (ig[i] != NUMBER) return -1;
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetRRTEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    err = SetRRTEGrid(cfac, ng, emin, emax);
  } else {
    return -1;
  }

  if (err < 0) return -1;

  return 0;
}

static int PSetSE(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  int c;

  if (argc != 1) return -1;
  c = atoi(argv[0]);

  SetSE(cfac, c);

  return 0;
}

static int PSetVP(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  int c;

  if (argc != 1) return -1;
  c = atoi(argv[0]);

  SetVP(cfac, c);

  return 0;
}

static int PSetBreit(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  int c;

  if (argc != 1) return -1;
  c = atoi(argv[0]);

  SetBreit(cfac, c);

  return 0;
}

static int PSetMS(int argc, char *argv[], int argt[],
                  ARRAY *variables) {
  int c1, c2;

  if (argc != 2) return -1;
  c1 = atoi(argv[0]);
  c2 = atoi(argv[1]);

  SetMS(cfac, c1, c2);

  return 0;
}

static int PSetScreening(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int n_screen;
  int *screened_n = NULL;
  double screened_charge;
  int i, kl;
  char *v[MAXNARGS];
  int t[MAXNARGS];

  n_screen = 0;
  screened_charge = 1.0;
  kl = 1;

  if (argc < 1 || argc > 3) return -1;
  if (argt[0] != LIST && argt[0] != TUPLE) return -1;
  n_screen = DecodeArgs(argv[0], v, t, variables);
  if (argc > 1) {
    screened_charge = atof(argv[1]);
    if (argc > 2) {
      kl = atoi(argv[2]);
    }
  }

  screened_n = malloc(sizeof(int)*n_screen);
  for (i = 0; i < n_screen; i++) {
    if (t[i] != NUMBER) return -1;
    screened_n[i] = atoi(v[i]);
    free(v[i]);
  }

  SetScreening(cfac, n_screen, screened_n, screened_charge, kl);

  return 0;
}

static int PSetTransitionCut(int argc, char *argv[], int argt[],
                             ARRAY *variables) {
  cfac_errmsg(cfac, "SetTransitionCut() is defunct\n");

  return 0;
}

static int PSetTransitionOptions(int argc, char *argv[], int argt[],
                                 ARRAY *variables) {
  int gauge, mode;

  if (argc != 2) return -1;

  gauge = atoi(argv[0]);
  mode = atoi(argv[1]);

  SetTransitionOptions(cfac, gauge, mode);

  return 0;
}

static int PSetUsrCEGrid(int argc, char *argv[], int argt[],
                         ARRAY *variables) {  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrCEEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetUsrCEEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCEEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetUsrCEGridType(int argc, char *argv[], int argt[],
                             ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrCEEGridType(atoi(argv[0]));
  return 0;
}

static int PSetUsrCIEGrid(int argc, char *argv[], int argt[],
                          ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrCIEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetUsrCIEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrCIEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetUsrCIEGridType(int argc, char *argv[], int argt[],
                              ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrCIEGridType(atoi(argv[0]));
  return 0;
}

static int PSetUsrPEGrid(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int n, ng, i, err;
  double xg[MAXNUSR];
  double emin, emax, eth;
  char *vg[MAXNARGS];
  int ig[MAXNARGS];

  eth = 0;
  n = argc;

  if (n == 1) {
    if (argt[0] == NUMBER) {
      ng = atoi(argv[0]);
      err = SetUsrPEGrid(cfac, ng, -1.0, -1.0, 0.0);
    } else if (argt[0] == LIST || argt[0] == TUPLE) {
      ng = DecodeArgs(argv[0], vg, ig, variables);
      for (i = 0; i < ng; i++) {
        xg[i] = atof(vg[i]);
        xg[i] /= HARTREE_EV;
        free(vg[i]);
      }
      err = SetUsrPEGridDetail(cfac, ng, xg);
    } else {
      return -1;
    }
  } else if (n == 3 || n == 4) {
    ng = atoi(argv[0]);
    emin = atof(argv[1]);
    emax = atof(argv[2]);
    if (n == 4) eth = atof(argv[3]);
    emin /= HARTREE_EV;
    emax /= HARTREE_EV;
    if (eth > 0) eth /= HARTREE_EV;
    err = SetUsrPEGrid(cfac, ng, emin, emax, eth);
  } else {
    return -1;
  }

  if (err < 0) return -1;
  return 0;
}

static int PSetUsrPEGridType(int argc, char *argv[], int argt[],
                             ARRAY *variables) {
  if (argc != 1 || argt[0] != NUMBER) return -1;
  SetUsrPEGridType(atoi(argv[0]));
  return 0;
}

static int PSolveBound(int argc, char *argv[], int argt[],
                       ARRAY *variables) {
  int n, kappa;
  ORBITAL *orb;
  int k;

  if (argc != 2 || argt[0] != NUMBER || argt[1] != NUMBER) return -1;
  n = atoi(argv[0]);
  kappa = atoi(argv[1]);

  if (n <= 0) {
    cfac_errmsg(cfac, "n must be greater than 0 for SolveBound\n");
    return -1;
  }

  k = OrbitalIndex(cfac, n, kappa, 0.0);
  if (k < 0) {
    cfac_errmsg(cfac, "Fatal error in solving Dirac equation\n");
    return -1;
  }

  orb = GetOrbital(cfac, k);
  printf("Energy = %16.8E\n", orb->energy);

  return 0;
}

static int PCutMixing(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int nlev, n, *ilev, *kg;
  double c;

  nlev = 0;
  n = 0;
  if (argc < 2 || argc > 3) return -1;
  nlev = SelectLevels(&ilev, argv[0], argt[0], variables);
  if (nlev <= 0) goto DONE;
  n = DecodeGroupArgs(&kg, 1, &(argv[1]), &(argt[1]), variables);
  if (n <= 0) goto DONE;
  if (argc == 3) c = atof(argv[2]);
  else c = 0.0;

  CutMixing(cfac, nlev, ilev, n, kg, c);

 DONE:
  if (nlev > 0) free(ilev);
  if (n > 0) free(kg);

  return 0;
}

static int PSetSymmetry(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  if (argc == 2 && argt[0] == NUMBER) {
    int nj, ip, *kg;
    ip = atoi(argv[0]);
    nj = IntFromList(argv[1], argt[1], variables, &kg);
    SetSymmetry(cfac, ip, nj, kg);
    free(kg);
    return 0;
  } else {
    return -1;
  }
}

static int PStructure(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int ng, ngp;
  int ip, nlevels;
  int *kg, *kgp;

  ng = 0;
  ngp = 0;
  kgp = NULL;
  ip = 0;

  if (argc < 1 || argc > 4) return -1;

  if (argc == 1) {
    if (argt[0] != STRING) return -1;
    ng = DecodeGroupArgs(&kg, 0, NULL, NULL, variables);
    if (ng < 0) return -1;
  } else
  if (argc == 2 && argt[0] == STRING && argt[1] == NUMBER) {
    int nele = atoi(argv[1]);
    ng = SelectNeleGroups(nele, &kg);
  } else {
    if (argc == 4) ip = atoi(argv[3]);
    if (argt[0] != STRING) return -1;
    if (argt[1] != LIST && argt[1] != TUPLE) return -1;
    ng = DecodeGroupArgs(&kg, 1, &(argv[1]), &(argt[1]), variables);
    if (ng < 0) return -1;
    if (argc >= 3) {
      if (argt[2] != LIST && argt[2] != TUPLE) return -1;
      ngp = DecodeGroupArgs(&kgp, 1, &(argv[2]), &(argt[2]), variables);
    }
  }

  if (ngp < 0) return -1;

  nlevels = cfac_get_num_levels(cfac);

  if (cfac_calculate_structure(cfac, ng, kg, ngp, kgp, ip) != 0) {
    return -1;
  }

  SaveLevels(cfac, argv[0], nlevels, -1);

  if (ng > 0) free(kg);
  if (ngp > 0) free(kgp);

  return 0;
}

static int PSetUTA(int argc, char *argv[], int argt[],
                   ARRAY *variables) {
  int m;

  if (argc == 1) {
    if (argt[0] != NUMBER) return -1;
    m = atoi(argv[0]);
  } else {
    return -1;
  }

  cfac_set_uta(cfac, m);

  return 0;
}

static int PPrepAngular(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;

  if (argc < 1 || argc > 2) return -1;
  nlow = SelectLevels(&low, argv[0], argt[0], variables);
  if (nlow <= 0) return -1;
  if (argc == 2) {
    nup = SelectLevels(&up, argv[1], argt[1], variables);
    if (nup <= 0) {
      free(low);
      return -1;
    }
  }
  PrepAngular(cfac, nlow, low, nup, up);

  return 0;
}

static int PTransitionTable(int argc, char *argv[], int argt[],
                            ARRAY *variables) {
  int m;
  int nlow, nup, *low, *up;

  nlow = 0;
  nup = 0;
  low = NULL;
  up = NULL;
  m = -1;

  if (argc < 3 || argc > 4) return -1;

  if (argt[0] != STRING) {
    return -1;
  }

  if (argc == 4) {
    if (argt[3] != NUMBER) {
      return -1;
    } else {
      m = atoi(argv[3]);
    }
  }

  if (argt[1] == NUMBER && argt[2] == NUMBER) {
    int nele = atoi(argv[1]);
    m = atoi(argv[2]);
    nlow = SelectNeleLevels(cfac, nele, &low);
    nup = nlow;
    up = low;
  } else {
    nlow = SelectLevels(&low, argv[1], argt[1], variables);
    nup  = SelectLevels(&up, argv[2], argt[2], variables);
  }
  if (nlow <= 0 || nup <= 0) {
    return -1;
  }

  if (m == 0) {
    cfac_errmsg(cfac, "m cannot be zero\n");
    return -1;
  }

  SaveTransition(cfac, nlow, (unsigned int*) low, nup, (unsigned int*) up,
    argv[0], m);

  if (low != up) {
    free(low);
  }
  if (up) {
    free(up);
  }

  return 0;
}

static int PWaveFuncTable(int argc, char *argv[], int argt[],
                          ARRAY *variables) {
  int k, n;
  double e;

  if (argc != 3 && argc != 4) return -1;

  n = atoi(argv[1]);
  k = atoi(argv[2]);
  if (argc == 4) {
    e = atof(argv[3]);
  } else {
    e = 0.0;
  }

  WaveFuncTable(cfac, argv[0], n, k, e);

  return 0;
}

static int PSetOptimizeMaxIter(int argc, char *argv[], int argt[],
                               ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetOptimizeMaxIter(cfac, m);
  return 0;
}

static int PSetOptimizeStabilizer(int argc, char *argv[], int argt[],
                                  ARRAY *variables) {
  double m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetOptimizeStabilizer(cfac, m);
  return 0;
}

static int PSetOptimizePrint(int argc, char *argv[], int argt[],
                             ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetOptimizePrint(cfac, m);
  return 0;
}

static int PSetOptimizeTolerance(int argc, char *argv[], int argt[],
                                 ARRAY *variables) {
  double m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atof(argv[0]);
  SetOptimizeTolerance(cfac, m);
  return 0;
}

static int PSetCELQR(int argc, char *argv[], int argt[],
                     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELQR(m);
  return 0;
}

static int PSetCELMax(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELMax(m);
  return 0;
}

static int PSetCELCB(int argc, char *argv[], int argt[],
                     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCELCB(m);
  return 0;
}

static int PSetCILQR(int argc, char *argv[], int argt[],
                     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILQR(m);
  return 0;
}

static int PSetCILMax(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILMax(m);
  return 0;
}

static int PSetCILMaxEject(int argc, char *argv[], int argt[],
                           ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILMaxEject(m);
  return 0;
}

static int PSetCILCB(int argc, char *argv[], int argt[],
                     ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetCILCB(m);
  return 0;
}

static int PSetTransitionMode(int argc, char *argv[], int argt[],
                              ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionMode(cfac, m);
  return 0;
}

static int PSetTransitionGauge(int argc, char *argv[], int argt[],
                               ARRAY *variables) {
  int m;

  if (argc != 1 || argt[0] != NUMBER) {
    return -1;
  }

  m = atoi(argv[0]);
  SetTransitionGauge(cfac, m);
  return 0;
}

static int PSetSlaterCut(int argc, char *argv[], int argt[],
                         ARRAY *variables) {
  int k0, k1;

  if (argc != 2) return -1;

  k0 = atoi(argv[0]);
  k1 = atoi(argv[1]);

  SetSlaterCut(cfac, k0, k1);

  return 0;
}


static int PAppendTable(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  if (argc != 1) return -1;
  AppendTable(cfac, argv[0]);

  return 0;
}

static int PJoinTable(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  if (argc != 3) return -1;

  JoinTable(cfac, argv[0], argv[1], argv[2]);

  return 0;
}

static int PSetFields(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int m;
  double b, e, a;

  if (argc < 3 || argc > 4) return -1;
  m = 0;
  b = atof(argv[0]);
  e = atof(argv[1]);
  a = atof(argv[2]);
  if (argc > 3) {
    m = atoi(argv[3]);
  }

  SetFields(cfac, b, e, a, m);

  return 0;
}

static int PStructureEB(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int n, *ilev, res;

  if (argc != 2 || argt[0] != STRING || argt[1] != LIST) return -1;

  n = SelectLevels(&ilev, argv[1], argt[1], variables);
  if (n <= 0) return -1;

  res = StructureEB(cfac, argv[0], n, ilev);
  free(ilev);

  if (res != 0) {
    return -1;
  } else {
    return 0;
  }
}

static int PTransitionTableEB(int argc, char *argv[], int argt[],
                              ARRAY *variables) {
  int m, nlow, *low, nup, *up;

  if (argc < 3 || argc > 4) return -1;

  if (argc == 4) m = atoi(argv[3]);
  else m = -1;

  if (m == 0) {
    cfac_errmsg(cfac, "m cannot be zero\n");
    return -1;
  }

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) {
    cfac_errmsg(cfac, "cannot determine levels in lower\n");
    return -1;
  }
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) {
    cfac_errmsg(cfac, "cannot determine levels in upper\n");
    return -1;
  }

  SaveTransitionEB(cfac, nlow, low, nup, up, argv[0], m);
  free(low);
  free(up);

  return 0;
}

static int PCETableEB(int argc, char *argv[], int argt[],
                      ARRAY *variables) {
  int nlow, nup, *low, *up, m;

  if (argc < 3 || argc > 4) return -1;
  if (argc == 4) m = atoi(argv[3]);
  else m = 0;

  nlow = SelectLevels(&low, argv[1], argt[1], variables);
  if (nlow <= 0) return -1;
  nup = SelectLevels(&up, argv[2], argt[2], variables);
  if (nup <= 0) return -1;

  if (m == 0) {
    SaveExcitationEB(cfac, nlow, low, nup, up, argv[0]);
  } else {
    SaveExcitationEBD(cfac, nlow, low, nup, up, argv[0]);
  }

  free(low);
  free(up);

  return 0;
}

static int PSlaterCoeff(int argc, char *argv[], int argt[],
                        ARRAY *variables) {
  int nlev, *ilev, na, nb, i, *n, *kappa;
  double *nq;
  SHELL *sa, *sb;

  if (argc != 4) return -1;
  if (argt[1] != LIST) return -1;
  if (argt[2] != STRING) return -1;
  if (argt[3] != STRING) return -1;

  nlev = SelectLevels(&ilev, argv[1], argt[1], variables);
  na = GetAverageConfigFromString(cfac, &n, &kappa, &nq, argv[2]);
  sa = malloc(sizeof(SHELL)*na);
  for (i = 0; i < na; i++) {
    sa[i].n = n[i];
    sa[i].kappa = kappa[i];
  }
  if (na > 0) {
    free(n);
    free(kappa);
    free(nq);
  }
  nb = GetAverageConfigFromString(cfac, &n, &kappa, &nq, argv[3]);
  sb = malloc(sizeof(SHELL)*nb);
  for (i = 0; i < nb; i++) {
    sb[i].n = n[i];
    sb[i].kappa = kappa[i];
  }
  if (nb > 0) {
    free(n);
    free(kappa);
    free(nq);
  }


  if (nlev > 0 && na > 0 && nb > 0) {
    SlaterCoeff(cfac, argv[0], nlev, ilev, na, sa, nb, sb);
    free(ilev);
    free(sa);
    free(sb);
  }

  return 0;
}

static int PGeneralizedMoment(int argc, char *argv[], int argt[],
                              ARRAY *variables) {
  int n0, k0, n1, k1, m;
  double e1;

  m = atoi(argv[1]);
  n0 = atoi(argv[2]);
  k0 = atoi(argv[3]);
  n1 = atoi(argv[4]);
  k1 = atoi(argv[5]);
  if (argc == 7) e1 = atof(argv[6]);
  else e1 = 0.0;
  PrintGeneralizedMoments(cfac, argv[0], m, n0, k0, n1, k1, e1);

  return 0;
}

static METHOD methods[] = {
  {"AITable", PAITable},
  {"AITableMSub", PAITableMSub},
  {"AppendTable", PAppendTable},
  {"AvgConfig", PAvgConfig},
  {"BasisTable", PBasisTable},
  {"CETable", PCETable},
  {"CETableEB", PCETableEB},
  {"CETableMSub", PCETableMSub},
  {"CITable", PCITable},
  {"CITableMSub", PCITableMSub},
  {"CheckEndian", PCheckEndian},
  {"Closed", PClosed},
  {"Config", PConfig},
  {"ConfigEnergy", PConfigEnergy},
  {"CorrectEnergy", PCorrectEnergy},
  {"CutMixing", PCutMixing},
  {"Exit", PExit},
  {"GeneralizedMoment", PGeneralizedMoment},
  {"GetConfigNR", PGetConfigNR},
  {"GetPotential", PGetPotential},
  {"Info", PInfo},
  {"JoinTable", PJoinTable},
  {"ListConfig", PListConfig},
  {"MemENTable", PMemENTable},
  {"OptimizeRadial", POptimizeRadial},
  {"Pause", PPause},
  {"PrepAngular", PPrepAngular},
  {"Print", PPrint},
  {"PrintTable", PPrintTable},
  {"RRMultipole", PRRMultipole},
  {"RRTable", PRRTable},
  {"RecStates", PRecStates},
  {"RefineRadial", PRefineRadial},
  {"SetAICut", PSetAICut},
  {"SetAngZCut", PSetAngZCut},
  {"SetAngZOptions", PSetAngZOptions},
  {"SetAngleGrid", PSetAngleGrid},
  {"SetAtom", PSetAtom},
  {"SetAvgConfig", PSetAvgConfig},
  {"SetBreit", PSetBreit},
  {"SetCEBorn", PSetCEBorn},
  {"SetCEGrid", PSetCEGrid},
  {"SetCEGridLimits", PSetCEGridLimits},
  {"SetCEGridType", PSetCEGridType},
  {"SetCELCB", PSetCELCB},
  {"SetCELMax", PSetCELMax},
  {"SetCELQR", PSetCELQR},
  {"SetCEPWGrid", PSetCEPWGrid},
  {"SetCEPWGridType", PSetCEPWGridType},
  {"SetCEPWOptions", PSetCEPWOptions},
  {"SetCIBorn", PSetCIBorn},
  {"SetCIEGrid", PSetCIEGrid},
  {"SetCIEGridLimits", PSetCIEGridLimits},
  {"SetCILCB", PSetCILCB},
  {"SetCILMax", PSetCILMax},
  {"SetCILMaxEject", PSetCILMaxEject},
  {"SetCILQR", PSetCILQR},
  {"SetCILevel", PSetCILevel},
  {"SetCIPWGrid", PSetCIPWGrid},
  {"SetCIQkMode", PSetCIQkMode},
  {"SetErrorOutput", PSetErrorOutput},
  {"SetFields", PSetFields},
  {"SetHydrogenicNL", PSetHydrogenicNL},
  {"SetIEGrid", PSetIEGrid},
  {"SetMS", PSetMS},
  {"SetMaxRank", PSetMaxRank},
  {"SetMixCut", PSetMixCut},
  {"SetOptimizeControl", PSetOptimizeControl},
  {"SetOptimizeMaxIter", PSetOptimizeMaxIter},
  {"SetOptimizePrint", PSetOptimizePrint},
  {"SetOptimizeStabilizer", PSetOptimizeStabilizer},
  {"SetOptimizeTolerance", PSetOptimizeTolerance},
  {"SetPEGrid", PSetPEGrid},
  {"SetPEGridLimits", PSetPEGridLimits},
  {"SetRRTEGrid", PSetRRTEGrid},
  {"SetRadialGrid", PSetRadialGrid},
  {"SetRecPWLimits", PSetRecPWLimits},
  {"SetRecPWOptions", PSetRecPWOptions},
  {"SetRecQkMode", PSetRecQkMode},
  {"SetRecSpectator", PSetRecSpectator},
  {"SetSE", PSetSE},
  {"SetScreening", PSetScreening},
  {"SetSlaterCut", PSetSlaterCut},
  {"SetSymmetry", PSetSymmetry},
  {"SetTEGrid", PSetTEGrid},
  {"SetTransitionCut", PSetTransitionCut},
  {"SetTransitionGauge", PSetTransitionGauge},
  {"SetTransitionMode", PSetTransitionMode},
  {"SetTransitionOptions", PSetTransitionOptions},
  {"SetUTA", PSetUTA},
  {"SetUsrCEGrid", PSetUsrCEGrid},
  {"SetUsrCEGridType", PSetUsrCEGridType},
  {"SetUsrCIEGrid", PSetUsrCIEGrid},
  {"SetUsrCIEGridType", PSetUsrCIEGridType},
  {"SetUsrPEGrid", PSetUsrPEGrid},
  {"SetUsrPEGridType", PSetUsrPEGridType},
  {"SetVP", PSetVP},
  {"SlaterCoeff", PSlaterCoeff},
  {"SolveBound", PSolveBound},
  {"StoreClose", PStoreClose},
  {"StoreInit", PStoreInit},
  {"StoreTable", PStoreTable},
  {"Structure", PStructure},
  {"StructureEB", PStructureEB},
  {"TRTable", PTransitionTable},
  {"TRTableEB", PTransitionTableEB},
  {"TransitionTable", PTransitionTable},
  {"WaveFuncTable", PWaveFuncTable},
  {"", NULL}
};

static void usage(FILE *fp, const char *progname) {
    fprintf(fp, "Usage:\n");
    fprintf(fp, "       %s            run in the interactive mode \n",
        progname);
    fprintf(fp, "       %s [FILE]...  execute commands in one or more files\n",
        progname);
    fprintf(fp, "       %s -V         print version info and exit\n",
        progname);
    fprintf(fp, "       %s -h         display this help and exit\n",
        progname);
}

static int InitFac() {
  gsl_ieee_env_setup();

  cfac = cfac_new();
  if (!cfac) {
    cfac_errmsg(cfac, "Initialization failed\n");
    return -1;
  }

  if (InitDBase() != 0) {
    return -1;
  }
  if (InitExcitation(cfac) != 0) {
    return -1;
  }
  if (InitRecombination() != 0) {
    return -1;
  }
  if (InitIonization(cfac) != 0) {
    return -1;
  }

  return 0;
}

int main(int argc, const char *argv[]) {
  int i;
  FILE *f;
  int cmdlen = 0;

/* fix non-standard number of exponent digits in the MSVC runtime */
#ifdef _WIN32
  _set_output_format(_TWO_DIGIT_EXPONENT);
#endif

  if (InitFac() < 0) {
    cfac_errmsg(cfac, "Initialization failed\n");
    exit(1);
  }

  for (i = 0; i < argc; i++) {
    cmdlen += strlen(argv[i]) + 1;
  }

  cmdline = malloc(cmdlen + 1);
  cmdline[0] = '\0';
  strcat(cmdline, argv[0]);
  for (i = 1; i < argc; i++) {
    strcat(cmdline, " ");
    strcat(cmdline, argv[i]);
  }

  if (argc == 1) {
    EvalFile(stdin, 1, methods, stderr);
  } else {
    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "-h")) {
        usage(stdout, argv[0]);
        exit(0);
      } else
      if (!strcmp(argv[i], "-V")) {
        cfac_verinfo();
        exit(0);
      } else {
        f = fopen(argv[i], "r");
        if (!f) {
          cfac_errmsg(cfac, "Cannot open file '%s'\n", argv[i]);
          exit(1);
        }
        EvalFile(f, 0, methods, stderr);
      }
    }
  }

  if (cfac->ncorrections > 0) {
    cfac_errmsg(cfac, "Warning: %d energy correction(s) have not been applied\n",
        cfac->ncorrections);
  }
  cfac_free(cfac);

  exit(0);
}
