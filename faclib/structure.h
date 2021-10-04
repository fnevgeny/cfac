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

#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "recouple.h"
#include "config.h"

typedef struct _HAMILTON_ {
  int     pj;        /* symmetry index                                   */

  int     dim;       /* dimension of the subset of basis of interest,
                        dim <= n_basis                                   */
  int     n_basis;   /* dimension of the basis                           */
  int     hsize;     /* size of Hamiltonian (.hamilton)                  */
  int     msize;     /* size of the mixing coefficients (.mixing),
                        n_basis*(dim+1)                                  */

  int    *basis;     /* basis                                            */
  double *hamilton;  /* matrix elements of H,
                        H1[dim*dim] &
                        H2[n_basis-dim] &
                        B[dim*(n_basis-dim)], Eq. (29) in structure.pdf  */
  double *mixing;    /* mixing coefficients                              */
} HAMILTON;

typedef struct _SHAMILTON_ {
  int pj;            /* symmetry index                                   */
  int nbasis, nlevs; /* equal to n_basis & dim of HAMILTON, respectively */
  STATE **basis;     /* n_basis pointers to STATE structures             */
  unsigned char closed[MBCLOSE];
} SHAMILTON;

typedef struct _LEVEL_ {
  int uta;             /* flag indicating it is an UTA level or not         */
  int ilev;            /* level # in corresponding pj Hamiltonian (non-UTA) */
  int pj;              /* parity & j encoded (only parity in UTA)           */
  int iham;            /* symmetry Hamiltonian # (non-UTA)                  */
  int n_basis;
  int pb;              /* principle (??) basis (non-UTA)                    */
  int kpb[NPRINCIPLE];
  int ibase;           /* index of the base level, if it can be determined  */
  int *basis;
  short *ibasis;
  double *mixing;
  double energy;       /* energy                                            */

  int uta_g;           /* used to be ilev + 1 in UTA                        */
  int uta_p;           /* used to be pj in UTA                              */
  int uta_cfg_g;       /* used to be iham in UTA                            */
  int uta_g_cfg;       /* used to be pb in UTA                              */

  char *name;          /* level name, unique within the same charge state   */
} LEVEL;

typedef struct _LEVEL_ION_ {
  int imin;
  int imax;
} LEVEL_ION;

typedef struct _ANGZ_DATUM_ {
  int ns;
  int *nz;
  void **angz;
  double **mk;
} ANGZ_DATUM;

typedef struct _ANGULAR_ZMIX_ {
  double coeff;
  short k;
  short k0;
  short k1;
} ANGULAR_ZMIX;

typedef struct _ANGULAR_ZFB_ {
  double coeff;
  short kb;
} ANGULAR_ZFB;

typedef struct _ANGULAR_ZxZMIX_ {
  double coeff;
  short k;
  short k0;
  short k1;
  short k2;
  short k3;
} ANGULAR_ZxZMIX;

typedef struct _ANGULAR_FROZEN_ {
  int nts, ncs;
  int *ts, *cs;
  int *nz, *nzfb, *nzxzfb;
  ANGULAR_ZMIX **z;
  ANGULAR_ZFB **zfb;
  ANGULAR_ZxZMIX **zxzfb;
} ANGULAR_FROZEN;

typedef struct _ECORRECTION_ {
  int nele;      /* number of electrons                                     */
  char *name;    /* level name                                              */
  double e;      /* corrected energy                                        */
  char *refname; /* name of the reference level                             */
  LEVEL *ref;    /* reference level, when resolved                          */
  int applied;   /* boolean flag indicating the correction has been applied */
} ECORRECTION;

HAMILTON *ConstructHamilton(cfac_t *cfac,
    int isym, int k, const int *kg, int kp, const int *kgp);
int ValidBasis(cfac_t *cfac, STATE *s, int k, int *kg, int n);
HAMILTON *ConstructHamiltonFrozen(cfac_t *cfac,
    int isym, int k, int *kg, int n, int nc, int *kc);
double HamiltonElement(cfac_t *cfac, int isym, int isi, int isj);
double HamiltonElementFrozen(cfac_t *cfac, int isym, int isi, int isj);
double HamiltonElementFB(cfac_t *cfac, int isym, int isi, int isj);
double Hamilton2E(cfac_t *cfac, int n_shells, SHELL_STATE *sbra,
                  SHELL_STATE *sket,INTERACT_SHELL *s);
double Hamilton1E(cfac_t *cfac, int n_shells, SHELL_STATE *sbra,
                  SHELL_STATE *sket,INTERACT_SHELL *s);
int DiagonalizeHamilton(const cfac_t *cfac, HAMILTON *h);
int AddToLevels(cfac_t *cfac, HAMILTON *h, int ng, const int *kg);

int AddECorrection(cfac_t *cfac, int nele, const char *name,
    const char *refname, double e);
void InitECorrectionData(void *p, int n);
void FreeECorrectionData(void *p);

LEVEL *GetLevel(const cfac_t *cfac, int k);
LEVEL *GetEBLevel(const cfac_t *cfac, int k);
int LevelTotalJ(cfac_t *cfac, int k);
int GetNumEBLevels(const cfac_t *cfac);
int GetLevNumElectrons(const cfac_t *cfac, const LEVEL *lev);
int GetNumElectrons(const cfac_t *cfac, int k);
int SortMixing(int start, int n, LEVEL *lev, SYMMETRY *sym);
int GetPrincipleBasis(double *mix, int d, int *kpb);
int SortLevels(cfac_t *cfac, int start, int n, int EB);
void SetSymmetry(cfac_t *cfac, int p, int n, int *j);
int GetBaseJ(cfac_t *cfac, STATE *s);
void AngularFrozen(cfac_t *cfac, int nts, int *ts, int ncs, int *cs);
void ClearAngularFrozen(cfac_t *cfac);
int PrepAngular(cfac_t *cfac, int n1, int *is1, int n2, int *is2);
int AngularZMix(cfac_t *cfac,
    ANGULAR_ZMIX **ang, int lower, int upper, int mink, int maxk);
int CompareAngularZMix(const void *c1, const void *c2);
int CompareAngularZxZMix(const void *c1, const void *c2);
int CompareAngularZFB(const void *c1, const void *c2);
int PackAngularZxZMix(int *n, ANGULAR_ZxZMIX **ang, int nz);
int PackAngularZMix(int *n, ANGULAR_ZMIX **ang, int nz);
int PackAngularZFB(int *n, ANGULAR_ZFB **ang, int nz);
int AngularZFreeBound(cfac_t *cfac, ANGULAR_ZFB **ang, int lower, int upper);
int AngularZMixStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2);
int AngZSwapBraKet(cfac_t *cfac, int nz, ANGULAR_ZMIX *ang, int p);
int AngularZFreeBoundStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2);
int AngularZxZFreeBoundStates(cfac_t *cfac, ANGZ_DATUM **ad, int ih1, int ih2);
int AddToAngularZxZ(cfac_t *cfac, int *n, int *nz, ANGULAR_ZxZMIX **ang,
                    int n_shells, int phase, SHELL_STATE *sbra,
                    SHELL_STATE *sket, INTERACT_SHELL *s, int m);
int AddToAngularZxZMix(const cfac_t *cfac,
    int *n, int *nz, ANGULAR_ZxZMIX **ang,
    int k, int k0, int k1, int k2, int k3, double coeff);
int AddToAngularZMix(const cfac_t *cfac, int *n, int *nz, ANGULAR_ZMIX **ang,
                     int k, int k0, int k1, double coeff);
int AddToAngularZFB(const cfac_t *cfac, int *n, int *nz, ANGULAR_ZFB **ang,
                    int kb, double coeff);
int AngularZxZFreeBound(cfac_t *cfac, ANGULAR_ZxZMIX **ang, int lower, int upper);
int GetBasisTable(cfac_t *cfac, char *fn, int m);
int ConstructLevelName(const cfac_t *cfac, const STATE *basis,
                       char *name, char *sname, char *nc, int *vnl);

int SaveEBLevels(cfac_t *cfac, char *fn, int m, int n);
int SetAngZOptions(cfac_t *cfac, int n, double mc, double c);
int SetAngZCut(cfac_t *cfac, double c);
int SetCILevel(cfac_t *cfac, int m);
int SetMixCut(cfac_t *cfac, double c, double c2);
int TestHamilton(cfac_t *cfac);
void CutMixing(cfac_t *cfac, int nlev, int *ilev, int n, int *kg, double c);
void SetFields(cfac_t *cfac, double b, double e, double a, int no_diamag);
void GetFields(const cfac_t *cfac, double *b, double *e, double *a);
int CodeBasisEB(int s, int m);
void DecodeBasisEB(int k, int *s, int *m);
HAMILTON *ConstructHamiltonEB(cfac_t *cfac, int n, int *ilev);
int StructureEB(cfac_t *cfac, char *fn, int n, int *ilev);
double HamiltonElementEB(const cfac_t *cfac, HAMILTON *h, int i0, int j0);

int SlaterCoeff(cfac_t *cfac,
    char *fn, int nlevs, int *ilevs, int na, SHELL *sa,
                int nb, SHELL *sb);
void AddSlaterCoeff(const cfac_t *cfac, double *c, double a, int n_shells,
                    SHELL_STATE *sbra, SHELL_STATE *sket,
                    INTERACT_SHELL *s, int na, SHELL *sa,
                    int nb, SHELL *sb);

#endif
