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

  int     dim0;      /* factually allocated dim, dim0 >= dim             */
  int     n_basis0;  /* same for n_basis                                 */
  int     hsize0;    /* same for hsize                                   */
  int     msize0;    /* same for msize                                   */

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
  int pj;
  int iham, ilev;
  int n_basis;
  int pb;
  int kpb[NPRINCIPLE];
  int ibase;
  int *basis;
  short *ibasis;
  double *mixing;
  double energy;
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
  int iref;
  int ilev;
  double e;
  int nmin;
  STATE *s;
} ECORRECTION;

typedef struct _TRANSITION_ {
    int nup;
    int nlo;
    LEVEL *lup;
    LEVEL *llo;
    double e;
} TRANSITION;

int ConstructHamilton(cfac_t *cfac,
    int isym, int k0, int k, int *kg, int kp, int *kgp, int md);
int ValidBasis(cfac_t *cfac, STATE *s, int k, int *kg, int n);
int ConstructHamiltonFrozen(cfac_t *cfac,
    int isym, int k, int *kg, int n, int nc, int *kc);
void HamiltonElement1E2E(cfac_t *cfac,
    int isym, int isi, int isj, double *r1, double *r2);
double HamiltonElement(cfac_t *cfac, int isym, int isi, int isj);
double HamiltonElementFrozen(cfac_t *cfac, int isym, int isi, int isj);
double HamiltonElementFB(cfac_t *cfac, int isym, int isi, int isj);
double Hamilton2E(cfac_t *cfac, int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);
double Hamilton1E(cfac_t *cfac, int n_shells, SHELL_STATE *sbra, 
		  SHELL_STATE *sket,INTERACT_SHELL *s);
int DiagonalizeHamilton(cfac_t *cfac);
int AddToLevels(cfac_t *cfac, int ng, int *kg);
int AddECorrection(cfac_t *cfac, int kref, int k, double e, int nmin);
LEVEL *GetLevel(const cfac_t *cfac, int k);
LEVEL *GetEBLevel(const cfac_t *cfac, int k);
int LevelTotalJ(cfac_t *cfac, int k);
int GetNumEBLevels(const cfac_t *cfac);
int GetNumLevels(const cfac_t *cfac);
int GetLevNumElectrons(cfac_t *cfac, const LEVEL *lev);
int GetNumElectrons(cfac_t *cfac, int k);
int SortMixing(int start, int n, LEVEL *lev, SYMMETRY *sym);
int GetPrincipleBasis(double *mix, int d, int *kpb);
int CompareLevels(cfac_t *cfac, LEVEL *lev1, LEVEL *lev2);
int SortLevels(cfac_t *cfac, int start, int n, int m);
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
int AddToAngularZxZMix(int *n, int *nz, ANGULAR_ZxZMIX **ang, 
		       int k, int k0, int k1, int k2, int k3, double coeff);
int AddToAngularZMix(int *n, int *nz, ANGULAR_ZMIX **ang,
		     int k, int k0, int k1, double coeff);
int AddToAngularZFB(int *n, int *nz, ANGULAR_ZFB **ang,
		    int kb, double coeff);
int AngularZxZFreeBound(cfac_t *cfac, ANGULAR_ZxZMIX **ang, int lower, int upper);
int GetBasisTable(cfac_t *cfac, char *fn, int m);
int ConstructLevelName(cfac_t *cfac, char *name, char *sname, char *nc, 
		       int *vnl, STATE *basis);
int GetTransition(const cfac_t *cfac,
    int nlo, int nup, TRANSITION *tr, int *swapped);
int SaveLevels(cfac_t *cfac, char *fn, int m, int n);
int SaveEBLevels(cfac_t *cfac, char *fn, int m, int n);
int SetAngZOptions(cfac_t *cfac, int n, double mc, double c);
int SetAngZCut(cfac_t *cfac, double c);
int SetCILevel(cfac_t *cfac, int m);
int SetMixCut(cfac_t *cfac, double c, double c2);
int TestHamilton(cfac_t *cfac);
void CutMixing(cfac_t *cfac, int nlev, int *ilev, int n, int *kg, double c);
int AllocHamMem(HAMILTON *h, int hdim, int nbasis);
void SetFields(cfac_t *cfac, double b, double e, double a, int m);
void GetFields(const cfac_t *cfac, double *b, double *e, double *a);
int CodeBasisEB(int s, int m);
void DecodeBasisEB(int k, int *s, int *m);
int ConstructHamiltonEB(cfac_t *cfac, int n, int *ilev);
void StructureEB(cfac_t *cfac, char *fn, int n, int *ilev);
double HamiltonElementEB(const cfac_t *cfac, int i, int j);

int SlaterCoeff(cfac_t *cfac,
    char *fn, int nlevs, int *ilevs, int na, SHELL *sa, 
		int nb, SHELL *sb);
void AddSlaterCoeff(const cfac_t *cfac, double *c, double a, int n_shells, 
		    SHELL_STATE *sbra, SHELL_STATE *sket, 
		    INTERACT_SHELL *s, int na, SHELL *sa, 
		    int nb, SHELL *sb);

#endif
