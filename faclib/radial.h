#ifndef _RADIAL_H_
#define _RADIAL_H_

#include "config.h"
#include "orbital.h"

typedef enum {
    INT_P1P2pQ1Q2 = 1,    /* P1*P2 + Q1*Q2 */
    INT_P1P2      = 2,    /* P1*P2         */
    INT_Q1Q2      = 3,    /* Q1*Q2         */ 
    INT_P1Q2pQ1P2 = 4,    /* P1*Q2 + Q1*P2 */
    INT_P1Q2mQ1P2 = 5,    /* P1*Q2 - Q1*P2 */
    INT_P1Q2      = 6,    /* P1*Q2         */
    INT_Q1P2      = 7     /* P1*Q2         */
} RadIntType;

typedef struct _SLATER_YK_ {
  short npts;
  float *yk;
  float coeff[2];
} SLATER_YK;

int SetBoundary(cfac_t *cfac, int nmax, double p, double bqp);
int RadialOverlaps(cfac_t *cfac, char *fn, int kappa);
void SetSlaterCut(int k0, int k1);
void SetPotentialMode(cfac_t *cfac, int m, double h);
void SetSE(int n);
void SetVP(int n);
void SetBreit(int n);
void SetMS(int nms, int sms);
int SetAWGrid(int n, double min, double max);
int GetAWGrid(double **a);
int SetRadialGrid(cfac_t *cfac,
    int maxrp, double ratio, double asymp, double rmin);
double SetPotential(POTENTIAL *potential, AVERAGE_CONFIG *acfg, int iter, double *vbuf);
int GetPotential(const cfac_t *cfac, char *s);
double GetResidualZ(const cfac_t *cfac);
double GetRMax(cfac_t *cfac);

/* solve the dirac equation for the given orbital */
int SolveDirac(POTENTIAL *potential, ORBITAL *orb);
int WaveFuncTable(cfac_t *cfac, char *s, int n, int kappa, double e);

/* get the index of the given orbital in the table */
int OrbitalIndex(const cfac_t *cfac, int n, int kappa, double energy);
int OrbitalExists(POTENTIAL *potential, int n, int kappa, double energy);
int AddOrbital(ORBITAL *orb);
ORBITAL *GetOrbital(int k);
ORBITAL *GetOrbitalSolved(POTENTIAL *potential, int k);
ORBITAL *GetNewOrbital(void);
int GetNumBounds(void);
int GetNumOrbitals(void);
int GetNumContinua(void);

double GetPhaseShift(cfac_t *cfac, int k);

/* radial optimization */
int SetAverageConfig(int nshells, int *n, int *kappa, double *nq);
void SetOptimizeMaxIter(int m);
void SetOptimizeStabilizer(double m);
void SetOptimizeTolerance(double c);
void SetOptimizePrint(int m);
void SetOptimizeControl(double tolerence, double stablizer, 
			int maxiter, int iprint);
void SetScreening(int n_screen, int *screened_n, 
		  double screened_harge, int kl);
int OptimizeRadial(cfac_t *cfac, int ng, int *kg, double *weight);
int RefineRadial(cfac_t *cfac, int maxfun, int msglvl);
int ConfigEnergy(cfac_t *cfac, int m, int mr, int ng, int *kg);
double TotalEnergyGroup(const cfac_t *cfac, int kg);
double ZerothEnergyConfig(cfac_t *cfac, CONFIG *cfg);
double ZerothResidualConfig(cfac_t *cfac, CONFIG *cfg);
double AverageEnergyConfig(const cfac_t *cfac, CONFIG *cfg);
double AverageEnergyAvgConfig(cfac_t *cfac, AVERAGE_CONFIG *cfg);

/* routines for radial integral calculations */
int GetYk(POTENTIAL *potential, int k, double *yk, ORBITAL *orb1, ORBITAL *orb2, 
	  int k1, int k2, RadIntType type);
int IntegrateF(POTENTIAL *potential, const double *f, const ORBITAL *orb1, const ORBITAL *orb2,
    RadIntType type, double x[], int id);
int IntegrateS(POTENTIAL *potential, const double *f, const ORBITAL *orb1, const ORBITAL *orb2, 
    RadIntType type, double *r, int id);
int IntegrateSubRegion(POTENTIAL *potential, int i0, int i1, 
		       const double *f,
                       const ORBITAL *orb1, const ORBITAL *orb2,
		       RadIntType type, double *r, int m, int last_only);
int IntegrateSinCos(POTENTIAL *potential, int j, double *x, double *y, 
		    double *phase, double *dphase, 
		    int i0, double *r, int last_only);
int SlaterTotal(cfac_t *cfac,
    double *sd, double *se, int *js, int *ks, int k, int mode);
double Vinti(POTENTIAL *potential, int k0, int k1);
double QED1E(cfac_t *cfac, int k0, int k1);
double SelfEnergyRatio(POTENTIAL *potential, ORBITAL *orb);
int Slater(const cfac_t *cfac, double *s, int k0, int k1, int k2, int k3, int k, int mode);
double BreitC(POTENTIAL *potential, int n, int m, int k, int k0, int k1, int k2, int k3);
double BreitS(POTENTIAL *potential, int k0, int k1, int k2, int k3, int k);
double BreitI(POTENTIAL *potential, int n, int k0, int k1, int k2, int k3, int m);
double Breit(POTENTIAL *potential, int k0, int k1, int k2, int k3, int k,
	     int kl0, int kl1, int kl2, int kl3);
void SortSlaterKey(int *kd);
int ResidualPotential(const cfac_t *cfac, double *s, int k0, int k1);
double MeanPotential(cfac_t *cfac, int k0, int k1);
int FreeResidualArray(void);
int FreeMultipoleArray(void);
int FreeSlaterArray(void);
int FreeSimpleArray(MULTI *ma);
int FreeMomentsArray(void);
int FreeGOSArray(void);

double RadialMoments(const cfac_t *cfac, int m, int k1, int k2);
double MultipoleRadialNR(cfac_t *cfac, int m, int k1, int k2, int guage);
int MultipoleRadialFRGrid(POTENTIAL *potential, double **p, int m, int k1, int k2, int guage);
double MultipoleRadialFR(cfac_t *cfac, double aw, int m, int k1, int k2, int guage);
double InterpolateMultipole(double aw2, int n, double *x, double *y);
double *GeneralizedMoments(cfac_t *cfac, int k0, int k1, int m);
void PrintGeneralizedMoments(cfac_t *cfac,
    char *fn, int m, int n0, int k0, int n1, int k1, double e1);
int SaveOrbital(int i);
int RestoreOrbital(int i); 
int FreeOrbital(int i);
int SaveAllContinua(int mode); 
int SaveContinua(double e, int mode);
int FreeAllContinua(void);
int FreeContinua(double e);
int ClearOrbitalTable(cfac_t *cfac, int m);
void LimitArrayRadial(int m, double n);
int InitRadial(cfac_t *cfac);
int ReinitRadial(cfac_t *cfac, int m);
int TestIntegrate(cfac_t *cfac);

#endif

