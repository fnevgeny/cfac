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
    INT_Q1P2      = 7     /* Q1*P2         */
} RadIntType;

typedef struct _SLATER_YK_ {
  short npts;
  float *yk;
  float coeff[2];
} SLATER_YK;

void SetSlaterCut(cfac_t *cfac, int k0, int k1);
void SetSE(cfac_t *cfac, int n);
void SetVP(cfac_t *cfac, int n);
void SetBreit(cfac_t *cfac, int n);
void SetMS(cfac_t *cfac, int nms, int sms);
int SetAWGrid(cfac_t *cfac, int n, double min, double max);
int SetRadialGrid(cfac_t *cfac,
    int maxrp, double ratio, double asymp, double rmin);
int GetPotential(const cfac_t *cfac, char *s);
double GetResidualZ(const cfac_t *cfac);
double GetRMax(cfac_t *cfac);

int WaveFuncTable(cfac_t *cfac, char *s, int n, int kappa, double e);

/* get the index of the given orbital in the table */
int OrbitalIndex(cfac_t *cfac, int n, int kappa, double energy);
int OrbitalExists(const cfac_t *cfac, int n, int kappa, double energy);
ORBITAL *GetOrbital(const cfac_t *cfac, int k);
int GetNumBounds(const cfac_t *cfac);
int GetNumOrbitals(const cfac_t *cfac);
int GetNumContinua(const cfac_t *cfac);

double GetPhaseShift(cfac_t *cfac, int k);

/* radial optimization */
int SetAverageConfig(cfac_t *cfac, int nshells, int *n, int *kappa, double *nq);
void SetOptimizeMaxIter(cfac_t *cfac, int m);
void SetOptimizeStabilizer(cfac_t *cfac, double m);
void SetOptimizeTolerance(cfac_t *cfac, double c);
void SetOptimizePrint(cfac_t *cfac, int m);
void SetOptimizeControl(cfac_t *cfac, double tolerence, double stablizer,
                        int maxiter, int iprint);
void SetScreening(cfac_t *cfac, int n_screen, int *screened_n,
                  double screened_harge, int kl);
int OptimizeRadial(cfac_t *cfac, int ng, int *kg, double *weight);
int RefineRadial(cfac_t *cfac, int maxfun, int msglvl);
int ConfigEnergy(cfac_t *cfac, int m, int mr, int ng, int *kg);
double TotalEnergyGroup(cfac_t *cfac, int kg);
double ZerothEnergyConfig(cfac_t *cfac, CONFIG *cfg);
double ZerothResidualConfig(cfac_t *cfac, CONFIG *cfg);
double AverageEnergyConfig(cfac_t *cfac, CONFIG *cfg);

/* routines for radial integral calculations */
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
double Vinti(cfac_t *cfac, int k0, int k1);
double QED1E(cfac_t *cfac, int k0, int k1);
double SelfEnergyRatio(POTENTIAL *potential, ORBITAL *orb);
int Slater(const cfac_t *cfac, double *s, int k0, int k1, int k2, int k3, int k, int mode);
double BreitC(cfac_t *cfac, int n, int m, int k, int k0, int k1, int k2, int k3);
double BreitS(cfac_t *cfac, int k0, int k1, int k2, int k3, int k);
double BreitI(cfac_t *cfac, int n, int k0, int k1, int k2, int k3, int m);
double Breit(cfac_t *cfac, int k0, int k1, int k2, int k3, int k,
             int kl0, int kl1, int kl2, int kl3);
void SortSlaterKey(int *kd);
int ResidualPotential(const cfac_t *cfac, double *s, int k0, int k1);
double MeanPotential(cfac_t *cfac, int k0, int k1);
int FreeMultipoleArray(cfac_t *cfac);
int FreeMomentsArray(cfac_t *cfac);
int FreeGOSArray(cfac_t *cfac);

double RadialMoments(const cfac_t *cfac, int m, int k1, int k2);
double MultipoleRadialNR(cfac_t *cfac, int m, int k1, int k2, int guage);
double MultipoleRadialFR(cfac_t *cfac, double aw, int m, int k1, int k2, int guage);
double InterpolateMultipole(double aw2, int n, double *x, double *y);
double *GeneralizedMoments(cfac_t *cfac, int k0, int k1, int m);
void PrintGeneralizedMoments(cfac_t *cfac,
    char *fn, int m, int n0, int k0, int n1, int k1, double e1);
int FreeOrbital(cfac_t *cfac, int i);
int SaveAllContinua(cfac_t *cfac, int mode);
int SaveContinua(cfac_t *cfac, double e, int mode);
int FreeAllContinua(cfac_t *cfac);
int FreeContinua(cfac_t *cfac, double e);
int ClearOrbitalTable(cfac_t *cfac, int m);
int InitRadial(cfac_t *cfac);
int ReinitRadial(cfac_t *cfac, int m);

#endif

