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

#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include <stdio.h>

#define MAX_COMPLEX 512
typedef struct _REC_COMPLEX_ {
  int n;
  ARRAY *rg;
  int s0;
  int s1;
} REC_COMPLEX;

int InitRecombination(void);
int ReinitRecombination(int m);
int FreeRecPk(void);
int FreeRecQk(void);
int FreeRecAngZ(void);
int SetAICut(double c);
int SetRRTEGrid(const cfac_t *cfac, int n, double emin, double emax);
int SetRRTEGridDetail(const cfac_t *cfac, int n, double *x);
int SetPEGrid(const cfac_t *cfac, int n, double emin, double emax, double eth);
int SetPEGridDetail(int n, double *x);
int SetPEGridLimits(double min, double max, int type);
int SetUsrPEGridType(int type);
int SetUsrPEGrid(const cfac_t *cfac, int n, double emin, double emax, double eth);
int SetUsrPEGridDetail(const cfac_t *cfac, int n, double *x);
int AddRecPW(int n, int step);
int SetRecQkMode(int m, double tol);
int SetRecPWOptions(int kl_interp, int max_kl);
int SetRecPWLimits(int m1, int m2);
int SetRecSpectator(int n_frozen, int n_spec);
int ConstructRecGroupName(char *rgn, char *gn, int n);
int RecStates(cfac_t *cfac, int n, int k, int *kg, char *fn);
int RecStatesFrozen(cfac_t *cfac, int n, int k, int *kg, char *fn);
int RRRadialMultipoleTable(cfac_t *cfac, double *qr, int k0, int k1, int m);
int RRRadialQkTable(cfac_t *cfac, double *qr, int k0, int k1, int m);
int RRRadialMultipole(cfac_t *cfac, double *rqc, double te, int k0, int k1, int m);
int RRRadialQk(cfac_t *cfac, double *rqc, double te, int k0, int k1, int m);
void RRRadialQkFromFit(int np, double *p, int n, double *x, double *logx,
                       double *y, double *dy, int ndy, void *extra);
void RRRadialQkHydrogenicParams(int np, double *p, double z, int n, int klb);
int BoundFreeMultipole(cfac_t *cfac, FILE *fp, int rec, int f, int m);
int BoundFreeOS(cfac_t *cfac, double *rqu, double *p,
                double *eb, int rec, int f, int m, int iuta);
int PrepRREGrids(const cfac_t *cfac, double eth, double emax0);
int SaveRRMultipole(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn, int m);
int SaveRecRR(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn, int m);
int SaveAI(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn,
           int msub);
int AIRadial1E(cfac_t *cfac, double *pk, int kb, int kappaf);
int AIRadialPk(cfac_t *cfac, double **pk, int k0, int k1, int kb, int kappaf, int k);
int AutoionizeRate(cfac_t *cfac, double *rate, double *e, int rec, int f, int msub);

#endif
