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

#ifndef _EXCITATION_H_
#define _EXCITATION_H_

#include "consts.h"
#include "transition.h"

int InitExcitation(void);

int SetCETEGrid(int n, double emin, double emax);
int SetCETEGridDetail(int n, double *x);
int SetAngleGrid(int m, int n, double xmin, double xmax);
int SetAngleGridDetail(int m, int n, double *xg);
int SetCEBorn(double e, double x, double x1, double x0);
void SetCELQR(int m);
void SetCELMax(int m);
void SetCELCB(int m);
int SetCEPWOptions(int qr, int max, int kl_cb);
int AddCEPW(int n, int step);
int SetCEFormat(int m);
int SetCEPWGrid(int ns, int *n, int *step);
int SetCEEGridLimits(double min, double max, int type);
int SetCEEGridType(int type);
int SetUsrCEEGridType(int type);
int SetCEPWGridType(int type);
int SetCEEGridDetail(int n, double *x);
int SetCEEGrid(int n, double emin, double emax, double eth);
int SetUsrCEEGridDetail(int n, double *x);
int SetUsrCEEGrid(int n, double emin, double emax, double eth);

int CERadialQkBorn(cfac_t *cfac, int k0, int k1, int k2, int k3, int k, 
		   double te, double e1, double *qk, int m);

int SaveExcitation(cfac_t *cfac, int nlow, int *low, int nup, int *up, int msub, char *fn);
int SaveExcitationEB(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn);
int SaveExcitationEBD(cfac_t *cfac, int nlow, int *low, int nup, int *up, char *fn);

#endif
