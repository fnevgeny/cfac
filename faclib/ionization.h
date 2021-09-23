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

#ifndef _IONIZATION_H_
#define _IONIZATION_H_

int InitIonization(cfac_t *cfac);

int SetIEGrid(int n, double emin, double emax);
int SetIEGridDetail(int n, double *x);
void SetCIBorn(int x);
void SetCILQR(int m);
void SetCILMax(int m);
void SetCILMaxEject(int m);
void SetCILCB(int m);
void SetCITol(double t);
int SetCIPWGrid(int ns, int *n, int *step);
int SetCIFormat(int m);
int SetCIEGrid(int n, double emin, double emax, double eth);
int SetCIEGridDetail(int n, double *x);
int SetCIFormat(int m);
int SetCIEGridLimits(double min, double max, int type);
int SetUsrCIEGridType(int type);
int SetUsrCIEGrid(int n, double emin, double emax, double eth);
int SetUsrCIEGridDetail(int n, double *x);
int SetCIQkMode(int m, double tol);
int SetCIMaxK(cfac_t *cfac, int k);

int SaveIonization(cfac_t *cfac, int nb, int *b, int nf, int *f, char *fn);
int SaveIonizationMSub(cfac_t *cfac, int nb, int *b, int nf, int *f, char *fn);

#endif
