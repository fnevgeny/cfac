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

#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

void SVDFit(int np, double *coeff, double tol,
	    int nd, double *x, double *logx, double *y, double *sig,
	    void Basis(int, double *, double, double));
int NLSQFit(int np, double *p, double tol,
	    double *fvec, double *fjac,
	    int n, double *x, double *logx, double *y, double *sigma,
	    void Func(int np, double *p, int n, double *x, double *logx, 
		     double *y, double *dy, int ndy, void *extra), 
	    void *extra);
double Simpson(const double *y, int ia, int ib);
int NewtonCotes(double r[], const double x[], int ilast,
                int last_only, int id);

#define UVIP3P uvip3p
void uvip3p(int nd, const double *xd, const double *yd,
		 int ni, const double *xi, double *yi);

#define UVIP3C uvip3c
void uvip3c(int nd, const double xd[], const double yd[],
		 double c1[], double c2[], double c3[]);

#endif
