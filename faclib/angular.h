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

#ifndef _ANGULAR_H_
#define _ANGULAR_H_ 1

/*************************************************************
  Header for module "angular".
  This module calculates Wigner 3j, 6j, 9j symbols,
  and related vector coupling coefficients.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

/*
** Public functions provided by *angular*
*/
double LnFactorial(unsigned int n);
double LnInteger(unsigned int n);

int    Triangle(int j1, int j2, int j3);
double W3j(int j1, int j2, int j3, int m1, int m2, int m3);
double W6j(int j1, int j2, int j3, int i1, int i2, int i3);
int    W6jTriangle(int j1, int j2, int j3, int i1, int i2, int i3);
double W9j(int j1, int j2, int j3,
           int i1, int i2, int i3,
           int k1, int k2, int k3);
int    W9jTriangle(int j1, int j2, int j3,
                   int i1, int i2, int i3,
                   int k1, int k2, int k3);
double WignerEckartFactor(int jf, int k, int ji,
                          int mf, int q, int mi);
double ClebschGordan(int j1, int m1, int j2, int m2, int jf, int mf);
double ReducedCL(int ja, int k, int jb);
double WignerDMatrix(double a, int j2, int m2, int n2);

#endif
