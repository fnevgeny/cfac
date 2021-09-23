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

#ifndef _COULOMB_H_
#define _COULOMB_H_ 1

/*************************************************************
  Header for module "coulomb".
  This module calculates quatities related to the H-like ions.

  Author: M. F. Gu, mfgu@stanford.edu
**************************************************************/

void    SetHydrogenicNL(cfac_t *cfac, int n, int kl, int nm, int klm);
void    GetHydrogenicNL(const cfac_t *cfac, int *n, int *kl, int *nm, int *klm);

double  HydrogenicDipole(const cfac_t *cfac, int n0, int kl0, int n1, int kl1);
double  HydrogenicExpectation(double z, int m, int n, int kl);
double  HydrogenicSelfEnergy(double z, int n, int k);
double  CoulombPhaseShift(double z, double e, int kappa);
double *GetCoulombBethe(const cfac_cbcache_t *cbcache,
    int ie2, int ite, int ie1, int t, int q);
double  GetCoulombBetheAsymptotic(double te, double e1);
int     CoulombBetheTail(int n, double *w, int nkl, double *kl, double *tcb);
int     PrepCoulombBethe(cfac_cbcache_t *cbcache,
                         int ne2, int nte, int ne1, double z,
                         double *e2, double *te, double *e1,
                         int nkl, double *kl, int mode);
#endif
