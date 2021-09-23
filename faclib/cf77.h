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

#ifndef _CF77_H_
#define _CF77_H_ 1

#include "sysdef.h"

#ifdef WITH_CPC_ACCEPTED

#include "cfortran.h"

     /* dirac coulomb function */
     PROTOCCALLSFSUB9(DCOUL, dcoul, DOUBLE, DOUBLE, INT, DOUBLE, DOUBLEV,\
                      DOUBLEV, DOUBLEV, DOUBLEV, INTV)
#define DCOUL(A1,A2,A3,A4,A5,A6,A7,A8,A9)                               \
     CCALLSFSUB9(DCOUL, dcoul, DOUBLE, DOUBLE, INT, DOUBLE, DOUBLEV,\
                 DOUBLEV, DOUBLEV, DOUBLEV, INTV,\
                 A1,A2,A3,A4,A5,A6,A7,A8,A9)

     /* coulomb multipole matrix elements */
     PROTOCCALLSFSUB10(CMULTIP, cmultip, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,\
                       INT, INT, INT, DOUBLEV, INTV)
#define CMULTIP(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)                   \
     CCALLSFSUB10(CMULTIP, cmultip, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,\
                  INT, INT, INT, DOUBLEV, INTV,                 \
                  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)

#else

#define DCOUL cfac_dummy_dcoul
#define CMULTIP cfac_dummy_cmultip

void cfac_dummy_dcoul(double, double, int, double, double *,\
                 double *, double *, double *, int *);
void cfac_dummy_cmultip(double, double, double, double, double,\
                  int, int, int, double *, int *);

#endif /* ACCEPT_CPC */

#endif
