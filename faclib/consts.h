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

#ifndef _CONSTS_H_
#define _CONSTS_H_ 1

#include <math.h>

/*
** MACRO:       IsOdd, IsEven
** PURPOSE:     determin if an integer is Odd or Even.
** INPUT:       {int x},
**              integer.
** RETURN:      {int},
**              0: false.
**              1: true.
** SIDE EFFECT:
** NOTE:
*/
#define IsOdd(x)  ((abs((x))&0x01)?1:0)
#define IsEven(x) ((abs((x))&0x01)?0:1)

/*
** MACRO:       Max, Min
** PURPOSE:     the larger and lesser of two numbers.
** INPUT:       {generic a},
**              number participating the comparison.
**              {generic b},
**              number participating the comparison.
** RETURN:      {generic},
**              the larger or lesser of a and b.
** SIDE EFFECT:
** NOTE:
*/
#define Max(a, b) (((a)>(b))?(a):(b))
#define Min(a, b) (((a)<(b))?(a):(b))

#define IsNan(a)  (!((a)>0) && !((a)<=0))

/*
** VARIABLE:    EPS1, ..., EPS30.
** TYPE:        macro constants.
** PURPOSE:     some small numbers.
** NOTE:
*/
#define EPS30 1E-30
#define EPS16 1E-16
#define EPS12 1E-12
#define EPS10 1E-10
#define EPS8  1E-08
#define EPS6  1E-06
#define EPS5  1E-05
#define EPS4  1E-04
#define EPS3  1E-03
#define EPS2  1E-02
#define EPS1  1E-01

/*
** VARIABLE:    TWO_PI
** TYPE:        macro constants, double
** PURPOSE:     2*PI
** NOTE:
*/
#define TWO_PI     (2*M_PI)

/*
** VARIABLE:    HARTREE_EV
** TYPE:        macro constant
** PURPOSE:     1 Hartree in eV.
** NOTE:
*/
#define HARTREE_EV 27.2113862

/*
** VARIABLE:    RATE_AU
** TYPE:        macro constants.
** PURPOSE:     atomic units of rate in 1/s
** NOTE:
*/
#define RATE_AU    4.13413733E16

/*
** VARIABLE:    AREA_AU20
** TYPE:        macro constants.
** PURPOSE:     atomic units of area in 10^-20 cm2
** NOTE:
*/
#define AREA_AU20  2.800285206E3

/*
** VARIABLE:    RBOHR
** TYPE:        macro constant.
** PURPOSE:     Bohr radius in cm.
** NOTE:
*/
#define RBOHR      5.29177211e-9

/*
** VARIABLE:    MBOHR
** TYPE:        macro constant.
** PURPOSE:     Bohr magneton in eV/Gauss.
** NOTE:
*/
#define MBOHR      5.78838183E-9

/*
** VARIABLE:    FINE_STRUCTURE_CONST, FINE_STRUCTURE_CONST2
** TYPE:        macro constants.
** PURPOSE:     fine structure constant and its square.
** NOTE:
*/
#define FINE_STRUCTURE_CONST  7.29735257E-3
#define FINE_STRUCTURE_CONST2 5.32513545E-5
#define AMU  1836.153

/* radial QK modes */
#define QK_DEFAULT    -1
#define QK_EXACT       0
#define QK_INTERPOLATE 1
#define QK_FIT         2
#define QK_CB          3
#define QK_DW          4
#define QK_BED         5

/* blocks for multi arrays */
#define MULTI_BLOCK2   128
#define MULTI_BLOCK3   64
#define MULTI_BLOCK4   25
#define MULTI_BLOCK5   15
#define MULTI_BLOCK6   10

/* orbital */
#define MAXRP      3000  /* maximum radial mesh */
#define DMAXRP     1200  /* default radial mesh points */
#define GRIDASYMP  36    /* no. points in one wavelength near infinity */
#define GRIDRATIO  1.1   /* ratio of successive mesh near origin */
#define GRIDRMIN   1E-6  /* starting point of the mesh is GRIDRMIN/Z */
#define ENERELERR  1E-5  /* relative energy error */
#define ENEABSERR  1E-8  /* absolute energy error */

/* config */
#define MCHSHELL           2048
#define MAX_SPEC_SYMBOLS   21
#define LEVEL_NAME_LEN     128
#define GROUP_NAME_LEN     64
#define MAX_GROUPS         600
#define MAX_SYMMETRIES     400
#define CONFIGS_BLOCK      1024
#define STATES_BLOCK       2048

/* radial */
#define ORBITALS_BLOCK     1024
#define OPTSTABLE          0.5
#define OPTTOL             1E-6
#define OPTNITER           128
#define OPTPRINT           0
#define QEDSE              5
#define QEDVP              2
#define QEDNMS             1
#define QEDSMS             1
#define QEDBREIT           5
#define RCOREMIN           50

/* structure */
#define MAX_HAMS           2000
#define MAX_HAMS2          (MAX_HAMS*MAX_HAMS)
#define LEVELS_BLOCK       1024
#define ANGZ_BLOCK         1024
#define ANGZxZ_BLOCK       8192
#define ANGZCUT            1E-5
#define MIXCUT             1E-5
#define MIXCUT2            1.0
#define NPRINCIPLE         2
#define MAXDN              3
#define MBCLOSE            8
#define MAXLEVEB           1000000

/* transition */
#define G_COULOMB          1
#define G_BABUSHKIN        2
#define M_FR               0
#define M_NR               1
#define DGAUGE             G_BABUSHKIN
#define DMODE              M_NR

/* recouple */
#define MAXRANK            20

/* coulomb */
#define NHYDROGEN          20
#define LHYDROGEN          10
#define NHYDROGENMAX       512
#define LHYDROGENMAX       20
#define DIPOLE_BLOCK       64
#define CBMULT             2
#define MAXNCB             ((CBMULT*(CBMULT+3))/2)
#define CBLMIN             15
#define CBLMAX             150

/* grid */
#define MAXNKL             50
#define MAXKL              512
#define MAXNUSR            30
#define MAXNE              20
#define MAXNTE             6
#define MAXNTHETA          30
#define MAXNPHI            60
#define TE_MAX_MIN         5.0

/* excitation */
#define NGOSK              256
#define EXCLQR             0
#define EXCLMAX            36
#define EXCLCB             36
#define XBORN              (-0.5)
#define XBORN1             (-1.0)
#define XBORN0             (0.25)
#define EBORN              100.0

/* ionization */
#define IONMAXK            6
#define IONLQR             0
#define IONLMAX            36
#define IONLEJEC           4
#define IONLCB             36

/* recombination */
#define RECNSPEC           8
#define RECNFROZEN         8
#define RECLMAX            12
#define AICUT              0.0

/* polarization */
#define MAXPOL             4 /* maximum multipol for polarization */

#endif
