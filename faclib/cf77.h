#ifndef _CF77_H_
#define _CF77_H_ 1

#include "cfortran.h"

     /* interpolation from TOMS */
     PROTOCCALLSFSUB7(UVIP3C, uvip3c, INT, INT, DOUBLEV, DOUBLEV,\
		      DOUBLEV, DOUBLEV, DOUBLEV)
#define UVIP3C(A1,A2,A3,A4,A5,A6,A7)\
     CCALLSFSUB7(UVIP3C, uvip3c, INT, INT, DOUBLEV, DOUBLEV,\
		 DOUBLEV, DOUBLEV, DOUBLEV, A1,A2,A3,A4,A5,A6,A7)

     PROTOCCALLSFSUB12(SUBPLX, subplx, ROUTINE, INT, DOUBLE,\
		       INT, INT, DOUBLEV, DOUBLEV, DOUBLEV,\
		       INTV, DOUBLEV, INTV, INTV)
#define SUBPLX(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)\
     CCALLSFSUB12(SUBPLX, subplx, ROUTINE, INT, DOUBLE,\
		  INT, INT, DOUBLEV, DOUBLEV, DOUBLEV, \
		  INTV, DOUBLEV, INTV, INTV,\
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)

     /* dirac coulomb function */
     PROTOCCALLSFSUB9(DCOUL, dcoul, DOUBLE, DOUBLE, INT, DOUBLE, DOUBLEV,\
		      DOUBLEV, DOUBLEV, DOUBLEV, INTV)
#define DCOUL(A1,A2,A3,A4,A5,A6,A7,A8,A9)				\
     CCALLSFSUB9(DCOUL, dcoul, DOUBLE, DOUBLE, INT, DOUBLE, DOUBLEV,\
		 DOUBLEV, DOUBLEV, DOUBLEV, INTV,\
                 A1,A2,A3,A4,A5,A6,A7,A8,A9)

     /* coulomb multipole matrix elements */
     PROTOCCALLSFSUB11(CMULTIP, cmultip, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,\
		       INT, INT, INT, DOUBLEV, INT, INTV)
#define CMULTIP(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)			  \
     CCALLSFSUB11(CMULTIP, cmultip, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,\
		  INT, INT, INT, DOUBLEV, INT, INTV,			\
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)

     PROTOCCALLSFSUB8(ACOFZ1, acofz1, DOUBLE, DOUBLE, INT, INT,\
		      DOUBLEV, DOUBLEV, INT, INT)
#define ACOFZ1(A1,A2,A3,A4,A5,A6,A7,A8)\
     CCALLSFSUB8(ACOFZ1, acofz1, DOUBLE, DOUBLE, INT, INT,\
		 DOUBLEV, DOUBLEV, INT, INT, A1,A2,A3,A4,A5,A6,A7,A8)

     PROTOCCALLSFSUB12(PIXZ1, pixz1, DOUBLE, DOUBLE, INT, INT,\
		       DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT,\
		       INT, INT, INT)
#define PIXZ1(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)\
     CCALLSFSUB12(PIXZ1, pixz1, DOUBLE, DOUBLE, INT, INT,\
		  DOUBLEV, DOUBLEV, DOUBLEV, DOUBLEV, INT,\
		  INT, INT, INT, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)

     /* minpack routines */
     PROTOCCALLSFSUB24(LMDER, lmder, ROUTINE, INT, INT, DOUBLEV,\
		       DOUBLEV, DOUBLEV, INT, DOUBLE, DOUBLE, DOUBLE,\
		       INT, DOUBLEV, INT, DOUBLE, INT, INTV, INTV,\
		       INTV, INTV, DOUBLEV, DOUBLEV, DOUBLEV,\
		       DOUBLEV, DOUBLEV)
#define LMDER(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,\
              B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12)\
     CCALLSFSUB24(LMDER, lmder, ROUTINE, INT, INT, DOUBLEV,\
		  DOUBLEV, DOUBLEV, INT, DOUBLE, DOUBLE, DOUBLE,\
		  INT, DOUBLEV, INT, DOUBLE, INT, INTV, INTV,\
		  INTV, INTV, DOUBLEV, DOUBLEV, DOUBLEV,\
		  DOUBLEV, DOUBLEV, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,\
		  B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12)

     /* BLAS, LAPACK */
     PROTOCCALLSFFUN5(DOUBLE, DDOT, ddot, INT, DOUBLEV, INT, DOUBLEV, INT)
#define DDOT(A1,A2,A3,A4,A5)\
     CCALLSFFUN5(DDOT, ddot, INT, DOUBLEV, INT, DOUBLEV, INT,\
		 A1,A2,A3,A4,A5)

     PROTOCCALLSFSUB4(DSCAL, dscal, INT, DOUBLE, DOUBLEV, INT)
#define DSCAL(A1,A2,A3,A4)\
     CCALLSFSUB4(DSCAL, dscal, INT, DOUBLE, DOUBLEV, INT,\
		 A1,A2,A3,A4)

     PROTOCCALLSFSUB11(DGEMV, dgemv, STRING, INT, INT, DOUBLE, DOUBLEV,\
		       INT, DOUBLEV, INT, DOUBLE, DOUBLEV, INT)
#define DGEMV(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)\
     CCALLSFSUB11(DGEMV, dgemv, STRING, INT, INT, DOUBLE, DOUBLEV,\
		  INT, DOUBLEV, INT, DOUBLE, DOUBLEV, INT,\
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11)

     PROTOCCALLSFSUB12(DSPEVD, dspevd, STRING, STRING, INT, DOUBLEV,\
		       DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, INTV, INT,\
		       INTV)
#define DSPEVD(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)\
     CCALLSFSUB12(DSPEVD, dspevd, STRING, STRING, INT, DOUBLEV,\
		  DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, INTV, INT,\
		  INTV, A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12)

     PROTOCCALLSFSUB9(DSPEV, dspev, STRING, STRING, INT, DOUBLEV,\
		      DOUBLEV, DOUBLEV, INT, DOUBLEV, INTV)
#define DSPEV(A1,A2,A3,A4,A5,A6,A7,A8,A9)\
     CCALLSFSUB9(DSPEV, dspev, STRING, STRING, INT, DOUBLEV,\
		 DOUBLEV, DOUBLEV, INT, DOUBLEV, INTV,\
		 A1,A2,A3,A4,A5,A6,A7,A8,A9)

     PROTOCCALLSFSUB14(DGEEV, dgeev, STRING, STRING, INT, DOUBLEV, INT,	\
		       DOUBLEV, DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, \
		       INT, INTV)
#define DGEEV(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)	\
     CCALLSFSUB14(DGEEV, dgeev, STRING, STRING, INT, DOUBLEV, INT,	\
		  DOUBLEV, DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV, \
		  INT, INTV, \
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)

     PROTOCCALLSFSUB8(DGESV, dgesv, INT, INT, DOUBLEV, INT,\
		      INTV, DOUBLEV, INT, INTV)
#define DGESV(A1,A2,A3,A4,A5,A6,A7,A8)\
     CCALLSFSUB8(DGESV, dgesv, INT, INT, DOUBLEV, INT,\
		 INTV, DOUBLEV, INT, INTV, A1,A2,A3,A4,A5,A6,A7,A8)

     PROTOCCALLSFSUB14(DGESDD, dgesdd, STRING, INT, INT, DOUBLEV, INT,\
		       DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV,\
		       INT, INTV, INTV)
#define DGESDD(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)\
     CCALLSFSUB14(DGESDD, dgesdd, STRING, INT, INT, DOUBLEV, INT,\
		  DOUBLEV, DOUBLEV, INT, DOUBLEV, INT, DOUBLEV,\
		  INT, INTV, INTV,\
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14)

     PROTOCCALLSFSUB10(DGBSV, dgbsv, INT, INT, INT, INT, DOUBLEV, INT,\
		       INTV, DOUBLEV, INT, INTV)
#define DGBSV(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)\
     CCALLSFSUB10(DGBSV, dgbsv, INT, INT, INT, INT, DOUBLEV, INT,\
		  INTV, DOUBLEV, INT, INTV,\
		  A1,A2,A3,A4,A5,A6,A7,A8,A9,A10)


     PROTOCCALLSFSUB5(NJFORM, njform, INT, INT, INTV, INTV, INTV)
#define NJFORM(A1,A2,A3,A4,A5)					\
     CCALLSFSUB5(NJFORM, njform, INT, INT, INTV, INTV, INTV, A1,A2,A3,A4,A5)

     PROTOCCALLSFSUB2(NJSUM, njsum, INTV, DOUBLEV)
#define NJSUM(A1,A2)					\
     CCALLSFSUB2(NJSUM, njsum, INTV, DOUBLEV, A1,A2)

     PROTOCCALLSFSUB3(CPYDAT, cpydat, INT, INTV, INT) 
#define CPYDAT(A1,A2,A3)\
     CCALLSFSUB3(CPYDAT, cpydat, INT, INTV, INT, A1,A2,A3)
     
     PROTOCCALLSFSUB0(FACTT, factt)
#define FACTT()        \
     CCALLSFSUB0(FACTT, factt)

#endif
