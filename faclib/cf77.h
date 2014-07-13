#ifndef _CF77_H_
#define _CF77_H_ 1

#include "sysdef.h"
#include "cfortran.h"

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
#endif
