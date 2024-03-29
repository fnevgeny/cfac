# Process this file with autoconf to produce a configure script.
AC_INIT

AC_CONFIG_AUX_DIR(ac-tools)

AC_CONFIG_HEADERS([sysdef.h:ac-tools/sysdef.h.in])
AC_PREFIX_DEFAULT("/usr/local")

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_F77(gfortran g77 f77 fort77 xlf f90 xlf90)
AC_F77_LIBRARY_LDFLAGS

# .exe file extension
AC_EXEEXT

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

AC_CHECK_PROG(AR, ar, ar, NONE)
if test "$AR" = "NONE"; then
  AC_MSG_ERROR(--> Can't find \`ar'!)
  AC_CACHE_SAVE
  exit 1
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_HEADER_STDC

#get optimization option for C and f77 compilers
AC_ARG_ENABLE(copt,
	[  --enable-copt=opt       optimization option for C compiler],
	[copt="$enableval"])
AC_ARG_ENABLE(fopt,
	[  --enable-fopt=opt       optimization option for F77 Compiler],
	[fopt="$enableval"])

# override the default CFLAGS and FFLAGS for GCC and G77
if test -n "$GCC"
then
  if test -z "$copt"
  then 
    copt="-O2"
  fi
fi
if test -n "$G77"
then
  if test -z "$fopt"
  then 
    fopt="-O2"
  fi
fi
CFLAGS="$copt"
FFLAGS=""

# Debugging option
AC_ARG_ENABLE(cdebug,
	[  --enable-cdebug=opt     Debugging option for C],
	[cdebug_opt=$enableval])
if test -n "$cdebug_opt"
then
  CFLAGS="$CFLAGS $cdebug_opt"
fi
AC_ARG_ENABLE(fdebug,
	[  --enable-fdebug=opt     Debugging option for F77],
	[fdebug_opt=$enableval])
if test -n "$fdebug_opt"
then
  FFLAGS="$FFLAGS $fdebug_opt"
fi
  
# gcc warnings
if test -n "$GCC"
then
    CFLAGS="$CFLAGS -fno-common -Wall -Wextra -pedantic"
    CFLAGS="$CFLAGS -Wno-sign-compare -Wno-unused-parameter"
fi
if test -n "$G77"
then
    FFLAGS="-std=legacy $FFLAGS -fno-common -Wall -Wextra"
    FFLAGS="$FFLAGS -Wno-compare-reals -Wno-unused-dummy-argument -Wno-unused-parameter"
fi

# use CPC-licensed code
AC_ARG_WITH(cpc_modules,
[  --with-cpc-modules      compile in CPC-licensed modules [[no]]],
[case "${withval}" in
  yes) cpc_modules=true ;;
  no)  cpc_modules=false;;
  *) AC_MSG_ERROR(bad value ${withval} for --with-cpc-modules) ;;
esac],[cpc_modules=false])

if test $cpc_modules = true
then
    echo ""
    cat coul/CPC-LICENSE | more
    answer="?"
    while test "x$answer" != "xyes" -a "x$answer" != "xno"
    do
        echo "Do you agree to these terms (type \"yes\" or \"no\")?"
        read answer <&AS_ORIGINAL_STDIN_FD
    done
    
    if test "x$answer" = "xyes"
    then
        echo ""
        echo "THE RESULTING BINARIES WILL NOT BE DISTRIBUTABLE!!!"
        echo ""
        echo "Press any key to continue"
        read answer <&AS_ORIGINAL_STDIN_FD
    else
        cpc_modules=false
    fi
fi

if test $cpc_modules = true
then
  AC_DEFINE(WITH_CPC_ACCEPTED)
fi

# which F77 known by cfortran.h is used
AC_ARG_WITH(cfortran, 
	[  --with-cfortran=Fortran Macro to be defined for cfortran.h],
	[xxFortran="$withval"])
if test -z "$xxFortran"
then
  if test -n "$G77"
  then 
    xxFortran=f2cFortran
  else
    case "$host_os" in
      *linux* )
        xxFortran=f2cFortran
      ;;
      *solaris*)
	xxFortran=sunFortran
      ;;
      *sunos*)
	xxFortran=sunFortran
      ;;
      *hpux*)
	xxFortran=hpuxFortran
      ;;
      *)
        echo "NOTE: F77 known by cfortran.h is not detected."
        echo "      using default f2cFortran."
        echo "      if not, change it using --with-cfortran=xxFortran"
        echo "      see faclib/cfortran.doc for available options."
        xxFortran=f2cFortran
      ;;
    esac
  fi
fi
AC_DEFINE_UNQUOTED([$xxFortran])

AC_ARG_WITH(extrainc,
	[  --with-extrainc         Extra compile flags],
	[extrainc=$withval])
AC_ARG_WITH(extralib,
	[  --with-extralib         Extra Libs],
	[extralib=$withval])

if test "x$extrainc" != "x"
then
  CPPFLAGS="$CPPFLAGS $extrainc"
fi
if test "x$extralib" != "x"
then
  LDFLAGS="$LDFLAGS $extralib"
fi

# libm - don't look for it again if $FLIBS already include it
fwithm=`echo "$FLIBS " | grep " -lm"`
if test "x$fwithm" = "x"
then
  AC_CHECK_LIB([m],[cos])
fi

AC_CHECK_DECLS([isfinite], [], [], [[#include <math.h>]])

# GSL libs 
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_odeiv2_driver_apply],[],
             [AC_MSG_ERROR(could not find required version of GSL)])

# Sqlite3
AC_CHECK_LIB([sqlite3],[sqlite3_open],[],
             [AC_MSG_ERROR(could not find SQLite3 library)])

# add fopt
FFLAGSNO="$FFLAGS"
AC_SUBST(FFLAGSNO)
FFLAGS="$FFLAGS $fopt"

AC_SUBST(EXEEXT)

AC_CONFIG_FILES([Make.conf:ac-tools/Make.conf.in])

AC_OUTPUT
