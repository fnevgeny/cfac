# Process this file with autoconf to produce a configure script.
AC_INIT

AC_CONFIG_HEADERS([sysdef.h:sysdef.h.in])
AC_PREFIX_DEFAULT("/usr/local")

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CC
AC_PROG_F77(gfortran g77 f77 fort77 xlf f90 xlf90)
AC_F77_LIBRARY_LDFLAGS

AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_RANLIB

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
    CFLAGS="$CFLAGS -fno-common -Wall -Wextra -pedantic -Wno-sign-compare -Wno-unused-parameter"
fi
if test -n "$G77"
then
    FFLAGS="$FFLAGS -fno-common -Wall -Wno-unused"
fi

# which F77 known by cfortran.h is used
AC_ARG_WITH(cfortran, 
	[  --with-cfortran=xxFortran    Macro to be defined for cfortran.h],
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

# setup libs 
LIBS="$LIBS -lgslodeiv2 -lgsl -lgslcblas -llapack -lblas -lsqlite3"
fwithm=`echo "$FLIBS " | grep " -lm"`
if test "x$fwithm" = "x"
then
  LIBS="$LIBS -lm"
fi

if test "x$extrainc" != "x"
then
  CPPFLAGS="$CPPFLAGS $extrainc"
fi
if test "x$extralib" != "x"
then
  LDFLAGS="$LDFLAGS $extralib"
fi

# add fopt
FFLAGSNO="$FFLAGS"
AC_SUBST(FFLAGSNO)
FFLAGS="$FFLAGS $fopt"
      
AC_CONFIG_FILES([Make.conf])

AC_OUTPUT
