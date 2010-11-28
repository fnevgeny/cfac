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
AC_CHECK_PROG([PYTHON], [python], [python])
AC_SUBST(PYTHON)

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

# get position independent code option for C and f77 compilers
AC_ARG_ENABLE(cpic,
	[  --enable-cpic=opt       PIC option for C compiler],
	[cpic="$enableval"])
AC_ARG_ENABLE(fpic,
	[  --enable-fpic=opt       PIC option for F77 Compiler],
	[fpic="$enableval"])

# override the default CFLAGS and FFLAGS for GCC and G77
if test -n "$GCC"
then
  if test -z "$copt"
  then 
    copt="-O2"
  fi
  if test -z "$cpic"
  then 
    cpic="-fPIC"
  fi
fi
if test -n "$G77"
then
  if test -z "$fopt"
  then 
    fopt="-O2"
  fi
  if test -z "$fpic"
  then
    fpic="-fPIC"
  fi
fi
CFLAGS="$copt $cpic"
FFLAGS="$fpic"

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
    CFLAGS="$CFLAGS -fno-common -Wall -Wextra -pedantic -Wno-sign-compare -Wno-unused"
fi
if test -n "$G77"
then
    FFLAGS="$FFLAGS -fno-common -Wall"
    FFLAGS=`echo $FFLAGS`
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

# MPI
AC_ARG_WITH(mpi,
	[  --with-mpi              Use MPI],
	[use_mpi=$withval])
AC_ARG_WITH(mpicompile,
	[  --with-mpicompile       MPI compile flags],
	[mpicompile=$withval])
AC_ARG_WITH(mpilink,
	[  --with-mpilink          MPI link flags],
	[mpilink=$withval])
AC_ARG_WITH(extrainc,
	[  --with-extrainc         Extra compile flags],
	[extrainc=$withval])
AC_ARG_WITH(extralib,
	[  --with-extralib         Extra Libs],
	[extralib=$withval])

# setup libs 
LIBS="$LIBS -lgsl -lgslcblas -llapack -lblas"
fwithm=`echo "$FLIBS " | grep " -lm"`
if test "x$fwithm" = "x"
then
  LIBS="$LIBS -lm"
fi

if test "x$use_mpi" != "x"
then
  AC_DEFINE_UNQUOTED([USE_MPI])
  if test "x$mpicompile" = "x"
  then 
    case "$use_mpi" in
      lam*)
        mpicompile=`mpicc -showme`
	mpicompile=[`echo $mpicompile | sed 's/^ *//'`]	
	mpicompile=[`echo $mpicompile | sed 's/^[^ ]* *//'`]
	mpicompile=[`echo $mpicompile | sed 's/-[Ll][^ ]*//g'`]
	mpicompile=[`echo $mpicompile | sed 's/ *$//'`]
      ;;
      *)
	echo "Supply MPI compile flags using --with-mpicompile"
      ;;
    esac
  fi
  if test "x$mpilink" = "x"
  then 
    case "$use_mpi" in
      lam*)
        mpilink=`mpicc -showme`
	mpilink=[`echo $mpilink | sed 's/^ *//'`]
	mpilink=[`echo $mpilink | sed 's/^[^ ]* *//'`]
	mpilink=[`echo $mpilink | sed 's/-[^Ll][^ ]*//g'`]
	mpilink=[`echo $mpilink | sed 's/ *$//'`]
      ;;
      *)
	echo "Supply MPI link flags using --with-mpilink"
      ;;
    esac
  fi
  CPPFLAGS="$CPPFLAGS $mpicompile"
  LIBS="$LIBS $mpilink"
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

if test -z "$PYTHON"
then
  echo "WARNING: NO PYTHON DETECTED ON THIS SYSTEM."
  echo "To install PFAC interface,"
  echo "download and install Python 2.0 or later"
  echo "You may still procede to install SFAC interface, though."
fi

AC_OUTPUT
