PROG =	fdemo

# Suffixes
.SUFFIXES : .c .f .sql .i

# Rules
.sql.i :
	./sql2cstr.sh $< > $@

LIBCFACDB = libcfacdb.a

LSRCS = cfacdb.c rates.c cfacdb_f.c
SQLS  = cfac_schema.sql cfac_schema_v1.sql cfac_schema_v2.sql

FSRCS = fdemo.f
CSRCS = 
CHDRS = 

LOBJS = ${LSRCS:.c=.o}
FOBJS =	${FSRCS:.f=.o}

SQLIS =	${SQLS:.sql=.i}

LIBS = -lsqlite3 -lgsl -lgslcblas

SRCS  = $(LSRCS) $(CSRCS)

include Make.conf

CFLAGS = $(DEBUG) $(LINT) $(OPTIMIZE) $(TCOVERAGE) $(PROFILING)
FFLAGS = $(DEBUG) $(OPTIMIZE)
LDFLAGS = $(DEBUG) 

all: $(PROG)

$(LIBCFACDB): $(LOBJS)
	$(RM) $@ && $(AR) $@ ${LOBJS} && $(RANLIB) $@

fdemo: $(FOBJS) $(LIBCFACDB)
	$(FC) $(LDFLAGS) -o $@ $(FOBJS) $(LIBCFACDB) $(LIBS)

parse_cfac.o : $(SQLIS)

include Make.dep

Make.dep: $(SRCS) $(SQLIS)
	@echo -n "Generating dependencies... "
	@echo "# Generated automatically by \`make depend'" > $@
	@$(CC) $(CFLAGS) -MM $(SRCS) >> $@
	@echo "done"

install: $(PROG)
	install $(PROG) $(DESTBIN)

clean:
	rm -f $(PROG) $(LIBCFACDB) \
	$(LOBJS) $(FOBJS) $(SQLIS) Make.dep \
	tags ChangeLog *.bak \
	*.bb *.bbg *.da *.gcda *.gcno *.gcov

depend: Make.dep

tags: $(SRCS) $(CHDRS)
	ctags -f $@ $(SRCS) $(CHDRS)
# DO NOT DELETE

rates.o: /usr/include/gsl/gsl_errno.h /usr/include/stdio.h
rates.o: /usr/include/features.h /usr/include/sys/cdefs.h
rates.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
rates.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
rates.o: /usr/include/bits/typesizes.h /usr/include/libio.h
rates.o: /usr/include/_G_config.h /usr/include/wchar.h
rates.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
rates.o: /usr/include/errno.h /usr/include/bits/errno.h
rates.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
rates.o: /usr/include/asm-generic/errno.h
rates.o: /usr/include/asm-generic/errno-base.h /usr/include/gsl/gsl_types.h
rates.o: /usr/include/gsl/gsl_integration.h /usr/include/stdlib.h
rates.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
rates.o: /usr/include/endian.h /usr/include/bits/endian.h
rates.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
rates.o: /usr/include/time.h /usr/include/sys/select.h
rates.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
rates.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
rates.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
rates.o: /usr/include/gsl/gsl_math.h /usr/include/math.h
rates.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
rates.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
rates.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
rates.o: /usr/include/bits/mathcalls.h /usr/include/gsl/gsl_sys.h
rates.o: /usr/include/gsl/gsl_inline.h /usr/include/gsl/gsl_machine.h
rates.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
rates.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
rates.o: /usr/include/bits/posix2_lim.h /usr/include/gsl/gsl_precision.h
rates.o: /usr/include/gsl/gsl_nan.h /usr/include/gsl/gsl_pow_int.h
rates.o: /usr/include/gsl/gsl_minmax.h /usr/include/gsl/gsl_const_num.h
rates.o: cfacdb.h /usr/include/sqlite3.h
