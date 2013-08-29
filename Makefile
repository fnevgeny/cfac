PROG =	fdemo

LIBCFACDB = libcfacdb.a

LSRCS = cfacdb.c

FSRCS = fdemo.f
CSRCS = 
CHDRS = 

SCHEMAS = 

LOBJS = ${LSRCS:.c=.o}
FOBJS =	${FSRCS:.f=.o}

LIBS = -lsqlite3 -lgsl -lgslcblas

SRCS  = $(LSRCS) $(CSRCS)

include Make.conf

CFLAGS = $(DEBUG) $(LINT) $(OPTIMIZE) $(TCOVERAGE) $(PROFILING)
FFLAGS = $(DEBUG) $(OPTIMIZE)
LDFLAGS = $(DEBUG) 

all: $(PROG) $(SCHEMAS)

$(LIBCFACDB): $(LOBJS)
	$(RM) $@ && $(AR) $@ ${LOBJS} && $(RANLIB) $@

fdemo: $(FOBJS) $(LIBCFACDB)
	$(FC) $(LDFLAGS) -o $@ $(FOBJS) $(LIBCFACDB) $(LIBS)

cfac_schema.i : cfac_schema.sql
	./sql2cstr.sh < $? > $@

parse_cfac.o : cfac_schema.i

include Make.dep

Make.dep: $(SRCS) cfac_schema.i
	@echo -n "Generating dependencies... "
	@echo "# Generated automatically by \`make depend'" > $@
	@$(CC) $(CFLAGS) -MM $(SRCS) >> $@
	@echo "done"

install: $(PROG)
	install $(PROG) $(DESTBIN)

clean:
	rm -f $(PROG) $(SCHEMAS) $(LIBCFACDB) \
	$(LOBJS) $(FOBJS) cfac_schema.i Make.dep \
	tags ChangeLog *.bak \
	*.bb *.bbg *.da *.gcda *.gcno *.gcov

depend: Make.dep

tags: $(SRCS) $(CHDRS)
	ctags -f $@ $(SRCS) $(CHDRS)
