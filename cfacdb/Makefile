TOP  = ..

include $(TOP)/Make.conf

ALL_CFLAGS = $(CPPFLAGS) -I$(TOP)/include $(CFLAGS)

.c.o: 
	$(CC) -c $(ALL_CFLAGS) $<

PROGS = cfacdbu$(EXE) fdemo$(EXE)

# Suffixes
.SUFFIXES : .c .f .sql .i

# Rules
.sql.i :
	./sql2cstr.sh $< > $@

LIBCFACDB = libcfacdb.a

LSRCS = cfacdb.c rates.c cache.c cfacdb_f.c
SQLS  = cfac_schema.sql cfac_schema_v1.sql cfac_schema_v2.sql cache_schema.sql

CSRCS = cfacdbu.c
CHDRS = cfacdb.h cfacdbP.h
FSRCS = fdemo.f

LOBJS = ${LSRCS:.c=.o}
COBJS =	${CSRCS:.c=.o}
FOBJS =	${FSRCS:.f=.o}

SQLIS =	${SQLS:.sql=.i}

SRCS  = $(LSRCS) $(CSRCS)

all: $(PROGS)

$(LIBCFACDB): $(LOBJS)
	$(RM) $@ && $(AR) cr $@ ${LOBJS} && $(RANLIB) $@

cfacdbu$(EXE): $(COBJS) $(LIBCFACDB)
	$(CC) $(LDFLAGS) -o $@ $(COBJS) $(LIBCFACDB) $(LIBS)

fdemo$(EXE): $(FOBJS) $(LIBCFACDB)
	$(FC) $(LDFLAGS) -o $@ $(FOBJS) $(LIBCFACDB) $(LIBS)

cfacdbu.1: cfacdbu$(EXE)
	help2man -s 1 -N ./cfacdbu$(EXE) -n "CFACDB utility" -S "cFAC" -o $@

include Make.dep

Make.dep: $(SRCS) $(SQLIS)
	@echo -n "Generating dependencies... "
	@echo "# Generated automatically by \`make depend'" > $@
	@$(CC) $(ALL_CFLAGS) -MM $(SRCS) >> $@
	@echo "done"

install: cfacdbu$(EXE) cfacdbu.1 $(TOP)/include/cfacdb.h $(LIBCFACDB)
	$(INSTALLDIR) $(BINDIR)
	$(INSTALL) cfacdbu$(EXE) $(BINDIR)
	$(INSTALLDIR) $(INCDIR)
	$(INSTALLDATA) $(TOP)/include/cfacdb.h $(INCDIR)
	$(INSTALLDIR) $(LIBDIR)
	$(INSTALLDATA) $(LIBCFACDB) $(LIBDIR)
	$(INSTALLDIR) $(MANDIR)/man1
	$(INSTALLDATA) cfacdbu.1 $(MANDIR)/man1

check :

clean:
	rm -f $(PROGS) $(LIBCFACDB) \
	$(LOBJS) $(COBJS) $(FOBJS) $(SQLIS) Make.dep \
	tags ChangeLog *.bak \
	*.bb *.bbg *.da *.gcda *.gcno *.gcov

depend: Make.dep

tags: $(SRCS) $(CHDRS)
	ctags -f $@ $(SRCS) $(CHDRS)

count:
	@sloccount . | \
	egrep --color 'Total ((Estimated Cost)|(Physical Source Lines)).*'
