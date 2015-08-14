TOP  = ..

include $(TOP)/Make.conf

ALL_CFLAGS = $(CPPFLAGS) -I$(TOP)/include $(CFLAGS)

.c.o: 
	$(CC) -c $(ALL_CFLAGS) $<

PROGS = cfacdbu	fdemo

# Suffixes
.SUFFIXES : .c .f .sql .i

# Rules
.sql.i :
	./sql2cstr.sh $< > $@

LIBCFACDB = libcfacdb.a

LSRCS = cfacdb.c rates.c cfacdb_f.c
SQLS  = cfac_schema.sql cfac_schema_v1.sql cfac_schema_v2.sql

CSRCS = cfacdbu.c
CHDRS = cfacdb.h cfacdbP.h
FSRCS = fdemo.f

LOBJS = ${LSRCS:.c=.o}
COBJS =	${CSRCS:.c=.o}
FOBJS =	${FSRCS:.f=.o}

SQLIS =	${SQLS:.sql=.i}

LIBS = -lsqlite3 -lgsl -lgslcblas -lm

SRCS  = $(LSRCS) $(CSRCS)

all: $(PROGS)

$(LIBCFACDB): $(LOBJS)
	$(RM) $@ && $(AR) cr $@ ${LOBJS} && $(RANLIB) $@

cfacdbu: $(COBJS) $(LIBCFACDB)
	$(CC) $(LDFLAGS) -o $@ $(COBJS) $(LIBCFACDB) $(LIBS)

fdemo: $(FOBJS) $(LIBCFACDB)
	$(FC) $(LDFLAGS) -o $@ $(FOBJS) $(LIBCFACDB) $(LIBS)

include Make.dep

Make.dep: $(SRCS) $(SQLIS)
	@echo -n "Generating dependencies... "
	@echo "# Generated automatically by \`make depend'" > $@
	@$(CC) $(ALL_CFLAGS) -MM $(SRCS) >> $@
	@echo "done"

install: cfacdbu $(TOP)/include/cfacdb.h $(LIBCFACDB)
	$(INSTALLDIR) $(BINDIR)
	$(INSTALL) cfacdbu $(BINDIR)
	$(INSTALLDIR) $(INCDIR)
	$(INSTALLDATA) $(TOP)/include/cfacdb.h $(INCDIR)
	$(INSTALLDIR) $(LIBDIR)
	$(INSTALLDATA) $(LIBCFACDB) $(LIBDIR)

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
