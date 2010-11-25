TOP = ..

include $(TOP)/Make.conf

ALL_CFLAGS = $(CPPFLAGS) -I$(TOP) -I$(TOP)/faclib $(CFLAGS) $(LDFLAGS) $(LIBS)

PROGS = sfac scrm spol

include Make.dep

all : $(PROGS)

Make.dep : $(SRCS)
	echo "# Generated automatically by \`make depend'" > $@
	$(CC) $(ALL_CFLAGS) -MM $(SRCS) \
		| sed 's/\.o:/\$$(O)\ :/g' >> $@


sfac : 	sfac.c stoken.c $(FACLIBS)
	$(CC) -o sfac sfac.c stoken.c $(FACLIBS) $(ALL_CFLAGS)
scrm : 	scrm.c stoken.c $(FACLIBS)
	$(CC) -o scrm scrm.c stoken.c $(FACLIBS) $(ALL_CFLAGS)
spol : 	spol.c stoken.c $(FACLIBS)
	$(CC) -o spol spol.c stoken.c $(FACLIBS) $(ALL_CFLAGS)

install :
	$(INSTALLDIR) $(BINDIR)
	$(INSTALL) $(PROGS) $(BINDIR)

clean :
	rm -f *.o *~ $(PROGS)

check : dummy

dummy :
