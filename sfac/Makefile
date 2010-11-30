TOP = ..

include $(TOP)/Make.conf

ALL_CFLAGS = $(CPPFLAGS) -I$(TOP) -I$(TOP)/faclib $(CFLAGS)

.c.o: 
	$(CC) -c $(ALL_CFLAGS) $<

PROGS = sfac scrm spol

SRCS = stoken.c sfac.c scrm.c spol.c

all : $(PROGS)

include Make.dep

Make.dep : $(SRCS)
	echo "# Generated automatically by \`make depend'" > $@
	$(CC) $(ALL_CFLAGS) -MM $(SRCS) >> $@

sfac : 	sfac.o stoken.o $(FACLIBS)
	$(CC) -o sfac sfac.o stoken.o $(FACLIBS) $(LDFLAGS) $(LIBS)
scrm : 	scrm.o stoken.o $(FACLIBS)
	$(CC) -o scrm scrm.o stoken.o $(FACLIBS) $(LDFLAGS) $(LIBS)
spol : 	spol.o stoken.o $(FACLIBS)
	$(CC) -o spol spol.o stoken.o $(FACLIBS) $(LDFLAGS) $(LIBS)

install :
	$(INSTALLDIR) $(BINDIR)
	$(INSTALL) $(PROGS) $(BINDIR)

clean :
	$(RM) *.o *~ $(PROGS)

check : dummy

dummy :
