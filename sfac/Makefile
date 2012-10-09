TOP = ..

include $(TOP)/Make.conf

ALL_CFLAGS = $(CPPFLAGS) -I$(TOP) -I$(TOP)/include -I$(TOP)/faclib $(CFLAGS)

.c.o: 
	$(CC) -c $(ALL_CFLAGS) $<

PROGS = sfac

SRCS = stoken.c sfac.c

all : $(PROGS)

include Make.dep

Make.dep : $(SRCS)
	echo "# Generated automatically by \`make depend'" > $@
	$(CC) $(ALL_CFLAGS) -MM $(SRCS) >> $@

sfac : 	sfac.o stoken.o $(FACLIBS)
	$(CC) -o sfac sfac.o stoken.o $(FACLIBS) $(LDFLAGS) $(LIBS)

install :
	$(INSTALLDIR) $(BINDIR)
	$(INSTALL) $(PROGS) $(BINDIR)

clean :
	$(RM) *.o *~ $(PROGS)

check : dummy

dummy :
