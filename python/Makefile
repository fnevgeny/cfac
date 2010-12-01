TOP = ..

include $(TOP)/Make.conf


all: mpy

mpy: 	mpy.c $(FACLIBS)
	$(CC) -o mpy -I$(TOP) -I$(TOP)/faclib $(CPPFLAGS) $(PYINC) $(CFLAGS) mpy.c \
	        $(LDFLAGS) $(FACLIBS) $(LIBS) $(PYLIB)

install:
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) mpy $(BINDIR)

clean:
	rm -rf *.o *~ mpy
