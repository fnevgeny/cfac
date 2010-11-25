#####################################################
# Top-level Makefile for cFAC                       #
#####################################################
# You should not change anything here.              #
#####################################################

TOP = .

include Make.conf

subdirs : configure Make.conf
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE)) || exit 1; done

all : subdirs


pfac:   $(FACLIBS)
	${PYTHON} setup.py build --force -extracomp="${CFLAGS}" -extralink="$(FACLIBS) ${LIBS}"
#mpy:
#	${PYTHON} setup.py build -mpy -extracomp="${CFLAGS}" -extralink="$(FACLIBS) ${LIBS}"

install : subdirs
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) install) || exit 1; done

check : subdirs
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) check) || exit 1; done

clean :
	@set -e; for i in $(SUBDIRS); do (cd $$i; $(MAKE) clean) || exit 1; done
	$(RM) -r build

distclean : clean
	$(RM) config.log config.status config.cache include/config.h Make.conf
	$(RM) -r autom4te.cache

devclean : distclean
	$(RM) configure NEWS ChangeLog

texts : ChangeLog

ChangeLog : dummy
	cvs2cl -F trunk

Make.conf : Make.conf.in configure
	@echo
	@echo 'Please re-run ./configure'
	@echo
	@exit 1

configure : configure.ac
	autoconf -o $@ configure.ac
	chmod +x $@

dummy :
