TOP = ..

DOCS = manual.pdf

CHAPTERS = overview.tex physics.tex sfaclang.tex fformat.tex cfacdb.tex \
	   faq.tex gpl-3.0.tex

INCLISTINGS = $(TOP)/demo/structure/test.sf \
	      $(TOP)/demo/structure/ref/lev.asc \
	      $(TOP)/demo/structure/ref/tr.asc \
	      $(TOP)/faclib/schema.sql

REFMAN = cfacdbdoc/refman.tex

all : $(DOCS)

include $(TOP)/Make.conf

manual.pdf : manual.tex $(CHAPTERS) facref.bib $(INCLISTINGS) $(REFMAN)
	TEXINPUTS=".:./cfacdbdoc:" pdflatex manual
	bibtex manual
	TEXINPUTS=".:./cfacdbdoc:" makeindex manual
	TEXINPUTS=".:./cfacdbdoc:" pdflatex manual
	TEXINPUTS=".:./cfacdbdoc:" pdflatex manual

$(REFMAN) : $(TOP)/include/cfacdb.h cfacdb-doxy.cfg
	doxygen cfacdb-doxy.cfg

clean :
	rm -f *~ *.aux *.bbl *.blg *.idx *.log *.out *.toc *.ilg *.ind *.lol

devclean : clean
	rm -f $(DOCS)
	rm -rf cfacdbdoc

check : dummy

install : $(DOCS)
	$(INSTALLDIR) $(DOCDIR)
	$(INSTALLDATA) $(DOCS) $(DOCDIR)

dummy : 
