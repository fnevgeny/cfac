TOP = ..

DOCS = manual.pdf

CHAPTERS = overview.tex physics.tex sfaclang.tex fformat.tex faq.tex gpl-3.0.tex

INCLISTINGS = $(TOP)/demo/structure/test.sf \
	      $(TOP)/demo/structure/ref/lev.asc \
	      $(TOP)/demo/structure/ref/tr.asc \
	      $(TOP)/faclib/schema.sql

all : $(DOCS)

include $(TOP)/Make.conf

manual.pdf : manual.tex $(CHAPTERS) facref.bib $(INCLISTINGS)
	pdflatex manual
	bibtex manual
	makeindex manual
	pdflatex manual
	pdflatex manual

clean :
	rm -f *~ *.aux *.bbl *.blg *.idx *.log *.out *.toc *.ilg *.ind *.lol

devclean : clean
	rm -f $(DOCS)

check : dummy

install : $(DOCS)
	$(INSTALLDIR) $(DOCDIR)
	$(INSTALLDATA) $(DOCS) $(DOCDIR)

dummy : 
