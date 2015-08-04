TOP = ..

DOCS = manual.pdf

all : $(DOCS)

include $(TOP)/Make.conf

manual.pdf : manual.tex facref.bib
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
	$(INSTALL) $(DOCS) $(DOCDIR)

dummy : 
