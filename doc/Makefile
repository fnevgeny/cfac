doc: manual.pdf

manual.pdf: manual.tex
	pdflatex manual
	bibtex manual
	makeindex manual
	pdflatex manual
	pdflatex manual

clean:
	rm -f *~ *.aux *.bbl *.blg *.idx *.log *.out *.toc *.ilg *.ind

devclean: clean
	rm -f manual.pdf

check: dummy

dummy : 
