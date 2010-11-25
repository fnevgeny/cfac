doc: manual.pdf

manual.pdf: manual.tex
	pdflatex manual
	bibtex manual
	makeindex manual
	pdflatex manual
	pdflatex manual

clean:
	rm -rf *~ *.aux *.bbl *.blg *.idx *.log *.out *.pdf *.toc *.ilg *.ind

check: dummy

dummy : 
