.PHONY: all
all:
	pdflatex nonlinear.tex
	biber nonlinear
	pdflatex nonlinear.tex
	pdflatex nonlinear.tex

.PHONY: clean
clean:
	rm -f *.aux
	rm -f *.bbl
	rm -f *.blg
	rm -f *.bcf
	rm -f *.ilg
	rm -f *.lof
	rm -f *.log
	rm -f *.lot
	rm -f *.nav
	rm -f *.out
	rm -f *.run.xml
	rm -f *.synctex*
	rm -f *.toc
