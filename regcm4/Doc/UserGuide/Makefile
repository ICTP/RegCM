

all: UserGuide.pdf

FIGURES = grid2.eps sigma_levels.eps xref1.fixed.eps xref2.fixed.eps \
          xref3.fixed.eps xref4.fixed.eps
TEXFILES = abstract.tex acronyms.tex description.tex intro.tex practice.tex  \
           ack.tex appendix.tex  preproc.tex
BIBFILES = references.bib 

UserGuide.pdf: UserGuide.tex $(TEXFILES) $(BIBFILES) $(FIGURES)
	latex UserGuide
	bibtex UserGuide
	latex UserGuide
	latex UserGuide
	dvipdf UserGuide

clean:
	rm -f *.aux *.dvi *.pdf *.log *.toc *.lof *.lot *.bbl *.blg
