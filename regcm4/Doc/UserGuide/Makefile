

all: UserGuide.pdf

FIGURES = grid2.eps sigma_levels.eps xref1.fixed.eps xref2.fixed.eps \
          xref3.fixed.eps xref4.fixed.eps
TEXFILES = abstract.tex acronyms.tex description.tex intro.tex practice.tex  \
           ack.tex appendix.tex  preproc.tex
BIBFILES = refs.sara.bib 

UserGuide.pdf: UserGuide.tex $(TEXFILES) $(BIBFILES) $(FIGURES)
	latex UserGuide.tex
	latex UserGuide.tex
	latex UserGuide.tex
	dvipdf UserGuide.dvi

clean:
	rm -f *.aux *.dvi *.pdf *.log *.toc *.lof *.lot
