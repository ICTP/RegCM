

all: UserGuide.pdf

FIGURES = grid2.eps sigma_levels.eps
TEXFILES = abstract.tex acronyms.tex description.tex intro.tex ack.tex
BIBFILES = references.bib 

UserGuide.pdf: UserGuide.tex $(TEXFILES) $(BIBFILES) $(FIGURES)
	latex UserGuide
	bibtex UserGuide
	latex UserGuide
	latex UserGuide
	dvipdf UserGuide

clean:
	rm -f *.aux *.dvi *.pdf *.log *.toc *.lof *.lot *.bbl *.blg
