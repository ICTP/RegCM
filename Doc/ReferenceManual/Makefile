

all: ReferenceMan.pdf

FIGURES = grid2.eps sigma_levels.eps
TEXFILES = abstract.tex acronyms.tex description.tex \
					 intro.tex ack.tex future.tex
BIBFILES = references.bib 

ReferenceMan.pdf: ReferenceMan.tex $(TEXFILES) $(BIBFILES) $(FIGURES)
	latex ReferenceMan
	bibtex ReferenceMan
	latex ReferenceMan
	latex ReferenceMan
	latex ReferenceMan
	dvipdf ReferenceMan

clean:
	rm -f *.aux *.dvi *.pdf *.log *.toc *.lof *.lot *.bbl *.blg
