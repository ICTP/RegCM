

all: DeveloperGuide.pdf

TEXFILES = MainDeveloperGuide.tex \
           DeveloperGuide.tex    TitleDeveloperGuide.tex
BIBFILLES = Biblio_Developer.bib

DeveloperGuide.pdf: DeveloperGuide.tex $(TEXFILES) $(BIBFILES) $(FIGURES)
	latex DeveloperGuide
	bibtex DeveloperGuide
	latex DeveloperGuide
	latex DeveloperGuide
	dvipdf DeveloperGuide

clean:
	  rm -f *.aux *.dvi *.pdf *.log *.toc *.lof *.lot *.bbl *.blg
