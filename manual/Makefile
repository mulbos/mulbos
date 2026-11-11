# Makefile for generating presentation $(TARGET).pdf from $(TARGET).tex

TARGET = Manual

Subfiles = $(TARGET).tex

PLOTS = 

all: $(TARGET).pdf $(Subfiles) MulBos.bib

$(TARGET).pdf: $(TARGET).tex MulBos.bib $(Subfiles)
	pdflatex $(TARGET).tex
	bibtex $(TARGET)
	pdflatex $(TARGET).tex
	bibtex $(TARGET)
	pdflatex $(TARGET).tex
#	dvips $(TARGET).dvi -o
#	ps2pdf $(TARGET).ps $(TARGET).pdf
#	xpdf $(TARGET).pdf &

clean:
	rm -f *.dvi *.log *.aux *~ *.bbl *.blg *.out *.toc

