#!/bin/sh

rm *.aux *.bbl *.blg *.log *.out *.toc
pdflatex EMCluster-guide.Rnw
bibtex EMCluster-guide
pdflatex EMCluster-guide.Rnw
pdflatex EMCluster-guide.Rnw
pdflatex EMCluster-guide.Rnw
rm *.aux *.bbl *.blg *.log *.out *.toc
