# SPDX-FileCopyrightText: 2025 Ovidio Garcia-Oliva
# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de

mkdir -p out
mkdir -p fig

Rscript isopleth_plots.R
Rscript physico_chemical_analysis.R
pdflatex Fig_3.tex
rm Fig_3.aux Fig_3.log Fig_3.out
mv Fig_3.pdf fig/Fig_3.pdf
