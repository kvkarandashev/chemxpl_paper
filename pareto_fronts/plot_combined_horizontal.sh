#!/bin/bash

fontsize=11
dl=2

tab_latex=$1
final_latex=$2

final_latex_full=$final_latex.tex

cat > $final_latex_full << EOF
\documentclass[preview]{standalone}%
\usepackage{amsmath,amsfonts,amsthm}
\usepackage{placeins}
\usepackage{tabularx,multirow,tikz}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage[english]{babel}
\usepackage{amsmath}
\usepackage{booktabs}
\usepackage{amsfonts}
\usepackage{amssymb}
\newcommand{\gap}{\Delta \epsilon}
\usepackage{xcolor}
\usetikzlibrary{calc}
\usetikzlibrary{positioning}
\usepackage[export]{adjustbox}
\textwidth=16.8cm
\begin{document}

\fontsize{$fontsize}{$dl}\selectfont
\begin{figure}[tbp]
\centering
EOF

sed -i 's/ \{1,\}/ /g' $tab_latex.tex
sed -i 's/& delete/ /g' $tab_latex.tex

sed "s/lll/ccc/g" $tab_latex.tex | grep -v "rule"  >> $final_latex_full
#grep -v "rule" $tab_latex.tex >> $final_latex_full

cat >> $final_latex_full << EOF
\end{figure}

\end{document}
EOF
pdflatex $final_latex_full

pdftoppm $final_latex.pdf $final_latex -png

mv $final_latex-1.png $final_latex.png

#png_cut $final_latex.png cropped_$final_latex.png "0.03:0.03:0:0.02"
#-sharpen 0x1.0
convert -density 250 $final_latex.png -trim -quality 100 $final_latex.png
