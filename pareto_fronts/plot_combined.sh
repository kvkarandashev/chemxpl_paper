#!/bin/bash

fontsize=11
dl=2

tab_latex=$1
final_latex=$2

cat > $final_latex << EOF
\documentclass[preview]{standalone}%
\usepackage{amsmath,amsfonts,amsthm}
\usepackage{placeins}
\usepackage{tabularx,multirow}
\usepackage[T1]{fontenc}
\usepackage{graphicx}
\usepackage[english]{babel}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\newcommand{\gap}{\Delta \epsilon}
\usepackage{xcolor}
\textwidth=16.8cm
\begin{document}

\fontsize{$fontsize}{$dl}\selectfont
\begin{figure}[tbp]
\centering
EOF

sed "s/llll/cccc/g" $tab_latex | grep -v "rule"  >> $final_latex

cat >> $final_latex << EOF
\end{figure}

\end{document}
EOF
pdflatex $final_latex