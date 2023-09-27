#!/bin/bash

inserted_png="opt_log_init_cond_comp.png"

tex_name=opt_log_init_cond_comp.tex
cat > $tex_name << EOF
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
\usepackage{xcolor}
\usetikzlibrary{calc}
\usetikzlibrary{positioning}
\usepackage[export]{adjustbox}
\newcommand{\dEsolv}{\Delta G^{\mathrm{cheap}}_{\mathrm{solv.}}}
\newcommand{\gap}{\Delta \epsilon^{\mathrm{cheap}}}
\textwidth=16.8cm
\fontsize{11}{2}\selectfont
\begin{document}

\begin{figure}[tbp]
\centering
\begin{tikzpicture}
\node[anchor=base, inner sep=0] (image) at (0.0ex,0.0ex) {\includegraphics[width=0.25\linewidth]{init_cond_comp_opt_log.png}};
\node[below right = 0.0ex and -14.0ex of image] (xaxis) {$ N^{\mathrm{step}}$};
\node[above left = -2.0ex and -2.0ex of image, rotate=90] (yaxis) {rel. improv.};
%\node[above left = 0.0ex and -28.0ex of image] (title) {$\dEsolv$/weak $\gap$ c.};
\node[below right = -14.5ex and -0.75ex of image] (legend) {\includegraphics[width=0.15\linewidth]{trimmed_init_cond_comp_legend_1.png}};
\end{tikzpicture}
\end{figure}

\end{document}

EOF
pdflatex $tex_name
prefix=$(echo $tex_name | cut -d'.' -f1)

pdftoppm -r 1000 $prefix.pdf $prefix -png
mv $prefix-1.png $prefix.png
convert -density 1000 $prefix.png -trim -quality 500 $prefix.png