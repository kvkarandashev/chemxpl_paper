#!/bin/bash

function create_for () {
    dataset=$1
    figure_dir="opt_log_figures"
    inserted_png="$figure_dir/opt_logs_"$dataset"_solvation_energy_weak.png"

    tex_name=opt_log_solv_en_weak_constr_$dataset.tes
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
\newcommand{\gap}{\Delta \epsilon}
\usepackage{xcolor}
\usetikzlibrary{calc}
\usetikzlibrary{positioning}
\usepackage[export]{adjustbox}
\newcommand{\dEsolv}{\Delta G_{\mathrm{solv.}}}
\textwidth=16.8cm
\fontsize{11}{2}\selectfont
\begin{document}

\begin{figure}[tbp]
\centering
\begin{tikzpicture}
\node[anchor=base, inner sep=0] (image) at (0.0ex,0.0ex) {\includegraphics[width=0.25\linewidth]
{$inserted_png}};
\node[below right = 0.0ex and -14.0ex of image] (xaxis) {$ N^{\mathrm{step}}$};
\node[above left = -2.0ex and -2.0ex of image, rotate=90] (yaxis) {rel. improv.};
\node[below right = -14.5ex and -0.75ex of image] (legend) {\includegraphics[width=0.1\linewidth]{$figure_dir/trimmed_legend_1.png}};
\end{tikzpicture}
\end{figure}

\end{document}

EOF
    pdflatex $tex_name
    prefix=$(echo $tex_name | cut -d'.' -f1)

    pdftoppm -r 1000 $prefix.pdf $prefix -png
    mv $prefix-1.png $prefix.png
    convert -density 1000 $prefix.png -trim -quality 500 $prefix.png
}

for dataset in QM9 EGP
do
    create_for $dataset
done