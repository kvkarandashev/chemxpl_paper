#!/bin/bash

picture_height=18.0

function create_for () {
    dataset=$1
    figure_dir="opt_log_figures"
    inserted_png="$figure_dir/opt_logs_"$dataset"_solvation_energy_weak.png"

    suffixes=("solvation_energy_strong" "dipole_weak" "dipole_strong")

    tex_name=opt_log_SI_$dataset.tex
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
EOF
    image_id=0
    image_y=0.0
    for suffix in ${suffixes[@]}
    do
        image_id=$((image_id+1))
        inserted_png="$figure_dir/opt_logs_"$dataset"_"$suffix".png"
        cat >> $tex_name << EOF
\node[anchor=base, inner sep=0] (image$image_id) at (0.0ex,${image_y}ex) {\includegraphics[width=0.25\linewidth]
{$inserted_png}};
\node[above left = -2.0ex and -2.0ex of image$image_id, rotate=90] (yaxis$imageid) {rel. improv.};
EOF
        image_y=$(echo $image_y $picture_height | awk '{print $1-$2}')
    done
    cat >> $tex_name << EOF
\node[below right = 0.0ex and -14.0ex of image$image_id] (xaxis) {$ N^{\mathrm{step}}$};
\node[below right = 2.5ex and -26.0ex of image$image_id] (legend) {\includegraphics[width=0.22\linewidth]{$figure_dir/trimmed_legend_2.png}};
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