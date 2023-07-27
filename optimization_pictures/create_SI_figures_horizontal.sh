#!/bin/bash

picture_width=27.0

rightmost_x_coord=$(echo $picture_width 3 | awk '{print $1*$2}')

function create_for () {
    dataset=$1
    figure_dir="opt_log_figures"
    inserted_png="$figure_dir/opt_logs_"$dataset"_solvation_energy_weak.png"

#    suffixes=("solvation_energy_strong" "dipole_weak" "dipole_strong")
    inv_suffixes=("dipole_strong" "dipole_weak" "solvation_energy_strong")
    inv_indices=("(c)" "(b)" "(a)")

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
\newcommand{\gap}{\Delta \epsilon^{\mathrm{cheap}}}
\usepackage{xcolor}
\usetikzlibrary{calc}
\usetikzlibrary{positioning}
\usepackage[export]{adjustbox}
\newcommand{\dEsolv}{\Delta G_{\mathrm{solv.}}^{\mathrm{cheap}}}
\newcommand{\dipole}{D^{\mathrm{cheap}}}
\textwidth=16.8cm
\fontsize{11}{2}\selectfont
\begin{document}

\begin{figure}[tbp]
\centering
\begin{tikzpicture}
EOF
    image_id=0
    image_x=$rightmost_x_coord
    for suffix in ${inv_suffixes[@]}
    do
        index=${inv_indices[$image_id]}
        image_id=$((image_id+1))
        inserted_png="$figure_dir/opt_logs_"$dataset"_"$suffix".png"
        cat >> $tex_name << EOF
\node[anchor=base, inner sep=0] (image$image_id) at (${image_x}ex,0.0ex) {\includegraphics[width=0.25\linewidth]
{$inserted_png}};
\node[above left = -1.0ex and -6.0ex of image$image_id] (index) {$index};
\node[below right = 0.0ex and -14.0ex of image$image_id] (xaxis) {$ N^{\mathrm{step}}$};
EOF
        if [ "$image_id" == "3" ]
        then
            echo "\node[above left = -2.0ex and -1.5ex of image$image_id, rotate=90] (yaxis$imageid) {rel. improv.};" >> $tex_name
        fi
        image_x=$(echo $image_x $picture_width | awk '{print $1-$2}')
    done
    cat >> $tex_name << EOF
\node[below right = 2.5ex and -26.0ex of image2] (legend) {\includegraphics[width=0.22\linewidth]{$figure_dir/trimmed_legend_2.png}};
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