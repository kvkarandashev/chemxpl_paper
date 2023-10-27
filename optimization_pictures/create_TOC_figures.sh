#!/bin/bash

function create_for () {
    dataset=$1
    figure_dir="opt_log_figures_TOC"
    inserted_image_prefix="$figure_dir/opt_logs_"$dataset"_solvation_energy_weak"
    inserted_image=$inserted_image_prefix.pdf
#    pdfcrop $inserted_image_prefix.pdf
#    inserted_image=$inserted_image_prefix-crop.pdf

    tex_name=TOC_opt_log_solv_en_weak_constr_$dataset.tex
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
\newcommand{\dEsolv}{\Delta G_{\mathrm{solv}}}
\newcommand{\gap}{\Delta \epsilon^{\mathrm{cheap}}}
\textwidth=16.8cm
\fontsize{11}{2}\selectfont
\begin{document}

\begin{figure}[tbp]
\centering
\begin{tikzpicture}
\node[anchor=base, inner sep=0] (image) at (0.0ex,0.0ex) {\includegraphics[width=0.3\linewidth]{$inserted_image}};
\node[below right = -0.5ex and -22.0ex of image] (xaxis) {num. MC steps};
\node[above left = -0.1ex and -2.0ex of image, rotate=90] (yaxis) {improv. over $dataset}; 
\node[above left = -0.75ex and -31.5ex of image] (title) {$\dEsolv$ optimization};
%\node[below right = -14.5ex and -0.75ex of image] (legend) {\includegraphics[width=0.1\linewidth]{$figure_dir/trimmed_legend_1.png}};
\end{tikzpicture}
\end{figure}

\end{document}

EOF
    pdflatex $tex_name
    prefix=$(echo $tex_name | cut -d'.' -f1)

    # Trim pdf file.
    pdf_version=$prefix.pdf
    pdfcrop $pdf_version
    cropped_pdf=$prefix-crop.pdf

    # Create png file.
    png_version=$prefix.png
    pdftoppm -r 1000 $pdf_version $prefix -png
    mv $prefix-1.png $png_version
    convert -density 1000 $png_version -trim -quality 500 $png_version

    # Create eps file.
    eps_version=$prefix.eps
#    inkscape $prefix-crop.pdf --export-filename=$eps_version
    gs -q -dNOCACHE -dNOPAUSE -dBATCH -dSAFER -sDEVICE=eps2write -sOutputFile=$eps_version $cropped_pdf
#    convert $cropped_pdf $eps_version
    svg_version=$prefix.svg
    pdf2svg $cropped_pdf $svg_version
}

for dataset in QM9 EGP
do
    create_for $dataset
done