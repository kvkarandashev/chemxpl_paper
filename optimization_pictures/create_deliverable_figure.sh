#!/bin/bash

python plot_all_optimization_logs_3.py deliverable

png_x_step=24.0
png_y_step=18.0

column_separator=4.0

figure_dir="opt_log_figures"
best_cand_dir="best_candidates"

tex_name=deliverable_opt_log_best_cand.tex

# The best candidate pictures are created separately.
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
EOF

cur_y=0.0
image_id=0

for quantity in solvation_energy dipole
do
    cur_x=0.0
    for png_type in opt_log candidate
    do
        for dataset in QM9 EGP
        do
            if [ "$quantity" == "solvation_energy" ]
            then
                axis_name="solv. en."
            else
                axis_name="\phantom{00}dipole\phantom{0}"
            fi

            if [ $png_type == opt_log ]
            then
                inserted_png="$figure_dir/opt_logs_"$dataset"_"$quantity"_weak.png"
            else
                inserted_png=$best_cand_dir"/"$dataset"_"$quantity"_weak.png"
            fi
            image_name=image$image_id
cat >>$tex_name << EOF
    \node[anchor=base, inner sep=0] ($image_name) at (${cur_x}ex,${cur_y}ex) {\includegraphics[width=0.225\linewidth]{$inserted_png}};
EOF
            image_id=$((image_id+1))
            if [ "$png_type" == "opt_log" ]
            then
                if [ "$dataset" == "QM9" ]
                then
                    cat >> $tex_name << EOF
    \node[above left = 0.0ex and -2.0ex of $image_name, rotate=90] (yaxis$image_name) {rel. improv.};
    \node[above left = -1.0ex and 1.0ex of $image_name, rotate=90] (yaxis$image_name) {{\\bf $axis_name}};
EOF
                fi
                if [ "$quantity" == "dipole" ]
                then
                    cat >> $tex_name << EOF
    \node[below right = 0.0ex and -14.0ex of $image_name] (xaxis$image_name) {$ N_{\mathrm{iter}}$};
EOF
                fi
            fi
            if [ "$quantity" == "solvation_energy" ]
            then
                title_x=$(echo $cur_x | awk '{print $1+1.0}')
                cat >> $tex_name << EOF
%    \node[above right = 0.0ex and -14.0ex of $image_name] (title$image_name) {{\\bf ${dataset}*}};
    \node[anchor=base, inner sep=0] (title$image_name) at (${title_x}ex,2.5) {{\\bf ${dataset}*}};
EOF

            fi
            cur_x=$(echo $cur_x $png_x_step | awk '{print $1+$2}')
        done
        cur_x=$(echo $cur_x $column_separator | awk '{print $1+$2}')
    done
    cur_y=$(echo $cur_y $png_y_step | awk '{print $1-$2}')
done
cat >> $tex_name << EOF
\end{tikzpicture}
\end{figure}

\end{document}

EOF

pdflatex $tex_name
prefix=$(echo $tex_name | cut -d'.' -f1)

pdftoppm -r 1000 $prefix.pdf $prefix -png
mv $prefix-1.png $prefix.png
convert -density 1000 $prefix.png -trim -quality 500 $prefix.png