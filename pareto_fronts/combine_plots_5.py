import glob, subprocess
import pandas as pd

biases = ["none", "weak", "stronger"]

bias_names = {
    "none": "$\\alpha_{b}=0.0$",
    "weak": "$\\alpha_{b}=0.2$",
    "stronger": "$\\alpha_{b}=0.4$",
}

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

datasets = ["QM9", "EGP"]

# gap_constr_name = {"strong": "strong $\gap$ constr.", "weak": "weak $\gap$ constr."}
yaxis_hspace = 35.0
cut_bottom = 10.0


def rotated_label(label, left=False, bottom=False):
    if left:
        indent = "\hspace{-34ex}"
    else:
        hspace = yaxis_hspace
        if not bottom:
            hspace += cut_bottom
        indent = "\hspace{-" + str(hspace) + "ex}"
    output = (
        "\\begin{tabular}{@{}c@{}}\\rotatebox[origin=c]{90}{"
        + label
        + indent
        + "}\end{tabular}"
    )
    if left:
        output = "\multicolumn{1}{c|}{" + output + "}"
    return output


def pareto_plot_for(dataset, quantity, gap_constraint, bias):
    return glob.glob(
        "from_Jan/plots/"
        + dataset
        + "_"
        + quantity
        + "_"
        + gap_constraint
        + "_"
        + bias
        + "_*_den.png"
    )[0]


def yaxis_label(dataset):
    return (
        "$\Delta \epsilon^{\mathrm{cheap}}/\mathrm{max}(\Delta\epsilon^{\mathrm{"
        + dataset
        + "}})$"
    )


def xaxis_label(dataset, quantity):
    if quantity == "dipole":
        return "$D^{\mathrm{cheap}}/\mathrm{max}(D^{\mathrm{" + dataset + "}})$"
    else:
        return (
            "$\Delta G_{\mathrm{solv.}}^{\mathrm{cheap}}/\mathrm{max}(|\Delta G_{\mathrm{solv.}}^{\mathrm{"
            + dataset
            + "}}|)$"
        )


def included_pareto_plot(dataset, quantity, gap_constraint, bias, bottom=False):
    axisy = 0.0  # 120.0
    if bottom:
        triml = "0"
        #        xaxis_label_line="\\node[below right] (xaxis) {"+xaxis_label(dataset, quantity)+"};"
        xaxis_label_line = ""
    else:
        #        triml="{0.08\height}"
        triml = "{" + str(cut_bottom) + "ex}"
        axisy += cut_bottom  # /2 # divide by 2 to account for center shifting
        xaxis_label_line = ""
    inline = (
        "\\begin{tikzpicture}\
        \\node[anchor=base, inner sep=0] (image) at (0.0ex,0.0ex) {\\adjincludegraphics[width=0.375\\linewidth, trim={0 "
        + triml
        + " 0 0}]\
{"
        + pareto_plot_for(dataset, quantity, gap_constraint, bias)
        + "}};\
\\node[rotate=90, above left = 0.0ex and 0.0ex of image] (yaxis) {"
        + yaxis_label(dataset)
        + "};\
"
        + xaxis_label_line
        + "\end{tikzpicture}"
    )
    return inline


#%\coordinate (text) at ($ (image) + (0.0ex,"+str(axisy)+"ex) $);\


phantom_label = "\phantom{\_}"


def legend():
    return (
        "\includegraphics[width=0.375\linewidth]{from_Jan/legends/legend_vertical.png}"
    )


def pareto_plot_table(dataset, quantity, gap_constraint):
    xaxis_label_line = xaxis_label(dataset, quantity) + " \\vspace{1.0ex}"
    max_colwidth = 0
    bias_col = []
    plot_col = []

    for bias in biases:
        is_bottom = bias == "stronger"
        bias_string = rotated_label(bias_names[bias], left=True)
        plot_string = included_pareto_plot(
            dataset, quantity, gap_constraint, bias, bottom=is_bottom
        )
        bias_col.append(bias_string)
        plot_col.append(plot_string)
        for s in [bias_string, plot_string]:
            max_colwidth = max(max_colwidth, len(s))

    for _ in range(2):
        bias_col.append(phantom_label)
    plot_col.append("\hspace{15.0ex}" + xaxis_label_line)
    plot_col.append(legend())

    table = {}
    for l in [bias_col, plot_col]:
        table[l[0]] = l[1:]
    return table, max_colwidth


def combine_pareto_plots(dataset, quantity, gap_constraint, file_prefix):
    table, max_colwidth = pareto_plot_table(dataset, quantity, gap_constraint)
    df = pd.DataFrame(table)
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            file_prefix + "_tab.tex", index=False, escape=False, multicolumn=True
        )
    subprocess.run(["./plot_combined_vertical.sh", file_prefix + "_tab", file_prefix])


def main():
    for dataset in datasets:
        for quantity in quantities:
            for gap_constraint in gap_constraints:
                combine_pareto_plots(
                    dataset,
                    quantity,
                    gap_constraint,
                    dataset + "_" + quantity + "_" + gap_constraint + "_pareto_plots",
                )


if __name__ == "__main__":
    main()
