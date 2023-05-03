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
cut_left = 10.0


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
        return (
            "\phantom{000}$D^{\mathrm{cheap}}/\mathrm{max}(D^{\mathrm{"
            + dataset
            + "}})$"
        )
    else:
        return (
            "$\Delta G_{\mathrm{solv.}}^{\mathrm{cheap}}/\mathrm{max}(|\Delta G_{\mathrm{solv.}}^{\mathrm{"
            + dataset
            + "}}|)$"
        )


def included_pareto_plot(dataset, quantity, gap_constraint, bias, left=False):
    if left:
        triml = "0"
        yaxis_label_line = (
            "\\node[above left = 0.0ex and 0.0ex of image, rotate=90] (yaxis) {"
            + yaxis_label(dataset)
            + "}; "
        )
    else:
        #        triml="{0.08\height}"
        triml = "{" + str(cut_left) + "ex}"
        yaxis_label_line = ""
    inline = (
        "\\begin{tikzpicture}\
        \\node[anchor=base, inner sep=0] (image) at (0.0ex,0.0ex) {\\includegraphics[scale=0.25, trim={"
        + triml
        + " 0 0 0}, clip]\
{"
        + pareto_plot_for(dataset, quantity, gap_constraint, bias)
        + "}}; \
\\node[below right = 0.0ex and -26.0ex of image] (xaxis) {"
        + xaxis_label(dataset, quantity)
        + "}; \
\\node[above right = 0.0ex and -20.0 ex of image] (bias) {"
        + bias_names[bias]
        + "}; \
"
        + yaxis_label_line
        + "\end{tikzpicture}"
    )
    return inline


#%\coordinate (text) at ($ (image) + (0.0ex,"+str(axisy)+"ex) $);\


phantom_label = "\phantom{\_}"


def legend():
    return "\includegraphics[width=0.6\linewidth]{from_Jan/legends/legend.pdf}"


def pareto_plot_table(dataset, quantity, gap_constraint):
    table = {}
    max_colwidth = 0
    for bias in biases:
        is_left = bias == "none"
        plot_string = included_pareto_plot(
            dataset, quantity, gap_constraint, bias, left=is_left
        )
        if is_left:
            col = "\multicolumn{3}{c}{" + legend() + "}"
        else:
            col = "delete"
        for s in [col, plot_string]:
            max_colwidth = max(max_colwidth, len(s))
        table[plot_string] = [col]

    return table, max_colwidth


def combine_pareto_plots(dataset, quantity, gap_constraint, file_prefix):
    table, max_colwidth = pareto_plot_table(dataset, quantity, gap_constraint)
    df = pd.DataFrame(table)
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            file_prefix + "_tab.tex", index=False, escape=False, multicolumn=True
        )
    subprocess.run(["./plot_combined_horizontal.sh", file_prefix + "_tab", file_prefix])


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
