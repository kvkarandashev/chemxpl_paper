import glob, subprocess
import pandas as pd

biases = ["none", "weak", "stronger"]

bias_names = {"none": "no bias", "weak": "weak bias", "stronger": "stronger bias"}

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

datasets = ["QM9", "EGP"]

gap_constr_name = {"strong": "strong $\gap$ constr.", "weak": "weak $\gap$ constr."}


def rotated_label(label):
    return (
        "\\begin{tabular}{@{}c@{}}\\rotatebox[origin=c]{90}{"
        + label
        + "\hspace{-32ex}}\end{tabular}"
    )


def pareto_plot_for(dataset, quantity, gap_constraint, bias):
    return glob.glob(
        "from_Jan/"
        + dataset
        + "/pictures/"
        + dataset
        + "_"
        + quantity
        + "_"
        + gap_constraint
        + "_"
        + bias
        + "_*_den.png"
    )[0]


def included_pareto_plot(dataset, quantity, gap_constraint, bias):
    return (
        "\includegraphics[width=0.3\linewidth]{"
        + pareto_plot_for(dataset, quantity, gap_constraint, bias)
        + "}"
    )


phantom_label = "\phantom{\_}"


def pareto_plot_table(dataset, quantity):
    max_colwidth = 0
    table = {phantom_label: []}
    for gap_constraint in gap_constraints:
        s = rotated_label(gap_constr_name[gap_constraint])
        max_colwidth = max(max_colwidth, len(s))
        table[phantom_label].append(s)
    table = {
        "\phantom{\_}": [
            rotated_label(gap_constr_name[gap_constraint])
            for gap_constraint in gap_constraints
        ]
    }
    for bias in biases:
        h = bias_names[bias]
        table[h] = []
        for gap_constraint in gap_constraints:
            s = included_pareto_plot(dataset, quantity, gap_constraint, bias)
            max_colwidth = max(len(s), max_colwidth)
            table[h].append(s)
    return table, max_colwidth


def combine_pareto_plots(dataset, quantity, file_prefix):
    table, max_colwidth = pareto_plot_table(dataset, quantity)
    df = pd.DataFrame(table)
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            file_prefix + "_tab.tex", index=False, escape=False, multicolumn=True
        )
    subprocess.run(
        ["./plot_combined.sh", file_prefix + "_tab.tex", file_prefix + ".tex"]
    )


def main():
    for dataset in datasets:
        for quantity in quantities:
            combine_pareto_plots(
                dataset, quantity, dataset + "_" + quantity + "_pareto_plots"
            )


if __name__ == "__main__":
    main()
