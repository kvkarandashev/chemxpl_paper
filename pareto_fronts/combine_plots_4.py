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

gap_constr_name = {"strong": "strong $\gap$ constr.", "weak": "weak $\gap$ constr."}


def rotated_label(label, right=False):
    output = (
        "\\begin{tabular}{@{}c@{}}\\rotatebox[origin=c]{90}{"
        + label
        + "\hspace{-27ex}}\end{tabular}"
    )
    if right:
        output = "\multicolumn{1}{c|}{" + output + "}"
    #    else:
    #        output+=" \hspace{-12ex}"
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


def included_pareto_plot(dataset, quantity, gap_constraint, bias, left=False):
    inline = (
        "\hspace{-2.0ex}\includegraphics[width=0.275\linewidth]{"
        + pareto_plot_for(dataset, quantity, gap_constraint, bias)
        + "}"
    )
    #    inline="\raisebox{-0.04ex}{"+inline+"}"
    #    inline="\begin{raisedtext}{-0.05ex}"+inline+"\end{raisedtext}"
    #    if left:
    #        cells="|c"
    #    else:
    #        cells="c"
    #    tab="\multicolumn{1}{"+cells+"}{" + inline + "}"
    #    return tab
    return inline


phantom_label = "\phantom{\_}"


def legend_multcol():
    return "\multicolumn{3}{c}{\includegraphics[width=0.6\linewidth]{from_Jan/legends/legend.pdf}}"


def pareto_plot_table(dataset, quantity):
    if quantity == "dipole":
        axis_label = "$D/\mathrm{max}(D^{\mathrm{" + dataset + "}})$"
    else:
        axis_label = (
            "$\Delta G_{\mathrm{solv.}}/\mathrm{max}(|\Delta G_{\mathrm{solv.}}^{\mathrm{"
            + dataset
            + "}}|)$"
        )
    axis_label += " \\vspace{1.5ex}"
    max_colwidth = 0
    lcol = []
    lcol1 = []
    for gap_constraint in gap_constraints:
        lcols = rotated_label(gap_constr_name[gap_constraint], right=True)
        lcols1 = rotated_label(
            "$\Delta \epsilon/\mathrm{max}(\Delta\epsilon^{\mathrm{" + dataset + "}})$"
        )
        lcol.append(lcols)
        lcol1.append(lcols1)
        max_colwidth = max(max_colwidth, len(lcols))
        max_colwidth = max(max_colwidth, len(lcols1))
    for _ in range(2):
        lcol.append(phantom_label)
        lcol1.append(phantom_label)
    #    lcol[0]="\cline{3-5}"+lcol[0]
    table = {phantom_label: lcol, "\phantom{1}": lcol1}
    for i, bias in enumerate(biases):
        h = "$\\vphantom{\prod}$" + bias_names[bias]
        table[h] = []
        for gap_constraint in gap_constraints:
            s = included_pareto_plot(
                dataset, quantity, gap_constraint, bias
            )  # , left=(gap_constraint=="weak"))
            max_colwidth = max(len(s), max_colwidth)
            table[h].append(s)
        table[h].append(axis_label)
        if i == 0:
            table[h].append(legend_multcol())
        else:
            table[h].append("delete")
    return table, max_colwidth


def combine_pareto_plots(dataset, quantity, file_prefix):
    table, max_colwidth = pareto_plot_table(dataset, quantity)
    df = pd.DataFrame(table)
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            file_prefix + "_tab.tex", index=False, escape=False, multicolumn=True
        )
    subprocess.run(["./plot_combined.sh", file_prefix + "_tab", file_prefix])


def main():
    for dataset in datasets:
        for quantity in quantities:
            combine_pareto_plots(
                dataset, quantity, dataset + "_" + quantity + "_pareto_plots"
            )


if __name__ == "__main__":
    main()
