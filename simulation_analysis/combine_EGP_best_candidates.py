import glob, subprocess
import pandas as pd


class PhantomId:
    def __init__(self):
        self.counter = 0

    def __call__(self):
        self.counter += 1
        return "\phantom{" + str(self.counter) + "}"


biases = ["none", "weak", "stronger"]

bias_names = {"none": "no bias", "weak": "weak bias", "stronger": "stronger bias"}

gap_constraints = ["weak", "strong"]

quantities = ["solvation_energy", "dipole"]

datasets = ["QM9", "EGP"]

gap_constr_name = {"strong": "strong $\gap$ constr.", "weak": "weak $\gap$ constr."}

opt_prob_name = {"solvation_energy": "min. $\dEsolv$", "dipole": "max. $\dipole$"}


def rotated_label(label):
    return (
        "\\begin{tabular}{@{}c@{}}\\rotatebox[origin=c]{90}{"
        + label
        + "\hspace{-32ex}}\end{tabular}"
    )


def included_best_plot(ref_folder, quantity, gap_constraint, bias):
    return (
        "\includegraphics[width=0.20\linewidth]{"
        + ref_folder
        + "/best_"
        + bias
        + "_"
        + gap_constraint
        + "_"
        + quantity
        + "_pngs/best_"
        + bias
        + "_"
        + gap_constraint
        + "_"
        + quantity
        + "_rot_0.png"
        + "}"
    )


phantom = "\phantom{\_}"


def best_candidate_table(ref_folder):
    phl = PhantomId()
    phl_quant = phl()
    max_colwidth = 0
    table = {
        (phl_quant, "\cline{2-3}\cline{5-6}"): [
            rotated_label(bias_names[bias]) for bias in biases
        ]
    }
    for quantity in quantities:
        if quantity == quantities[-1]:
            table[(phl(), phantom)] = [phantom for _ in biases]
        for gap_constraint in gap_constraints:
            new_row = []
            for bias in biases:
                new_str = included_best_plot(ref_folder, quantity, gap_constraint, bias)
                max_colwidth = max(max_colwidth, len(new_str))
                new_row.append(new_str)
            table[(opt_prob_name[quantity], gap_constr_name[gap_constraint])] = new_row
    return table, max_colwidth


def main():
    ref_folder = "summary_xTB_dipsolv_opt_egp_cheap_3"
    output_prefix = "EGP_best_candidates"
    table, max_colwidth = best_candidate_table(ref_folder)
    df = pd.DataFrame(table)
    with pd.option_context("max_colwidth", max_colwidth):
        df.to_latex(
            output_prefix + "_tab.tex", index=False, escape=False, multicolumn=True
        )
    subprocess.run(
        ["./plot_combined.sh", output_prefix + "_tab.tex", output_prefix + ".tex"]
    )


if __name__ == "__main__":
    main()
