import matplotlib.pyplot as plt
from glob import glob
from subprocess import run
import numpy as np

bias_values = ["none", "weak", "stronger"]

bias_linecolor = {"none": "blue", "weak": "green", "stronger": "red"}

solv_en = "solvation_energy"
dipole = "dipole"
QM9 = "QM9"
EGP = "EGP"

quantity_names = [solv_en, dipole]
datasets = [QM9, EGP]

gap_constraints = ["weak", "strong"]


data_dir = "./running_minima"

output_dir = "./opt_log_figures"

plot_averages_werrors = True

linewidth = 2.0
err_linewidth = 1.5
err_capsize = 6.0
markersize = 16

minor_tick_length_coeff = 2.0
major_tick_length_coeff = 3.5
markeredge_coeff = 0.75

title_fontsize = 40.0

legend_fontsize = 40.0

ticks_fontsize = 32.0

leftover_mult_x = 2.0
leftover_mult_y = 0.14

bias_linestyle = {"none": "dotted", "weak": "dashdot", "stronger": "dashed"}

negligible_step_float = 1.0e-5
plotted_averages_number = 16

xtick_labels = {
    0.2: r"$2\times10^{-1}$",
    0.4: r"$4\times10^{-1}$",
    0.3: r"$3\times10^{-1}$",
    0.6: r"$6\times10^{-1}$",
    0.8: r"$8\times10^{-1}$",
    0.05: r"$5\times10^{-2}$",
}

for i in range(-1, 8):
    xtick_labels[10.0**i] = r"$10^{" + str(i) + "}$"


def latex_scientific(number_in, nfigures=0):
    def_sci = ("{:0." + str(nfigures) + "e}").format(number_in)
    parts = def_sci.split("e")
    output = parts[0]
    if len(parts) > 1:
        extra = parts[1]
        if extra != "+00":
            output += "\cdot 10^{" + str(int(extra)) + "}"
    return output


class latex_format:
    def __init__(self, scientific=True, nfigures=0):
        self.scientific = scientific
        self.nfigures = nfigures

    def __call__(self, number_in):
        if not self.scientific:
            return str(int(number_in))
        if number_in == 0.0:
            return "0.0"
        output = latex_scientific(number_in, nfigures=self.nfigures)
        output = r"$" + output + "$"
        return output


ytick_positions = {
    QM9: {
        solv_en: [1e-2, 0.0, -1e-2, -2e-2, -3e-2, -4e-2],
        dipole: [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    },
    EGP: {
        solv_en: [-0.0, -1.0e-1, -2.0e-1, -3.0e-1, -4.0e-1, -5.0e-1],
        dipole: list(range(0, 50, 10)),
    },
}

ytick_label_format = {
    QM9: {
        solv_en: {"scientific": True, "nfigures": 0},
        dipole: {"scientific": False, "nfigures": 0},
    },
    EGP: {
        solv_en: {"scientific": True, "nfigures": 0},
        dipole: {"scientific": True, "nfigures": 0},
    },
}

# xtick_positions=[1.0, 1.e+1, 1.e+2, 1.e+3, 1.e+4, 1.e+5]
xtick_positions = [1.0, 1.0e2, 1.0e4]

xtick_labels = [xtick_labels[xtick_pos] for xtick_pos in xtick_positions]

left_indent = 0.25
bottom_indent = 0.1

right_indent = 0.05
top_indent = 0.01

position = [
    left_indent,
    bottom_indent,
    1.0 - right_indent - left_indent,
    1.0 - top_indent - bottom_indent,
]


def line2step_quant_SMILES(line):
    lsplit = line.split()
    return int(lsplit[0]), float(lsplit[1]), lsplit[2]


def log_spaced_values(lines):
    step_start = 1.0
    step_finish, _, _ = line2step_quant_SMILES(lines[-1])
    log2start = np.log2(step_start)
    log2finish = np.log2(float(step_finish))

    points_x = []
    points_y = []

    unrounded_step_vals = np.logspace(
        log2start, log2finish, num=plotted_averages_number, base=2
    )
    for unrounded_step_val in unrounded_step_vals:
        step_val = int(unrounded_step_val + negligible_step_float)
        if (len(points_x) != 0) and (points_x[-1] == step_val):
            continue
        cur_step_val, cur_val, _ = line2step_quant_SMILES(lines[step_val])
        assert cur_step_val == step_val
        points_x.append(step_val)
        points_y.append(cur_val)
    return points_x, points_y


def optimization_log_points(filename):
    input = open(filename, "r")
    lines = input.readlines()
    input.close()

    if plot_averages_werrors:
        return log_spaced_values(lines)

    step0, cur_quant, cur_SMILES = line2step_quant_SMILES(lines[0])
    assert step0 == 0

    cur_start_step = 1
    points_x = [cur_start_step]
    points_y = [cur_quant]

    for l in lines:
        step, quant, SMILES = line2step_quant_SMILES(l)
        if SMILES == cur_SMILES:
            continue
        if step != cur_start_step:
            points_x.append(step)
            points_y.append(cur_quant)
        points_x.append(step)
        points_y.append(quant)
        cur_quant = quant
        cur_start_step = step
    last_step, last_quant, _ = line2step_quant_SMILES(lines[-1])
    if last_step != cur_start_step:
        points_x.append(last_step)
        points_y.append(last_quant)
    return points_x, points_y


def comp_wNone(val, prev_ext, ext_func):
    if prev_ext is None:
        return val
    else:
        return ext_func(val, prev_ext)


def minmax_vals(val, prev_min, prev_max):
    return comp_wNone(val, prev_min, min), comp_wNone(val, prev_max, max)


def list_minmax_vals(val_list, cur_list_min=None, cur_list_max=None):
    list_min = cur_list_min
    list_max = cur_list_max
    for val in val_list:
        list_min, list_max = minmax_vals(val, list_min, list_max)
    return list_min, list_max


def list_minmax_vals_werrs(mean_list, std_list, cur_list_min=None, cur_list_max=None):
    list_min = cur_list_min
    list_max = cur_list_max
    for mean, err in zip(mean_list, std_list):
        list_min = comp_wNone(mean - err, list_min, min)
        list_max = comp_wNone(mean + err, list_max, max)
    return list_min, list_max


def plot_ax_bias(ax, points_x, points_y, bias, yerr=None):
    color = bias_linecolor[bias]
    linestyle = bias_linestyle[bias]
    plot_kwargs = {"color": color, "linestyle": linestyle, "linewidth": linewidth}
    ax.plot(points_x, points_y, **plot_kwargs)
    if yerr is None:
        return
    ax.errorbar(
        points_x,
        points_y,
        yerr=yerr,
        **plot_kwargs,
        elinewidth=err_linewidth,
        capsize=err_capsize
    )


def plot_opt_log_filenames(
    input_filenames, bias_vals, ytick_positions, ytick_label_format, output_filename
):
    fig = plt.figure(constrained_layout=True)
    ax = fig.add_subplot()
    plt.setp(ax.spines.values(), linewidth=linewidth)

    if plot_averages_werrors:
        all_points_xy_dict = {}
        for bias in bias_values:
            all_points_xy_dict[bias] = []

    min_x = None
    max_x = None
    min_y = None
    max_y = None

    for input_filename, bias in zip(input_filenames, bias_vals):
        points_x, points_y = optimization_log_points(input_filename)
        min_x, max_x = list_minmax_vals(
            points_x, cur_list_min=min_x, cur_list_max=max_x
        )
        if plot_averages_werrors:
            all_points_xy_dict[bias].append(points_y)
        else:
            plot_ax_bias(ax, points_x, points_y, bias)
            min_y, max_y = list_minmax_vals(
                points_y, cur_list_min=min_y, cur_list_max=max_y
            )

    if plot_averages_werrors:
        for bias in bias_vals:
            y_means = []
            y_stds = []
            for y_tuple in zip(*all_points_xy_dict[bias]):
                y_array = np.array(y_tuple)
                y_means.append(np.mean(y_array))
                y_stds.append(np.std(y_array))
            min_y, max_y = list_minmax_vals_werrs(
                y_means, y_stds, cur_list_min=min_y, cur_list_max=max_y
            )
            plot_ax_bias(ax, points_x, y_means, bias, yerr=y_stds)
    #                ax.plot(x_vals, y_vals, label=legend, color=color, marker=marker, markersize=markersize,
    #                            markerfacecolor=facecolors, linestyle=linestyle, markeredgewidth=markeredge_coeff*linewidth,
    #                            linewidth=linewidth)

    final_min_x = min_x / leftover_mult_x
    final_max_x = max_x * leftover_mult_x

    y_spread = max_y - min_y
    extra_y_spread = y_spread * leftover_mult_y

    final_min_y = min_y - extra_y_spread
    final_max_y = max_y + extra_y_spread

    ax.set_xlim(final_min_x, final_max_x)
    ax.set_ylim(final_min_y, final_max_y)

    #    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.tick_params(
        which="minor",
        width=linewidth,
        length=minor_tick_length_coeff * linewidth,
        labelsize=0,
    )
    ax.tick_params(
        which="major",
        width=linewidth,
        length=major_tick_length_coeff * linewidth,
        labelsize=0,
    )

    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=ticks_fontsize)

    format_func = latex_format(**ytick_label_format)

    cur_ytick_positions = []
    cur_ytick_labels = []
    for pos in ytick_positions:
        if (pos < final_min_y) or (pos > final_max_y):
            continue
        cur_ytick_positions.append(pos)
        cur_ytick_labels.append(format_func(pos))

    ax.set_yticks(cur_ytick_positions)
    ax.set_yticklabels(cur_ytick_labels, fontsize=ticks_fontsize)

    fig.set_figwidth(8.0)
    fig.set_figheight(6.0)

    ax.set_position(position)

    fig.savefig(output_filename)


def plot_opt_log_diff_bias(
    data_dir, dataset, quantity_name, gap_constraint, figure_dir
):
    batch_name = dataset + "_" + quantity_name + "_" + gap_constraint
    all_input_filenames = []
    all_biases = []
    for bias in bias_values:
        input_filename_family_name = (
            data_dir + "/running_min_" + batch_name + "_" + bias
        )
        input_filenames = glob(input_filename_family_name + "_*.txt")
        all_input_filenames += input_filenames
        all_biases += [bias for _ in input_filenames]
    figure_file = figure_dir + "/opt_logs_" + batch_name + ".png"
    plot_opt_log_filenames(
        all_input_filenames,
        all_biases,
        ytick_positions[dataset][quantity_name],
        ytick_label_format[dataset][quantity_name],
        figure_file,
    )


def main():
    run(["mkdir", "-p", output_dir])
    for dataset in datasets:
        for quantity_name in quantity_names:
            for gap_constraint in gap_constraints:
                plot_opt_log_diff_bias(
                    data_dir, dataset, quantity_name, gap_constraint, output_dir
                )


if __name__ == "__main__":
    main()
