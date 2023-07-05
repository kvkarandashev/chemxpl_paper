import matplotlib.pyplot as plt
from glob import glob
from subprocess import run
import numpy as np
from matplotlib.ticker import MultipleLocator

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
leftover_mult_y = 0.05

fig_width = 8.0
fig_height = 5.5

bias_linestyle = {"none": "dotted", "weak": "dashdot", "stronger": "dashed"}

init_val_linestyle = "dashed"
init_val_color = "black"

negligible_step_float = 1.0e-5
plotted_averages_number = 16

# Quantities for relative improvement calculation.
best_ref_vals = {
    "QM9": {
        "weak": {
            "dipole": (
                5.2624256900482225,
                8.190772053238254e-05,
                "8@7:7#2@8:7@5@6@8:7@6@7:7@7@8:6#3:6#1:6:6",
            ),
            "solvation_energy": (
                -0.036118152765830525,
                9.940944965402545e-05,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "atomization_energy": (
                -5.489771247476437,
                4.704659650874689e-05,
                "6#3@7:6#3@7:6#3@8:6#3@8:6#3@8:6#2@6@7:6#2@8:6#1:6",
            ),
        },
        "strong": {
            "dipole": (
                3.3568648712900764,
                0.06473921973977188,
                "8@7:8@8:7#1@7@8:7#1@5@7:7#1@6@8:6#2@6:6#2:6:6",
            ),
            "solvation_energy": (
                -0.021467490751841414,
                0.0005614388897913243,
                "8@7:8@8:7#2@7:7#2@8:7#1@5@7:6#2@6:6#2@8:6:6",
            ),
            "atomization_energy": (
                -5.489771247476437,
                4.704659650874689e-05,
                "6#3@7:6#3@7:6#3@8:6#3@8:6#3@8:6#2@6@7:6#2@8:6#1:6",
            ),
        },
    },
    "EGP": {
        "weak": {
            "dipole": (
                5.217191799780868,
                0.0022688571698954913,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "solvation_energy": (
                -0.03608624382949954,
                8.482093682841212e-05,
                "8@6:7#2@7:7#2@8:7#1@7@8:7@6@7:7@6@8:6:6:6",
            ),
            "atomization_energy": (
                -8.130926853176323,
                0.0015182550242717613,
                "6#3@8:6#3@8:6#3@9:6#3@9:6#3@10:6#3@10:6#3@11:6#3@11:6#1@12:6#1@12:6#1@13:6#1@13:6@13:6",
            ),
        },
        "strong": {
            "dipole": (
                3.1670498606166886,
                0.07476372139614328,
                "16@2@4@5@8:16@3@6@7@9:9:9:8:8:8:8:6#2@9:6#2",
            ),
            "solvation_energy": (
                -0.02976335641089712,
                0.0002888209390425601,
                "16@2@3@6@11:16@4@5@7@12:8:8:8:8:7#2:7#2:6#2@10@11:6#2@10@12:6#2:6#2:6#2",
            ),
            "atomization_energy": (
                -8.065834964640352,
                0.0007362348666310831,
                "7#2@12:7#2@13:6#2@4@12:6#2@5@13:6#2@6:6#2@7:6#2@8:6#2@9:6#2@10:6#2@11:6#2@11:6#2:6#2:6#2",
            ),
        },
    },
}


quant_STD = {
    "QM9": {
        "weak": {
            "dipole": 0.7202490150375122,
            "solvation_energy": 0.003800219983701787,
        },
        "strong": {
            "dipole": 0.5330517652144172,
            "solvation_energy": 0.0027043460958190825,
        },
    },
    "EGP": {
        "weak": {
            "dipole": 0.740092558470457,
            "solvation_energy": 0.0035497259204171077,
        },
        "strong": {
            "dipole": 0.5868571986460336,
            "solvation_energy": 0.003148923583866458,
        },
    },
}


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
            if self.nfigures == 0:
                return str(int(number_in))
            else:
                return "{:1.1f}".format(number_in)
        if number_in == 0.0:
            return "0"
        output = latex_scientific(number_in, nfigures=self.nfigures)
        output = r"$" + output + "$"
        return output


ytick_positions_tuples = {
    QM9: {
        "weak": {solv_en: (-100.0, 20.0, 2.0), dipole: (-100.0, 20.0, 2.0)},
        "strong": {solv_en: (-100.0, 20.0, 2.0), dipole: (-100.0, 20.0, 2.0)},
    },
    EGP: {
        "weak": {solv_en: (-20.0, 140.0, 20.0), dipole: (-20.0, 140.0, 10.0)},
        "strong": {solv_en: (-10.0, 140.0, 10.0), dipole: (-20.0, 140.0, 10.0)},
    },
}

ytick_label_format = {
    QM9: {
        solv_en: {"scientific": True, "nfigures": 1},
        dipole: {"scientific": False, "nfigures": 1},
    },
    EGP: {
        solv_en: {"scientific": True, "nfigures": 0},
        dipole: {"scientific": True, "nfigures": 0},
    },
}

# xtick_positions=[1.0, 1.e+1, 1.e+2, 1.e+3, 1.e+4, 1.e+5]
xtick_position_powers = [0, 2, 4]
xtick_positions = [10.0**p for p in xtick_position_powers]
xtick_labels = [r"$10^{" + str(p) + "}$" for p in xtick_position_powers]

display_xtick_labels = {
    QM9: {
        "weak": {solv_en: True, dipole: False},
        "strong": {solv_en: False, dipole: True},
    },
    EGP: {
        "weak": {solv_en: True, dipole: False},
        "strong": {solv_en: False, dipole: True},
    },
}

num_minor_ticks = 2

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

# For the legend.
bias_coeff_symbol = "\\alpha_{b}"
bias_true_val = {"none": 0.0, "weak": 0.2, "stronger": 0.4}


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


def plot_found_initial_value(ax, xlim, init_improvement, **add_plot_kwargs):
    ax.plot(
        list(xlim),
        [init_improvement for _ in xlim],
        linestyle=init_val_linestyle,
        linewidth=linewidth,
        color=init_val_color,
        **add_plot_kwargs
    )


def plot_initial_value(
    ax, input_filename, best_val_ref, STD_val_coeff, **add_plot_kwargs
):
    input = open(input_filename, "r")
    first_line = input.readline()
    input.close()
    _, init_val, _ = line2step_quant_SMILES(first_line)
    init_improvement = (init_val - best_val_ref) / STD_val_coeff
    xlim = ax.get_xlim()
    plot_found_initial_value(ax, xlim, init_improvement, **add_plot_kwargs)


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


def plot_ax_bias(ax, points_x, points_y, bias, yerr=None, **add_plot_kwargs):
    color = bias_linecolor[bias]
    linestyle = bias_linestyle[bias]
    plot_kwargs = {
        "color": color,
        "linestyle": linestyle,
        "linewidth": linewidth,
        **add_plot_kwargs,
    }
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
    input_filenames,
    bias_vals,
    ytick_positions_tuple,
    ytick_label_format,
    cur_best_ref_val,
    cur_val_STD_coeff,
    output_filename,
    cur_display_xtick_labels=True,
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
        points_x, unscaled_points_y = optimization_log_points(input_filename)
        points_y = [
            (unscaled_point_y - cur_best_ref_val) / cur_val_STD_coeff
            for unscaled_point_y in unscaled_points_y
        ]
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

    # Plot initial value.
    plot_initial_value(ax, input_filenames[0], cur_best_ref_val, cur_val_STD_coeff)

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
    if cur_display_xtick_labels:
        true_xtick_labels = xtick_labels
    else:
        true_xtick_labels = ["" for _ in xtick_labels]
    ax.set_xticklabels(true_xtick_labels, fontsize=ticks_fontsize)
    ax.yaxis.set_minor_locator(
        MultipleLocator(ytick_positions_tuple[2] / float(num_minor_ticks))
    )

    format_func = latex_format(**ytick_label_format)

    cur_ytick_positions = []
    cur_ytick_labels = []
    ytick_upper_bound = ytick_positions_tuple[1] + negligible_step_float
    pos = ytick_positions_tuple[0]
    tick_step = ytick_positions_tuple[2]
    while pos < ytick_upper_bound:
        if pos < final_min_y:
            pos += tick_step
            continue
        if pos > final_max_y:
            break
        cur_ytick_positions.append(pos)
        cur_ytick_labels.append(format_func(pos))
        pos += tick_step

    ax.set_yticks(cur_ytick_positions)
    ax.set_yticklabels(cur_ytick_labels, fontsize=ticks_fontsize)

    fig.set_figwidth(fig_width)
    fig.set_figheight(fig_height)

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

    # Specific values.
    cur_ytick_positions_tuple = ytick_positions_tuples[dataset][gap_constraint][
        quantity_name
    ]
    cur_ytick_label_format = ytick_label_format[dataset][quantity_name]
    cur_best_ref_val = best_ref_vals[dataset][gap_constraint][quantity_name][0]
    cur_val_STD_coeff = quant_STD[dataset][gap_constraint][quantity_name]
    if quantity_name == solv_en:
        cur_val_STD_coeff *= -1
    plot_opt_log_filenames(
        all_input_filenames,
        all_biases,
        cur_ytick_positions_tuple,
        cur_ytick_label_format,
        cur_best_ref_val,
        cur_val_STD_coeff,
        figure_file,
        cur_display_xtick_labels=display_xtick_labels[dataset][gap_constraint][
            quantity_name
        ],
    )


def plot_legend(output_dir, ncol):
    output_filename = output_dir + "/legend_" + str(ncol) + ".png"
    fig = plt.figure(constrained_layout=True)
    fig.set_figwidth(fig_width * 0.75)
    fig.set_figheight(fig_height * 0.75)
    legend_plot = fig.add_subplot()
    for bias in bias_values:
        bias_label = (
            r"$" + bias_coeff_symbol + "=" + "{:1.1f}".format(bias_true_val[bias]) + "$"
        )
        plot_ax_bias(legend_plot, [], [], bias, yerr=None, label=bias_label)
    plot_found_initial_value(legend_plot, [], 0.0, label="init. val.")
    # Re-assign linewidth.
    leg = legend_plot.legend()
    for line in leg.get_lines():
        line.set_linewidth(linewidth * 20.0)
    handles, labels = legend_plot.get_legend_handles_labels()
    legend_plot.axis("off")
    legend_plot.get_xaxis().set_visible(False)
    legend_plot.get_yaxis().set_visible(False)
    legend_plot.legend(
        handles,
        labels,
        numpoints=1,
        ncol=ncol,
        fontsize=legend_fontsize * 0.6,
        frameon=False,
        handlelength=2.0,
    )
    fig.savefig(output_filename)
    run(
        [
            "convert",
            output_filename,
            "-trim",
            output_dir + "/trimmed_legend_" + str(ncol) + ".png",
        ]
    )


def main():
    run(["mkdir", "-p", output_dir])
    for dataset in datasets:
        for quantity_name in quantity_names:
            for gap_constraint in gap_constraints:
                plot_opt_log_diff_bias(
                    data_dir, dataset, quantity_name, gap_constraint, output_dir
                )
    for ncol in [1, 2]:
        plot_legend(output_dir, ncol)


if __name__ == "__main__":
    main()
