import sys
from bmapqml.chemxpl.rdkit_draw_utils import draw_chemgraph_to_file
from bmapqml.chemxpl.rdkit_utils import SMILES_to_egc
from bmapqml.utils import mkdir
import subprocess

SMILES_filename = sys.argv[1]
output_filename = sys.argv[2]

sim_dict_SMILES = {}

with open(SMILES_filename, "r") as f:
    SMILES_lines = f.readlines()
    for l in SMILES_lines:
        try:
            cur_step = int(l)
            sim_dict_SMILES[cur_step] = []
        except:
            sim_dict_SMILES[cur_step].append(l[:-1])

rdkit_draw_kwargs = {
    "size": (300, 200),
    "bw_palette": True,
}

sim_dict_pictures = {}

pict_dir = "pictures"
mkdir(pict_dir)

MolPict_width = 2.0
Inf_width = 1.0
pict_height = 1.0
arrow_space = 0.5

arrow_vertical_shift = 0.5


class MolPict:
    def __init__(self, pict_path):
        self.base_init(MolPict_width)
        self.pict_path = pict_path
        self.x_shift = 0.0
        self.y_shift = 0.0

    def base_init(self, width):
        self.width = width
        self.height = pict_height
        self.position = None
        self.left_handle = None
        self.right_handle = None

    def tikz_node_int_string(self):
        return "\includegraphics[width=0.10\\textwidth]{" + self.pict_path + "}"

    def tikz_node(self, name, x_placement, y_placement):
        return (
            "\\node[anchor=north west, inner sep=0] ("
            + name
            + ") at ("
            + str(x_placement + self.x_shift)
            + ","
            + str(y_placement + self.y_shift)
            + ") {"
            + self.tikz_node_int_string()
            + "};\n"
        ), y_placement - self.height


class InfPict(MolPict):
    def __init__(self):
        self.base_init(Inf_width)
        self.x_shift = 0.0
        self.y_shift = -0.4

    def tikz_node_int_string(self):
        return "\Large{$\\mathbf{\\cdot\\cdot\\cdot}$}"


def arrow_contact_y(i):
    return -arrow_vertical_shift - pict_height * i


def mutarr_node(contact_tuple, x_placement, col_width):
    x0 = x_placement + col_width
    x1 = x0 + arrow_space
    y0 = arrow_contact_y(contact_tuple[0])
    y1 = arrow_contact_y(contact_tuple[1])
    return (
        "\\mutarr{" + str(x0) + "," + str(y0) + "}{" + str(x1) + "," + str(y1) + "};\n"
    )


def pict_list(SMILES_list, step_id):
    pict_list = []
    for s_id, s in enumerate(SMILES_list):
        cg = SMILES_to_egc(s).chemgraph
        pict_name = pict_dir + "/comp_" + str(step_id) + "_" + str(s_id) + ".png"
        draw_chemgraph_to_file(cg, pict_name, **rdkit_draw_kwargs)
        pict_list.append(MolPict(pict_name))
    return pict_list


for step_id, SMILES_list in sim_dict_SMILES.items():
    sim_dict_pictures[step_id] = pict_list(SMILES_list, step_id)

contact_arr_mutate = [(0, 0), (1, 1), (2, 2), (3, 3)]
contact_arr_crossover = [(0, 1), (1, 0), (2, 3), (3, 2)]

contact_arrs = [
    contact_arr_mutate,
    contact_arr_mutate,
    contact_arr_crossover,
    contact_arr_mutate,
]


def column_int():
    return [InfPict() for _ in range(4)]


columns = [column_int()]
for step_picts in sim_dict_pictures.values():
    columns.append(step_picts)
columns.append(column_int())

tex_output = """
\\documentclass[preview]{standalone}%
\\usepackage{amsmath,amsfonts,amsthm}
\\usepackage{placeins}
\\usepackage{tabularx,multirow}
\\usepackage[T1]{fontenc}
\\usepackage{graphicx, tikz}
\\usepackage[english]{babel}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\newcommand{\\gap}{\\Delta \\epsilon}
\\newcommand{\\dipole}{D}
\\newcommand{\\dEsolv}{\\Delta G_{\\mathrm{solv.}}}
\\usepackage{xcolor}
\\textwidth=16.8cm

\\usetikzlibrary{arrows.meta}
\\usetikzlibrary{shapes.misc}
\\usetikzlibrary{shapes.arrows}
\\usetikzlibrary{positioning}

% Thanks to https://tex.stackexchange.com/questions/197793/how-to-draw-gradient-arrows-with-tikz
\\makeatletter
\\def\\createshadingfromlist#1#2#3{%
  \\pgfutil@tempcnta=0\\relax
  \\pgfutil@for\\pgf@tmp:={#3}\\do{\\advance\\pgfutil@tempcnta by1}%
  \\ifnum\\pgfutil@tempcnta=1\\relax%
    \\edef\\pgf@spec{color(0)=(#3);color(100)=(#3)}%
  \\else%
    \\pgfmathparse{50/(\\pgfutil@tempcnta-1)}\\let\\pgf@step=\\pgfmathresult%
    %
    \\pgfutil@tempcntb=1\\relax%
    \\pgfutil@for\\pgf@tmp:={#3}\\do{%
      \\ifnum\\pgfutil@tempcntb=1\\relax%
        \\edef\\pgf@spec{color(0)=(\\pgf@tmp);color(25)=(\\pgf@tmp)}%
      \\else%
        \\ifnum\\pgfutil@tempcntb<\\pgfutil@tempcnta\\relax%
          \\pgfmathparse{25+\\pgf@step/4+(\\pgfutil@tempcntb-1)*\\pgf@step}%
          \\edef\\pgf@spec{\\pgf@spec;color(\\pgfmathresult)=(\\pgf@tmp)}%
        \\else%
          \\edef\\pgf@spec{\\pgf@spec;color(75)=(\\pgf@tmp);color(100)=(\\pgf@tmp)}%
        \\fi%
      \\fi%
      \\advance\\pgfutil@tempcntb by1\\relax%
    }%
  \\fi%
  \\csname pgfdeclare#2shading\\endcsname{#1}{100}\\pgf@spec%
}
\\createshadingfromlist{temparrshading}{vertical}{red,blue}


\\newcommand{\\mutarrend}[1]{Triangle[#1,fill=#1,scale=.5]}
\\newcommand{\\mutarr}[2]{\\draw [blue, arrows={-\\mutarrend{blue}}, line width=1.5mm] (#1) -- (#2)}

\\newcommand{\\steparr}[2]{\\draw [red, arrows={-\\mutarrend{red}}, line width=2.5mm] (#1) -- (#2)}
\\newcommand{\\temparr}[1]{\\node [shading=temparrshading,
  shape=double arrow,
  double arrow head extend=0.125cm, 
  shape border rotate=90, 
  minimum height=3.5cm, line width=1.5mm] (temparr) at (#1) {}}

\\begin{document}

\\fontsize{11}{2}\\selectfont
\\begin{figure}[tbp]
\\centering\\begin{tikzpicture}
"""

x_placement = 0.0
for col_id, column in enumerate(columns):
    y_placement = 0.0
    for pict_id, pict in enumerate(column):
        line, y_placement = pict.tikz_node(
            str(col_id) + "_" + str(pict_id), x_placement, y_placement
        )
        tex_output += line
    col_width = column[-1].width
    if col_id != len(columns) - 1:
        contact_arr = contact_arrs[col_id]
        for contact_tuple in contact_arr:
            tex_output += mutarr_node(contact_tuple, x_placement, col_width)
    x_placement += col_width + arrow_space

steparr_y = y_placement - 0.5 * pict_height
steparr_x = x_placement - arrow_space - 0.4
tex_output += (
    "\\steparr{"
    + str(0.0)
    + ","
    + str(steparr_y)
    + "}{"
    + str(steparr_x)
    + ","
    + str(steparr_y)
    + "};\n"
)

steparr_label_x = steparr_x / 2 - 1.25
steparr_label_y = steparr_y - 0.7

temparr_x = -0.75
temparr_y = y_placement

tex_output += (
    """\\node[anchor=south west] (steparrlabel) at ("""
    + str(steparr_label_x)
    + ","
    + str(steparr_label_y)
    + """) {optimization};
\\temparr{"""
    + str(temparr_x)
    + """,-2.1};
\\node[above left = 1.5cm and -1.5cm of temparr] (exploration) {exploration};
\\node[below left = 1.5cm and -1.5cm of temparr] (exploitation) {exploitation};
"""
)

tex_output += """
\\end{tikzpicture}
\\end{figure}
\\end{document}
"""

with open(output_filename, "w") as f:
    f.write(tex_output)

subprocess.run(["pdflatex", output_filename])

output_filename_prefix = output_filename[:-4]

pdf = output_filename_prefix + ".pdf"

final_png = output_filename_prefix + ".png"

subprocess.run(["pdftoppm", pdf, output_filename_prefix, "-png"])
subprocess.run(["mv", output_filename_prefix + "-1.png", final_png])
subprocess.run(
    [
        "convert",
        "-density",
        "250",
        final_png,
        "-trim",
        "-quality",
        "100",
        "final_" + final_png,
    ]
)
