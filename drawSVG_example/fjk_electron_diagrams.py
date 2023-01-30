# Example of how the arrow was drawn in the orbital representation paper.
import drawSvg as draw


class simple_arrow:
    def __init__(
        self,
        arrow_end_width,
        arrow_end_height,
        arrow_width,
        arrow_length,
        color="black",
    ):
        self.arrow_end_width = arrow_end_width
        self.arrow_end_height = arrow_end_height
        self.arrow_width = arrow_width
        self.arrow_length = arrow_length
        self.color = color
        self.arrow = draw.Marker(
            0.0,
            -arrow_end_width / 2,
            arrow_end_height,
            arrow_end_width / 2,
            scale=1.0,
            orient="auto",
        )  # TO-DO what is scale?
        self.arrow.append(
            draw.Lines(
                0.0,
                -arrow_end_width / 2,
                0.0,
                arrow_end_width / 2,
                arrow_end_height,
                0,
                fill=self.color,
                close=True,
            )
        )

    def vertical_arrow(self, pos_x, pos_y, D, upside=True):
        if upside:
            s = 1.0
        else:
            s = -1.0
        self.add_arrow(
            pos_x,
            pos_y - s * self.arrow_length / 2,
            pos_x,
            pos_y + s * (self.arrow_length / 2 - self.arrow_end_height * 2),
            D,
        )

    def add_arrow(self, x_beg, y_beg, x_end, y_end, D):
        D.d.append(
            draw.Line(
                x_beg,
                y_beg,
                x_end,
                y_end,
                stroke=self.color,
                stroke_width=self.arrow_width,
                fill="none",
                marker_end=self.arrow,
            )
        )

    def horizontal_arrow(self, pos_x, pos_y, D, to_the_right=True):
        if to_the_right:
            s = 1.0
        else:
            s = -1.0
        self.add_arrow(
            pos_x - s * self.arrow_length / 2,
            pos_y,
            pos_x + s * (self.arrow_length / 2 - self.arrow_end_height * 2),
            pos_y,
            D,
        )


class Drawing:
    def __init__(
        self,
        width=200,
        height=200,
        origin="center",
        displayInline=False,
        spin_arrow=None,
    ):
        self.d = draw.Drawing(width, height, origin="center", displayInline=False)
        self.spin_arrow = spin_arrow

    def add_state(self, pos_x, pos_y, occ=0, spin=0):
        self.d.append(
            draw.Line(
                pos_x - rect_width / 2,
                pos_y,
                pos_x + rect_width / 2,
                pos_y,
                stroke="black",
                stroke_width=plank_width,
            )
        )
        if occ != 0:
            if spin == 0:
                self.spin_arrow.vertical_arrow(
                    pos_x - occ_state_spin_spacing / 2, pos_y, self, upside=True
                )
                self.spin_arrow.vertical_arrow(
                    pos_x + occ_state_spin_spacing / 2, pos_y, self, upside=False
                )
            else:
                upside = spin == 1
                self.spin_arrow.vertical_arrow(pos_x, pos_y, self, upside)

    def add_diagram(self, pos_x, pos_ys, occs, spins):
        for pos_y, occ, spin in zip(pos_ys, occs, spins):
            self.add_state(pos_x, pos_y, occ=occ, spin=spin)

    def Text(self, text, fontsize, pos_x, pos_y):
        self.d.append(draw.Text(text, fontsize, pos_x, pos_y, center=0.6, fill="black"))

    def save(self, name):
        self.d.savePng(name)


rect_height = 20
rect_width = 20


spin_arrow_end_width = 3.0
spin_arrow_end_height = 2.5
spin_arrow_width = 3.0
plank_width = 3.0

axis_arrow_end_width = 3.0
axis_arrow_end_height = 6.0
axis_width = 3.0
axis_length = 185

big_arrow_end_width = 2.0
big_arrow_end_height = 2.0
big_arrow_width = 10.0
big_arrow_length = 25

fontsize = 25

occ_state_spin_spacing = 10.0

energies = [-75, -46, -23, 30, 57, 85]

occs = [2, 2, 2, 2, 0, 0]
occs_HOMO = [2, 2, 2, 1, 0, 0]
occs_LUMO = [2, 2, 2, 2, 1, 0]
occs_gap = [2, 2, 2, 1, 1, 0]

spins = [0, 0, 0, 0, 0, 0]
spins_HOMO = [0, 0, 0, 1, 0, 0]
spins_LUMO = [0, 0, 0, 0, 1, 0]
spins_gap = [0, 0, 0, 1, -1, 0]

energies_oxidation_1 = [-79, -51, -27, 26, 52, 80]
energies_oxidation_2 = [-70, -40, -18, 33, 62, 89]
spins_oxidation_1 = [1, 1, 1, 1, 0, 0]
spins_oxidation_2 = [-1, -1, -1, 0, 0, 0]
occs_oxidation_1 = [1, 1, 1, 1, 0, 0]
occs_oxidation_2 = [1, 1, 1, 0, 0, 0]

energies_reduction_1 = [-80, -51, -28, 27, 53, 81]
energies_reduction_2 = [-70, -42, -19, 35, 61, 88]
spins_reduction_1 = [1, 1, 1, 1, 1, 0]
spins_reduction_2 = [-1, -1, -1, -1, 0, 0]
occs_reduction_1 = [1, 1, 1, 1, 1, 0]
occs_reduction_2 = [1, 1, 1, 1, 0, 0]


energies_spin_change_1 = [-79, -51, -27, 28, 53, 82]
energies_spin_change_2 = [-71, -43, -20, 35, 61, 88]
spins_spin_change_1 = [1, 1, 1, 1, 1, 0]
spins_spin_change_2 = [-1, -1, -1, -1, 0, 0]
occs_spin_change_1 = [1, 1, 1, 1, 1, 0]
occs_spin_change_2 = [1, 1, 1, 0, 0, 0]


spin_arrow = simple_arrow(
    spin_arrow_end_width, spin_arrow_end_height, spin_arrow_width, rect_height
)

axis = simple_arrow(
    axis_arrow_end_width, axis_arrow_end_height, axis_width, axis_length
)

big_arrow = simple_arrow(
    big_arrow_end_width,
    big_arrow_end_height,
    big_arrow_width,
    big_arrow_length,
    color="blue",
)

plot_names = [
    "HOMO_change",
    "LUMO_change",
    "first_excitation",
    "oxidation",
    "reduction",
    "spin_change",
]

for plot_name in plot_names:
    D = Drawing(spin_arrow=spin_arrow)
    axis.vertical_arrow(-70, 0, D)
    D.add_diagram(-45, energies, occs, spins)
    D.Text("E", fontsize, -85, 80)
    big_arrow.horizontal_arrow(0, 0, D)
    if plot_name == "HOMO_change":
        D.add_diagram(60, energies, occs_HOMO, spins_HOMO)
    elif plot_name == "LUMO_change":
        D.add_diagram(60, energies, occs_LUMO, spins_LUMO)
    elif plot_name == "first_excitation":
        D.add_diagram(60, energies, occs_gap, spins_gap)
    elif plot_name == "oxidation":
        D.add_diagram(55, energies_oxidation_1, occs_oxidation_1, spins_oxidation_1)
        D.add_diagram(85, energies_oxidation_2, occs_oxidation_2, spins_oxidation_2)
    elif plot_name == "reduction":
        D.add_diagram(55, energies_reduction_1, occs_reduction_1, spins_reduction_1)
        D.add_diagram(85, energies_reduction_2, occs_reduction_2, spins_reduction_2)
    elif plot_name == "spin_change":
        D.add_diagram(
            55, energies_spin_change_1, occs_spin_change_1, spins_spin_change_1
        )
        D.add_diagram(
            85, energies_spin_change_2, occs_spin_change_2, spins_spin_change_2
        )
    D.save(plot_name + ".png")
# d.saveSvg('example.svg')
