choice_order = {
    "add node (\elementarychange{1}a)": [
        "type of heavy element to be added",
        "node to which the new node will be connected with a covalent bond",
        "order of the new covalent bond",
    ],
    "remove node (\elementarychange{1}b)": [
        "type of heavy element to be removed",
        "removed node",
    ],
    "change bond order (\elementarychange{2})": [
        "step by which bond order is changed",
        "altered pair of nodes",
        "resonance structure",
    ],
    "replace heavy atom (\elementarychange{3})": [
        "type of heavy atom to be created",
        "changed node",
    ],
    "change valence / change hydrogen number (\elementarychange{4})": [
        "node whose valence is to be changed",
        "new number of hydrogens connected to the node",
    ],
    "change valence / add heavy atoms (\elementarychange{5}a)": [
        "element of the added nodes",
        "node to which the created nodes will be connected via covalent bonds",
        "order of the new covalent bonds",
    ],
    "change valence / remove heavy atoms (\elementarychange{5}b)": [
        "removed nodes element",
        "node whose neighbors will be removed",
        "order of covalent bonds connecting removed nodes to the molecule",
    ],
    "change valence / change bond order (\elementarychange{6})": [
        "step by which bond order is changed",
        "altered pair of nodes",
        "resonance structure",
    ],
}

output_name = "mutation_choice_order.tex"

output = open(output_name, "w")

print(
    """
\\begin{tabular}{p{0.25\\textwidth}|p{0.7\\textwidth}}
\\toprule""",
    file=output,
)


skip = "2.0ex"
skip_after = "3.0ex"


def begenum(empty_label=False):
    if empty_label:
        int_s = ",labelwidth=0pt,leftmargin=0pt"
    else:
        int_s = ""
    return (
        "\\vspace{-" + skip + "}\\begin{enumerate}[noitemsep,topsep=0pt" + int_s + "]"
    )


endenum = "\end{enumerate}\n\n\\phantom{0}\\vspace{-" + skip_after + "}"


def inenum(s, empty_label=False):
    return begenum(empty_label=empty_label) + s + endenum


headers = ["Procedure", "Choice order"]

for i, h in enumerate(headers):
    output.write(inenum("\item[]" + " " + h + " ", empty_label=(i == 0)) + " ")
    if i == len(headers) - 1:
        output.write("\\\\\n ")
    else:
        output.write(" & ")

for k, vals in choice_order.items():
    output.write("\midrule " + inenum(" \item[] " + k, empty_label=True) + " & ")
    l = ""
    for v in vals:
        l += "\item " + v + " "
    output.write(inenum(l, empty_label=False) + "\\\\\n ")

output.write("\\bottomrule ")

output.write("\end{tabular} ")
output.close()
