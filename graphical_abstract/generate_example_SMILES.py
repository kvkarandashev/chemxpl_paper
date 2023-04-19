from bmapqml.chemxpl.random_walk import (
    randomized_change,
    randomized_cross_coupling,
    TrajectoryPoint,
    ChemGraph,
    add_heavy_atom_chain,
    replace_heavy_atom,
    remove_heavy_atom,
    change_bond_order,
)
from bmapqml.chemxpl.rdkit_utils import chemgraph_to_canonical_rdkit  # SMILES_to_egc
from bmapqml.chemxpl.minimized_functions.morfeus_quantity_estimates import (
    morfeus_FF_xTB_code_quants,
)
from bmapqml.utils import loadpkl, mkdir
import random, sys
import numpy as np


class BadChange(Exception):
    pass


def solv_en(tp):
    if isinstance(tp, ChemGraph):
        true_tp = TrajectoryPoint(cg=tp)
    else:
        true_tp = tp
    quant = "solvation_energy"
    true_tp.calculated_data = {}
    output = morfeus_FF_xTB_code_quants(
        true_tp,
        num_conformers=16,
        num_attempts=16,
        ff_type="MMFF94",
        quantities=[quant],
        solvent="water",
        remaining_rho=0.9,
    )["mean"][quant]
    if output is None:
        raise BadChange
    return output


change_kwargs = {
    "nhatoms_range": [1, 12],
    "forbidden_bonds": [(7, 7), (7, 8), (7, 9), (8, 8), (8, 9), (9, 9)],
    "possible_elements": ["C", "N", "O", "F"],
    "bond_order_changes": [-1, 1],
}

nhatom_change_list = [add_heavy_atom_chain, remove_heavy_atom]

mutation_lists = [
    nhatom_change_list,
    [replace_heavy_atom],
    nhatom_change_list,
    [change_bond_order],
]

ind_step_tries = 20


def mutations_single_change(tp, mutation_list, i):
    if i > 1:
        init_en = solv_en(tp)
    tp.possibility_dict = None
    for _ in range(ind_step_tries):
        try:
            trial_tp, prob_balance = randomized_change(
                tp,
                change_prob_dict=mutation_list,
                delete_chosen_mod_path=True,
                **change_kwargs
            )
        except:
            continue
        if prob_balance is None:
            continue
        trial_en = solv_en(trial_tp)
        if trial_en is None:
            continue
        if i > 1:
            if trial_en > init_en:
                continue
        return trial_tp
    raise BadChange


def mutations_change(tp_list):
    output = []
    for i, (init_tp, mutation_list) in enumerate(zip(tp_list, mutation_lists)):
        output.append(mutations_single_change(init_tp, mutation_list, i))

    if output[0].chemgraph().nhatoms() == output[2].chemgraph().nhatoms():
        raise BadChange
    return output


def ind_crossover_change(cg_tuple, en, i):
    for _ in range(ind_step_tries):
        new_pair, prob_balance = randomized_cross_coupling(cg_tuple, **change_kwargs)
        if prob_balance is None:
            continue
        if solv_en(new_pair[1]) > solv_en(new_pair[0]):
            new_pair = new_pair[::-1]
        if i == 1:
            if en < solv_en(new_pair[1]):
                continue
        return [TrajectoryPoint(cg=cg) for cg in new_pair]
    raise BadChange


def crossover_change(tp_list):
    output = []
    for i in range(len(tp_list) // 2):
        cg1 = tp_list[i * 2].chemgraph()
        cg2 = tp_list[i * 2 + 1].chemgraph()
        output += ind_crossover_change(
            (cg1, cg2), min(solv_en(tp_list[i * 2]), solv_en(tp_list[i * 2 + 1])), i
        )
    return output


random.seed(1)
np.random.seed(1)

num_tries = 1000

SMILES_dir = "example_SMILES"
mkdir(SMILES_dir)

change_list = [mutations_change, crossover_change]


def make_step_dict(init_tp_list):
    init_tp_tuples = []
    for tp in init_tp_list:
        tp.possibility_dict = None
        en = solv_en(tp)
        if en is None:
            raise BadChange
        init_tp_tuples.append((tp, en))
    init_tp_tuples.sort(key=lambda x: x[1], reverse=True)

    cur_tp_list = [t[0] for t in init_tp_tuples]

    step_dict = {0: cur_tp_list}
    for i, change in enumerate(change_list):
        print("Changing step:", i)
        cur_tp_list = change(cur_tp_list)
        step_dict[i + 1] = cur_tp_list
    return step_dict


pkl_file = sys.argv[1]
restart_object = loadpkl(pkl_file)
histogram = restart_object["histogram"]

for try_id in range(num_tries):
    print(try_id)
    init_tp_list = random.sample(histogram, len(mutation_lists))
    print("sampled")
    try:
        step_dict = make_step_dict(init_tp_list)
    except BadChange:
        continue
    with open(SMILES_dir + "/SMILES_example_" + str(try_id) + ".txt", "w") as f:
        for step_id, tp_list in step_dict.items():
            print(step_id, file=f)
            for tp in tp_list:
                print(
                    chemgraph_to_canonical_rdkit(tp.chemgraph(), SMILES_only=True),
                    file=f,
                )
