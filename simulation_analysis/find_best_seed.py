import sys, os
from misc_procedures import find_seed_with_best_candidate
import pdb

simulation_names = ["QM9"]


for name in simulation_names:

    paths = os.listdir("/data/jan/konstantin/{}".format(name))
    print(paths)

    best_seeds = open("best_seeds_{}.txt".format(name), "w")
    for p in paths:
        #pdb.set_trace()
        best = find_seed_with_best_candidate("/data/jan/konstantin/{}/{}".format(name,p))
        print(p, best)
        best_seeds.write("{}\n".format(best))

    best_seeds.close()


