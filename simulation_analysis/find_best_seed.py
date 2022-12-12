import sys, os
from misc_procedures import find_seed_with_best_candidate


paths = os.listdir("/data/jan/konstantin")
print(paths)

best_seeds = open("best_seeds.txt", "w")
for p in paths:
    best = find_seed_with_best_candidate("/data/jan/konstantin/{}".format(p))
    print(p, best)
    best_seeds.write("{}\n".format(best))

best_seeds.close()


