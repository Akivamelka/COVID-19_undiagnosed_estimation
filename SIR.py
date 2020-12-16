import numpy as np
import math
import random

num_leaves = (int)(math.pow(2, 17)) # Total population
gamma_list = [0.3, 0.5, 0.7] # detection rate
gamma_slope = 0.0
num_initial_infected = 10
max_num_step = 1000000
traget_removed = 1000
averaging_length = 1000
beta = 1
beta_step = 1000
num_simul = 5
file = f"SIR.csv"
beta_fix = True

def compute_betas(slope):

    distribution = range(1, beta_step + 1)
    distribution = [1 / np.power(x, slope) for x in distribution]
    sumdist = sum(distribution)
    distribution = [x / sumdist for x in distribution]

    cumsum = 0
    cumdist = []
    i = 0
    while i < len(distribution):
        cumsum += distribution[i]
        cumdist.append(cumsum)
        i += 1

    betas = []
    k = 0
    while k < num_leaves:
        var = random.random()
        if var > 0 and var < cumdist[0]:
            betas.append(1)

        else:
            h = 1
            while h < len(cumdist):
                if var > cumdist[h - 1] and var < cumdist[h]:
                    betas.append(h + 1)
                h += 1
        k += 1

    ss = sum(betas)
    betas = [x / ss for x in betas]
    return betas

def create_population(betas):

    tree_infected = [0] * (2 * num_leaves - 1)  # Tree of infected
    tree_susceptible = [0] * (2 * num_leaves - 1)  # Tree of betas

    tree_susceptible[num_leaves - 1:] = betas # [1] * num_leaves
    tree_infected[num_leaves - 1: num_leaves - 1 + num_initial_infected] = [1] * num_initial_infected
    tree_susceptible[num_leaves - 1: num_leaves - 1 + num_initial_infected] = [0] * num_initial_infected

    i = num_leaves - 2
    while i >= 0:
        tree_infected[i] = tree_infected[2 * i + 1] + tree_infected[2 * i + 2]
        tree_susceptible[i] = tree_susceptible[2 * i + 1] + tree_susceptible[2 * i + 2]
        i -= 1

    repertoire = {x + num_leaves - 1: num_leaves - 1 for x in range(num_initial_infected)}

    return tree_infected, tree_susceptible, repertoire

def correct_tree(leaf, tree):

    i = (leaf - 1) // 2
    while i >= 0:
        tree[i] = tree[2 * i + 1] + tree[2 * i + 2]
        i -= 1
        i //= 2

def choose_leaf(tree):

    var = random.random() * tree[0]
    if var == 0:
        var += 1
    leaf = 0
    while leaf < len(tree):
        if var > tree[leaf]:
            var -= tree[leaf]
            leaf += 1
        leaf = 2 * leaf + 1

    return (leaf - 1) // 2

def process():

    graph = []

    if beta_fix:
        betas = [1] * num_leaves
    else:
        betas = compute_betas(2)

    for gamma in gamma_list:

        for simul in range(num_simul):

            tree_infected, tree_susceptible, repertoire = create_population(betas)
            scale = tree_susceptible[0]
            print("gamma = ", gamma)
            num_infected = num_initial_infected
            num_removed = 0
            num_removed_detected = 0

            graph1 = []
            graph2 = []
            graph3 = []
            graph4 = []

            step = 0
            while num_removed <= traget_removed and step < max_num_step and tree_infected[0] > 0: # num_removed <= traget_removed and

                a = tree_susceptible[0] / scale
                gamma_eff = gamma + gamma_slope * step
                proba = beta * a  / (beta * a + gamma_eff)
                rand_var = random.random()

                if rand_var < proba:
                    # Infection
                    leaf = choose_leaf(tree_susceptible)
                    tree_susceptible[leaf] = 0
                    correct_tree(leaf, tree_susceptible)

                    # Choose contaminator
                    leaf_contaminator = choose_leaf(tree_infected)
                    repertoire[leaf] = leaf_contaminator

                    # Update tree
                    tree_infected[leaf] = 1
                    correct_tree(leaf, tree_infected)
                    num_infected += 1

                else:
                    # Detection
                    leaf = choose_leaf(tree_infected)
                    tree_infected[leaf] = 0
                    correct_tree(leaf, tree_infected)
                    num_removed += 1

                    # Checking
                    leaf_contaminator = repertoire[leaf]
                    if tree_infected[leaf_contaminator] == 0:
                        num_removed_detected += 1
                    del repertoire[leaf]

                if num_removed_detected > 0: # and num_infected > 100:
                    graph1.append(num_infected)
                    graph2.append(num_removed)
                    graph3.append(num_removed_detected)

                if step % 1000 == 0:
                    print(step)
                step +=1

            print("final step is ", step)
            print("number of removed is ", num_removed)
            print("number of detected is ", num_removed_detected)
            print("number of infected is ", num_infected)
            print("number of infected remaining is ", tree_infected[0])

            mean_ccf = 0
            mean_ksf = 0

            if len(graph1) > averaging_length:
                infected = np.array(graph1[-averaging_length:])
                removed = np.array(graph2[-averaging_length:])
                detected = np.array(graph3[-averaging_length:])
                mean_ccf = np.mean(infected / removed)
                mean_ksf = np.mean(removed / detected - 1)

            graph.append([simul, gamma, mean_ksf, mean_ccf])

    graph_save = np.asmatrix(graph)
    np.savetxt(file, graph_save, delimiter=",")

process()
