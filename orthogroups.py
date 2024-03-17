import pathlib
import os 
import pandas as pd 
import numpy as np 
from collections import deque, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import functools
import itertools


def plot_heatmap(symmetric_matrix, cols, gc="unassined", title=None):
    symmetric_df = pd.DataFrame(symmetric_matrix, index=cols, columns=cols)
    plt.figure(figsize=(20, 14))
    ax = sns.heatmap(symmetric_df, annot=True, fmt=".0f", annot_kws={"size": 15}, cmap="viridis")
    ax.tick_params(axis='both', which='both', labelsize=15)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    cbar = ax.collections[0].colorbar
    
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.tick_params(labelsize=15)
    cbar.update_ticks()
    if not title:
        plt.title(f"Orthogroups intersection {gc} (%)", fontsize = 20)
    else:
        plt.title(f"Orthogroups: {title} (%)", fontsize = 20)
    output_path = pathlib.Path.cwd() / "orthogroups_compare" / f"{gc}0.png"
    plt.savefig(output_path)
    plt.close()


def compare_diff(left_sets, right_sets):

    distinct_orthogroups = set()

    if left_sets == right_sets:
        return left_sets, 0, 0

    left_diffs = left_sets - right_sets 
    right_diffs = right_sets - left_sets
    
    occurred = set()
    gen_intersection_dict = defaultdict(lambda: defaultdict(dict))
    prev = set()

    for left, right in itertools.product(left_diffs, right_diffs):

        left_gens = len(left)
        right_gens = len(right)

        intersection = left & right

        if len(intersection) > 0:
            occurred.add(right)

        if left not in prev:
            max_intersection = 0
            max_intersection_union = 0
            max_intersection_left = 0
            max_intersection_right = 0
            prev.add(left)

        max_intersection = np.amax((max_intersection, len(intersection)))
        max_intersection_union = np.amax((max_intersection_union, max_intersection / (left_gens + right_gens - len(intersection))))
        max_intersection_left = np.amax((max_intersection_left, max_intersection / left_gens))
        max_intersection_right = np.amax((max_intersection_right, max_intersection / right_gens))


        gen_intersection_dict[left]["max_intersection"] = max_intersection
        gen_intersection_dict[left]["max_intersection_divide_union"] = max_intersection_union
        gen_intersection_dict[left]["max_intersection_divide_left"]  = max_intersection_left
        gen_intersection_dict[left]["max_intersection_divide_right"]  = max_intersection_right
    
    distinct_orthogroups = left_sets | (right_diffs - occurred)       
    distinct_right, orthogroup_overlap = right_diffs - occurred, occurred

    return distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict


orthofinder_outputs = pathlib.Path.cwd() / "ExampleData0" / "OrthoFinder"
species_fasta = pathlib.Path.cwd() / "ExampleData0"
species_gens = {file.name.split(".")[0]: {} for file in species_fasta.iterdir() if file.is_file() }

orthogroup_file_dict = {}
unassigned_gens_file_dict = {}
for output_dir in orthofinder_outputs.iterdir():
    output_dir_name = output_dir.name 
    output_matrix_result = output_dir_name.split("_")[-1]
    
    for category_dir in output_dir.iterdir():
        if category_dir.name == "Orthogroups":
            orthogroup_file_dict[output_matrix_result]  = category_dir / "Orthogroups.tsv"
            unassigned_gens_file_dict[output_matrix_result] = category_dir / "Orthogroups_UnassignedGenes.tsv"

def join_rows(row):
    row_items = set()
    for label, item in row.items():
        if item != "nan":
            for gen in item.split(", "):
                row_items.add((label, gen))
    return frozenset(row_items)
        

def create_ogs_sets(file_dict):
    orthogroup_sets_dict = {}
    for matrix_name, file in file_dict.items():
        matrix_og_sets = set()

        with open(file) as reader:
            for i, line in enumerate(reader):
                line = line.split("\t")[1:]
                line[-1] = line[-1].replace("\n", "")

                if i == 0:
                    species = line
                else:
                    gens = set()
                    for i in range(len(line)):
                        if len(line[i]) != 0:
                            for l in line[i].split(", "):
                                if len(l) != 0:
                                    gens.add((species[i], l))
                    
                    gens = frozenset(gens)
                    matrix_og_sets.add(gens)

        if all(len(item) == 0 for item in matrix_og_sets):
            print(f"The implementation of the {matrix_name} did not generate any Orthogroups!")
        else:
            orthogroup_sets_dict[matrix_name] = matrix_og_sets
    return orthogroup_sets_dict


orthogroup_sets_dict = create_ogs_sets(orthogroup_file_dict)
matrices = [*orthogroup_sets_dict.keys()]

intersection_list0 = []
left_diff_list0 = []
right_diff_list0 = []

distinct_right_list0 = []
overlap_right_list0 = []

updated_left_orthogroups = {}
orthogroups_similarity_dict = {}
orthogroups_diff_per_dict = {}

for i in range(len(matrices)):
    left_set = orthogroup_sets_dict[matrices[i]]
    
    intersection_list1 = deque()
    left_diff_list1 = deque()
    right_diff_list1 = deque()

    distinct_right_list1 = deque()
    overlap_right_list1 = deque()
    
    for j in range(i+1, len(matrices)):
        print(matrices[i], matrices[j])
        right_set = orthogroup_sets_dict[matrices[j]]
        
        intersection = left_set & right_set
        left_diff = left_set - right_set
        right_diff = right_set - left_set

        intersection_num = len(intersection) 
        left_diff_num = len(left_diff)
        right_diff_num = len(right_diff)

        intersection_list1.append(intersection_num)
        left_diff_list1.append(left_diff_num)
        right_diff_list1.append(right_diff_num)

        distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict = compare_diff(left_diff.copy(), right_diff.copy())
        distinct_right_num, orthogroup_overlap_num = len(distinct_right), len(orthogroup_overlap)
        updated_left_orthogroups["-".join((matrices[i],matrices[j]))] = left_set | distinct_orthogroups
        orthogroups_similarity_dict["-".join((matrices[i],matrices[j]))] = gen_intersection_dict
       
        print(len(left_diff), len( right_diff), len(gen_intersection_dict))
        print(len(left_set), len(right_set), intersection_num, left_diff_num, right_diff_num, distinct_right_num, orthogroup_overlap_num)
        
        distinct_right_list1.append(100 * distinct_right_num / right_diff_num)
        overlap_right_list1.append(100 * orthogroup_overlap_num / right_diff_num)

    while len(intersection_list1) < len(matrices):
        intersection_list1.appendleft(0)
    
    while len(left_diff_list1) < len(matrices):
        left_diff_list1.appendleft(0)
    
    while len(right_diff_list1) < len(matrices):
        right_diff_list1.appendleft(0)

    while len(distinct_right_list1) < len(matrices):
        distinct_right_list1.appendleft(0)
    
    while len(overlap_right_list1) < len(matrices):
        overlap_right_list1.appendleft(0)
    
    intersection_list0.append(np.array(intersection_list1))
    left_diff_list0.append(np.array(left_diff_list1))
    right_diff_list0.append(np.array(right_diff_list1))

    distinct_right_list0.append(np.array(distinct_right_list1))
    overlap_right_list0.append(np.array(overlap_right_list1))


def create_symmetric_matrix(a_list):

    arr = np.stack(a_list)
    symmetric_matrix = arr + arr.T - np.diag(arr.diagonal())

    return symmetric_matrix


intersection_arr = create_symmetric_matrix(intersection_list0)
left_diff_arr = create_symmetric_matrix(left_diff_list0)
right_diff_arr = create_symmetric_matrix(right_diff_list0)
distinct_right_arr = create_symmetric_matrix(distinct_right_list0)
overlap_right_arr = create_symmetric_matrix(overlap_right_list0)

intersection_left = 100 * intersection_arr / (intersection_arr + left_diff_arr) 
intersection_right = 100 * intersection_arr / (intersection_arr + right_diff_arr)
intersection_union = 100 * intersection_arr / (intersection_arr + left_diff_arr + right_diff_arr)


plot_heatmap(intersection_left, matrices, gc="divided by left")
plot_heatmap(intersection_right, matrices, gc="divided by right")
plot_heatmap(intersection_union, matrices, gc="divided by both")
# plot_heatmap(distinct_right_arr, matrices, 
#              gc="distinct_right_ratio", title="distinct gens in the right group")
# plot_heatmap(overlap_right_arr, matrices, 
#              gc="overlap_right_ratio", title="shared gens in the right group")
                

def create_unassgined_gens_sets(file_dict):
    unassigned_gens_sets = {}

    for matrix_name, file in file_dict.items():
        with open(file) as reader:
            for i, line in enumerate(reader):
                line = line.split("\t")[1:]
                line[-1] = line[-1].replace("\n", "")

                if i == 0:
                    species = line
                    for l in species:
                        if l not in unassigned_gens_sets:
                            unassigned_gens_sets[l] = {}

                        unassigned_gens_sets[l][matrix_name] = set()

                else:
                    for i in range(len(line)):
                        if len(line[i]) != 0:
                            unassigned_gens_sets[species[i]][matrix_name].add(line[i])

    return unassigned_gens_sets 

print()
print("--------------------- Unassigned -------------------------------")
print()
unassigned_gens_sets = create_unassgined_gens_sets(unassigned_gens_file_dict)


# unassigned_gens_arr = {species: {} for species in unassigned_gens_sets}
for species, unassigned_set in unassigned_gens_sets.items():
    matrices = [*unassigned_set.keys()]

    unassigned_intersection_list0 = []
    unassigned_left_diff_list0 = []
    unassigned_right_diff_list0 = []
   
    for i in range(len(matrices)):
        unassigned_left_set = unassigned_set[matrices[i]]
        unassigned_intersection_list1 = deque()
        unassigned_left_diff_list1 = deque()
        unassigned_right_diff_list1 = deque()

        for j in range(i+1, len(matrices)):
            print(matrices[i], matrices[j])
            unassigned_right_set = unassigned_set[matrices[j]]

            unassigned_intersection = unassigned_left_set & unassigned_right_set
            unassigned_left_diff = unassigned_left_set - unassigned_right_set
            unassigned_right_diff = unassigned_right_set - unassigned_left_set

            unassigned_intersection_num = len(unassigned_intersection) 
            unassigned_left_diff_num = len(unassigned_left_diff)
            unassigned_right_diff_num = len(unassigned_right_diff)

            unassigned_intersection_list1.append(unassigned_intersection_num)
            unassigned_left_diff_list1.append(unassigned_left_diff_num)
            unassigned_right_diff_list1.append(unassigned_right_diff_num)

        while len(unassigned_intersection_list1) < len(matrices):
            unassigned_intersection_list1.appendleft(0)
        
        while len(unassigned_left_diff_list1) < len(matrices):
            unassigned_left_diff_list1.appendleft(0)
        
        while len(unassigned_right_diff_list1) < len(matrices):
            unassigned_right_diff_list1.appendleft(0)

        unassigned_intersection_list0.append(np.array(unassigned_intersection_list1))
        unassigned_left_diff_list0.append(np.array(unassigned_left_diff_list1))
        unassigned_right_diff_list0.append(np.array(unassigned_right_diff_list1))

    unassigned_intersection_arr = create_symmetric_matrix(unassigned_intersection_list0)
    unassigned_left_diff_arr = create_symmetric_matrix(unassigned_left_diff_list0)
    unassigned_right_diff_arr = create_symmetric_matrix(unassigned_right_diff_list0)

    unassigned_intersection_left = 100 * unassigned_intersection_arr / (unassigned_intersection_arr + unassigned_left_diff_arr) 
    unassigned_intersection_right = 100 * unassigned_intersection_arr / (unassigned_intersection_arr + unassigned_right_diff_arr)
    unassigned_intersection_union = 100 * unassigned_intersection_arr / (unassigned_intersection_arr + unassigned_left_diff_arr + unassigned_right_diff_arr)

    plot_heatmap(unassigned_intersection_left, matrices, gc=f"{species}_(unassigned gens) divided by left")
    plot_heatmap(unassigned_intersection_right, matrices, gc=f"{species}_(unassigned gens) divided by right")
    plot_heatmap(unassigned_intersection_union, matrices, gc=f"{species}_(unassigned gens) divided by both")
