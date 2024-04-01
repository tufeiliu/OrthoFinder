import pathlib
import os 
import pandas as pd 
import numpy as np 
from collections import deque, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import functools
import itertools
import sys
import multiprocessing
import copy


def plot_heatmap(symmetric_matrix, cols, data_dirname, output_filename,
                 gc="orthogroup", analysis_type="intersection",
                 save_to="orthogroups_compare", title=None):

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
        plt.title(f"{gc}: {analysis_type} %", fontsize = 20)
    else:
        plt.title(title, fontsize = 20)
    output_dir = pathlib.Path.cwd() / "_".join((data_dirname, save_to))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_path = output_dir / f"{output_filename}.png"
    plt.savefig(output_path)
    plt.close()


def compare_diff(left_sets, right_sets):

    distinct_orthogroups = set()

    if left_sets == right_sets:
        return left_sets, set(), set(), {}

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


def create_ogs_sets(matrix_name, file):

    matrix_og_sets = set()
    if file is not None:
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
    

    if all(len(item) == 0 for item in matrix_og_sets) or len(matrix_og_sets) == 0:
        print(f"The implementation of the {matrix_name} did not generate anything!")
        return
    else:
        return matrix_og_sets


def compute_expected(left_diff_num, intersection_num, right_diff_num, gen_intersection_dict):
    
    left_total = left_diff_num + intersection_num
    total = left_diff_num + intersection_num + right_diff_num
    union_ratio_total = intersection_num
    # left_ratio_total = intersection_num
    # right_ratio_total = intersection_num
    
    for left, val in gen_intersection_dict.items():
        union_ratio_total += val["max_intersection_divide_union"]
        # left_ratio_total += val["max_intersection_divide_left"]
        # right_ratio_total += val["max_intersection_divide_right"]

    total_ratio = union_ratio_total / total
    # total_ratio = left_ratio_total / total
    expected = total_ratio * left_total 
    return total_ratio * 100, expected


def create_symmetric_matrix(a_list):

    arr = np.stack(a_list)
    symmetric_matrix = arr + arr.T - np.diag(arr.diagonal())

    return symmetric_matrix


def create_unassgined_gens_matrix_sets(file_dict):

    unassigned_gens_matrix_sets = {}

    for matrix_name, file in file_dict.items():
        unassigned_gens_matrix_sets[matrix_name] = set()
        with open(file) as reader:
            for i, line in enumerate(reader):
                line = line.split("\t")[1:]
                line[-1] = line[-1].replace("\n", "")

                if i == 0:
                    species = line
                else:
                    for i in range(len(line)):
                        if len(line[i]) != 0:
                            unassigned_gens_matrix_sets[matrix_name].add((species[i], line[i]))
   
    return unassigned_gens_matrix_sets


def obtain_ogs_datafiles(datadir):

    orthofinder_outputs = pathlib.Path.cwd() / datadir / "OrthoFinder"
    species_fasta = pathlib.Path.cwd() / datadir

    orthogroup_file_dict = {}
    unassigned_gens_file_dict = {}
    orthologue_file_dict = {}
    for output_dir in orthofinder_outputs.iterdir():
        output_dir_name = output_dir.name 
        output_matrix_result = output_dir_name.split("_")[-1]
        orthologue_file_dict[output_matrix_result] = []
        orthogroup_file_dict[output_matrix_result] = None
        unassigned_gens_file_dict[output_matrix_result] = None

        for category_dir in output_dir.iterdir():
            if category_dir.name == "Orthogroups":
                orthogroup_file_dict[output_matrix_result]  = category_dir / "Orthogroups.tsv"
                unassigned_gens_file_dict[output_matrix_result] = category_dir / "Orthogroups_UnassignedGenes.tsv"
                
            if category_dir.name == "Orthologues":
                for species_dir in category_dir.iterdir():
                    if not species_dir.is_file():
                        for file in os.listdir(species_dir):
                            orthologue_file_dict[output_matrix_result].append(species_dir / file)

    return  orthogroup_file_dict, unassigned_gens_file_dict, orthologue_file_dict


def orthogroups_sets_compare(left_set, right_set, left_matrix, right_matrix):

    intersection = left_set & right_set
    left_diff = left_set - right_set
    right_diff = right_set - left_set

    intersection_num = len(intersection) 
    left_diff_num = len(left_diff)
    right_diff_num = len(right_diff)

    distinct_orthogroups, distinct_right, orthogroup_overlap, gen_intersection_dict = compare_diff(left_diff.copy(), right_diff.copy())

    distinct_right_num, orthogroup_overlap_num = len(distinct_right), len(orthogroup_overlap)
    updated_left_orthogroups = left_set | distinct_orthogroups
    similarity, expected = compute_expected(left_diff_num, intersection_num, right_diff_num, gen_intersection_dict)

    # print(left_matrix, right_matrix)
    # print(len(left_diff), len( right_diff), len(gen_intersection_dict), similarity, expected)
    # print(len(left_set), len(right_set), intersection_num, left_diff_num, right_diff_num, distinct_right_num, orthogroup_overlap_num)

    results = (intersection_num, left_diff_num, right_diff_num, \
               similarity, expected, left_matrix, right_matrix)

    return results

def get_arr_division(numerator_arr, denominator_arr):

    if np.any(denominator_arr == 0):
        arr = np.empty(numerator_arr.shape)
        arr.fill(np.nan)
        non_zero_indices = np.where(denominator_arr != 0)
        arr[non_zero_indices] = 100 * numerator_arr[non_zero_indices] / denominator_arr[non_zero_indices]
    else:
        arr = 100 * numerator_arr / denominator_arr

    return arr

        
def gen_sets_analysis(matrices, gens_matrix_sets):

    intersection_list0 = []
    left_diff_list0 = []
    right_diff_list0 = []

    for i in range(len(matrices)):
        left_set = gens_matrix_sets[matrices[i]]

        intersection_list1 = deque()
        left_diff_list1 = deque()
        right_diff_list1 = deque()

        for j in range(i+1, len(matrices)):
            right_set = gens_matrix_sets[matrices[j]]

            intersection = left_set & right_set 
            left_diff = left_set - right_set 
            right_diff = right_set - left_set 

            intersection_num = len(intersection)
            left_diff_num = len(left_diff)
            right_diff_num = len(right_diff)
            
            intersection_list1.append(intersection_num)
            left_diff_list1.append(left_diff_num)
            right_diff_list1.append(right_diff_num)

        while len(intersection_list1) < len(matrices):
            intersection_list1.appendleft(0)
        
        while len(left_diff_list1) < len(matrices):
            left_diff_list1.appendleft(0)
        
        while len(right_diff_list1) < len(matrices):
            right_diff_list1.appendleft(0)
  

        intersection_list0.append(np.array(intersection_list1))
        left_diff_list0.append(np.array(left_diff_list1))
        right_diff_list0.append(np.array(right_diff_list1))

    intersection_arr = create_symmetric_matrix(intersection_list0)
    left_diff_arr = create_symmetric_matrix(left_diff_list0)
    right_diff_arr = create_symmetric_matrix(right_diff_list0)

    intersection_left = get_arr_division(intersection_arr, left_diff_arr)
    intersection_right = get_arr_division(intersection_arr, right_diff_arr)
    intersection_union = get_arr_division(intersection_arr, intersection_arr + left_diff_arr + right_diff_arr)
        
    return intersection_left, intersection_right, intersection_union

def sort_output(output_dict, ordered_matrix):
        
        ordered_output = {}
        for left_matrix, val in output_dict.items():
            ordered_val = sorted(val, key=lambda x: ordered_matrix[left_matrix].index(x[0]))
            val_list = deque([x2 for x1, x2 in ordered_val])
            while len(val_list) < len(ordered_matrix):
                val_list.appendleft(0)
            ordered_output[left_matrix] = np.array(val_list)

        return ordered_output

def create_ndarray(ordered_output_dict, matrices):
    output_tuple = list(ordered_output_dict.items())
    ordered_tuple = sorted(output_tuple, key=lambda x: matrices.index(x[0]))
    val_arr = create_symmetric_matrix(np.stack([x2 for x1, x2 in ordered_tuple]))

    return val_arr

def orthogroups_analysis(datadir):

    print()
    print("--------------------------- Orthgroups -------------------------------")
    print()

    orthogroup_file_dict, _, _ = obtain_ogs_datafiles(datadir)

    orthogroup_sets_dict = {}
    with multiprocessing.Pool() as pool1:
        items1 = [( matrix, file) for matrix, file in orthogroup_file_dict.items()]
        results1 = pool1.starmap(create_ogs_sets, items1)
        for i in range(len(results1)):
            if results1[i]:
                orthogroup_sets_dict[items1[i][0]] = results1[i]

    matrices = [*orthogroup_sets_dict.keys()]
    matrice_with_ogs_num = []
    items2 = []

    track_left_matrix = {matrices[0]: []}

    for i in range(len(matrices)):
        left_set = orthogroup_sets_dict[matrices[i]]
        matrice_with_ogs_num.append("-".join((matrices[i], str(len(left_set)))))
        
        for j in range(i+1, len(matrices)):
            right_set = orthogroup_sets_dict[matrices[j]]
            items2.append((left_set, right_set, matrices[i], matrices[j]))

    comparison_stats = {}

    with multiprocessing.Pool() as pool2:
        results2 = pool2.starmap(orthogroups_sets_compare, items2)

        intersection_num, left_diff_num, right_diff_num, similarity, \
        expected, left_matrix, right_matrix = results2[0]

        comparison_stats["intersection"] = intersection_num
        comparison_stats["left_diff"] = left_diff_num
        comparison_stats["right_diff"] = right_diff_num
        comparison_stats["similarity"] = similarity
        comparison_stats["expectation"] = expected
        comparison_stats["labels"] = (left_matrix, right_matrix)

    return  comparison_stats


def unassigned_gens_analysis(datadir):    

    print()
    print("--------------------- Unassigned gens -------------------------------")
    print()
    
    _, unassigned_gens_file_dict, _ = obtain_ogs_datafiles(datadir)
    unassgined_gens_matrix_sets = create_unassgined_gens_matrix_sets(unassigned_gens_file_dict)
    matrices = [*unassgined_gens_matrix_sets.keys()]
    unassigned_intersection_left, unassigned_intersection_right, \
    unassigned_intersection_union = gen_sets_analysis(matrices, unassgined_gens_matrix_sets)
    

def orthologues_analysis(datadir):


    print()
    print("--------------------- Orthologues -------------------------------")
    print()
    _, _, orthologue_file_dict = obtain_ogs_datafiles(datadir)
    items3 = []
    for matrix, files in orthologue_file_dict.items():
        if len(files) != 0:
            for file in files:
                items3.append((matrix, file))
        else:
            items3.append((matrix, None))

    orthologue_sets_dict = {}
    with multiprocessing.Pool() as pool3:
        
        results3 = pool3.starmap(create_ogs_sets, items3)
        for i in range(len(results3)):
            if results3[i]:
                orthologue_sets_dict[items3[i][0]] = results3[i]
            else:
                orthologue_sets_dict[items3[i][0]] = set()

    matrices = [*orthologue_sets_dict.keys()]

    orthologue_intersection_left, orthologue_intersection_right, \
    orthologue_intersection_union = gen_sets_analysis(matrices, orthologue_sets_dict)
 

if __name__ == "__main__":
    
    args = sys.argv[1:]
    datadir = "ExampleData0" if len(args) == 0 else args[-1]    
    comparison_stats = orthogroups_analysis(datadir)
    # unassigned_gens_analysis(datadir)
    # orthologues_analysis(datadir)
