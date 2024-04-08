import pathlib
import pandas as pd 
import numpy as np 
import os
from collections import deque 
import matplotlib.pyplot as plt
import seaborn as sns

orthofinder_outputs = pathlib.Path.cwd() / "ExampleData0" / "OrthoFinder"
species_fasta = pathlib.Path.cwd() / "ExampleData0"
species_gens = {file.name.split(".")[0]: {} for file in species_fasta.iterdir() if file.is_file() }

orthologue_file_dict = {}
for output_dir in orthofinder_outputs.iterdir():
    output_dir_name = output_dir.name 
    output_matrix_result = output_dir_name.split("_")[-1]
    orthologue_file_dict[output_matrix_result] = {}
    
    for category_dir in output_dir.iterdir():
        if category_dir.name == "Orthologues":
            for species_dir in category_dir.iterdir():
                if not species_dir.is_file():

                    left_species = species_dir.name.replace("Orthologues_", "")
                    orthologue_file_dict[output_matrix_result][left_species] = []
                    for file in os.listdir(species_dir):
                        orthologue_file_dict[output_matrix_result][left_species].append(species_dir / file)
                        right_species = file.split(".")[0].split("__v__")[1]
                        species_gens[left_species][right_species] = {}


def plot_heatmap(symmetric_matrix, cols, left_species, right_species, gc="unassined"):
    symmetric_df = pd.DataFrame(symmetric_matrix, index=cols, columns=cols)
    plt.figure(figsize=(20, 14))
    ax = sns.heatmap(symmetric_df, annot=True, fmt=".0f", annot_kws={"size": 15}, cmap="viridis")
    ax.tick_params(axis='both', which='both', labelsize=15)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    cbar = ax.collections[0].colorbar
    
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.tick_params(labelsize=15)
    cbar.update_ticks()
    plt.title(f"Orthologues intersection {gc} (%) ({left_species}, {right_species})", fontsize = 20)
    output_path = pathlib.Path.cwd() / "orthologues_compare" / f'{left_species}_vs_{right_species}_{gc}.png'
    plt.savefig(output_path)
    plt.close()


# for matrix_name in orthologue_file_dict.keys():
#     # print(f"********************* {matrix_name} ****************************")
#     for left_species, files in orthologue_file_dict[matrix_name].items():
        
#         for file in files:
#             right_species = file.name.split(".")[0].split("__v__")[1]
#             orthologue_df = pd.read_csv(file, sep="\t")
#             orthologue_df = orthologue_df. \
#                             map(lambda x: x.split(", ") if "," in x else x)

#             compared_species = orthologue_df.columns[1:]
#             for col in compared_species:
#                 orthologue_df = orthologue_df.explode(col)

#             orthologue_df.reset_index(drop=True, inplace=True)
#             orthologue_df["_v_".join(compared_species)] = orthologue_df.iloc[:, 1] + ", " + orthologue_df.iloc[:, 2]
#             orthologue_df.drop(columns=compared_species, axis=1, inplace=True)
            
#             if matrix_name not in species_gens[left_species][right_species]:
#                 species_gens[left_species][right_species][matrix_name] = orthologue_df


for matrix_name in orthologue_file_dict.keys():
    for left_species, files in orthologue_file_dict[matrix_name].items():
        
        for file in files:
            right_species = file.name.split(".")[0].split("__v__")[1]
            matrix_olg_sets = set()
            with open(file) as reader:
                for i, line in enumerate(reader):
                    line = line.split("\t")[1:]
                    line[-1] = line[-1].replace("\n", "")

                    if i > 0:
                        gens = set()

                        for l in line[0].split(", "):
                            if len(l) != 0:
                                gens.add((left_species, l))
                        for r in line[1].split(", "):
                            if len(r) != 0:
                                gens.add((right_species, l))
                        
                        gens = frozenset(gens)
                        matrix_olg_sets.add(gens)

            species_gens[left_species][right_species][matrix_name] =  matrix_olg_sets  

intersection_dict = {}
left_diff_dict = {}
right_diff_dict = {}
for left_species in species_gens:
    intersection_dict[left_species] = {}
    left_diff_dict[left_species] = {}
    right_diff_dict[left_species] = {}

    for right_species in species_gens[left_species]:
        matrices = [*species_gens[left_species][right_species].keys()]

        intersection_list0 = []
        left_diff_list0 = []
        right_diff_list0 = []

        for i in range(len(matrices)):
            left_set = species_gens[left_species][right_species][matrices[i]]

            intersection_list1 = deque() 
            left_diff_list1 = deque()
            right_diff_list1 = deque()

            for j in range(i+1, len(matrices)):
                right_set = species_gens[left_species][right_species][matrices[j]]

                intersection = len(left_set & right_set)
                left_diff = len(left_set - right_set)
                right_diff = len(right_set - left_set)

                intersection_list1.append(intersection)
                left_diff_list1.append(left_diff)
                right_diff_list1.append(right_diff)
        
            while len(intersection_list1) < len(matrices):
                intersection_list1.appendleft(0)
            
            while len(left_diff_list1) < len(matrices):
                left_diff_list1.appendleft(0)
            
            while len(right_diff_list1) < len(matrices):
                right_diff_list1.appendleft(0)
            

            intersection_list0.append(np.array(intersection_list1))
            left_diff_list0.append(np.array(left_diff_list1))
            right_diff_list0.append(np.array(right_diff_list1))
        
        intersection_arr = np.stack(intersection_list0)
        left_diff_arr = np.stack(left_diff_list0)
        right_diff_arr = np.stack(right_diff_list0)

        intersection_arr = intersection_arr + intersection_arr.T - np.diag(intersection_arr.diagonal())
        left_diff_arr = left_diff_arr + left_diff_arr.T - np.diag(left_diff_arr.diagonal())
        right_diff_arr = right_diff_arr + right_diff_arr.T - np.diag(right_diff_arr.diagonal())

        intersection_left = 100 * intersection_arr / (intersection_arr + left_diff_arr) 
        intersection_right = 100 * intersection_arr / (intersection_arr + right_diff_arr)
        intersection_union = 100 * intersection_arr / (intersection_arr + left_diff_arr + right_diff_arr)
        
        plot_heatmap(intersection_left, matrices, left_species, right_species, gc=f"divided by left species")
        plot_heatmap(intersection_right, matrices, left_species, right_species, gc=f"divided by right species")
        plot_heatmap(intersection_union, matrices, left_species, right_species, gc=f"divided by both")
                



                