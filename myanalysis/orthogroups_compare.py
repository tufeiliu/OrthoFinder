import pathlib
import os 
import pandas as pd 
import numpy as np 
from collections import deque 
import seaborn as sns
import matplotlib.pyplot as plt
import functools
from multiprocessing.pool import Pool


# colnames = ["BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70", "PAM250"]
colnames = ["BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM30", "PAM70", "PAM250"]
print(colnames )

orthofinder_outputs = pathlib.Path.cwd() / "ExampleData" / "OrthoFinder"
species_fasta = pathlib.Path.cwd() / "ExampleData"
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



def get_intersection(gdf):

    intersection_per_dict = {}

    for i in range(gdf.shape[1]):
        coli = gdf.columns[i]
        intersection_per_dict[coli] = deque()

        for j in range(i, gdf.shape[1]):
            colj = gdf.columns[j]
            
            # intersection_num_left = gdf[coli].isin(gdf[colj]).value_counts()

            # print(f"Number of {coli} in {colj}: {intersection_num}")

            intersection_per = gdf[coli].isin(gdf[colj]).value_counts(normalize=True) * 100
            # print(f"Percentage of {coli} in {colj}: {intersection_per}")
            # intersection_per
            intersection_per_dict[coli].append(intersection_per[True])
    
    extended_ins_per_dict = {}
    for key, val in intersection_per_dict.items():
        
        while len(val) < gdf.shape[1]:
            val.appendleft(0.0)
        extended_ins_per_dict[key] = np.array(val)
    
    
    cols = extended_ins_per_dict.keys()
    vals = np.stack(list(extended_ins_per_dict.values()))
    symmetric_val = vals + vals.T - np.diag(vals.diagonal())
    
    return symmetric_val
    

def unassigned_comparison(gens_file_dict, species_gens):

    for matrix_name, file in gens_file_dict.items():
        # print(f"***************** {matrix_name} **********************")
        df =  pd.read_csv(file, sep="\t")

        for col in df.columns[1:]:
            
            gens = df[col].dropna().sort_values(key=lambda x: x.str.lower()).to_frame(name=matrix_name)
            gens.reset_index(drop=True, inplace=True)
        
            if matrix_name not in species_gens[col]:
                species_gens[col][matrix_name] = gens

    for species in species_gens.keys():
        # print(species)
        species_method_df = pd.concat(species_gens[species].values(),
                                        join="outer", axis=1)
        species_method_df = species_method_df[colnames]
        # print(species_method_df.info())
        val_matrix = get_intersection(species_method_df)
        plot_heatmap(val_matrix, colnames, species, gc="unassigned")


def plot_heatmap(symmetric_matrix, cols, species, gc="unassined"):
    symmetric_df = pd.DataFrame(symmetric_matrix, index=cols, columns=cols)
    plt.figure(figsize=(16, 12))
    ax = sns.heatmap(symmetric_df, annot=True, fmt=".0f", annot_kws={"size": 15}, cmap="viridis")
    ax.tick_params(axis='both', which='both', labelsize=15)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    cbar = ax.collections[0].colorbar
    
    cbar.formatter.set_powerlimits((0, 0))
    cbar.ax.tick_params(labelsize=15)
    cbar.update_ticks()
    plt.savefig(f'{species}_{gc}.png')

def merge_dfs(dfs):
    merged_df = functools.reduce(lambda df1, df2: 
                                 pd.merge(df1, df2, how="outer", on="Orthogroup"),
                                 dfs)
    return merged_df
    

def get_species_methods_dict(gens_file_dict, species_gens):

    species_dfs = {species: [] for species, _ in species_gens.items()}

    for matrix_name, file in gens_file_dict.items():
        # print(f"***************** {matrix_name} **********************")
    
        df = pd.read_csv(file, sep="\t")

        for col in df.columns[1:]:
            gens_df = df[["Orthogroup", col]].dropna() \
                            .map(lambda x: x.split(",") if "," in x else x)

            gens_expanded = gens_df.explode(col)
            gens_expanded.sort_values(by=["Orthogroup", col], ignore_index=True, inplace=True)
            gens_expanded.reset_index(drop=True, inplace=True)
            gens_expanded.rename(columns={col: matrix_name}, inplace=True)

            # gens = gens_expanded.groupby("Orthogroup")[matrix_name].count()         
            # if matrix_name not in species_dfs[col]:
                # species_gens[col][matrix_name] = gens_expanded
            species_dfs[col].append(gens_expanded)
    
    
    # species_method_df = {}
    # for species, _ in species_gens.items():
    #     species_method_df[species] = list(species_gens[species].values()) #
    # print(species_method_df)
    return species_dfs


def compare_ortho(df1, df2):

    colnames = [df1.columns[1], df2.columns[1], 
                df1.columns[1] + "_num", df2.columns[1] + "_num",
                "intersection", df1.columns[1] + "-" + df2.columns[1], 
                "intersection/" + df1.columns[1], "intersection/" + df2.columns[1], 
                "intersection/(" + df1.columns[1] + "+" + df2.columns[1] + ")",
                "expected"]
    
    similarity_freq = {}
    similarity_freq_dfs = {}
    for labell, dfl in df1.groupby("Orthogroup"):
        similarity_freq[labell] = []
        for labelr, dfr in df2.groupby("Orthogroup"):

            df1_in_df2 = dfl.iloc[:, 1].isin(dfr.iloc[:, 1]).value_counts()
            
            if len(df1_in_df2.index) > 1:
                df1_in_df2_num = df1_in_df2[True]
                df1_notin_df2_num = df1_in_df2[False]
                df1_in_df2_per_df1 = df1_in_df2_num / dfl.shape[0]
                df1_in_df2_per_df2 = df1_in_df2_num / dfr.shape[0]
                df1_in_df2_per_df12 = df1_in_df2_num / (dfl.shape[0] + dfr.shape[0] - df1_in_df2_num)
                df1_in_df2_expected = df1_in_df2_per_df12 * dfl.shape[0]
                similarity_freq[labell].append((
                                        labell,
                                        labelr,
                                        dfl.shape[0],
                                        dfr.shape[0],
                                        df1_in_df2_num, 
                                        df1_notin_df2_num,
                                        df1_in_df2_per_df1,
                                        df1_in_df2_per_df2,
                                        df1_in_df2_per_df12,
                                        df1_in_df2_expected
                                        ))

            elif df1_in_df2.index == True:
                df1_in_df2_num = df1_in_df2[True]
                df1_notin_df2_num = 0
                df1_in_df2_per_df1 = df1_in_df2_num / dfl.shape[0]
                df1_in_df2_per_df2 = df1_in_df2_num / dfr.shape[0]
                df1_in_df2_per_df12 = df1_in_df2_num / (dfl.shape[0] + dfr.shape[0] - df1_in_df2_num)
                df1_in_df2_expected = df1_in_df2_per_df12 * (dfl.shape[0] + dfr.shape[0])
                similarity_freq[labell].append((
                                        labell,
                                        labelr,
                                        dfl.shape[0],
                                        dfr.shape[0],
                                        df1_in_df2_num, 
                                        df1_notin_df2_num,
                                        df1_in_df2_per_df1,
                                        df1_in_df2_per_df2,
                                        df1_in_df2_per_df12,
                                        df1_in_df2_expected
                                        ))

            elif df1_in_df2.index == False:

                df1_in_df2_num = 0
                df1_notin_df2_num = df1_in_df2[False]
                df1_in_df2_per_df1 = 0.0
                df1_in_df2_per_df2 = 0.0
                df1_in_df2_per_df12 = 0.0
                df1_in_df2_expected = 0.0
                similarity_freq[labell].append((
                                        labell,
                                        labelr,
                                        dfl.shape[0],
                                        dfr.shape[0],
                                        df1_in_df2_num, 
                                        df1_notin_df2_num,
                                        df1_in_df2_per_df1,
                                        df1_in_df2_per_df2,
                                        df1_in_df2_per_df12,
                                        df1_in_df2_expected
                                        ))

        similarity_freq_df = pd.DataFrame(similarity_freq[labell], 
                                          columns=colnames
                                        #   index=df1.groupby("Orthogroup").groups.keys(),
                                        )

        similarity_freq_dfs[labell] = similarity_freq_df

    return similarity_freq_dfs

def orthogroup_comparison(species_ortho_dict, dfs):
    
    for n in range(len(dfs)):

        df1 = dfs[n]
        for m in range(n+1, len(dfs)):
            df2 = dfs[m]
            similarity_freq1 = compare_ortho(df1, df2)
            similarity_freq2 = compare_ortho(df2, df1)
            species_ortho_dict[species].append((similarity_freq1, similarity_freq2))
    

    return species_ortho_dict   

def ortho_comm(freqs, criteria="expected"):
    
    ortho_max = {}

    for label, df in freqs.items():
        if not df.empty:
            df_max = df[df[criteria] == df[criteria].max()].copy()
            df_max.drop_duplicates(subset=df_max.columns[2:], inplace=True)
            if df_max.shape[0] > 1:
                df_max.drop_duplicates(subset=df_max.columns[0], inplace=True)
            
            ortho_max[label] = df_max
            
    if len(ortho_max) != 0:
        ortho_max_df = pd.concat(ortho_max.values(), ignore_index=True)
        ortho_intersec = ortho_max_df[ortho_max_df.intersection != 0].copy()
        print(ortho_intersec)
        return ortho_intersec
    else:
        return pd.DataFrame()

# unassigned_comparison(unassigned_gens_file_dict, species_gens)
# species_gens_dict = get_species_methods_dict(unassigned_gens_file_dict, species_gens)
species_gens_dict = get_species_methods_dict(orthogroup_file_dict, species_gens)
# print(species_gens_dict)
ortho_similarity = {}
for species, dfs in species_gens_dict.items():
    species_ortho_dict = {species: []}
    species_ortho_dict = orthogroup_comparison(species_ortho_dict, dfs)
    print(f"****************************** {species} *****************************")
    for smf1, smf2 in species_ortho_dict[species]:
        smf1_filtered = ortho_comm(smf1)
        smf2_filtered = ortho_comm(smf2)
        # if not smf1_filtered.empty or smf2_filtered:
        ortho_similarity[species] = (smf1_filtered, smf2_filtered)

for species, vals in ortho_similarity.items():
    print(f"****************************** {species} *****************************")
    if not vals[0].empty:
        smf1_total_similarity = vals[0].expected.sum()
        print(smf1_total_similarity)
    if not vals[1].empty:
        smf2_total_similarity = vals[1].expected.sum()
        print(smf2_total_similarity)


# def custom_callback(result_iterable):
# 	for result in result_iterable:
# 		print(f'Got result: {result}', flush=True)
 
# def custom_error_callback(error):
# 	print(f'Got error: {error}')

# if __name__ == '__main__':
#     items = [(species, dfs) for species, dfs in species_gens_dict.items()]
#     with Pool() as pool:
#         result = pool.starmap(orthogroup_comparison, 
#                                 items,
#                                 # chunksize=100,
#                                 # callback=custom_callback,
#                                 # error_callback=custom_error_callback
#                                 )
#         # print(result.get())                 
#         # iterate results
#         for result in result.get():
#             print(f'Got result: {result}', flush=True)


# if __name__ == "__main__":
    # 
             
         

        
        

