import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from collections import deque 
import itertools 
import os
import numpy.typing as npt
import scipy as spy
import seaborn as sns
import json
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_samples, silhouette_score


def str_compare(s1: str, s2: str):
    diff_list = []
    non_aa = "XOJUZB-"
    for i, c in enumerate(itertools.zip_longest(s1, s2)):
        if (c[0] != c[1]) or (c[0] in non_aa) or (c[1] in non_aa):
            diff_list.append((i, c))
    return diff_list


def write_to_file(filename: str, rownames: str, colnames: str, arr: npt.NDArray[np.float64]):
    

    BZX_row = np.array([[-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, 
                     -3, -3, -2, 0, -1, -4, -3, -3,  4,  1, -1],
                    [-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, 
                     -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1],
                    [0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1,
                     -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1]])


    BZX_col = np.array([[-2, -1, 0],
                        [-1, 0, -1],
                        [3, 0, -1],
                        [4, 1, -1],
                        [-3, -3, -2],
                        [0, 3, -1],
                        [1, 4, -1],
                        [-1, -2, -1],
                        [0,  0, -1],
                        [-3, -3, -1],
                        [-4, -3, -1],
                        [0, 1, -1],
                        [-3, -1, -1],
                        [-3, -3, -1],
                        [-2, -1, -2],
                        [0, 0, 0],
                        [-1, -1, 0],
                        [-4, -3, -2],
                        [-3, -2, -1],
                        [-3, -2, -1],
                    ])

    # row_index = list(rownames + "BZX")
    # col_index = list(colnames + "BZX")

    row_index = list(rownames)
    col_index = list(colnames)
    
    arr = np.round(arr).astype(int)
    # arr = np.hstack((arr, BZX_col))
    # arr = np.vstack((arr, BZX_row))
    df = pd.DataFrame(arr, index=row_index, columns=col_index)

    filedir = os.path.join(os.getcwd(), "substitution_matrix")
    if not os.path.exists(filedir):
        os.makedirs(filedir, exist_ok=False)

    filepath = os.path.join(filedir, filename)
    df.to_csv(filepath, sep='\t')
    # with open(filepath, 'a') as f:
    #     f.write(df.to_string())
    


def spearman_corr(data_arr: npt.NDArray[np.float64]):

    corr_matrix, p_matrix = spy.stats.spearmanr(data_arr, axis=1)
    corr_dist = np.sqrt(2 - 2 * corr_matrix)

    return corr_matrix, corr_dist

def pearson_corr(sdata_arr: npt.NDArray[np.float64]):
    
    corr_matrix = np.corrcoef(sdata_arr)
    corr_dist = np.sqrt(2 - 2 * corr_matrix)

    return corr_matrix, corr_dist


def fancy_dendrogram(*args, **kwargs):

    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = spy.cluster.hierarchy.dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):

        plt.title('Hierarchical Clustering Dendrogram')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')

        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')

    return ddata


def hierarchy_clustering(distm: npt.NDArray[np.float64], labels: list[str]):

    # dist_condensed = distm[np.tril_indices(distm.shape[0], -1)]
    metrics = ['euclidean', 'minkowski', 'cityblock', 
               'seuclidean', 'sqeuclidean', 'cosine',
               'correlation', 'chebyshev', 'canberra']
    methods = ["single", "complete", "ward", "average", "weighted", "centroid", "median"]
    
    max_c = 0
    print("Cophenetic correlation coefficients:")
    for metric in metrics:
        dist_condensed = spy.spatial.distance.pdist(distm, metric=metric)
        for method in methods:
            Z = spy.cluster.hierarchy.linkage(dist_condensed, method=method)
            c, coph_dists = spy.cluster.hierarchy.cophenet(Z, dist_condensed)
            print(f"{metric} - {method}: {c}")
            if c > max_c:
                print(f"*** Best method: {metric} - {method} ***", end="\n"*2)
                max_c = c 
                max_Z = Z

    return max_Z



def silhouette_scores(X: npt.NDArray[np.float64], 
                      cluster_labels: npt.NDArray[np.int_]):
    
    silhouette_avg = silhouette_score(X, cluster_labels)
    n_clusters = len(set(cluster_labels))
    print(
        "For n_clusters =",
        n_clusters,
        "The average silhouette_score is :",
        silhouette_avg,
    )
    
    # Compute the silhouette scores for each sample
    sample_silhouette_values = silhouette_samples(X, cluster_labels)

    fig, ax = plt.subplots()
    # fig.set_size_inches(18, 7)

    # The 1st subplot is the silhouette plot
    # The silhouette coefficient can range from -1, 1 but in this example all
    # lie within [-0.1, 1]
    ax.set_xlim([-0.1, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, len(X) + (n_clusters + 1) * 10])

    y_lower = 10
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / n_clusters)
        ax.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([-0.1, 0, 0.2, 0.4, 0.6, 0.8, 1])
    plt.show()



# def get_cluster_classes(den, label='ivl'):
#     cluster_idxs = defaultdict(list)
#     for c, pi in zip(den['color_list'], den['icoord']):
#         for leg in pi[1:3]:
#             i = (leg - 5.0) / 10.0
#             if abs(i - int(i)) < 1e-5:
#                 cluster_idxs[c].append(int(i))

#     cluster_classes = {}
#     for c, l in cluster_idxs.items():
#         i_l = [den[label][i] for i in l]
#         cluster_classes[c] = i_l

#     return cluster_classes



def clustermap(smatrix: npt.NDArray[np.float64], labels: list[str]):
    
    df =  pd.DataFrame(smatrix, index=labels, columns=labels)
    g = sns.clustermap(df, cmap=sns.cubehelix_palette(as_cmap=True))
    plt.show()

    row_linkage = spy.cluster.hierarchy.linkage(
        spy.spatial.distance.pdist(smatrix), method='average')

    col_linkage = spy.cluster.hierarchy.linkage(
        spy.spatial.distance.pdist(smatrix.T), method='average')

    sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, method="average")
    plt.show()

    den = spy.cluster.hierarchy.dendrogram(row_linkage,
                                         labels = labels,
                                         leaf_font_size=8.,
                                         color_threshold=0.60) 
    plt.show()

    return row_linkage



with open("aaindex2.txt") as reader:
    mname = {}
    id_list = deque()
    row_label_list = deque()
    col_label_list = deque()
    data_list = []
    data_lists = deque()
    missing_val_matrix = []

    for line in reader:
        if "H " in line:
            id_list.append(line.rstrip("\n"))
            data_lists.append(data_list)

        if "M rows =" in line:
            aa_list = line.rstrip("\n").split("=")
            rows = aa_list[1].split(",")[0].lstrip()
            cols = aa_list[2].lstrip()
            ncols = len(cols)
            row_label_list.append(rows)
            col_label_list.append(cols)
            row_vals = data_lists.popleft()
            if len(row_vals) != 0: 
                prev_rows = row_label_list.popleft()
                prev_cols = col_label_list.popleft()
                idf = id_list.popleft() + f" rows = {prev_rows} cols = {prev_cols}"
                mname[idf] = row_vals

            data_list = []

        try:
            row_val_list = line.rstrip("\n").split()
            row_vals = []
            for v in row_val_list:
                if v.strip() == "-":
                    print(f"Missing value in {id_list[-1]}")
                    v = 0.0
                    missing_val_matrix.append(id_list[-1].split()[1])
                row_vals.append(float(v))

            if len(row_vals) < ncols:
                zeros_str = "0.0 " * (ncols - len(row_vals))
                zeros = [float(z) for z in zeros_str.rstrip().split()]
                row_vals.extend(zeros)

            data_list.append(row_vals)
            if len(row_vals) == 0:
                data_lists[0] = data_list
            
        except:
            continue

BLOSUM = {"HENS920101": "BLOSUM45",
          "HENS920102": "BLOSUM62",
          "HENS920103": "BLOSUM80",
          "HENS920104": "BLOSUM50",
}

sub_matrix_val_list = []
ssub_matrix_val_list = [] # standarlised
sub_matrix_ids = []
for name, vals in mname.items():
    # print(name)
    matrix_id = name.split()[1]
    val_array = np.array(vals)
    nrows = val_array.shape[0]
    ncols = val_array.shape[1]
    rownames = name.split()[4]
    colnames = name.split()[-1]
    nrownames = len(rownames)
    ncolnames = len(colnames)
    
    diff_list = str_compare(rownames, colnames)

    if nrownames != 20 or ncolnames != 20:
        count = 0
        for i, (r, c) in diff_list:
            if count != 0:
                i -= 1
                
            if r and i < val_array.shape[0]:
                val_array = np.delete(val_array, (i), axis=0)
                rownames = rownames.replace(r, "")
                
            if c and i < val_array.shape[1]:
                val_array = np.delete(val_array, (i), axis=1)
                colnames = colnames.replace(c, "")
            
            count += 1
    

    if not np.allclose(val_array, val_array.T):
        if matrix_id != "GIAG010101":
            val_array = val_array + val_array.T - np.diag(np.diag(val_array))
        
    # print(name)
    # print(nrows, ncols)
    # print(val_array)


    if matrix_id not in missing_val_matrix:
        filename = ".".join((matrix_id, "txt"))
        
        write_to_file(filename, rownames, colnames, val_array)
    
        # val_tril_index = np.tril_indices(val_array.shape[0])
        # lower_tri = val_array[val_tril_index]
        # val_vector = np.ravel(lower_tri)
        val_vector = np.ravel(val_array)
        sub_matrix_val_list.append(val_vector)
        ssub_matrix_val_list.append((val_vector - np.mean(val_vector)) / np.std(val_vector, ddof=1))
        sub_matrix_ids.append(matrix_id)


sub_val_arr = np.stack(sub_matrix_val_list)
ssub_val_arr = np.stack(ssub_matrix_val_list)

df = pd.DataFrame(sub_val_arr, index=sub_matrix_ids)
sdf = pd.DataFrame(ssub_val_arr, index=sub_matrix_ids)

# -------- Spearman --------
spearman_corr_matrix, spearman_corr_dist = spearman_corr(sub_val_arr)
Z = hierarchy_clustering(spearman_corr_matrix, sub_matrix_ids)
sns_Z = clustermap(spearman_corr_matrix, sub_matrix_ids)

fancy_dendrogram(
    Z,
    labels=sub_matrix_ids,
    # truncate_mode='lastp',
    # p=12,
    # orientation="left", 
    leaf_rotation=90.,
    leaf_font_size=8.,
    # show_contracted=True,
    annotate_above=0.8,  # useful in small plots so annotations don't overlap
)
plt.show()

sns_cluster_labels = spy.cluster.hierarchy.fcluster(sns_Z, 1.5, depth=3)
cluster_labels_dict = {str(key): [] for key in np.unique(sns_cluster_labels)}

for k, l in zip(sns_cluster_labels, sub_matrix_ids):
    k = str(k)
    if k in cluster_labels_dict:
        cluster_labels_dict[k].append(l)

print(json.dumps(cluster_labels_dict, sort_keys=True, indent=4))


# silhouette_scores(sub_val_arr, cluster_labels)

# -------- Pearson --------
# pearson_corr_matrix, pearson_corr_dist = pearson_corr(ssub_val_arr)
# Z = hierarchy_clustering(pearson_corr_matrix, sub_matrix_ids)
# sns_Z =clustermap(pearson_corr_matrix, sub_matrix_ids)
# fancy_dendrogram(
#     Z,
#     labels=sub_matrix_ids,
#     leaf_rotation=90.,
#     leaf_font_size=8.,
#     annotate_above=0.5,  # useful in small plots so annotations don't overlap
# )
# plt.show()



## ---- plot cluster map -----
# corr_matrix_df =  pd.DataFrame(corr_matrix, index=sub_matrix_ids, columns=sub_matrix_ids)
# sns.clustermap(corr_matrix_df, cmap=sns.cubehelix_palette(as_cmap=True))
# plt.show()


# corr_dist = np.sqrt(2 - 2 * corr_matrix)
# # st = True if np.allclose(corr_dist, corr_dist.T) else False
# # corr_dist_condensed = spy.spatial.distance.squareform(corr_dist[:5, :5])

# corr_dist_condensed = corr_dist[np.tril_indices(corr_dist.shape[0], -1)]

# methods = ["single", "complete", "ward", "average", "weighted", "centroid", "median"]

# max_c = 0
# print("Cophenetic correlation coefficients:")
# for method in methods:
#     Z = spy.cluster.hierarchy.linkage(corr_dist_condensed, method=method)
#     c, coph_dists = spy.cluster.hierarchy.cophenet(Z, corr_dist_condensed)
#     print(f"{method}: {c}")
#     if c > max_c:
#         max_c = c 
#         max_Z = Z

# fig, ax = plt.subplots()
# dn = spy.cluster.hierarchy.dendrogram(max_Z, labels=sub_matrix_ids,
#                                     #   orientation="left", 
#                                       truncate_mode='lastp',  # show only the last p merged clusters
#                                       p=20,  # show only the last p merged clusters
#                                       show_contracted=True,  # to get a distribution impression in truncated branches
#                                      )
# ax.set_title('Hierarchical Clustering Dendrogram')
# ax.set_ylabel('sample index')
# ax.yaxis.set_label_position("right")
# ax.yaxis.tick_right()
# ax.set_xlabel('distance')
# plt.show()

# # move ticks
# plt.tick_params(axis='y', which='both', labelleft=False, labelright=True)

# # move label
# plt.ylabel('Your label here', labelpad=-725, fontsize=18)







            

        


            

        
