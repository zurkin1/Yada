import numpy as np
import pandas as pd
from scipy.special import comb
#from tslearn import metrics


# Given a mixture and a sorted list of genes, check if expression levels are in the same order.
# Provide a score for each gene.
def spearman_rank(set_1, set_2):
    # order the sets
    set_1_ord = sorted(set_1)
    set_2_ord = sorted(set_2)

    # append relevant rank to each value in set
    set_1_ranked = []
    set_2_ranked = []
    d = []

    for i in range(len(set_1)):
        set_1_ranked.append([set_1[i], set_1_ord.index(set_1[i]) + 1])
        set_2_ranked.append([set_2[i], set_2_ord.index(set_2[i]) + 1])
        d.append(set_1_ranked[i][1] - set_2_ranked[i][1])

    d_sq = [i**2 for i in d]
    # Indices of smallest elements.
    indices_small = np.argsort(d_sq)[0:round(0.9 * len(set_2))]
    sum_d_sq = sum(d_sq)
    # calculate n^3 - n
    n_cu_min_n = len(set_1)**3 - len(set_1)
    r = 1 - ((6.0 * sum_d_sq) / n_cu_min_n)  # Spearman rank.
    return indices_small, round(r, 2)


# Set_1 should be the pure expression values.
def kendall_rank(set_1, set_2):
    dist = np.zeros(len(set_1))
    for i in range(len(set_1)):
        for j in range(i, len(set_1)):
            if ((set_1[i] <= set_1[j]) and (set_2[i] > set_2[j])) or ((set_1[i] > set_1[j]) and (set_2[i] <= set_2[j])):
                dist[j] += 1
                dist[i] += 1
    indices_small = np.argsort(dist)[0:round(0.8 * len(set_2))]
    rank = dist.sum() / comb(len(set_1), 2, exact=True)
    return indices_small, round(rank, 2)


def dtw_rank(set_1, set_2):
    dist = np.zeros(len(set_1))
    for i in range(len(set_1)):
        set_1_temp = np.delete(set_1.values, i)
        set_2_temp = np.delete(set_2, i)
        P = np.array([set_1_temp, np.arange(0, len(set_1_temp))])
        Q = np.array([set_1_temp, np.arange(0, len(set_2_temp))])
        dist[i] = metrics.dtw(P, Q)
    indices_large = np.argsort(dist)[len(dist):-round(0.9 * len(dist)) - 1:-1]
    return indices_large, 0


# This function does a gene differentiation analysis. It finds for each cell types its majorodizing genes.
# We also get a factor of how much (strongly majordizing) each gene is.
def gene_diff(pure, mix):
    num_genes_per_cell = 80
    gene_list_per_cell = {}
    # Loop on all cells
    for cell_type in pure:
        gene_list = []
        data = pure.copy()
        for factor in range(100, 1, -1):  # 2
            data2 = data.copy()
            data2['max_val'] = data.drop(cell_type, axis=1).max(axis=1)
            data2['med_val'] = data.drop(cell_type, axis=1).median(axis=1)
            data2 = data2.loc[(data2[cell_type] > 0.5 * factor * data2.max_val) & (data2.max_val > 0)]  # Ignore cases where all other mixes are zeros.
            # data2 = data2.loc[(data2[cell_type] > 0.5 * factor * data2.max_val) & (data2.med_val > 0)]  # Ignore cases where all other mixes are zeros.
            if (len(data2) == 0):
                continue
            data2.sort_values(cell_type, inplace=True, ascending=False)
            #gene_list += list(zip(data2.index, [factor for i in range(len(data2.index))]))
            gene_list += list(data2.index)
            data.drop(data2.index, inplace=True)
            if (len(gene_list) > num_genes_per_cell):
                gene_list = gene_list[0:num_genes_per_cell]
                break

        # Select maximum correlated genes (max(mix columns) agrees with pure).
        #max_col = mix.loc[gene_list].sum().idxmax()
        #mix_column = mix[max_col].loc[gene_list]
        #indices_small, score = kendall_rank(pure[cell_type].loc[gene_list], mix_column.values)
        #gene_list = [gene_list[i] for i in indices_small]

        #print(f'{cell_type} : {gene_list}')
        if (len(gene_list) < 8):
            print(f'WARNING: {cell_type} short differential gene list, length {len(gene_list)}.')
        gene_list_per_cell[cell_type] = gene_list

    gene_list_df = pd.DataFrame(dict([(k, pd.Series([i for i in v])) for k, v in gene_list_per_cell.items()]), columns=pure.columns)
    #factor_list = pd.DataFrame(dict([(k, pd.Series([j * 0.5 for (i, j) in v])) for k, v in gene_list_per_cell.items()]), columns=pure.columns)

    #import IPython
    # IPython.embed()
    # for cell in cells:
    #  print(cell, ',', *[x for x in gene_list_df[cell].values.tolist() if str(x) != 'nan'])
    return gene_list_df
