import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
from scipy.optimize import nnls
from diffexp import gene_diff
import os

pure = pd.DataFrame()
mix = pd.DataFrame()
gene_list_df = pd.DataFrame()
factor_list = pd.DataFrame()
real_weight = pd.DataFrame()
mixtures_sum = []
other_res = 100
num_mixes = 0
num_cells = 0
test_cases = ['EPIC']

#This function calculate a deconvolution algorithm.
def deconv(useSimpleDist):
    global pure, mix, gene_list_df, factor_list, real_weight, mixtures_sum

    O_array = np.zeros((num_cells, num_mixes))  # The per cell sorted array - how far each mix is from the maximal mix.
    i = 0
    #Loop on all cell types.
    for cell_type in gene_list_df:
        cell_vals = []
        cell_genelist = gene_list_df[cell_type].dropna().sample(frac=0.35)
        mix_temp = mix.loc[cell_genelist]
        max_ind = mix_temp.sum().idxmax() #Mix with maximum sum of gene expression.
        max_column = mix_temp[max_ind].copy()
        max_column.sort_values(ascending=False, inplace=True)
        #print(f'Max ind: {max_ind} len(max_column): {len(max_column)}')

        #Loop on all mixes.
        k = 0
        for mix_col in mix_temp:
            column = mix_temp[mix_col].copy()
            column = column[max_column.index] #Sort according to the maximum column index.

            if(useSimpleDist):
                dist = 0
                j = 0
                len_cell_genes = len(cell_genelist)
                for gene in cell_genelist:
                    # Subtract the average of other cells on this gene to compensate for over estimation of ratio.
                    delta = ((column[gene] + 0.000001) / (max_column[gene] + 0.000001)) * \
                            ((factor_list[cell_type][j] + 0.000001) / (factor_list[cell_type].sum() + 0.000001))
                    #Handle missing values
                    if np.isnan(delta):
                        len_cell_genes -= 1
                    else:
                        dist += delta
                    j += 1
                dist = dist/len_cell_genes
            #print(f'Cell type: {cell_type}, mix col: {mix_col}, dist: {dist}')
            cell_vals.append(dist)
            k += 1
        O_array[i] = cell_vals
        i += 1
    #We have to scale to [0,1] in order for the sum of proportions to be meaningful. Scale along each cell.
    if(useSimpleDist):
        O_scaled = np.transpose(O_array)
    else:
        O_scaled = np.transpose(1 - minmax_scale(O_array, axis=1))
    #print(O_scaled)
    solution = nnls(O_scaled, mixtures_sum)[0]
    solution_mat = np.diag(solution)
    estimate_wt = np.matmul(solution_mat, O_scaled.T)
    return(estimate_wt)


if __name__ == '__main__':
    desc_file = pd.read_csv('/input/input.csv')
    #cells = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells'] #, 'neutrophils', 'monocytic.lineage', 'fibroblasts', 'endothelial.cells']
    k = 0
    output_df = pd.DataFrame(columns=['dataset.name', 'sample.id', 'cell.type', 'prediction'])

    for index, row in desc_file.iterrows():
        pure = pd.read_csv('/pure.csv', index_col=0)
        mix = pd.read_csv('/input/' + row['hugo.expr.file'], index_col=0)
        both_genes = list(set(mix.index) & set(pure.index))
        pure = pure.reindex(both_genes)  # Drop genes that don't appear in mix.
        mix = mix.reindex(both_genes)
        num_mixes = len(mix.columns)
        num_cells = len(pure.columns)
        print(f'Test case: {row["dataset.name"]}, num_cells: {num_cells}, num_mixes: {num_mixes}')
        result=0
        mixtures_sum=[1 for i in range(num_mixes)]
        gene_list_df, factor_list = gene_diff(pure)

        for method in [[0, deconv, True, 10]]:
            ens_estimate_wt = np.zeros((num_cells, num_mixes))
            for ens_i in range(method[3]):
                estimate_wt = method[1](method[2])
                ens_estimate_wt += estimate_wt
            ens_estimate_wt /= method[3]

        dist_from_mix = np.round(np.sum(np.sum(np.abs(mix - np.dot(pure, ens_estimate_wt)) / mix.size)), 2)
        print(f'Method distance from mix {dist_from_mix}')

        for i in range(num_cells): #Loop on cell-types.
            for j in range(num_mixes): #Loop on samples.
                output_df.loc[k] = [row['dataset.name'], mix.columns[j], pure.columns[i], ens_estimate_wt[i,j]]
                k += 1

    if not os.path.exists('/output'):
        os.makedirs('/output')
    output_df.to_csv('/output/predictions.csv', index=None)
    print(output_df)