import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
from scipy.optimize import nnls
from tslearn import metrics
from diffexp import gene_diff
from dataLoader import load_data
from dsection import dsection

pure = pd.DataFrame()
mix = pd.DataFrame()
gene_list_df = pd.DataFrame()
factor_list = pd.DataFrame()
real_weight = pd.DataFrame()
mixtures_sum = []
other_res = 100
num_mixes = 0
num_cells = 0
test_cases = ['RatBrain', 'Abbas', 'TIMER', 'DSA', 'DeconRNASeq', 'BreastBlood', '10x', 'EPIC', 'CIBERSORT', 'PertU']

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

        #Loop on all mixes.
        k = 0
        for mix_col in mix_temp:
            column = mix_temp[mix_col].copy()
            column = column[max_column.index] #Sort according to the maximum column index.

            P = np.array([column, np.arange(0,len(column.index))])
            Q = np.array([max_column, np.arange(0, len(max_column.index))])
            dist = metrics.dtw(Q,P) #0.342
            if (len(cell_genelist) <= 5):
                dist = max_column.mean() - column.mean()
                #print('.', end="")
            if(useSimpleDist):
                dist = 0
                j = 0
                for gene in cell_genelist:
                    # Subtract the average of other cells on this gene to compensate for over estimation of ratio.
                    dist += mixtures_sum[k] * \
                            ((column[gene] + 0.000001)) / (max_column[gene] + 0.000001) * \
                            ((factor_list[cell_type][j] + 0.000001) / ((factor_list[cell_type].sum() + 0.000001)))
                    j += 1

            cell_vals.append(dist)
            k += 1
        O_array[i] = cell_vals
        i += 1
    #We have to scale to [0,1] in order for the sum of proportions to be meaningful. Scale along each cell.
    if(useSimpleDist):
        O_scaled = np.transpose(O_array)
    else:
        O_scaled = np.transpose(1 - minmax_scale(O_array, axis=1))
    solution = nnls(O_scaled, mixtures_sum)[0]
    solution_mat = np.diag(solution)
    estimate_wt = np.matmul(solution_mat, O_scaled.T)
    return(estimate_wt)

"""This method calculate the score (dist) of an algorithm using L1-norm of distances of proportions.
result: Calculate how far pure*proportions is from mix.
dist requires knowledge of the true values. But result does not."""
def calc_scores(a_prop, b_prop, method_i):
    dist_from_real_value = np.round(np.sum(np.linalg.norm(np.abs(a_prop - b_prop), axis=0)), 2)
    #all_marker_genes = pd.concat([gene_list_df.iloc[:, i] for i in range(len(gene_list_df.columns))]).reset_index(drop=True).dropna()
    #Deviding by mixtures_sum since guessing the right proportions on 0.7 mixtures is harder than on mixtures that are 1.
    mixtures_sum_quotient = [x**-2 for x in mixtures_sum]
    dist_from_mix = np.round(np.sum(np.dot(np.abs(mix - np.dot(pure,b_prop.T)), mixtures_sum_quotient)/mix.size), 2)
    if (other_res < 5):
      dist_from_mix = np.round(np.square(mix - np.dot(pure,b_prop.T)).values.mean(), 2)
    print(f'Method {method_i}, distance from real is: {dist_from_real_value}, distance from mix {dist_from_mix}')
    return dist_from_real_value, dist_from_mix

if __name__ == '__main__':
    for test in test_cases:
        mix, pure, real_weight, other_result, mixtures_sum = load_data(test)
        num_mixes = len(mix.columns)
        num_cells = len(pure.columns)
        print(f'\nTest case: {test}, num_cells: {num_cells}, num_mixes: {num_mixes}')
        result=[0,0,0]
        dist=[0,0,0]
        if(len(mixtures_sum) == 0):
            mixtures_sum=[1 for i in range(num_mixes)]
        gene_list_df, factor_list = gene_diff(pure)
        other_dist, _ = calc_scores(real_weight, other_result, 'CIBERSORT')
        other_res = np.round(np.sum(np.dot(np.abs(mix - np.dot(pure,other_result.T)), mixtures_sum)/mix.size), 2)

        for method in [[0, deconv, True, 10],
                       [1, deconv, False, 400],
                       [2, dsection, {'mix':mix, 'pure':pure, 'gene_list_df':gene_list_df}, 1]]: #[index, method, method parameters, number of ensembles]
            ens_estimate_wt = np.zeros((num_cells, num_mixes))
            for ens_i in range(method[3]):
                estimate_wt = method[1](method[2])
                ens_estimate_wt += estimate_wt
            ens_estimate_wt /= method[3]
            dist[method[0]], result[method[0]] = calc_scores(real_weight, ens_estimate_wt.T, method[0] + 1)

        ind = np.argmin([result[0], result[1], result[2]])
        print(f'Best result: method {ind + 1}.')
        assert dist[ind] <= other_dist, "Failed test " + test