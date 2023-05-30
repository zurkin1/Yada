import numpy as np
import pandas as pd
from diffexp import gene_diff
import scipy.stats as st
import time
import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
from scipy.optimize import nnls, minimize
from scipy import stats
from scipy.spatial import minkowski_distance
from tslearn import metrics
import logging
import warnings
import random
from random import choice
#mport similaritymeasures
#from similaritymeasures import pcm
#from tqdm import tqdm
import multiprocessing as mp


warnings.filterwarnings("ignore")
logger = logging.getLogger("pymc3")
logger.propagate = False
logger.setLevel(logging.ERROR)


# This function calculate a deconvolution algorithm using basic ratios algorithm.
def basic_deconv(P, Q):
    dist = 0
    for index, value in P.items():
        # Subtract the average of other cells on this gene to compensate for over estimation of ratio.
        dist += ((value + 0.000001) / (Q.loc[index] + 0.000001))  # * \
        #((factor_list[cell_type][j] + 0.000001) / (factor_list[cell_type].sum() + 0.000001))
    return(dist)


def dtw_metric(P, Q):
    #fns = [metrics.dtw] #, similaritymeasures.dtw]
    #P = P.sample(frac = 1)
    #Q = Q[P.index]
    P1 = np.array([P.values, np.linspace(0, np.max(P), len(P))]) # np.arange(0, len(P))
    Q1 = np.array([Q.values, np.linspace(0, np.max(P), len(Q))])
    if (len(P) <= 5):
        return(P.mean() - Q.mean())
    #factor = max(0.01, np.corrcoef(P, Q)[0][1])
    #return similaritymeasures.dtw(Q1.T,P1.T)[0] #0.342
    return(metrics.dtw(P1, Q1))  # 0.342
    #return(pcm(P1, Q1))
    # return choice(fns)(P1, Q1)


# This function calculate a deconvolution algorithm using DTW ditance and DSA algorithm.
def dtw_deconv(pure, mix, gene_list_df, metric='dtw'):
    #mix[mix < 0] = 0
    gene_list_df.replace(["NaN", 'NaN', 'nan'], np.nan, inplace=True)
    num_cells = len(pure.columns)
    num_mixes = len(mix.columns)
    #mixtures_sum = [round(1 - abs(random.gauss(0, 0.045)), 2) for i in range(num_mixes)]
    mixtures_sum = [1 for i in range(num_mixes)]
    O_array = np.zeros((num_cells, num_mixes))  # The per cell sorted array - how far each mix is from the maximal mix.
    i = 0
    # Loop on all cell types.
    for cell_type in pure:
        cell_vals = []
        #gene_list_df[cell_type] = gene_list_df[cell_type].map(str.lower)
        cell_genelist = gene_list_df[cell_type].dropna().sample(frac=0.35) # round(random.gauss(0.4, 0.03), 2))  # 0.35
        #If marker list or sample list is short, don't sample.
        if (len(gene_list_df[cell_type].dropna()) < 8):  # or (len(cell_genelist) < 5)
            cell_genelist = gene_list_df[cell_type].dropna()
        #Make sure mix has all the genes in the list.
        cell_genelist = list(set(mix.index) & set(cell_genelist))
        mix_temp = mix.loc[cell_genelist]
        max_ind = mix_temp.mean().idxmax()  # Mix with maximum mean of gene expression.
        pure_temp = pure.loc[cell_genelist]
        pure_column = pure_temp[cell_type].copy()
        #if mix_temp.max().max() > 2*pure_temp.max().max():
        #    max_column = mix_temp[max_ind].copy()
            #print('.', end="")
        #else:
        max_column = pure_column
        max_column.sort_values(ascending=False, inplace=True)

        # Loop on all mixes.
        k = 0
        for mix_col in mix_temp:
            column = mix_temp[mix_col].copy()
            column = column[max_column.index]  # Sort according to the maximum column index.
            if metric == 'dtw':
                dist = dtw_metric(max_column, column)  # 0.342
            elif metric == 'avg':
                dist = np.mean(max_column - column)
            elif metric == 'abs':
                dist = np.mean(abs(max_column.values - column.values))
            elif metric == 'basic':
                dist = basic_deconv(max_column, column)
            elif metric == 'ks':
                dist, pval = stats.ks_2samp(max_column.values, column.values) # Kolmogarov Smirnoff distance.
            elif metric == 'euclid':
                dist = minkowski_distance(max_column, column, p=2)
            elif metric == 'taxi':
                dist = minkowski_distance(max_column, column, p=1)
            cell_vals.append(dist)
            k += 1
        O_array[i] = cell_vals
        i += 1

    #Scale to [0,1] in order for the sum of proportions to be meaningful. Scale along each cell.
    O_scaled = np.transpose(1 - minmax_scale(O_array, axis=1))
    solution = nnls(O_scaled, mixtures_sum)[0]
    solution_mat = np.diag(solution)
    estimate_wt = np.matmul(solution_mat, O_scaled.T)
    #estimate_wt = pd.DataFrame(data=estimate_wt, index=pure.columns).T
    #estimate_wt.index = mix.columns
    return estimate_wt


def run_dtw_deconv_ensemble(pure, mix, gene_list_df):
    num_loops = 400
    pool = mp.Pool()
    num_mixes = len(mix.columns)
    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    
    results = [pool.apply_async(dtw_deconv, args=(pure, mix, gene_list_df)) for i in range(num_loops)]
    #results = [dtw_deconv(pure, mix, gene_list_df) for i in range(num_loops)]
    for ens_i in range(num_loops):
        print('\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i].get()
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    ens_estimate_wt.index = mix.columns
    #ens_estimate_wt.to_csv('./data/results.csv')
    pool.close()
    return(ens_estimate_wt)


def calc_corr(prop, ens_estimate_wt_2):
    real_weight = pd.read_csv(f'./data/{prop}/labels.csv', index_col=0)
    both_mixes = list(set(ens_estimate_wt_2.index) & set(real_weight.index))
    ens_estimate_wt_2 = ens_estimate_wt_2.reindex(both_mixes)
    real_weight = real_weight.reindex(both_mixes)
    result = []
    for col in real_weight:
        if col in ens_estimate_wt_2.columns:
            pearson = np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]
            spearman = st.spearmanr(real_weight[col], ens_estimate_wt_2[col])
            result.append([prop, col, pearson, spearman[0], spearman[1]])
    return result


def preprocess(pure, mix):
    mix = pd.read_csv(mix, index_col=0)
    mix.index = mix.index.map(str.lower)
    mix.index = mix.index.map(str.strip)
    mix = mix.groupby(mix.index).first()
    mix.fillna(0, inplace=True)
    if mix.max().max() < 20:
        mix = 2 ** mix
    num_mixes = len(mix.columns)
    pure = pd.read_csv(pure, index_col=0)
    pure.index = pure.index.map(str.lower)
    pure.index = pure.index.map(str.strip)
    pure = pure.groupby(pure.index).first()
    pure.fillna(0, inplace=True)
    if pure.max().max() < 20:
        pure = 2 ** pure

    #Drop genes that are not shared by mix and pure.
    both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
    pure = pure.reindex(both_genes)
    mix = mix.reindex(both_genes)

    #Standardize.
    mix = (mix-mix.min())/mix.mean() #mix.apply(lambda x: (x-x.mean())/(x.std()), axis=1) #+1e-16
    pure = (pure-pure.min())/pure.mean() #pure.apply(lambda x: (x-x.mean())/(x.std()), axis=1)

    return(pure, mix)


if __name__ == '__main__':
    #Benchmark of running deconvolution of few datasets.
    result = pd.DataFrame()
    for file in ['xCell', 'ABIS', 'GSE123604', '10x', 'CIBERSORT', 'EPIC', 'TIMER', 'Abbas', 'BreastBlood', 'DeconRNASeq', 'DSA', 'RatBrain']:
        pure, mix = preprocess(f'./data/{file}/pure.csv', f'./data/{file}/mix.csv')
        gene_list_df = gene_diff(pure, mix)
        print(file)
        res = run_dtw_deconv_ensemble(pure, mix, gene_list_df)
        file_res = pd.DataFrame(calc_corr(file, res))
        result = pd.concat([result, file_res])
    result.columns = ['dataset', 'celltype', 'pearson', 'spearman', 'p']
    print(result) #.set_index(2).iloc[:,1:5])
    result.to_csv('./data/result.csv')