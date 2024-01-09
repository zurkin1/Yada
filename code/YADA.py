import numpy as np
import pandas as pd
from diffexp import gene_diff
import scipy.stats as st
import time
from sklearn.preprocessing import minmax_scale
from scipy.optimize import nnls, minimize
from scipy import stats
from scipy.spatial import minkowski_distance
#from tslearn import metrics
import logging
import warnings
import random
from random import choice
#mport similaritymeasures
#from similaritymeasures import pcm
from tqdm import tqdm
import multiprocessing as mp
import sys
import cProfile


warnings.filterwarnings("ignore")
logger = logging.getLogger("pymc3")
logger.propagate = False
logger.setLevel(logging.ERROR)

FULL_PATH = 'c:\\Users\\danili\\Documents\\GitHub\\Yada\\data\\'

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


# This function calculates a deconvolution algorithm using DTW distance and YADA algorithm.
def dtw_deconv(pure, mix, gene_list_df, metric='avg'):
    #mix[mix < 0] = 0
    gene_list_df.replace(["NaN", 'NaN', 'nan', ''], np.nan, inplace=True)
    num_cells = len(pure.columns)
    num_mixes = len(mix.columns)
    #mixtures_sum = [round(1 - abs(random.gauss(0, 0.045)), 2) for i in range(num_mixes)]
    mixtures_sum = [1 for i in range(num_mixes)]
    O_array = np.zeros((num_cells, num_mixes))  #The proportion matrix approximation.
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
        max_column = mix_temp[max_ind].copy()
        #pure_temp = pure.loc[cell_genelist]
        #pure_column = pure_temp[cell_type].copy()
        #max_column = pure_column
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
        logger.log(logging.DEBUG, '\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i].get()
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, columns=mix.columns, index=pure.columns).T
    ens_estimate_wt.index = mix.columns
    #ens_estimate_wt.to_csv('../data/results.csv')
    pool.close()
    return(ens_estimate_wt)


def calc_corr(prop, ens_estimate_wt_2):
    real_weight = pd.read_csv(FULL_PATH + f'{prop}/labels.csv', index_col=0)
    results = pd.read_csv(FULL_PATH + f'{prop}/results.csv')
    both_mixes = list(set(ens_estimate_wt_2.index) & set(real_weight.index))
    ens_estimate_wt_2 = ens_estimate_wt_2.reindex(both_mixes)
    real_weight = real_weight.reindex(both_mixes)
    for col in real_weight:
        if col in ens_estimate_wt_2.columns:
            pearson = np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]
            spearman = st.spearmanr(real_weight[col], ens_estimate_wt_2[col])
            results.loc[len(results)] = ([prop, col, pearson, spearman[0], spearman[1]])
    results.to_csv(FULL_PATH + f'{prop}/results.csv', index=None)

def preprocess_only_marker(pure, mix):
    mix = pd.read_csv(mix, index_col=0)
    mix.index = mix.index.map(str.lower)
    mix.index = mix.index.map(str.strip)
    mix = mix.groupby(mix.index).first()
    mix.fillna(0, inplace=True)
    if mix.max().max() < 20:
        mix = 2 ** mix
    num_mixes = len(mix.columns)
	
    pure = pd.read_csv(pure, index_col=0)
    for col in pure.columns:
       pure[col] = [str.lower(word) if word is not np.nan else '' for word in pure[col]]

    #Standardize.
    mix = (mix-mix.min())/mix.mean() #mix.apply(lambda x: (x-x.mean())/(x.std()), axis=1) #+1e-16
	
    pure_genes = set()
    for i in range(pure.values.shape[0]):
      pure_genes = pure_genes.union(set(pure.values[i]))
    both_genes = list(set(mix.index) & pure_genes)
    mix = mix.reindex(both_genes)
    for col in pure.columns:
       pure[col] = [word if word in both_genes else '' for word in pure[col]]
    return(pure, mix)


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


def run_benchmark():
    #Benchmark of running deconvolution of few datasets.
    result = pd.DataFrame()
    for file in ['ABIS', 'GSE123604', '10x', 'CIBERSORT', 'EPIC', 'TIMER', 'Abbas', 'BreastBlood', 'DeconRNASeq', 'DSA', 'RatBrain']:
        pure, mix = preprocess(f'../data/{file}/pure.csv', f'../data/{file}/mix.csv')
        gene_list_df = gene_diff(pure, mix)
        logger.log(logging.DEBUG, file)
        res = run_dtw_deconv_ensemble(pure, mix, gene_list_df)
        file_res = pd.DataFrame(calc_corr(file, res))
        result = pd.concat([result, file_res])
    result.columns = ['dataset', 'celltype', 'pearson', 'spearman', 'p']
    logger.log(logging.DEBUG, result)
    result.to_csv('../data/result.csv')


#Randomly select a cell from each group to create the pure matrix.
def get_random_indices(df):
    groups = ['Group1', 'Group2', 'Group3']  #Ensure the order of groups.
    indices = []

    for group in sorted(groups):  #Sort groups in ascending order.
        group_indices = df[df['Group'] == group].index.tolist() #Assuming 'df' is pandas DataFrame with a 'Group' column.
        indices.append(random.choice(group_indices))  #Add a random index for each group.

    return indices

#@profile kernprof -l -v yada.py
def create_simulator():
    #Pure
    logger.log(logging.DEBUG, 'Creating pure matrix.')
    tpm = pd.read_csv(FULL_PATH + 'sim/splater.csv', index_col=0)
    cells = pd.read_csv(FULL_PATH + 'sim/celldata.csv', index_col=0)
    random_indices = get_random_indices(cells)
    pure = tpm.T.loc[random_indices].T.copy()
    
    #Mixture
    logger.log(logging.DEBUG, 'Creating mixture.')
    num_mixes = 50
    mix = pure.copy()
    prop_data = np.empty((0,len(pure.columns)))
    # Loop on mixes.
    for i in range(0,num_mixes):
        prop = np.random.dirichlet(np.ones(len(pure.columns)), size=1)
        j = 0
        #Create a temporary DataFrame to hold intermediate results.
        temp_mix = pd.DataFrame()
        #Loop on cell types.
        for col in pure.columns:
            noise = np.random.normal(1, 0.4, len(mix[col]))
            temp_mix[col] = mix[col] * prop[0][j] * noise
            j+=1
        # Sum along the columns to get the updated mix[f'mix{i}']
        mix[f'mix{i}'] = temp_mix.sum(axis=1)
        prop_data = np.append(prop_data, prop, axis=0)

    columns_to_drop = list(pure.columns)
    mix.drop(columns_to_drop, axis=1, inplace=True)
    labels = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,num_mixes)], index=pure.columns).T
    labels.to_csv(FULL_PATH + 'sim/labels.csv')
    mix = mix[[f'mix{i}' for i in range(0,num_mixes)]]
    
    #Markers
    logger.log('Creating markers.')
    de = pd.read_csv(FULL_PATH + 'sim/dedata.csv', index_col=0)
    de2 = de.loc[(de.DEFacGroup1 > 1.8) | (de.DEFacGroup2 > 1.8) | (de.DEFacGroup3 > 1.8)].copy()
    de2.drop(['Gene','BaseGeneMean','OutlierFactor','GeneMean', 'Length'], axis=1, inplace=True)
    #Convert columns to numeric values.
    df = de2
    df = df.apply(pd.to_numeric, errors='coerce')

    #Initialize empty lists for each group.
    group1, group2, group3 = [], [], []

    #Iterate through rows to assign genes to respective groups.
    for gene, row in df.iterrows():
        max_val = row.max()
        group = row[row == max_val].index[0]
        if group == 'DEFacGroup1':
            group1.append(gene)
        elif group == 'DEFacGroup2':
            group2.append(gene)
        elif group == 'DEFacGroup3':
            group3.append(gene)

    #Fill empty lists with None to equalize lengths.
    max_length = max(len(group1), len(group2), len(group3))
    group1.extend([None] * (max_length - len(group1)))
    group2.extend([None] * (max_length - len(group2)))
    group3.extend([None] * (max_length - len(group3)))

    #Create a new DataFrame with genes split into respective groups.
    markers = pd.DataFrame({
        'DEFacGroup1': group1,
        'DEFacGroup2': group2,
        'DEFacGroup3': group3
    })
    markers.columns = pure.columns
    return pure, mix, labels, markers

def run_simulator():
    pure, mix, labels, markers = create_simulator()
    #pure = pd.read_csv(FULL_PATH + 'sim/pure.csv', index_col=0)
    #mix = pd.read_csv(FULL_PATH + 'sim/mix.csv', index_col=0)
    #markers = pd.read_csv(FULL_PATH + 'sim/markers.csv', index_col=0)
    result = run_dtw_deconv_ensemble(pure, mix, markers)
    calc_corr('sim', result)
    

if __name__ == '__main__':
    #cProfile.run('create_simulator()', sort='cumulative')
    #run_benchmark()
    run_simulator()