import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import *
import scipy.stats as st
import time


def calc_corr(metric, prop, ens_estimate_wt_2):
    real_weight = pd.read_csv(f'./data/{prop}/labels.csv', index_col=0)
    both_mixes = list(set(ens_estimate_wt_2.index) & set(real_weight.index))
    ens_estimate_wt_2 = ens_estimate_wt_2.reindex(both_mixes)
    real_weight = real_weight.reindex(both_mixes)
    result = []
    for col in real_weight:
        if col in ens_estimate_wt_2.columns:
            pearson = np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]
            spearman = st.spearmanr(real_weight[col], ens_estimate_wt_2[col])
            result.append([metric, prop, col, pearson, spearman[0], spearman[1]])
    return result


def preprocess(mix, pure, GO = False):
    mix = pd.read_csv(mix, index_col=0)
    mix.index = mix.index.map(str.lower)
    mix.index = mix.index.map(str.strip)
    mix = mix.groupby(mix.index).first()
    mix.fillna(0, inplace=True)
    if mix.max().max() < 20:
        mix = 2 ** mix
    num_mixes = len(mix.columns)
    pure = pd.read_csv(pure, index_col=0).iloc[:,:6]
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

    # Gene differentiation algorithm.
    gene_list_df = gene_diff(pure, mix)

    #Standardize.
    mix = (mix-mix.min())/mix.mean() #mix.apply(lambda x: (x-x.mean())/(x.std()), axis=1) #+1e-16
    pure = (pure-pure.min())/pure.mean() #pure.apply(lambda x: (x-x.mean())/(x.std()), axis=1)

    return(mix, pure, gene_list_df)


if __name__ == '__main__':
    result = pd.DataFrame()
    for file in ['xCell', 'Mysort', 'ABIS', 'GSE123604', '10x', 'CIBERSORT', 'EPIC', 'TIMER', 'Abbas', 'BreastBlood', 'DeconRNASeq', 'DSA', 'RatBrain']:
        mix, pure, gene_list_df = preprocess(f'./data/{file}/mix.csv', f'./data/{file}/pure.csv')
        for method in [run_dtw_deconv_ensemble, cibersort]:
            print(file)
            res = method(mix, pure, gene_list_df)
            file_res = pd.DataFrame(calc_corr(method.__name__, file, res))
            result = pd.concat([result, file_res])
    #print(result.set_index(2).iloc[:,1:5])
    result.to_csv('./data/result.csv')