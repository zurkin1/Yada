import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import *
# import matplotlib
# import matplotlib.pyplot as plt
import scipy.stats as st
import multiprocessing as mp
import time


def calc_corr(metric, prop, platform, ens_estimate_wt_2):
    #real_weight = pd.read_csv('./data/Challenge/prop-' + prop, index_col=0).T
    real_weight = pd.read_csv(f'./data/{prop}/labels.csv', index_col=0)
    # matplotlib.style.use('ggplot')
    result = []
    for col in real_weight:
        # plt.title(f'{test}' + col)
        # plt.scatter(real_weight[col], ens_estimate_wt_2[col])
        # plt.show()
        # plt.savefig(f'c:/work/data/GEO/RNASeq/images/{test}{col}.png', bbox_inches='tight')
        # plt.close()
        if col in ens_estimate_wt_2.columns:
            pearson = np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]
            spearman = st.spearmanr(real_weight[col], ens_estimate_wt_2[col])
            result.append([metric, prop, platform, col, pearson, spearman[0], spearman[1]])
    return result


def run_deconv(mix, pure, test, metric):
    num_loops = 1
    pool = mp.Pool()
    mix = pd.read_csv(mix, index_col=0)
    mix.index = mix.index.map(str.lower)
    mix.index = mix.index.map(str.strip)
    mix = mix.groupby(mix.index).first()
    if mix.max().max() < 20:
        mix = 2 ** mix
    num_mixes = len(mix.columns)
    pure = pd.read_csv(pure, index_col=0)
    pure.index = pure.index.map(str.lower)
    pure.index = pure.index.map(str.strip)
    pure = pure.groupby(pure.index).first()
    if pure.max().max() < 20:
        pure = 2 ** pure

     # Drop genes that are not shared by mix and pure.
    both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
    pure = pure.reindex(both_genes)
    mix = mix.reindex(both_genes)

    # Gene differentiation algorithm.
    gene_list_df = gene_diff(pure, mix)

    print(f'Deconvolution, num_cells: {len(pure.columns)}, num_mixes: {num_mixes}')
    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    
    results = [pool.apply_async(dtw_deconv, args=(mix, pure, gene_list_df, metric)) for i in range(num_loops)]
    #results = [dtw_deconv(mix, pure, gene_list_df) for i in range(num_loops)]
    for ens_i in range(num_loops):
        print('\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i].get()
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    #ens_estimate_wt.to_csv('./data/results.csv')
    pool.close()
    return(ens_estimate_wt)


"""
if __name__ == '__main__':
    infiles = pd.read_csv('./data/Challenge/input.csv')
    result = pd.DataFrame([['a', 'b', 'c', 0.0, 0.0, 0.0]]*2)
    for line in infiles.iterrows():
        file = line[1]['dataset.name'] + '.csv'
        print('\n' + file)
        platform = line[1]['type']
        num_cells = str(line[1]['cells'])
        pure = './data/' + platform + '_' + num_cells + '.csv'
        #pure = './data/Challenge/pure-' + file
        res = run_dtw_deconv('./data/Challenge/mix-' + file, pure)
        res = pd.DataFrame(calc_corr(file, platform, res))
        result = pd.concat([result, res])
    result.to_csv('./data/result.csv')
"""

if __name__ == '__main__':
    result = pd.DataFrame([['a', 'b', 'c', 'd', 0.0, 0.0, 0.0]]*2)
    for metric in ['dtw', 'avg', 'abs', 'basic', 'ks', 'euclid', 'taxi']:
        for file in ['10x', 'Abbas', 'BreastBlood', 'xCell', 'CIBERSORT', 'DeconRNASeq', 'DSA', 'EPIC', 'RatBrain', 'TIMER']:
            print(metric, file, ' ', end="")
            res = run_deconv(f'./data/{file}/mix.csv', f'./data/{file}/pure.csv', file, metric)
            file_res = pd.DataFrame(calc_corr(metric, file, 'Microarray', res))
            print(file_res.describe())
            result = pd.concat([result, file_res])
    result[2:].to_csv('./data/result.csv')