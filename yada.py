import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import *
# import matplotlib
# import matplotlib.pyplot as plt
import scipy.stats as st
import multiprocessing as mp
import time


def calc_corr(test, ens_estimate_wt_2):
    real_weight = pd.read_csv(test, index_col=0)
    # matplotlib.style.use('ggplot')
    for col in real_weight:
        # plt.title(f'{test}' + col)
        # plt.scatter(real_weight[col], ens_estimate_wt_2[col])
        # plt.show()
        # plt.savefig(f'c:/work/data/GEO/RNASeq/images/{test}{col}.png', bbox_inches='tight')
        # plt.close()
        if col in ens_estimate_wt_2.columns:
            print(f'{col}, {np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]}, {st.spearmanr(real_weight[col], ens_estimate_wt_2[col])}')


def run_dtw_deconv(mix, pure):
    num_loops = 100
    pool = mp.Pool()
    mix = pd.read_csv(mix, index_col=0)
    mix.index = mix.index.map(str.lower)
    mix.index = mix.index.map(str.strip)
    mix = mix.groupby(mix.index).first()
    num_mixes = len(mix.columns)
    pure = pd.read_csv(pure, index_col=0)
    pure.index = pure.index.map(str.lower)
    pure.index = pure.index.map(str.strip)
    pure = pure.groupby(pure.index).first()
    print(f'Deconvolution, num_cells: {len(pure.columns)}, num_mixes: {num_mixes}')

     # Drop genes that are not shared by mix and pure.
    both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
    pure = pure.reindex(both_genes)
    mix = mix.reindex(both_genes)

    # Gene differentiation algorithm.
    gene_list_df = gene_diff(pure, mix)

    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    estimate_wt = np.zeros((num_cells, num_mixes))

    #results = [pool.apply_async(dtw_deconv, args=(mix, pure, gene_list_df)) for i in range(num_loops)]
    results = [dtw_deconv(mix, pure, gene_list_df) for i in range(num_loops)]
    for ens_i in range(num_loops):
        print('\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i] #.get()  # estimate_wt
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    ens_estimate_wt.to_csv('./data/results.csv')
    pool.close()
    return(ens_estimate_wt)