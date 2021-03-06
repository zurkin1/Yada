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


def run_deconv(num, mix, pure, method):
    pool = mp.Pool()
    mix = pd.read_csv(mix, index_col=0)
    mix = mix.groupby(mix.index).first()
    mix.index = mix.index.str.lower()
    num_mixes = len(mix.columns)
    print(f'Deconvolution, num_mixes: {num_mixes}')
    print(f'{method.__code__.co_name}')
    pure = pd.read_csv(pure, index_col=0)
    gene_list_df = pure

    # Gene differentiation algorithm.
    # Drop genes that are not shared by mix and pure.
    both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
    pure = pure.reindex(both_genes)
    mix_loop = mix.copy()
    mix_loop = mix_loop.reindex(both_genes)
    gene_list_df = gene_diff(pure, mix_loop)

    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    estimate_wt = np.zeros((num_cells, num_mixes))

    results = [pool.apply_async(method, args=(mix_loop, pure, gene_list_df)) for i in range(num)]
    for ens_i in range(num):
        print('\r', f"{ens_i / num * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i].get()  # estimate_wt
    ens_estimate_wt /= num
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    ens_estimate_wt.to_csv('./data/results.csv')
    pool.close()
    return(ens_estimate_wt)
