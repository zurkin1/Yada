import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import *
# import matplotlib
# import matplotlib.pyplot as plt
import scipy.stats as st
import multiprocessing as mp
import time


def poincare_norm(a):
    #ret = np.arccosh(1 + 2 * (np.linalg.norm(a) ** 2) / (1 - np.linalg.norm(a) ** 2))
    ret = np.linalg.norm(a)
    return(ret)


def go_norm(expr):
    # Go function gene scaling.
    hig2vec = pd.read_csv('hig2vec.csv', index_col = 0)
    hig2vec.index = hig2vec.index.map(str.lower)
    
    #both_genes = list(set(expr.index) & set(hig2vec.index))
    #expr = expr.reindex(both_genes)
    #hig2vec = hig2vec.reindex(both_genes)
    hig2vec['norm'] = hig2vec.apply(poincare_norm, axis = 0)
    hig2vec3 = expr.join(hig2vec)['norm'].fillna(1)
    return(expr.mul(hig2vec3, axis = 0))


def calc_corr(metric, prop, platform, ens_estimate_wt_2):
    #real_weight = pd.read_csv('./data/Challenge/prop-' + prop, index_col=0).T
    real_weight = pd.read_csv(f'./data/{prop}/labels.csv', index_col=0)
    both_mixes = list(set(ens_estimate_wt_2.index) & set(real_weight.index))
    ens_estimate_wt_2 = ens_estimate_wt_2.reindex(both_mixes)
    real_weight = real_weight.reindex(both_mixes)
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


def preprocess(mix, pure, GO = False):
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

    #Drop genes that are not shared by mix and pure.
    both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
    pure = pure.reindex(both_genes)
    mix = mix.reindex(both_genes)

    #Change to UDP.
    #merged_df = pd.concat([mix, pure], axis=1)
    #df = calc_udp_multi_process(merged_df, True)
    #mix = df.iloc[:,0:num_mixes].copy()
    #pure = df.iloc[:,num_mixes:].copy()

    # Gene differentiation algorithm.
    gene_list_df = gene_diff(pure, mix)

    #Standardize.
    mix = (mix-mix.min())/mix.mean() #mix.apply(lambda x: (x-x.mean())/(x.std()), axis=1) #+1e-16
    pure = (pure-pure.min())/pure.mean() #pure.apply(lambda x: (x-x.mean())/(x.std()), axis=1)

    if GO:
        mix = go_norm(mix)
        #pure = go_norm(pure)

    return(mix, pure, gene_list_df)


def run_dtw_deconv(mix, pure):
    #metrics = ['dtw', 'avg' , 'abs', 'ks', 'euclid', 'taxi'] #'basic',
    metric = 'dtw'
    num_loops = 1
    #pool = mp.Pool()
    num_mixes = len(mix.columns)
    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    
    #results = [pool.apply_async(dtw_deconv, args=(mix, pure, gene_list_df, metric)) for i in range(num_loops)]
    results = [dtw_deconv(mix, pure, gene_list_df, metric) for i in range(num_loops)]
    for ens_i in range(num_loops):
        print('\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i] #.get()
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    ens_estimate_wt.index = mix.columns
    #ens_estimate_wt.to_csv('./data/results.csv')
    #pool.close()
    return(ens_estimate_wt)


"""
if __name__ == '__main__':
    infiles = pd.read_csv('./data/Challenge/input1.csv')
    result = pd.DataFrame([['a', 'b', 'c', 0.0, 0.0, 0.0]]*2)
    for line in infiles.iterrows():
        file = line[1]['dataset.name'] + '.csv'
        platform = line[1]['type']
        num_cells = str(line[1]['cells'])
        #pure = './data/' + platform + '_' + num_cells + '.csv'
        mix, pure, gene_list_df = preprocess(file)
        for metric in ['dtw']: #, 'avg', 'abs', 'basic', 'ks', 'euclid', 'taxi']:
            print(metric, file, ' ', end="")
            res = run_deconv(mix, pure, gene_list_df, metric)
            res = pd.DataFrame(calc_corr(metric, file, platform, res))
            result = pd.concat([result, res])
        #print(result.loc[result[0] == metric].describe())
            result.to_csv('./data/result.csv')
"""


if __name__ == '__main__':
    result = pd.DataFrame()
    for file in ['xCell', 'Mysort', 'ABIS']: #, 'GSE123604', '10x', 'CIBERSORT', 'EPIC', 'TIMER', 'Abbas', 'BreastBlood', 'DeconRNASeq', 'DSA', 'RatBrain']:
        for GO in [False, True]:
            mix, pure, gene_list_df = preprocess(f'./data/{file}/mix.csv', f'./data/{file}/pure.csv', GO)
            print(file, GO) #f'num_cells: {len(pure.columns)}, num_mixes: {len(mix.columns)}')
            for method in [cibersort, nnls_deconv_constrained]:
                res = method(mix, pure)
                file_res = pd.DataFrame(calc_corr(method.__name__, file, GO, res))
                result = pd.concat([result, file_res])

    result.to_csv('./data/result.csv')