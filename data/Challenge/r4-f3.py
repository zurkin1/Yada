import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import dtw_deconv, nnls_deconv_constrained
# import matplotlib
# import matplotlib.pyplot as plt
import scipy.stats as st
import multiprocessing as mp
import time


pd.options.display.width = 0


def calc_corr(test, ens_estimate_wt_2):
    real_weight = pd.read_csv(f'c:/input/prop-{test["dataset.name"]}.csv', index_col=0).T
    # matplotlib.style.use('ggplot')
    for col in real_weight:
        # plt.title(f'{test}' + col)
        # plt.scatter(real_weight[col], ens_estimate_wt_2[col])
        # plt.show()
        # plt.savefig(f'c:/work/data/GEO/RNASeq/images/{test}{col}.png', bbox_inches='tight')
        # plt.close()
        if col in ens_estimate_wt_2.columns:
            print(f'{test["dataset.name"]}, {col}, {np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]}, {st.spearmanr(real_weight[col], ens_estimate_wt_2[col])}')


if __name__ == '__main__':
    desc_file = pd.read_csv('/input/input.csv')
    output_df = pd.DataFrame(columns=['dataset.name', 'sample.id', 'cell.type', 'prediction'])
    cells_lm14 = ['memory.B.cells', 'naive.B.cells', 'memory.CD4.T.cells', 'naive.CD4.T.cells', 'regulatory.T.cells', 'memory.CD8.T.cells', 'naive.CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytes', 'myeloid.dendritic.cells', 'macrophages', 'fibroblasts', 'endothelial.cells']
    cells_lm8 = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytic.lineage', 'fibroblasts', 'endothelial.cells']
    # mix = pd.read_csv('c:/data/GEO/Microarray/mix-72642-110085/mix.csv', index_col=0)

    pool = mp.Pool()
    k = 0
    for test_index, test in desc_file.iterrows():
        mix = pd.read_csv('/input/' + test['hugo.expr.file'], index_col=0)
        num_mixes = len(mix.columns)

        methods = [
            [10, dtw_deconv, ['memory.CD4.T.cells', 'fibroblasts', 'regulatory.T.cells'], 'pure_orig8v2_14', 10],  # 400
            [10, dtw_deconv, ['naive.CD4.T.cells', 'naive.CD8.T.cells'], 'pure_quanTISeq', 10],
            [10, dtw_deconv, ['NK.cells', 'macrophages', 'myeloid.dendritic.cells', 'neutrophils'], 'LM14', 10],  # 400
            [10, dtw_deconv, ['naive.B.cells'], 'ABIS', 10],
            [10, nnls_deconv_constrained, ['endothelial.cells', 'monocytes', 'memory.B.cells'], 'pure_orig8v2_14', 10],
            [10, dtw_deconv, ['memory.CD8.T.cells'], 'pure_orig8v3_14', 10]
        ]

        print(f'{test["dataset.name"]}, {test_index}, num_mixes: {num_mixes}, {test["cancer.type"]}, {test["platform"]}, {test["scale"]}, {test["normalization"]}, {test["symbol.compression.function"]}')
        ens_estimate_wt_2 = pd.DataFrame(data=np.zeros((num_mixes, len(cells_lm14))), columns=cells_lm14)
        for method in methods:
            print(f'{method[1].__code__.co_name}, {method[3]}')
            pure = pd.read_csv('./data/Challenge/' + method[3] + '.csv', index_col=0)
            # Drop genes that are not shared by mix and pure.
            both_genes = list(set(mix.index) & set(pure.index))  # - set(BRCA)
            pure = pure.reindex(both_genes)
            mix_loop = mix.copy()
            mix_loop = mix_loop.reindex(both_genes)
            gene_list_df = gene_diff(pure, mix_loop)
            # gene_list_df.to_csv('gene_list_df.csv')
            #gene_list_df = pd.read_csv('gene_list_df.csv', index_col=0)

            num_cells = len(pure.columns)
            ens_estimate_wt = np.zeros((num_cells, num_mixes))
            estimate_wt = np.zeros((num_cells, num_mixes))
            # print(time.ctime())
            results = [pool.apply_async(method[1], args=(mix_loop, pure, gene_list_df)) for i in range(method[0])]
            for ens_i in range(method[0]):
                print('\r', f"{ens_i / method[0] * 100:.0f}%", end='')
                #estimate_wt = method[1](mix_loop, pure, gene_list_df)
                ens_estimate_wt += results[ens_i].get()  # estimate_wt
            ens_estimate_wt /= method[0]
            ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
            #ens_estimate_wt_2 = pd.concat([ens_estimate_wt_2, ens_estimate_wt[method[2]]], axis=1)
            ens_estimate_wt_2[method[2]] += ens_estimate_wt[method[2]]
            # print(time.ctime())

        #ens_estimate_wt_2 = ens_estimate_wt_2.groupby(ens_estimate_wt_2.index).mean()
        for col in ens_estimate_wt_2:  # Loop on cell-types.
            for j in range(num_mixes):  # Loop on samples.
                output_df.loc[k] = [test['dataset.name'], mix.columns[j], col, ens_estimate_wt_2[col][j]]
                k += 1

        calc_corr(test, ens_estimate_wt_2)
    output_df.to_csv('/input/predictions.csv', index=None)
    pool.close()
