#Used in round 3 submit 4.
import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import cibersort, dtw_deconv, nnls_deconv_constrained
#import matplotlib
#import matplotlib.pyplot as plt
#import scipy.stats as st


def calc_corr(ens_estimate_wt_2, test):
    real_weight = pd.read_csv(f'c:/input/prop-{test["dataset.name"]}.csv', index_col=0).T #-Affymetrix HG-U133 Plus 2.0, MAS5
    ens_estimate_wt_2 = ens_estimate_wt_2[real_weight.columns]
    #matplotlib.style.use('ggplot')
    for col in real_weight:
        #plt.title(f'{test}' + col)
        #plt.scatter(real_weight[col], ens_estimate_wt_2[col])
        # plt.show()
        # plt.savefig(f'c:/work/data/GEO/RNASeq/images/{test}{col}.png', bbox_inches='tight')
        #plt.close()
        print(f'Pearson, Spearman correlation of {col}: {np.corrcoef(real_weight[col], ens_estimate_wt_2[col])[0][1]}, {st.spearmanr(real_weight[col], ens_estimate_wt_2[col])}')


if __name__ == '__main__':
    desc_file = pd.read_csv('/input/input.csv')
    output_df = pd.DataFrame(columns=['dataset.name', 'sample.id', 'cell.type', 'prediction'])
    cells_lm14 = ['memory.B.cells','naive.B.cells','memory.CD4.T.cells','naive.CD4.T.cells','regulatory.T.cells','memory.CD8.T.cells','naive.CD8.T.cells','NK.cells','neutrophils','monocytes','myeloid.dendritic.cells','macrophages','fibroblasts','endothelial.cells']
    cells_lm8 = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytic.lineage', 'fibroblasts', 'endothelial.cells']
    #mix = pd.read_csv('c:/data/GEO/Microarray/mix-72642-110085/mix.csv', index_col=0)

    k = 0
    for test_index, test in desc_file.iterrows():
        mix = pd.read_csv('/input/' + test['hugo.expr.file'], index_col=0)
        # Anti-log if max < 50 in mixture file
        if 'Log2' in test["scale"]:
            mix = 2 ** mix
        # If Microarray, quantile normalize data (only in some cases like MAS5, RMA already does so).
        if 'MAS5' in test['platform']:
            rank_mean = mix.stack().groupby(mix.rank(method='first').stack().astype(int)).mean()
            mix = mix.rank(method='min').stack().astype(int).map(rank_mean).unstack()

        num_mixes = len(mix.columns)
        num_cells = 14
        methods = [
                   [400, dtw_deconv, ['naive.B.cells', 'regulatory.T.cells'], 'pure_orig8v2_14'],
                   [400, dtw_deconv, ['NK.cells'], 'ImmunoState'],
                   [400, dtw_deconv, ['fibroblasts', 'endothelial.cells'], 'LM8'],  # 400
                   [1, cibersort, ['naive.CD4.T.cells'], 'LM14'],
                   [1, cibersort, ['naive.CD8.T.cells'], 'ImmunoState'],
                   [400, dtw_deconv, ['myeloid.dendritic.cells', 'monocytes', 'memory.CD8.T.cells'], 'pure_orig8v2_14'],  # 400
                   [400, dtw_deconv, ['memory.CD4.T.cells', 'macrophages'], 'LM14'],  # 40
                   [1, cibersort, ['memory.B.cells'], 'pure_orig8v2_14']
                  ]
        if 'HG-U133 Plus 2.0' in test["platform"] and 'MAS5' in test["normalization"]:
            methods += [[1, cibersort, ['neutrophils'], 'pure_orig8v2_14']]
        else:
            methods += [[400, dtw_deconv, ['neutrophils'], 'pure_orig8v2_8']]
        if 'CPM' in test["normalization"] or 'TPM' in test["normalization"]:
            methods = [[50, nnls_deconv_constrained, cells_lm14, 'pure_orig8v2_14']] #50
        if 'Illumina HiSeq 4000' in test["platform"]:
            methods = [[50, nnls_deconv_constrained, ['naive.B.cells', 'memory.CD4.T.cells', 'naive.CD4.T.cells', 'NK.cells', 'neutrophils', 'monocytes', 'myeloid.dendritic.cells', 'macrophages', 'fibroblasts', 'endothelial.cells'], 'pure_orig8v2_14'],
                       [1, cibersort, ['memory.B.cells'], 'LM14'],
                       [400, dtw_deconv, ['naive.CD8.T.cells', 'memory.CD8.T.cells', 'regulatory.T.cells'], 'pure_orig8v2_14']
                      ]

        ens_estimate_wt_2 = pd.DataFrame()
        for method in methods:
            print(f'{method[1].__code__.co_name} {test["dataset.name"]}, {test_index}, num_mixes: {num_mixes}, {test["cancer.type"]}, {test["platform"]}, {test["scale"]}, {test["normalization"]}, {test["symbol.compression.function"]}, {method[3]}')
            pure = pd.read_csv('/' + method[3] + '.csv', index_col=0)

            # Drop genes that are not shared by mix and pure.
            both_genes = list(set(mix.index) & set(pure.index))
            pure = pure.reindex(both_genes)
            mix_loop = mix.copy()
            mix_loop = mix_loop.reindex(both_genes)
            gene_list_df, factor_list = gene_diff(pure, test["platform"], test["normalization"])

            num_cells = len(pure.columns)
            ens_estimate_wt = np.zeros((num_cells, num_mixes))
            estimate_wt = np.zeros((num_cells, num_mixes))
            for ens_i in range(method[0]):
                print('\r', f"{ens_i / method[0] * 100:.0f}%", end='')
                estimate_wt = method[1](mix_loop, pure, gene_list_df)
                ens_estimate_wt += estimate_wt
            ens_estimate_wt /= method[0]
            ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
            ens_estimate_wt_2 = pd.concat([ens_estimate_wt_2, ens_estimate_wt[method[2]]], axis=1)

        for col in ens_estimate_wt_2: #Loop on cell-types.
          for j in range(num_mixes): #Loop on samples.
              output_df.loc[k] = [test['dataset.name'], mix.columns[j], col, ens_estimate_wt_2[col][j]]
              k += 1

    output_df.to_csv('/output/predictions.csv', index=None)
    #calc_corr(ens_estimate_wt_2, test)