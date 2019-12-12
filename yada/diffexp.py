import numpy as np
import pandas as pd

#This function does a gene differentiation analysis. It finds for each cell types its majorodizing genes.
#We also get a factor of how much (strongly majordizing) each gene is.
def gene_diff(pure, platform, normalization):
    gene_list_per_cell = {}
    # Loop on all cells
    for cell_type in pure:
        gene_list = []
        data = pure.copy()
        for factor in range(100, 1, -1):
            data2 = data.copy()
            data2['max_val'] = data.drop(cell_type, axis=1).max(axis=1)
            data2['min_val'] = data.drop(cell_type, axis=1).min(axis=1)
            data2 = data2.loc[(data2[cell_type] > 0.5 * factor * data2.max_val) & (data2.max_val > 0)] #Ignore cases where all other mixes are zeros.
            if (len(data2) == 0):
                continue
            data2.sort_values(cell_type, inplace=True, ascending=False)
            gene_list += list(zip(data2.index, [factor for i in range(len(data2.index))]))
            data.drop(data2.index, inplace=True)
            if (len(gene_list) > 80):
                gene_list = gene_list[0:80]
                break
        #print(f'{cell_type} : {gene_list}')
        gene_list_per_cell[cell_type] = gene_list

    gene_list_df = pd.DataFrame(dict([(k, pd.Series([i for (i, j) in v])) for k, v in gene_list_per_cell.items()]), columns=pure.columns)
    factor_list = pd.DataFrame(dict([(k, pd.Series([j * 0.5 for (i, j) in v])) for k, v in gene_list_per_cell.items()]), columns=pure.columns)

    if 'Affymetrix Human Gene 1.1 ST' in platform:
        gene_list_df['neutrophils'][0:80] = ['TMEM55A', 'MEGF9', 'WDFY3', 'BEST1', 'STX3', 'GCA', 'TLR8', 'UBXN2B', 'MXD1', 'GLT1D1', 'TMEM88', 'ITPRIP', 'ROPN1L', 'RNF149', 'DCP2',
                                             'TLR6', 'NFIL3', 'CEACAM4', 'HAL', 'MNDA', 'PFKFB4', 'RP2', 'TET2', 'RNF24', 'LRRC4', 'CMTM2', 'REPS2', 'KIAA0232', 'LIN7A', 'AQP9',
                                             'LRRK2', 'CMTM6', 'NRBF2', 'LILRB3', 'NFAM1', 'TLR2', 'ALPK1', 'MOSC1', 'DOCK5', 'CLEC4E', 'SLC22A4', 'ACSL1', 'KCNE3', 'QPCT', 'SULT1B1',
                                             'PLXNC1', 'DENND3', 'FRAT2', 'LRP10', 'TMEM49', 'NAMPT', 'RASGRP4', 'TNFRSF10C', 'TRIM25', 'USP32', 'BST1', 'DKFZp761E198', 'CYB5R4',
                                             'MTMR3', 'CPD', 'RRAGD', 'BCL6', 'ZBTB34', 'CPEB4', 'FCGR2A', 'EGLN1', 'IL1R2', 'RBM47', 'CSF3R', 'NPL', 'PPP1R3B', 'FPR1', 'TLR4', 'KCNJ15',
                                             'PPP4R1', 'CEACAM3', 'SPOPL', 'RTN3', 'LAMP2', 'HSPA6']
        gene_list_df['monocytes'][0:80] = ['NAGA', 'CD68', 'ADAP2', 'ANXA2', 'SLC43A3', 'ZNF385A', 'RASSF4', 'RIN2', 'CPVL', 'CD33', 'LOC284837', 'PLXNB2', 'CRTAP', 'CECR1',
                                                   'KLF4', 'ANXA2P2', 'NPC2', 'CD86', 'TTYH2', 'CD300C', 'KCTD12', 'DPYSL2', 'PEA15', 'CYBB', 'CST3', 'CTNND1', 'PLD3', 'C16orf70', 'GPBAR1',
                                                   'MS4A7', 'FLVCR2', 'MYCL1', 'ADAM15', 'CTTNBP2NL', 'KCNMB1', 'CCR2', 'SLC46A2', 'MAN2B1', 'RPS6KA4', 'SLC37A2',
                                                   'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan',
                                                   'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'
                                                   ]
    if 'FARMS' in normalization:
        gene_list_df['neutrophils'][0:80] = ['GAB2', 'BEST1', 'ADAMTSL4-AS1', 'LOC729603', 'RBM47', 'ZDHHC18', 'DOCK5', 'RASGRP4', 'PFKFB4', 'CHST15', 'STX3', 'PELI2', 'TNFRSF1A',
                                           'PHC2', 'ABTB1', 'PREX1', 'NDST1', 'RN7SL600P', 'NFAM1', 'CEACAM3', 'NOTCH1', 'EPHB1', 'DYSF', 'NINJ1', 'C5AR1', 'KDM6B', 'C19orf35',
                                           'TMEM127', 'LILRB3', 'MYO1F', 'CXCR2', 'ARAP1', 'MBOAT7', 'RN7SL473P', 'RNF24', 'SEC14L1', 'SIGLEC9', 'REPS2', 'SIRPD', 'ITPRIP', 'IL6R',
                                           'IMPDH1', 'CEACAM4', 'CSF3R', 'NCF4', 'PADI4', 'SLED1', 'LRP10', 'HAL', 'NUMB', 'DENND3', 'ABHD5', 'DEF8', 'ATG2A', 'ITPK1', 'STK40', 'MXD1',
                                           'NHSL2', 'KLHL21', 'MEGF9', 'IGF2R', 'IL1R1', 'LAT2', 'PLAUR', 'ERV3-1', 'NADK', 'WAS', 'OR52K2', 'SLC45A4', 'ADCY4', 'PRRG4', 'SLC6A6',
                                           'LINC00999', 'PPP1R3B', 'MGAM', 'FCGRT', 'CXCR1', 'GLT1D1', 'DHX34', 'MAST3']
        gene_list_df['monocytes'][0:80] = ['LGALS1', 'CD300C', 'TTYH2', 'GSTP1', 'PHPT1', 'CLCN5', 'CECR1', 'CD68', 'RHOC', 'CST3', 'NLN', 'PLXNB2', 'KIAA0930', 'CD300E', 'CNP', 'TMEM205',
                                                   'CALHM2', 'GPBAR1', 'DPYSL2', 'CD1D', 'GLB1L', 'SH3TC1', 'UBXN11', 'AIMP2', 'ACP2', 'PLXND1', 'TBC1D8', 'CPVL', 'ZNF385A', 'ABI3', 'COMT',
                                                   'TSPAN4', 'SLC43A3', 'FGD2', 'RNH1', 'TMEM150B', 'MS4A14', 'CSF1R', 'ANAPC2', 'ATP6V1F',
                                                   'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan',
                                                   'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan'
                                                   ]
    return gene_list_df, factor_list