import numpy as np
import pandas as pd
import GEOparse
#from diffexp import gene_diff
from random import randrange


# cells_lm8 = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytic.lineage', 'fibroblasts', 'endothelial.cells']
# cells_lm14 = ['memory.B.cells','naive.B.cells','memory.CD4.T.cells','naive.CD4.T.cells','naive.CD8.T.cells','NK.cells','neutrophils','monocytes','macrophages','fibroblasts','endothelial.cells','regulatory.T.cells','memory.CD8.T.cells','myeloid.dendritic.cells']
# quantiseq_cells = ['B.cells', 'Macrophages M1', 'Macrophages M2', 'monocytes', 'neutrophils', 'NK.cells', 'CD4.T.cells', 'CD8.T.cells', 'regulatory.T.cells', 'dendritic.cells']
def create_microarray_lm():
    studies = pd.read_csv(f'c:/data/Geo/Microarray/samples-cg.csv')
    global_table = pd.DataFrame()
    gses = {}

    # Loop on all cell types
    for index, row in studies.iterrows():
        if row['Data set'] not in gses:
            gses[row['Data set']] = GEOparse.get_GEO(geo=row['Data set'], destdir="c:/data/Geo/Microarray/Gse")
            gse = gses[row['Data set']]
            # Parse platform data.
            for gpl_name, gpl in gse.gpls.items():
                # print("Metadata:",)
                # for key, value in gpl.metadata.items():
                #    print(" - %s : %s" % (key, ", ".join(value)))
                gpl_table = gpl.table[['ID', 'Gene Symbol']]
                # Split rows that have probes mapped to multiple genes: probex,geneA///geneB///geneC
                gpl_table = gpl_table.set_index(['ID'])['Gene Symbol'].str.split('\/\/\/', expand=True).stack().reset_index()
                gpl_table.drop(['level_1'], inplace=True, axis=1)
                gpl_table.columns = ['ID_REF', 'Symbol']
                gpl_table['Symbol'] = gpl_table.Symbol.str.strip()
                break
        else:
            gse = gses[row['Data set']]

        # Count how many samples from a given cell type to be able to calculate the mean.
        count = len(studies.loc[studies.celltype == row['celltype']])
        # Parse GSM examples.
        for gsm_name, gsm in gse.gsms.items():
            if gsm_name == row['Sample ID']:
                print("Name: ", gsm_name)
                # print("Metadata:",)
                # for key, value in gsm.metadata.items():
                #    print(" - %s : %s" % (key, ", ".join(value)))
                table = gsm.table
                table = pd.merge(table, gpl_table, on='ID_REF', how='left')
                table.drop('ID_REF', axis=1, inplace=True)
                table.dropna(inplace=True)
                table = table.groupby('Symbol').VALUE.min() / count
                table.rename(row['celltype'], inplace=True)

        print(f'{row[0]}.')
        global_table = pd.concat([global_table, table], axis=1)

    # Combine all cell samples from one group of cell type to one column.
    temp = global_table.T
    global_table = temp.groupby(temp.index).sum().T

    if('NK cells resting' in list(studies.celltype.unique())):
        global_table['NK.cells'] = global_table[['NK cells resting', 'NK cells activated']].max(axis=1)
        global_table.drop(['NK cells resting', 'NK cells activated'], inplace=True, axis=1)

    if('T cells CD4 memory resting' in list(studies.celltype.unique())):
        global_table['memory.CD4.T.cells'] = global_table[['T cells CD4 memory resting', 'T cells CD4 memory activated']].max(axis=1)
        global_table.drop(['T cells CD4 memory resting', 'T cells CD4 memory activated'], inplace=True, axis=1)

    if('Dendritic cells resting' in list(studies.celltype.unique())):
        global_table['myeloid.dendritic.cells'] = global_table[['Dendritic cells resting', 'Dendritic cells activated']].max(axis=1)
        global_table.drop(['Dendritic cells resting', 'Dendritic cells activated'], inplace=True, axis=1)

    if('Macrophages M0' in list(studies.celltype.unique())):
        global_table['macrophages'] = global_table[['Macrophages M0', 'Macrophages M1', 'Macrophages M2']].max(axis=1)
        global_table.drop(['Macrophages M0', 'Macrophages M1', 'Macrophages M2'], inplace=True, axis=1)

    # Anti-log if needed.
    for col in global_table:
        if global_table[col].max() < 20:
            global_table[col] = 2 ** global_table[col]

    # Selecting genes for cell types that are not in LM22 (Endothelial, Fibroblasts and CD8) is better done before quantile normalization.
    gene_list_df, _ = gene_diff(global_table)

    # Quantile normalize data.
    rank_mean = global_table.stack().groupby(global_table.rank(method='first').stack().astype(int)).mean()
    global_table = global_table.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    # Gene selection for cell types that don't show in LM22.
    lm22 = pd.read_csv('c:/work/src/Yada/yada/data/CIBERSORT/pure.csv', index_col=0)
    global_table = global_table.reindex(set(lm22.index) |
                                        set(gene_list_df.fibroblasts))  # |
                                        # set(gene_list_df['memory.CD8.T.cells']) |
                                        # set(gene_list_df['Macrophages M1']) |
                                        # set(gene_list_df['Macrophages M2']))
    global_table.dropna(inplace=True)
    global_table.to_csv(f'c:/data/Geo/Microarray/LM8.csv')


# Create RNASeq matrix using quanTIseq data. We were not able to reproduce the matrix in the paper.
def create_lm_rnaseq():
    gse39652 = pd.read_csv('c:/data/Geo/RNASeq/GSE36952.csv', index_col=0)
    gse39652 = gse39652 * 1000000 / gse39652.sum()
    gse39652['Macrophages M1'] = gse39652[['M1_1', 'M1_2', 'M1_3']].min(axis=1)
    gse39652['Macrophages M2'] = gse39652[['M2_1', 'M2_2', 'M2_3']].min(axis=1)
    gse39652.drop(['M1_1', 'M1_2', 'M1_3', 'M2_1', 'M2_2', 'M2_3'], axis=1, inplace=True)
    gse39652 = gse39652.reset_index()
    gse39652.columns = ['genes', 'Macrophages M1', 'Macrophages M2']

    gse64655 = pd.read_csv('c:/data/Geo/RNASeq/GSE64655.csv')
    dict = gse64655[['Gene ID', 'genes']].copy()
    # dict.set_index(dict['Gene ID'], inplace=True)
    # dict.drop_duplicates(inplace=True)
    # dict = dict.reset_index()

    gse64655.drop('Gene ID', inplace=True, axis=1)
    gse64655.set_index('genes', inplace=True)
    gse64655 = gse64655 * 1000000 / gse64655.sum()
    gse64655['dendritic.cells'] = gse64655[['dendritic.cells1', 'dendritic.cells2', 'dendritic.cells3', 'dendritic.cells4',
                                            'dendritic.cells5', 'dendritic.cells6', 'dendritic.cells7', 'dendritic.cells8']].min(axis=1)
    gse64655.drop(['dendritic.cells1', 'dendritic.cells2', 'dendritic.cells3', 'dendritic.cells4',
                   'dendritic.cells5', 'dendritic.cells6', 'dendritic.cells7', 'dendritic.cells8'], inplace=True, axis=1)
    gse64655 = gse64655.reset_index()
    gse64655 = pd.merge(gse64655, gse39652, on='genes')

    gse60424 = pd.read_csv('c:/data/Geo/RNASeq/gse60424.csv')
    gse60424 = pd.merge(gse60424, dict, on='Gene ID')
    gse60424.drop('Gene ID', inplace=True, axis=1)
    gse60424.set_index('genes', inplace=True)
    gse60424 = gse60424 * 1000000 / gse60424.sum()
    gse60424['CD4.T.cells'] = gse60424[['CD4.T.cells1', 'CD4.T.cells2', 'CD4.T.cells3', 'CD4.T.cells4']].min(axis=1)
    gse60424.drop(['CD4.T.cells1', 'CD4.T.cells2', 'CD4.T.cells3', 'CD4.T.cells4'], inplace=True, axis=1)
    gse60424['CD8.T.cells'] = gse60424[['CD8.T.cells1', 'CD8.T.cells2', 'CD8.T.cells3', 'CD8.T.cells4']].min(axis=1)
    gse60424.drop(['CD8.T.cells1', 'CD8.T.cells2', 'CD8.T.cells3', 'CD8.T.cells4'], inplace=True, axis=1)

    gse64655 = pd.merge(gse64655, gse60424, how='left', on='genes')
    gse64655['monocytes'] = gse64655[['monocytes1', 'monocytes2', 'monocytes3', 'monocytes4', 'monocytes5', 'monocytes6']].min(axis=1)
    gse64655.drop(['monocytes1', 'monocytes2', 'monocytes3', 'monocytes4', 'monocytes5', 'monocytes6'], inplace=True, axis=1)
    gse64655['neutrophils'] = gse64655[['neutrophils1', 'neutrophils2', 'neutrophils3', 'neutrophils4', 'neutrophils5', 'neutrophils6']].min(axis=1)
    gse64655.drop(['neutrophils1', 'neutrophils2', 'neutrophils3', 'neutrophils4', 'neutrophils5', 'neutrophils6'], inplace=True, axis=1)
    gse64655['B.cells'] = gse64655[['B.cells1', 'B.cells2', 'B.cells3', 'B.cells4', 'B.cells5', 'B.cells6']].min(axis=1)
    gse64655.drop(['B.cells1', 'B.cells2', 'B.cells3', 'B.cells4', 'B.cells5', 'B.cells6'], inplace=True, axis=1)
    gse64655['NK.cells'] = gse64655[['NK.cells1', 'NK.cells2', 'NK.cells3', 'NK.cells4', 'NK.cells5', 'NK.cells6']].min(axis=1)
    gse64655.drop(['NK.cells1', 'NK.cells2', 'NK.cells3', 'NK.cells4', 'NK.cells5', 'NK.cells6'], inplace=True, axis=1)
    # gse64655.columns = ['genes', 'dendritic.cells', 'monocytes1', 'neutrophils1', 'B.cells1', 'NK.cells1']
    gse64655 = gse64655.groupby('genes').median().reset_index()

    gse64655 = gse64655.loc[gse64655.max(axis=1) > 2]
    gse64655.to_csv('c:/data/GEO/RNASeq/rnaseq.csv')


# Create RNASeq matrix. Data is mostly from ABIS paper.
def create_lm_107011():
    base = 'c:/Users/Admin/Google Drive/src/GEO/RNASeq/mix-107011-107019-TPM/'
    data = pd.read_csv(base+'GSE107011_Processed_data_TPM.txt', delimiter='\t')
    genes = pd.read_csv(base+'genes2.csv')
    id_to_type = pd.read_csv(base+'107011-IDtoType.csv')
    id_to_type_dict = pd.read_csv(base+'107011-IDtoType.csv', index_col=0)['cell.type'].to_dict()
    id_to_name = pd.read_csv(base+'107011-IDtoName.csv', index_col=0)['sampleID'].to_dict()
    data = pd.merge(data, genes, on='geneID', how='left')
    data.drop('geneID', inplace=True, axis=1)
    data.dropna(inplace=True)
    data = data.groupby('gene').mean()
    data.columns = [id_to_name[i] for i in data.columns]
    data = data[list(id_to_type.sampleID)] #.T
    data = data.T.reset_index()
    data['cell.type'] = data['index'].apply(lambda x: id_to_type_dict[x])
    data = data.set_index('index')
    data = data.convert_objects(convert_numeric=True)
    data = data.groupby('cell.type').mean()
    data = data.T
    # data = data*1000000/data.sum()
    # data.to_csv('c:/data/GEO/RNAseq/pure_orig.csv')

    # Handle RNASeq fibroblast macrophages and endothelial
    # data = pd.read_csv('c:/data/GEO/RNASeq/pure_orig.csv', index_col=0).reset_index()
    # fib = pd.read_csv('c:/data/GEO/RNASeq/mix-115978/GSE115978_tpm.csv', index_col=0)[['cy82_CD45_neg_2_H12_S288_comb',
    #                                                                                   'Cy81_Bulk_CD45_G02_S170_comb',
    #                                                                                   'cy79_p3_CD45_neg_PDL1_neg_E09_S249_comb',
    #                                                                                   'cy53_1_CD45_neg_A03_S291_comb',
    #                                                                                   'cy82_CD45_pos_2_H12_S576_comb',
    #                                                                                   'cy80_CD45_neg_H12_S480_comb',
    #                                                                                   'Cy81_FNA_CD45_H12_S288_comb',
    #                                                                                   'CY88CD45POS_2_D12_S432_comb',
    #                                                                                   'cy81_Bulk_CD45_neg_E02_S146_comb',
    #                                                                                   'cy79_p5_CD45_neg_PDL1_neg_C04_S796_comb',
    #                                                                                   'CY89NEG_C11_S35_comb',
    #                                                                                   'cy94_cd45neg_cd90pos_D08_S332_comb']]
    fib = pd.read_csv(base+'GSE115978_tpm_reduced.csv', index_col=0)
    fib.columns = ['fib1', 'fib2', 'fib3', 'fib4', 'mac1', 'mac2', 'mac3', 'mac4', 'end1', 'end2', 'end3', 'end4']
    fib = (2**fib - 1) * 10 #Convert to TPM
    fib['endothelial'] = fib[['end1', 'end2', 'end3', 'end4']].mean(axis=1)
    fib['fibroblast'] = fib[['fib1', 'fib2', 'fib3', 'fib4']].mean(axis=1)
    fib['macrophages'] = fib[['mac1', 'mac2', 'mac3', 'mac4']].mean(axis=1)
    fib.drop(['end1', 'end2', 'end3', 'end4', 'fib1', 'fib2', 'fib3', 'fib4', 'mac1', 'mac2', 'mac3', 'mac4'], inplace=True, axis=1)
    fib = fib.reset_index()
    fib.columns = ['gene', 'endothelial', 'fibroblasts', 'macrophages']
    fib = pd.merge(fib, data, on='gene', how='left')
    # fib.to_csv('c:/data/GEO/RNASeq/pure_orig2.csv')

    # Handle RNASeq memory.CD4.T.cells.
    # cd4 = pd.read_csv('c:/data/GEO/RNASeq/73213-memoryCD4/memoryCD4.csv', index_col=0)
    # Change to TPM. Original data is single ended FPKM (RPKM).
    # cd4 = cd4 * 1000000 / cd4.sum()
    # cd4['memory.CD4.T.cells'] = cd4.mean(axis=1)
    # cd4 = cd4.reset_index()
    # cd4.drop(['4659', '5053', '5131', '5291', '4659_2', '5053_2', '5131_2', '5291_2'], inplace=True, axis=1)
    # fib = pd.merge(fib, cd4, on='gene', how='left')
    fib.dropna(inplace=True)
    fib = fib.set_index('gene')

    # Gene selection, we use quanTISeq gene list.
    # Gene diff.
    # gene_list_df, factor_list = gene_diff(fib)
    # gene_list = pd.concat([gene_list_df[celltype] for celltype in gene_list_df]).values
    # gene_list = list(set([x for x in gene_list if str(x) != 'nan']))
    # gene_df_1 = pd.read_csv('c:/data/Geo/RNASeq/quanTIseq/223180-4.csv') #quanTIseq
    gene_df_1 = pd.read_csv('c:/pure_orig8v2_14.csv')
    gene_df_2 = pd.read_csv('c:/ABIS.csv')
    gene_df_3 = pd.read_csv(base+'genes_memory_CD4.csv') #Xcell.
    # gene_df_3 = pd.read_csv('c:/data/Geo/RNASeq/Signatures/Kong_pone.0215987.s007.csv') #Kong
    # temp = pd.read_csv('c:/data/Geo/RNASeq/Signatures/xcell_orig.csv', index_col=0)
    # temp = np.concatenate(temp.values)
    # temp = set(list(temp))
    # gene_df_4 = {x for x in temp if x==x}
    #Understanding memory CD8+ T cells/Samji, Khanna
    memCD8 = set(['Eomes', 'ID3', 'BCL-6', 'STAT3', 'FOXO1', 'Serpina3g', 'TCF-7', 'Fos', 'ATF-2', 'BHLHE40', 'Fosb', 'Hopx', 'ID2', 'Klf4', 'Runx2', 'Cebpb', 'Myb', 'FasL', 'KLF2', 'Prdm1', 'ASCL2', 'RUNX1', 'LEF1', 'FOXP1', 'KLF12', 'ZNF511', 'IKZF5', 'GTF2H2', 'RERE', 'KLF7', 'ZNF318', 'TFDP2', 'BBX', 'HOXB7', 'ZNF274', 'ZNF638', 'ZNF365', 'EWSR1', 'ZNF18', 'CNOT7', 'ZNF22', 'ZNF428', 'FOXN3', 'ZNF83', 'GTF3A', 'Litaf', 'Nr4a1', 'Nr4a2', 'KLF2', 'Zfp683', 'ETS1', 'ZFP36L1', 'GPBP1', 'NFE2L2', 'RELA', 'FOS', 'NFAT5', 'ZC3H12A', 'NOTCH1', 'MYB', 'KLF5', 'CITED2', 'NR3C1', 'REL', 'ZNF32', 'AHR', 'TWIST1', 'EPAS1', 'KLF6', 'NFKB1', 'CREM', 'NR4A3', 'NFKB2', 'BATF', 'IRF4', 'RBPJ', 'PBX4', 'ATF3', 'FOSB', 'BACH2', 'RUNX3', 'EGR2'])
    ##ssGsea
    #CD4 = set(['ANP32B', 'ASF1A', 'ATF2', 'BATF', 'C13orf34', 'CD28', 'DDX50', 'FAM111A', 'FRYL', 'GOLGA8A', 'ICOS', 'ITM2A', 'LRBA', 'NAP1L4', 'NUP107', 'PHF10', 'PPP2R5C', 'RPA1', 'SEC24C', 'SLC25A12', 'SRSF10', 'TRA', 'UBE2L3', 'YME1L1', 'AQP3', 'ATF7IP', 'ATM', 'CASP8', 'CDC14A', 'CEP68', 'CLUAP1', 'CREBZF', 'CYLD', 'DOCK9', 'FAM153B', 'FOXP1', 'FYB', 'HNRPH1', 'INPP4B', 'KLF12', 'LOC441155', 'MAP3K1', 'MLL', 'N4BP2L2-IT2', 'NEFL', 'NFATC3', 'PCM1', 'PCNX', 'PDXDC2', 'PHC3', 'POLR2J2', 'PSPC1', 'REPS1', 'RPP38', 'SLC7A6', 'SNRPN', 'ST3GAL1', 'STX16', 'TIMM8A', 'TRAF3IP3', 'TXK', 'TXLNGY', 'USP9Y', 'AKT3', 'C7orf54', 'CCR2', 'DDX17', 'EWSR1', 'FLI1', 'GDPD5', 'LTK', 'MEFV', 'NFATC4', 'PRKY', 'TBC1D5', 'TBCD', 'TRA', 'VIL2', 'B3GAT1', 'BLR1', 'C18orf1', 'CDK5R1', 'CHGB', 'CHI3L2', 'CXCL13', 'HEY1', 'HIST1H4K', 'ICA1', 'KCNK5', 'KIAA1324', 'MAF', 'MAGEH1', 'MKL2', 'MYO6', 'MYO7A', 'PASK', 'PDCD1', 'POMT1', 'PTPN13', 'PVALB', 'SH3TC1', 'SIRPG', 'SLC7A10', 'SMAD1', 'ST8SIA1', 'STK39', 'THADA', 'TOX', 'TSHR', 'ZNF764', 'C1orf61', 'CD160', 'FEZ1', 'TARP', 'TRD', 'TRGV9', 'APBB2', 'APOD', 'ATP9A', 'BST2', 'BTG3', 'CCL4', 'CD38', 'CD70', 'CMAH', 'CSF2', 'CTLA4', 'DGKI', 'DOK5', 'DPP4', 'DUSP5', 'EGFL6', 'GGT1', 'HBEGF', 'IFNG', 'IL12RB2', 'IL22', 'LRP8', 'LRRN3', 'LTA', 'SGCB', 'SYNGR3', 'ZBTB32', 'IL17A', 'IL17RA', 'RORC', 'ADCY1', 'AHI1', 'ANK1', 'BIRC5', 'CDC25C', 'CDC7', 'CENPF', 'CXCR6', 'DHFR', 'EVI5', 'GATA3', 'GSTA4', 'HELLS', 'IL26', 'LAIR2', 'LIMA1', 'MB', 'MICAL2', 'NEIL3', 'PHEX', 'PMCH', 'PTGIS', 'SLC39A14', 'SMAD2', 'SNRPD1', 'WDHD1', 'FOXP3'])
    CD4 = set(['Abracl','Il2rb','Aldoa','Ablim1','Ccnd2','Btg1','Arl6ip1','Acot7','Itga4','Angptl2','Acot7','Eif4a2','Capg','AW112010','Actb','Itgb2','Asap1','Acp5','Gm10073','Ccnd2','Bhlhe40','Actg1','Itgb7','Batf','Actn1','Il7r','Cd2','Ccl3','Actr3','Klrc1','BC021614','Adgre5','Malat1','Cd74','Ccl4','Adgre5','Krtcap2','Bcl2a1b','Adk','Rabac1','Ctla4','Ccl5','Agpat4','Lfng','Borcs8','Ahnak','Rpl12','Foxp3','Ccr5','Ahnak','Lgals1','Cd160','Anp32a','Rpl13','Gbp7','Cd160','Anxa1','Lgals3','Cd200','Anxa2','Rpl13-ps3','Gimap7','Cd27','Anxa2','Lsp1','Cd3g','Arl5c','Rpl18a','Gpx4','Cd3g','Anxa5','Ly6c2','Cd82','Atp1b3','Rpl36','H2-T22','Cd7','Anxa6','Mrpl33','Cox14','Atp2b1','Rpl36a','Ifngr1','Cst7','Ap2s1','Ms4a4b','Ctsb','Bcl2','Rpl37a','Ikzf2','Ctla2a','Aprt','Ms4a6b','Cxcr5','Bin2','Rpl38','Il7r','Ctla4','Arhgdib','Mtpn','Ddit4','Ccdc28b','Rpl9','Izumo1r','Ctsb','Arl6ip5','Myl12a','Dennd2d','Ccr7','Rpl9-ps6','Ltb','Ctsw','Arpc5','Myl6','Eea1','Cd9','Rplp1','Mbnl1','Cxcr6','Atp1a1','Myo1f','Fam162a','Cdc25b','Rplp2','Peli1','Cyba','Atp5c1','Myo1g','Fyn','Cdc42se1','Rps14','Rgs1','Dnaja1','Atp5h','Ndufb7','Gapdh','Cdkn2d','Rps15','Samsn1','Efhd2','Atp5l','Nkg7','Gdi2','Crip1','Rps15a','Sdf4','Fasl','AW112010','Nptn','Gimap5','Crip2','Rps16','Sell','Fyn','Bhlhe40','Ostf1','Gna13','Dock2','Rps19','Serinc3','Gadd45b','Calm1','Pglyrp1','Gng2','Emb','Rps21','Shisa5','Gimap7','Capzb','Plac8','Hif1a','Ezr','Rps23','Tnfrsf18','Gpr65','Ccl5','Plek','Hmgb1','Fam177a','Rps24','Tnfrsf4','Gzmk','Ccr2','Podnl1','Icos','Fam65b','Rps27','H2-K1','Cd2','Ppib','Ifi27l2a','Glud1','Rps27rt','Hspa1a','Cd47','Ppp1ca','Isg15','Hint1','Rps28','Hspa8','Cd48','Ppp3ca','Izumo1r','Id3','Rps29','Hsph1','Cd52','Prelid1','Limd2','Il6ra','Rps5','Id2','Cdc42','Psmb3','Lpp','Iqgap1','Rps9','Lag3','Cdc42ep3','Pycard','Lrmp','Itgb1','Lilr4b','Cfl1','Rap1b','Maf','Itpkb','Nkg7','Clic1','Rasgrp2','Marcksl1','Jund','Nr4a2','Clta','Rbm3','Matk','Klf2','Pdcd1','Coro1a','Reep5','Mif','Klf3','Plac8','Cox17','Rnf138','Mmd','Klf6','Psmb10','Cox5a','Rnf166','Nsg2','mt-Co1','Psmb8','Cox5b','Rora','Nt5e','Pbxip1','Rgs1','Crip1','Rpa2','P2rx7','Pdlim1','Rnaset2b','Crot','Runx3','Pdcd1','Prr13','Samsn1','Ctla2a','S100a10','Pfkl','Psma6','Serpina3g','Ctsc','S100a11','Pfkp','Raf1','Sh2d2a','Ctsd','S100a13','Pgam1','Rasa3','Shisa5','Ctsw','S100a4','Pkm','Rasgrp2','Slamf7','Cxcr6','S100a6','Ppp1r14b','Rora','Sub1','Cyba','S1pr4','Prkca','Rpa1','Tapbpl','Cyth4','Sec61b','Ptp4a2','S100a10','Tnfrsf1b','Dnajc15','Sec61g','Ptpn11','S100a11','Traf4','Dok2','Selplg','Ptprcap','S100a4','Ubb','Ech1','Sept11','Ptrh1','S100a6','Emp3','Serpinb6b','Rab37','Sfr1','Eno1','Serpinb9','Rgs10','Slamf6','Epsti1','Sh3bgrl3','Rnaset2a','Spn','Esyt1','Slamf1','Rnaset2b','Srpk1','Fam107b','Sms','Rpsa','Stk24','Gabarapl2','Spn','Scd2','Stk38','Gapdh','St3gal6','Sept7','Tagap','Gbp7','Stmn1','Sh2d1a','Tagln2','Ggh','Sub1','Smco4','Tcf7','Gimap7','Tagln2','Sostdc1','Tcp11l2','Glipr2','Taldo1','Tbc1d4','Tspan13','Glrx','Tbx21','Tnfaip8','Tuba1a','Gm4950','Tceb2','Tnfsf8','Vim','Gmfg','Thy1','Tox','Vsir','Gna15','Tmed2','Tox2','Xrn2','Gramd3','Tmsb4x','Tpi1','Gzmb','Tpm4','Trim8','H2afy','Tpst2','Zap70','H2afz','Tspo','Zfp36l1','Hmgb2','Txn1','Hsp90b1','Txndc5','Id2','Ube2g2','Idh3a','Vim','Ifngr1','Zyx','Il18r1'])
    macrophages = set(['10-SEP', 'ABCD1', 'ABI1', 'ABTB2', 'ACADVL', 'ACP2', 'ACSM5', 'ACTR10', 'ACTR2', 'ACTR3', 'ADAMDEC1', 'ADCK2',
                    'ADCY3', 'ADO', 'ADRA2B', 'AFG3L2', 'AGGF1', 'AGPS', 'AKR7A2', 'ALCAM', 'ALDH9A1', 'ALG9', 'ALK', 'ANGPT4',
                    'ANKFY1', 'ANXA11', 'ANXA2', 'AP1B1', 'AP1M2', 'AQP8', 'ARFGEF2', 'ARHGEF11', 'ARL8B', 'ARPC4', 'ARSB', 'ATOX1',
                    'ATP2A2', 'ATP2C1', 'ATP6AP2', 'ATP6V0A1', 'ATP6V0C', 'ATP6V0D1', 'ATP6V0E1', 'ATP6V1A', 'ATP6V1C1', 'ATP6V1D',
                    'ATP6V1E1', 'ATP6V1F', 'ATP6V1H', 'BAG3', 'BAIAP2', 'BCAP31', 'BCKDK', 'BLVRA', 'BPI', 'BTBD1', 'C10orf76',
                    'C12orf4', 'C12orf49', 'C16orf62', 'C1QA', 'C1QB', 'C3AR1', 'C7orf25', 'CALR', 'CAMP', 'CANX', 'CARD14', 'CCDC47',
                    'CCDC85C', 'CCDC88A', 'CCL1', 'CCL18', 'CCL19', 'CCL22', 'CCL24', 'CCL7', 'CCL8', 'CCR1', 'CD163', 'CD164',
                    'CD300C', 'CD48', 'CD52', 'CD63', 'CD80', 'CD81', 'CD84', 'CD9', 'CDS2', 'CECR5', 'CEPT1', 'CETN2', 'CHD9',
                    'CHIT1', 'CIAO1', 'CIR1', 'CLCN7', 'CLEC4E', 'CLEC5A', 'CLIP1', 'CLPB', 'CLTC', 'CMKLR1', 'CNIH4', 'COL4A3BP',
                    'COMMD8', 'COMMD9', 'COQ2', 'CORO7', 'COX15', 'COX5A', 'COX5B', 'COX7B', 'COX8A', 'CPNE6', 'CRYBB1', 'CSF1',
                    'CSF1R', 'CXCL9', 'CYBA', 'CYBB', 'CYC1', 'CYFIP1', 'CYP19A1', 'DAGLA', 'DBI', 'DERA', 'DHX57', 'DLAT', 'DNAJC13',
                    'DNASE1L3', 'DNASE2B', 'DOT1L', 'ECHS1', 'EFR3A', 'ELK1', 'ELOVL1', 'EMILIN1', 'ERP29', 'EXOC1', 'EXOC5',
                    'FAM32A', 'FANCE', 'FCER1G', 'FDX1', 'FEZ2', 'FGR', 'FH', 'FKBP15', 'FLT1', 'FOLR2', 'FPR2', 'FPR3', 'FTL',
                    'G6PC3', 'GABARAP', 'GGA1', 'GLB1', 'GLRX2', 'GORASP1', 'GP1BA', 'GPD1', 'GRB2', 'GSTO1', 'GUCA1A', 'GUF1',
                    'HADHB', 'HAMP', 'HAUS2', 'HCCS', 'HEXA', 'HEXB', 'HIGD2A', 'HK3', 'HMGCL', 'HPS1', 'HS3ST2', 'HSD17B12', 'HSPB7',
                    'HSPH1', 'HTT', 'HYAL2', 'IARS2', 'IBTK', 'IFNAR1', 'IGSF6', 'IL10', 'IL12B', 'IL17RA', 'IPPK', 'ITGAE', 'ITGAX',
                    'ITGB1BP1', 'KCMF1', 'KCNJ1', 'KCNJ5', 'KCNK13', 'KCNMB1', 'KCTD5', 'KIAA0196', 'KIFC3', 'KLHL12', 'LAIR1',
                    'LAMP1', 'LDHAL6B', 'LILRA2', 'LILRB1', 'LILRB4', 'LILRB5', 'LIMD2', 'LONP1', 'LONRF3', 'LY86', 'M6PR', 'MAPK13',
                    'MAPKAP1', 'MARCO', 'MDH1', 'MFN1', 'MFSD7', 'MGST3', 'MKL2', 'MLX', 'MMP19', 'MMP8', 'MRM1', 'MRPL12', 'MRPL40',
                    'MRS2', 'MS4A4A', 'MS4A6A', 'MSR1', 'MT2A', 'MTHFR', 'MTMR14', 'MUL1', 'MYBPH', 'MYH11', 'MYO15A', 'MYO7A',
                    'MYO9B', 'MYOF', 'MYOZ1', 'NAGPA', 'NARS', 'NCAPH', 'NCKAP1L', 'NDUFA8', 'NDUFAF1', 'NDUFB1', 'NDUFB3', 'NDUFB6',
                    'NDUFS2', 'NDUFS3', 'NDUFS6', 'NDUFS8', 'NECAP2', 'NFS1', 'NOP10', 'NPR1', 'NPTN', 'NRBP1', 'NSMAF', 'NSUN3',
                    'NUBP1', 'NUDT9', 'NUMB', 'OGFR', 'ORMDL2', 'OS9', 'OSBPL11', 'OTUD4', 'P2RX7', 'PABPC4', 'PANK3', 'PCDHB11',
                    'PCMT1', 'PDCD6IP', 'PDCL', 'PDE1B', 'PEX14', 'PEX19', 'PHLDB1', 'PICK1', 'PKD2L1', 'PLEKHB2', 'PLEKHM2',
                    'PMFBP1', 'POGK', 'PPCS', 'PQLC2', 'PRDX1', 'PRDX3', 'PSMD10', 'PSME1', 'PTGIR', 'PTPN12', 'PTPRA', 'QDPR',
                    'RAB1A', 'RAB3IL1', 'RAB5C', 'RABGGTA', 'RAC1', 'RALA', 'RB1', 'RELA', 'RENBP', 'RIN2', 'RNH1', 'RRP1', 'RTN4',
                    'S100A11', 'S100A6', 'S1PR2', 'SCAMP2', 'SDCBP', 'SDHB', 'SDHD', 'SDS', 'SETD3', 'SH3GLB1', 'SIGLEC1', 'SIGLEC7',
                    'SIGLEC9', 'SLAMF8', 'SLC11A1', 'SLC1A2', 'SLC25A24', 'SLC25A46', 'SLC30A5', 'SLC31A1', 'SLC38A7', 'SLC39A1',
                    'SLC6A12', 'SLC6A7', 'SLC9A6', 'SMG5', 'SNAPC2', 'SNUPN', 'SNX1', 'SNX2', 'SNX3', 'SNX4', 'SNX5', 'SPG21',
                    'SPIN1', 'SPR', 'SRC', 'STAB1', 'STAM2', 'STIP1', 'STX12', 'STX18', 'STX4', 'STYXL1', 'SUMO3', 'TAF10',
                    'TBC1D16', 'TBC1D9B', 'TCEAL4', 'TCEB1', 'TDRD7', 'TFEC', 'TFRC', 'TGOLN2', 'TIE1', 'TM2D1', 'TMBIM4', 'TMED5',
                    'TMEM115', 'TMEM126B', 'TMEM127', 'TMEM147', 'TMEM184C', 'TMEM33', 'TMEM70', 'TMEM9B', 'TMX1', 'TNFRSF12A',
                    'TNFSF14', 'TNPO1', 'TPP1', 'TRAF3', 'TRAPPC2L', 'TRAPPC3', 'TREM2', 'TRIP4', 'TSPO', 'TTLL4', 'TULP4',
                    'TYROBP', 'UBE2D4', 'UBXN6', 'UCP3', 'UGP2', 'UNC50', 'UQCR10', 'UQCR11', 'UQCRC2', 'USF2', 'USP14', 'UTP3',
                    'VAMP8', 'VIM', 'VPS33A', 'VPS35', 'VPS53', 'VSIG4', 'VTI1B', 'WDFY3', 'WDR11', 'WSB2', 'WTAP', 'XPNPEP2',
                    'YIF1B', 'ZC3H15', 'ZC3H3', 'ZCCHC4', 'ZDHHC24', 'ZDHHC3', 'ZMPSTE24', 'ZNF219', 'ZZZ3'])
    macrophages2 = set(['CD14', 'CD16', 'CD64', 'CD68', 'CD71', 'CD86', 'CD80', 'CD68', 'MHCII', 'IL', '1R', 'TLR2', 'TLR4', 'iNOS', 'SOCS3', 'CD163', 'MHCII', 'SR', 'MMR', 'CD206', 'CD200R', 'TGM2', 'DecoyR', 'IL', '1R ', 'I', 'Mouse ', 'nly', ' ', 'm1', '2', 'Fizz1', 'Arg', '1', 'CD86', 'MHCII', 'CD163', 'TLR1', 'TLR8', 'VEGF', 'CCL2', 'CD3', 'CD31', 'CD68', 'CD163', 'HLA', 'DR', 'IL', '10', 'iNOS', 'PCNA', 'VEGF', 'CD68', 'Fyn', 'LAT', 'Lck', 'ZAP70', 'CD11b', 'CD11c', 'CD68', 'F4', '80', 'MHCII'])
    #fib = fib.reindex(list(set(gene_df_1.gene.values) | set(gene_df_2.GeneSymbol) | memCD8)) #| set(gene_df_3.gene.values) 
    #fib = fib.reindex(list(set(gene_df_1.gene.values) | set(gene_df_2.GeneSymbol) | macrophages2))
    fib.dropna(inplace=True)
    fib.to_csv(base+'pure_orig8v6.csv')
    
    # Handle RNASeq Macrophages.
    # gse39652 = pd.read_csv('c:/data/Geo/RNASeq/quanTISeq/GSE36952.csv', index_col=0)
    # gse39652 = gse39652 * 1000000 / gse39652.sum()
    # gse39652['Macrophages'] = gse39652.mean(axis=1)
    # gse39652.drop(['M1_1', 'M1_2', 'M1_3', 'M2_1', 'M2_2', 'M2_3'], axis=1, inplace=True)
    # gse39652 = gse39652.reset_index()
    # gse39652.columns = ['gene', 'macrophages']
    # data = pd.merge(data, gse39652, on='gene', how='left')
    # data.dropna(inplace=True)
    # data.to_csv('c:/data/GEO/RNASeq/pure_orig4.csv')

# Create an RNASeq mix from a pure matrix.
def create_rnaseq_mix_115978():
    # data = pd.read_csv('c:/data/GEO/RNASeq/pure_orig.csv', index_col=0)
    annotations = pd.read_csv('c:/data/GEO/RNASeq/mix-115978/GSE115978_cell.annotations.csv')
    cell_type = ['Endo.', 'T.CD4', 'CAF', 'T.CD8', 'NK', 'B.cell']
    cell_type2 = ['endothelial.cells', 'CD4.T.cells', 'fibroblasts', 'CD8.T.cells', 'NK.cells', 'B.cells']
    # cell_type = ['T.CD4', 'T.CD8', 'NK', 'B.cell', 'Macrophage']
    cols = []
    for i in cell_type:
        annot_temp = annotations.loc[annotations['cell.types'] == i]
        row = annot_temp.iloc[0] #randrange(len(annot_temp))
        cols += [row['cells']]
    data = pd.read_csv('c:/data/GEO/RNASeq/mix-115978/GSE115978_tpm.csv', index_col=0)[cols]
    data = (2**data-1)*10 #Convert to TPM
    data.columns = cell_type2
    data2 = data.copy()
    prop_data = np.empty((0, len(data.columns)))
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = 1 #np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2.to_csv('c:/data/GEO/RNASeq/mix-115978/mix.csv')
    data.to_csv('c:/data/GEO/RNASeq/mix-115978/pure.csv')
    prop_df.index = cell_type2
    prop_df.to_csv('c:/data/GEO/RNASeq/mix-115978/prop.csv')


# Create an RNASeq mix from a pure matrix.
def create_rnaseq_mix_115736():
    cell_type = ['dendritic', 'macrophages', 'CD4.T.cells', 'CD8.T.cells', 'naive.B.cells', 'NK.cells'] #, 'memory.B.cells']
    data = pd.read_csv('c:/data/Geo/RNASeq/mix-115736/GSE115736_Haemopedia-Human-RNASeq_rpkm.txt', delimiter='\t', index_col = 0)[cell_type]
    data.reset_index(inplace=True)
    data.rename({'index':'geneID'}, inplace=True, axis=1)
    genes = pd.read_csv('c:/data/GEO/RNASeq/genes2.csv')
    genes.geneID = genes.geneID.apply(lambda x: x.split('.')[0])
    genes = pd.merge(genes, data, on='geneID', how='left')
    genes.drop('geneID', inplace=True, axis=1)
    genes.dropna(inplace=True)
    genes = genes.groupby('gene').mean()
    # Tranfer to TPM from RPKM.
    genes = genes * 1000000 / genes.sum()
    data = genes
    data2 = data.copy()
    prop_data = np.empty((0, len(data.columns)))
    # Loop on mixes.
    for i in range(0,20):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = 1 #np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,20)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,20)]]
    data2.to_csv('c:/data/GEO/RNASeq/mix-115736/mix2.csv')
    data.to_csv('c:/data/GEO/RNASeq/mix-115736/pure2.csv')
    prop_df.index = cell_type
    prop_df.to_csv('c:/data/GEO/RNASeq/mix-115736/prop2.csv')


# Create an RNASeq mix from a pure matrix.
def create_rnaseq_mix_118165():
    cell_type = ['memory.CD8.T.cells', 'naive.CD8.T.cells', 'myeloid.dendritic.cells', 'monocytes', 'naive.B.cells', 'NK.cells', 'memory.B.cells', 'regulatory.T.cells']
    data = pd.read_csv('c:/Work/data/GEO/RNASeq/mix-118165-tximport/GSE118165_RNA_gene_abundance.csv')
    gencode = pd.read_csv('c:/Work/data/GEO/RNASeq/mix-118165-tximport/gencode.csv')
    gencode['ENS_ID'] = gencode.ENSID.apply(lambda x: x.split('.')[0])
    gencode = gencode[['ENS_ID', 'ID', 'length']]
    gencode['length'] = gencode.length/1000
    data = pd.merge(data, gencode, on='ENS_ID', how='left')
    data.dropna(inplace=True)
    data.drop(['ENS_ID'], axis=1, inplace=True)
    data.drop_duplicates('ID', keep='first', inplace=True)
    data = data.set_index('ID')
    data = data.div(data.length, axis=0)
    data.drop(['length'], axis=1, inplace=True)
    # Convert to TPM
    data = (data*1000000).div(data.sum(axis=0).values, axis=1)
    data.to_csv('c:/Work/data/GEO/RNASeq/mix-118165-tximport/abundance/pure.csv')

    data2 = data.copy()
    prop_data = np.empty((0, len(data.columns)))
    # Loop on mixes.
    for i in range(0,20):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = 1 #np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,20)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,20)]]
    data2.to_csv('c:/Work/data/GEO/RNASeq/mix-118165-tximport/abundance/mix.csv')
    prop_df.index = cell_type
    prop_df.to_csv('c:/Work/data/GEO/RNASeq/mix-118165-tximport/abundance/prop.csv')


# Create microarray mix from GSE22886, the file Tom sent.
def create_microarray_mix():
    cells_ = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytic.lineage']
    data = pd.read_csv('c:/data/GEO/mix1-22886/GSE22886_Microarray_with_GSM.csv', index_col=0)
    data2 = data[['CD4.T.cells', 'CD4.T.cells.1']].copy()
    prop_data = np.empty((0, len(cells_)))
    # Loop on mixes.
    for i in range(0,20):
        prop = np.random.dirichlet(np.ones(len(cells_)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in cells_:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,20)],
                           index=cells_)
    data2 = data2[[f'mix{i}' for i in range(0,20)]]
    data2 = data2.reset_index()
    data2.rename(columns={'index':'ID_REF'}, inplace=True)
    u133a = pd.read_csv('c:/data/GEO/U133A.csv', index_col=0)
    data2 = pd.merge(data2, u133a, on='ID_REF', how='left')
    data2.dropna(inplace=True)
    data2.drop('ID_REF', inplace=True, axis=1)
    data2 = data2.groupby('Symbol').mean()
    data2.to_csv('c:/data/GEO/mix1-22886/mix1_array.csv')
    prop_df.to_csv('c:/data/GEO/mix1-22886/prop1_array.csv')


def create_rnaseq_mix_107019():
    cells_ = ['naive.B.cells', 'memory.B.cells', 'naive.CD4.T.cells', 'naive.CD8.T.cells', 'memory.CD8.T.cells', 'regulatory.T.cells', 'monocytes', 'NK.cells', 'myeloid.dendritic.cells', 'neutrophils']
    data = pd.read_csv('c:/data/Geo/RNASeq/mix-107019/GSE107011_Processed_data_TPM.csv')[['geneID'] + cells_]
    genes = pd.read_csv('c:/data/Geo/RNASeq/mix-107019/genes2.csv')
    data = pd.merge(data, genes, on='geneID', how='left')
    data.drop('geneID', inplace=True, axis=1)
    data.dropna(inplace=True)
    data = data.groupby('gene').mean()
    prop_data = np.empty((0, len(cells_)))
    data2 = data[['naive.CD4.T.cells', 'naive.CD8.T.cells']].copy()
    # Loop on mixes.
    for i in range(0,20):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,20)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,20)]]
    data2.to_csv('c:/data/GEO/RNASeq/mix-107019/mix2/mix2.csv')
    data.to_csv('c:/data/GEO/RNASeq/mix-107019/mix2/pure2.csv')
    prop_df.to_csv('c:/data/GEO/RNASeq/mix-107019/mix2/prop2.csv')


def create_rnaseq_mix_112101():
    cells_ = ['B.cells', 'CD4.T.cells', 'monocytic.lineage', 'endothelial.cells', 'neutrophils', 'fibroblasts']
    data = pd.read_csv('c:/data/Geo/RNASeq/mix-112101/GSE112101_NormalizedReadCounts.csv', index_col=0)
    prop_data = np.empty((0, len(cells_)))
    data2 = data.copy()
    # Loop on mixes.
    for i in range(0,20):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,20)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,20)]]
    data2.to_csv('c:/data/GEO/RNASeq/mix-112101/mix.csv')
    data.to_csv('c:/data/GEO/RNASeq/mix-112101/pure.csv')
    prop_df.to_csv('c:/data/GEO/RNASeq/mix-112101/prop.csv')


def create_microarray_mix_109348_3982_28490_28492_58173_72642_93777():
    global_table = pd.DataFrame()

    # Download GSE.
    gse = GEOparse.get_GEO(geo='GSE93777', destdir="c:/data/Geo/Microarray/GSE")
    gse2 = GEOparse.get_GEO(geo='GSE82159', destdir="c:/data/Geo/Microarray/GSE")
    # Parse platform data.
    for gpl_name, gpl in gse.gpls.items():
        # print("Metadata:",)
        # for key, value in gpl.metadata.items():
        #    print(" - %s : %s" % (key, ", ".join(value)))
        gpl_table = gpl.table[['ID', 'Gene Symbol']]
        # Split rows that have probes mapped to multiple genes: probex,geneA///geneB///geneC
        gpl_table = gpl_table.set_index(['ID'])['Gene Symbol'].str.split('\/\/\/', expand=True).stack().reset_index()
        gpl_table.drop(['level_1'], inplace=True, axis=1)
        gpl_table.columns = ['ID_REF', 'Symbol']
        gpl_table['Symbol'] = gpl_table.Symbol.str.strip()
        break

    # Loop on all cell types
    # celltypes = {'GSM2940351':'memory.CD4.T.cells',
    #              'GSM2940345': 'myeloid.dendritic.cells',
    #              'GSM2940349': 'naive.CD4.T.cells',
    #              'GSM2940358': 'neutrophils',
    #              'GSM2940348': 'NK.cells',
    #              'GSM2940355': 'naive.B.cells',
    #              'GSM2940353': 'naive.CD8.T.cells',
    #              'GSM2940341': 'monocytes'}
    # celltypes = {'GSM90838': 'macrophages', 'GSM90854': 'memory.CD4.T.cells', 'GSM90666': 'myeloid.dendritic.cells', 'GSM90845': 'B.cells', 'GSM90844': 'neutrophils', 'GSM90851': 'NK.cells'}
    celltypes_28490 = {'GSM705326': 'neutrophils', 'GSM705299': 'B.cells', 'GSM705302': 'CD4.T.cells', 'GSM705312': 'CD8.T.cells', 'GSM705287': 'monocytic.lineage', 'GSM705307': 'NK.cells'} #'GSM705321': 'myeloid.dendritic.cells',
    celltypes_28492 = {'GSM705297': 'B.cells', 'GSM705302': 'CD4.T.cells', 'GSM705312': 'CD8.T.cells', 'GSM705287': 'monocytic.lineage', 'GSM705326': 'neutrophils', 'GSM705309': 'NK.cells'} #'GSM705321': 'myeloid.dendritic.cells' 'GSM705420'
    celltypes_28492 = {'GSM705299': 'B.cells', 'GSM705302': 'CD4.T.cells', 'GSM705312': 'CD8.T.cells', 'GSM705287': 'monocytic.lineage', 'GSM705326': 'neutrophils', 'GSM705307': 'NK.cells'} #'GSM705321': 'myeloid.dendritic.cells' 'GSM705420'
    celltypes_58173 = {'GSM1402765': 'B.cells', 'GSM1402752': 'CD4.T.cells', 'GSM1402755': 'CD8.T.cells', 'GSM1402758': 'monocytic.lineage', 'GSM1402761': 'neutrophils', 'GSM1402768': 'NK.cells'}

    celltypes_72642 = {'GSM1867065': 'B.cells', 'GSM1867066': 'CD4.T.cells', 'GSM1867067': 'CD8.T.cells', 'GSM1867068': 'monocytic.lineage', 'GSM1867070': 'neutrophils', 'GSM1867069': 'NK.cells'}
    celltypes_72642c = {'GSM1867065': 'B.cells', 'GSM1867067': 'CD8.T.cells', 'GSM1867068': 'monocytic.lineage', 'GSM1867070': 'neutrophils', 'GSM1867069': 'NK.cells'}
    celltypes_110085 = {'GSM2977380': 'endothelial.cells', 'GSM2977390': 'fibroblasts'}

    celltypes_93777 = {'GSM2461874': 'CD4.T.cells', 'GSM2461949': 'CD8.T.cells', 'GSM2461958': 'monocytic.lineage', 'GSM2462007': 'neutrophils', 'GSM2462022': 'NK.cells'}
    celltypes_82159 = {'GSM2185120' : 'B.cells', 'GSM2185129': 'fibroblasts'}
    celltypes = celltypes_93777

    # Parse GSM examples.
    for gsm_name, gsm in gse.gsms.items():
        if gsm_name in celltypes.keys():
            print("Name: ", gsm_name)
            # print("Metadata:",)
            # for key, value in gsm.metadata.items():
            #    print(" - %s : %s" % (key, ", ".join(value)))
            table = gsm.table
            table = pd.merge(table, gpl_table, on='ID_REF', how='left')
            table.drop('ID_REF', axis=1, inplace=True)
            table.dropna(inplace=True)
            table = table.groupby('Symbol').VALUE.min()
            table.rename(celltypes[gsm_name], inplace=True)
            global_table = pd.concat([global_table, table], axis=1)

    for gsm_name, gsm in gse2.gsms.items():
        if gsm_name in celltypes_82159.keys():
            print("Name: ", gsm_name)
            # print("Metadata:",)
            # for key, value in gsm.metadata.items():
            #    print(" - %s : %s" % (key, ", ".join(value)))
            table = gsm.table
            table = pd.merge(table, gpl_table, on='ID_REF', how='left')
            table.drop('ID_REF', axis=1, inplace=True)
            table.dropna(inplace=True)
            table = table.groupby('Symbol').VALUE.min()
            table.rename(celltypes_82159[gsm_name], inplace=True)
            global_table = pd.concat([global_table, table], axis=1)

    # Anti-log if needed.
    for col in global_table:
        if global_table[col].max() < 20:
            global_table[col] = 2 ** global_table[col]

    global_table.dropna(inplace=True)
    prop_data = np.empty((0, len(celltypes)+len(celltypes_82159)))
    data2 = global_table.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(len(celltypes)+len(celltypes_82159)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in global_table:
            noise = np.random.normal(1, 0.25, len(global_table[col]))
            data2[f'mix{i}'] += global_table[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=global_table.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2.to_csv('c:/data/Geo/Microarray/mix-93777/mix.csv')
    prop_df.to_csv('c:/data/Geo/Microarray/mix-93777/prop.csv')
    global_table.to_csv(f'c:/data/Geo/Microarray/mix-93777/pure.csv')


def create_microarray_mix_40240():
    data = pd.read_csv('c:/data/Geo/Microarray/mix-40240/GSE40240_series_matrix.csv')
    gpl = pd.read_csv('c:/data/Geo/Microarray/mix-40240/GPL6244-17930.csv')
    gpl = gpl.set_index(['ID'])['gene_assignment'].str.split('\/\/', expand=True)[1].reset_index()
    gpl.columns = ['ID_REF', 'gene']
    gpl.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl.dropna(inplace=True)
    data['ID_REF'] = data.ID_REF.astype(np.int64)
    data = pd.merge(data, gpl, on='ID_REF', how='left')
    data.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    data.dropna(inplace=True)
    data.drop('ID_REF', axis=1, inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    data.to_csv('c:/data/Geo/Microarray/mix-40240/mix.csv')


def create_microarray_mix_77343():
    global_table = pd.DataFrame()

    # Download GSE.
    gse = GEOparse.get_GEO(geo='GSE77343', destdir="c:/data/Geo/Microarray/GSE")

    # Parse platform data.
    gpl = gse.gpls['GPL11532']
    gpl_table = gpl.table[['ID', 'gene_assignment']].copy()
    gpl_table['gene_assignment'] = gpl_table.gene_assignment.str.split('\/\/', expand=True)[1] #.reset_index()
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    prop_data = pd.DataFrame(columns=['neutrophils', 'monocyte'])
    k = 0
    for gsm_name, gsm in gse.gsms.items():
        # test = gse.gsms['GSM2049464']
        prop_data.loc[k] = [float(gsm.metadata['characteristics_ch1'][3][23:]), float(gsm.metadata['characteristics_ch1'][4][22:])]
        k += 1
        temp = gsm.table.set_index('ID_REF')
        global_table = pd.concat([global_table, temp], axis=1)

    global_table.columns = [f'mix-{i}' for i in range(1,198)]
    global_table = global_table.reset_index()
    global_table = pd.merge(global_table, gpl_table, on='ID_REF', how='left')
    global_table.drop('ID_REF', axis=1, inplace=True)
    global_table.dropna(inplace=True)
    global_table['gene'] = global_table.gene.str.strip()
    global_table = global_table.groupby('gene').agg(np.mean)
    prop_data.T.to_csv('c:/data/Geo/Microarray/mix-77343/prop.csv')
    global_table.to_csv('c:/data/Geo/Microarray/mix-77343/mix.csv')


def create_microarray_mix_44621():
    # Download GSE.
    gse1 = GEOparse.get_GEO(geo='GSE44621', destdir="c:/data/Geo/Microarray/GSE") #B.cells, monocytes.
    gse2 = GEOparse.get_GEO(geo='GSE102693', destdir="c:/data/Geo/Microarray/GSE") #CD4.T.cells.
    gse3 = GEOparse.get_GEO(geo='GSE68003', destdir="c:/data/Geo/Microarray/GSE") #CD8.T.cells

    # Parse platform data.
    gpl = gse1.gpls['GPL6244']
    gpl_table = gpl.table[['ID', 'gene_assignment']].copy()
    gpl_table['gene_assignment'] = gpl_table.gene_assignment.str.split('\/\/', expand=True)[1] #.reset_index()
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    data = pd.DataFrame()
    gsm1 = gse1.gsms['GSM1087965']
    gsm1 = gsm1.table.set_index('ID_REF')
    data = pd.concat([data, gsm1], axis=1)
    gsm2 = gse1.gsms['GSM1087971']
    gsm2 = gsm2.table.set_index('ID_REF')
    data = pd.concat([data, gsm2], axis=1)
    gsm3 = gse2.gsms['GSM2743091']
    gsm3 = gsm3.table.set_index('ID_REF')
    data = pd.concat([data, gsm3], axis=1)
    gsm4 = gse3.gsms['GSM1660727']
    gsm4 = gsm4.table.set_index('ID_REF')
    data = pd.concat([data, gsm4], axis=1)
    data.dropna(inplace=True)
    data.columns = ['B.cells', 'monocytic.lineage', 'CD4.T.cells', 'CD8.T.cells']
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)

    prop_data = np.empty((0, 4))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(4), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2.to_csv('c:/data/Geo/Microarray/mix-44621/mix.csv')
    prop_df.to_csv('c:/data/Geo/Microarray/mix-44621/prop.csv')
    data.to_csv(f'c:/data/Geo/Microarray/mix-44621/pure.csv')


def create_microarray_mix_107019():
    data = pd.read_csv('c:/data/Geo/Microarray/mix-107019/mix-orig.csv')
    gpl = pd.read_csv('c:/data/Geo/Microarray/mix-107019/GPL.csv')
    data = pd.merge(data, gpl, on='ID_REF', how='left')
    data.dropna(inplace=True)
    data.drop('ID_REF', axis=1, inplace=True)
    data['Symbol'] = data.Symbol.str.strip()
    data = data.groupby('Symbol').agg(np.mean)
    data.to_csv('c:/data/Geo/Microarray/mix-107019/mix.csv')


def find_correlated_genes_77343():
    mix = pd.read_csv('c:/data/Geo/Microarray/mix-77343/mix.csv', index_col=0).T
    prop = pd.read_csv('c:/data/Geo/Microarray/mix-77343/prop.csv', index_col=0).T
    # test = np.corrcoef(mix.T.loc['A1BG'], prop.T.loc['neutrophils'])[0][1]
    df1 = prop['neutrophils'].T
    df1.index = mix.index
    df2 = mix.corrwith(df1, axis=0)
    df3 = df2.sort_values(ascending=False)
    # df3.index[0:80] = ['TMEM55A', 'MEGF9', 'WDFY3', 'BEST1', 'STX3', 'GCA', 'TLR8', 'UBXN2B',
       # 'MXD1', 'GLT1D1', 'TMEM88', 'ITPRIP', 'ROPN1L', 'RNF149', 'DCP2',
       # 'TLR6', 'NFIL3', 'CEACAM4', 'HAL', 'MNDA', 'PFKFB4', 'RP2', 'TET2',
       # 'RNF24', 'LRRC4', 'CMTM2', 'REPS2', 'KIAA0232', 'LIN7A', 'AQP9',
       # 'LRRK2', 'CMTM6', 'NRBF2', 'LILRB3', 'NFAM1', 'TLR2', 'ALPK1', 'MOSC1',
       # 'DOCK5', 'CLEC4E', 'SLC22A4', 'ACSL1', 'KCNE3', 'QPCT', 'SULT1B1',
       # 'PLXNC1', 'DENND3', 'FRAT2', 'LRP10', 'TMEM49', 'NAMPT', 'RASGRP4',
       # 'TNFRSF10C', 'TRIM25', 'USP32', 'BST1', 'DKFZp761E198', 'CYB5R4',
       # 'MTMR3', 'CPD', 'RRAGD', 'BCL6', 'ZBTB34', 'CPEB4', 'FCGR2A', 'EGLN1',
       # 'IL1R2', 'RBM47', 'CSF3R', 'NPL', 'PPP1R3B', 'FPR1', 'TLR4', 'KCNJ15',
       # 'PPP4R1', 'CEACAM3', 'SPOPL', 'RTN3', 'LAMP2', 'HSPA6']
    df4 = prop['monocytic.lineage'].T
    df4.index = mix.index
    df5 = mix.corrwith(df4, axis = 0)
    df6 = df5.sort_values(ascending=False)
    # df6.index[0:40] = ['NAGA', 'CD68', 'ADAP2', 'ANXA2', 'SLC43A3', 'ZNF385A', 'RASSF4',
       # 'RIN2', 'CPVL', 'CD33', 'LOC284837', 'PLXNB2', 'CRTAP', 'CECR1', 'KLF4',
       # 'ANXA2P2', 'NPC2', 'CD86', 'TTYH2', 'CD300C', 'KCTD12', 'DPYSL2',
       # 'PEA15', 'CYBB', 'CST3', 'CTNND1', 'PLD3', 'C16orf70', 'GPBAR1',
       # 'MS4A7', 'FLVCR2', 'MYCL1', 'ADAM15', 'CTTNBP2NL', 'KCNMB1', 'CCR2',
       # 'SLC46A2', 'MAN2B1', 'RPS6KA4', 'SLC37A2']


def find_correlated_genes_40240():
    mix = pd.read_csv('c:/data/Geo/Microarray/mix-40240/mix.csv', index_col=0).T
    prop = pd.read_csv('c:/data/Geo/Microarray/mix-40240/prop.csv', index_col=0).T
    # test = np.corrcoef(mix.T.loc['A1BG'], prop.T.loc['neutrophils'])[0][1]
    df1 = prop['neutrophils'].T
    df1.index = mix.index
    df2 = mix.corrwith(df1, axis=0)
    df3 = df2.sort_values(ascending=False)
    # df3.index[0:80] = ['GAB2', 'BEST1', 'ADAMTSL4-AS1', 'LOC729603', 'RBM47', 'ZDHHC18',
       # 'DOCK5', 'RASGRP4', 'PFKFB4', 'CHST15', 'STX3', 'PELI2', 'TNFRSF1A',
       # 'PHC2', 'ABTB1', 'PREX1', 'NDST1', 'RN7SL600P', 'NFAM1', 'CEACAM3',
       # 'NOTCH1', 'EPHB1', 'DYSF', 'NINJ1', 'C5AR1', 'KDM6B', 'C19orf35',
       # 'TMEM127', 'LILRB3', 'MYO1F', 'CXCR2', 'ARAP1', 'MBOAT7', 'RN7SL473P',
       # 'RNF24', 'SEC14L1', 'SIGLEC9', 'REPS2', 'SIRPD', 'ITPRIP', 'IL6R',
       # 'IMPDH1', 'CEACAM4', 'CSF3R', 'NCF4', 'PADI4', 'SLED1', 'LRP10', 'HAL',
       # 'NUMB', 'DENND3', 'ABHD5', 'DEF8', 'ATG2A', 'ITPK1', 'STK40', 'MXD1',
       # 'NHSL2', 'KLHL21', 'MEGF9', 'IGF2R', 'IL1R1', 'LAT2', 'PLAUR', 'ERV3-1',
       # 'NADK', 'WAS', 'OR52K2', 'SLC45A4', 'ADCY4', 'PRRG4', 'SLC6A6',
       # 'LINC00999', 'PPP1R3B', 'MGAM', 'FCGRT', 'CXCR1', 'GLT1D1', 'DHX34',
       # 'MAST3']
    df4 = prop['monocytic.lineage'].T
    df4.index = mix.index
    df5 = mix.corrwith(df4, axis = 0)
    df6 = df5.sort_values(ascending=False)
    # df6.index[0:40] = ['LGALS1', 'CD300C', 'TTYH2', 'GSTP1', 'PHPT1', 'CLCN5', 'CECR1', 'CD68',
       # 'RHOC', 'CST3', 'NLN', 'PLXNB2', 'KIAA0930', 'CD300E', 'CNP', 'TMEM205',
       # 'CALHM2', 'GPBAR1', 'DPYSL2', 'CD1D', 'GLB1L', 'SH3TC1', 'UBXN11',
       # 'AIMP2', 'ACP2', 'PLXND1', 'TBC1D8', 'CPVL', 'ZNF385A', 'ABI3', 'COMT',
       # 'TSPAN4', 'SLC43A3', 'FGD2', 'RNH1', 'TMEM150B', 'MS4A14', 'CSF1R',
       # 'ANAPC2', 'ATP6V1F']


def create_microarray_mix_gpl6244():
    # Download GSE.
    gse1 = GEOparse.get_GEO(geo='GSE23321', destdir="c:/data/Geo/Microarray/GSE") #naive.CD8.T.cells, memory.CD8.T.cells. Log2.
    gse2 = GEOparse.get_GEO(geo='GSE53455', destdir="c:/data/Geo/Microarray/GSE") #naive.CD4.T.cells, memory.CD4.T.cells. Linear.
    gse3 = GEOparse.get_GEO(geo='GSE102693', destdir="c:/data/Geo/Microarray/GSE") #regulatory.T.cells. Log2.
    gse4 = GEOparse.get_GEO(geo='GSE42724', destdir="c:data/Geo/Microarray/GSE") #naive.B.cells, memory.B.cells. Log2.

    # Parse platform data.
    gpl = gse1.gpls['GPL6244']
    gpl_table = gpl.table[['ID', 'gene_assignment']].copy()
    gpl_table['gene_assignment'] = gpl_table.gene_assignment.str.split('\/\/', expand=True)[1] #.reset_index()
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    data = pd.DataFrame()
    gsm1 = gse1.gsms['GSM572423'] #memory.CD8.T.cells
    gsm1 = gsm1.table.set_index('ID_REF')
    gsm1 = 2 ** gsm1
    data = pd.concat([data, gsm1], axis=1)
    gsm2 = gse1.gsms['GSM572425'] #naive.CD8.T.cells
    gsm2 = gsm2.table.set_index('ID_REF')
    gsm2 = 2 ** gsm2
    data = pd.concat([data, gsm2], axis=1)

    gsm3 = gse2.gsms['GSM1293938'] #naive.CD4.T.cells
    gsm3 = gsm3.table.set_index('ID_REF')
    data = pd.concat([data, gsm3], axis=1)
    gsm4 = gse2.gsms['GSM1293949'] #memory.CD4.T.cells
    gsm4 = gsm4.table.set_index('ID_REF')
    data = pd.concat([data, gsm4], axis=1)

    gsm5 = gse3.gsms['GSM2743108'] #regulatory.T.cells
    gsm5 = gsm5.table.set_index('ID_REF')
    gsm5 = 2 ** gsm5
    data = pd.concat([data, gsm5], axis=1)

    gsm6 = gse4.gsms['GSM1048787'] #naive.B.cells
    gsm6 = gsm6.table.set_index('ID_REF')
    gsm6 = 2 ** gsm6
    data = pd.concat([data, gsm6], axis=1)
    gsm7 = gse4.gsms['GSM1048789'] #memory.B.cells
    gsm7 = gsm7.table.set_index('ID_REF')
    gsm7 = 2 ** gsm7
    data = pd.concat([data, gsm7], axis=1)

    data.dropna(inplace=True)
    data.columns = ['memory.CD8.T.cells', 'naive.CD8.T.cells', 'naive.CD4.T.cells', 'memory.CD4.T.cells', 'regulatory.T.cells',
                    'naive.B.cells', 'memory.B.cells']
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    # Quantile normalize data.
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    data = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    prop_data = np.empty((0, 7))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(7), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2.to_csv('c:/data/Geo/Microarray/mix-6244/mix.csv')
    prop_df.to_csv('c:/data/Geo/Microarray/mix-6244/prop.csv')
    data.to_csv(f'c:/data/Geo/Microarray/mix-6244/pure.csv')


def create_microarray_mix_gpl6244_2():
    # Download GSE.
    gse1 = GEOparse.get_GEO(geo='GSE23321', destdir="c:/data/Geo/Microarray/GSE") #naive.CD8.T.cells, memory.CD8.T.cells. Log2.
    gse2 = GEOparse.get_GEO(geo='GSE32901', destdir="c:/data/Geo/Microarray/GSE") #naive.CD4.T.cells, memory.CD4.T.cells. Log2.
    gse3 = GEOparse.get_GEO(geo='GSE18893', destdir="c:/data/Geo/Microarray/GSE") #regulatory.T.cells. Log2.
    gse4 = GEOparse.get_GEO(geo='GSE51528', destdir="c:data/Geo/Microarray/GSE") #naive.B.cells, memory.B.cells. Log2.

    # Parse platform data.
    gpl = gse1.gpls['GPL6244']
    gpl_table = gpl.table[['ID', 'gene_assignment']].copy()
    gpl_table['gene_assignment'] = gpl_table.gene_assignment.str.split('\/\/', expand=True)[1] #.reset_index()
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    data = pd.DataFrame()
    gsm1 = gse1.gsms['GSM572427'] #memory.CD8.T.cells
    gsm1 = gsm1.table.set_index('ID_REF')
    gsm1 = 2 ** gsm1
    data = pd.concat([data, gsm1], axis=1)
    gsm2 = gse1.gsms['GSM572429'] #naive.CD8.T.cells
    gsm2 = gsm2.table.set_index('ID_REF')
    gsm2 = 2 * gsm2
    data = pd.concat([data, gsm2], axis=1)

    gsm3 = gse2.gsms['GSM814478'] #naive.CD4.T.cells
    gsm3 = gsm3.table.set_index('ID_REF')
    gsm3 = 2 ** gsm3
    data = pd.concat([data, gsm3], axis=1)
    gsm4 = gse2.gsms['GSM814482'] #memory.CD4.T.cells
    gsm4 = gsm4.table.set_index('ID_REF')
    gsm4 = 2 ** gsm4
    data = pd.concat([data, gsm4], axis=1)

    gsm5 = gse3.gsms['GSM468272'].table[['ID_REF', 'VALUE']] #regulatory.T.cells
    gsm5 = gsm5.set_index('ID_REF')
    gsm5 = 2 ** gsm5
    data = pd.concat([data, gsm5], axis=1)

    gsm6 = gse4.gsms['GSM1247420'] #naive.B.cells
    gsm6 = gsm6.table.set_index('ID_REF')
    gsm6 = 2 ** gsm6
    data = pd.concat([data, gsm6], axis=1)
    gsm7 = gse4.gsms['GSM1247426'] #memory.B.cells
    gsm7 = gsm7.table.set_index('ID_REF')
    gsm7 = 2 ** gsm7
    data = pd.concat([data, gsm7], axis=1)

    data.dropna(inplace=True)
    data.columns = ['memory.CD8.T.cells', 'naive.CD8.T.cells', 'naive.CD4.T.cells', 'memory.CD4.T.cells', 'regulatory.T.cells',
                    'naive.B.cells', 'memory.B.cells']
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    # Quantile normalize data.
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    data = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    prop_data = np.empty((0, 7))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(7), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2.to_csv('c:/data/Geo/Microarray/mix-6244v2/mix.csv')
    prop_df.to_csv('c:/data/Geo/Microarray/mix-6244v2/prop.csv')
    data.to_csv(f'c:/data/Geo/Microarray/mix-6244v2/pure.csv')


def create_microarray_mix_gpl6244_3():
    # Download GSE.
    gse1 = GEOparse.get_GEO(geo='GSE23321', destdir="c:/Temp/GSE") #naive.CD8.T.cells, memory.CD8.T.cells. Log2.
    gse2 = GEOparse.get_GEO(geo='GSE52129', destdir="c:/Temp/GSE") #naive.CD4.T.cells. Log2.
    gse6 = GEOparse.get_GEO(geo='GSE73968', destdir="c:/Temp/GSE") #memory.CD4.T.cells. Log2.
    gse3 = GEOparse.get_GEO(geo='GSE64176', destdir="c:/Temp/GSE") #regulatory.T.cells. Linear.
    gse4 = GEOparse.get_GEO(geo='GSE104373', destdir="c:/Temp/GSE") #monocytes. Log2.

    # Parse platform data.
    gpl = gse1.gpls['GPL6244']
    gpl_table = gpl.table[['ID', 'gene_assignment']].copy()
    gpl_table['gene_assignment'] = gpl_table.gene_assignment.str.split('\/\/', expand=True)[1] #.reset_index()
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    data = pd.DataFrame()
    gsm1 = gse1.gsms['GSM572431'] #memory.CD8.T.cells
    gsm1 = gsm1.table.set_index('ID_REF')
    gsm1 = 2 ** gsm1
    data = pd.concat([data, gsm1], axis=1)
    gsm2 = gse1.gsms['GSM572433'] #naive.CD8.T.cells
    gsm2 = gsm2.table.set_index('ID_REF')
    gsm2 = 2 * gsm2
    data = pd.concat([data, gsm2], axis=1)

    gsm3 = gse2.gsms['GSM1260083'] #naive.CD4.T.cells
    gsm3 = gsm3.table.set_index('ID_REF')
    gsm3 = 2 ** gsm3
    data = pd.concat([data, gsm3], axis=1)

    gsm4 = gse6.gsms['GSM1906885'] #memory.CD4.T.cells
    gsm4 = gsm4.table.set_index('ID_REF')
    gsm4 = 2 ** gsm4
    data = pd.concat([data, gsm4], axis=1)

    gsm5 = gse3.gsms['GSM1565824'].table[['ID_REF', 'VALUE']] #regulatory.T.cells
    gsm5 = gsm5.set_index('ID_REF')
    # gsm5 = gsm5.apply(lambda x: np.log2(x))
    data = pd.concat([data, gsm5], axis=1)

    gsm6 = gse4.gsms['GSM2796198'] #monocytes
    gsm6 = gsm6.table.set_index('ID_REF')
    gsm6 = 2 ** gsm6
    data = pd.concat([data, gsm6], axis=1)

    data.dropna(inplace=True)
    data.columns = ['memory.CD8.T.cells', 'naive.CD8.T.cells', 'naive.CD4.T.cells', 'memory.CD4.T.cells', 'regulatory.T.cells',
                    'monocytes']
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    # Quantile normalize data.
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    data = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()


    prop_data = np.empty((0, 6))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,50):
        prop = np.random.dirichlet(np.ones(6), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,50)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,50)]]
    data2 = data2.apply(lambda x: np.log2(x))
    data2.to_csv('c:/data/Geo/Microarray/mix-6244v3/mix.csv')
    prop_df.to_csv('c:/data/Geo/Microarray/mix-6244v3/prop.csv')
    data = data.apply(lambda x: np.log2(x))
    data.to_csv(f'c:/data/Geo/Microarray/mix-6244v3/pure.csv')


def create_microarray_mix_gse25638():
    # Download GSE.
    gse = GEOparse.get_GEO(geo='GSE25638', destdir="c:/Temp/GSE")

    # Parse platform data.
    gpl = gse.gpls['GPL570']
    gpl_table = gpl.table[['ID', 'Gene Symbol']].copy()
    gpl_table['gene_assignment'] = gpl_table['Gene Symbol'].str.split('\/\/\/', expand=True)[0] #.reset_index()
    gpl_table.drop('Gene Symbol', axis=1, inplace=True)
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    data = pd.DataFrame()
    gsm1 = gse.gsms['GSM629990'] #T.cells
    gsm1 = gsm1.table.set_index('ID_REF')
    gsm1 = 2 ** gsm1
    data = pd.concat([data, gsm1], axis=1)

    gsm2 = gse.gsms['GSM629986'] #memory.B.cells
    gsm2 = gsm2.table.set_index('ID_REF')
    gsm2 = 2 * gsm2
    data = pd.concat([data, gsm2], axis=1)

    gsm3 = gse.gsms['GSM629988'] #naive.B.cells
    gsm3 = gsm3.table.set_index('ID_REF')
    gsm3 = 2 ** gsm3
    data = pd.concat([data, gsm3], axis=1)

    data.dropna(inplace=True)
    data.columns = ['T.cells', 'memory.B.cells', 'naive.B.cells']
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    # Quantile normalize data.
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    data = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    num_cells = 3
    num_mixes = 20
    prop_data = np.empty((0, num_cells))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,num_mixes):
        prop = np.random.dirichlet(np.ones(num_cells), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,num_mixes)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,num_mixes)]]
    data2.to_csv('c:/work/data/Geo/Microarray/mix-25638/mix.csv')
    prop_df.to_csv('c:/work/data/Geo/Microarray/mix-25638/prop.csv')
    data.to_csv(f'c:/work/data/Geo/Microarray/mix-25638/pure.csv')


def create_microarray_mix_gse():
    gse_name = '93683' #'15659'
    # Download GSE.
    gse = GEOparse.get_GEO(geo='GSE'+gse_name, destdir="c:/Temp/GSE")

    # Parse platform data.
    gpl = gse.gpls['GPL570']
    gpl_table = gpl.table[['ID', 'Gene Symbol']].copy()
    gpl_table['gene_assignment'] = gpl_table['Gene Symbol'].str.split('\/\/\/', expand=True)[0] #.reset_index()
    gpl_table.drop('Gene Symbol', axis=1, inplace=True)
    gpl_table.replace(["NaN", 'NaT', 'nan', 'None'], np.nan, inplace=True)
    gpl_table.dropna(inplace=True)
    gpl_table.columns = ['ID_REF', 'gene']

    # Parse GSM examples.
    # gse_samples = ['392071', '392072', '392073']
    gse_samples = ['2460456', '2460433']
    # gse_samples_names = ['naive.T.cells', 'memory.CD4.T.cells', 'regulatory.T.cells']
    gse_samples_names = ['memory.CD8.T.cells', 'naive.CD8.T.cells']
    data = pd.DataFrame()
    gsm1 = gse.gsms['GSM'+gse_samples[0]]
    gsm1 = gsm1.table.set_index('ID_REF')
    # gsm1.drop(['ABS_CALL', 'DETECTION P-VALUE'], axis=1, inplace=True)
    data = pd.concat([data, gsm1], axis=1)

    gsm2 = gse.gsms['GSM'+gse_samples[1]]
    gsm2 = gsm2.table.set_index('ID_REF')
    # gsm2.drop(['ABS_CALL', 'DETECTION P-VALUE'], axis=1, inplace=True)
    data = pd.concat([data, gsm2], axis=1)

    gsm3 = gse.gsms['GSM'+gse_samples[2]]
    gsm3 = gsm3.table.set_index('ID_REF')
    gsm3.drop(['ABS_CALL', 'DETECTION P-VALUE'], axis=1, inplace=True)
    data = pd.concat([data, gsm3], axis=1)

    data = 2 ** data
    data.dropna(inplace=True)
    data.columns = gse_samples_names
    data = data.reset_index()
    data = pd.merge(data, gpl_table, on='ID_REF', how='left')
    data.drop('ID_REF', axis=1, inplace=True)
    data.dropna(inplace=True)
    data['gene'] = data.gene.str.strip()
    data = data.groupby('gene').agg(np.mean)
    # Quantile normalize data.
    rank_mean = data.stack().groupby(data.rank(method='first').stack().astype(int)).mean()
    data = data.rank(method='min').stack().astype(int).map(rank_mean).unstack()

    num_cells = 2
    num_mixes = 20
    prop_data = np.empty((0, num_cells))
    data2 = data.iloc[:,0:1].copy()
    # Loop on mixes.
    for i in range(0,num_mixes):
        prop = np.random.dirichlet(np.ones(num_cells), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,num_mixes)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,num_mixes)]]
    data2.to_csv(f'c:/work/data/Geo/Microarray/mix-{gse_name}/mix.csv')
    prop_df.to_csv(f'c:/work/data/Geo/Microarray/mix-{gse_name}/prop.csv')
    data.to_csv(f'c:/work/data/Geo/Microarray/mix-{gse_name}/pure.csv')


def create_mix_rnaseq_89134():
    base = 'c:/work/data/Geo/RNASeq/mix-89134/'
    data = pd.read_csv(base+'GSE89134_Log2NormalizedCPM.csv', index_col=0)
    cell_type = ['memory.CD8.T.cells', 'naive.CD8.T.cells']
    # cols = []
    # for i in cell_type:
    #    annot_temp = annotations.loc[annotations['cell.types'] == i]
    #    row = annot_temp.iloc[0] #randrange(len(annot_temp))
    #    cols += [row['cells']]
    # data = pd.read_csv('c:/data/GEO/RNASeq/mix-115978/GSE115978_tpm.csv', index_col=0)[cols]
    data = 2**data
    data.columns = cell_type
    create_props(base, data)


def create_props(base, data):
    data2 = data.copy()
    prop_data = np.empty((0, len(data.columns)))
    # Loop on mixes.
    num_mixes = 50
    for i in range(0,num_mixes):
        prop = np.random.dirichlet(np.ones(len(data.columns)), size=1)
        j = 0
        data2[f'mix{i}'] = 0
        # Loop on cell types.
        for col in data:
            noise = np.random.normal(1, 0.25, len(data[col]))
            data2[f'mix{i}'] += data[col] * prop[0][j] * noise
            j+=1
        prop_data = np.append(prop_data, prop, axis=0)

    prop_df = pd.DataFrame(data=prop_data.T, columns=[f'mix{i}' for i in range(0,num_mixes)], index=data.columns)
    data2 = data2[[f'mix{i}' for i in range(0,num_mixes)]]
    data2.to_csv(base+'mix10.csv')
    prop_df.index = data.columns
    prop_df.to_csv(base+'prop10.csv')
    data.to_csv(base+'pure10.csv')


def create_mix_rnaseq_97861():
	base = 'c:/Users/Admin/Google Drive/src/Geo/RNASeq/mix-97861-TPM/'
	data = pd.read_csv(base+'GSM3321078_Tregs.csv', index_col=0)['Tregs']
	gsm1 = pd.read_csv(base+'GSM3005222_Memory_B_cells.csv', index_col=0)['memory.B.cells']
	gsm2 = pd.read_csv(base+'GSM3005216_Naive_B_cells.csv', index_col=0)['naive.B.cells']
	gsm3 = pd.read_csv(base+'GSE97861_Memory_CD4_T_cells_Naive_CD4_T_cells.csv', index_col=0)[['GS0272_NAIVE_RNASeq', 'GS0272_TCM_RNASeq']]
	gsm3.columns = ['naive.CD4.T.cells', 'memory.CD4.T.cells']
	gsm6 = pd.read_csv(base+'GSE114407_B_cells_CD4_T_cells_CD8_T_cells_Monocytes.csv', index_col=0)[['E_80_Postsort_CD14_M', 'E_80_Postsort_CD8_T']]
	gsm6.columns = ['monocytes', 'naive.CD8.T.cells']
	gsm4 = pd.read_csv(base+'GSM1888829_4659-Memory-T24.csv', index_col=0)['memory.CD4.T.cells']
	gsm4.dropna(inplace=True)
	gsm5 = pd.read_csv(base+'GSE111907_Fibroblasts_Endothelial.csv', index_col=0)[['fibroblasts', 'endothelial.cells']]
	gsm5.columns = ['fibroblasts', 'endothelial.cells']
	# gsm5 = pd.read_csv('c:/Users/Admin/Google Drive/src/Geo/RNASeq/mix-121922/GSE121922_All_Processed_Data.csv', index_col=0, encoding='latin-1')['endothelial.cells']
	# gsm7 = pd.read_csv(base+'log2FPKM_CD8TCell_CRC.csv', index_col=0)['memory.CD8.T.cells']
	genes = pd.read_csv(base+'genes2.csv')
	# gsm7 = pd.concat([gsm7, genes], axis=1).set_index('gene')
	# gsm7.dropna(inplace=True)
	# gsm7 = gsm7.groupby(gsm7.index).mean()
	# gsm7 = pd.read_csv(base+'GSE90728_All_Counts.lung_cancer.csv', index_col=0, encoding='latin-1')['memory.CD8.T.cells']
	# gsm7.drop_duplicates(inplace=True)
	# gsm7 = pd.read_csv(base+'GSE63144_GEO_Schultze_RNA_seq_normalizedData.csv', index_col=0)['memory.CD8.T.cells']
	gsm7 = pd.read_csv(base+'GSE85527_AbATE_single_cell_profile_raw_counts.csv', index_col=0)['memory.CD8.T.cells']
	gsm7.dropna(inplace=True)
	genes['symbol'] = genes.gene_name.apply(lambda x: str(x).split('.')[0])
	genes = genes.set_index('symbol')
	genes.drop('gene_name', inplace=True, axis=1)
	genes.drop_duplicates(inplace=True)
	gsm7 = pd.concat([gsm7, genes], axis=1).set_index('gene')
	gsm7 = gsm7.groupby(gsm7.index).mean()
	#Mix8.
	#gsm8 = pd.read_csv(base+'GSE100382_RNAseq_gene_expression-Macrophages.txt', index_col=0, delimiter='\t')['macrophages']
	#Mix9.
	gsm8 = pd.read_csv(base+'GSE104174_TPM_SSc_study-Macrophages.txt', index_col=0, delimiter='\t')['macrophages']
	gsm8.dropna(inplace=True)
	gsm8 = pd.concat([gsm8, genes], axis=1).set_index('gene')
	gsm8 = gsm8.groupby(gsm8.index).mean()
	#Mix10.
	gsm9 = pd.read_csv(base+'GSE84865_genes_expression_tpm.txt', index_col=0, delimiter='\t')['myeloid.dendritic.cells']
	data = pd.concat([data, gsm1, gsm2, gsm3['naive.CD4.T.cells'], gsm4, gsm5, gsm6, gsm7, gsm8, gsm9], axis=1)
	data.fillna(0, inplace=True)
	create_props(base, data)
 

# Create Microarray mixes using the CSV file from the challenge. Didn't work due to platforms differences.
# cells_lm8 = ['B.cells', 'CD4.T.cells', 'CD8.T.cells', 'NK.cells', 'neutrophils', 'monocytic.lineage', 'fibroblasts', 'endothelial.cells']
# studies = pd.read_csv(f'c:/Work/src/Yada/yada/data/Challenge/geo-rnaseq-immune-cells.csv')
# platforms = studies.platform_id.unique()
# platform = platforms[randrange(len(platforms))]
# studies2 = studies.loc[studies.platform_id == platform].copy()
# cells = list(set(studies2['cell.type'].unique()) & set(cells_lm8))

# gene_list = list(pure.index)
# sample_size = int(0.35*len(pure))
# sorted_sample = [
#     gene_list[i] for i in sorted(random.sample(range(len(gene_list)), sample_size))
# ]
# mix = mix.loc[gene_list]


if __name__ == '__main__':
    create_mix_rnaseq_97861()
    #create_lm_107011()

    """
    # Create test cases.
    fractions = pd.read_csv('c:/input/quanTIseq_SimRNAseq_read_fractions.csv', index_col=0)
    # fractions = fractions.loc[fractions.Tumor == 0]
    mix = pd.read_csv('c:/input/quanTIseq_SimRNAseq_mixture_orig.csv', index_col=0)
    for i in range(100):
        temp_frac = fractions.sample(frac=0.01)
        mix_temp = mix[temp_frac.index]
        temp_frac.to_csv(f'c:/input/frac_{i}.csv')
        mix_temp.to_csv(f'c:/input/mix_{i}.csv')
    """