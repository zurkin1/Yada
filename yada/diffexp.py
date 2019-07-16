import pandas as pd

#This function does a gene differentiation analysis by finding majoroizing genes, with the factor.
def gene_diff(pure):
    print('Running gene differentiation')
    gene_list_per_cell = {}
    # Loop on all cells
    for cell_type in pure:
        gene_list = []
        data = pure.copy()
        for factor in range(100, 2, -1):
            data2 = data.copy()
            data2['max_val'] = data.drop(cell_type, axis=1).max(axis=1)
            data2['min_val'] = data.drop(cell_type, axis=1).min(axis=1)
            data2 = data2.loc[(data2[cell_type] > 0.5 * factor * data2.max_val)] #Ignore cases where all other mixes are zeros.
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
    return gene_list_df, factor_list