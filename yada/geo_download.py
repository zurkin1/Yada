import pandas as pd
import GEOparse

studies = pd.read_csv('c:/data/Geo/geo-expression-array-immune-cells-studies.csv')
studies = studies.loc[studies['cell.types'] == 'fibroblast']
studies = studies.loc[studies.title.str.contains('HG-')]
global_table = pd.read_csv('c:/work/src/Yada/yada/global_table.csv', index_col=0)
counter = 1148
temp_studies = studies.copy()

for index, row in studies.iterrows():
#row = studies.iloc[0]
    gse = GEOparse.get_GEO(geo=row[0], destdir="c:/data/Geo/Microarray/fibroblast")
    #Parse platform data.
    for gpl_name, gpl in gse.gpls.items():
        #print("Metadata:",)
        #for key, value in gpl.metadata.items():
        #    print(" - %s : %s" % (key, ", ".join(value)))
        gpl_table = gpl.table[['ID', 'Gene Symbol']]
        gpl_table.columns = ['ID_REF', 'Symbol']
        break

    #Parse GSM examples.
    for gsm_name, gsm in gse.gsms.items():
        #print("Name: ", gsm_name)
        #print("Metadata:",)
        #for key, value in gsm.metadata.items():
        #    print(" - %s : %s" % (key, ", ".join(value)))
        table = gsm.table
        table = pd.merge(table, gpl_table, on='ID_REF', how='left')
        table.drop('ID_REF', axis=1, inplace=True)
        table.dropna(inplace=True)
        table = table.groupby('Symbol').VALUE.mean()
        global_table = pd.concat([global_table, table], axis=1).sum(axis=1)
        counter +=1
    print(f'{row[0]}, counter: {counter}.')
    global_table.to_csv('global_table.csv')
    temp_studies.drop(index, inplace=True)
    temp_studies.to_csv('geo-expression-array-immune-cells-studies-temp.csv', index=False)

global_table /= counter
global_table.to_csv('global_table.csv')