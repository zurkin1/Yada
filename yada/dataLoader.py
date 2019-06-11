import pandas as pd

def load_data(dataset_name):
    # In some cases datasets are missing (cell proportions do not sum to 1)
    # This parameter allows us to calibrate this information into the model.
    mixtures_sum = []

    if(dataset_name == 'DSA'):
        #DSA (3,11) ,liver, brain, lung
        mix = pd.read_csv('data/DSA/mix.csv', index_col=0)
        pure = pd.read_csv('data/DSA/pure.csv', index_col=0)
        real_weight = pd.read_csv('data/DSA/weight.csv')
        real_weight = real_weight[3:]/100
        real_weight = real_weight.reset_index(drop=True)
        mix = 2 ** mix
        pure = 2 ** pure
        other_result = pd.read_csv('data/DSA/cibersort_result.csv', index_col=0)
        other_result.columns = real_weight.columns
        other_result.reset_index(drop=True, inplace=True)
    elif(dataset_name == 'Abbas'):
        #Abbas et. al (4,12), jurkat T-cels, im9, raji B-cell, thp1 Monocytic cells
        mix = pd.read_csv('data/Abbas/GSE11103_matrix_mixtures.csv', index_col=0)
        real_weight = pd.read_csv('data/Abbas/labels.csv', index_col=0)
        other_result = pd.read_csv('data/Abbas/CIBERSORT.Output_Job1.csv', index_col=0)
        pure = pd.read_csv('data/Abbas/GSE11103_matrix_classes.GSE11103_matrix_pure.bm.K999.0-signature gense.csv', index_col=0)
    elif(dataset_name == 'CIBERSORT'):
        #CIBERSORT Tutorial (22,5)
        mix = pd.read_csv('data/CIBERSORT/TutorialExampleMixtures-GEPs.csv', index_col=0)
        pure = pd.read_csv('data/CIBERSORT\LM22.csv', index_col=0)
        real_weight = pd.read_csv('data/CIBERSORT/TutorialExampleMixtures-GroundTruth.csv')
        other_result = pd.read_csv('data/CIBERSORT/cibersort_tutorial_results.csv')
    elif(dataset_name == '10x'):
        #Single Cell from 10x (5,10)
        mix = pd.read_csv('data/10x/mix.csv', index_col=0)
        real_weight = pd.read_csv('data/10x/weights.csv').T/100
        pure = pd.read_csv('data/10x/pure.csv', index_col=0)
        other_result = pd.read_csv('data/10x/cibersort.csv', index_col=0)
        other_result.index = real_weight.index
        other_result.columns = real_weight.columns
    elif(dataset_name == 'DeconRNASeq'):
        #DeconRNASeq (5,10)
        mix = pd.read_csv('data/DeconRNASeq/mix.csv', index_col = 0)
        pure = pd.read_csv('data/DeconRNASeq/pure.csv', index_col=0)
        real_weight = pd.read_csv('data/DeconRNASeq/real_weights.csv')
        other_result = pd.read_csv('data/DeconRNASeq/cibersort_result.csv')
        other_result.columns = real_weight.columns
    elif(dataset_name == 'EPIC'):
        #EPIC (4,4)
        mix = pd.read_csv('data/EPIC/mix.csv', index_col=0)
        pure = pd.read_csv('data/EPIC/TRefProfiles.csv', index_col=0)
        real_weight = pd.read_csv('data/EPIC/ground_truth.csv')
        other_result = pd.read_csv('data/EPIC/cibersort_results.csv')
        other_result.columns = real_weight.columns
        mixtures_sum = [0.19, 0.13, 0.7, 0.62]
    elif(dataset_name == 'TIMER'):
        #TIMER (6,5)
        mix = pd.read_csv('data/TIMER/mix.csv', index_col=0)
        real_weight = pd.read_csv('data/TIMER/real_values.csv', index_col=0)
        pure = pd.read_csv('data/TIMER/pure.csv', index_col = 0)[['B_cell','CD4_Tcell','CD8_Tcell','Neutrophil','Macrophage','Dendritic']]
        pure = pure.groupby(pure.index).mean()
        other_result = pd.read_csv('data/TIMER/cibersort_result.csv', index_col=0)
        mixtures_sum = [0.8, 1.3, 1.18, 1.7, 0.85]
    elif(dataset_name == 'BreastBlood'):
        #Breast blood (2,9)
        mix = pd.read_csv('data/BreastBlood_GSE29832/mix.txt', sep='\t', index_col=0)
        real_weight = pd.read_csv('data/BreastBlood_GSE29832/coef.txt', sep='\t', index_col=0).T
        pure = pd.read_csv('data/BreastBlood_GSE29832/pure2.txt', sep='\t', index_col=0)
        other_result = pd.read_csv('data/BreastBlood_GSE29832/cibersort_result.csv')
        other_result.columns = real_weight.columns
        other_result.index = real_weight.index
    elif(dataset_name == 'RatBrain'):
        #RatBrain (4,10)
        mix = pd.read_csv('data/RatBrain_GSE19380/mix1.txt', sep='\t', index_col=0)
        real_weight = pd.read_csv('data/RatBrain_GSE19380/coef.txt', sep='\t', index_col=0).T
        pure = pd.read_csv('data/RatBrain_GSE19380/pure1.txt', sep='\t', index_col=0)
        other_result = pd.read_csv('data/RatBrain_GSE19380/cibersort_result.csv')
        mix = 2 ** mix
        pure = 2 ** pure
        other_result.columns = real_weight.columns
        other_result.index = real_weight.index
    elif(dataset_name == 'PertU'):
        #Pert uncultured (11,4)
        mix = pd.read_csv('data/PERT_uncultured_GSE40830/mix.txt', sep='\t', index_col=0)
        real_weight = pd.read_csv('data/PERT_uncultured_GSE40830/coef.txt', sep='\t', index_col=0).T/100
        pure = pd.read_csv('data/PERT_uncultured_GSE40830/sig.txt', sep='\t', index_col=0)
        other_result = pd.read_csv('data/PERT_uncultured_GSE40830/ciber_result.csv', index_col=0)
    elif (dataset_name == 'Retina'):
        # (2,8)
        mix = pd.read_csv('data/Retina_GSE33076/mix.txt', sep='\t', index_col=0)
        real_weight = pd.read_csv('data/Retina_GSE33076/coef.txt', sep='\t', index_col=0).T
        pure = pd.read_csv('data/Retina_GSE33076/pure.txt', sep='\t', index_col=0)
        other_result = pd.read_csv('data/Retina_GSE33076/cibersort_result.csv', index_col=0)
        mix = 2 ** mix
        pure = 2 ** pure
        other_result.columns = real_weight.columns
        other_result.index = real_weight.index
    elif(dataset_name == 'PertC'):
        #Pert uncultured (11,4)
        mix = pd.read_csv('data/PERT_Cultured_GSE16589/mix.txt', sep='\t', index_col=0)
        real_weight = pd.read_csv('data/PERT_Cultured_GSE16589/coef.txt', sep='\t', index_col=0).T/100
        pure = pd.read_csv('data/PERT_Cultured_GSE16589/sig.txt', sep='\t', index_col=0)
        other_result = pd.read_csv('data/PERT_Cultured_GSE16589/ciber_result.csv', index_col=0)
    elif(dataset_name == 'MySort'):
        #Pert uncultured (11,4)
        mix = pd.read_csv('data/MySort/mix.csv', index_col=0)
        real_weight = pd.read_csv('data/MySort/real_values.csv', index_col=0)/100
        pure = pd.read_csv('data/MySort/pure.csv', index_col=0)[real_weight.columns]
        other_result = pd.read_csv('data/MySort/cibersort_result.csv', index_col=0)

    both_genes = list(set(mix.index) & set(pure.index))
    pure = pure.reindex(both_genes)  # Drop genes that don't appear in mix.
    mix = mix.reindex(both_genes)
    return mix, pure, real_weight, other_result, mixtures_sum