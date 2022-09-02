import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
from scipy.optimize import nnls, minimize
from scipy import stats
from scipy.spatial import minkowski_distance
from tslearn import metrics
from sklearn.svm import LinearSVR  # NuSVR
from sklearn.linear_model import Lasso
from sklearn.decomposition import FastICA
# import pymc3 as pm
import logging
import warnings
import random
from random import choice
#mport similaritymeasures
#from similaritymeasures import pcm
#import gseapy as gp
#from singscore import *
#from tqdm import tqdm
import multiprocessing as mp


warnings.filterwarnings("ignore")
logger = logging.getLogger("pymc3")
logger.propagate = False
logger.setLevel(logging.ERROR)


# This function calculate a deconvolution algorithm using basic ratios algorithm.
def basic_deconv(P, Q):
    dist = 0
    for index, value in P.items():
        # Subtract the average of other cells on this gene to compensate for over estimation of ratio.
        dist += ((value + 0.000001) / (Q.loc[index] + 0.000001))  # * \
        #((factor_list[cell_type][j] + 0.000001) / (factor_list[cell_type].sum() + 0.000001))
    return(dist)


def dtw_metric(P, Q):
    #fns = [metrics.dtw] #, similaritymeasures.dtw]
    #P = P.sample(frac = 1)
    #Q = Q[P.index]
    P1 = np.array([P.values, np.linspace(0, np.max(P), len(P))]) # np.arange(0, len(P))
    Q1 = np.array([Q.values, np.linspace(0, np.max(P), len(Q))])
    if (len(P) <= 5):
        return(P.mean() - Q.mean())
    #factor = max(0.01, np.corrcoef(P, Q)[0][1])
    #return similaritymeasures.dtw(Q1.T,P1.T)[0] #0.342
    return(metrics.dtw(P1, Q1))  # 0.342
    #return(pcm(P1, Q1))
    # return choice(fns)(P1, Q1)


# This function calculate a deconvolution algorithm using DTW ditance and DSA algorithm.
def dtw_deconv(mix, pure, gene_list_df, metric='dtw'):
    #mix[mix < 0] = 0
    gene_list_df.replace(["NaN", 'NaN', 'nan'], np.nan, inplace=True)
    num_cells = len(pure.columns)
    num_mixes = len(mix.columns)
    #mixtures_sum = [round(1 - abs(random.gauss(0, 0.045)), 2) for i in range(num_mixes)]
    mixtures_sum = [1 for i in range(num_mixes)]
    O_array = np.zeros((num_cells, num_mixes))  # The per cell sorted array - how far each mix is from the maximal mix.
    i = 0
    # Loop on all cell types.
    for cell_type in pure:
        cell_vals = []
        #gene_list_df[cell_type] = gene_list_df[cell_type].map(str.lower)
        cell_genelist = gene_list_df[cell_type].dropna().sample(frac=0.35) # round(random.gauss(0.4, 0.03), 2))  # 0.35
        #If marker list or sample list is short, don't sample.
        if (len(gene_list_df[cell_type].dropna()) < 8):  # or (len(cell_genelist) < 5)
            cell_genelist = gene_list_df[cell_type].dropna()
        #Make sure mix has all the genes in the list.
        cell_genelist = list(set(mix.index) & set(cell_genelist))
        mix_temp = mix.loc[cell_genelist]
        max_ind = mix_temp.mean().idxmax()  # Mix with maximum mean of gene expression.
        pure_temp = pure.loc[cell_genelist]
        pure_column = pure_temp[cell_type].copy()
        #if mix_temp.max().max() > 2*pure_temp.max().max():
        #    max_column = mix_temp[max_ind].copy()
            #print('.', end="")
        #else:
        max_column = pure_column
        max_column.sort_values(ascending=False, inplace=True)

        # Loop on all mixes.
        k = 0
        for mix_col in mix_temp:
            column = mix_temp[mix_col].copy()
            column = column[max_column.index]  # Sort according to the maximum column index.
            if metric == 'dtw':
                dist = dtw_metric(max_column, column)  # 0.342
            elif metric == 'avg':
                dist = np.mean(max_column - column)
            elif metric == 'abs':
                dist = np.mean(abs(max_column.values - column.values))
            elif metric == 'basic':
                dist = basic_deconv(max_column, column)
            elif metric == 'ks':
                dist, pval = stats.ks_2samp(max_column.values, column.values) # Kolmogarov Smirnoff distance.
            elif metric == 'euclid':
                dist = minkowski_distance(max_column, column, p=2)
            elif metric == 'taxi':
                dist = minkowski_distance(max_column, column, p=1)
            cell_vals.append(dist)
            k += 1
        O_array[i] = cell_vals
        i += 1

    #Scale to [0,1] in order for the sum of proportions to be meaningful. Scale along each cell.
    O_scaled = np.transpose(1 - minmax_scale(O_array, axis=1))
    solution = nnls(O_scaled, mixtures_sum)[0]
    solution_mat = np.diag(solution)
    estimate_wt = np.matmul(solution_mat, O_scaled.T)
    #estimate_wt = pd.DataFrame(data=estimate_wt, index=pure.columns).T
    #estimate_wt.index = mix.columns
    return estimate_wt


def run_dtw_deconv_ensemble(mix, pure, gene_list_df, metric = 'dtw'):
    #metrics = ['dtw', 'avg' , 'abs', 'ks', 'euclid', 'taxi'] #'basic',
    num_loops = 400
    pool = mp.Pool()
    num_mixes = len(mix.columns)
    num_cells = len(pure.columns)
    ens_estimate_wt = np.zeros((num_cells, num_mixes))
    
    results = [pool.apply_async(dtw_deconv, args=(mix, pure, gene_list_df, metric)) for i in range(num_loops)]
    #results = [dtw_deconv(mix, pure, gene_list_df, metric) for i in range(num_loops)]
    for ens_i in range(num_loops):
        print('\r', f"{ens_i / num_loops * 100:.0f}%", end='')
        ens_estimate_wt += results[ens_i].get()
    ens_estimate_wt /= num_loops
    ens_estimate_wt = pd.DataFrame(data=ens_estimate_wt, index=pure.columns).T
    ens_estimate_wt.index = mix.columns
    #ens_estimate_wt.to_csv('./data/results.csv')
    pool.close()
    return(ens_estimate_wt)


#SVR discovers a hyperplane that fits the maximal possible number of points within a constant distance ϵ, thus performing a regression.
def cibersort(mix, pure, gene_list_df):
    epsilon = 0.25
    # Standardize data.
    #mix = (mix - mix.mean()) / mix.std()
    #pure = (pure - np.mean(pure.values)) / np.std(pure.values)

    p = np.zeros([len(mix.columns), len(pure.columns)])
    i = 0
    for col in mix:
        mix_col = mix[col].copy()
        # model = NuSVR(nu=nu, kernel='linear')
        model = LinearSVR(epsilon=epsilon, fit_intercept=False, max_iter=2000)
        model.fit(pure, mix_col)
        # print(model.support_vectors_)
        w = model.coef_
        # w = w.clip(min=0) #Limit to positive values only in case we look for proportions (not in case we look only for correlation).
        # w = w/sum(sum(w))
        #w = w / sum(w)
        p[i] = w
        i += 1
    return pd.DataFrame(data = p, index = mix.columns, columns = pure.columns)


"""
# This function calculates a deconvolution algorithm using Bayes algorithm.
def dsection(mix, pure, params):
    # Use only differential meaningful genes.
    # gene_list = pd.concat([gene_list_df.iloc[:, i] for i in range(len(gene_list_df.columns))]).values
    # gene_list = list(set([x for x in gene_list if str(x) != 'nan']))
    # mix = mix.loc[gene_list]
    # pure = pure.loc[gene_list]

    # Standardize data.
    mix = (mix - mix.mean()) / mix.std()
    pure = (pure - np.mean(pure.values)) / np.std(pure.values)

    p_start = np.ones([len(mix.columns), len(pure.columns)])
    # Variance calculation
    joined_data = pd.concat([pure, mix], axis=1)
    joined_std = pd.DataFrame([joined_data.std(axis=1)] * len(pure.columns)).T

    basic_model = pm.Model()
    with basic_model:
        # Priors
        gamma = pm.Gamma('gamma', alpha=1, beta=1)
        epsilon = pm.Deterministic('ϵ', gamma + joined_std.T)
        pure_ = pm.Normal('pure_', mu=pure.T, shape=pure.T.shape, sd=epsilon)
        p_ = pm.Dirichlet('p_', a=p_start, shape=p_start.shape)

        # Likelihood (sampling distribution) of observations
        gamma2 = pm.Gamma('gamma2', alpha=1, beta=1)
        mu = pm.math.matrix_dot(p_, pure_)  # np.dot(pure_, p_) # + epsilon
        Y_obs = pm.Normal('Y_obs', mu=mu, observed=mix.T, sd=gamma2)

        start = pm.find_MAP(model=basic_model)
        step = pm.Metropolis()
        trace = pm.sample(1000, tune=50, step=step, cores=1, start=start, chains=4, progressbar=False)

    return(np.round(trace['p_'].mean(axis=0), 2).T)
"""

#Lasso regression adds a regularization term to residual sum of squares (called L1 norm), it may set coefficients to 0 and therefore perform feature selection.
def lasso(mix, pure, params):
    p = np.zeros([len(mix.columns), len(pure.columns)])
    i = 0
    for col in mix:
        mix_col = mix[col].copy()
        model = Lasso(alpha=0, positive=True, fit_intercept=False, selection='random')
        model.fit(pure, mix_col)
        w = model.coef_
        w = w / sum(w)
        p[i] = w
        i += 1
    return np.transpose(p)


# This function calculate a deconvolution algorithm using non-negative least square optimization of direct equation Prop * Pure = Mix
# under constraints of Prop == 1 (similar to deconRNASeq). Cellmix doesn't use pure matrix hence NNMF is needed.
# https://stackoverflow.com/questions/33385898/how-to-include-constraint-to-scipy-nnls-function-solution-so-that-it-sums-to-1?lq=1
def nnls_deconv_constrained(mix, pure, params = None):
    # if params is not None:
    #    mix = params[0].copy()
    #    pure = params[1].copy()
    num_cells = len(pure.columns)
    num_mixes = len(mix.columns)
    # gene_list = list(pure.index)
    # sample_size = int(0.35*len(pure))
    # sorted_sample = [gene_list[i] for i in sorted(random.sample(range(len(gene_list)), sample_size))]
    # mix = mix.loc[sorted_sample]
    # pure = pure.loc[sorted_sample]

    #pure = pure.sample(frac=0.75)
    #mix = mix.loc[pure.index]

    result = np.zeros((num_mixes, num_cells))

    # Loop on all mixes.
    i = 0
    for mix_col in mix:
        # Use nnls to get initial guess
        x0, rnorm = nnls(pure, mix[mix_col])
        # result[i] = x0
        # i += 1
        # return (result.T)

        # Define minimisation function
        def fn(x, A, b):
            return np.linalg.norm(A.dot(x) - b)

        # Define constraints and bounds
        cons = {'type': 'eq', 'fun': lambda x: np.sum(x) - 1}
        bounds = [[0., None]] * num_cells

        # Call minimisation subject to these values
        minout = minimize(fn, x0, args=(pure, mix[mix_col]), method='SLSQP', bounds=bounds, constraints=cons)
        x = minout.x

        result[i] = x
        i += 1
    return pd.DataFrame(data = result, index = mix.columns, columns = pure.columns)

#Independent Components Analysis is written as in Eq. (2.1) maximizing independence and non-Gaussianity of columns. It was first formulated by Herault and Jutten (1986).
#The independence can be measured with entropy, kurtosis, mutual information or negentropy measure.
def ica_deconv(mix, pure, gene_list_df):
    # Use only differential meaningful genes.
    gene_list = pd.concat([gene_list_df.iloc[:, i] for i in range(len(gene_list_df.columns))]).values
    gene_list = list(set([x for x in gene_list if str(x) != 'nan']))
    # mix = mix.loc[gene_list]
    pure = pure.loc[gene_list]
    pure = pure.sample(frac=0.75)
    mix = mix.loc[pure.index]

    # Standardize data.
    # mix = (mix - mix.mean()) / mix.std()
    # pure = (pure - np.mean(pure.values)) / np.std(pure.values)

    model = FastICA(n_components=len(pure.columns))  # , fun='cube'
    S_ = model.fit_transform(mix)  # Reconstruct signals.
    A_ = model.mixing_  # Get estimated mixing matrix.
    return A_.T


def pxcell(expr):
    expr = pd.read_csv(expr, index_col=0)
    # Reduce the expression dataset to contain only the required genes
    genes = pd.read_csv('./data/xCell/genes.csv', index_col=0)
    genes['x'] = genes.x.map(str.lower)
    expr = expr.loc[expr.index.intersection(list(genes.x))]
    #if (len(expr.index) < 5000):
    #   raise("ERROR: not enough genes")

    # Transform the expression to rank
    for col in expr:
        expr[col] = expr[col].rank()

    # Run ssGSEA analysis for the ranked gene expression dataset
    # txt, gct file input
    gene_sets = {}

    with open('./data/xCell/signatures.txt') as f:
        for line in f:
          key, val = line.split('\t', 1)
          val = val.split()
          gene_sets[key] = [str.lower(i) for i in val]

    print(f'Number of samples: {len(expr.columns)}, number of gene sets: {len(gene_sets)}')

    ssg = gp.ssgsea(data=expr,
                gene_sets=gene_sets, #gene_sets={'A':['gene1', 'gene2',...], 'B':['gene2', 'gene4',...],  ...}
                #outdir='test/ssgsea_report',
                sample_norm_method='custom', # choose 'custom' for your own rank list
                permutation_num=0, # skip permutation procedure, because you don't need it
                no_plot=True, # skip plotting, because you don't need these figures
                processes=20,
                seed=9,
                scale=False)

    scores = pd.DataFrame(ssg.resultsOnSamples)

    # Rescale on gene sets.
    scores = scores.subtract(scores.min(axis='columns'), axis='rows')

    # Combine signatures for same cell types
    scores['cell_type'] = scores.index
    scores['cell_type'] = scores.cell_type.apply(lambda x: x.split('%', 1)[0])
    scores = scores.groupby('cell_type').mean()

    fv = pd.read_csv('./data/xCell/spill_fv.csv', index_col=0)
    fv = fv.reindex(scores.index)

    #Rescale on cell types.
    scores = scores.subtract(scores.min(axis='columns'), axis='rows')
    scores = scores/5000

    #Use fv formula tscores <- (tscores^fit.vals[A,2])/(fit.vals[A,3]*2).
    scores = scores.pow(fv.V2, axis='rows')
    scores = scores.divide(fv.V3*2, axis='rows')

    #scores = spillOver(transformed.scores,xCell.data$spill.array$K)
    K = pd.read_csv('./data/xCell/spill_K.csv', index_col=0)
    K = K.reindex(scores.index)
    K = K[K.index]
    alpha = 0.5
    K = K * alpha

    def spillOver(sample):
      #Apply correction on samples: scores <- apply(transformedScores[rows, ], 2, function(x) pracma::lsqlincon(K[rows,rows], x, lb = 0)).
      x, rnorm = nnls(K, sample)
      return x

    scores = scores.apply(lambda x: spillOver(x), axis='rows')
    return scores


def process_expr(expr, gene_sets):
    for col in expr:
        expr[col] = expr[col].rank()

    ssg = gp.ssgsea(data=expr,
            gene_sets=gene_sets, #gene_sets={'A':['gene1', 'gene2',...], 'B':['gene2', 'gene4',...],  ...}
            #outdir='test/ssgsea_report',
            sample_norm_method='custom', # choose 'custom' for your own rank list
            permutation_num=0, # skip permutation procedure, because you don't need it
            no_plot=True, # skip plotting, because you don't need these figures
            processes=20,
            seed=9,
            scale=False)
    return pd.DataFrame(ssg.resultsOnSamples)

def msigdb(mix, pure):
    score_mix = {}
    score_pure = {}
    gene_sets = {}

    with open('./data/c8.all.v7.5.1.symbols.gmt') as f:
        for line in tqdm(f):
            key, val = line.split('\t', 1)
            val = val.split()
            val = [str.lower(i) for i in val]
            score_mix[key] = score(val, mix).total_score.values
            score_pure[key] = score(val, pure).total_score.values
            gene_sets[key] = val

    score_mix = pd.DataFrame(score_mix).T
    score_mix.columns = mix.columns
    score_mix.dropna(inplace=True)
    #score_mix = process_expr(mix, gene_sets)
    print('Score pure.')
    score_pure = pd.DataFrame(score_pure).T
    score_pure.columns = pure.columns
    score_pure.dropna(inplace=True)
    #score_pure = process_expr(pure, gene_sets)
    return(cibersort(score_mix, score_pure))