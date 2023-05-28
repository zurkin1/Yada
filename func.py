import numpy as np
import pandas as pd
from diffexp import gene_diff
from func import *
import scipy.stats as st
import time
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
import gseapy as gp
#from singscore import *
#from tqdm import tqdm
import multiprocessing as mp


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
    #expr = pd.read_csv(expr, index_col=0)
    # Reduce the expression dataset to contain only the required genes
    expr.index = expr.index.map(str.lower)
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


if __name__ == '__main__':
    mix = pd.read_csv('./data/xCell/chicago.csv', index_col=0)
    result = pxcell(mix)
    result.to_csv('./data/xCell/results.csv')