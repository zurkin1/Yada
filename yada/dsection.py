import numpy as np
import pandas as pd
import pymc3 as pm

#This function calculates a deconvolution algorithm using Bayes algorithm.
def dsection(params):
    #Use only differential meaningful genes.
    gene_list = pd.concat([params['gene_list_df'].iloc[:, i] for i in range(len(params['gene_list_df'].columns))]).values
    gene_list = [x for x in gene_list if str(x) != 'nan']
    mix = params['mix'].loc[gene_list]
    pure = params['pure'].loc[gene_list]
    p_start = np.ones([len(pure.columns), len(mix.columns)])
    # Variance calculation
    joined_data = pd.concat([pure, mix], axis=1)
    joined_std = pd.DataFrame([joined_data.std(axis=1)] * len(pure.columns)).T

    basic_model = pm.Model()
    with basic_model:
        # Priors
        gamma = pm.Gamma('gamma', alpha=1, beta=1)
        ϵ = pm.Deterministic('ϵ', gamma + joined_std.T)
        pure_ = pm.Normal('pure_', mu=pure.T, shape=pure.T.shape, sd=ϵ)
        p_ = pm.Dirichlet('p_', a=p_start.T, shape=p_start.T.shape)

        # Likelihood (sampling distribution) of observations
        gamma2 = pm.Gamma('gamma2', alpha=1, beta=1)
        mu = pm.math.matrix_dot(p_, pure_)  # np.dot(pure_, p_) # + epsilon
        Y_obs = pm.Normal('Y_obs', mu=mu, observed=mix.T, sd=gamma2)

        start = pm.find_MAP(model=basic_model)
        step = pm.Metropolis()
        trace = pm.sample(1000, tune=50, step=step, cores=1, start=start, chains=4)

    return(np.round(trace['p_'].mean(axis=0), 2).T)