import argparse
import os, sys
from scsim import scsim
import numpy as np
import pandas as pd
import time



def save_df(obj, filename):
    '''Save pandas dataframe in compressed format'''
    #np.savez_compressed(filename, data=obj.values, index=obj.index.values, columns=obj.columns.values)
    obj.to_csv(filename + '.csv')


def load_df(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj
    
    
def parse_args():
    parser = argparse.ArgumentParser(description='Run scsim with specified input arguments')
    parser.add_argument('--outdir', type=str, default='./', #scsim-%s-%s-%s-%s-%s-%s-%s-%s',
                        help='Output directory base')
    parser.add_argument('--seed', type=int, help='simulation seed')
    parser.add_argument('--numsims', type=int, help='number of sims to run',
                        default=20)
    parser.add_argument('--deloc', type=float,
                        help='devalue',
                        default=1.)
    parser.add_argument('--K', type=int,
                        help='Number of identity programs',
                        default=10)
    parser.add_argument('--nproggoups', type=int,
                        help='Number of groups expressing activity program. Default is 1/3 of K rounded down',
                        default=None)
    parser.add_argument('--ncells', type=int,
                        help='Total number of cells',
                        default=100)
    parser.add_argument('--doubletfrac', type=float,
                        help='Percentage of doublet cells',
                        default=0.)
    a = parser.parse_args()
    return(a.outdir, a.seed, a.numsims, a.deloc, a.K, a.nproggoups, a.ncells, a.doubletfrac)





def main(test):
    (outdir, randseed, numsims, deval, K, nproggroups, ncells, doubletfrac) = parse_args()
    outdir = f'./data/{test}-'
    ngenes=1000
    nproggenes = 400
    ndoublets=int(doubletfrac*ncells)
    
    deloc=deval
    progdeloc=deval
    descale=1.0
    progcellfrac = .35
    deprob = test*5/100 #.0025 #Of all genes in group what percentage are differentially expressed.

    if nproggroups is None:
        nproggroups = int(K/3)
        
    proggroups = list(range(1, nproggroups+1))

    simulator = scsim(ngenes=ngenes, ncells=ncells, ngroups=K, libloc=7.64, libscale=0.78,
                 mean_rate=7.68,mean_shape=0.34, expoutprob=0.00286,
                 expoutloc=6.15, expoutscale=0.49,
                 diffexpprob=deprob, diffexpdownprob=0., diffexploc=deloc, diffexpscale=descale,
                 bcv_dispersion=0.448, bcv_dof=22.087, ndoublets=ndoublets,
                 nproggenes=nproggenes, progdownprob=0., progdeloc=progdeloc,
                 progdescale=descale, progcellfrac=progcellfrac, proggoups=proggroups,
                 minprogusage=.1, maxprogusage=.7, seed=randseed)


    start_time = time.time()
    simulator.simulate()
    end_time = time.time()
    #print('%.3f minutes elapsed for seed %d' % ((end_time-start_time)/60, randseed))

    save_df(simulator.cellparams, outdir + 'cellparams')
    save_df(simulator.geneparams, outdir + 'geneparams')
    save_df(simulator.counts.T, outdir + 'pure')

    #Create mixes.
    prop_data = np.empty((0, 6))
    data = simulator.counts.T.iloc[:,:6]
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
    #data2 = data2.apply(lambda x: np.log2(x))
    data2.to_csv(f'./data/{test}-mix.csv')
    prop_df.to_csv(f'./data/{test}-prop.csv')
    #data = data.apply(lambda x: np.log2(x))
    #data.to_csv(f'c:/data/Geo/Microarray/mix-6244v3/pure.csv')


if __name__ == '__main__':
    for test in range(20):
        main(test)