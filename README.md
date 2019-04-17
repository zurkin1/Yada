# Yada Package

- This is the Yada deconvolution package.
- Run it using the main function in run.py.
- run.py runs an ensemble of methods to find the best approximation to cell proportions.
- Ten data sets are available under the data folder, together with the true proportions and CIBERSORT result.
	- TIMER
	- PertU
	- DSA
	- DeconRNASeq
	- Abbas
	- BreastBlood
	- RatBrain
	- 10x
	- EPIC
	- CIBERSORT
- Use the first line in run.py to select the datasets to run.
- The module will run 3 ensemble algorithms and present the one that achieves the best score on each dataset.

# Requirements on datasets:
- Two files:
	- pure.csv: pure cell genes expression file. n(genes) x k(cell types)
	- mix.csv: mixtures genes expression file. n(genes) x m(mixtures). Gene symbols in column 1; Mixture labels in row 1.
- Tabular format with no missing entries.
- It is OK if some genes are missing from the either file since we run a joined list operation.
- Data should be in non-log space. If the dataset has maximum expression value is <50 it is
 given in log space, and we run anti-log on all expression values by 2^x.
- Yada performs a marker gene selection algorithm and therefore typically does not use all genes in the signature matrix.