# Yada Deconvolution Package

- Ten benchmark data sets are available in the data folder, together with the true proportions (TIMER, PertU, DSA, DeconRNASeq, Abbas, BreastBlood, RatBrain, EPIC, CIBERSORT).
- Yada supports ensemble of algorithms. Follow the instructions for its configuration.

# Requirements on Dataset:
- Two files:
	- pure.csv: pure cell genes expression file. n(genes) x k(cell types)
	- mix.csv: mixtures genes expression file. n(genes) x m(mixtures). Gene symbols in column 1; Mixture labels in row 1.
- Tabular format with no missing entries.
- It is OK if some genes are missing from the either file.
- Data is assumed to be in non-log space. If the dataset maximum expression value is less than 50, we run anti-log on all expression values.
- Yada performs a marker gene selection algorithm and therefore typically does not use all genes in the signature matrix.

# Running on Anaconda
- pip install tslearn
- git clone https://github.com/zurkin1/Yada.git
- Using Jupyter notebook open Yada.ipynb and follow the instructions.