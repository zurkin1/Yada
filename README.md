# Yada Deconvolution Package.
![Yada Flow](/data/Yada.jpg)

Yada is an Python library for biological cell types deconvolution. Given gene expression data, a deconvolution algorithm is capable of estimating cell type proportions in mixture of cells. Yada is capable of deconvoluting either by using a list of marker genes or by using a complete pure gene expression matrix. Yada offers the following novelties:

- Performance: Yada’s results on benchmark datasets reached top results on a recent Dream challenge.
- Flexibility: Can be used with pure gene expression matrix or with marker-genes list only. Its core algorithm supports different sequencing platforms.
- Speed: Yada is very fast compared to other methods due to its parralel design.
- Yada is one of the few deconvolution tools that are based on Python.

## Resources.
- Matrix decomposition https://en.wikipedia.org/wiki/Matrix_decomposition
- Literature overview. https://urszulaczerwinska.github.io/UCzPhDThesis/methods.html#literature-overview

## Sample Datasets.
- Benchmark data sets are available in the data folder.

## Requirements on Input Datasets.
- Two files:
	- pure.csv: pure cell genes expression file. (n genes) x (k cell types)
	- mix.csv: mixtures genes expression file. (n genes) x (m mixtures).
	- Gene symbols in column 1; Mixture labels in row 1.
	- Tabular format with no missing entries.
	- It is OK if some genes are missing from the either file.
	- Data is assumed to be in non-log space (scale). If the dataset maximum expression value is less than 50, we run anti-log on all expression values.
- Yada performs a marker gene selection algorithm and therefore typically does not use all genes in the signature matrix. If this step is not needed a simple code change should comment the relevant lines.

## Sample Notebook.
- Use Yada.ipynb for sample run.