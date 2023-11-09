# Yada Deconvolution Package
![Yada Flow](/data/Yada.jpg)

Yada is a Python library designed for deconvoluting mixed gene expression data to estimate cell type proportions. Deconvolution algorithms leverage the expression signature of pure cell populations to quantify their relative abundances in complex, heterogeneous samples. Yada implements two approaches:

1) Marker-based deconvolution using curated gene signatures representative of each cell type. 

2) Reference-based deconvolution utilizing a complete gene expression matrix from pure cell populations.

Key features and innovations of Yada include:
- High performance - Yada achieves state-of-the-art accuracy on benchmark datasets, as evidenced by top results in a recent challenge.
- Flexibility - Works with either marker genes or full reference profiles.
- Broad applicability - Core algorithm supports gene expression data from various sequencing platforms. 
- Computational efficiency - Optimized parallel implementation allows rapid deconvolution of large datasets.
- User-friendly Python API - One of the few deconvolution libraries natively implemented in Python.

In summary, Yada provides an accessible, high-speed Python toolkit for accurate and flexible deconvolution of bulk gene expression data, facilitating estimation of immune and other cell type mixtures from sample transcriptomes. Its user-friendly interface, performance, and multifaceted capabilities make Yada ideally suited for computational immunology and systems biology applications.

## Resources
- Matrix decomposition https://en.wikipedia.org/wiki/Matrix_decomposition
- Literature overview https://urszulaczerwinska.github.io/UCzPhDThesis/methods.html#literature-overview
- Collection of papers https://github.com/changwn/Deconvolution_paper

## Dataset For Training
- Benchmark data sets are available in the data folder. This is the dataset I used for training to the DREAM challenge. It is collected from open datasets as well as some synthetic datasets that I have created.

## Requirements on Input Files
- pure.csv: A gene expression matrix for pure cell populations, with dimensions of (n genes) x (k cell types). For marker-based deconvolution, this file contains only gene symbols. Refer to the sample notebook for formatting details. 
- mix.csv: A gene expression matrix for mixed cell samples, with dimensions of (n genes) x (m mixtures). Row 1 contains mixture labels. Additional input guidelines:
- Gene symbols in column 1 of both files.
- It is acceptable for some genes to be missing in either file. 
- Data should be in non-log scale. If maximum expression value is <50, anti-log transform is applied. 
- Yada implements internal marker gene selection - not all signature genes are used.

## Sample Notebooks
- [Using reference matrix](Yada.ipynb)
- [Using only marker list](Yada-only_markers.ipynb)