# YADA Deconvolution Package
![Yada Flow](/data/Yada.jpg)

YADA is an innovative biological deconvolution algorithm developed by myself and my collaborator as part of my doctoral research in the Systems Biomedicine Lab under the supervision of Professor Efroni. Its purpose is to estimate the proportions of distinct cell types within complex, heterogeneous gene expression samples.

The fundamental premise behind YADA is that the transcriptomic signatures of pure cell populations can be leveraged to deconvolute mixed expression profiles and quantify the relative abundance of each constituent cell type. By analyzing gene expression patterns, deconvolution algorithms computationally unravel the complexities inherent to these cellular mixtures.

YADA implements a robust approach to perform this deconvolution task, accurately estimating immune and other cell type fractions from bulk transcriptomic data. It represents a novel contribution stemming from my doctoral studies focused on advancing computational methods for dissecting intricate systems-level biomedical data. Under the guidance of Professor Efroni, we were able to design, validate, and optimize YADA to address a crucial need in the field of computational immunology.

YADA implements two approaches:

1) Marker-based deconvolution using curated gene signatures representative of each cell type. 

2) Reference-based deconvolution utilizing a complete gene expression matrix from pure cell populations.

Key features and innovations of YADA include:
- High performance - YADA achieves state-of-the-art accuracy on benchmark datasets, as evidenced by top results in a recent challenge.
- Flexibility - Works with either marker genes or full reference profiles.
- Broad applicability - Core algorithm supports gene expression data from various sequencing platforms. 
- Computational efficiency - Optimized parallel implementation allows rapid deconvolution of large datasets.
- User-friendly Python API - One of the few deconvolution libraries natively implemented in Python.

In summary, YADA provides an accessible, high-speed Python toolkit for accurate and flexible deconvolution of bulk gene expression data, facilitating estimation of immune and other cell type mixtures from sample transcriptomes. Its user-friendly interface, performance, and multifaceted capabilities make YADA ideally suited for computational immunology and systems biology applications.

## Resources
- Matrix decomposition https://en.wikipedia.org/wiki/Matrix_decomposition
- Collection of papers https://github.com/changwn/Deconvolution_paper

## Training Dataset
- The training dataset for YADA comprises benchmark datasets available in the "data" folder. This comprehensive collection includes data from publicly available sources as well as synthetically generated datasets I created. These datasets were utilized for training and validation purposes during the DREAM challenge, a community-based deconvolution benchmarking effort.

## Input File Requirements
- The input data for YADA consists of two files:
    - pure.csv: A gene expression matrix for purified cell populations, with dimensions (n_genes) x (k_cell_types). For marker-based deconvolution, this file should contain only gene symbols. Refer to the sample notebook for formatting details.
    - mix.csv: A gene expression matrix for mixed cell samples, with dimensions (n_genes) x (m_mixtures). The first row must contain mixture labels.
    
Additional guidelines for input files:

- Gene symbols should be provided in the first column for both files.
It is acceptable for some genes to be missing in either the pure or mixed file.
- Expression data should be in non-log scale. If the maximum expression value is <50, an anti-log transformation is automatically applied.
- YADA performs internal marker gene selection; therefore, not all provided signature genes may be utilized.

By following these specifications, users can ensure their input data is properly formatted for YADA to perform accurate deconvolution.

## Sample Notebooks
- [Using reference matrix](/code/YADA-gene-diff.ipynb)
- [Using marker list](/code/YADA.ipynb)
- [YADA challenge result](/data/Challenge/challenge.ipynb)