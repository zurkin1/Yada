EPIC package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Description
-----------

Package implementing EPIC method to estimate the proportion of immune, stromal, endothelial and cancer or other cells from bulk gene expression data. It is based on reference gene expression profiles for the main non-malignant cell types and it predicts the proportion of these cells and of the remaining "*other cells*" (that are mostly cancer cells) for which no reference profile is given.

This method is described in the publication from *Racle et al., 2017* available at <https://elifesciences.org/articles/26476>.

Usage
-----

The main function in this package is `EPIC`. It needs as input a matrix of the TPM (or RPKM) gene expression from the samples for which to estimate cell proportions. One can also define the reference cells to use

``` r
# library(EPIC) ## If the package isn't loaded (or use EPIC::EPIC and so on).
out <- EPIC(bulk = bulkSamplesMatrix)
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList)
```

`out` is a list containing the various mRNA and cell fractions in each samples as well as some *data.frame* of the goodness of fit.

Values of mRNA per cell and signature genes to use can also be changed:

``` r
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell = mRNA_cell_vector, sigGenes = sigGenes_vector)
out <- EPIC(bulk = bulkSamplesMatrix, reference = referenceCellsList, mRNA_cell_sub = mRNA_cell_sub_vector)
```

Various other options are available and are well documented in the help pages from EPIC:

``` r
?EPIC::EPIC
?EPIC::EPIC.package
```

Installation
------------

``` r
install.packages("devtools")
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
```

Python wrapper
--------------

A pyhton wrapper has been written by Stephen C. Van Nostrand von MIT and is available at <https://github.com/scvannost/epicpy>.

License
-------

EPIC can be used freely by academic groups for non-commercial purposes. The product is provided free of charge, and, therefore, on an "*as is*" basis, without warranty of any kind. Please read the file "*LICENSE*" for details.

If you plan to use EPIC (version 1.1) in any for-profit application, you are required to obtain a separate license. To do so, please contact Ece Auffarth (<eauffarth@licr.org>) at the Ludwig Institute for Cancer Research Ltd.

Contact information
-------------------

Julien Racle (<julien.racle@unil.ch>), and David Gfeller (<david.gfeller@unil.ch>).

FAQ
---

##### What do the "*other cells*" represent?

-   EPIC predicts the proportions of the various cell types for which we have gene expression reference profiles (and corresponding gene signatures). But, depending on the bulk sample, it is possible that some other cell types are present for which we don't have any reference profile. EPIC returns the proportion of these remaining cells under the name "*other cells*". In the case of tumor samples, most of these other cells would certainly correspond to the cancer cells, but it could be that there are also some stromal cells or epithelial cells for example.

##### I receive an error message "*attempt to set 'colnames' on an object with less than two dimensions*". What can I do?

-   This is certainly that some of your data is a vector instead of a matrix. Please make sure that your bulk data is in the form of a matrix (and also your reference gene expression profiles if using custom ones).

##### I receive a warning message that "*the optimization didn't fully converge for some samples*". What does it mean?

-   When estimating the cell proportions EPIC performs a least square regression between the observed expression of the signature genes and the expression of these genes predicted based on the estimated proportions and gene expression reference profiles of the various cell types.

    When such a warning message appears, it means that the optimization didn’t manage to fully converge for this regression, for some of the samples. You can then check the "*fit.gof$convergeCode*" (and possibly also "*fit.gof$convergeMessage*") that is outputted by EPIC alongside the cell proportions. This will tell you which samples had issue with the convergence (a value of 0 means it converged ok, while other values are errors/warnings, their meaning can be found in the help of "*optim*" (or "*constrOptim*") function from R (from "*stats*" package) which is used during the optimization and we simply forward the message it returns).

    The error code that usually comes is a "1" which means that the maximum number of iterations has been reached in the optimization. This could mean there is an issue with the bulk gene expression data that maybe don’t completely follow the assumption of equation (1) from our manuscript. In practice, I’ve observed that even when there was such an error message the proportions were predicted well, it is maybe that the optimization just wants to be *too precise*, or maybe few of the signature genes didn’t match well but the rest of signature genes could be used to have a good estimate of the proportions.

    If you have some samples that seem to have strange results, it could however be useful to check that the issue is not that these samples didn’t converge well. To be more conservative you could also remove all the samples that didn't converge well as these are maybe outliers, if it is only a small fraction from your original samples. Another possibility would be to change the parameters of the optim/constrOptim function to allow for more iterations or maybe a weaker tolerance for the convergence, but for this you would need to tweak it directly in the code of EPIC, I didn't implement such option for EPIC.

##### Who should I contact in case of a technical or other issue?

-   Julien Racle (<julien.racle@unil.ch>). Please provide as much details as possible and ideally send also an example input file (and/or reference profiles) that is causing the issue.
