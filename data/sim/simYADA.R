# Set multiple parameters at once (using a list)
#https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html#41_Getting_and_setting
#https://www.bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatPop.html
suppressPackageStartupMessages({
  library(splatter)
  library(scater)
})
#params.batches <- newSplatPopParams(similarity.scale = 8,
#                                    batchCells = c(1000, 1000, 1000),
#                                    batch.size = 3,
#                                    batch.facLoc = c(0.5, 0.1,  0.1),
#                                    batch.facScale = c(0.6, 0.15, 0.15))

params <- newSplatParams(lib.loc = 12, lib.scale = 0.6)
sim <- splatSimulate(params, nGenes = 10000, group.prob = c(0.3, 0.3, 0.4), method = "groups", batchCells = 3000)
#sim <- logNormCounts(sim)
sim
#counts(sim)[1:5, 1:5]
#head(rowData(sim))
#head(colData(sim))
#names(assays(sim))
#assays(sim)$CellMeans[1:5, 1:5]
sim <- logNormCounts(sim)
sim <- runPCA(sim)
plotPCA(sim, colour_by = "Group")
#sim <- simpleSimulate(verbose = FALSE)
sim <- addGeneLengths(sim)
#head(colData(sim))
write.csv(colData(sim), "dedata.csv")
tpm(sim) <- calculateTPM(sim, rowData(sim)$Length)
tpm(sim)[1:5, 1:5]
write.csv(tpm(sim), "splater.csv")
write.csv(rowData(sim), "dedata.csv")
write.csv(assays(sim), "assays.csv")