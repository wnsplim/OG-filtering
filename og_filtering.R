library(dplyr)

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

og_filter <- dget("./functions/og_filter.R")
## Usage: og_filter(og, threshold, n, outgroup, statistics, quiet)
# og: OrthologousMatrix.txt output file of OMA
# threshold: Threshold value
# n: Optional, consider OG with at least n number of taxa, default = minimum OG size of og
# outgroup: Optional, vector of an outgroup species to prune before calculation, default = NULL
# statistics: Optional, if TRUE return some statistics that can be used to determine a threshold value, default = FALSE
# quiet: Optional, if FALSE print summary statistics, default = FALSE

threshold <- dget("./functions/threshold.R")
## Usage: threshold(og, maxthres, diff, outgroup)
# og: OrthologousMatrix.txt output file of OMA
# maxthres: Maximum number of threshold value to consider
# diff: Increment of the threshold sequence
# outgroup: Optional, vector of an outgroup species to prune before calculation, default = NULL


og_Bibionomorpha <- read.delim("./data/Bibionomorpha_OrthologousMatrix.txt", skip = 4, header = TRUE)

## Bibionomorpha OGset retaining outgroup

threshold_og1 <- threshold(og_Bibionomorpha, 30, 0.1)
og_filtered1 <- og_filter(og_Bibionomorpha, 5.3)
og_bibionomorpha <- rownames(og_filtered1[[1]])
writeLines(og_bibionomorpha, "./result/og_Bibionomorpha.txt")

## Bibionomorpha OGset excluding outgroup (used in the main study)

Brachycera <- c("Drosophila_melanogaster_AA")

threshold_og2 <- threshold(og_Bibionomorpha, 30, 0.1, outgroup = Brachycera)
og_filtered2 <- og_filter(og = og_Bibionomorpha, threshold = 8.8, outgroup = Brachycera)
og_bibionomorpha_noOut <- rownames(og_filtered2[[1]])
writeLines(og_bibionomorpha_noOut, "./result/og_Bibionomorpha_noOut.txt")
