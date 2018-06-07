#!/usr/bin/env Rscript
library(methods)
library(Matrix)
library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(Rtsne)
library(RANN)
library(glmnet)
devtools::load_all('../team13')


# Read raw counts
counts = readMM("../data/ica_bone_marrow_ref_set.mtx")
genes = read.table("../data/ica_bone_marrow_ref_set.genes.tsv")
cells = read.table("../data/ica_bone_marrow_ref_set.barcodes.tsv")
rownames(counts) = genes[,1]
colnames(counts) = cells[,1]
counts = as(counts, "dgCMatrix")


# Perform QC, CLR normalization and scran gene selection
preprocess = density_wrapper(counts,
                             qc_method="simple",
                             qc_params=list(mito.frac=0.05, max.umis=7500),
                             normalisation_method="clr",
                             gene_selection_method="scran",
                             dimred_method="pca",
                             dimred_params=list(dims=50),
                             density_method="none")

norm_reads = preprocess$retained_counts
dimred = preprocess$dimred


# Create 2D projection using tsne
tsne = Rtsne(dimred, pca=FALSE, check_duplicates=FALSE)
projection = tibble(x=tsne$Y[,1], y=tsne$Y[,2])

ggp <- ggplot(projection, aes(x=x, y=y)) +
    geom_point(size=0.3) +
    theme_bw()

ggsave('tsne_projection.png', ggp, width=5, height=5, units='in')


# Calculate for each cell 50 nn using euclidean distance
max_k <- 50
knn <- nn2(dimred, k=max_k+1)
k50 <- knn$nn.dists[,50+1]


# Plot distance of 50th nn over tsne
ggdata <- projection %>%
          mutate(k50=k50)

ggp <- ggplot(ggdata, aes(x=x, y=y, col=log2(k50))) +
    geom_point(size=0.3) +
    scale_color_gradient2(low='#1a9850', high='#d73027', mid='#d9ef8b', midpoint=3, space='Lab') +
    theme_bw()

ggsave('tsne_log_k50.png', ggp, width=6, height=5, units='in')


# Subsample half the cell
# Calculate improvement in distance of 50th nn of the full sample compared to the subsample
repeats <- 1000

diffs <- list()
for (i in 1:repeats) {
    ncell <- nrow(dimred)
    idxs <- sample(1:ncell, floor(ncell/2), replace=FALSE)
    subsample <- dimred[idxs,]
    before <- nn2(subsample, k=max_k+1)

    diff_matrix <- array(NA, dim=c(ncell, max_k+1))
    diff_matrix[idxs,] <- before$nn.dists - knn$nn.dists[idxs,]
    diffs <- c(diffs, list(diff_matrix))
}


# Calculate knn distance improvement for different valus of k
d1  <- sapply(1:repeats, function(i) {diffs[[i]][,1+1]}) %>% rowMeans(na.rm=TRUE)
d10 <- sapply(1:repeats, function(i) {diffs[[i]][,10+1]}) %>% rowMeans(na.rm=TRUE)
d25 <- sapply(1:repeats, function(i) {diffs[[i]][,25+1]}) %>% rowMeans(na.rm=TRUE)
d50 <- sapply(1:repeats, function(i) {diffs[[i]][,50+1]}) %>% rowMeans(na.rm=TRUE)


# Plot histograms knn distance imrovements
ggdata <- tibble(d1=d1, d10=d10, d25=d25, d50=d50)
ggp1 <- ggplot(ggdata, aes(x=d1)) +
    geom_histogram(binwidth=0.05) +
    theme_bw()
ggp10 <- ggplot(ggdata, aes(x=d10)) +
    geom_histogram(binwidth=0.05) +
    theme_bw()
ggp25 <- ggplot(ggdata, aes(x=d25)) +
    geom_histogram(binwidth=0.05) +
    theme_bw()
ggp50 <- ggplot(ggdata, aes(x=d50)) +
    geom_histogram(binwidth=0.05) +
    theme_bw()
ggp <- plot_grid(ggp1, ggp10, ggp25, ggp50, ncol=2)
ggsave('knn_diff_histogram.png', ggp, width=10, height=10, units='in')


# Plot improvement of 50th nn distance over tsne
ggdata <- projection %>%
          mutate(d50=d50)

ggp <- ggplot(ggdata, aes(x=x, y=y, col=pmin(d50, 2.5))) +
    geom_point(size=0.3) +
    scale_color_gradient2(low='#1a9850', high='#d73027', mid='#d9ef8b', midpoint=1, space='Lab') +
    theme_bw()

ggsave('tsne_diff50.png', ggp, width=6, height=5, units='in')


# Use logistic regression to identify "dense" regions
d50_thresh <- 0.8
Y <- 2*(d50 < d50_thresh) - 1
X <- t(norm_counts)

fit <- glmnet(X, Y, family='binomial', alpha=1)
