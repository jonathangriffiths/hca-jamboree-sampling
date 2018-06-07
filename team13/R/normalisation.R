# Functions should take a counts matrix and return it normalised and possibly variance stabilised (e.g. log)

#log2 normalised count + 1 (sparse)
normalise_scran = function(count_matrix){
  sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts = count_matrix))
  qclust= scran::quickCluster(x = sce, method = "igraph", max.size = 3000, pc.approx = TRUE)

  sce = scran::computeSumFactors(sce, clusters = qclust)
  sce = scater::normalise(sce)

  return(SingleCellExperiment::logcounts(sce))

}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#log2( (count + pseudocount) / geom_mean(across cell) )
normalise_clr = function(count_matrix, eps=1) {
    count_matrix = as.matrix(count_matrix) + eps
    geom_mean = apply(count_matrix, 2, gm_mean)
    norm_counts = count_matrix / matrix(geom_mean, nrow=nrow(count_matrix), ncol=ncol(count_matrix), byrow=TRUE)
    return(log2(norm_counts))
}

#log2 cpm
normalise_cpm = function(count_matrix){
  libs = Matrix::colSums(count_matrix)
  swept = sweep(count_matrix, 2, libs, "/") * 1e6
  return(log2(swept+1))
}
