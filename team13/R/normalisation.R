#log2 normalised count + 1 (sparse)
normalise_scran = function(count_matrix){
  require(SingleCellExperiment)
  require(scran)
  require(scater)
  
  sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts = count_matrix))
  qclust= scran::quickCluster(x = sce, method = "igraph", max.size = 3000, pc.approx = TRUE)
  
  sce = scran::computeSumFactors(sce, clusters = qclust)
  sce = scater::normalise(sce)
  
  return(SingleCellExperiment::logcounts(sce))
    
}

#log2( (count + pseudocount) / geom_mean(across cell) )
normalise_clr = function(count_matrix){
  geom_mean = apply(count_matrix, 2, function(x){exp(sum(log(x+1))/length(x))})
  swept = sweep(count_matrix + 1, 2, geom_mean, "/")
  return(log2(swept))
}

#log2 cpm
normalise_cpm = function(count_matrix){
  libs = Matrix::colSums(count_matrix)
  swept = sweep(count_matrix, 2, libs, "/") * 1e6
  return(log2(swept+1))
}
