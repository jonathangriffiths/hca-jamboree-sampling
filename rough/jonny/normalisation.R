#log2 normalised count + 1 (sparse)
normalise_scran = function(count_matrix){
  require(SingleCellExperiment)
  require(scran)
  
  sce = SingleCellExperiment(assays = list(counts = count_matrix))
  qclust= quickCluster(x = sce, method = "igraph", max.size = 3000, pc.approx = TRUE)
  
  sce = computeSumFactors(sce, clusters = qclust)
  sce = normalise(sce)
  
  return(logcounts(sce))
    
}

#log2( (count + pseudocount) / geom_mean(across cell) )
normalise_clr = function(count_matrix){
  
  geom_mean = apply(count_matrix, 2, function(x){product(x+1)^(1/length(x))})
  swept = sweep(count_matrix + 1, 2, geom_mean, "/")
  return(log2(swept))
}

#log2 cpm
normalise_cpm = function(count_matrix){
  swept = sweep(count_matrix + 1, 2, Matrix::colSums(count_matrix) + 1, "/") * 1e6
  return(log2(swept))
}
