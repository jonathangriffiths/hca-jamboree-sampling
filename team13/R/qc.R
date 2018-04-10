qc_simple = function(count_matrix, mito.frac = 0.03, min.umis = 1000, gene_df = genes){
  count_matrix = count_matrix[, Matrix::colSums(count_matrix) >= min.umis]
  mt.genes = grepl("^MT-", gene_df[,1])
  mito_frac = Matrix::colSums(count_matrix[mt.genes,])/Matrix::colSums(count_matrix)
  count_matrix = count_matrix[,mito_frac < mito.frac]
  return(count_matrix)
}
