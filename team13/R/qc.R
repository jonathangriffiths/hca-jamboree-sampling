#Functions should take a count matrix and return a column-subsetted counts matrix
#after QC purposes.

qc_simple = function(count_matrix, mito.frac = 0.03, min.umis = 1000, max.umis=NULL) {
  count_matrix = count_matrix[, Matrix::colSums(count_matrix) >= min.umis]
  if (!is.null(max.umis)) {
    count_matrix = count_matrix[, Matrix::colSums(count_matrix) <= max.umis]
  }
  mt.genes = grepl("^MT-", rownames(count_matrix))
  mito_frac = Matrix::colSums(count_matrix[mt.genes,])/Matrix::colSums(count_matrix)
  count_matrix = count_matrix[,mito_frac < mito.frac]
  return(count_matrix)
}
