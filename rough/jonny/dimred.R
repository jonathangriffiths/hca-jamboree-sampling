dimred_pca = function(count_matrix, dims, scale = FALSE){
  require(irlba)
  pca = prcomp_irlba(t(count_matrix), n = dims, scale. = scale)
  return(pca$x)
}