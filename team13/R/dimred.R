#Functions should take a log normalised count matrix
#and return coordinates in a dimension-reduced space

dimred_pca = function(count_matrix, dims = 20, scale = FALSE){
  count_matrix = t(count_matrix)
  pca = irlba::prcomp_irlba(count_matrix, n = dims, scale. = scale)
  return(pca$x)
}

dimred_dmap = function(count_matrix, dims, scale = FALSE){
  require(destiny)
  count_matrix = t(count_matrix)
  if(scale) {
    count_matrix = scale(count_matrix)
  }
  dmap = destiny::DiffusionMap(count_matrix, n_eigs = dims)
  return(as.matrix(as.data.frame(dmap)))
}
