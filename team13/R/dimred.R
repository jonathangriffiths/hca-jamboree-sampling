#Functions should take a log normalised count matrix
#and return coordinates in a dimension-reduced space

dimred_pca = function(count_matrix, dims, scale = FALSE){
  require(irlba)
  pca = prcomp_irlba(Matrix::t(count_matrix), n = dims, scale. = scale)
  return(pca$x)
}

dimred_dmap = function(count_matrix, dims, scale = FALSE){
  require(destiny)
  if(scale)
    count_matrix = Matrix::t(scale(Matrix::t(count_matrix)))
  dmap = DiffusionMap(t(as.matrix(count_matrix)), n_eigs = dims)
  return(as.matrix(as.data.frame(dmap)))
}