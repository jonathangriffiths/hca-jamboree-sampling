get_kth_distances = function(dim_red, k = 5){
  require(FNN)
  knns = get.knn(data = dim_red, k = k+1)
  return(knns$nn.dist[,ncol(knns$nn.dist)])
}

get_cell_densities = function(dim_red, k = 5){
  distances = get_kth_distances(dim_red, k=k)
  logvolumes = log(distances) * ncol(dim_red)
  return(-logvolumes)
}