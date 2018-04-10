get_kth_distances = function(dim_red, k = 5){
  require(FNN)
  knns = get.knn(data = t(dim_red), k = k+1)
  return(knns$nn.dist[,ncol(nn_dist)])
}

get_cell_densities = function(dim_red, k = 5){
  distances = get_kth_distances(dim_red, k=k)
  volumes = distances ^ ncol(dim_red)
  return(1/volumes)
}