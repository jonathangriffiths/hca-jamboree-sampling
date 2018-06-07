#TODO: local dimension estimation

get_k_distances = function(dim_red, k = 5){
  require(FNN)
  knns = get.knn(data = dim_red, k = k+1)
  return(knns$nn.dist[,-1])
}

get_kth_distances = function(dim_red, k = 5){
  dist = get_k_distances(dim_red)
  return(dist[,ncol(knns$dist)])
}

get_knn_densities = function(dim_red, k = 5){
  distances = get_kth_distances(dim_red, k=k)
  logvolumes = log(distances) * ncol(dim_red)
  return(-logvolumes)
}

normal_kernel_density = function(vec){
  #constants excluded, as these are always relative
  return(sum(exp(-1/2 * vec ^ 2)))
}

get_kernel_densities = function(dim_red, k = 50, kernel_function = normal_kernel_density){
  distances = get_k_distances(dim_red, k=k)
  densities = apply(distances, 1, kernel_function)
  return(densities)
}
