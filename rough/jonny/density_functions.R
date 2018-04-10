get_kth_distances = function(count_matrix, k = 5){
  require(FNN)
  knns = get.knn(data = t(count_matrix), k = k+1)
  return(knns$nn.dist[,ncol(nn_dist)])
}

