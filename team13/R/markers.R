select_cells = function(wrapper_out, 
                        method = c("quantile", "mixmod"),
                        quantile = 0.5){
  #note that hits refers to a logical vector for low density cells
  if(method == "quantile"){
    hits = wrapper_out$densities <= quantile(wrapper_out$densities, quantile)
  } else if(method == "mixmod"){
    require(mixtools)
    mixmod = normalmixEM(wrapper_out$densities, mu = c(min(wrapper_out$densities), max(wrapper_out$densities)))
    small_dist = which.min(mixmod$mu)
    hits = mixmod$posterior[,small_dist] > mixmod$posterior[,-small_dist]
  } else{
    "No suitable cell labelling method provided."
  }
  
  return(hits)
  
}

find_genes = function(hits, model_counts, sample_fraction = 0.7, runs = 10, n_genes = 6){
  require(glmnet)
  
  list = lapply(1:runs, function(x){
    retain = sample(length(hits), round(length(hits)*sample_fraction))
    model = glmnet(x= t(model_counts[,retain]), y = as.numeric(hits[retain]), family = "binomial", alpha = 1, dfmax = 100)
    hits = which(model$beta[, ncol(model$beta)-1]>0)
    return(names(hits))
  })
  
  select = data.frame(gene = unique(unlist(list)),
                      freq = sapply(unique(unlist(list)), function(x) sum(unlist(list)==x)))
  
  select$mean.position = sapply(select$gene, function(x){
    posns = sapply(list, function(part){
      if(x %in% part){
        return(which(part == x))
      } else {
        return(NA)
      }
    }
    )
    
    
    return(mean(posns, na.rm = TRUE))
  })
  
  select$avg = sapply(1:nrow(select), function(i) mean(c(select$freq[i], 50-select$mean.position[i])))
  select = select[order(-select$avg),]

  return(select$gene[1:n_genes])
  
}



make_gate = function(genes, densities, count_matrix){
  
}
