### SELECT CELLS FUNCTIONS
# INPUT cell density values
# OUTPUT a list:
#   hits = logical vector of cell selection
#   plot = diagnostic plot of selected cells.

select_cells_lowdensity_quantile = function(densities, quantile = 0.5, quantile_type = c("upper", "lower")){
  hits = densities <= quantile(densities, quantile)
  if(quantile_type == "upper")
    hits = !hits
  
  plot = ggplot(mapping = aes(x = densities, fill = hits)) +
    geom_histogram(bins = 50) +
    scale_fill_brewer(palette = "Paired", name = "Selected\ncell") +
    labs(x = "Log-density", y = "Count")
  
  return(list(hits = hits, plot = plot))
}

select_cells_mixmod = function(densities, quantile_type = c("lower", "upper")){
  require(mixtools)
  mixmod = normalmixEM(densities, mu = c(min(densities), max(densities)))
  small_dist = which.min(mixmod$mu)
  hits = mixmod$posterior[,small_dist] > mixmod$posterior[,-small_dist]
  if(quantile_type == "upper")
    hits = !hits
  
  plot = ggplot(mapping = aes(x = densities, fill = hits)) +
    geom_histogram(bins = 50) +
    scale_fill_brewer(palette = "Paired", name = "Selected\ncell") +
    labs(x = "Log-density", y = "Count")
  
  return(list(hits = hits, plot = plot))
  
}

### FIND GENES FUNCTIONS
# INPUT:
#   hits - logical vector of low density cells; 
#   model_counts - expression matrix (e.g. logcounts) for genes usable for sorting e.g. surface proteins
# OUTPUT

find_genes_lasso = function(hits, model_counts, sample_fraction = 0.7, runs = 10, n_genes = 6){
  require(glmnet)
  
  main_model = glmnet(x= t(model_counts), y = as.numeric(hits), family = "binomial", alpha = 1, dfmax = n_genes)
  main_hits = rownames(model_counts)[which(main_model$beta[, ncol(main_model$beta)-1]!=0)]
  
  #here we repear the glmnet run on randomly subsampled data
  #this is a sort of robustness estimate on the lasso selections
  subsamp = lapply(1:runs, function(x){
    retain = sample(length(hits), round(length(hits)*sample_fraction))
    model = glmnet(x= t(model_counts[,retain]), y = as.numeric(hits)[retain], family = "binomial", alpha = 1, dfmax = n_genes*2)
    hits = which(model$beta[, ncol(model$beta)-1]>0)
    return(names(hits))
  })
  
  

  select = data.frame(gene = unique(c(unlist(subsamp), main_hits)))
  select$subsample.freq = sapply(select$gene, function(x) sum(unlist(subsamp)==x))

  select$mean.position = sapply(select$gene, function(x){
    posns = sapply(subsamp, function(part){
      if(x %in% part){
        return(which(part == x))
      } else {
        return(NA)
      }
    }
    )
    return(mean(posns, na.rm = TRUE))
  })
  
  select$main.hits = select$gene %in% main_hits
  
  select = select[order(-select$subsample.freq, select$mean.position),]

  return(list(genes = select$gene[1:n_genes], model = main_model))
  
}



make_gate = function(genes, densities, count_matrix){
  
}
