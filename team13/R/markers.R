select_cells = function(wrapper_out, 
                        method = c("quantile"),
                        quantile = 0.5){
  
  if(method == "quantile"){
    hits = wrapper_out$densities >= quantile(wrapper_out$densities, quantile)
  } else{
    "No suitable cell labelling method provided."
  }
  
  return(hits)
  
}

find_genes = function(hits, model_counts){
  require(glmnet)
  
  model = glmnet(x= t(model_counts), y = as.numeric(hits), family = "binomial", alpha = 1, dfmax = 6)
  
}