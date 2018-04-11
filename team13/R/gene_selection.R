# Function should take a normalised log counts matrix and return a row-subsetted matrix

scran_select = function(count_matrix, max.mean = 1, min.mean = 1e-3, pval = 0.01){
  if(as.numeric(R.version$minor) >= 5){ #TODO: FIX JANKINESS
    trend = scran::trendVar(count_matrix, loess.args = list(span = 0.05))
  } else {
    trend = scran::trendVar(count_matrix, span = 0.05)
  }
  decomp = scran::decomposeVar(count_matrix, fit = trend)  
  decomp = decomp[decomp$mean < max.mean & decomp$mean > min.mean, ]
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  genes = rownames(decomp)[decomp$FDR <= pval]
  return(count_matrix[rownames(count_matrix) %in% genes,])
}