density_wrapper = function(count_matrix,
                   qc_method = c("simple"),
                   qc_params = list(),
                   normalisation_method = c("scran", "clr", "cpm", "none"),
                   normalisation_params = list(),
                   gene_selection_method = "scran",
                   dimred_method = "pca",
                   density_k = 10){

  #QC
  qc_funcs = list()
  qc_funcs$simple = qc_simple

  qc_func = qc_funcs[[qc_method[1]]]
  if (is.null(qc_func)) {
      stop("Unknwon QC method chosen:", qc_method[1])
  }
  qc_params = c(list(count_matrix), qc_params)
  count_matrix = do.call(qc_func, qc_params)

  #normalise
  norm_funcs = list()
  norm_funcs$scran = normalise_scran
  norm_funcs$clr   = normalise_clr
  norm_funcs$cpm   = normalise_cpm
  norm_funcs$none  = function(counts, ...) counts

  norm_func = norm_funcs[[normalisation_method[1]]]
  if (is.null(norm_func)) {
      stop("Unknwon normalisation method chosen:", normalisation_method[1])
  }
  normalisation_params = c(list(count_matrix), normalisation_params)
  norm_counts = do.call(norm_func, normalisation_params)

  #gene selection
  if(gene_selection_method == "scran"){
    norm_counts = scran_select(norm_counts)
  } else {
    stop("No suitable gene selection method chosen.")
  }

  #dimension reduction
  if(dimred_method == "pca"){
    dimred = dimred_pca(as.matrix(norm_counts), dims = 20)
  } else{
    stop("No suitable dimension reduction method chosen.")
  }
}

select_markers = function(){
  #TODO: Bo.
}
