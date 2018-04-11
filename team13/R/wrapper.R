density_wrapper = function(count_matrix,
                   qc_function = c("simple"),
                   qc_params = list(),
                   normalisation_method = c("scran", "clr", "cpm", "none"),
                   gene_selection_method = "scran",
                   dimred_method = "pca",
                   density_k = 10){

  qc_funcs <- list()
  qc_funcs$simple = qc_simple

  qc_function = qc_function[1]
  if (!(qc_function %in% names(qc_funcs))) {
      stop("Unknwon QC method chosen:", qc_function)
  }
  qc_params = c(list(count_matrix), qc_params)
  do.call(qc_funcs[[qc_function]], qc_params)

  #normalise
  if(normalisation_method[1] == "scran"){
    norm_counts = normalise_scran(count_matrix)
  } else if(normalisation_method[1] == "clr"){
    norm_counts = normalise_clr(count_matrix)
  } else if(normalisation_method[1] == "cpm"){
    norm_counts = normalise_cpm(count_matrix)
  } else if(normalisation_method[1] == "none"){
    print("Hmmmmm.....")
    norm_counts = log2(count_matrix + 1)
  } else {
    stop("No suitable normalisation method chosen.")
  }

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
