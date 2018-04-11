density_wrapper = function(count_matrix,
                   qc_method = "simple",
                   qc_params = list(),
                   normalisation_method = c("scran", "clr", "cpm", "none"),
                   normalisation_params = list(),
                   gene_selection_method = "scran",
                   gene_selection_params = list(),
                   dimred_method = "pca",
                   dimred_params = list(),
                   density_method = c("knn", "none"),
                   density_params = list())

  #QC
  qc_funcs = list()
  qc_funcs$simple = qc_simple

  qc_func = qc_funcs[[qc_method[1]]]
  if (is.null(qc_func)) {
      stop("Unknwon QC method chosen: '", qc_method[1], "'")
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
      stop("Unknwon normalisation method chosen: '", normalisation_method[1], "'")
  }
  normalisation_params = c(list(count_matrix), normalisation_params)
  norm_counts = do.call(norm_func, normalisation_params)

  #gene selection
  gene_selection_funcs = list()
  gene_selection_funcs$scran = scran_select

  gene_selection_func = gene_selection_funcs[[gene_selection_method[1]]]
  if (is.null(gene_selection_func)) {
      stop("Unknwon gene selection method chosen: '", gene_selection_method[1], "'")
  }
  gene_selection_params = c(list(norm_counts), gene_selection_params)
  norm_counts = do.call(gene_selection_func, gene_selection_params)

  #dimension reduction
  dimred_funcs = list()
  dimred_funcs$pca = dimred_pca

  dimred_func = dimred_funcs[[dimred_method[1]]]
  if (is.null(dimred_func)) {
      stop("Unknwon dimension reduction method chosen: '", dimred_method[1], "'")
  }
  dimred_params = c(list(norm_counts), dimred_params)
  dimred = do.call(dimred_func, dimred_params)

  #densities
  density_funcs = list()
  density_funcs$knn = get_knn_densities
  density_funcs$none = function(dimred, ...) { return NULL }

  density_func = density_funcs[[density_method[1]]]
  if (is.null(density_func)) {
      stop("Unknwon density estimation method chosen: '", density_method[1], "'")
  }
  density_params = c(list(dimred), density_params)
  densities = do.call(density_func, density_params)

  return(list(densities = densities, dimred = dimred, retained_counts = norm_counts))
}

select_markers = function(){
  #TODO: Bo.
}
