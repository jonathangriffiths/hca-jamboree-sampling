density_wrapper = function(count_matrix,
                   qc_method = "simple",
                   qc_params = list(),
                   normalisation_method = c("scran", "clr", "cpm", "none"),
                   normalisation_params = list(),
                   gene_selection_method = "scran",
                   gene_selection_params = list(),
                   dimred_method = "pca",
                   dimred_params = list(),
                   density_k = 10){

  #QC
  qc_funcs = list()
  qc_funcs$simple = qc_simple

  qc_func = qc_funcs[[qc_method[1]]]
  if (is.null(qc_func)) {
      stop("Unknown QC method chosen: '", qc_method[1], "'")
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
      stop("Unknown normalisation method chosen: '", normalisation_method[1], "'")
  }
  normalisation_params = c(list(count_matrix), normalisation_params)
  norm_counts = do.call(norm_func, normalisation_params)
  
  all_gene_counts = norm_counts

  #gene selection
  gene_selection_funcs = list()
  gene_selection_funcs$scran = scran_select

  gene_selection_func = gene_selection_funcs[[gene_selection_method[1]]]
  if (is.null(gene_selection_func)) {
      stop("Unknown gene selection method chosen: '", gene_selection_method[1], "'")
  }
  gene_selection_params = c(list(norm_counts), gene_selection_params)
  norm_counts = do.call(gene_selection_func, gene_selection_params)

  #dimension reduction
  dimred_funcs = list()
  dimred_funcs$pca = dimred_pca

  dimred_func = dimred_funcs[[dimred_method[1]]]
  if (is.null(dimred_func)) {
      stop("Unknown dimension reduction method chosen: '", dimred_method[1], "'")
  }
  dimred_params = c(list(norm_counts), dimred_params)
  dimred = do.call(dimred_func, dimred_params)

  #densities
  densities = get_cell_densities(dimred, k = density_k)

  return(list(densities = densities, dimred = dimred, all_expr = all_gene_counts))
}

select_markers = function(density_wrapper_out,
                          gene_selection_method = c("quantile", "mixture_model"),
                          cell_selection_params = list(),
                          gate_gene_method = c("lasso")
                          ){
  #cell selection
  cell_selection_funcs = list()
  cell_selection_funcs$quantile = select_cells_quantile
  cell_selection_funcs$mixmod = select_cells_mixmod
  
  cell_selection_func = cell_selection_funcs[[cell_selection_method[1]]]
  if (is.null(cell_selection_func)) {
    stop("Unknown dimension reduction method chosen: '", cell_selection_method[1], "'")
  }
  cell_selection_params = c(list(density_wrapper_out$densities), cell_selection_params)
  cell_selection = do.call(cell_selection_func, cell_selection_params)
  
  #identify key genes
  gate_gene_funcs = list()
  gate_gene_funcs$lasso = find_genes_lasso
  
  gate_gene_func = gate_gene_funcs[[gate_gene_method[1]]]
  if (is.null(gate_gene_func)) {
    stop("Unknown dimension reduction method chosen: '", gate_gene_method[1], "'")
  }
  gate_gene_params = c(list(cell_selection$hits, density_wrapper_out$all_expr), gate_gene_params)
  gate_gene = do.call(gate_gene_func, gate_gene_params)
}
