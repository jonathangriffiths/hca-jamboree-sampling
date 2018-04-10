density_wrapper = function(count_matrix,
                   qc_function = c("simple"),
                   min.umis = 1000, max.mito = 0.03,
                   normalisation_method = c("scran", "clr", "cpm", "none"),
                   gene_selection_method = "scran",
                   dimred_method = "pca",
                   density_k = 10){
  #qc
  if(qc_function[1] == "simple"){
    count_matrix = qc_simple(count_matrix, mito.frac = max.mito, min.umis = min.umis, gene_df = genes)
  } else {
    stop("No suitable QC method chosen.")
  }

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
    dimred = dimred_pca(norm_counts, dims = 20)
  } else{
    stop("No suitable dimension reduction method chosen.")
  }

  #densities
  densities = get_cell_densities(dimred, k = density_k)

  #return
  return(list(densities = densities, dimred = dimred, retained_counts = norm_counts))
}

select_markers = function(){
  #TODO: Bo.
}
