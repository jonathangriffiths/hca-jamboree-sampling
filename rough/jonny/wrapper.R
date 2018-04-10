source("~/Documents/hca-jamboree-sampling/rough/jonny/density_functions.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/dimred.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/gene_selection.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/normalisation.R")

wrapper = function(count_matrix){
  #normalise
  norm_counts = normalise_scran(count_matrix)
  #gene selection
  norm_counts = scran_select(norm_counts)
  #dimension reduction
  dimred = dimred_pca(norm_counts, dims = 20)
  #densities
  densities = get_cell_densities(dimred)
  #return
  return(densities)
}

test = wrapper(count_matrix)
