source("~/Documents/hca-jamboree-sampling/rough/jonny/density_functions.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/dimred.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/gene_selection.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/normalisation.R")
source("~/Documents/hca-jamboree-sampling/rough/jonny/qc.R")


counts = readMM("~/Documents/hca-jamboree-sampling/data/ica_bone_marrow_ref_set.mtx")
genes = read.table("~/Documents/hca-jamboree-sampling/data/ica_bone_marrow_ref_set.genes.tsv")
cells = read.table("~/Documents/hca-jamboree-sampling/data/ica_bone_marrow_ref_set.barcodes.tsv")
rownames(counts) = genes[,2]
rownames(cells) = cells[,1]
counts = as(counts, "dgCMatrix")

wrapper = function(count_matrix){
  #qc
  count_matrix = qc_simple(count_matrix)
  #normalise
  norm_counts = normalise_scran(count_matrix)
  #gene selection
  norm_counts = scran_select(norm_counts)
  #dimension reduction
  dimred = dimred_pca(norm_counts, dims = 20)
  #densities
  densities = get_cell_densities(dimred)
  #return
  return(list(densities = densities, dimred = dimred, retained_counts = norm_counts))
}

test = wrapper(counts)

tsne = Rtsne(test$dimred, pca = FALSE)

p1 = ggplot(as.data.frame(tsne$Y), aes(x = V1, y = V2, col = test$densities)) +
  geom_point() +
  scale_color_viridis()

p2 = ggplot(as.data.frame(tsne$Y), aes(x = V1, y = V2, col = test$retained_counts["ENSG00000198851",])) +
  geom_point() +
  scale_color_viridis(name=  "CD3E")

plot_grid(p1, p2)


