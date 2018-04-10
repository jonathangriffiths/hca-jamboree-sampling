library(Matrix)
library(Rtsne)
library(viridis)
library(cowplot)
devtools::load_all('team13')
counts = readMM("data/ica_bone_marrow_ref_set.mtx")
genes = read.table("data/ica_bone_marrow_ref_set.genes.tsv")
cells = read.table("data/ica_bone_marrow_ref_set.barcodes.tsv")
rownames(counts) = genes[,2]
colnames(counts) = cells[,1]
counts = as(counts, "dgCMatrix")




test = density_wrapper(counts, normalisation_method = "scran")

tsne = Rtsne(test$dimred, pca = FALSE, check_duplicates = FALSE)

p1 = ggplot(as.data.frame(tsne$Y), aes(x = V1, y = V2, col = test$densities)) +
  geom_point() +
  scale_color_viridis()

p2 = ggplot(as.data.frame(tsne$Y), aes(x = V1, y = V2, col = test$retained_counts["ENSG00000198851",])) +
  geom_point() +
  scale_color_viridis(name=  "CD3E")

p2 = ggplot(as.data.frame(tsne$Y), aes(x = V1, y = V2, col = log2(counts["ENSG00000198851", colnames(counts) %in% colnames(test$retained_counts)]+1))) +
  geom_point() +
  scale_color_viridis(name=  "CD3E")

plot_grid(p1, p2)

hits = select_cells(test, quantile = 0.8)
