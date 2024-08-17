library(SCISSORS)
library(Seurat)

setwd('~/data/PBMC/')

print('load data')

sp.list <- readRDS('./rds/sp.list.quality_filter.rds')
MC.seu.list <- readRDS('./rds/supercell.seu.list.rds')

ComputeSilhouetteScores <- function(seurat.obj = NULL,
                                    dist.metric = "cosine",
                                    avg = TRUE,
                                    group.by = 'seurat_clusters',
                                    reduction = "pca") {
  # check inputs
  if (is.null(seurat.obj)) { stop("You didn't supply a Seurat object to ComputeSilhouetteScores().") }
  # prepare input matrix
  pca_df <- data.frame(Seurat::Embeddings(seurat.obj, reduction = reduction))
  pca_mat <- as.matrix(pca_df)
  # calculate distance matrix -- default is cosine dissimilarity
  if (dist.metric == "cosine") {
    pc_dists <- CosineDist(input.mat = pca_mat)
  } else {
    pc_dists <- stats::dist(x = pca_mat, method = dist.metric)
  }
  clust_list <- as.integer(unlist(seurat.obj[[group.by]])) - 1L
  res <- cluster::silhouette(dist = pc_dists, x = clust_list)
  # prepare results & return
  if (avg) {
    avg_widths <- summary(res)$clus.avg.widths
    avg_widths <- unlist(avg_widths)
    val <- avg_widths
  } else {
    val <- data.frame(Cluster = as.factor(res[, 1]), Score = res[, 3])
  }
  return(val)
}

for(i in names(MC.seu.list)){
  sp.list[[i]]$multi_annotation <- as.factor(sp.list[[i]]$multi_annotation)
  MC.seu.list[[i]]$cell_line <- as.factor(MC.seu.list[[i]]$cell_line)
}


ComputeSilhouetteScores(sp.list[[1]], group.by = 'multi_annotation', reduction = 'harmony') %>% as.data.frame()

ComputeSilhouetteScores(MC.seu.list[[1]], group.by = 'cell_line', reduction = 'pca') %>% as.data.frame()

print('silhouette of single cell')

silhouette1 <- ComputeSilhouetteScores(sp.list[[1]], group.by = 'multi_annotation', reduction = 'harmony') %>% as.data.frame()
colnames(silhouette1) <- names(sp.list[1])
silhouette1$cluster <- rownames(silhouette1)
for(i in 2:length(sp.list)){
  silhouettei <- ComputeSilhouetteScores(sp.list[[i]], group.by = 'multi_annotation', reduction = 'harmony') %>% as.data.frame()
  colnames(silhouettei) <- names(sp.list[i])
  silhouettei$cluster <- rownames(silhouettei)
  silhouette1 <- merge(silhouette1, silhouettei, by = 'cluster', all = T)
}

print('silhouette of meta cell')

silhouette2 <- ComputeSilhouetteScores(MC.seu.list[[1]], group.by = 'cell_line', reduction = 'pca') %>% as.data.frame()
colnames(silhouette2) <- names(MC.seu.list[1])
silhouette2$cluster <- rownames(silhouette2)
for(i in 2:length(MC.seu.list)){
  silhouettei <- ComputeSilhouetteScores(MC.seu.list[[i]], group.by = 'cell_line', reduction = 'pca') %>% as.data.frame()
  colnames(silhouettei) <- names(MC.seu.list[i])
  silhouettei$cluster <- rownames(silhouettei)
  silhouette2 <- merge(silhouette2, silhouettei, by = 'cluster', all = T)
}

silhouette1 <- silhouette1[,-1]
silhouette1 <- tidyr::gather(silhouette1, 'Species', 'Score')
silhouette1$Cell_type <- 'Single_cell'

silhouette2 <- silhouette2[,-1]
silhouette2 <- tidyr::gather(silhouette2, 'Species', 'Score')
silhouette2$Cell_type <- 'Meta_cell'

sil_score <- rbind(silhouette1,silhouette2)
sil_score$Species <- factor(sil_score$Species, levels = c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
                                                          'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human'))
print('plot of silhouette score')

library(ggplot2)

ggplot(sil_score, aes(x = Species, y = Score, fill = Cell_type)) +
  geom_boxplot() +
  labs(x = "", y = "Silhouette Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./plot/supercell/silhouette_score.pdf')

print('statistic of gene number')

gene_number1 <- sp.list[[1]][['nFeature_SCT']]
gene_number1$Species <- names(sp.list[1])
gene_number1$Cell_type <- 'Single_cell'


for(i in 2:length(sp.list)){
  gene_numberi <- sp.list[[i]][['nFeature_SCT']]
  gene_numberi$Species <- names(sp.list[i])
  gene_numberi$Cell_type <- 'Single_cell'
  gene_number1 <- rbind(gene_number1, gene_numberi)
}


gene_number2 <- MC.seu.list[[1]][['nFeature_RNA']]
gene_number2$Species <- names(MC.seu.list[1])
gene_number2$Cell_type <- 'Meta_cell'


for(i in 2:length(MC.seu.list)){
  gene_numberi <- MC.seu.list[[i]][['nFeature_RNA']]
  gene_numberi$Species <- names(MC.seu.list[i])
  gene_numberi$Cell_type <- 'Meta_cell'
  
  gene_number2 <- rbind(gene_number2, gene_numberi)
}

colnames(gene_number1)[1] <- 'Gene_number'
colnames(gene_number2)[1] <- 'Gene_number'
gene_number <- rbind(gene_number1, gene_number2)

gene_number$Species <- factor(gene_number$Species, levels = c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
                                                          'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human'))
print('plot of gene number')

ggplot(gene_number, aes(x = Species, y = Gene_number, fill = Cell_type)) +
  geom_boxplot() +
  labs(x = "", y = "Number of genes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./plot/supercell/number_of_genes.pdf')
