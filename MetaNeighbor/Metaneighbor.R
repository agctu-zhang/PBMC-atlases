library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)
library(stringr)
library(ggplot2)
library(gghalves)
library(ggsignif)
library(ggpubr)
library(Seurat)
library(colorRamp2)
#### Unsupervised metaneighbor
setwd('~/data/PBMC/')
getwd()

MC.seu.list <- readRDS('./rds/supercell.seu.list.rds')
############ filter cells ############
mn.list <- list()
for (i in names(MC.seu.list)){
  mn.list[[i]] <- subset(MC.seu.list[[i]], subset = purity ==1)
}


cellname <- list()
for (i in names(mn.list)){
  cellname[[i]] <- mn.list[[i]]@meta.data
  cellname[[i]] <- cellname[[i]][ ,c('size','cell_line')]
  cellname[[i]]$cell_line <- gsub(' ', '_', cellname[[i]]$cell_line)
  cellname[[i]]$cell_name <- paste(cellname[[i]]$cell_line, names(mn.list[i]), rownames(cellname[[i]]), sep = '_')
}

for (i in names(mn.list)){
  mn.list[[i]] <- GetAssayData(mn.list[[i]],slot = 'data') %>% as.matrix()
}

for (i in names(mn.list)){
  colnames(mn.list[[i]]) <- cellname[[i]]$cell_name
}

##############  gene convert file ######
feature_all <- list()
for (i in species){
  feature_all[[i]] = rownames(mn.list[[i]]) %>% as.data.frame()
  colnames(feature_all[[i]])[1] = paste0(i, '.gene.name')
}

names(feature_all)

for (i in names(feature_all)) {
  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_human.xlsx", sheet = i)
  column_name <- paste0(i, '.gene.name')
  feature_all[[i]][["Human.gene.name"]] <- all_ortholog$Human.gene.name[match(feature_all[[i]][[column_name]], all_ortholog[[column_name]])]
}

saveRDS(feature_all, './rds/orthologs_to_human.rds')
# feature_all <- readRDS('./rds/orthologs_to_mouse.rds')
######## gene name change #######
for (i in names(mn.list)){
  feature_all[[i]][is.na(feature_all[[i]])] <- "Unknown"  
  rownames(mn.list[[i]]) <- feature_all[[i]]$Human.gene.name
  idx <- which(rownames(mn.list[[i]]) == 'Unknown')
  mn.list[[i]] <- mn.list[[i]][-idx, ]
}

for (i in names(mn.list)){
  mn.list[[i]] <- mn.list[[i]] %>% as.data.frame()
  mn.list[[i]]$gene <- rownames(mn.list[[i]])
}

all_cell <- mn.list[[1]]
all_cell$gene <- rownames(all_cell)
for(i in 2:length(mn.list)){
  all_celli <- mn.list[[i]]
  all_celli$gene <- rownames(all_celli)
  all_cell <- merge(all_cell, all_celli, by = 'gene', all = T)
}

all_cell[4,4]
all_cell[is.na(all_cell)] <- 0

rownames(all_cell) <- all_cell$gene
all_cell[3,3]
all_cell <- all_cell[,-1] # remove gene column
all_cell <- all_cell %>% as.matrix()
all_cell[3,3]
saveRDS(all_cell,'./rds/metaneighbor_cell.rds')

# all_cell <- readRDS('./rds/metaneighbor_cell.rds')

source("/home/data/t090402/download/MetaNeighbor/2017-08-28-runMN-US.R")
load("/home/data/t090402/download/MetaNeighbor/MetaNeighbor_US_data.Rdata")
var.genes=get_variable_genes(data, pheno)
celltype.NV=run_MetaNeighbor_US(var.genes, data, celltypes, pheno)
get_top_hits(celltype.NV, pheno, threshold=0.9, filename="filename.txt")
cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.6,cexCol=0.6)

############## all cells ############
species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')

all_name <- c('tfd','sslg','pss','gga','ray','ssc',
              'bta','mmu','rno','mcc','ptr','hsa')

all_pheno <- colnames(all_cell) %>% as.data.frame()
colnames(all_pheno) <- 'Sample_ID'
matched_species <- sapply(species, function(x) ifelse(grepl(x, all_pheno$Sample_ID), x, ""))
matched_species <- matched_species[matched_species != ""]
all_pheno$Study_ID <- matched_species
all_pheno$Celltype <- gsub("_\\d+$", "", all_pheno$Sample_ID)

original_species <- sub(".*?_", "", all_pheno$Study_ID)
levels(as.factor(original_species))

for (i in seq_along(species)) {
  all_pheno$Study_ID <- gsub(species[i], all_name[i], all_pheno$Study_ID)
  all_pheno$Celltype<- gsub(paste0('_',species[i]), paste0('.',all_name[i]), all_pheno$Celltype)
}

all_celltypes <- levels(as.factor(all_pheno$Celltype))

var.genes=get_variable_genes(all_cell, all_pheno)
celltype.NV=run_MetaNeighbor_US(var.genes, all_cell, all_celltypes, all_pheno)
get_top_hits(celltype.NV, all_pheno, threshold=0.85, filename="/home/data/t090402/data/PBMC/plot/metaneighbor/filename.txt")
cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
pdf('./plot/metaneighbor/heatmap.2.pdf', width = 10, height = 10)
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.6,cexCol=0.6)
dev.off()

# pdf('./plot/metaneighbor/heatmap.pdf', width = 15, height = 15)
# Heatmap(celltype.NV,
#         row_dend_width = unit(20, "mm"),
#         column_dend_height = unit(20, "mm"),
#         clustering_method_rows = 'centroid',
#         clustering_method_columns = 'centroid',
#         show_heatmap_legend = F,
#         col = colorRamp2(c(0,1), c('white', "red")))
# dev.off()


celltype.NV <- as.matrix(celltype.NV)
all_name <- c('tfd','sslg','pss','gga','ray','ssc',
              'bta','mmu','rno','mcc','ptr','hsa')
species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')
cell_type <- sub("\\..*$", "", rownames(celltype.NV)) %>% stringr::str_replace_all("_", " ")
species_abb <- sub(".*\\.", "", colnames(celltype.NV))
species_abb <- factor(species_abb, levels = all_name)
# c4a_gui()

species_col <- cols4all::c4a("renoir",12)
names(species_col) <- species

species_abb <- factor(species[match(species_abb, all_name)], levels = species)


hb = HeatmapAnnotation(
  Species = species_abb,
  Celltypes = cell_type,
  col = list(
    Species = species_col,
    Celltypes = cell.colors1
  )
)

hr = rowAnnotation(
  Species = species_abb,
  Celltypes = cell_type,
  col = list(
    Species = species_col,
    Celltypes = cell.colors1
  ),
  show_legend = FALSE, 
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

pdf('./plot/metaneighbor/heatmap_legend.pdf', width = 8.5, height = 8)
Heatmap(celltype.NV,
        show_row_names = F,
        show_column_names = F,
        row_dend_width = unit(20, "mm"),
        column_dend_height = unit(20, "mm"),
        row_names_max_width = unit(10, "cm"),
        column_names_max_height = unit(10, "cm"),
        show_heatmap_legend = T,
        bottom_annotation = hb,
        right_annotation = hr,
        width = unit(10, "cm"), height = unit(10, "cm"),
        heatmap_legend_param = list(title = 'AUROC score'),
        col = colorRamp2(c(0,1), c('white', "red")))
dev.off()

write.csv(as.data.frame(celltype.NV),'/home/data/t090402/data/PBMC/plot/metaneighbor/celltype.NV.csv')
celltype.NV <- read.csv("/home/data/t090402/data/PBMC/plot/metaneighbor/celltype.NV.csv", row.names = 1)


sp.list <- readRDS('./rds/sp.list.rds')
############ filter cells ############
cellname <- list()
for (i in names(sp.list)){
  cellname[[i]] <- sp.list[[i]]@meta.data
  cellname[[i]] <- cellname[[i]][ ,c('orig.ident','multi_annotation')]
  cellname[[i]]$multi_annotation <- gsub(' ', '_', cellname[[i]]$multi_annotation)
  cellname[[i]]$cell_name <- paste(cellname[[i]]$multi_annotation, names(sp.list[i]), rownames(cellname[[i]]), sep = '_')
}

for (i in names(sp.list)){
  sp.list[[i]] <- GetAssayData(sp.list[[i]],slot = 'data') %>% as.matrix()
}

for (i in names(sp.list)){
  colnames(sp.list[[i]]) <- cellname[[i]]$cell_name
}

##############  gene convert file ######
feature_all <- list()
for (i in species){
  feature_all[[i]] = rownames(sp.list[[i]]) %>% as.data.frame()
  colnames(feature_all[[i]])[1] = paste0(i, '.gene.name')
}

names(feature_all)

for (i in names(feature_all)) {
  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_human.xlsx", sheet = i)
  column_name <- paste0(i, '.gene.name')
  feature_all[[i]][["Human.gene.name"]] <- all_ortholog$Human.gene.name[match(feature_all[[i]][[column_name]], all_ortholog[[column_name]])]
}

saveRDS(feature_all, './rds/all_cell_orthologs_to_human.rds')
# feature_all <- readRDS('./rds/orthologs_to_mouse.rds')
######## gene name change #######
for (i in names(sp.list)){
  feature_all[[i]][is.na(feature_all[[i]])] <- "Unknown" 
  rownames(sp.list[[i]]) <- feature_all[[i]]$Human.gene.name
  idx <- which(rownames(sp.list[[i]]) == 'Unknown')
  sp.list[[i]] <- sp.list[[i]][-idx, ]
}

for (i in names(sp.list)){
  sp.list[[i]] <- sp.list[[i]] %>% as.data.frame()
  sp.list[[i]]$gene <- rownames(sp.list[[i]])
}

all_cell <- sp.list[[1]]
all_cell$gene <- rownames(all_cell)
for(i in 2:length(sp.list)){
  all_celli <- sp.list[[i]]
  all_celli$gene <- rownames(all_celli)
  all_cell <- merge(all_cell, all_celli, by = 'gene', all = T)
}

all_cell[3,3]

rownames(all_cell) <- all_cell$gene
all_cell[3,3]
all_cell <- all_cell[,-1] # remove gene column
all_cell <- all_cell %>% as.matrix()

all_cell[is.na(all_cell)] <- 0
all_cell[3,3]

saveRDS(all_cell,'./rds/metaneighbor_all_cell.rds')

# all_cell <- readRDS('./rds/metaneighbor_cell.rds')

species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')

all_name <- c('tfd','sslg','pss','gga','ray','ssc',
              'bta','mmu','rno','mcc','ptr','hsa')

all_pheno <- colnames(all_cell) %>% as.data.frame()
colnames(all_pheno) <- 'Sample_ID'
matched_species <- sapply(species, function(x) ifelse(grepl(x, all_pheno$Sample_ID), x, ""))
matched_species <- matched_species[matched_species != ""]
all_pheno$Study_ID <- matched_species
all_pheno <- all_pheno %>%
  mutate(Celltype = str_extract(Sample_ID, "^[^_]+"))

original_species <- sub(".*?_", "", all_pheno$Study_ID)
levels(as.factor(original_species))
for (i in seq_along(species)) {
  all_pheno$Study_ID <- gsub(species[i], all_name[i], all_pheno$Study_ID)
  all_pheno$Celltype<- gsub(paste0('_',species[i]), paste0('.',all_name[i]), all_pheno$Celltype)
}

all_celltypes <- levels(as.factor(all_pheno$Celltype))

var.genes=get_variable_genes(all_cell, all_pheno)
celltype.NV=run_MetaNeighbor_US(var.genes, all_cell, all_celltypes, all_pheno)
saveRDS(celltype.NV, './plot/metaneighbor/all_cell.celltype.NV.rds')
get_top_hits(celltype.NV, all_pheno, threshold=0.85, filename="/home/data/t090402/data/PBMC/plot/metaneighbor/filename_all_cell.txt")
cols=rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(100))
breaks=seq(0,1,length=101)
# pdf('./plot/metaneighbor/heatmap_cell.2.pdf', width = 10, height = 10)
heatmap.2(celltype.NV,trace="none",density.info="none",col=cols,breaks=breaks,cexRow=0.6,cexCol=0.6)
# dev.off()

pdf('./plot/metaneighbor/heatmap_cell.pdf', width = 15, height = 15)
Heatmap(celltype.NV,
        row_dend_width = unit(20, "mm"),
        column_dend_height = unit(20, "mm"),
        clustering_method_rows = 'centroid',
        clustering_method_columns = 'centroid',
        show_heatmap_legend = F,
        col = colorRamp2(c(0,1), c('white', "red")))
dev.off()
pdf('./plot/metaneighbor/heatmap_cell_legend.pdf', width = 15, height = 15)
Heatmap(celltype.NV,
        row_dend_width = unit(40, "mm"),
        column_dend_height = unit(40, "mm"),
        row_names_max_width = unit(10, "cm"),
        column_names_max_height = unit(10, "cm"),
        show_heatmap_legend = T,
        heatmap_legend_param = list(title = 'AUROC score'),
        col = colorRamp2(c(0,1), c('white', "red")))
dev.off()

write.csv(as.data.frame(celltype.NV),'/home/data/t090402/data/PBMC/plot/metaneighbor/celltype_cell.NV.csv')


#BiocManager::install('factoextra') 
library(factoextra) 
library(spatstat.geom)
sampleDists <- dist(celltype.NV)

res1 <- hcut(sampleDists, k = 1, stand = FALSE,hc_method ="average" ) 
res2 <- hclust(sampleDists)
fviz_dend(res2,
          rect_fill = T,
          k_colors = 'black',
          cex = 0.8,
          color_labels_by_k=T,
          horiz=T)
ggsave('./plot/metaneighbor/cluster_celltype.pdf', height = 15, width = 5)

###### similar cell type pairs among major categories #######
celltype.NV <- read.csv('/home/data/t090402/data/PBMC/plot/metaneighbor/celltype.NV.csv')

head(celltype.NV)

melted_data <- reshape2::melt(celltype.NV, id.vars = "X")
melted_data$variable <- as.character(melted_data$variable)
head(melted_data)

cell_types <- sub("\\..*", "", rownames(celltype.NV))

melted_data$interaction <- paste0(sub("\\..*", "", melted_data$X),'-',sub("\\..*", "", melted_data$variable))
melted_data$cross <- ifelse(sub("\\..*", "", melted_data$X) == sub("\\..*", "", melted_data$variable), 'Same-category', 'Cross-category')
melted_data$species <- ifelse(substring(melted_data$X, nchar(melted_data$X) - 2) == 
                                substring(melted_data$variable, nchar(melted_data$variable) - 2), "Same-specie", "Cross-specie")
melted_data$category <- ifelse(melted_data$species == "Same-specie", "Inner-specie",
                              ifelse(melted_data$cross == 'Same-category' & melted_data$species == "Cross-specie", 'Inner-category', 'Filter'))

head(melted_data)

openxlsx::write.xlsx(melted_data,'./plot/metaneighbor/melted_data.xlsx')

ggplot(melted_data[melted_data$X != melted_data$variable, ],
       aes(x = interaction, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "", y = "AUROC score") +
  theme_classic() +
  # scale_fill_manual(values = fill_color) +
  NoLegend() +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

########### same or different category #########
ggviolin(melted_data, x="cross", y="value", 
         color = "cross",
         fill="cross",
         #palette =c("#B3CDE3","#DECBE4"),
         add = "boxplot",
         add.params = list(color="white"),
         xlab = F
         )+
  NoLegend()+
  labs(x = "", y = "AUROC score") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_signif(comparisons = list(levels(as.factor(melted_data$cross)))
              ,test=t.test,map_signif_level = F)

ggsave('./plot/metaneighbor/category_same_cross.pdf', width = 4.5, height = 4.5)


########### same or different category #########
ggviolin(melted_data, x="species", y="value", 
         color = "species",
         fill="species",
         #palette =c("#B3CDE3","#DECBE4"),
         add = "boxplot",
         add.params = list(color="white"),
         xlab = F
)+
  NoLegend()+
  labs(x = "", y = "AUROC score") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_signif(comparisons = list(levels(as.factor(melted_data$species)))
              ,test=t.test,map_signif_level = F)

ggsave('./plot/metaneighbor/species_same_cross.pdf', width = 4.5, height = 4.5)


########### same or different category #########
ggviolin(melted_data[melted_data$category != 'Filter',], x="category", y="value", 
         color = "category",
         fill="category",
         #palette =c("#B3CDE3","#DECBE4"),
         add = "boxplot",
         add.params = list(color="white"),
         xlab = F
)+
  NoLegend()+
  labs(x = "", y = "AUROC score") +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  geom_signif(comparisons = list(levels(as.factor(melted_data[melted_data$category != 'Filter',]$category)))
              ,test=t.test,map_signif_level = F)

ggsave('./plot/metaneighbor/inner_specie_category.pdf', width = 4.5, height = 4.5)

