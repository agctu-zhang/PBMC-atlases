lapply(c("dplyr","Seurat"), library, character.only = T)
library(ggplot2)
library(patchwork)
library(stringr)
library(DoubletFinder)
library(cols4all)
library(harmony)
library(dplyr)
library(tidyr)

setwd('~/data/PBMC/')

plan('multisession', workers = 8)
plan()

species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')

sp.list <- list(
  readRDS('./rds/Catfish_raw.rds'),
  readRDS('./rds/Jacopever_raw.rds'),
  readRDS('./rds/Turtle_raw.rds'),
  readRDS('./rds/Chicken_raw.rds'),
  readRDS('./rds/Bat_raw.rds'),
  readRDS('./rds/Pig_raw.rds'),
  readRDS('./rds/Cattle_raw.rds'),
  readRDS('./rds/Mouse_raw.rds'),
  readRDS('./rds/Rat_raw.rds'),
  readRDS('./rds/Monkey_raw.rds'),
  readRDS('./rds/Chimpanzee_raw.rds'),
  readRDS('./rds/Human_raw.rds')
)


names(sp.list) <- species


for(i in names(sp.list)){
  DimPlot(sp.list[[i]],group.by = 'seurat_clusters',label = T,repel = T) +
    umap.theme +
    ggtitle('')
  filename = paste('./plot/cluster/seurat_clusters_',i,'.pdf',sep = '')
  ggsave(filename = filename)
}


for(i in names(sp.list)){
  DimPlot(sp.list[[i]],reduction = 'harmony',group.by = 'sample') +
    umap.theme +
    ggtitle('')
  filename = paste('./plot/quality/harmony_',i,'.pdf',sep = '')
  ggsave(filename = filename)
}

job::job({
  plan('multisession', workers = 8)
  for(i in names(sp.list)){
    Idents(sp.list[[i]]) <- 'seurat_clusters'
    marker <- FindAllMarkers(sp.list[[i]])
    filename = paste('./file/markers_seurat_clusters_',i,'.xlsx',sep = '')
    openxlsx::write.xlsx(marker,file=filename)
  }
},title = 'FindMarker')

for(i in names(sp.list)){
  filename <- paste0('./file/markers_seurat_clusters_',i,'.xlsx')
  markers <- openxlsx::read.xlsx(filename)
  if(i == 'Catfish'){
    catfish_anno <- read.csv('./ncbi/catfish_annotation.csv')
    catfish_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Catfish')
    markers$human <- catfish_ortho$Human.gene.name[match(markers$gene,catfish_ortho$Catfish.gene.name)]
    markers$description <- catfish_anno$product[match(markers$gene,catfish_anno$gene)]
    openxlsx::write.xlsx(markers,filename)
  }else if(i == 'Jacopever'){
    jacopever_anno <- read.csv('./ncbi/jacopever_annotation.csv')
    jacopever_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Jacopever')
    markers$human <- jacopever_ortho$Human.gene.name[match(markers$gene,jacopever_ortho$Jacopever.gene.name)]
    markers$description <- jacopever_anno$product[match(markers$gene,jacopever_anno$gene)]
    openxlsx::write.xlsx(markers,filename)
  }else{
    markers_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = i)
    markers$human <- markers_ortho$Human.gene.name[match(markers$gene,markers_ortho[[paste0(i,'.gene.name')]])]
    openxlsx::write.xlsx(markers,filename)
  }
}


for (i in names(sp.list)){
  print(i)
  print(table(sp.list[[i]]$multi_annotation))
}

for (i in names(sp.list)){
  print(i)
  print(ncol(sp.list[[i]]))
}

celltype <- sp.list[[1]]$multi_annotation %>% as.factor() %>% levels()

for (i in 2:length(sp.list)){
  celltypei <- sp.list[[i]]$multi_annotation %>% as.factor() %>% levels()
  celltype <- c(celltype,celltypei)
}

celltype <- unique(celltype)
celltype

for(i in names(sp.list)){
  sp.list[[i]] <- sp.list[[i]][,!sp.list[[i]]$multi_annotation %in% c("Doublets","Low quality cells",'Cycling cells','Erythrocytes')]
}

cell_num <- ncol(sp.list[[1]])
for (i in 2:length(sp.list)){
  cell_num <- cell_num + ncol(sp.list[[i]])
}
cell_num


for (i in names(sp.list)){
  print(i)
  print(ncol(sp.list[[i]]))
}


for(i in names(sp.list)){
  sp.list[[i]] <- sp.list[[i]] %>%
    SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = 'sample')
  
  n_pcs <- npcs(object = sp.list[[i]], reduction = 'harmony')
  
  sp.list[[i]] <- sp.list[[i]] %>%
    RunUMAP(reduction = "harmony", dims = 1:n_pcs)
}

for(i in names(sp.list)){
  DimPlot(sp.list[[i]],group.by = 'multi_annotation',cols = cell.colors1) +
    umap.theme +
    ggtitle('') +
    NoLegend()
  filename = paste('./plot/picture/dimplot_',i,'.pdf',sep = '')
  ggsave(filename = filename,width = 6,height = 6)
}

for(i in names(sp.list)){
  DimPlot(sp.list[[i]],group.by = 'multi_annotation',cols = cell.colors1,label = T,repel = T) +
    umap.theme +
    ggtitle('') +
    NoLegend()
  filename = paste('./plot/picture/dimplot_label_',i,'.pdf',sep = '')
  ggsave(filename = filename,width = 6,height = 6)
}

DimPlot(sp.list[[10]],group.by = 'multi_annotation',cols = cell.colors1,label = T,repel = T) +
  umap.theme +
  ggtitle('') +
  NoLegend()



cell.types <- c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','Plasma cells',
                'Platelets','Neutrophils','Low quality cells','Doublets','Prolymphocytes')


cell.data <- data.frame(cell.types = cell.types,
                        x = 1:length(cell.types),
                        y = rep(1, length(cell.types)))

cell.data$cell.types <- factor(cell.data$cell.types, levels = cell.types)

ggplot(cell.data, aes(x = x, y = y, fill = cell.types)) +
  geom_point(size = 5, shape = 21) +
  geom_text(aes(label = cell.types), vjust = 0.5, size = 3, angle = -90, hjust = 0) +
  scale_fill_manual(values = cell.colors1) +
  scale_colour_manual(values = cell.colors1) +
  theme_void()


ggsave('./plot/picture/all_legend.pdf')

job::job({
  plan('multisession', workers = 6)
  for(i in names(sp.list)){
    Idents(sp.list[[i]]) <- 'multi_annotation'
    marker <- FindAllMarkers(sp.list[[i]])
    filename = paste('./file/markers_multi_annotation_',i,'.xlsx',sep = '')
    openxlsx::write.xlsx(marker,file=filename)
  }
},title = 'FindMarker')

for(i in names(sp.list)){
  filename <- paste0('./file/markers_multi_annotation_',i,'.xlsx')
  markers <- openxlsx::read.xlsx(filename)
  if(i == 'Catfish'){
    catfish_anno <- read.csv('./ncbi/catfish_annotation.csv')
    catfish_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Catfish')
    markers$human <- catfish_ortho$Human.gene.name[match(markers$gene,catfish_ortho$Catfish.gene.name)]
    markers$description <- catfish_anno$product[match(markers$gene,catfish_anno$gene)]
    openxlsx::write.xlsx(markers,filename)
  }else if(i == 'Jacopever'){
    jacopever_anno <- read.csv('./ncbi/jacopever_annotation.csv')
    jacopever_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Jacopever')
    markers$human <- jacopever_ortho$Human.gene.name[match(markers$gene,jacopever_ortho$Jacopever.gene.name)]
    markers$description <- jacopever_anno$product[match(markers$gene,jacopever_anno$gene)]
    openxlsx::write.xlsx(markers,filename)
  }else{
    markers_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = i)
    markers$human <- markers_ortho$Human.gene.name[match(markers$gene,markers_ortho[[paste0(i,'.gene.name')]])]
    openxlsx::write.xlsx(markers,filename)
  }
}


################# dotplot ###################
cell.types <- c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','DCs','Platelets','Neutrophils','Prolymphocytes')

for (i in names(sp.list)){
  sp.list[[i]]$multi_annotation <- as.factor(sp.list[[i]]$multi_annotation)
  cell_order <- cell.types[cell.types %in% levels(sp.list[[i]]$multi_annotation)]
  sp.list[[i]]$multi_annotation <- factor(sp.list[[i]]$multi_annotation, levels = cell_order)
  Idents(sp.list[[i]]) <- 'multi_annotation'
}

gene.list <- list(
  'Catfish'    = c('cd79a','tcf7','csf1ra','flt3','gp1bb','ncf1','ikzf2'),
  'Jacopever'  = c('cd79a','tcf7','nkl.4','marco','flt3','gp1bb','ncf1'),
  'Turtle'     = c('JCHAIN','CD3E','TGFBI','MEIS1'),
  'Chicken'    = c('CD79B','TCF7','GNLY','CSF1R','NCF1C','GP1BB','IKZF2'),
  'Bat'        = c('CD79A','CD3E','GNLY','CSF1R','FCER1A','GP1BB'),
  'Pig'        = c('CD79A','CD3E','GNLY','CSF1R','FCER1A'),
  'Cattle'     = c('CD79A','CD3E','GNLY','CSF1R','FCER1A','FLT3','NCF1'),
  'Mouse'      = c('Cd79a','Cd3e','Ncr1','Csf1r','Siglech','S100a8'),
  'Rat'        = c('Cd79a','Cd3e','Ncr1','Csf1r','Clec9a','Siglech','Gp1bb','S100a9'),
  'Monkey'     = c('CD79A','CD3E','GNLY','CSF1R','FLT3','GP1BB'),
  'Chimpanzee' = c('CD79A','CD3E','GNLY','CSF1R','FCER1A'),
  'Human'      = c('CD79A','CD3E','GNLY','CD36','FCER1A','CLEC4C','TUBB1')
)

DotPlot(sp.list[[9]], features = gene.list[[9]]) +
  labs(x = '', y = '', title = names(sp.list[9]))

for(i in names(sp.list)){
  DotPlot(sp.list[[i]], features = gene.list[[i]]) +
    labs(x = '', y = '', title = i)
  
  filename = paste('./plot/picture/dotplot_',i,'.pdf',sep = '')
  ggsave(filename = filename)
}

table(sp.list[[1]]$cell_annotation)


# jjDotPlot(object = sp.list[[1]],
#           gene= gene.list[[1]],
#           id = 'multi_annotation',
#           dot.col = c('white','#023049'),
#           ytree = F,rescale = T)

# for(i in names(sp.list)){
#   jjDotPlot(object = sp.list[[i]],
#             gene= gene.list[[i]],
#             id = 'multi_annotation',dot.col = c('white','#023049'),ytree = F,rescale = T)
#   filename = paste('./plot/picture/dotplot_',i,'.pdf',sep = '')
#   ggsave(filename = filename)
# }
calculate_cell_stats <- function(seurat_obj) {
  cell_counts <- table(seurat_obj$multi_annotation)
  percentages <- prop.table(cell_counts) * 100
  return(list(counts = cell_counts, percentages = percentages))
}


species_names <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
                   'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')
all_stats <- lapply(sp.list, calculate_cell_stats)


all_cell_types <- unique(unlist(lapply(all_stats, function(x) names(x$counts))))
count_df <- data.frame(row.names = all_cell_types)
percent_df <- data.frame(row.names = all_cell_types)


for (i in seq_along(all_stats)) {
  species <- species_names[i]
  counts <- all_stats[[i]]$counts
  percentages <- all_stats[[i]]$percentages
  
  count_df[[species]] <- counts[match(rownames(count_df), names(counts))]
  percent_df[[species]] <- percentages[match(rownames(percent_df), names(percentages))]
}


count_df[is.na(count_df)] <- 0
percent_df[is.na(percent_df)] <- 0


print("Cell Counts:")
print(count_df)
print("\nCell Percentages:")
print(percent_df)

write.csv(count_df, "./file/cell_type_counts.csv")
write.csv(percent_df, "./file/cell_type_percentages.csv")

saveRDS(sp.list, './rds/sp.list.rds')

############### filter #################
celltype <- sp.list[[1]]$multi_annotation %>% as.factor() %>% levels()

for (i in 2:length(sp.list)){
  celltypei <- sp.list[[i]]$multi_annotation %>% as.factor() %>% levels()
  celltype <- c(celltype,celltypei)
}

celltype <- unique(celltype)
celltype


#### celltype number ####
cell_proportion <- table(sp.list[[1]]$multi_annotation) %>% as.data.frame()
colnames(cell_proportion) <- c('Cell_type', names(sp.list[1]))

for (i in 2:length(sp.list)){
  cell_proportioni <- table(sp.list[[i]]$multi_annotation) %>% as.data.frame()
  colnames(cell_proportioni) <- c('Cell_type', names(sp.list[i]))
  cell_proportion <- merge(cell_proportion, cell_proportioni, by = 'Cell_type', all = T)
}
rownames(cell_proportion) <- cell_proportion$Cell_type
cell_proportion <- cell_proportion[ ,-1]
cell_proportion[is.na(cell_proportion)] <- 0




Cellratio <- cell_proportion
Cellratio <- as.data.frame(Cellratio)
Cellratio$Cell_type <- rownames(Cellratio)
Cellratio <- tidyr::gather(Cellratio, 'Species', 'Proportion', -Cell_type)
Cellratio$Cell_type <- factor(Cellratio$Cell_type, levels = c('T cells','B cells','NK cells','Monocytes','pDCs','mDCs','DCs',
                                                              'Neutrophils','Platelets','Prolymphocytes'))
Cellratio$Species <- factor(Cellratio$Species,levels = rev(c('Catfish','Jacopever','Turtle','Chicken',
                                                         'Bat','Cattle','Pig','Mouse','Rat','Monkey','Chimpanzee','Human')))


ggplot(Cellratio) + 
  geom_bar(aes(x = Proportion, y = Species, fill = Cell_type),stat = "identity",width = 0.6,linewidth = 0,colour = '#222222')+ 
  theme_classic() +
  labs(y='Species', x = 'Cell Number', fill = 'Cell type')+
  scale_fill_manual(values = cell.colors1)+
  theme()
ggsave('./plot/picture/cell_proportion.pdf', width = 8)



meta_data <- list()
for (i in names(sp.list)){
  meta_data[[i]] <- sp.list[[i]]@meta.data
}

for (i in names(meta_data)){
  meta_data[[i]]$percentage <- NA
  for(j in 1:nrow(meta_data[[i]])){
    if (meta_data[[i]]$multi_annotation[j] %in% c('T cells','B cells','NK cells','pDCs','Prolymphocytes')){
      meta_data[[i]]$percentage[j] <- 'Lymphocytes'
    }else if (meta_data[[i]]$multi_annotation[j] %in% c('Monocytes','mDCs')){
      meta_data[[i]]$percentage[j] <- 'Myelocytes'
    }else {
      meta_data[[i]]$percentage[j] <- 'Others'
    }
  }
}

cell_proportion <- table(meta_data[[1]]$percentage) %>% as.data.frame()
colnames(cell_proportion) <- c('Cell_type', names(meta_data[1]))

for (i in 2:length(meta_data)){
  cell_proportioni <- table(meta_data[[i]]$percentage) %>% as.data.frame()
  colnames(cell_proportioni) <- c('Cell_type', names(meta_data[i]))
  cell_proportion <- merge(cell_proportion, cell_proportioni, by = 'Cell_type', all = T)
}
rownames(cell_proportion) <- cell_proportion$Cell_type
cell_proportion <- cell_proportion[ ,-1]
cell_proportion[is.na(cell_proportion)] <- 0
cell_proportion <- cell_proportion[-3, ]

Cellratio <- prop.table(as.matrix(cell_proportion), margin = 2)
Cellratio <- cell_proportion
Cellratio <- as.data.frame(Cellratio)
Cellratio$Cell_type <- rownames(Cellratio)
Cellratio <- tidyr::gather(Cellratio, 'Species', 'Proportion', -Cell_type)
Cellratio$Cell_type <- factor(Cellratio$Cell_type, levels = c('T cells','B cells','NK cells','Monocytes','pDCs','mDCs','DCs',
                                                              'Neutrophils','Platelets','Prolymphocytes'))
Cellratio$Species <- factor(Cellratio$Species,levels = rev(c('Catfish','Jacopever','Turtle','Chicken',
                                                             'Bat','Cattle','Pig','Mouse','Rat','Monkey','Chimpanzee','Human')))


ggplot(Cellratio) + 
  geom_bar(aes(x = Proportion, y = Species, fill = Cell_type),stat = "identity",width = 0.6,linewidth = 0,colour = '#222222')+ 
  theme_classic() +
  labs(y='Species', x = 'Cell Number', fill = 'Cell type')+
  scale_fill_manual(values = cell.colors1)+
  theme()
ggsave('./plot/picture/cell_proportion.pdf', width = 8)
