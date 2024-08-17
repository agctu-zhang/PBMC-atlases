library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(webshot)
library(networkD3)
library(scFunctions)
library(ComplexHeatmap)
library(BiocParallel)
library(cols4all)
library(ggpubr)
library(cowplot)
library(ggalluvial)

setwd('~/data/PBMC/')
# dir.create('./SCENIC')

# scenicOptions <- initializeScenic(org="hgnc",#mouse填'mgi', human填'hgnc',fly填'dmel') 
#                                   dbDir="./../../download/database/SCENIC/",
#                                   dbs = list('500bp' = 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather',
#                                              '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'),
                                  # nCores=8)#这里可以设置并行计算


############# export expression matrix and metadata #######
species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')
sp.list <- readRDS('./rds/sp.list.quality_filter.rds')

# exprMat <- as.matrix(MC.seu.list[['Human']]@assays$RNA@data)
# dim(exprMat)
# cellInfo <- MC.seu.list[['Human']]@meta.data[,c(5,2,3)]
# colnames(cellInfo) <- c('CellType', 'nGene', 'nUMI')
# rownames(cellInfo) <- paste(names(MC.seu.list['Human']), rownames(cellInfo), sep = '_')
# colnames(exprMat) <- paste(names(MC.seu.list['Human']), colnames(exprMat), sep = '_')
# exprMat[1:4,1:4]
# 
# write.csv(t(exprMat), file = './SCENIC/human.csv')
# write.table(cellInfo,'./SCENIC/human_metadata.xls',sep='\t',quote=F)
# 
# exprMat <- read.csv('./SCENIC/human.csv')
# 
# exprMat <- as.matrix(MC.seu.list[['Mouse']]@assays$RNA@data)
# dim(exprMat)
# cellInfo <- MC.seu.list[['Mouse']]@meta.data[,c(5,2,3)]
# colnames(cellInfo) <- c('CellType', 'nGene', 'nUMI')
# rownames(cellInfo) <- paste(names(MC.seu.list['Mouse']), rownames(cellInfo), sep = '_')
# colnames(exprMat) <- paste(names(MC.seu.list['Mouse']), colnames(exprMat), sep = '_')
# exprMat[1:4,1:4]
# 
# write.csv(t(exprMat), file = './SCENIC/mouse.csv')
# write.table(cellInfo,'./SCENIC/mouse_metadata.xls',sep='\t',quote=F)
# 
# exprMat <- read.csv('./SCENIC/mouse.csv')
# 
# ########## plot ################
# rss <- readRDS('./SCENIC/human/human_rss.rds')
# regulon <- load('./SCENIC/human/human_regulon_RSS.Rdata')
# loom <- open_loom('./SCENIC/human/human.aucell.loom')
# 
# Heatmap(rss)

##############  #############
exprMat.list <- list()

for (i in names(sp.list)){
  exprMat.list[[i]] <- sp.list[[i]][,sp.list[[i]]$multi_annotation %in% c('Monocytes','B cells','T cells','NK cells','pDCs','mDCs','DCs')]
  
  exprMat.list[[i]] <- GetAssayData(exprMat.list[[i]], assay = 'RNA', slot = 'data')
}


feature_all <- list()
for (i in species){
  feature_all[[i]] = rownames(exprMat.list[[i]]) %>% as.data.frame()
  colnames(feature_all[[i]])[1] = i
}


names(feature_all)

# 同源转换
for (i in names(feature_all)) {
  # 读取对应的工作表数据
  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_human.xlsx", sheet = i)
  
  
  # 匹配并赋值给相应的列
  feature_all[[i]]$Human.gene.name <- all_ortholog$Human.gene.name[match(feature_all[[i]][[i]], all_ortholog[[paste0(i,'.gene.name')]])]
}


#######filter matrix
for (i in names(exprMat.list)){
  feature_all[[i]][is.na(feature_all[[i]])] <- "Unknown"  # 将没有匹配到的行名设置为 "Unknown"
  rownames(exprMat.list[[i]]) <- feature_all[[i]]$Human.gene.name
  # 找到行名为'Unknown'的行的下标
  idx <- which(rownames(exprMat.list[[i]]) == 'Unknown')
  
  # 从data.in中删除行名为'Unknown'的行
  exprMat.list[[i]] <- exprMat.list[[i]][-idx, ]
}

for (i in names(exprMat.list)){
  filename <- paste0('./SCENIC/others/',i,'.csv')
  colnames(exprMat.list[[i]]) <- paste(i, colnames(exprMat.list[[i]]), sep = '_')
  exprMat.list[[i]] <- as.matrix(exprMat.list[[i]])
  exprMat.list[[i]] <- t(exprMat.list[[i]])
  write.csv(exprMat.list[[i]],filename)
}

for (i in names(exprMat.list)){
  filename <- paste0('./SCENIC/others/',i,'_metadata.xls')
  cellInfo <- sp.list[[i]][,sp.list[[i]]$multi_annotation %in% 
                             c('Monocytes','B cells','T cells','NK cells','pDCs','mDCs','DCs')]@meta.data[,c('multi_annotation','nCount_RNA','nFeature_RNA')]
  colnames(cellInfo) <- c('CellType', 'nGene', 'nUMI')
  rownames(cellInfo) <- paste(i, rownames(cellInfo), sep = '_')
  message(head(rownames(cellInfo[[i]])))
  
  write.table(cellInfo,filename,sep='\t',quote=F)
}


############# HOW MANY ORTHOLUES OF TF #########
regulon.list <- list()
for (i in species){
  filename <- paste0('./SCENIC/others/',i,'.AUCell.loom')
  
  loom <- open_loom(filename)
  
  regulon.list[[i]] <- get_regulons(loom, column.attr.name="Regulons")
  
  close_loom(loom)
  
  regulon.list[[i]] <- rownames(regulon.list[[i]])
}

regulon_num <- data.frame(specie = species,
                          number = numeric(12),
                          orthologue = numeric(12),
                          row.names = species)

for (i in names(regulon.list)){
  regulon_num[i,"number"] <- length(regulon.list[[i]])
  regulon_num[i,"orthologue"] <- length(intersect(regulon.list[[i]], regulon.list[[12]]))
}

regulon_num$specie <- factor(regulon_num$specie, levels = species)
head(regulon_num)
regulon_num <-  gather(regulon_num, key = "variable", value = "value", -specie)

ggplot(regulon_num, aes(x = specie, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Number of TF", fill = '') +  # Set axis labels
  theme_classic() +
  scale_fill_manual(
    values = c("number" = "grey", "orthologue" = "darkblue"),
    labels = c("number" = "Identified by SCENIC", "orthologue" = "Intersect with human")
  ) +
  geom_text(aes(label = value), vjust = -0.3, size = 3, 
            position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5)) +
  theme(text = element_text(size = 12),
        axis.text = element_text(colour = "black")) +
  theme(legend.position = "top", legend.direction = "horizontal")

ggsave('./plot/SCENIC/TF_number.pdf', width = 6, height = 4)

regulon_ratio <- data.frame(specie = species,
                            ratio = (regulon_num[regulon_num$variable == "orthologue",]$value/regulon_num[regulon_num$variable == "number",]$value)*100)
regulon_ratio$specie <- factor(regulon_ratio$specie, levels = species)

ggplot(regulon_ratio, aes(x = specie, y = ratio)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Intersect with human %", fill = '') +  # Set axis labels
  theme_classic() +
  scale_fill_manual(
    values = c("ratio" = "grey"))+
  geom_text(aes(label = sprintf("%.2f", ratio)), vjust = -0.3, size = 3, 
            position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(plot.title = element_text(size = 12, 
                                  face = "bold", hjust = 0.5)) +
  theme(text = element_text(size = 12),
        axis.text = element_text(colour = "black")) +
  theme(legend.position = "top", legend.direction = "horizontal")
ggsave('./plot/SCENIC/TF_ratio.pdf', width = 4.5, height = 4)


#### conserved TF #######
regulon.conserved <- list()
for (i in names(regulon.list)){
  regulon.conserved[[i]] <- intersect(regulon.list[[i]], regulon.list[[12]])
}


# ########## plot ################
loom <- open_loom('./SCENIC/others/Human.AUCell.loom')

regulon_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulon_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
close_loom(loom)
head(regulonAucThresholds)

regulonAucThresholds <- setNames(names(regulonAucThresholds), regulonAucThresholds)
head(regulonAucThresholds)

# This function will be included in the next version of AUCell
binarizeAUC <- function(auc, thresholds)
{
  thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
  regulonsCells <- setNames(lapply(names(thresholds), 
                                   function(x) {
                                     trh <- thresholds[x]
                                     names(which(getAUC(auc)[x,]>trh))
                                   }),names(thresholds))
  
  regulonActivity <- reshape2::melt(regulonsCells)
  binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
  class(binaryRegulonActivity) <- "matrix"  
  
  return(binaryRegulonActivity)
}

binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)

kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
head(kmeans_thresholds)
binary_regulons <- binarize_regulons(regulonAUC,kmeans_thresholds)

meta <- read.table('./SCENIC/others/Human_metadata.xls',sep='\t',header=T,stringsAsFactor=F)

cellinfo <- meta[,c('CellType','nGene', 'nUMI')]
colnames(cellinfo)=c('CellType', 'nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'CellType'))
selectedResolution <- "CellType"

sub_regulonAUC <- regulonAUC
sub_regulonAUC@assays@data@listData$AUC <- binaryRegulonActivity

rss <- calcRSS(AUC=getAUC(sub_regulonAUC), # binaryRegulonActivity
               cellAnnotation=cellTypes[colnames(regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)

rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
ggsave('./plot/SCENIC/rss_human.pdf', width = 4, height = 9)


options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
plotRSS_oneSet(rss, setName = "B cells") # cluster ID
rss_df <- rssPlot$df
colnames(rss_df) <- c('regulon','cell_type','RSS','Z')

head(rss_df)
plot_rrs_ranking(rss_df,
                 cell_type = "B cells",
                 ggrepel_force = 1,
                 ggrepel_point_padding = 0.2,
                 top_genes = 4,
                 plot_extended = FALSE)


plot.list <- list()
for (i in levels(rss_df$cell_type)){
  rrs_df_sub <- rss_df %>% subset(cell_type == i) %>% 
    arrange(desc(RSS))
  rrs_df_sub <- rrs_df_sub %>% mutate(rank = as.numeric(rownames(rrs_df_sub)))
  
  plot.list[[i]] <- 
  ggplot(rrs_df_sub, aes(rank, RSS, 
                         label = regulon)) +
    geom_point(color = "grey20",size = 2) +
    geom_point(data = subset(rrs_df_sub, 
                             rank < 4), color = "red", size = 2) +
    geom_text_repel(data = subset(rrs_df_sub, rank < 4), force = 1, point.padding = 0.2) + 
    labs(x = "Regulons", y = "Regulon Specificity Score", title = i) +
    theme_minimal() +  # 使用最小化主题
    theme(
      panel.border = element_rect(
        color = "black",  # 边框颜色
        fill = NA,  # 填充颜色为空（透明）
        linewidth = 1  # 边框宽度
      ),
      panel.grid.major = element_blank(),  # 移除主要网格线
      panel.grid.minor = element_blank(),  # 移除次要网格线
      axis.line = element_line(color = "black"),  # 坐标轴线颜色
      axis.text = element_text(size = 12, color = "black"),  # 坐标轴文本颜色和大小
      axis.title = element_text(size = 14, color = "black"),  # 坐标轴标题颜色和大小
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)  # 图表标题颜色、大小和位置
    )
}
plot.list[[1]]
ggarrange(
  plot.list[[5]] + labs(x = '', y = 'Regulon Specificity Score'),
  plot.list[[1]] + labs(x = '', y = ''),
  plot.list[[3]] + labs(x = '', y = ''),
  plot.list[[2]] + labs(x = '', y = ''),
  plot.list[[4]] + labs(x = '', y = ''),
  plot.list[[6]] + labs(x = '', y = ''),
  ncol = 6,
  align = c("v")
)

ggsave('./plot/SCENIC/ranking_tf_human.pdf', width = 12, height = 4)

regulon_cell <- binaryRegulonActivity %>% t() %>% as.data.frame()
regulon_cell$celltype <- cellinfo$CellType

regulon_cell <- lapply(X = unique(x = regulon_cell$celltype), FUN = function(ident) {
  data.use <- regulon_cell[regulon_cell$celltype == ident, 1:(ncol(x = regulon_cell) - 1), drop = FALSE]
  avg.AUC <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
    return(mean(x = x))
  })
  pct.AUC <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  res <- data.frame(id = ident, avg.AUC = avg.AUC, pct.AUC = pct.AUC * 100)
  res$gene <- rownames(res)
  return(res)
}) %>% do.call("rbind", .) %>% data.frame()

# Scale the average expression values
regulon_cell <- purrr::map_df(unique(regulon_cell$gene), function(x) {
  tmp <- regulon_cell %>% dplyr::filter(gene == x)
  avg.AUC.scale <- tmp %>% dplyr::select(avg.AUC) %>% scale(.) %>% scales::rescale(., to = c(0, 1))
  tmp$avg.AUC.scaled <- as.numeric(avg.AUC.scale)
  return(tmp)
})

head(rss_df)
head(regulon_cell)

merged_df <- rss_df %>%
  left_join(regulon_cell, by = c("regulon" = "gene", "cell_type" = "id"))

merged_df <- merged_df %>%
  group_by(cell_type) %>%
  mutate(RSSZ = scale(RSS))

head(merged_df)

ggplot(filter(merged_df,cell_type == 'B cells'), aes(x = RSSZ, y = avg.AUC)) +
  geom_point(color = "darkblue", size = 2.5, shape = 21) +  # 设置散点的颜色和大小
  labs(x = "RSSZ", y = "AUCell score") +  # 设置坐标轴标签
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      color = "black",  # 边框颜色
      fill = NA,  # 填充颜色为空（透明）
      linewidth = 1  # 边框宽度
    ),
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    axis.line = element_line(color = "black"),  # 坐标轴线颜色
    axis.text = element_text(size = 12, color = "black"),  # 坐标轴文本颜色和大小
    axis.title = element_text(size = 14, color = "black"),  # 坐标轴标题颜色和大小
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # 图表标题颜色、大小和位置
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "darkgreen")



library(ggridges)
rrs_df_nona <- subset(rss_df,RSS > 0)
ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 5, alpha = 0.75) +
  geom_vline(xintercept = 0.1) +
  theme(legend.position = "none")

regulon_csi <- calculate_csi(sub_regulonAUC,
                              calc_extended = FALSE)

head(regulon_csi)

plot_csi_modules(regulon_csi,
                 nclust = 10,
                 font_size_regulons = 8)

csi_matrix <- pivot_wider(regulon_csi, names_from = regulon_2, values_from = CSI, values_fill = NA) %>% as.data.frame()

# Set the row names as the first column (regulon_1)
rownames(csi_matrix) <- csi_matrix$regulon_1
head(csi_matrix)
csi_matrix <- dplyr::select(csi_matrix, -regulon_1)%>% as.matrix()
head(csi_matrix)
# csi_matrix <- csi_matrix[rowSums(csi_matrix) > 5.347594e-03,]
# csi_matrix <- csi_matrix[,colSums(csi_matrix) > 5.347594e-03]
head(csi_matrix)

factoextra::fviz_nbclust(csi_matrix, kmeans, method = "wss",k.max = 20)
km_csi <- kmeans(csi_matrix,7)

km_csi <- data.frame(TF = names(km_csi$cluster),
                     Module = km_csi$cluster)
km_csi$Module <- paste('Module',km_csi$Module,sep = ' ')
km_csi$Module <- factor(km_csi$Module, levels = c(paste(rep('Module'),seq(1,7,1),sep = ' ')))

openxlsx::write.xlsx(km_csi,'./plot/SCENIC/CSI_module_human.xlsx')
# c4a_gui()
col_use <- c4a('paired',7)
names(col_use) <- levels(km_csi$Module)

top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col_use, col = col_use), # 设置填充色
                       labels = levels(km_csi$Module),
                       labels_gp = gpar(cex = 0.3, col = col_use),
                       height = unit(0.2, "cm")
  ))
left_anno <- rowAnnotation(
  cluster = anno_block(gp = gpar(fill = col_use, col = col_use), # 设置填充色
                       labels = levels(km_csi$Module),
                       labels_gp = gpar(cex = 0.3, col = col_use),
                       width = unit(0.2, "cm")
  ))

pdf('./plot/SCENIC/heatmap_human_CSI.pdf',width = 7.65, height = 5.29)
Heatmap(csi_matrix,col = colorRampPalette(c("#F7F6D6","#242162"))(100),
        row_split = km_csi$Module,
        column_split = km_csi$Module,
        top_annotation = top_anno,
        left_annotation = left_anno,
        column_dend_height = unit(8, "mm"),
        row_dend_width = unit(8, "mm"),
        show_column_names = F,
        show_row_names = F,
        row_title_side ="right",
        row_title_gp = gpar(col = col_use),
        row_title_rot = 0,
        column_title = NULL,
        heatmap_legend_param = list(title = 'CSI',
                                    title_gp = gpar(fontsize = 10, 
                                                    fontface = "plain"),
                                    title_position = "topcenter", 
                                    border = NA,
                                    # legend_height = unit(40, "mm"),
                                    legend_width = unit(30, "mm"),
                                    labels_gp = gpar(fontsize = 9),
                                    legend_direction = "horizontal",
                                    grid_width = unit(4, "mm")))

dev.off()

######## go enrichment analysis #######
module.ego <- list()
for (i in levels(km_csi$Module)){
  filename <- paste0('./plot/SCENIC/GO_module/',i,'.csv')
  module.ego[[i]] <- read.csv(filename)
  module.ego[[i]]$Module <- i
}

module_ego <- do.call(rbind, module.ego)
module_ego$LogP <- -module_ego$LogP
module_ego <- filter(module_ego, LogP > 6)
head(module_ego)

module_ego <- slice_max(module_ego, by = Module, order_by = LogP, n=10)

module_matrix <- module_ego %>% #[,c('Description','Module','LogP')]
  pivot_wider(names_from = Module, values_from = LogP)

module_matrix <- module_matrix[,c(4,18:24)]

module_matrix <- module_matrix %>%
  group_by(Description) %>%
  summarize_all(~ sum(., na.rm = TRUE)) %>%
  ungroup() 

module_description <- module_matrix$Description

module_matrix <- module_matrix %>%
  as.data.frame() %>%
  select(-Description) %>%
  as.matrix()

rownames(module_matrix) <- module_description

wrapped_rownames <- str_wrap(rownames(module_matrix), width = 40)

head(module_matrix)
ht <- Heatmap(t(module_matrix),
              width = unit(28, "cm"), height = unit(4.5, "cm"),
              # column_labels = wrapped_rownames,
              row_names_gp = gpar(fontsize = 13),
              column_names_gp = gpar(fontsize = 13),
              heatmap_legend_param = list(
                title = '-LogP',
                direction = "horizontal",  # 图例方向设置为水平
                legend_width = unit(4, "cm")  # 调整图例宽度
              ))


pdf('./plot/SCENIC/heatmap_module_ego.pdf',width = 16,height = 10)

draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 20), "mm"))
dev.off()

####### regulons ######
regulon_filter <- filter(regulon_csi, CSI > 0.9)
regulon_filter$module <- km_csi$Module[match(regulon_filter$regulon_1,km_csi$TF)]
# regulon_filter <- dplyr::select(regulon_filter, -CSI)
regulon_filter$regulon_1 <- gsub("\\(\\+\\)", "", regulon_filter$regulon_1)
regulon_filter$regulon_2 <- gsub("\\(\\+\\)", "", regulon_filter$regulon_2)
regulon_filter$module <- as.character(regulon_filter$module)

head(regulon_filter)

colnames(regulon_filter) <- c('node1','node2','CSI','module')

write.csv(regulon_filter,'./plot/SCENIC/Human.grn.filter.csv',row.names = F)

#### heatmap ####
regulon_cell <- regulon_cell %>%
  left_join(km_csi, by = c("gene" = "TF"))

write.csv(regulon_cell, './plot/SCENIC/regulon_auc.csv')

regulon_module <- regulon_cell %>%
  group_by(Module) %>%
  filter(avg.AUC > 0.5)
  # slice_max(avg.AUC,n = 20)
  
regulon_module <- regulon_module %>%
  group_by(id,Module) %>%
  mutate(module.AUC = mean(avg.AUC))

regulon_module <- regulon_module[,c('id','Module','module.AUC')]

regulon_module <- regulon_module %>%
  distinct(id, Module, .keep_all = TRUE) %>%
  group_by(id) %>%
  mutate(module.AUC.scaled = scales::rescale(module.AUC))

regulon_matrix <- regulon_module %>%
  select(Module, id, module.AUC.scaled) %>%
  pivot_wider(names_from = Module, values_from = module.AUC.scaled) %>%
  as.data.frame()

rownames(regulon_matrix) <- regulon_matrix$id
regulon_matrix <- regulon_matrix[,-1]

# for(i in colnames(regulon_matrix)){
#   regulon_matrix[[i]] <- scales::rescale(regulon_matrix[[i]])
# }

regulon_matrix[is.na(regulon_matrix)] = 0

head(regulon_matrix)

# pdf('./plot/SCENIC/heatmap_avg_AUC.pdf', width = 7.65, height = 3.5)
Heatmap(regulon_matrix,cluster_columns = F,
        show_row_dend = F,
        col = colorRampPalette(c("#F7F6D6","darkred"))(100),
        heatmap_legend_param = list(title = 'AVG. Regulon Activity',
                                    title_gp = gpar(fontsize = 10, 
                                                    fontface = "plain"), title_position = "leftcenter-rot", 
                                    border = NA, legend_height = unit(40, "mm"), labels_gp = gpar(fontsize = 9),
                                    legend_direction = "vertical",
                                    grid_width = unit(4, "mm")))
# dev.off()

openxlsx::write.xlsx(regulon_module, './plot/SCENIC/module_select_regulon.xlsx')

binaryRegulonActivity[1:3,1:3]

col_use <- c4a('paired',7)
names(col_use) <- levels(km_csi$Module)
km_csi$ModuleColor <- col_use[km_csi$Module]
col_module <- as.vector(km_csi$ModuleColor)
names(col_module) <- km_csi$Module


col_cell <- c("B cells" = "#fc8002", 
               "Monocytes" = "#369f2d",
               "mDCs" = "#846e89", 
               "NK cells" = "#fac7b3",
               "pDCs" = "#4995c6", 
               "T cells" = "#ee4431")
cellinfo$CellColor <- col_cell[cellinfo$CellType]
col_type <- as.vector(cellinfo$CellColor)
names(col_type) <- rownames(cellinfo)


top_anno <- HeatmapAnnotation(Celltype = anno_simple(rownames(cellinfo),
                                                    col = col_type,
                              height = unit(0.2, "cm")))

left_anno <- rowAnnotation(Module = anno_simple(as.vector(km_csi$Module),
                                                  col = col_module,
                                                  width = unit(0.2, "cm")))

legend.list <- list(
  Legend(labels = names(col_cell), legend_gp = gpar(fill = col_cell, col = col_cell)),
  Legend(labels = names(col_use), legend_gp = gpar(fill = col_use, col = col_use))
)

binary_dist <- dist(binaryRegulonActivity)
binary_cluster <- hclust(binary_dist)
plot(binary_cluster,hang = -1,cex=0.6,axes=FALSE,ann=FALSE)
annotation_genes <- data.frame(genes = c('BCL11A(+)','SPIB(+)','IRF8(+)',
                                         'TCF7(+)','FOXO1(+)','ETS1(+)',
                                         'EOMES(+)','TBX21(+)','RUNX3(+)',
                                         'RXRA(+)','CEBPA(+)','NFE2(+)',
                                         'VDR(+)','SPI1(+)','ARID3A(+)',
                                         'IRF7(+)','MYCN(+)',
                                         'FOS(+)','FOSL2(+)','CEBPD(+)',
                                         'KLF4(+)','CUX1(+)','BACH1(+)','TCF7L2(+)',
                                         'SMAD3(+)','MYB(+)','RARG(+)','HOXA1(+)','SOX10(+)','ETV1(+)',
                                         'SOX4(+)','GATA6(+)','JUN(+)','SREBF1(+)'))

heat <- Heatmap(binaryRegulonActivity,
                col = c("white","black"),
                # row_split = km_csi$Module,
                # column_split = cellinfo$CellType,
                gap = unit(0.5,'mm'),
                # cluster_rows = binary_cluster,
                row_gap = unit(0.5,'mm'),
                show_column_names = F,
                show_row_names = T,
                show_column_dend = F,
                row_title = NULL,
                top_annotation = top_anno,
                left_annotation = left_anno,
                heatmap_legend_param = list(title = '',
                                            title_gp = gpar(fontsize = 10,fontface = "plain"),
                                            title_position = "topleft",
                                            border = 1, legend_height = unit(40, "mm"),
                                            at = c(0, 1),
                                            labels = c('off', 'on'),
                                            labels_gp = gpar(fontsize = 9),
                                            legend_gp = gpar(fill = 1:2),
                                            legend_direction = "vertical",
                                            grid_width = unit(4, "mm"),
                                            ncol = 2)) 
  # rowAnnotation(link = anno_mark(at = which(rownames(binaryRegulonActivity) %in% annotation_genes$genes),
  #                                labels = annotation_genes$genes, labels_gp = gpar(fontsize = 10)))

# pdf('./plot/SCENIC/heatmap_human_activity_rownames.pdf',width = 7.65, height = 22)

draw(heat, annotation_legend_list = legend.list)

# dev.off()

############# module of different species ##########
head(km_csi)
head(regulon.list)
# Create an empty data frame
regulon_df <- km_csi

for (i in names(regulon.list)){
  regulon_animal <- data.frame(TF = regulon.list[[i]], regulon = regulon.list[[i]])
  colnames(regulon_animal)[2] <- i
  
  regulon_df <- left_join(regulon_df, regulon_animal, by = 'TF')
}

head(regulon_df)

regulon_ortho <- cbind(regulon_df[ ,1:2],regulon_df[ ,3:14] %>% mutate_all(~ ifelse(is.na(.), 0, 1)))

regulon_ortho <- regulon_ortho %>%
    group_by(Module) %>%
    summarize(
      Catfish = sum(Catfish),
      Jacopever = sum(Jacopever),
      Turtle = sum(Turtle),
      Chicken = sum(Chicken),
      Bat = sum(Bat),
      Pig = sum(Pig),
      Cattle = sum(Cattle),
      Mouse = sum(Mouse),
      Rat = sum(Rat),
      Monkey = sum(Monkey),
      Chimpanzee = sum(Chimpanzee),
      Human = sum(Human)
    )

head(regulon_ortho)

regulon_ortho <- regulon_ortho %>%
  pivot_longer(
    cols = -Module,
    names_to = "Specie",
    values_to = "Count"
  )

regulon_ortho$Specie <- factor(regulon_ortho$Specie, levels = rev(species))

ggplot(regulon_ortho, aes(x = Module, y = Specie)) +
  geom_tile(aes(fill = Count), color = "white") +
  geom_text(aes(label = Count), vjust = 1, size = 3) +
  scale_fill_gradient(low = "white", high = "red") +  # Adjust colors as needed
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave('./')


data <- regulon_ortho

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$Module), ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$Module <- rep(levels(data$Module), each=empty_bar)
data <- rbind(data, to_add)
data <- data %>% arrange(Module)
data$id <- seq(1, nrow(data))

# Get the name and the y position of each label
label_data <- data
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(Module) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))


module_angle <- label_data %>%
  group_by(Module) %>%
  summarize(angle = mean(angle, na.rm = TRUE)) %>%
  as.data.frame()

module_angle$angle <- module_angle$angle - c(90,90,90,180,270,270,270) + c(5,5,5,10,5,5,5)

base_data <- left_join(x = base_data,
                       y = module_angle)  
  
module_ratio <- regulon_ortho %>%
  group_by(Module) %>%
  summarize(ratio = ((sum(Count) - max(Count))/11)/max(Count) * 100) %>% # 计算比对到人的平均比例
  as.data.frame()

module_ratio$ratio <- round(module_ratio$ratio, digits = 2)
module_ratio$ratio <- paste0(module_ratio$ratio,'%')

base_data <- left_join(x = base_data,
                       y = module_ratio)

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=Count, fill=Module)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=Count, fill=Module), stat="identity", alpha=0.5) +
  scale_fill_manual(values = c4a('paired',7)) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 10, xend = start, yend = 10), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data$id),4), y = c(0,10,20,40), label = c("0", "10", "20", "40") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=Count, fill=Module), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  geom_text(data=label_data, aes(x=id, y=Count+10, label=Specie, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # 
  # # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=Module), hjust=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5),
            vjust = c(0.5,0.5,0.5,0.5,0.5,0.5,0.5), colour = "black",
            alpha=0.8, size=3, angle = base_data$angle, fontface="bold", inherit.aes = FALSE) +
  geom_text(data=base_data, aes(x = title, y = -18, label=ratio), hjust=c(0.5,0.5,0.5,0.5,0.5,0.5,0.5),
            vjust = c(2,2,2,2,2,2,2), colour = "darkred",
            alpha=0.8, size=3, angle = base_data$angle, fontface="bold", inherit.aes = FALSE)


p

ggsave('./plot/SCENIC/module_ortholog_human.pdf', width = 9, height = 9)

### 查看一下module ####
head(km_csi)
head(regulon.list)
# Create an empty data frame
regulon_df <- km_csi

for (i in names(regulon.list)){
  regulon_animal <- data.frame(TF = regulon.list[[i]], regulon = regulon.list[[i]])
  colnames(regulon_animal)[2] <- i
  
  regulon_df <- left_join(regulon_df, regulon_animal, by = 'TF')
}

head(regulon_df)
head(regulon_df)

regulon_df$Count <- rowSums(!is.na(regulon_df[, 3:14]))

openxlsx::write.xlsx(regulon_df, './plot/SCENIC/module_all_species.xlsx')

head(regulon_df)

module_conserved <- filter(regulon_df, Module == 'Module 4')
module_conserved$TF <- factor(module_conserved$TF, levels = module_conserved[order(module_conserved$Count), ]$TF)

ggplot(module_conserved, aes(x = TF, y = Count)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  geom_text(aes(label = Count), colour = 'darkred', hjust = 0.5, vjust = 0, size = 4) +
  labs(x = '', y = '') +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(), # 去掉图形背景
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave('./plot/SCENIC/mudole_conserved_regulon.pdf', width = 6.5, height = 4.5)


########### DIMPLOT ###########
# binaryRegulonActivity来自human
regulon_conserved <- binaryRegulonActivity[c("MYB(+)","ERG(+)","HLTF(+)","HMGA1(+)"),] %>% t() %>% as.data.frame()

regulon_conserved <- binaryRegulonActivity[c("TCF7L2(+)","SPI1(+)","FOSL2(+)","FOS(+)"),] %>% t() %>% as.data.frame()


# sp.list <- readRDS('./rds/sp.list.quality_filter.rds')

human.umap <- SeuratObject::Embeddings(sp.list[['Human']], reduction = 'umap') %>% as.data.frame()

# rm(sp.list)

head(regulon_conserved)
head(human.umap)
regulon_conserved$cellid <- rownames(regulon_conserved)
human.umap$cellid <- paste(rep('Human',nrow(human.umap)),rownames(human.umap),sep = '_')

regulon_conserved <- left_join(regulon_conserved, human.umap)

gg1 <- ggplot() +
  geom_point(data = regulon_conserved, aes(x = UMAP_1, y = UMAP_2, color = factor(`MYB(+)`)), size = 0.5) +
  scale_color_manual(values = c('0' = 'grey', '1' = 'darkred'), labels = c('0' = 'OFF', '1' = 'ON')) +
  labs(title = colnames(regulon_conserved)[1]) +
  labs(color = NULL) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom',
    legend.direction = "horizontal",
  )+ NoLegend()
gg1
# 创建图例图
gg2 <- ggplot() +
  geom_point(data = regulon_conserved, aes(x = UMAP_1, y = UMAP_2, color = factor(`ERG(+)`)), size = 0.5) +
  scale_color_manual(values = c('0' = 'grey', '1' = 'darkred'), labels = c('0' = 'OFF', '1' = 'ON')) +
  labs(title = colnames(regulon_conserved)[2]) +
  labs(color = NULL) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom',
    legend.direction = "horizontal",
  )+ NoLegend()

gg3 <- ggplot() +
  geom_point(data = regulon_conserved, aes(x = UMAP_1, y = UMAP_2, color = factor(`HLTF(+)`)), size = 0.5) +
  scale_color_manual(values = c('0' = 'grey', '1' = 'darkred'), labels = c('0' = 'OFF', '1' = 'ON')) +
  labs(title = colnames(regulon_conserved)[3]) +
  labs(color = NULL) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'bottom',
    legend.direction = "horizontal",
  )+ NoLegend()

gg4 <- ggplot() +
  geom_point(data = regulon_conserved, aes(x = UMAP_1, y = UMAP_2, color = factor(`HMGA1(+)`)), size = 0.5) +
  scale_color_manual(values = c('0' = 'grey', '1' = 'darkred'), labels = c('0' = 'OFF', '1' = 'ON')) +
  labs(title = colnames(regulon_conserved)[4]) +
  labs(color = NULL) + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'right'
  )

ggarrange(gg1,gg2,gg3,gg4,
          ncol = 4,
          widths = c(1,1,1,1.3))

ggsave('./plot/SCENIC/module1_regulon_dimplot.pdf', width = 14, height = 3.5)



# sp.list <- readRDS('./rds/sp.list.quality_filter.rds')
human.umap <- SeuratObject::Embeddings(sp.list[['Human']], reduction = 'umap') %>% as.data.frame()
human.umap$cellid <- paste(rep('Human',nrow(human.umap)),rownames(human.umap),sep = '_')

head(human.umap)

plot.list <- list()
for (i in levels(as.factor(km_csi$Module))){
  regulon_module <- binaryRegulonActivity[filter(km_csi, Module == i)$TF,] %>% t() %>% as.data.frame()
  
  regulon_module[1:3,1:3]
  
  regulon_module <- regulon_module %>% mutate(avg.RAS = rowMeans(.))
  
  # rm(sp.list)
  regulon_module$cellid <- rownames(regulon_module)
  
  regulon_module <- left_join(regulon_module, human.umap)
  
  plot.list[[i]]<- 
    ggplot() +
    geom_point(data = regulon_module, aes(x = UMAP_1, y = UMAP_2, color = avg.RAS), size = 0.5) +
    scale_color_gradient(low = "grey", high = "#3C0D03") +
    labs(title = i) +
    labs(color = NULL) + 
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = 'right'
    )
  
}
plot.list[[1]]+plot.list[[2]]+plot.list[[3]]+plot.list[[4]]+plot.list[[5]]+plot.list[[6]]+plot.list[[7]]
ggsave('./plot/SCENIC/module_umap.pdf',width = 7.64, height = 6)

######## target gene ########
regulon_target <- read.table('./SCENIC/others/Human.grn.tsv', header = T, sep = '\t')

regulon_target_module <- filter(regulon_target, TF %in% c("MYB","ERG","HLTF","HMGA1"))

openxlsx::write.xlsx(regulon_target_module, './plot/SCENIC/grn_conserved_module1.xlsx')

regulon_target_module <- filter(regulon_target, TF %in% c("TCF7L2","SPI1","FOSL2","FOS"))

openxlsx::write.xlsx(regulon_target_module, './plot/SCENIC/grn_conserved_module4.xlsx')


############## snkey plot ###########

for (i in species){
  loomfile <- paste0('./SCENIC/others/',i,'.AUCell.loom')
  loom <- open_loom(loomfile)
  
  regulon_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulon_incidMat)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  regulonAucThresholds <- get_regulon_thresholds(loom)
  close_loom(loom)
  head(regulonAucThresholds)
  
  regulonAucThresholds <- setNames(names(regulonAucThresholds), regulonAucThresholds)
  head(regulonAucThresholds)
  
  # This function will be included in the next version of AUCell
  binarizeAUC <- function(auc, thresholds)
  {
    thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
    regulonsCells <- setNames(lapply(names(thresholds), 
                                     function(x) {
                                       trh <- thresholds[x]
                                       names(which(getAUC(auc)[x,]>trh))
                                     }),names(thresholds))
    
    regulonActivity <- reshape2::melt(regulonsCells)
    binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
    class(binaryRegulonActivity) <- "matrix"  
    
    return(binaryRegulonActivity)
  }
  
  binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
  
  # kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
  # head(kmeans_thresholds)
  # binary_regulons <- binarize_regulons(regulonAUC,kmeans_thresholds)
  
  metafile <- paste0('./SCENIC/others/',i,'_metadata.xls')
  meta <- read.table(metafile,sep='\t',header=T,stringsAsFactor=F)
  
  cellinfo <- meta[,c('CellType','nGene', 'nUMI')]
  colnames(cellinfo)=c('CellType', 'nGene' ,'nUMI')
  cellTypes <-  as.data.frame(subset(cellinfo,select = 'CellType'))
  selectedResolution <- "CellType"

  bi_sankey <- binaryRegulonActivity %>%
    as.data.frame()
  
  bi_sankey$TF <- rownames(bi_sankey)
  
  km_csi <- openxlsx::read.xlsx('./plot/SCENIC/CSI_module_human.xlsx')
  
  bi_sankey <- left_join(bi_sankey, km_csi)
  bi_sankey <- na.omit(bi_sankey)
  
  head(bi_sankey)
  
  for (k in 1:nrow(bi_sankey)) {
    for (j in 1:(ncol(bi_sankey)-2)) {
      bi_sankey[k, j] <- ifelse(bi_sankey[k, j] == 1, bi_sankey[k, 'Module'], NA)
    }
  }
  
  bi_sankey <- select(bi_sankey, -c(TF, Module)) %>% as.data.frame()
  
  bi_sankey <- bi_sankey[, !sapply(bi_sankey, function(x) all(is.na(x)))] # 去除全部为 NA 的列
  
  head(bi_sankey)
  
  bi_sankey_melted <- data.frame()
  
  # 遍历 bi_sankey 数据框的列
  for (col in colnames(bi_sankey)) {
    # 提取非缺失值的 Module 列
    module_values <- bi_sankey[, col][!is.na(bi_sankey[, col])]
    
    # 创建一个数据框，其中包含两列：CellID 和 Module
    cell_id <- col
    data <- data.frame(CellID = cell_id, Module = module_values)
    
    # 将数据添加到结果数据框中
    bi_sankey_melted <- rbind(bi_sankey_melted, data)
  }
  
  head(bi_sankey_melted)
  
  cellinfo$CellID <- rownames(cellinfo)
  bi_sankey_melted <- left_join(bi_sankey_melted, cellinfo[,c('CellType','CellID')])
  
  bi_sankey_melted <- na.omit(bi_sankey_melted)
  # bi_sankey$Module <- paste(rep('Module',nrow(bi_sankey)), bi_sankey$Module, sep = ' ')
  
  bi_sankey_melted <- bi_sankey_melted %>% mutate_if(is.character, as.factor) # convert to factor
  nodes <- data.frame(node = c(bi_sankey_melted$Module, bi_sankey_melted$CellType)) %>% unique() #node label
  bi_sankey_melted$IDsource <- match(bi_sankey_melted$Module, nodes$node) - 1 # start from 0
  bi_sankey_melted$IDtarget <- match(bi_sankey_melted$CellType, nodes$node) - 1 # start from 0
  bi_sankey_melted$value <- 1
  head(bi_sankey_melted)
  
  # Create the Sankey diagram
  sankey_plot <- sankeyNetwork(Links = bi_sankey_melted, Nodes = nodes, 
                               colourScale = 'd3.scaleOrdinal().domain(["Module 1","Module 2","Module 3","Module 4","Module 5","Module 6","Module 7","B cells","Monocytes","mDCs","NK cells","pDCs","T cells"]).range(["#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#fc8002","#369f2d","#846e89","#fac7b3","#4995c6","#ee4431"])',
                               Source = "IDsource", Target = "IDtarget", Value = "value", NodeID = "node", 
                               fontSize = 18, nodeWidth = 30)
  
  
  # sankey_plot
  
  sankeyplot <- paste0('./plot/SCENIC/sankeyplot_',i,'.html')
  sankeyplot.pdf <- paste0('./plot/SCENIC/sankeyplot_',i,'.pdf')
  saveNetwork(sankey_plot, file = sankeyplot)
  webshot::webshot(sankeyplot, file =  sankeyplot.pdf, vwidth = 450, vheight = 600)
}


for (i in species){
sankeyplot.pdf <- paste0('./plot/SCENIC/sankeyplot_',i,'.pdf')
sankeyplot.png <- paste0('./plot/SCENIC/sankeyplot_',i,'.png')
pdftools::pdf_convert(sankeyplot.pdf, format = 'png', pages = 1, sankeyplot.png,  dpi = 300)
}


rss.list <- list()
for (i in species){
  loomfile <- paste0('./SCENIC/others/',i,'.AUCell.loom')
  loom <- open_loom(loomfile)
  
  regulon_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulon_incidMat)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  regulonAucThresholds <- get_regulon_thresholds(loom)
  close_loom(loom)
  head(regulonAucThresholds)
  
  regulonAucThresholds <- setNames(names(regulonAucThresholds), regulonAucThresholds)
  head(regulonAucThresholds)
  
  # This function will be included in the next version of AUCell
  binarizeAUC <- function(auc, thresholds)
  {
    thresholds <- thresholds[intersect(names(thresholds), rownames(auc))]
    regulonsCells <- setNames(lapply(names(thresholds), 
                                     function(x) {
                                       trh <- thresholds[x]
                                       names(which(getAUC(auc)[x,]>trh))
                                     }),names(thresholds))
    
    regulonActivity <- reshape2::melt(regulonsCells)
    binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
    class(binaryRegulonActivity) <- "matrix"  
    
    return(binaryRegulonActivity)
  }
  
  binaryRegulonActivity <- binarizeAUC(regulonAUC, regulonAucThresholds)
  
  metafile <- paste0('./SCENIC/others/',i,'_metadata.xls')
  meta <- read.table(metafile,sep='\t',header=T,stringsAsFactor=F)
  
  cellinfo <- meta[,c('CellType','nGene', 'nUMI')]
  colnames(cellinfo)=c('CellType', 'nGene' ,'nUMI')
  cellTypes <-  as.data.frame(subset(cellinfo,select = 'CellType'))
  selectedResolution <- "CellType"
  
  sub_regulonAUC <- regulonAUC[,colnames(binaryRegulonActivity)]
  sub_regulonAUC@assays@data@listData$AUC <- binaryRegulonActivity
  
  rss <- calcRSS(AUC=getAUC(sub_regulonAUC), # binaryRegulonActivity
                 cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                          selectedResolution])
  rss=na.omit(rss)
  
  rss.list[[i]] <- scale(rss)
  
}


plot.list <- list()
for(i in species){
  annotation_genes <- data.frame(genes = rownames(rss.list[[i]])[apply(rss.list[[i]], 1, function(row) any(row > 3))])
  plot.list[[i]] <- Heatmap(rss.list[[i]],
                            circlize::colorRamp2(c(-4, 0, 4), c("white", "white", "red")),
                            heatmap_legend_param = list(title = 'RSSZ',
                                                        title_gp = gpar(fontsize = 10,fontface = "plain"),
                                                        title_position = "topleft")) +
    rowAnnotation(link = anno_mark(at = which(rownames(rss.list[[i]]) %in% annotation_genes$genes), 
                                   labels = annotation_genes$genes, labels_gp = gpar(fontsize = 10)))
}

plot.list[[2]]

for(i in species){
  filename <- paste0('./plot/SCENIC/RSSZ_heatmap_', i, '.pdf')
  pdf(filename, width = 3.5, height = 6)
  print(plot.list[[i]])
  dev.off()
}

all_cell <- c('Monocytes','T cells','NK cells','B cells','mDCs','pDCs','DCs')
all_name <- c('tfd','sslg','pss','gga','ray','ssc',
              'bta','mmu','rno','mcc','ptr','hsa')

names(all_name) <- species

rss.merge <- list()
for(i in rev(species)){
  rss.merge[[i]] <- rss.list[[i]] %>% as.data.frame()
  colnames(rss.merge[[i]]) <- paste(colnames(rss.merge[[i]]),rep(all_name[i],ncol(rss.merge[[i]])),sep = '.')
  rss.merge[[i]]$TF <- rownames(rss.merge[[i]])
}

head(rss.merge[[1]])

rss.merge <- rss.merge %>%
  reduce(left_join, by = 'TF')

rownames(rss.merge) <- rss.merge$TF

rss.merge <- select(rss.merge, -TF)

rss.merge[is.na(rss.merge)] <- -4

all_column <- list()
for(i in all_cell){
  all_column[[i]] <- paste(i, all_name, sep = '.')
}

all_column <- unlist(all_column)
all_column <- all_column[all_column %in% colnames(rss.merge)]

rss.merge <- rss.merge[,all_column]

rss.merge <- as.matrix(rss.merge)

species_count <- apply(rss.merge, 1, function(row) sum(row > 1, na.rm = TRUE))
selected_genes <- names(species_count[species_count > 6])

annotation_genes <- data.frame(genes = rownames(rss.merge)[apply(rss.merge, 1, function(row) any(row > 4))])

all_name <- c('tfd','sslg','pss','gga','ray','ssc',
              'bta','mmu','rno','mcc','ptr','hsa')
species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')
cell_type <- sub("\\..*$", "", colnames(rss.merge))
species_abb <- sub(".*\\.", "", colnames(rss.merge))
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


pdf('./plot/SCENIC/RSSZ_heatmap_merged.pdf', width = 12, height = 8)
Heatmap(rss.merge,
        cluster_columns = F,
        bottom_annotation = hb,
        show_column_names = F,
        circlize::colorRamp2(c(-4, 0, 4), c("white", "white", "red")),
        width = unit(20, "cm"), height = unit(20, "cm"),
        heatmap_legend_param = list(title = 'RSSZ',
                                    legend_height = unit(50, "mm"),
                                    title_gp = gpar(fontsize = 10,fontface = "plain"),
                                    title_position = "topleft")) +
  rowAnnotation(link = anno_mark(at = which(rownames(rss.merge) %in% selected_genes),
                                 labels = selected_genes, labels_gp = gpar(fontsize = 10)))
dev.off()

pdf('./plot/SCENIC/RSSZ_heatmap_merged_trans.pdf', width = 22, height = 12)
Heatmap(t(rss.merge),
        cluster_rows = F,
        circlize::colorRamp2(c(-4, 0, 4), c("white", "white", "red")),
        heatmap_legend_param = list(title = 'RSSZ',
                                    legend_height = unit(50, "mm"),
                                    title_gp = gpar(fontsize = 10,fontface = "plain"),
                                    title_position = "topleft"))
dev.off()

########### regulons ############
regulons.list <- list()
for (i in species){
  loomfile <- paste0('./SCENIC/others/',i,'.AUCell.loom')
  loom <- open_loom(loomfile)
  regulon_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons.list[[i]] <- regulonsToGeneLists(regulon_incidMat)
  close_loom(loom)
}

for(i in names(regulons.list)){
  regulons.list[[i]] <- regulons.list[[i]][c("FOSL2(+)","FOS(+)","TCF7L2(+)","SPI1(+)")]
}


regulons.data <- list()
for(i in names(regulons.list[[1]])){
  # regulons.data[[i]] <- data.frame(
  #   Catfish = character(),
  #   Jacopever = character(),
  #   Turtle = character(),
  #   Chicken = character(),
  #   Bat = character(),
  #   Pig = character(),
  #   Cattle = character(),
  #   Mouse = character(),
  #   Rat = character(),
  #   Monkey = character(),
  #   Chimpanzee = character(),
  #   Human = character()
  # )
  regulons.data[[i]] <- list()
   for(j in names(regulons.list)){
     regulons.data[[i]][[j]] <- regulons.list[[j]][[i]]
   }
}


head(regulons.data[[1]])


for(i in names(regulons.data)){
  # Find the maximum length of subsets
  max_length <- max(sapply(regulons.data[[i]], length))
  
  # Create a list with the same structure as regulons.data but filled with NA
  regulons.data[[i]] <- lapply(regulons.data[[i]], function(subset) {
    if (length(subset) < max_length) {
      c(subset, rep(NA, max_length - length(subset)))
    } else {
      subset
    }
  })
}

head(regulons.data[[1]])

for(i in names(regulons.data)){
  regulons.data[[i]] <- as.data.frame(regulons.data[[i]])
}

for(i in names(regulons.data)){
  filename <- paste0('./plot/SCENIC/conserved_regulon_',i,'.xlsx')
  openxlsx::write.xlsx(regulons.data[[i]], filename)
  filename <- paste0('./plot/SCENIC/conserved_regulon_',i,'.csv')
  write.csv(regulons.data[[i]], filename, row.names = F)
}

View(regulons.data[[1]])


####### 提取转录因子靶基因 #######
regulons.conserved <- list()
for(i in names(regulons.list[[1]])){
  regulons.conserved[[i]] <- list()
  for(j in names(regulons.list)){
    regulons.conserved[[i]][[j]] <- regulons.list[[j]][[i]]
  }
}


head(regulons.conserved[[1]])


for(i in names(regulons.conserved)){
  for(j in names(regulons.conserved[[i]])){
    regulons.conserved[[i]][[j]] <- as.data.frame(regulons.conserved[[i]][[j]])
    colnames(regulons.conserved[[i]][[j]]) <- j
    regulons.conserved[[i]][[j]]$target <- regulons.conserved[[i]][[j]][[j]]
  } 
}

regulons_conserved <- list()

for(i in names(regulons.conserved)){
  regulons_conserved[[i]] <- regulons.conserved[[i]][[1]]
}

for(i in names(regulons.conserved)){
  for(j in 2:12){
    regulons_conserved[[i]] <- merge(regulons_conserved[[i]], regulons.conserved[[i]][[j]], by = 'target', all = T)
  }
}

for(i in names(regulons.conserved)){
  regulons_conserved[[i]]$count <- rowSums(!is.na(regulons_conserved[[i]]))
}

for(i in names(regulons.conserved)){
  regulons_conserved[[i]] <- filter(regulons_conserved[[i]], count > 6)
}

regulon.target <- list()
for (i in species){
  filename <- paste0('./SCENIC/others/',i,'.grn.tsv')
  regulon.target[[i]] <- read.table(filename, header = T, sep = '\t')
  regulon.target[[i]] <- filter(regulon.target[[i]], TF %in% c("TCF7L2","SPI1","FOSL2","FOS"))
}

regulon.target <- do.call(rbind, regulon.target)

regulon.target$pair <- paste(regulon.target$TF, regulon.target$target, sep = '_')

head(regulon.target)

regulon.target <- regulon.target[duplicated(regulon.target$pair) == 'FALSE',]

regulon.target <- dplyr::select(regulon.target, -pair)

write.csv(regulon.target, './plot/SCENIC/conserved_regulon_monocyte.csv')
openxlsx::write.xlsx(regulon.target, './plot/SCENIC/conserved_regulon_monocyte.xlsx')




