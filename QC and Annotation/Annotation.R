lapply(c("dplyr","Seurat"), library, character.only = T)
library(ggplot2)
library(patchwork)
library(future)
library(stringr)
library(cols4all)
library(dplyr)
library(tidyr)
library(harmony)

########## Calculate the number of PCs that contain some proportion (95%) of the variance #######
npcs <- function(
    object,
    var.toal=0.95,
    reduction="pca"
){
  if(is.null(object@reductions[[reduction]])){
    cat("Reduction", reduction, "not found!")
    return(NULL)
  }
  
  tmp.var <- (object@reductions[[reduction]]@stdev)^2
  var.cut <- var.toal*sum(tmp.var)
  n.pcs=0
  var.sum = 0
  while(var.sum < var.cut){
    n.pcs = n.pcs + 1
    var.sum <- var.sum + tmp.var[n.pcs]
  }
  
  return(n.pcs)
}

setwd('~/data/PBMC/')

plan('multisession', workers = 8)
plan()

sp.list <- readRDS('./rds/sp.list.rds')
################## catfish ###############
catfish <- merge(sp.list[['Yellow_catfish1']],sp.list[['Yellow_catfish2']],
               add.cell.ids = c('Yellow_catfish1','Yellow_catfish2'),project = 'Yellow_catfish')

meta <- catfish@meta.data
meta[,'orig.ident'] = 'Yellow_catfish'
meta$sample <- grepl("^Yellow_catfish1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Yellow_catfish1', 'Yellow_catfish2')
catfish@meta.data <- meta

DefaultAssay(catfish) <- 'RNA'
genes <- rownames(catfish) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
catfish_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_mouse.xlsx',sheet = 'Catfish')
genes$mouse <- catfish_ortho$Mouse.gene.name[match(genes$gene,catfish_ortho$Catfish.gene.name)]
genes$newgene <- NA
for (i in 1:nrow(genes)){
  if(is.na(genes[i,2])=='TRUE'){
    genes[i,3] = genes[i,1]
  }else{
    genes[i,3] = genes[i,2]
  }
}

genes$LOC <- grepl('^LOC[[:digit:]]',genes$gene)
table(genes$LOC)
table(is.na(genes$mouse))

genes$rp <- grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",genes$newgene)
genes$mt <- grepl("^mt-|^Mt-|^MT-",genes$newgene)
table(genes$rp)

genes$dup <- duplicated(genes)
table(genes$dup)

table(grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",rownames(catfish[["RNA"]])))
table(grepl("^mt-|^Mt-|^MT-",rownames(catfish[["RNA"]])))
catfish[["percent.rp"]]  = PercentageFeatureSet(catfish, pattern = "Rp[sl][[:digit:]]|^rp[sl][[:digit:]]")
catfish[["percent.mt"]]  = PercentageFeatureSet(catfish, pattern = "^mt-|^Mt-|^MT-")
dim(catfish)
summary(catfish[["nCount_RNA"]]  )
summary(catfish[["nFeature_RNA"]])
summary(catfish[["percent.mt"]]  )
summary(catfish[["percent.rp"]]  )
VlnPlot(catfish,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(catfish, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(catfish, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(catfish, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3

catfish$quality <- ifelse(catfish$nFeature_RNA < 300, "Low quality", "High quality")

genes <- rownames(catfish) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",rownames(genes))
genes$mt <- grepl("^mt-|^Mt-|^MT-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')

catfish <- catfish[genes$gene,]

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(catfish))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(catfish))
catfish <- CellCycleScoring(catfish,g2m.features = g2m_genes,s.features = s_genes)
table(catfish$Phase)

catfish <- catfish %>%
  SCTransform(method = 'glmGamPoi',vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

catfish <- catfish %>% RunHarmony(group.by.vars = 'sample')

DimPlot(catfish,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_catfish.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = catfish, reduction = 'harmony')

catfish <- catfish %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

job::job({
  plan('multisession', workers = 4)
  Idents(catfish) <- 'seurat_clusters'
  catfish.marker <- FindAllMarkers(catfish)
  catfish_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Catfish')
  catfish_anno <- read.csv('./ncbi/catfish_annotation.csv')
  catfish.marker$human <- catfish_ortho$Human.gene.name[match(catfish.marker$gene,catfish_ortho$Catfish.gene.name)]
  catfish.marker$description <- catfish_anno$product[match(catfish.marker$gene,catfish_anno$gene)]
  openxlsx::write.xlsx(catfish.marker,'./file/markers_seurat_clusters_Catfish.xlsx')
},title = 'catfish_FindMarker')

DimPlot(catfish, label = T, repel = T)
ggsave('./plot/cluster/seurat_clusters_Catfish.pdf')

library(SingleR)
library(celldex)
mouseRNA <- MouseRNAseqData()

sct_matrix <- GetAssayData(catfish,slot = 'data')
clusters <- catfish$seurat_clusters
genes <- rownames(sct_matrix) %>% as.data.frame()
colnames(genes) <- 'gene'
genes$mouse <- catfish_ortho$Mouse.gene.name[match(genes$gene,catfish_ortho$Yellow.catfish.gene.name)]

rownames(sct_matrix) <- genes$mouse

catfish_singler <- SingleR(test = sct_matrix, ref = mouseRNA,
                           labels = mouseRNA$label.main,
                           clusters = clusters,
                           assay.type.test = 'logcounts',
                           assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(catfish$seurat_clusters),
                       cell_annotation=catfish_singler$labels)

celltype

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_costume_mouse.xlsx',"Peripheral blood")

scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- catfish$seurat_clusters[match(rownames(cellmax),rownames(catfish@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)
celltype

table(catfish$seurat_clusters)

catfish.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Catfish.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)

markers <- catfish.marker %>% filter(cluster %in% c('1')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$mouse,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=20)
head(ego,n=10L)

kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Mm.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)


VlnPlot(catfish,features = 'cd79a') #B cells   B cell receptor CD22-like
VlnPlot(catfish,features = 'LOC113663309') #B cells   B cell receptor CD22-like
VlnPlot(catfish,features = 'flt3') #T cells   M1-specific T cell receptor alpha chain-like
VlnPlot(catfish,features = 'fn1a') #Monocytes
VlnPlot(catfish,features = 'marco') #Macrophages
VlnPlot(catfish,features = 'ncf1') #Neutrophils
VlnPlot(catfish,features = 'LOC113646159') #Platelets ## platelet glycoprotein V
VlnPlot(catfish,features = 'LOC113656766') #Erythrocytes ## hemoglobin subunit alpha-like
VlnPlot(catfish,features = 'mafb')

FeaturePlot(catfish,features = 'cd79a') #B cells   B cell receptor CD22-like
FeaturePlot(catfish,features = 'cd79b') #B cells   B cell receptor CD22-like
FeaturePlot(catfish,features = 'ikzf2') #T cells ## Cd3e
FeaturePlot(catfish,features = 'LOC113643771') #mDCs   CD209
FeaturePlot(catfish,features = 'ikzf2') #HSCs
FeaturePlot(catfish,features = 'mafb') # pDCs
FeaturePlot(catfish,features = 'tcf7') #T cells
FeaturePlot(catfish,features = 'csf1ra') #Monocytes
FeaturePlot(catfish,features = 'marco') #Macrophages
FeaturePlot(catfish,features = 'ncf1') #Neutrophils mpx
FeaturePlot(catfish,features = 'LOC113646159') #Platelets ## platelet glycoprotein V
FeaturePlot(catfish,features = 'LOC113656766') #Erythrocytes ## hemoglobin subunit alpha-like

FeaturePlot(catfish,features = 'mef2b')
FeaturePlot(catfish,features = 'nkl.4')
FeaturePlot(catfish,features = 'fn1a')
FeaturePlot(catfish,features = 'flt3')
FeaturePlot(catfish,features = 'gp1bb') #Platelets

Idents(catfish) <- 'seurat_clusters'
catfish <- RenameIdents(catfish,
                        '17'='B cells','19'='B cells',
                        '1'='HSCs','10'='HSCs','14'='HSCs',
                        '15'='T cells',
                        '20'='DCs',
                        '13'='Monocytes',
                        '4'='Monocytes','16'='Monocytes','18'='Monocytes',
                        '2'='Neutrophils','3'='Neutrophils','5'='Neutrophils','7'='Neutrophils','8'='Neutrophils',
                        '9'='Neutrophils','11'='Neutrophils',
                        '6'='Platelets',
                        '12'='Erythrocytes') 

cluster_letters <- Idents(catfish)
names(cluster_letters) <- colnames(catfish)
catfish <- AddMetaData(object = catfish,
                         metadata = cluster_letters,
                         col.name = "cell_annotation")

# add low quality and doublet 
metadata <- catfish@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

catfish$multi_annotation <- metadata$multi_annotation
table(catfish$multi_annotation)

# c4a_gui()
DimPlot(catfish, group.by = 'multi_annotation')
ggsave('./plot/cluster/catfish_cell_cluster.pdf')
saveRDS(catfish,'./rds/Catfish_raw.rds')

catfish <- readRDS('./rds/Catfish_raw.rds')

################## jacopever #######################
jacopever <- merge(sp.list[['Jacopever1']],sp.list[['Jacopever2']],
                   add.cell.ids = c('Jacopever1','Jacopever2'),project = 'Jacopever')

meta <- jacopever@meta.data
meta[,'orig.ident'] = 'Jacopever'
meta$sample <- grepl("^Jacopever1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Jacopever1', 'Jacopever2')
jacopever@meta.data <- meta

DefaultAssay(jacopever) <- 'RNA'
genes <- rownames(jacopever) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
jacopever_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_mouse.xlsx',sheet = 'Jacopever')
genes$mouse <- jacopever_ortho$Mouse.gene.name[match(genes$gene,jacopever_ortho$Jacopever.gene.name)]
genes$newgene <- NA
for (i in 1:nrow(genes)){
  if(is.na(genes[i,2])=='TRUE'){
    genes[i,3] = genes[i,1]
  }else{
    genes[i,3] = genes[i,2]
  }
}

genes$LOC <- grepl('^LOC[[:digit:]]',genes$gene)
table(genes$LOC)

genes$rp <- grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",genes$newgene)
genes$mt <- grepl("^mt-",genes$newgene)


genes$dup <- duplicated(genes)


table(grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",rownames(jacopever[["RNA"]])))
table(grepl("^mt-",rownames(jacopever[["RNA"]])))
jacopever[["percent.rp"]]  = PercentageFeatureSet(jacopever, pattern = "Rp[sl][[:digit:]]|^rp[sl][[:digit:]]")
jacopever[["percent.mt"]]  = PercentageFeatureSet(jacopever, pattern = "^mt-")
dim(jacopever)
summary(jacopever[["nCount_RNA"]]  )
summary(jacopever[["nFeature_RNA"]])
summary(jacopever[["percent.mt"]]  )
summary(jacopever[["percent.rp"]]  )
VlnPlot(jacopever,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(jacopever, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(jacopever, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(jacopever, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3


jacopever$quality <- ifelse(jacopever$nFeature_RNA < 300, "Low quality", "High quality")

genes <- rownames(jacopever) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]",rownames(genes))
genes$mt <- grepl("^mt-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')

jacopever <- jacopever[genes$gene,]

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(jacopever))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(jacopever))
jacopever <- CellCycleScoring(jacopever,g2m.features = g2m_genes,s.features = s_genes)
table(jacopever$Phase)

jacopever <- jacopever %>%
  SCTransform(method = 'glmGamPoi',vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

jacopever <- jacopever %>% RunHarmony(group.by.vars = 'sample')

DimPlot(jacopever,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_jacopever.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = jacopever, reduction = 'harmony')

jacopever <- jacopever %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

job::job({
  plan('multisession', workers = 4)
  Idents(jacopever) <- 'seurat_clusters'
  jacopever.marker <- FindAllMarkers(jacopever,)
  openxlsx::write.xlsx(jacopever.marker,'./file/markers_seurat_clusters_Jacopever.xlsx')
},title = 'jacopever_FindMarker')

jacopever_anno <- read.csv('./ncbi/jacopever_annotation.csv')
jacopever_ortho <- openxlsx::read.xlsx('./protein/ortholog_to_human.xlsx',sheet = 'Jacopever')
jacopever.marker$human <- jacopever_ortho$Human.gene.name[match(jacopever.marker$gene,jacopever_ortho$Jacopever.gene.name)]
jacopever.marker$description <- jacopever_anno$product[match(jacopever.marker$gene,jacopever_anno$gene)]

openxlsx::write.xlsx(jacopever.marker,'./file/markers_seurat_clusters_Jacopever.xlsx')

DimPlot(jacopever, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Jacopever.pdf')


jacopever.marker <- openxlsx::read.xlsx('./file/markers_Jacopever.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)

markers <- jacopever.marker %>% filter(cluster %in% c('17')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$mouse,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=20)
head(ego,n=10L)

kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Mm.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)

VlnPlot(jacopever,features = 'LOC119484985') #B cells   B cell receptor CD22-like
VlnPlot(jacopever,features = 'LOC119498499') #T cells   M1-specific T cell receptor alpha chain-like
VlnPlot(jacopever,features = 'nkl.4') #NK cells
VlnPlot(jacopever,features = 'fgl2a') #Monocytes
VlnPlot(jacopever,features = 'marco') #Macrophages
VlnPlot(jacopever,features = 'flt3') #DCs
VlnPlot(jacopever,features = 'mctp1a') #Neutrophils

FeaturePlot(jacopever,features = 'cd79a')
FeaturePlot(jacopever,features = 'cd79b')
FeaturePlot(jacopever,features = 'LOC119484985') #B cells   B cell receptor CD22-like
FeaturePlot(jacopever,features = 'tcf7') #T cells   cd3g
# FeaturePlot(jacopever,features = 'nkl.4') #NK cells
FeaturePlot(jacopever,features = 'marco') #Monocytes
FeaturePlot(jacopever,features = 'marco') #Macrophages
FeaturePlot(jacopever,features = 'flt3') #mDCs
FeaturePlot(jacopever,features = 'ncf1') #Neutrophils
FeaturePlot(jacopever,features = 'top2a') #Cycling cells
FeaturePlot(jacopever,features = 'gp1bb') #Platelets

FeaturePlot(jacopever,features = 'LOC119490960') # cd3g
FeaturePlot(jacopever,features = 'LOC119481023')
FeaturePlot(jacopever,features = 'mef2b')

Idents(jacopever) <- 'seurat_clusters'
jacopever <- RenameIdents(jacopever,
                       '2'='B cells','3'='B cells','4'='B cells','5'='B cells','8'='B cells','12'='B cells',
                       '15'='B cells',
                       '16'='Cycling cells',
                       '6'='T cells','11'='T cells',
                       '10'='NK cells',
                       '13'='Monocytes',
                       '9'='Neutrophils','19'='Neutrophils',
                       '1'='Platelets','7'='Platelets','14'='Platelets',
                       '17'='DCs','18'='DCs') 

cluster_letters <- Idents(jacopever)
names(cluster_letters) <- colnames(jacopever)
jacopever <- AddMetaData(object = jacopever,
                      metadata = cluster_letters,
                      col.name = "cell_annotation")

metadata <- jacopever@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

jacopever$multi_annotation <- metadata$multi_annotation


#c4a_gui()
DimPlot(jacopever, group.by = 'multi_annotation', cols = c4a('redon',12), label = T, repel = T)
ggsave('./plot/cluster/jacopever_cell_cluster.pdf')

saveRDS(jacopever,'./rds/Jacopever_raw.rds')

jacopever <- readRDS('./rds/Jacopever_raw.rds')

################## turtle ##############
turtle  <- merge(sp.list[['Chinese_softshell_turtle1']],sp.list[['Chinese_softshell_turtle2']],
               add.cell.ids = c('Chinese_softshell_turtle1','Chinese_softshell_turtle2'),project = 'Chinese_softshell_turtle')

meta <- turtle@meta.data
meta[,'orig.ident'] = 'Chinese_softshell_turtle'
meta$sample <- grepl("^Chinese_softshell_turtle1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Chinese_softshell_turtle1', 'Chinese_softshell_turtle2')
turtle@meta.data <- meta


turtle_gtf <- read.table("./ncbi/genome/GCF_000230535.1.gtf/ncbi_dataset/data/GCF_000230535.1/genomic.gtf", sep="\t", header=FALSE, stringsAsFactors=FALSE)

turtle_gtf <- cbind(turtle_gtf[1:8], t(sapply(strsplit(turtle_gtf$V9, "; "), function(x) {
  sapply(c("gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "exon_number"), function(y) {
    tmp <- grep(y, x, fixed=TRUE)
    if (length(tmp) == 0) {
      return("")
    } else {
      return(sub(paste0("^", y, " "), "", x[tmp]))
    }
  })
})))


colnames(turtle_gtf) <- c("chromsome", "source", "feature", "start", "end", "score", "strand", "frame",
                       "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "exon_number")

turtle_gtf <- filter(turtle_gtf,feature == 'gene') %>% filter(gene_biotype == 'protein_coding;')

mt_genes <- filter(turtle_gtf,chromsome == 'MT')$gene_name

DefaultAssay(turtle) <- 'RNA'
genes <- rownames(turtle) %>% as.data.frame()
colnames(genes) <- 'gene'
rownames(genes) <- genes$gene
genes$newgene <- NA
for (i in 1:nrow(genes)){
  if (genes[i, 'gene'] %in% mt_genes){
    genes$newgene[i] <- paste('MT', genes[i, 'gene'], sep = '-')
  }else{
    genes$newgene[i] = genes$gene[i]
  }
}

genes$mt <- grepl("^MT-",genes$newgene)
genes$rp <- grepl("^RP[SL][[:digit:]]",genes$newgene)
table(genes$mt)
table(genes$rp)
genes$LOC <- grepl('^LOC[[:digit:]]',genes$gene)
table(genes$LOC)

table(grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]|^RP[SL][[:digit:]]",rownames(turtle[["RNA"]])))
table(grepl("^MT-",rownames(turtle[["RNA"]])))
turtle[["percent.rp"]]  = PercentageFeatureSet(turtle, pattern = "^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]|^RP[SL][[:digit:]]")
turtle[["percent.mt"]]  = PercentageFeatureSet(turtle, pattern = "^MT-")
dim(turtle)
summary(turtle[["nCount_RNA"]]  )
summary(turtle[["nFeature_RNA"]])
summary(turtle[["percent.mt"]]  )
summary(turtle[["percent.rp"]]  )

VlnPlot(turtle,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(turtle, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(turtle, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(turtle, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3


turtle$quality <- ifelse(turtle$nFeature_RNA < 300, "Low quality", "High quality")

genes <- rownames(turtle) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^Rp[sl][[:digit:]]|^rp[sl][[:digit:]]|^RP[SL][[:digit:]]",rownames(genes))
genes$mt <- grepl("^mt-|^Mt-|^MT-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')

turtle <- turtle[genes$gene,]

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(turtle))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(turtle))
turtle <- CellCycleScoring(turtle,g2m.features = g2m_genes,s.features = s_genes)
table(turtle$Phase)

turtle <- turtle %>%
  SCTransform(method = 'glmGamPoi',vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

turtle <- turtle %>% RunHarmony(group.by.vars = 'sample')

DimPlot(turtle,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_turtle.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = turtle, reduction = 'harmony')

turtle <- turtle %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

job::job({
  plan('multisession', workers = 4)
  Idents(turtle) <- 'seurat_clusters'
  turtle.marker <- FindAllMarkers(turtle)
  turtle.marker$mouse <- turtle_ortho$Mouse.gene.name[match(turtle.marker$gene,turtle_ortho$Chinese.softshell.turtle.symbol)]
  openxlsx::write.xlsx(turtle.marker,'./file/markers_Turtle.xlsx')
},title = 'turtle_FindMarker')


DimPlot(turtle, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Turtle.pdf')


library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(turtle,slot = 'data')
clusters <- turtle$seurat_clusters

human_singler <- SingleR(test = sct_matrix, ref = immRNA,
                         labels = immRNA$label.main,
                         clusters = clusters,
                         assay.type.test = 'logcounts',
                         assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(turtle$seurat_clusters),
                       cell_annotation=human_singler$labels)


lapply(c("dplyr","Seurat","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")

scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- turtle$seurat_clusters[match(rownames(cellmax),rownames(turtle@meta.data))]

sum_df <- cellmax %>%
  group_by(seurat_clusters, cell_annotation) %>%
  summarise(n = n()) %>%
  ungroup()

count_df <- cellmax %>%
  group_by(seurat_clusters) %>%
  summarise(total_cells = n()) %>%
  ungroup()

result_df <- sum_df %>%
  left_join(count_df, by = "seurat_clusters") %>%
  mutate(prop = n / total_cells)

result_df <- result_df %>%
  group_by(seurat_clusters) %>%
  slice_max(prop,n=1)

result_df <- left_join(result_df,celltype,by='seurat_clusters')
colnames(result_df)[2] <- 'sctype'
colnames(result_df)[6] <- 'singleR'

result_df$sctype <- paste(result_df$sctype,'s',sep = '')

result_df$match <- result_df$sctype == result_df$singleR

result_df

table(turtle$seurat_clusters)

turtle.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Turtle.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)

markers <- turtle.marker %>% filter(cluster %in% c('1','2','3','7','9','13','18')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$mouse,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=20)
head(ego,n=10L)
openxlsx::write.xlsx(ego,'./file/enrichment/turtle_seu_1.ego.xlsx')

kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Mm.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)

openxlsx::write.xlsx(kk,'./file/enrichment/turtle_seu_1.ekegg.xlsx')

VlnPlot(turtle,features = 'LOC112545017') #B cells
VlnPlot(turtle,features = 'CD3E') #T cells
VlnPlot(turtle,features = 'JCHAIN') #Plasma cells
VlnPlot(turtle,features = 'CD44') #Macrophages
VlnPlot(turtle,features = 'F10') #Monocytes  ###Fc gamma R-mediated phagocytosis

FeaturePlot(turtle,features = 'LOC112545017') #B cells
FeaturePlot(turtle,features = 'CD79B') #B cells
FeaturePlot(turtle,features = 'CD3E') #T cells
FeaturePlot(turtle,features = 'JCHAIN') #Plasma cells
FeaturePlot(turtle,features = 'CSF1R') #Monocytes
FeaturePlot(turtle,features = 'MEIS1') #HSCs
FeaturePlot(turtle,features = 'TGFBI') 

FeaturePlot(turtle,features = '')
Idents(turtle) <- 'seurat_clusters'
turtle <- RenameIdents(turtle,
                      '6'='B cells',
                      '17'='B cells',
                      '14'='T cells',
                      '4'='Monocytes','5'='Monocytes','8'='Monocytes','10'='Monocytes','11'='Monocytes',
                      '12'='Monocytes','15'='Monocytes','16'='Monocytes',
                      '1'='HSCs','2'='HSCs','3'='HSCs','7'='HSCs','9'='HSCs',
                      '13'='HSCs','18'='HSCs') 

cluster_letters <- Idents(turtle)
names(cluster_letters) <- colnames(turtle)
turtle <- AddMetaData(object = turtle,
                       metadata = cluster_letters,
                       col.name = "cell_annotation")

metadata <- turtle@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

turtle$multi_annotation <- metadata$multi_annotation


#c4a_gui()
DimPlot(turtle, group.by = 'multi_annotation')
ggsave('./plot/cluster/turtle_cell_cluster.pdf')

saveRDS(turtle,'./rds/Turtle_raw.rds')

turtle <- readRDS('./rds/Turtle_raw.rds')

################## chicken ######################
chicken <- merge(sp.list[['Chicken1']],sp.list[['Chicken2']],
               add.cell.ids = c('Chicken1','Chicken2'),project = 'Chicken')

meta <- chicken@meta.data
meta[,'orig.ident'] = 'Chicken'
meta$sample <- grepl("^Chicken1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Chicken1', 'Chicken2')
chicken@meta.data <- meta

DefaultAssay(chicken) <- 'RNA'
table(grepl("^RP[SL][[:digit:]]",rownames(chicken[["RNA"]])))
table(grepl("^J6367-",rownames(chicken[["RNA"]])))
chicken[["percent.rp"]]  = PercentageFeatureSet(chicken, pattern = "^RP[SL][[:digit:]]")
chicken[["percent.mt"]]  = PercentageFeatureSet(chicken, pattern = "^J6367-")
dim(chicken)
summary(chicken[["nCount_RNA"]]  )
summary(chicken[["nFeature_RNA"]])
summary(chicken[["percent.mt"]]  )
summary(chicken[["percent.rp"]]  )

VlnPlot(chicken,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(chicken, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(chicken, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(chicken, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3

chicken$quality <- ifelse(chicken$percent.mt > 20 | chicken$nFeature_RNA < 300, "Low quality", "High quality")


genes <- rownames(chicken) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^RP[SL][[:digit:]]",rownames(genes))
genes$mt <- grepl("^J6367-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')

chicken <- chicken[genes$gene,]

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(chicken))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(chicken))
chicken <- NormalizeData(chicken) #
chicken <- CellCycleScoring(chicken,g2m.features = g2m_genes,s.features = s_genes)
table(chicken$Phase)

chicken <- chicken %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

chicken <- chicken %>% RunHarmony(group.by.vars = 'sample')

DimPlot(chicken,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_chicken.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = chicken, reduction = 'harmony')

chicken <- chicken %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

DimPlot(chicken,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  chicken.marker <- FindAllMarkers(chicken)
  openxlsx::write.xlsx(chicken.marker,'./file/markers_seurat_clusters_Chicken.xlsx')
},title = 'chicken_FindMarker')
chicken.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Chicken.xlsx',sheet = 'Sheet 1')

DimPlot(chicken,group.by = 'sample')
DimPlot(chicken, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Chicken.pdf')

library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(chicken,slot = 'data')
clusters <- chicken$seurat_clusters

chicken_singler <- SingleR(test = sct_matrix, ref = immRNA,
                         labels = immRNA$label.main,
                         clusters = clusters,
                         assay.type.test = 'logcounts',
                         assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(chicken$seurat_clusters),
                       cell_annotation=chicken_singler$labels)

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")

scRNAseqData = sct_matrix %>% as.matrix() 
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

View(es.max)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- chicken$seurat_clusters[match(rownames(cellmax),rownames(chicken@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

table(chicken$seurat_clusters)
celltype

chicken.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Chicken.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Gg.eg.db)
markers <- chicken.marker %>% filter(cluster %in% c('21')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Gg.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Gg.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=40)
head(ego,n=10L)

kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'gga',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Gg.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)

FeaturePlot(chicken,features = 'HBA1')

DimPlot(chicken, label = T, repel = T)
ggsave('./plot/cluster/seurat_clusters_chicken.pdf')

VlnPlot(chicken,features = 'CD79B') #B cells
VlnPlot(chicken,features = 'TCF7') #T cells
VlnPlot(chicken,features = 'GNLY') #NK cells
VlnPlot(chicken,features = 'EXFABP') #Monocytes
VlnPlot(chicken,features = 'NCF1C') #
VlnPlot(chicken,features = 'IL5RA') #
VlnPlot(chicken,features = 'TUBB1') #Platelets
VlnPlot(chicken,features = 'JCHAIN') #Plasma cells
VlnPlot(chicken,features = 'HBA1') #Erythrocytes
VlnPlot(chicken,features = 'ERC2') #Lymphoid progenitor cellss

FeaturePlot(chicken,features = 'CD79B') #B cells
FeaturePlot(chicken,features = 'TCF7') #T cells
FeaturePlot(chicken,features = 'GNLY') #NK cells
FeaturePlot(chicken,features = 'CSF1R') #Monocytes
FeaturePlot(chicken,features = 'CSF3R') # Neutrophils
FeaturePlot(chicken,features = 'CSF2RB') #pDCs
FeaturePlot(chicken,features = 'GP1BB') #Platelets
FeaturePlot(chicken,features = 'JCHAIN') #Plasma cells
FeaturePlot(chicken,features = 'HBA1') #Erythrocytes
FeaturePlot(chicken,features = 'MAF') #T cells

FeaturePlot(chicken,features = 'EXFABP') #
FeaturePlot(chicken,features = 'JCHAIN') #Plasma cells

FeaturePlot(chicken,features = 'CTSG')

Idents(chicken) <- 'seurat_clusters'
chicken <- RenameIdents(chicken,
                        '1'='T cells','6'='T cells','10'='T cells','11'='T cells','12'='T cells','15'='T cells',
                        '5'='T cells',
                        '9'='NK cells',
                        '14'='B cells','16'='B cells',
                        '3'='Monocytes','19'='Monocytes',
                        '17'='HSCs',
                        '21'='pDCs',
                        '2'='Platelets','4'='Platelets','7'='Platelets','8'='Platelets',
                        '18'='Platelets','20'='Platelets',
                        '13'='B cells',
                        '23'='Erythrocytes')

cluster_letters <- Idents(chicken)
names(cluster_letters) <- colnames(chicken)
chicken <- AddMetaData(object = chicken,
                   metadata = cluster_letters,
                   col.name = "cell_annotation")

# add low quality and doublet 
metadata <- chicken@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '22'){
    metadata[i,4] = 'Doublets'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

chicken$multi_annotation <- metadata$multi_annotation
table(chicken$multi_annotation)

#c4a_gui()
DimPlot(chicken, group.by = 'multi_annotation')
ggsave('./plot/cluster/chicken_cell_cluster.pdf')

saveRDS(chicken,'./rds/Chicken_raw.rds')

chicken <- readRDS('./rds/Chicken_raw.rds')

################## bat #################
bat <- merge(sp.list[['Bat1']], y = c(sp.list[['Bat2']], sp.list[['Bat3']], sp.list[['Bat4']]),
                add.cell.ids = c('Bat1','Bat2','Bat3','Bat4'),project = 'Bat')

meta <- bat@meta.data
meta[,'orig.ident'] = 'Bat'
head(str_extract(rownames(meta),'Bat\\d+'))
meta$sample <- str_extract(rownames(meta),'Bat\\d+')
bat@meta.data <- meta

## percent.mt
DefaultAssay(bat) <- 'RNA'
table(grepl("^RP[SL][[:digit:]]",rownames(bat[["RNA"]])))
table(grepl("^MT-",rownames(bat[["RNA"]])))
bat[["percent.mt"]]  = PercentageFeatureSet(bat, pattern = "^MT-")
bat[["percent.rp"]]  = PercentageFeatureSet(bat, pattern = "^RP[SL][[:digit:]]")
dim(bat)
summary(bat[["nCount_RNA"]]  )
summary(bat[["nFeature_RNA"]])
summary(bat[["percent.rp"]])

VlnPlot(bat,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.rp"), ncol = 3)

plot1 <- FeatureScatter(bat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(bat, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2

# low quality cell 
bat$quality <- ifelse(bat$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp and mt genes
genes <- rownames(bat) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^RP[SL][[:digit:]]",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE')

## remove rp and mt genes
bat <- bat[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(bat))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(bat))
bat <- NormalizeData(bat) #
bat <- CellCycleScoring(bat,g2m.features = g2m_genes,s.features = s_genes)
table(bat$Phase)

bat <- bat %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()


bat <- bat %>% RunHarmony(group.by.vars = 'sample')

DimPlot(bat,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_bat.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = bat, reduction = 'harmony')

bat <- bat %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()
DimPlot(bat,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  bat.marker <- FindAllMarkers(bat)
  openxlsx::write.xlsx(bat.marker,'./file/markers_seurat_clusters_Bat.xlsx')
},title = 'bat_FindMarker')

bat_ortho <- openxlsx::read.xlsx('./protein/ortholog_all_species.xlsx',sheet = 'Egyptian.fruit.bat')
bat.marker$mouse <- bat_ortho$Mouse.gene.name[match(bat.marker$gene,bat_ortho$Egyptian.fruit.bat.gene.name)]
openxlsx::write.xlsx(bat.marker,'./file/markers_seurat_clusters_Bat.xlsx')

# annotation 
DimPlot(bat, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Bat.pdf')

# singleR annotation 
library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(bat,slot = 'data')
clusters <- bat$seurat_clusters

bat_singler <- SingleR(test = sct_matrix, ref = immRNA,
                       labels = immRNA$label.main,
                       clusters = clusters,
                       assay.type.test = 'logcounts',
                       assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(bat$seurat_clusters),
                       cell_annotation=bat_singler$labels)


lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")

scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)


cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- bat$seurat_clusters[match(rownames(cellmax),rownames(bat@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

celltype

table(bat$seurat_clusters)

bat.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Bat.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)
markers <- bat.marker %>% filter(cluster %in% c('16')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$mouse,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Mm.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=20)
head(ego,n=10L)
openxlsx::write.xlsx(ego,'./file/enrichment/bat_seu_28.ego.xlsx')

kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Mm.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)

openxlsx::write.xlsx(kk,'./file/enrichment/bat_seu_28.ekegg.xlsx')

VlnPlot(bat,features = 'HBB') #Erythrocytes

VlnPlot(bat,features = 'MS4A1') #B cells
VlnPlot(bat,features = 'CD3E') #T cells
VlnPlot(bat,features = 'GNLY') #NK cells
VlnPlot(bat,features = 'CSF1R') #Monocytes
VlnPlot(bat,features = 'FCER1A') #DCs
VlnPlot(bat,features = 'CD59') #Hematopoietic stem cells
VlnPlot(bat,features = 'TYMS') #Proliferative T cells
VlnPlot(bat,features = 'HOXA9') #Progenitor cells

FeaturePlot(bat,features = 'CD79A') # B cells
FeaturePlot(bat,features = 'CD3E') # T cells
FeaturePlot(bat,features = 'GNLY') # NK cells
FeaturePlot(bat,features = 'CSF1R') # Monocytes
FeaturePlot(bat,features = 'FCER1A') # mDCs
FeaturePlot(bat,features = 'GP1BB') # Platelets
FeaturePlot(bat,features = 'HBM') # Erythrocytes

FeaturePlot(bat,features = 'CD8B')
FeaturePlot(bat,features = 'FLT3')
FeaturePlot(bat,features = 'FCER1A')
FeaturePlot(bat,features = 'JCHAIN')

Idents(bat) <- 'seurat_clusters'
bat <- RenameIdents(bat,
                    '1'='B cells','4'='B cells','6'='B cells','10'='B cells','12'='B cells','13'='B cells','15'='B cells',
                    '17'='B cells','20'='B cells','22'='B cells','24'='B cells',
                    '7'='T cells','8'='T cells','9'='T cells','19'='T cells','29'='T cells',
                    '3'='NK cells','11'='NK cells','25'='NK cells',
                    '2'='Monocytes','5'='Monocytes','14'='Monocytes','16'='Monocytes','21'='Monocytes','23'='Monocytes',
                    '26'='Monocytes','27'='Monocytes','31'='Monocytes',
                    '30'='mDCs',
                    '28'='Platelets',
                    '18'='Erythrocytes','32'='Erythrocytes'
                    ) 

cluster_letters <- Idents(bat)
names(cluster_letters) <- colnames(bat)
bat <- AddMetaData(object = bat,
                   metadata = cluster_letters,
                   col.name = "cell_annotation")

metadata <- bat@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}


bat$multi_annotation <- metadata$multi_annotation
table(bat$multi_annotation)

#c4a_gui()
DimPlot(bat, group.by = 'multi_annotation', cols = c4a('redon',12))
ggsave('./plot/cluster/bat_cell_cluster.pdf')

saveRDS(bat,'./rds/Bat_raw.rds')

bat <- readRDS('./rds/Bat_raw.rds')

################## cattle #################
cattle <- merge(sp.list[['Cattle1']],sp.list[['Cattle2']],
                 add.cell.ids = c('Cattle1','Cattle2'),project = 'Cattle')

meta <- cattle@meta.data
meta[,'orig.ident'] = 'Cattle'
meta$sample <- grepl("^Cattle1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Cattle1', 'Cattle2')
cattle@meta.data <- meta

## percent.mt
DefaultAssay(cattle) <- 'RNA'
table(grepl("^RP[SL][[:digit:]]",rownames(cattle[["RNA"]])))
table(grepl("^MT-",rownames(cattle[["RNA"]])))
cattle[["percent.rp"]]  = PercentageFeatureSet(cattle, pattern = "^RP[SL][[:digit:]]")
cattle[["percent.mt"]]  = PercentageFeatureSet(cattle, pattern = "^MT-")
dim(cattle)
summary(cattle[["nCount_RNA"]]  )
summary(cattle[["nFeature_RNA"]])
summary(cattle[["percent.mt"]]  )
summary(cattle[["percent.rp"]]  )

VlnPlot(cattle,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(cattle, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(cattle, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(cattle, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 + plot3
# low quality cell
cattle$quality <- ifelse(cattle$percent.mt > 20 | cattle$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp and mt genes
genes <- rownames(cattle) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^RP[SL][[:digit:]]",rownames(genes))
genes$mt <- grepl("^MT-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')
## remove rp and mt genes
cattle <- cattle[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(cattle))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(cattle))
cattle <- NormalizeData(cattle) #
cattle <- CellCycleScoring(cattle,g2m.features = g2m_genes,s.features = s_genes)
table(cattle$Phase)

cattle <- cattle %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

cattle <- cattle %>% RunHarmony(group.by.vars = 'sample')

DimPlot(cattle,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_cattle.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = cattle, reduction = 'harmony')

cattle <- cattle %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>% # resolution = c(seq(0.1,1,0.1))
  identity()

DimPlot(cattle,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  cattle.marker <- FindAllMarkers(cattle)
  openxlsx::write.xlsx(cattle.marker,'./file/markers_seurat_clusters_Cattle.xlsx')
},title = 'cattle_FindMarker')
cattle.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Cattle.xlsx',sheet = 'Sheet 1')

# annotation
DimPlot(cattle, group.by = 'sample')
DimPlot(cattle, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Cattle.pdf')

# singleR annotation
library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(cattle,slot = 'data')
clusters <- cattle$seurat_clusters

# auto annotation
cattle_singler <- SingleR(test = sct_matrix, ref = immRNA,
                           labels = immRNA$label.main,
                           clusters = clusters,
                           assay.type.test = 'logcounts',
                           assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(cattle$seurat_clusters),
                       cell_annotation=cattle_singler$labels)

# sctype
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")
# assign cell types
scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

View(es.max)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- cattle$seurat_clusters[match(rownames(cellmax),rownames(cattle@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

DimPlot(cattle, label = T, repel = T)
ggsave('./plot/cluster/seurat_clusters_Cattle.pdf')

cattle.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Cattle.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Bt.eg.db)
markers <- cattle.marker %>% filter(cluster %in% c('20')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Bt.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Bt.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=40)
head(ego,n=10L)



kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'bta',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Bt.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)

VlnPlot(cattle,features = 'CD79B') #B cells
VlnPlot(cattle,features = 'CD3E') #T cells
VlnPlot(cattle,features = 'GNLY') #NK cells
VlnPlot(cattle,features = 'CD14') #Monocytes
VlnPlot(cattle,features = 'FCER1A') #mDCs
VlnPlot(cattle,features = 'IL5RA') #pDCs
VlnPlot(cattle,features = 'TUBB1') #Platelets
VlnPlot(cattle,features = 'JCHAIN') #Plasma cells
VlnPlot(cattle,features = 'HBA1') #Erythrocytes
VlnPlot(cattle,features = 'SLC4A4') #innate HSCs
VlnPlot(cattle,features = 'HBB') #Erythrocytes

FeaturePlot(cattle,features = 'CD79B') #B cells
FeaturePlot(cattle,features = 'CD3E') #T cells
FeaturePlot(cattle,features = 'GNLY') #NK cells
FeaturePlot(cattle,features = 'CD14') #Monocytes
FeaturePlot(cattle,features = 'FCER1A')#mDCs
FeaturePlot(cattle,features = 'FLT3')#pCDs
FeaturePlot(cattle,features = 'CXCR1') #Neutrophils
FeaturePlot(cattle,features = 'JCHAIN') #Plasma cells
FeaturePlot(cattle,features = 'HBB') #Erythrocytes

FeaturePlot(cattle,features = 'SLC4A4') #innate HSCs

FeaturePlot(cattle,features = 'FLT3') 
FeaturePlot(cattle,features = 'AXL') 
FeaturePlot(cattle,features = 'CLEC4C')

Idents(cattle) <- 'seurat_clusters'
cattle <- RenameIdents(cattle,
                        '1'='T cells','2'='T cells','3'='T cells','6'='T cells','8'='T cells','9'='T cells',
                        '18'='T cells','22'='T cells',
                        '12'='NK cells','14'='NK cells','23'='NK cells',
                        '5'='B cells','7'='B cells','10'='B cells','13'='B cells','16'='B cells','17'='B cells',
                        '4'='Monocytes','15'='Monocytes','22'='Monocytes',
                        '21'='mDCs',
                        '20'='pDCs',
                        '11'='Neutrophils',
                        '19'='Erythrocytes')

cluster_letters <- Idents(cattle)
names(cluster_letters) <- colnames(cattle)
cattle <- AddMetaData(object = cattle,
                       metadata = cluster_letters,
                       col.name = "cell_annotation")

metadata <- cattle@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '22'){
    metadata[i,4] = 'Doublets'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

cattle$multi_annotation <- metadata$multi_annotation
table(cattle$multi_annotation)

# c4a_gui()
DimPlot(cattle, group.by = 'multi_annotation')
ggsave('./plot/cluster/cattle_cell_cluster.pdf')

saveRDS(cattle,'./rds/Cattle_raw.rds')

cattle <- readRDS('./rds/Cattle_raw.rds')

################## pig ##################
pig <- merge(sp.list[['Pig1']],sp.list[['Pig2']],
                add.cell.ids = c('Pig1','Pig2'),project = 'Pig')

meta <- pig@meta.data
meta[,'orig.ident'] = 'Pig'
meta$sample <- grepl("^Pig1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Pig1', 'Pig2')
pig@meta.data <- meta


pig_gtf <- read.table("./ensembl/Sus_scrofa.Sscrofa11.1.109.chr.gtf", sep="\t", header=FALSE, stringsAsFactors=FALSE)

pig_gtf <- cbind(pig_gtf[1:8], t(sapply(strsplit(pig_gtf$V9, "; "), function(x) {
  sapply(c("gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "exon_number"), function(y) {
    tmp <- grep(y, x, fixed=TRUE)
    if (length(tmp) == 0) {
      return("")
    } else {
      return(sub(paste0("^", y, " "), "", x[tmp]))
    }
  })
})))

colnames(pig_gtf) <- c("chromsome", "source", "feature", "start", "end", "score", "strand", "frame",
                       "gene_id", "gene_name", "gene_biotype", "transcript_id", "transcript_name", "transcript_biotype", "exon_number")

pig_gtf <- filter(pig_gtf,feature == 'gene') %>% filter(gene_biotype == 'protein_coding;')

mt_genes <- filter(pig_gtf,chromsome == 'MT')$gene_name


DefaultAssay(pig) <- 'RNA'
genes <- rownames(pig) %>% as.data.frame()
colnames(genes) <- 'gene'
rownames(genes) <- genes$gene
genes$newgene <- NA
for (i in 1:nrow(genes)){
  if (genes[i, 'gene'] %in% mt_genes){
    genes$newgene[i] <- paste('MT', genes[i, 'gene'], sep = '-')
  }else{
    genes$newgene[i] = genes$gene[i]
  }
}

genes$mt <- grepl("^MT-",genes$newgene)
genes$rp <- grepl("^RP[SL][[:digit:]]",genes$newgene)
table(genes$mt)

pig <- RenameGenesSeurat(pig, genes$newgene)

table(grepl("^RP[SL][[:digit:]]",rownames(pig[["RNA"]])))
table(grepl("^MT-",rownames(pig[["RNA"]])))
pig[["percent.rp"]]  = PercentageFeatureSet(pig, pattern = "^RP[SL][[:digit:]]")
pig[["percent.mt"]]  = PercentageFeatureSet(pig, pattern = "^MT-")
dim(pig)
summary(pig[["nCount_RNA"]]  )
summary(pig[["nFeature_RNA"]])
summary(pig[["percent.mt"]]  )
summary(pig[["percent.rp"]]  )
VlnPlot(pig,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(pig, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(pig, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(pig, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 + plot3
# filter
# low quality cell
pig$quality <- ifelse(pig$percent.mt > 20 | pig$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp & mt
genes <- filter(genes,mt=='FALSE') %>% filter(rp=='FALSE')

pig <- pig[genes$newgene,]
# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(pig))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(pig))
pig <- NormalizeData(pig)
pig <- CellCycleScoring(pig,g2m.features = g2m_genes,s.features = s_genes)
table(pig$Phase)

pig <- pig %>%
  SCTransform(method = 'glmGamPoi',vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

# determine number of PCs to keep based on variance explained
pig <- pig %>% RunHarmony(group.by.vars = 'sample')

DimPlot(pig,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_pig.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = pig, reduction = 'harmony')

pig <- pig %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>% # resolution = c(seq(0.1,1,0.1))
  identity()

DimPlot(pig,label = T,repel = T)

# annotation
job::job({
  plan('multisession', workers = 5)
  pig.marker <- FindAllMarkers(pig,)
  openxlsx::write.xlsx(pig.marker,'./file/markers_seurat_clusters_Pig.xlsx')
},title = 'pig_FindMarker')


# annotation
DimPlot(pig,group.by = 'sample')
DimPlot(pig, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Pig.pdf')

# singleR annotation
library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(pig,slot = 'data')
clusters <- pig$seurat_clusters

# auto annotation
pig_singler <- SingleR(test = sct_matrix, ref = immRNA,
                       labels = immRNA$label.main,
                       clusters = clusters,
                       assay.type.test = 'logcounts',
                       assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(pig$seurat_clusters),
                       cell_annotation=pig_singler$labels)

# sctype ###
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")
# assign cell types
scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)
# Select max score annotation
cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- pig$seurat_clusters[match(rownames(cellmax),rownames(pig@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)


table(pig$seurat_clusters)

celltype


pig.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Pig.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Ss.eg.db)
markers <- pig.marker %>% filter(cluster %in% c('26')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Ss.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Ss.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=40)
head(ego,n=10L)


kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'ssc',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Ss.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)


DimPlot(pig, group.by = 'seurat_clusters', label = T, repel = T)
# filter  #erythrocyte 
FeaturePlot(pig,features = 'HBB')
# cluster none
VlnPlot(pig,features = 'CD79B') #B cells
VlnPlot(pig,features = 'CD3E') #T cells
VlnPlot(pig,features = 'GNLY') #NK cells
VlnPlot(pig,features = 'CD14') #Monocytes
VlnPlot(pig,features = 'TUBB1') #Platelets
VlnPlot(pig,features = 'FLT3') #mDCs
VlnPlot(pig,features = 'CSF3R') #Neutrophils
VlnPlot(pig,features = 'HBB') #Erythrocytes


FeaturePlot(pig,features = 'CD79B') #B cells
FeaturePlot(pig,features = 'CD3E') #T cells
FeaturePlot(pig,features = 'GNLY') #NK cells
FeaturePlot(pig,features = 'CD14') #Monocytes
FeaturePlot(pig,features = 'TUBB1') #Platelets
FeaturePlot(pig,features = 'FLT3') #mDCs
FeaturePlot(pig,features = 'FCERIA') #mDCs
FeaturePlot(pig,features = 'LGALS3') #Neutrophils
FeaturePlot(pig,features = 'HBB') #Erythrocytes

FeaturePlot(pig,features = 'LOC106510122')

Idents(pig) <- 'seurat_clusters'

pig <- RenameIdents(pig,
                    '8'='B cells','10'='B cells',
                    '1'='T cells','2'='T cells','3'='T cells','4'='T cells','5'='T cells','7'='T cells','9'='T cells',
                    '11'='T cells','12'='T cells','15'='T cells','16'='T cells','20'='T cells','21'='T cells',
                    '23'='T cells','24'='T cells','25'='T cells',
                    '18'='NK cells',
                    '6'='Monocytes','17'='Monocytes','22'='Monocytes',
                    '13'='Erythrocytes','14'='Erythrocytes','19'='Erythrocytes',
                    '26'='mDCs') 

cluster_letters <- Idents(pig)
names(cluster_letters) <- colnames(pig)
pig <- AddMetaData(object = pig,
                   metadata = cluster_letters,
                   col.name = "cell_annotation")

# add low quality and doublet 
metadata <- pig@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

for(i in 1:nrow(metadata)){
  if(metadata[i,3] == 'Erythrocytes'){
    metadata[i,4] = 'Erythrocytes'
  }
}

cell_id <- colnames(subset(pig[,pig$cell_annotation %in% 'Erythrocytes'], percent.mt > 20))
for(i in 1:length(cell_id)){
  metadata[cell_id[i],4] = 'Low quality cells'
}


pig$multi_annotation <- metadata$multi_annotation


# c4a_gui()
DimPlot(pig, group.by = 'multi_annotation')
ggsave('./plot/cluster/pig_cell_cluster.pdf')

saveRDS(pig,'./rds/Pig_raw.rds')

pig <- readRDS('./rds/Pig_raw.rds')


################## mouse ################
mouse <- merge(sp.list[['Mouse_balbc']],sp.list[['Mouse_c57bl6']],
               add.cell.ids = c('Mouse_balbc','Mouse_c57bl6'),project = 'Mouse')

meta <- mouse@meta.data
meta[,'orig.ident'] = 'Mouse'
meta$sample <- grepl("^Mouse_balbc", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Mouse_balbc', 'Mouse_c57bl')
mouse@meta.data <- meta

# percent.mt
DefaultAssay(mouse) <- 'RNA'
table(grepl("^Rp[sl][[:digit:]]",rownames(mouse[["RNA"]])))
table(grepl("^mt-",rownames(mouse[["RNA"]])))
mouse[["percent.rp"]]  = PercentageFeatureSet(mouse, pattern = "^Rp[sl][[:digit:]]")
mouse[["percent.mt"]]  = PercentageFeatureSet(mouse, pattern = "^mt-")
dim(mouse)
summary(mouse[["nCount_RNA"]]  )
summary(mouse[["nFeature_RNA"]])
summary(mouse[["percent.mt"]]  )
summary(mouse[["percent.rp"]]  )

VlnPlot(mouse,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3

# low quality cell
mouse$quality <- ifelse(mouse$percent.mt > 10 | mouse$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp and mt genes
genes <- rownames(mouse) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^Rp[sl][[:digit:]]",rownames(genes))
genes$mt <- grepl("^mt-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')
# remove rp and mt genes
mouse <- mouse[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(mouse))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(mouse))
mouse <- CellCycleScoring(mouse,g2m.features = g2m_genes,s.features = s_genes)
table(mouse$Phase)

mouse <- mouse %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

mouse <- mouse %>% RunHarmony(group.by.vars = 'sample')

DimPlot(mouse,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_mouse.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = mouse, reduction = 'harmony')

mouse <- mouse %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

DimPlot(mouse,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  mouse.marker <- FindAllMarkers(mouse)
  openxlsx::write.xlsx(mouse.marker,'./file/markers_seurat_clusters_Mouse.xlsx')
},title = 'mouse_FindMarker')

# singleR annotation
library(SingleR)
mouseRNA <- MouseRNAseqData()

sct_matrix <- GetAssayData(mouse,slot = 'data')
clusters <- mouse$seurat_clusters

# auto annotation
mouse_singler <- SingleR(test = sct_matrix, ref = mouseRNA,
                         labels = mouseRNA$label.main,
                         clusters = clusters,
                         assay.type.test = 'logcounts',
                         assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(mouse$seurat_clusters),
                       cell_annotation=mouse_singler$labels)

# sctype ###
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_costume_mouse.xlsx','Peripheral blood')
# assign cell types
scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)
# Select max score annotation
cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- mouse$seurat_clusters[match(rownames(cellmax),rownames(mouse@meta.data))]
table(cellmax$cell_annotation)

# manul
# annotation
DimPlot(mouse,group.by = 'sample')
DimPlot(mouse, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Mouse.pdf')


# filter ########## whether > 10 cells
table(mouse$seurat_clusters)

# filter  # Erythrocytes
VlnPlot(mouse,features = 'Hbb-bs')

VlnPlot(mouse,features = 'Cd79a') #B cells
VlnPlot(mouse,features = 'Cd3d') #T cells
VlnPlot(mouse,features = 'Ncr1') #NK cells
VlnPlot(mouse,features = 'Csf1r') #Monocytes
VlnPlot(mouse,features = 'S100a8') #Neutrophils
VlnPlot(mouse,features = 'Siglech') #Dendritic cells
VlnPlot(mouse,features = 'Cd34') #
VlnPlot(mouse,features = 'Iglv1')
VlnPlot(mouse,features = 'Iglv2')
VlnPlot(mouse,features = 'Iglc1')
VlnPlot(mouse,features = 'Iglc2')
VlnPlot(mouse,features = 'Igkv1-117')
VlnPlot(mouse,features = 'Igkv10-96')
VlnPlot(mouse,features = 'Igkv1-110')
VlnPlot(mouse,features = 'Igkv1-135')
VlnPlot(mouse,features = 'Igkv14-126')
VlnPlot(mouse,features = 'Hspa1a')

VlnPlot(mouse,features = 'Myof')

FeaturePlot(mouse,features = 'Cd79a') #B cells
FeaturePlot(mouse,features = 'Cd3d') #T cells
FeaturePlot(mouse,features = 'Ncr1') #NK cells
FeaturePlot(mouse,features = 'Csf1r') #Monocytes
FeaturePlot(mouse,features = 'S100a8') #Neutrophils
FeaturePlot(mouse,features = 'Siglech') #Dendritic cells
FeaturePlot(mouse,features = 'Cd34') #

FeaturePlot(mouse,features = 'Iglv1')
FeaturePlot(mouse,features = 'Iglv2')
FeaturePlot(mouse,features = 'Ighv9-3')
FeaturePlot(mouse,features = 'Igkv1-117')
FeaturePlot(mouse,features = 'Igkv10-96')
FeaturePlot(mouse,features = 'Igkv1-110')
FeaturePlot(mouse,features = 'Igkv1-135')
FeaturePlot(mouse,features = 'Igkv14-126')
FeaturePlot(mouse,features = 'Hspa1a')
VlnPlot(mouse,features = 'Camk1')

FeaturePlot(mouse,features = 'Siglech')
VlnPlot(mouse,features = 'Flt3')

Idents(mouse) <- 'seurat_clusters'
mouse <- RenameIdents(mouse,
                      '1'='B cells','4'='B cells','5'='B cells','9'='B cells','17'='B cells','24'='B cells',
                      '8'='B cells','10'='B cells','11'='B cells','12'='B cells','13'='B cells','18'='B cells',
                      '20'='B cells','21'='B cells','23'='B cells','27'='B cells','31'='B cells',
                      '2'='T cells','6'='T cells','7'='T cells','25'='T cells','28'='T cells',
                      '19'='NK cells',
                      '15'='Monocytes','22'='Monocytes',
                      '3'='Monocytes','26'='Monocytes','29'='Monocytes','32'='Monocytes',
                      '16'='Neutrophils',
                      '30'='pDCs',
                      '14'='Erythrocytes') 

cluster_letters <- Idents(mouse)
names(cluster_letters) <- colnames(mouse)
mouse <- AddMetaData(object = mouse,
                     metadata = cluster_letters,
                     col.name = "cell_annotation")


# add low quality and doublet 
metadata <- mouse@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

for(i in 1:nrow(metadata)){
  if(metadata[i,3] == 'Erythrocytes'){
    metadata[i,4] = 'Erythrocytes'
  }
}

cell_id <- colnames(subset(mouse[,mouse$cell_annotation %in% 'Erythrocytes'], percent.mt > 10))
for(i in 1:length(cell_id)){
  metadata[cell_id[i],4] = 'Low quality cells'
}


mouse$multi_annotation <- metadata$multi_annotation


# c4a_gui()
DimPlot(mouse, group.by = 'multi_annotation', cols = c4a('redon',11))
ggsave('./plot/cluster/mouse_cell_cluster.pdf')

saveRDS(mouse,'./rds/Mouse_raw.rds')

mouse <- readRDS('./rds/Mouse_raw.rds')

################## rat ################
rat <- merge(sp.list[['Rat1']],sp.list[['Rat2']],
             add.cell.ids = c('Rat1','Rat2'),project = 'Rat')

meta <- rat@meta.data
meta[,'orig.ident'] = 'Rat'
meta$sample <- grepl("^Rat1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Rat1', 'Rat2')
rat@meta.data <- meta


DefaultAssay(rat) <- 'RNA'
table(grepl("^Rp[sl][[:digit:]]",rownames(rat[["RNA"]])))
table(grepl("^Mt-",rownames(rat[["RNA"]])))
rat[["percent.rp"]]  = PercentageFeatureSet(rat, pattern = "^Rp[sl][[:digit:]]")
rat[["percent.mt"]]  = PercentageFeatureSet(rat, pattern = "^Mt-")
dim(rat)
summary(rat[["nCount_RNA"]]  )
summary(rat[["nFeature_RNA"]])
summary(rat[["percent.mt"]]  )
summary(rat[["percent.rp"]]  )
VlnPlot(rat,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(rat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(rat, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(rat, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3

# low quality cell
rat$quality <- ifelse(rat$percent.mt > 10 | rat$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp and mt genes
genes <- rownames(rat) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^Rp[sl][[:digit:]]",rownames(genes))
genes$mt <- grepl("^Mt-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')
## remove rp and mt genes
rat <- rat[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(rat))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(rat))
rat <- CellCycleScoring(rat,g2m.features = g2m_genes,s.features = s_genes)
table(rat$Phase)

rat <- rat %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

rat <- rat %>% RunHarmony(group.by.vars = 'sample')

DimPlot(rat,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_rat.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = rat, reduction = 'harmony')

rat <- rat %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

DimPlot(rat,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  rat.marker <- FindAllMarkers(rat)
  openxlsx::write.xlsx(rat.marker,'./file/markers_seurat_clusters_Rat.xlsx')
},title = 'rat_FindMarker')


# annotation
DimPlot(rat,group.by = 'sample')
DimPlot(rat, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Rat.pdf')



# singleR annotation
library(SingleR)
mouseRNA <- MouseRNAseqData()

sct_matrix <- GetAssayData(rat,slot = 'data')
clusters <- rat$seurat_clusters

# auto annotation
rat_singler <- SingleR(test = sct_matrix, ref = mouseRNA,
                       labels = mouseRNA$label.main,
                       clusters = clusters,
                       assay.type.test = 'logcounts',
                       assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(rat$seurat_clusters),
                       cell_annotation=rat_singler$labels)


lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_costume_mouse.xlsx','Peripheral blood')
# assign cell types
scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- rat$seurat_clusters[match(rownames(cellmax),rownames(rat@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)


table(rat$seurat_clusters)

# filter  # Erythrocytes
VlnPlot(rat,features = 'Hbb-bs')
FeatuPlot(rat,features = 'Hbb-bs')

VlnPlot(rat,features = 'Cd79a') #B cells
VlnPlot(rat,features = 'Cd3d') #T cells
VlnPlot(rat,features = 'Ncr1') #NK cells
VlnPlot(rat,features = 'Csf1r') #Monocytes
VlnPlot(rat,features = 'Siglech') #pDendritic cells
VlnPlot(rat,features = 'Flt3') #mDendritic cells
VlnPlot(rat,features = 'S100a9') #Neutrophils
VlnPlot(rat,features = 'Pf4') #Megakaryocytes
VlnPlot(rat,features = 'Car1') #Hematopoietic stem cells

VlnPlot(rat,features = 'Iglv1')
VlnPlot(rat,features = 'Iglv2')
VlnPlot(rat,features = 'Iglc1')
VlnPlot(rat,features = 'Iglc2')
VlnPlot(rat,features = 'Igkv1-117')
VlnPlot(rat,features = 'Igkv10-96')
VlnPlot(rat,features = 'Igkv1-110')
VlnPlot(rat,features = 'Igkv1-135')
VlnPlot(rat,features = 'Igkv14-126')
VlnPlot(rat,features = 'Hspa1a')



FeaturePlot(rat,features = 'Cd79a') #B cells
FeaturePlot(rat,features = 'Cd3d') #T cells
FeaturePlot(rat,features = 'Ncr1') #NK cells
FeaturePlot(rat,features = 'Fcnb') #Monocytes
FeaturePlot(rat,features = 'Csf1r') #Macrophages
FeaturePlot(rat,features = 'Siglech') #p Dendritic cells
FeaturePlot(rat,features = 'Flt3') #m Dendritic cells
FeaturePlot(rat,features = 'S100a9') #Neutrophils
FeaturePlot(rat,features = 'Gp1bb') #Platelets
# FeaturePlot(rat,features = 'Car1') #Hematopoietic stem cells

FeaturePlot(rat,features = 'Iglv1')
FeaturePlot(rat,features = 'Iglv2')
FeaturePlot(rat,features = 'Iglc1')
FeaturePlot(rat,features = 'Iglc2')
FeaturePlot(rat,features = 'Igkv1-117')
FeaturePlot(rat,features = 'Igkv10-96')
FeaturePlot(rat,features = 'Igkv1-110')
FeaturePlot(rat,features = 'Igkv1-135')
FeaturePlot(rat,features = 'Igkv14-126')
FeaturePlot(rat,features = 'Siglech')
FeaturePlot(rat,features = 'Flt3')

Idents(rat) <- 'seurat_clusters'
rat <- RenameIdents(rat,
                    '4'='B cells','17'='B cells','19'='B cells','20'='B cells',
                    '3'='T cells','8'='T cells','10'='T cells','13'='T cells','15'='T cells','16'='T cells','23'='T cells',
                    '2'='NK cells','6'='NK cells','11'='NK cells','21'='NK cells',
                    '1'='Monocytes','7'='Monocytes','12'='Monocytes',
                    '5'='Monocytes','9'='Monocytes',
                    '22'='mDCs',
                    '14'='pDCs',
                    '24'='Neutrophils',
                    '25'='Platelets',
                    '18'='Erythrocytes') 

cluster_letters <- Idents(rat)
names(cluster_letters) <- colnames(rat)
rat <- AddMetaData(object = rat,
                   metadata = cluster_letters,
                   col.name = "cell_annotation")


metadata <- rat@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '26'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '27'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

for(i in 1:nrow(metadata)){
  if(metadata[i,3] == 'Erythrocytes'){
    metadata[i,4] = 'Erythrocytes'
  }
}


rat$multi_annotation <- metadata$multi_annotation


# c4a_gui()
DimPlot(rat, group.by = 'multi_annotation')
ggsave('./plot/cluster/rat_cell_cluster.pdf')

saveRDS(rat,'./rds/Rat_raw.rds')

rat <- readRDS('./rds/Rat_raw.rds')

################## monkey ##################
monkey <- merge(sp.list[['Monkey1']], y = c(sp.list[['Monkey2']], sp.list[['Monkey3']], sp.list[['Monkey4']]),
             add.cell.ids = c('Monkey1','Monkey2','Monkey3','Monkey4'),project = 'Monkey')

meta <- monkey@meta.data
meta[,'orig.ident'] = 'Monkey'
head(str_extract(rownames(meta),'Monkey\\d+'))
meta$sample <- str_extract(rownames(meta),'Monkey\\d+')
monkey@meta.data <- meta

monkey_gtf <- read.table("./ncbi/genome/GCF_003339765.1.gtf/ncbi_dataset/data/GCF_003339765.1/genomic.gtf", sep="\t", header=FALSE, stringsAsFactors=FALSE)

monkey_gtf <- cbind(monkey_gtf[1:8], t(sapply(strsplit(monkey_gtf$V9, "; "), function(x) {
  sapply(c("gene_id", "gene_biotype"), function(y) {
    tmp <- grep(y, x, fixed=TRUE)
    if (length(tmp) == 0) {
      return("")
    } else {
      return(sub(paste0("^", y, " "), "", x[tmp]))
    }
  })
})))

colnames(monkey_gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame",
                       "gene_id", "gene_biotype")

monkey_gtf <- filter(monkey_gtf, feature == 'gene') %>% filter(gene_biotype == 'protein_coding')

mt_genes <- filter(monkey_gtf, seqname == 'NC_005943.1')$gene_id

## percent.mt
DefaultAssay(monkey) <- 'RNA'
genes <- rownames(monkey) %>% as.data.frame()
colnames(genes) <- 'gene'
rownames(genes) <- genes$gene
genes$newgene <- genes$gene

genes$mt <- grepl("^KEG06-",genes$newgene)
genes$rp <- grepl("^RP[SL][[:digit:]]",genes$newgene)
table(genes$mt)
table(genes$rp)

# monkey <- RenameGenesSeurat(monkey, genes$newgene)

table(grepl("^RP[SL][[:digit:]]",rownames(monkey[["RNA"]])))
table(grepl("^KEG06-",rownames(monkey[["RNA"]])))
monkey[["percent.rp"]]  = PercentageFeatureSet(monkey, pattern = "^RP[SL][[:digit:]]")
monkey[["percent.mt"]]  = PercentageFeatureSet(monkey, pattern = "^KEG06-")
dim(monkey)
summary(monkey[["nCount_RNA"]]  )
summary(monkey[["nFeature_RNA"]])
summary(monkey[["percent.mt"]]  )
summary(monkey[["percent.rp"]]  )

VlnPlot(monkey,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(monkey, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(monkey, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(monkey, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 +plot3

# low quality cell
monkey$quality <- ifelse(monkey$nFeature_RNA < 300, "Low quality", "High quality")

# filter rp and mt genes
genes <- genes %>% filter(rp == 'FALSE') %>% filter(mt == 'FALSE')

## remove rp and mt genes
monkey <- monkey[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(monkey))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(monkey))
monkey <- NormalizeData(monkey) #
monkey <- CellCycleScoring(monkey,g2m.features = g2m_genes,s.features = s_genes)
table(monkey$Phase)

monkey <- monkey %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()


monkey <- monkey %>% RunHarmony(group.by.vars = 'sample')

DimPlot(monkey, reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_monkey.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = monkey, reduction = 'harmony')

monkey <- monkey %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>% # resolution = c(seq(0.1,1,0.1))
  identity()


DimPlot(monkey,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  monkey.marker <- FindAllMarkers(monkey)
  openxlsx::write.xlsx(monkey.marker,'./file/markers_seurat_clusters_Monkey.xlsx')
},title = 'monkey_FindMarker')


# annotation
DimPlot(monkey, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Monkey.pdf')

# singleR annotation
library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(monkey,slot = 'data')
clusters <- monkey$seurat_clusters
genes <- rownames(sct_matrix) %>% as.data.frame()
colnames(genes) <- 'gene'
monkey_human <- read.csv('./protein/macaque_human.csv')
genes$human <- monkey_human$Human.gene.name[match(genes$gene,monkey_human$Macaque.gene.name)]

genes <- na.omit(genes)
sct_matrix <- sct_matrix[genes$gene,]
rownames(sct_matrix) <- genes$human


# auto annotation
monkey_singler <- SingleR(test = sct_matrix, ref = immRNA,
                          labels = immRNA$label.main,
                          clusters = clusters,
                          assay.type.test = 'logcounts',
                          assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(monkey$seurat_clusters),
                       cell_annotation=monkey_singler$labels)

# sctype ###
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")

scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)


View(es.max)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- monkey$seurat_clusters[match(rownames(cellmax),rownames(monkey@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

celltype


table(monkey$seurat_clusters)

## cluster none
VlnPlot(monkey,features = 'HBB') #Erythrocytes
VlnPlot(monkey,features = 'CD79B') #B cells
VlnPlot(monkey,features = 'CD3E') #T cells
VlnPlot(monkey,features = 'GNLY') #NK cells
VlnPlot(monkey,features = 'MAFB') #Monocytes
VlnPlot(monkey,features = 'TUBB1') #Platelets
VlnPlot(monkey,features = 'FCER1A') #mDCs

FeaturePlot(monkey,features = 'CD79B') #B cells
FeaturePlot(monkey,features = 'CD3E') #T cells
FeaturePlot(monkey,features = 'GNLY') #NK cells
FeaturePlot(monkey,features = 'MAFB') #Monocytes
FeaturePlot(monkey,features = 'TUBB1') #Platelets
FeaturePlot(monkey,features = 'FLT3') #mDCs
FeaturePlot(monkey,features = 'HBB') #Erythrocytes

FeaturePlot(monkey,features = 'FCER1A')

Idents(monkey) <- 'seurat_clusters'
monkey <- RenameIdents(monkey,
                       '1'='B cells','3'='B cells','4'='B cells','5'='B cells','7'='B cells','8'='B cells','10'='B cells',
                       '13'='B cells','15'='B cells','21'='B cells',
                       '2'='T cells','14'='T cells','8'='T cells','17'='T cells','24'='T cells',
                       '6'='NK cells','9'='NK cells','12'='NK cells','16'='NK cells','20'='NK cells',
                       '11'='Monocytes','19'='Monocytes',
                       '18'='pDCs',
                       '23'='Platelets',
                       '22'='Erythrocytes') 

cluster_letters <- Idents(monkey)
names(cluster_letters) <- colnames(monkey)
monkey <- AddMetaData(object = monkey,
                      metadata = cluster_letters,
                      col.name = "cell_annotation")

# add low quality and doublet 
metadata <- monkey@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '22'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

for (i in 1:nrow(metadata)){
  if(metadata[i,3] == 'Monocytes'){
    metadata[i,4] = 'Monocytes'
  }else if(metadata[i,3] == 'Erythrocytes'){
    metadata[i,4] = 'Erythrocytes'
  }
}

monkey$multi_annotation <- metadata$multi_annotation
table(monkey$multi_annotation)

# c4a_gui()
DimPlot(monkey, group.by = 'multi_annotation', cols = c4a('redon',12))
ggsave('./plot/cluster/monkey_cell_cluster.pdf')

saveRDS(monkey,'./rds/Monkey_raw.rds')

monkey <- readRDS('./rds/Monkey_raw.rds')

################## chimpanzee ##############
chimpanzee <- merge(sp.list[['Chimpanzee1']], y = c(sp.list[['Chimpanzee2']], sp.list[['Chimpanzee3']], sp.list[['Chimpanzee4']]),
                add.cell.ids = c('Chimpanzee1','Chimpanzee2','Chimpanzee3','Chimpanzee4'),project = 'Chimpanzee')

meta <- chimpanzee@meta.data
meta[,'orig.ident'] = 'Chimpanzee'
head(str_extract(rownames(meta),'Chimpanzee\\d+'))
meta$sample <- str_extract(rownames(meta),'Chimpanzee\\d+')
chimpanzee@meta.data <- meta

DefaultAssay(chimpanzee) <- 'RNA'
table(grepl("^RP[SL][[:digit:]]",rownames(chimpanzee[["RNA"]])))
table(grepl("^MT-",rownames(chimpanzee[["RNA"]])))
chimpanzee[["percent.mt"]]  = PercentageFeatureSet(chimpanzee, pattern = "^MT-")
chimpanzee[["percent.rp"]]  = PercentageFeatureSet(chimpanzee, pattern = "^RP[SL][[:digit:]]")
dim(chimpanzee)
summary(chimpanzee[["nCount_RNA"]]  )
summary(chimpanzee[["nFeature_RNA"]])
summary(chimpanzee[["percent.rp"]])

VlnPlot(chimpanzee,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.rp"), ncol = 3)

plot1 <- FeatureScatter(chimpanzee, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(chimpanzee, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(chimpanzee, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 + plot3


chimpanzee$quality <- ifelse(chimpanzee$nFeature_RNA < 300, "Low quality", "High quality")


genes <- rownames(chimpanzee) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^RP[SL][[:digit:]]",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE')

## remove rp and mt genes
chimpanzee <- chimpanzee[genes$gene,]

# normalize and cluster
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(chimpanzee))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(chimpanzee))
chimpanzee <- NormalizeData(chimpanzee) #
chimpanzee <- CellCycleScoring(chimpanzee,g2m.features = g2m_genes,s.features = s_genes)
table(chimpanzee$Phase)

chimpanzee <- chimpanzee %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()


# determine number of PCs to keep based on variance explained
chimpanzee <- chimpanzee %>% RunHarmony(group.by.vars = 'sample')

DimPlot(chimpanzee, reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_chimpanzee.pdf', height = 5.64, width = 7)

# determine number of PCs to keep based on variance explained
n_pcs <- npcs(object = chimpanzee, reduction = 'harmony')

chimpanzee <- chimpanzee %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>% # resolution = c(seq(0.1,1,0.1))
  identity()


DimPlot(chimpanzee,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  chimpanzee.marker <- FindAllMarkers(chimpanzee)
  openxlsx::write.xlsx(chimpanzee.marker,'./file/markers_seurat_clusters_Chimpanzee.xlsx')
},title = 'chimpanzee_FindMarker')


# annotation
DimPlot(chimpanzee, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Chimpanzee.pdf')

# singleR annotation
library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(chimpanzee,slot = 'data')
clusters <- chimpanzee$seurat_clusters
genes <- rownames(sct_matrix) %>% as.data.frame()
colnames(genes) <- 'gene'
chimpanzee_human <- read.csv('./protein/chimpanzee_human.csv')
genes$human <- chimpanzee_human$Human.gene.name[match(genes$gene,chimpanzee_human$Chimpanzee.gene.name)]

genes <- na.omit(genes)
sct_matrix <- sct_matrix[genes$gene,]
rownames(sct_matrix) <- genes$human


# auto annotation
chimpanzee_singler <- SingleR(test = sct_matrix, ref = immRNA,
                              labels = immRNA$label.main,
                              clusters = clusters,
                              assay.type.test = 'logcounts',
                              assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(chimpanzee$seurat_clusters),
                       cell_annotation=chimpanzee_singler$labels)

### sctype ###
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")
# assign cell types
scRNAseqData = sct_matrix %>% as.matrix() #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- chimpanzee$seurat_clusters[match(rownames(cellmax),rownames(chimpanzee@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

celltype

table(chimpanzee$seurat_clusters)

VlnPlot(chimpanzee,features = 'HBB') #Erythrocytes
VlnPlot(chimpanzee,features = 'CD79A') #B cells
VlnPlot(chimpanzee,features = 'JCHAIN') #B cells
VlnPlot(chimpanzee,features = 'CD3E') #T cells
VlnPlot(chimpanzee,features = 'CD8A') #T cells
VlnPlot(chimpanzee,features = 'CD8B') #T cells
VlnPlot(chimpanzee,features = 'GNLY') # NK cells
VlnPlot(chimpanzee,features = 'KLRC1') #NK cells
VlnPlot(chimpanzee,features = 'MAFB') #Monocytes
VlnPlot(chimpanzee,features = 'CLEC4C') #pDCs

FeaturePlot(chimpanzee,features = 'CD79B') #B cells
FeaturePlot(chimpanzee,features = 'JCHAIN') #B cells
FeaturePlot(chimpanzee,features = 'CD3E') #T cells
FeaturePlot(chimpanzee,features = 'GNLY') #NK cells
FeaturePlot(chimpanzee,features = 'MS4A7') #Monocytes
FeaturePlot(chimpanzee,features = 'CD68') #Macrophages
FeaturePlot(chimpanzee,features = 'CLEC4C') #pDCs
FeaturePlot(chimpanzee,features = 'HBB') #Erythrocytes
FeaturePlot(chimpanzee,features = 'FCER1A') #mDCs
FeaturePlot(chimpanzee,features = 'CD34') #Circulating angiogenic cells

FeaturePlot(chimpanzee,features = 'MAFB')
FeaturePlot(chimpanzee,features = 'CD8B') ### T cells
FeaturePlot(chimpanzee,features = 'CD33')

job::job({
  plan('multisession', workers = 5)
  markers <- FindMarkers(chimpanzee,ident.1 = 28)
},title = 'Marker')


Idents(chimpanzee) <- 'seurat_clusters'
chimpanzee <- RenameIdents(chimpanzee,
                           '7'='B cells','13'='B cells','16'='B cells','18'='B cells','22'='B cells',
                           '3'='NK cells','12'='NK cells','17'='NK cells','15'='NK cells',
                           '1'='T cells','2'='T cells','4'='T cells','6'='T cells','9'='T cells',
                           '10'='T cells','11'='T cells','14'='T cells',
                           '5'='Monocytes','8'='Monocytes','19'='Monocytes',
                           '20'='mDCs') 

cluster_letters <- Idents(chimpanzee)
names(cluster_letters) <- colnames(chimpanzee)
chimpanzee <- AddMetaData(object = chimpanzee,
                          metadata = cluster_letters,
                          col.name = "cell_annotation")

metadata <- chimpanzee@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else if(metadata[i,3] == '21'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

chimpanzee$multi_annotation <- metadata$multi_annotation
table(chimpanzee$multi_annotation)

#c4a_gui()
DimPlot(chimpanzee, group.by = 'multi_annotation', cols = c4a('redon',12))
ggsave('./plot/cluster/chimpanzee_cell_cluster.pdf')

saveRDS(chimpanzee,'./rds/Chimpanzee_raw.rds')

chimpanzee <- readRDS('./rds/Chimpanzee_raw.rds')

################## human ##################
human <- merge(sp.list[['Human1']],sp.list[['Human2']],
               add.cell.ids = c('Human1','Human2'),project = 'Human')

meta <- human@meta.data
meta[,'orig.ident'] = 'Human'
meta$sample <- grepl("^Human1", rownames(meta))
meta$sample <- ifelse(meta$sample == 'TRUE', 'Human1', 'Human2')
human@meta.data <- meta


DefaultAssay(human) <- 'RNA'
table(grepl("^RP[SL][[:digit:]]",rownames(human[["RNA"]])))
table(grepl("^MT-",rownames(human[["RNA"]])))
human[["percent.rp"]]  = PercentageFeatureSet(human, pattern = "^RP[SL][[:digit:]]")
human[["percent.mt"]]  = PercentageFeatureSet(human, pattern = "^MT-")
dim(human)
summary(human[["nCount_RNA"]]  )
summary(human[["nFeature_RNA"]])
summary(human[["percent.mt"]]  )
summary(human[["percent.rp"]]  )

VlnPlot(human,group.by = 'sample',features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rp"), ncol = 4)

plot1 <- FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = 'sample')+NoLegend()
plot2 <- FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = 'sample')+NoLegend()
plot3 <- FeatureScatter(human, feature1 = "nCount_RNA", feature2 = "percent.rp",group.by = 'sample')+NoLegend()
plot1 + plot2 + plot3

human$quality <- ifelse(human$percent.mt > 20 | human$nFeature_RNA < 300, "Low quality", "High quality")


genes <- rownames(human) %>% as.data.frame()
colnames(genes)[1] <- 'gene'
rownames(genes) <- genes$gene
genes$rp <- grepl("^RP[SL][[:digit:]]",rownames(genes))
genes$mt <- grepl("^MT-",rownames(genes))

genes <- genes %>%
  filter(rp == 'FALSE') %>%
  filter(mt == 'FALSE')

human <- human[genes$gene,]

g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes,match = rownames(human))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes,match = rownames(human))
human <- NormalizeData(human) #
human <- CellCycleScoring(human,g2m.features = g2m_genes,s.features = s_genes)
table(human$Phase)

human <- human %>%
  SCTransform(method = 'glmGamPoi', vars.to.regress = c('G2M.Score','S.Score')) %>%
  RunPCA()

human <- human %>% RunHarmony(group.by.vars = 'sample')

DimPlot(human,reduction = 'harmony',group.by = 'sample')
ggsave('./plot/quality/harmony_human.pdf', height = 5.64, width = 7)

n_pcs <- npcs(object = human, reduction = 'harmony')

human <- human %>%
  RunUMAP(reduction = "harmony", dims = 1:n_pcs) %>%
  FindNeighbors(reduction = "harmony", dims = 1:n_pcs) %>%
  FindClusters(algorithm = 4) %>%
  identity()

DimPlot(human,label = T,repel = T)

job::job({
  plan('multisession', workers = 4)
  human.marker <- FindAllMarkers(human)
  openxlsx::write.xlsx(human.marker,'./file/markers_seurat_clusters_Human.xlsx')
},title = 'human_FindMarker')


DimPlot(human,group.by = 'sample')
DimPlot(human, label = T, repel = T, group.by = 'seurat_clusters')
ggsave('./plot/cluster/seurat_clusters_Human.pdf')


library(SingleR)
library(celldex)
immRNA <- DatabaseImmuneCellExpressionData()

sct_matrix <- GetAssayData(human,slot = 'data')
clusters <- human$seurat_clusters

human_singler <- SingleR(test = sct_matrix, ref = immRNA,
                         labels = immRNA$label.main,
                         clusters = clusters,
                         assay.type.test = 'logcounts',
                         assay.type.ref = 'logcounts')

celltype <- data.frame(seurat_clusters=levels(human$seurat_clusters),
                       cell_annotation=human_singler$labels)


lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source('../../download/sc-type-1.0/R/gene_sets_prepare.R')
source('../../download/sc-type-1.0/R/sctype_score_.R')
gs_list <- gene_sets_prepare('../../download/sc-type-1.0/ScTypeDB_short.xlsx',"Immune system")

scRNAseqData = sct_matrix %>% as.matrix() 
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive)

View(es.max)

cellmax <- apply(es.max, 2, function(x) rownames(es.max)[which.max(x)]) %>% as.data.frame()
colnames(cellmax) <- 'cell_annotation'
cellmax$seurat_clusters <- human$seurat_clusters[match(rownames(cellmax),rownames(human@meta.data))]
table(cellmax$cell_annotation,cellmax$seurat_clusters)

table(human$seurat_clusters)

human.marker <- openxlsx::read.xlsx('./file/markers_seurat_clusters_Human.xlsx', sheet = 'Sheet 1')

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

markers <- human.marker %>% filter(cluster %in% c('22')) %>% filter(avg_log2FC > 0.25)

genes <- bitr(markers$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')

ego <- enrichGO(gene          = genes$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego@result$Description,n=20)
head(ego,n=10L)


kk <- enrichKEGG(gene         = genes$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk <- setReadable(kk,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID')
head(kk@result$Description,n=30)
head(kk,n=12)



FeaturePlot(human,features = 'HBB')

VlnPlot(human,features = 'CD79B') #B cells
VlnPlot(human,features = 'CD27') #T cells
VlnPlot(human,features = 'GNLY') #NK cells
VlnPlot(human,features = 'S100A8') #Monocytes
VlnPlot(human,features = 'PPBP') #Megakaryocyte
VlnPlot(human,features = 'JCHAIN')#Plasma cells
VlnPlot(human,features = 'FCER1A')#mDCs
VlnPlot(human,features = 'CLEC4C')#pCDs
VlnPlot(human,features = 'MS4A7')#Macrophages


FeaturePlot(human,features = 'CD79B') #B cells
FeaturePlot(human,features = 'CD27') #T cells
FeaturePlot(human,features = 'GNLY') #NK cells
FeaturePlot(human,features = 'S100A8') #Monocytes
FeaturePlot(human,features = 'PPBP') #Megakaryocyte
FeaturePlot(human,features = 'JCHAIN')#Plasma cells
FeaturePlot(human,features = 'FCER1A')#mDCs
FeaturePlot(human,features = 'CLEC4C')#pCDs
FeaturePlot(human,features = 'MS4A7')#Macrophages

FeaturePlot(human,features = 'TCF7L2')
FeaturePlot(human,features = 'C1QB')

celltype

Idents(human) <- 'seurat_clusters'

human <- RenameIdents(human,
                      '7'='B cells','11'='B cells',
                      '1'='T cells','2'='T cells','4'='T cells','8'='T cells','10'='T cells','16'='T cells','17'='T cells','23'='T cells',
                      '18'='B cells',
                      '5'='NK cells','9'='NK cells','12'='NK cells','13'='NK cells',
                      '3'='Monocytes','6'='Monocytes','14'='Monocytes',
                      '19'='Monocytes',
                      '20'='mDCs',
                      '15'='Platelets','21'='Platelets',
                      '22'='pDCs')

cluster_letters <- Idents(human)
names(cluster_letters) <- colnames(human)
human <- AddMetaData(object = human,
                     metadata = cluster_letters,
                     col.name = "cell_annotation")

############ add low quality and doublet 
metadata <- human@meta.data
metadata$cell_annotation <- as.character(metadata$cell_annotation)
metadata <- metadata[,c('DF_hi.lo','quality','cell_annotation')]
metadata$multi_annotation <- NA
for (i in 1:nrow(metadata)){
  if(metadata[i,1] == 'Doublet_hi'){
    metadata[i,4] = 'Doublets'
  }else if(metadata[i,2] == 'Low quality'){
    metadata[i,4] = 'Low quality cells'
  }else{
    metadata[i,4] = metadata[i,3]
  }
}

human$multi_annotation <- metadata$multi_annotation


#c4a_gui()
DimPlot(human, group.by = 'multi_annotation', cols = c4a('redon',12))
ggsave('./plot/cluster/human_cell_cluster.pdf')

saveRDS(human,'./rds/Human_raw.rds')

human <- readRDS('./rds/Human_raw.rds')
