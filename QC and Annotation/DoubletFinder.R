lapply(c("dplyr","Seurat"), library, character.only = T)
library(ggplot2)
library(patchwork)
library(stringr)
library(DoubletFinder)
library(future)
library(cols4all)
library(dplyr)
library(tidyr)
setwd('~/PBMC/')

plan('multisession', workers = 10)
plan()
# Calculate the number of PCs that contain some proportion (95%) of the variance
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

species <- c('Catfish1','Catfish2',
             'Jacopever1','Jacopever2',
             'Turtle1','Turtle2',
             'Chicken1','Chicken2',
             'Bat1','Bat2','Bat3','Bat4',
             "Cattle1",'Cattle2',
             'Pig1','Pig2',
             'Mouse1','Mouse2',
             'Rat1','Rat2',
             'Monkey1','Monkey2','Monkey3','Monkey4',
             'Chimpanzee1','Chimpanzee2','Chimpanzee3','Chimpanzee4',
             'Human1','Human2')

sp.list <- list(
Yellow_catfish1 = Read10X('./BMK/Pelteobagrus fulvidraco sample1/filtered_feature_bc_matrix/'),
Yellow_catfish2 = Read10X('./BMK/Pelteobagrus fulvidraco sample2/filtered_feature_bc_matrix/'),

Jacopever1 = Read10X('./BMK/Sebastes schlegelii sample1/filtered_feature_bc_matrix/'),
Jacopever2 = Read10X('./BMK/Sebastes schlegelii sample2/filtered_feature_bc_matrix/'),

Turtle1 = Read10X('./BMK/Pelodiscus sinensis sample1/filtered_feature_bc_matrix/'),
Turtle2 = Read10X('./BMK/Pelodiscus sinensis sample2/filtered_feature_bc_matrix/'),

Chicken1 = Read10X('./BMK/Gallus gallus domesticus sample1/filtered_feature_bc_matrix/'),
Chicken2 = Read10X('./BMK/Gallus gallus domesticus sample2/filtered_feature_bc_matrix/'),

Bat1 = Read10X_h5('./monkey/GSM6736397_200221_Dr_Sato_RA_mock_filtered_feature_bc_matrix.h5'),
Bat2 = Read10X_h5('./monkey/GSM6736398_3_koumori_filtered_feature_bc_matrix.h5'),
Bat3 = Read10X_h5('./monkey/GSM6736399_9_RA_SEB_filtered_feature_bc_matrix.h5'),
Bat4 = Read10X_h5('./monkey/GSM6736400_6_200408_Dr_Sato_RA_LPS_filtered_feature_bc_matrix.h5'),

Cattle1 = Read10X('./cattle/GSM5066754/'),
Cattle2 = Read10X('./cattle/GSM5066755/'),

Pig1 = Read10X('./pig/GSM5005383/'),
Pig2 = Read10X('./pig/GSM5005384/'),

Mouse1 = Read10X('./10x/balbc/filtered_feature_bc_matrix/'),
Mouse2 = Read10X('./10x/c57bl6/filtered_feature_bc_matrix/'),

Rat1 = Read10X('./BMK/Rattus norvegicus sample1/filtered_feature_bc_matrix/'),
Rat2 = Read10X('./BMK/Rattus norvegicus sample2/filtered_feature_bc_matrix/'),

Monkey1 = Read10X_h5('./monkey/GSM6736393_3_akagezaru_191004_1120_PBL_filtered_feature_bc_matrix.h5'),
Monkey2 = Read10X_h5('./monkey/GSM6736394_2_akage_filtered_feature_bc_matrix.h5'),
Monkey3 = Read10X_h5('./monkey/GSM6736395_8_MM_SeV_filtered_feature_bc_matrix.h5'),
Monkey4 = Read10X_h5('./monkey/GSM6736396_4_200408_Dr_Sato_MM_LPS_filtered_feature_bc_matrix.h5'),

Chimpanzee1 = Read10X_h5('./monkey/GSM6736389_2_chimpanzee_190724_paru_PBL_filtered_feature_bc_matrix.h5'),
Chimpanzee2 = Read10X_h5('./monkey/GSM6736390_4_chimp_HSV_filtered_feature_bc_matrix.h5'),
Chimpanzee3 = Read10X_h5('./monkey/GSM6736391_7_PT_SeV_filtered_feature_bc_matrix.h5'),
Chimpanzee4 = Read10X_h5('./monkey/GSM6736392_6_chimp_LPS_filtered_feature_bc_matrix.h5'),

Human1 = Read10X('./BMK/Homo sapiens sample1/filtered_feature_bc_matrix/'),
Human2 = Read10X('./BMK/Homo sapiens sample2/filtered_feature_bc_matrix/')
)

sp.list[['Bat3']] = Read10X_h5('./monkey/GSM6736399_9_RA_SEB_filtered_feature_bc_matrix.h5')
sp.list[['Bat4']] =  Read10X_h5('./monkey/GSM6736400_6_200408_Dr_Sato_RA_LPS_filtered_feature_bc_matrix.h5')
sp.list[['Monkey3']] = Read10X_h5('./monkey/GSM6736395_8_MM_SeV_filtered_feature_bc_matrix.h5')
sp.list[['Monkey4']] = Read10X_h5('./monkey/GSM6736396_4_200408_Dr_Sato_MM_LPS_filtered_feature_bc_matrix.h5')
sp.list[['Chimpanzee3']] = Read10X_h5('./monkey/GSM6736391_7_PT_SeV_filtered_feature_bc_matrix.h5')
sp.list[['Chimpanzee4']] = Read10X_h5('./monkey/GSM6736392_6_chimp_LPS_filtered_feature_bc_matrix.h5')

for (i in names(sp.list)){
  sp.list[[i]] <- CreateSeuratObject(sp.list[[i]])
}


for (i in names(sp.list)){
  print(i)
  sp.list[[i]] <- sp.list[[i]] %>%
    SCTransform(method = 'glmGamPoi') %>%
    RunPCA(verbose = F)
  
  # determine number of PCs to keep based on variance explained
  n_pcs <- npcs(object = sp.list[[i]], reduction = 'pca')
  
  sp.list[[i]] <- sp.list[[i]] %>%
    RunUMAP(reduction = "pca", dims = 1:n_pcs) %>%
    FindNeighbors(reduction = "pca", dims = 1:n_pcs) %>%
    FindClusters(algorithm = 4) %>%
    identity()
}

for(i in names(sp.list)){
  print(i)
  n_pcs <- npcs(object = sp.list[[i]], reduction = 'pca')
  sweep.res.list <- paramSweep_v3(sp.list[[i]], PCs = 1:n_pcs, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn_nc <- find.pK(sweep.stats)
  
  pK_bcmvn <- bcmvn_nc$pK[which.max(bcmvn_nc$BCmetric)] %>% as.character() %>% as.numeric()
  ## Homotypic Doublet Proportion Estimate ----------------------------------------------
  annotations <- sp.list[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sp.list[[i]]@meta.data$ClusteringResults
  #nExp_poi <- round(DoubletRate*ncol(sp.list[[i]]))
  nExp_poi <- round(0.075*nrow(sp.list[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying csp.list[[i]]ssification stringencies --------------------------
  sp.list[[i]] <- doubletFinder_v3(sp.list[[i]], PCs = 1:n_pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  sp.list[[i]] <- doubletFinder_v3(sp.list[[i]], PCs = 1:n_pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  
  gc()
}

for(i in names(sp.list)){
  print(i)
  num_meta <- ncol(sp.list[[i]]@meta.data)
  colnames(sp.list[[i]]@meta.data)[(num_meta-3):num_meta] <- c('pANN_nExp_poi','DF.classifications_nExp_poi',
                                            'pANN_nExp_poi.adj','DF.classifications_nExp_poi.adj')
}

for(i in names(sp.list)){
  print(i)
  sp.list[[i]]@meta.data[,"DF_hi.lo"] <- sp.list[[i]]@meta.data$DF.classifications_nExp_poi
  sp.list[[i]]@meta.data$DF_hi.lo[which(sp.list[[i]]@meta.data$DF_hi.lo == "Doublet" & sp.list[[i]]@meta.data$DF.classifications_nExp_poi.adj == "Singlet")] <- "Doublet_lo"
  sp.list[[i]]@meta.data$DF_hi.lo[which(sp.list[[i]]@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
}
table(sp.list[[1]]@meta.data$DF_hi.lo)
DimPlot(sp.list[[1]],group.by = "DF_hi.lo",cols = c("red","gold","black"))

saveRDS(sp.list,'./rds/sp.list.rds')


