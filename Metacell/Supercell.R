library(SuperCell)
library(Seurat)
library(patchwork)

setwd('~/data/PBMC/')

sp.list <- readRDS('./rds/sp.list.rds')

species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig','Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')

MC.list <- list()
MC.seu.list <- list()
gamma <- 20 # Graining level

for (i in species){
  n.pcs <- npcs(sp.list[[i]])
  hvg <- VariableFeatures(sp.list[[i]])
  # Compute metacells using SuperCell package
  MC.list[[i]] <- SCimplify(
    X = GetAssayData(sp.list[[i]]), # single-cell log-normalized gene expression data
    gamma = gamma,
    genes.use = hvg,
    n.pc = n.pcs,
    sc.cell.annotation. = sp.list[[i]]$multi_annotation
  )
  # Compute gene expression of metacells by simply averaging gene expression within each metacell
  MC.ge <- supercell_GE(
    ge = GetAssayData(sp.list[[i]]),
    groups = MC.list[[i]]$membership
  )
  # Alternatively, counts can be averaged (summed up) followed by a lognormalization step (this approach is used in the MetaCell and SEACell algorithms)
  if(0){
    MC.counts <- supercell_GE(
      ge = GetAssayData(sp.list[[i]], slot = "counts"),
      mode = "sum", # summing counts instead of the default averaging
      groups = MC$membership
    )
    
    MC.ge <- Seurat::LogNormalize(MC.counts, verbose = FALSE)
  }
  
  # Annotate metacells to cells line
  MC.list[[i]]$cell_line <- supercell_assign(
    cluster = sp.list[[i]]$multi_annotation,          # single-cell assignment to cell lines 
    supercell_membership = MC.list[[i]]$membership,  # single-cell assignment to metacells
    method = "absolute" # available methods are c("jaccard", "relative", "absolute"), function's help() for explanation
  )
  
  # Compute purity of metacells as :
  #  * a proportion of the most abundant cell type withing metacells (`method = `"max_proportion)
  #  * an entropy of cell type within metacells (`method = "entropy"`)
  method_purity <- c("max_proportion", "entropy")[1]
  MC.list[[i]]$purity <- supercell_purity(
    clusters = sp.list[[i]]$multi_annotation,
    supercell_membership = MC.list[[i]]$membership, 
    method = method_purity
  )
  
  # Metacell purity distribution
  summary(MC.list[[i]]$purity)
  
  ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  ##       1       1       1       1       1       1
  #### save to picture 
  filename = paste('./plot/supercell/','purity_of_metacells_',names(MC.list[i]),'.pdf',sep = '')
  pdf(filename)
  hist(MC.list[[i]]$purity, main = paste0("Purity of metacells \nin terms of cell line composition (", method_purity,")"))
  dev.off()
  #### covert to seurat object 
  MC.seu.list[[i]] <- supercell_2_Seurat(
    SC.GE = MC.ge, 
    SC = MC.list[[i]], 
    fields = c("cell_line", "purity"),
    var.genes = MC.list[[i]]$genes.use,
    N.comp = 10
  )
  
  Idents(MC.seu.list[[i]]) <- "cell_line"
  
  MC.seu.list[[i]] <- RunUMAP(MC.seu.list[[i]], dims = 1:npcs(MC.seu.list[[i]]))
}

for(i in names(MC.list)){
  filename = paste('./plot/supercell/','supercell_',i,'.pdf',sep = '')
  pdf(filename)
  supercell_plot(
    MC.list[[i]]$graph.supercells, 
    group = MC.list[[i]]$cell_line, 
    color.use = cell.colors1[levels(as.factor(MC.list[[i]]$cell_line))],
    seed = 1, 
    alpha = -pi/2,
    main  = "Metacells colored by cell line assignment"
  )
  dev.off()
}

for(i in names(MC.seu.list)){
  filename = paste('./plot/supercell/','supercell_dimplot_',i,'.pdf',sep = '')
  DimPlot(MC.seu.list[[i]], cols = cell.colors1, label=T, repel = T)+
    umap.theme+
    NoLegend()
  ggsave(filename, width = 6, height = 6)

}

supercell_plot(
  MC.list[[1]]$graph.supercells, 
  group = MC.list[[1]]$cell_line, 
  color.use = cell.colors1[levels(as.factor(MC.list[[1]]$cell_line))],
  seed = 1, 
  alpha = -pi/2,
  main  = "Metacells colored by cell line assignment"
)


library(patchwork)
library(stringr)

plot1 <- DimPlot(MC.seu.list[[3]], reduction = "umap")

plot2 <- FeaturePlot(MC.seu.list[[3]], features = 'purity')

plot1 + plot2

plot3 <- DimPlot(subset(MC.seu.list[[1]],subset = purity == 1))

plot4 <- FeaturePlot(subset(MC.seu.list[[1]],subset = purity == 1), features = 'purity')

plot3 + plot4


