lapply(c("dplyr","Seurat"), library, character.only = T)
library(ggplot2)
library(cowplot)
library(patchwork)
library(stringr)
library(DoubletFinder)
library(future)
library(cols4all)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(COSG)
setwd('~/data/PBMC/')

plan('multisession', workers = 5)
plan()

sp.list <- readRDS('./rds/sp.list.quality_filter.rds')

species <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
             'Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')

pbmc.list <- list()
for (i in names(sp.list)){
  pbmc.list[[i]] <- sp.list[[i]][ ,sp.list[[i]]$multi_annotation %in% c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','DCs')]
}

for (i in names(sp.list)){
  pbmc.list[[i]][['pbmc_annotation']] <- pbmc.list[[i]][['multi_annotation']]
  pbmc.list[[i]]$pbmc_annotation <- pbmc.list[[i]]$pbmc_annotation %>%
    as.character() %>%
    str_replace_all('mDCs|pDCs','DCs') %>%
    as.factor()
}

for (i in names(pbmc.list)){
  print(i)
  print(levels(pbmc.list[[i]]$pbmc_annotation))
}

mouse.marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Mouse.xlsx',sheet = 'Sheet 1')
de_marker <- unique(mouse.marker$gene)

mouse_num <- pbmc.list[['Mouse']][de_marker,]
mouse_num <- CreateSeuratObject(GetAssayData(mouse_num, assay = 'RNA', slot = 'counts'),
                                meta.data = mouse_num@meta.data[,c('multi_annotation','pbmc_annotation')])
mouse_num <- mouse_num@meta.data[,c('nFeature_RNA','pbmc_annotation')]

mouse_sum <- mouse_num %>%
  group_by(pbmc_annotation) %>%
  summarize(
    Mouse = mean(nFeature_RNA, na.rm = TRUE)
  )

human.marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Human.xlsx',sheet = 'Sheet 1')
de_marker <- unique(human.marker$gene)

human_num <- pbmc.list[['Human']][de_marker,]
human_num <- CreateSeuratObject(GetAssayData(human_num, assay = 'RNA', slot = 'counts'),
                                meta.data = human_num@meta.data[,c('multi_annotation','pbmc_annotation')])
human_num <- human_num@meta.data[,c('nFeature_RNA','pbmc_annotation')]

human_sum <- human_num %>%
  group_by(pbmc_annotation) %>%
  summarize(
    Human = mean(nFeature_RNA, na.rm = TRUE)
  )

orders <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
            'Cattle','Rat','Monkey','Chimpanzee','Human')

matrix.list <- list()
for (i in names(pbmc.list)) {
  matrix.list[[i]] <- GetAssayData(pbmc.list[[i]], assay = 'RNA', slot = 'counts')
}

mouse.marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Mouse.xlsx',sheet = 'Sheet 1')
de_marker <- unique(mouse.marker$gene)
mouse.list <- list()
# gene to mouse orthologue
for (i in names(matrix.list)[names(matrix.list) != 'Mouse']) {

  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_mouse.xlsx", sheet = i)
  gene_name <- rownames(matrix.list[[i]]) %>% as.data.frame()
  colnames(gene_name) <- paste0(i, '.gene.name')
  column_name <- colnames(gene_name)[1]

  gene_name[["Mouse.gene.name"]] <- all_ortholog$Mouse.gene.name[match(gene_name[[column_name]], all_ortholog[[column_name]])]
  gene_name[gene_name == ''] <- NA
  gene_name <- na.omit(gene_name) 
  mouse.list[[i]] <- matrix.list[[i]][gene_name[[column_name]],]
  rownames(mouse.list[[i]]) <- gene_name$Mouse.gene.name
  mouse.list[[i]] <- mouse.list[[i]][rownames(mouse.list[[i]]) %in% de_marker,]
}

ortho.list <- list()
for(i in names(mouse.list)){
  ortho.list[[i]] <- CreateSeuratObject(mouse.list[[i]],
                                        meta.data = pbmc.list[[i]]@meta.data[,c('multi_annotation','pbmc_annotation')])
}

ortho_num <- list()
ortho_sum <- list()
for (i in names(ortho.list)){
  print(i)
  ortho_num[[i]] <- ortho.list[[i]]@meta.data[,c('nFeature_RNA','pbmc_annotation')]
  ortho_sum[[i]] <- ortho_num[[i]] %>%
    group_by(pbmc_annotation) %>%
    summarize(
      Avg_nFeature = mean(nFeature_RNA, na.rm = TRUE)
    )
  colnames(ortho_sum[[i]])[2] <- i
}

for (i in names(ortho_sum)){
  mouse_sum <- merge(mouse_sum, ortho_sum[[i]], by = 'pbmc_annotation', all = T)
}


human.marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Human.xlsx',sheet = 'Sheet 1')
de_marker <- unique(human.marker$gene)
human.list <- list()
# gene to human orthologue
for (i in names(matrix.list)[names(matrix.list) != 'Human']) {

  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_human.xlsx", sheet = i)
  gene_name <- rownames(matrix.list[[i]]) %>% as.data.frame()
  colnames(gene_name) <- paste0(i, '.gene.name')
  column_name <- colnames(gene_name)[1]

  gene_name[["Human.gene.name"]] <- all_ortholog$Human.gene.name[match(gene_name[[column_name]], all_ortholog[[column_name]])]
  gene_name[gene_name == ''] <- NA
  gene_name <- na.omit(gene_name)
  human.list[[i]] <- matrix.list[[i]][gene_name[[column_name]],]
  rownames(human.list[[i]]) <- gene_name$Human.gene.name
  human.list[[i]] <- human.list[[i]][rownames(human.list[[i]]) %in% de_marker,]
}


ortho.list <- list()
for(i in names(human.list)){
  ortho.list[[i]] <- CreateSeuratObject(human.list[[i]],
                                        meta.data = pbmc.list[[i]]@meta.data[,c('multi_annotation','pbmc_annotation')])
}

ortho_num <- list()
ortho_sum <- list()
for (i in names(ortho.list)){
  print(i)
  ortho_num[[i]] <- ortho.list[[i]]@meta.data[,c('nFeature_RNA','pbmc_annotation')]
  ortho_sum[[i]] <- ortho_num[[i]] %>%
    group_by(pbmc_annotation) %>%
    summarize(
      Avg_nFeature = mean(nFeature_RNA, na.rm = TRUE)
    )
  colnames(ortho_sum[[i]])[2] <- i
}

for (i in names(ortho_sum)){
  human_sum <- merge(human_sum, ortho_sum[[i]], by = 'pbmc_annotation', all = T)
}

head(mouse_sum)
head(human_sum)
mouse_propo <- mouse_sum %>%
  mutate(
    Catfish = Catfish / Mouse,
    Jacopever = Jacopever / Mouse,
    Turtle = Turtle / Mouse,
    Chicken = Chicken / Mouse,
    Bat = Bat / Mouse,
    Pig = Pig / Mouse,
    Cattle = Cattle / Mouse,
    Rat = Rat / Mouse,
    Monkey = Monkey / Mouse,
    Chimpanzee = Chimpanzee / Mouse,
    Human = Human / Mouse
  )

human_propo <- human_sum %>%
  mutate(
    Catfish = Catfish / Human,
    Jacopever = Jacopever / Human,
    Turtle = Turtle / Human,
    Chicken = Chicken / Human,
    Bat = Bat / Human,
    Pig = Pig / Human,
    Cattle = Cattle / Human,
    Mouse = Mouse / Human,
    Rat = Rat / Human,
    Monkey = Monkey / Human,
    Chimpanzee = Chimpanzee / Human
  )

head(mouse_propo)

mouse_propo <- mouse_propo %>%
  dplyr::select(-Mouse)


mouse_propo <- tidyr::pivot_longer(
  mouse_propo,
  cols = -pbmc_annotation,
  names_to = "Species",
  values_to = "Proportion"
)

mouse_propo$Species <- factor(mouse_propo$Species, levels = rev(c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
                                                        'Cattle','Rat','Monkey','Chimpanzee','Human')))

ggplot(mouse_propo, aes(x = Species, y = Proportion, group = pbmc_annotation, color = pbmc_annotation)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cell.colors1) +
  labs(x = "", y = "Proportion of orthologs with mouse",color = 'Cell category') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave('./plot/signature/Proportion of orthologs with mouse.pdf')

head(human_propo)

human_propo <- human_propo %>%
  dplyr::select(-Human)

human_propo <- tidyr::pivot_longer(
  human_propo,
  cols = -pbmc_annotation,
  names_to = "Species",
  values_to = "Proportion"
)

human_propo$Species <- factor(human_propo$Species, levels = rev(c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig',
                                                                  'Cattle','Mouse','Rat','Monkey','Chimpanzee')))

ggplot(human_propo, aes(x = Species, y = Proportion, group = pbmc_annotation, color = pbmc_annotation)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = cell.colors1) +
  labs(x = "", y = "Proportion of orthologs with human",color = 'Cell category') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

ggsave('./plot/signature/Proportion of orthologs with human.pdf')
######## use cosg to test markers #########

# Set cell identities for the 'Human' dataset
Idents(sp.list[['Human']]) <- 'multi_annotation'

# Perform cosg analysis on the 'Human' dataset
human_cosg <- cosg(sp.list[['Human']], assay = 'SCT', n_genes_user = 400)

# Extract gene scores and annotations from cosg results
human_cosg.gene <- as.data.frame(human_cosg[[1]])
human_cosg.score <- as.data.frame(human_cosg[[2]])

# Create a new dataframe to store cosg results
human_cosg_new <- cbind(human_cosg.gene[,1], human_cosg.score[,1]) %>% as.data.frame()
colnames(human_cosg_new) <- c('gene', 'score')
human_cosg_new$multi_annotation <- colnames(human_cosg[[1]])[1]

# Iterate over remaining columns of cosg results and append to the new dataframe
for (i in 2:ncol(human_cosg[[1]])) {
  human_cosg_newi <- cbind(human_cosg.gene[,i], human_cosg.score[,i]) %>% as.data.frame()
  colnames(human_cosg_newi) <- c('gene', 'score')
  human_cosg_newi$multi_annotation <- colnames(human_cosg.gene)[i]
  human_cosg_new <- rbind(human_cosg_new, human_cosg_newi)
}

############ read markers of human ############

# Read markers of human from an Excel file
human_marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Human.xlsx', sheet = 'Sheet 1')

# Perform a left join to merge cosg results with human markers
human_cosg_new <- left_join(human_cosg_new, human_marker, relationship = "many-to-many")

# Remove rows with NA values
human_cosg_new <- na.omit(human_cosg_new)

# Check if cell annotation matches cluster and store the result in 'overlay' column
for (i in 1:nrow(human_cosg_new)) {
  human_cosg_new$overlay[i] = human_cosg_new$multi_annotation[i] == human_cosg_new$cluster[i]
}

# Filter rows where the cell annotation matches the cluster
human_cosg_new <- filter(human_cosg_new, overlay == 'TRUE')
human_cosg_new <- human_cosg_new[,1:8]

############### scale gene expression ##############

# Fetch gene expression data for the selected genes from the 'Human' dataset
scale.data <- Seurat::FetchData(object = sp.list[['Human']], vars = unique(human_cosg_new$gene), slot = 'data')
scale.data$multi_annotation <- sp.list[['Human']]@meta.data[['multi_annotation']]

# Scale the gene expression data
scale.data <- lapply(X = unique(x = scale.data$multi_annotation), FUN = function(ident) {
  data.use <- scale.data[scale.data$multi_annotation == ident, 1:(ncol(x = scale.data) - 1), drop = FALSE]
  avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
    return(mean(x = expm1(x = x)))
  })
  pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
  res$gene <- rownames(res)
  return(res)
}) %>% do.call("rbind", .) %>% data.frame()

# Scale the average expression values
scale.data <- purrr::map_df(unique(scale.data$gene), function(x) {
  tmp <- scale.data %>% dplyr::filter(gene == x)
  avg.exp.scale <- tmp %>% dplyr::select(avg.exp) %>% scale(.) %>% scales::rescale(., to = c(0, 1))
  tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
  return(tmp)
})

# Update the 'human_cosg_new' dataframe with scaled gene expression values
for (i in 1:nrow(scale.data)) {
  gene_name <- scale.data$gene[i]
  id <- as.character(scale.data$id[i])
  matching_row <- which(human_cosg_new$gene == gene_name)
  if (length(matching_row) > 0) {
    human_cosg_new[matching_row, id] <- scale.data$avg.exp.scaled[i]
  }
}

# Identify cells to exclude based on thresholds
human_cosg_new$score <- as.numeric(human_cosg_new$score)

cell_use <- colnames(human_cosg_new)[9:15]
for (i in 1:nrow(human_cosg_new)) {
  exclude_cell <- human_cosg_new[i, 'multi_annotation']
  exclude_cell <- setdiff(cell_use, exclude_cell)
  exclude_cell <- na.omit(exclude_cell)
  mean_score <- mean(filter(human_cosg_new,multi_annotation == human_cosg_new[i, 'multi_annotation'])$score)
  human_cosg_new[i, 'remove'] <- any(human_cosg_new[i, exclude_cell] > 0.5) | 
    human_cosg_new[i, 'pct.2'] > 0.3 | 
    human_cosg_new[i, 'score'] < mean_score
}

# Filter out rows to be removed
human_cosg_new <- filter(human_cosg_new, remove == 'FALSE')

# Display the frequency of cell annotations
table(human_cosg_new$multi_annotation)

# Save the results to an Excel file
filename <- paste('./file/cosg.markers_', 'human', '_multi_annotation.xlsx', sep = '')
openxlsx::write.xlsx(human_cosg_new, file = filename)

######## mouse cosg ##############

# Set cell identities for the 'Mouse' dataset
Idents(sp.list[['Mouse']]) <- 'multi_annotation'

# Perform cosg analysis on the 'Mouse' dataset
mouse_cosg <- cosg(sp.list[['Mouse']], assay = 'SCT', n_genes_user = 400)

# Extract gene scores and annotations from cosg results
mouse_cosg.gene <- as.data.frame(mouse_cosg[[1]])
mouse_cosg.score <- as.data.frame(mouse_cosg[[2]])

# Create a new dataframe to store cosg results
mouse_cosg_new <- cbind(mouse_cosg.gene[,1], mouse_cosg.score[,1]) %>% as.data.frame()
colnames(mouse_cosg_new) <- c('gene', 'score')
mouse_cosg_new$multi_annotation <- colnames(mouse_cosg[[1]])[1]

# Iterate over remaining columns of cosg results and append to the new dataframe
for (i in 2:ncol(mouse_cosg[[1]])) {
  mouse_cosg_newi <- cbind(mouse_cosg.gene[,i], mouse_cosg.score[,i]) %>% as.data.frame()
  colnames(mouse_cosg_newi) <- c('gene', 'score')
  mouse_cosg_newi$multi_annotation <- colnames(mouse_cosg.gene)[i]
  mouse_cosg_new <- rbind(mouse_cosg_new, mouse_cosg_newi)
}


############ read markers of mouse ############

# Read markers of mouse from an Excel file
mouse_marker <- openxlsx::read.xlsx('./file/markers_multi_annotation_Mouse.xlsx', sheet = 'Sheet 1')

# Perform a left join to merge cosg results with mouse markers
mouse_cosg_new <- left_join(mouse_cosg_new, mouse_marker, relationship = "many-to-many")

# Remove rows with NA values
mouse_cosg_new <- na.omit(mouse_cosg_new)

# Check if cell annotation matches cluster and store the result in 'overlay' column
for (i in 1:nrow(mouse_cosg_new)) {
  mouse_cosg_new$overlay[i] = mouse_cosg_new$multi_annotation[i] == mouse_cosg_new$cluster[i]
}

# Filter rows where the cell annotation matches the cluster
mouse_cosg_new <- filter(mouse_cosg_new, overlay == 'TRUE')
mouse_cosg_new <- mouse_cosg_new[,1:8]


############### scale gene expression ##############

# Get common genes between cosg results and the 'Mouse' dataset
common_genes <- intersect(mouse_cosg_new$gene, rownames(sp.list[['Mouse']]@assays$SCT))

# Fetch gene expression data for the common genes from the 'Mouse' dataset
scale.data <- Seurat::FetchData(object = sp.list[['Mouse']], vars = common_genes, slot = 'data')
scale.data$multi_annotation <- sp.list[['Mouse']]@meta.data[['multi_annotation']]

# Scale the gene expression data
scale.data <- lapply(X = unique(x = scale.data$multi_annotation), FUN = function(ident) {
  data.use <- scale.data[scale.data$multi_annotation == ident, 1:(ncol(x = scale.data) - 1), drop = FALSE]
  avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
    return(mean(x = expm1(x = x)))
  })
  pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
  res$gene <- rownames(res)
  return(res)
}) %>% do.call("rbind", .) %>% data.frame()

# Scale the average expression values
scale.data <- purrr::map_df(unique(scale.data$gene), function(x) {
  tmp <- scale.data %>% dplyr::filter(gene == x)
  avg.exp.scale <- tmp %>% dplyr::select(avg.exp) %>% scale(.) %>% scales::rescale(., to = c(0, 1))
  tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
  return(tmp)
})

# Update the 'mouse_cosg_new' dataframe with scaled gene expression values
for (i in 1:nrow(scale.data)) {
  gene_name <- scale.data$gene[i]
  id <- as.character(scale.data$id[i])
  matching_row <- which(mouse_cosg_new$gene == gene_name)
  if (length(matching_row) > 0) {
    mouse_cosg_new[matching_row, id] <- scale.data$avg.exp.scaled[i]
  }
}

# Identify cells to exclude based on thresholds
mouse_cosg_new$score <- as.numeric(mouse_cosg_new$score)

cell_use <- colnames(mouse_cosg_new)[9:14]
for (i in 1:nrow(mouse_cosg_new)) {
  exclude_cell <- mouse_cosg_new[i, 'multi_annotation']
  exclude_cell <- setdiff(cell_use, exclude_cell)
  exclude_cell <- na.omit(exclude_cell)
  mean_score <- mean(filter(mouse_cosg_new,multi_annotation == mouse_cosg_new[i, 'multi_annotation'])$score)
  mouse_cosg_new[i, 'remove'] <- any(mouse_cosg_new[i, exclude_cell] > 0.5) | 
    mouse_cosg_new[i, 'pct.2'] > 0.3 | 
    mouse_cosg_new[i, 'score'] < mean_score
}
# Filter out rows to be removed
mouse_cosg_new <- filter(mouse_cosg_new, remove == 'FALSE')

# Display the frequency of cell annotations
table(mouse_cosg_new$multi_annotation)

# Save the results to an Excel file
filename <- paste('./file/cosg.markers_', 'mouse', '_multi_annotation.xlsx', sep = '')
openxlsx::write.xlsx(mouse_cosg_new, file = filename)

############# add ortholog symbol ###############

# Read the ortholog data between human and mouse from an Excel file
human_mouse <- openxlsx::read.xlsx('./protein/ortholog_to_mouse.xlsx', sheet = 'Human')

# Add mouse symbol to the human cosg dataframe
human_cosg_new$Mouse.gene.name <- human_mouse$Mouse.gene.name[match(human_cosg_new$gene, human_mouse$Human.gene.name)]

# Add human symbol to the mouse cosg dataframe
mouse_cosg_new$Human.gene.name <- human_mouse$Human.gene.name[match(mouse_cosg_new$gene, human_mouse$Mouse.gene.name)]

# Remove rows with NA values
human_cosg_new <- na.omit(human_cosg_new)
mouse_cosg_new <- na.omit(mouse_cosg_new)

# # Get highly variable genes from human and mouse datasets
# human_HVG <- VariableFeatures(sp.list[['Human']])
# mouse_HVG <- VariableFeatures(sp.list[['Mouse']])

# Mark rows for removal based on the presence of HVGs in the opposite species
# human_cosg_new$remove <- human_cosg_new$Mouse.gene.name %in% mouse_HVG
# mouse_cosg_new$remove <- mouse_cosg_new$Human.gene.name %in% human_HVG
# 
# # Filter out rows marked for removal
# human_cosg_new <- human_cosg_new %>% filter(remove == 'TRUE')
# 
# mouse_cosg_new <- mouse_cosg_new %>% filter(remove == 'TRUE')

table(human_cosg_new$multi_annotation)
table(mouse_cosg_new$multi_annotation)

# Filter cell types of interest
celltype_use <- c('B cells', 'T cells', 'NK cells', 'Monocytes', 'pDCs', 'mDCs','Neutrophils','Platelets')
human_cosg_new <- human_cosg_new[human_cosg_new$multi_annotation %in% celltype_use, ]
mouse_cosg_new <- mouse_cosg_new[mouse_cosg_new$multi_annotation %in% celltype_use, ]


############## cell type scoring ###############

# Select top markers for each cell type in human and mouse datasets
top_human <- human_cosg_new[, c(1, 2, 3, 17)]
colnames(top_human)[1] <- 'Human.gene.name'
top_mouse <- mouse_cosg_new[, c(1, 2, 3, 16)]
colnames(top_mouse)[1] <- 'Mouse.gene.name'

# Combine top markers from human and mouse datasets
top_marker <- rbind(top_human, top_mouse)

# Group markers by cell annotation and keep only non-duplicated markers
top_marker <- top_marker %>%
  group_by(multi_annotation) %>%
  mutate(dup = duplicated(Human.gene.name) | duplicated(Mouse.gene.name))

# Filter out duplicated markers
top_marker <- top_marker %>% filter(dup == 'FALSE')

# Perform cell type scoring by assigning ranks based on marker scores
top_marker <- top_marker %>%
  mutate(score = as.numeric(score)) %>%
  group_by(multi_annotation)

# Display the frequency of cell annotations
table(top_marker$multi_annotation)

############ all ortholog genes ########
for (i in species) {

  all_ortholog <- openxlsx::read.xlsx("./protein/ortholog_to_human.xlsx", sheet = i)
  
  top_marker[[paste0(i, ".gene.name")]] <- all_ortholog[[paste0(i, ".gene.name")]][match(top_marker$Human.gene.name, all_ortholog$Human.gene.name)]
}

top_marker <- top_marker[ ,c(2:3,6:12,4,13:15,1)]
top_marker <- as.data.frame(top_marker)

########## save markers that expressed in cells ##########
animal_name <- c('Catfish','Jacopever','Turtle','Chicken','Bat','Pig','Cattle','Rat','Monkey','Chimpanzee')
animal_name <- species
for (i in animal_name){
  symbol_name <- paste0(i, ".gene.name")
  common_genes <- intersect(top_marker[,symbol_name],rownames(sp.list[[i]]))
  top_marker[,symbol_name] <- ifelse(top_marker[,symbol_name] %in% common_genes, top_marker[,symbol_name], NA)
}

openxlsx::write.xlsx(top_marker, './file/conserved_marker_cross_species_orig.xlsx')

##### Fetch gene expression data for the common genes from the  dataset ####
top_marker <- openxlsx::read.xlsx('./file/conserved_marker_cross_species_orig.xlsx')

used_genes <- na.omit(top_marker[,'Catfish.gene.name']) 

scale.data <- GetAssayData(object = sp.list[['Catfish']], slot = 'count') %>% as.matrix()

total_counts <- colSums(scale.data)

relative_expression <- t(t(scale.data) / total_counts)

scaling_factor <- 1e6

relative_expression <- relative_expression * scaling_factor

cpm <- log10(relative_expression+1)

cpm <- cpm[used_genes, ] %>% t() %>% as.data.frame()

cpm$multi_annotation <- sp.list[['Catfish']]@meta.data[['multi_annotation']]

# Scale the gene expression data
cpm <- lapply(X = unique(x = cpm$multi_annotation), FUN = function(ident) {
  data.use <- cpm[cpm$multi_annotation == ident, 1:(ncol(x = cpm) - 1), drop = FALSE]
  avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
    return(mean(x = expm1(x = x)))
  })
  pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
  res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
  res$gene <- rownames(res)
  return(res)
}) %>% do.call("rbind", .) %>% data.frame()

tail(cpm)
cpm$gene <- gsub("\\.\\d+$", "", cpm$gene)

cpm$multi_annotation <- top_marker$multi_annotation[match(cpm$gene, top_marker$Catfish.gene.name)]

cpm$Species <- 'tfd' ### Abbreviations

################# all species #########
animal_name <- c('Jacopever','Turtle','Chicken','Bat','Pig','Cattle','Mouse','Rat','Monkey','Chimpanzee','Human')
abb_name <- c('sslg','pss','gga','ray','ssc','bta','mmu','rno','mcc','ptr','hsa')

for (i in 1:length(animal_name)){
  used_animal <- paste0(animal_name[i], '.gene.name')
  used_genes <- na.omit(top_marker[ ,used_animal])
  scale.data <- GetAssayData(object = sp.list[[animal_name[i]]], slot = 'count') %>% as.matrix()

  total_counts <- colSums(scale.data)
  
  relative_expression <- t(t(scale.data) / total_counts)
  
  scaling_factor <- 1e6
  
  relative_expression <- relative_expression * scaling_factor
  
  cpmi <- log10(relative_expression+1)
  
  cpmi <- cpmi[used_genes, ] %>% t() %>% as.data.frame()
  
  cpmi$multi_annotation <- sp.list[[animal_name[i]]]@meta.data[['multi_annotation']]
  
  # Scale the gene expression data
  cpmi <- lapply(X = unique(x = cpmi$multi_annotation), FUN = function(ident) {
    data.use <- cpmi[cpmi$multi_annotation == ident, 1:(ncol(x = cpmi) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
    res <- data.frame(id = ident, avg.exp = avg.exp, pct.exp = pct.exp * 100)
    res$gene <- rownames(res)
    return(res)
  }) %>% do.call("rbind", .) %>% data.frame()
  
  tail(cpmi)
  cpmi$gene <- gsub("\\.\\d+$", "", cpmi$gene)
  
  cpmi$multi_annotation <- top_marker$multi_annotation[match(cpmi$gene, top_marker[ ,used_animal])]
  
  cpmi$Species <- abb_name[i]
  
  cpm <- rbind(cpm, cpmi)
}

saveRDS(cpm,'./rds/Mean_expression.rds')
###### plotting #########

cell.colors1 <- c('HSCs' = '#be9457',
                  "Neutrophils" = "#cedfef",
                  "Cycling cells" = "#fabb6e",
                  "Platelets" = "#1663a9",
                  "Erythrocytes" = "#b9181a",
                  "B cells" = "#fc8002", 
                  "Monocytes" = "#369f2d",
                  "mDCs" = "#846e89", 
                  "NK cells" = "#fac7b3",
                  "pDCs" = "#4995c6", 
                  "T cells" = "#ee4431",
                  'DCs' = '#c6d182'
)
suffix_order <- c('tfd','sslg','pss','gga','ray','ssc','bta','mmu','rno','mcc','ptr','hsa')

animal_name <- species

celltype_use <- c('B cells', 'T cells', 'NK cells', 'Monocytes', 'pDCs', 'mDCs','Neutrophils','Platelets')

for(i in 1:length(celltype_use)){
  used_marker <- filter(top_marker,multi_annotation == celltype_use[i])
  used_genes <- used_marker$Catfish.gene.name
  for(j in 2:length(animal_name)){
    used_genes <- c(used_genes, used_marker[[paste0(animal_name[j],'.gene.name')]])
  }
  used_genes <- used_genes %>% na.omit() %>% unique()
  
  cpmp <- cpm[cpm$gene %in% used_genes,]
  
  cpmp$group <- paste(cpmp$id, cpmp$Species, sep = '.')
  
  # Scale the average expression values
  cpmp <- purrr::map_df(unique(cpmp$gene), function(x) {
    tmp <- cpmp %>% dplyr::filter(gene == x)
    avg.exp.scale <- tmp %>% dplyr::select(avg.exp) %>% scale(.) %>% scales::rescale(., to = c(0, 10))
    tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
    return(tmp)
  })
  

  grouped_levels <- split(levels(as.factor(cpmp$group)), sapply(strsplit(levels(as.factor(cpmp$group)), "\\."), tail, 1))
  
  sorted_levels <- unlist(lapply(suffix_order, function(suffix) {
    if (suffix %in% names(grouped_levels)) {
      sorted_group <- sort(grouped_levels[[suffix]])
    }
  }))

  cpmp$group <- factor(cpmp$group, levels = sorted_levels)
  
  levels(cpmp$group)
  
  fill_color <- ifelse(grepl(celltype_use[i], levels(cpmp$group)), cell.colors1[celltype_use[i]], "white")
  
  ggplot(cpmp, aes(x = group, y = avg.exp.scaled, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Mean expression, log10(CPM + 1)", title = celltype_use[i]) +
    theme_classic() +
    scale_fill_manual(values = fill_color) +
    NoLegend() +
    scale_y_continuous(limits = c(0, 10)) +
    geom_vline(xintercept = c(7,14,18,25,31,36,43,49,57,63,68)+0.5, 
               linetype = "dashed", color = "red", na.rm = TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank())
  filename = paste0('./plot/signature/boxplot_', gsub(" ", "_", celltype_use[i]), '_signature.pdf')
  
  ggsave(filename, height = 4)
}


plot.list <- list()
for(i in celltype_use){
  used_marker <- filter(top_marker,multi_annotation == i)
  used_genes <- used_marker$Catfish.gene.name
  for(j in 2:length(animal_name)){
    used_genes <- c(used_genes, used_marker[[paste0(animal_name[j],'.gene.name')]])
  }
  used_genes <- used_genes %>% na.omit() %>% unique()
  
  cpmp <- cpm[cpm$gene %in% used_genes,]
  
  cpmp$group <- paste(cpmp$id, cpmp$Species, sep = '.')
  
  # Scale the average expression values
  cpmp <- purrr::map_df(unique(cpmp$gene), function(x) {
    tmp <- cpmp %>% dplyr::filter(gene == x)
    avg.exp.scale <- tmp %>% dplyr::select(avg.exp) %>% scale(.) %>% scales::rescale(., to = c(0, 10))
    tmp$avg.exp.scaled <- as.numeric(avg.exp.scale)
    return(tmp)
  })

  grouped_levels <- split(levels(as.factor(cpmp$group)), sapply(strsplit(levels(as.factor(cpmp$group)), "\\."), tail, 1))
  

  sorted_levels <- unlist(lapply(suffix_order, function(suffix) {
    if (suffix %in% names(grouped_levels)) {
      sorted_group <- sort(grouped_levels[[suffix]])
    }
  }))

  cpmp$group <- factor(cpmp$group, levels = sorted_levels)
  
  levels(cpmp$group)
  
  fill_color <- ifelse(grepl(i, levels(cpmp$group)), cell.colors1[i], "white")
  common_plot <- 
    ggplot(cpmp, aes(x = group, y = avg.exp.scaled, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = "Mean expression, log10(CPM + 1)", title = i) +
    theme_classic() +
    scale_fill_manual(values = fill_color) +
    NoLegend() +
    scale_y_continuous(limits = c(0, 10), 
                       breaks = seq(0, 10, length.out = 3),
                       labels = function(x) sprintf("%.2f", x)) +
    geom_vline(xintercept = c(7,14,18,25,31,36,43,49,57,63,68)+0.5,
               linetype = "dashed", color = "red", na.rm = TRUE) +
    theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14),
          axis.title.x = element_blank())
  
  if (i %in% c(celltype_use[8])) {
    plot.list[[i]] <- common_plot +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12))
  } else {
    plot.list[[i]] <- common_plot +
      theme(axis.text.x = element_blank())
  }
}


for (i in 1:8) {
  if (i != 4) { 
    plot.list[[i]] <- plot.list[[i]] + 
      theme(axis.title.y = element_blank())
  }
}

plots <- plot.list[[1]]/plot.list[[2]]/plot.list[[3]]/plot.list[[4]]/plot.list[[5]]/plot.list[[6]]/plot.list[[7]]/plot.list[[8]]
plots <- plots + plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 1, 1, 1, 1))

plots

ggsave('./plot/signature/boxplot_total_signature.pdf', width = 11, height = 12)

###################### cell scoring #############
library(UCell)
library(GSEABase)
library(BiocParallel)
set.seed(1)

############### human  mouse ###########
human <- sp.list[['Human']]

# Select top markers for each cell type based on scores
top_marker <- top_marker %>% group_by(multi_annotation) %>% slice_max(score, n = 20)

signatures <- list()
for(i in 1:length(celltype_use)){
  signatures[[i]] = filter(top_marker,multi_annotation == celltype_use[i])$Human.gene.name
}

names(signatures) <- paste(gsub(' ','_',celltype_use),'signature',sep = '_')
  

human <- AddModuleScore_UCell(human,features = signatures,ncores = 3)

human_metadata <- human@meta.data[,c(22:29)]
human_metadata <- cbind(human_metadata,Embeddings(human,reduction = 'umap'))
human_metadata$id <- rownames(human_metadata)
human_metadata = human_metadata[,c(11,9,10,1:8)]

DimPlot(human,label = T,repel = T)

##  mouse
mouse <- sp.list[['Mouse']]

signatures <- list()
for(i in 1:length(celltype_use)){
  signatures[[i]] = filter(top_marker,multi_annotation == celltype_use[i])$Mouse.gene.name %>% na.omit()
}

names(signatures) <- paste(gsub(' ','_',celltype_use),'signature',sep = '_')


mouse <- AddModuleScore_UCell(mouse,features = signatures,ncores = 3)

mouse_metadata <- mouse@meta.data[,c(22:29)]
mouse_metadata <- cbind(mouse_metadata,Embeddings(mouse,reduction = 'umap'))
mouse_metadata$id <- rownames(mouse_metadata)
mouse_metadata = mouse_metadata[,c(11,9,10,1:8)]

DimPlot(mouse,label = T,repel = T)

######################### other species ###################
celltype_use <- c('B cells', 'T cells', 'NK cells', 'Monocytes', 'pDCs', 'mDCs','Neutrophils','Platelets')
animal_name <- species

animal_metadata <- list()

for (i in animal_name){
  animal <- sp.list[[i]]
  used_animal <- paste0(i, '.gene.name')
  
  signatures <- list()
  for(j in celltype_use){
    signatures[[j]] = filter(top_marker,multi_annotation == j)[[used_animal]] %>%
      na.omit()
  }
  
  names(signatures) <- paste(gsub(' ','_',celltype_use),'signature',sep = '_')
  
  animal <- AddModuleScore_UCell(animal,features = signatures,ncores = 3)
  
  animal_meta <- animal@meta.data[,c((ncol(animal@meta.data)-8):ncol(animal@meta.data))]
  animal_meta <- cbind(animal_meta,Embeddings(animal,reduction = 'umap'))
  animal_meta$id <- rownames(animal_meta)
  animal_meta = animal_meta[,c(ncol(animal_meta),1:(ncol(animal_meta)-1))]
  
  animal_metadata[[i]] <- animal_meta
}

for (i in names(animal_metadata)){
  filename = paste0('./file/','UCell_Score_',names(animal_metadata[i]),'.xlsx')
  openxlsx::write.xlsx(animal_metadata[[i]],filename)
}

suffix_order <- c('tfd','sslg','pss','gga','ray','ssc','bta','mmu','rno','mcc','ptr','hsa')

animal_name <- species

animal_meta <- animal_metadata[[1]]
animal_meta$Species <- suffix_order[1]
for(i in 2:length(animal_metadata)){
  animal_metai <- animal_metadata[[i]]
  animal_metai$Species <- suffix_order[i]
  animal_meta <- rbind(animal_meta, animal_metai)
}

animal_meta$group <- paste(animal_meta$multi_annotation, animal_meta$Species, sep = '.')

grouped_levels <- split(levels(as.factor(animal_meta$group)), sapply(strsplit(levels(as.factor(animal_meta$group)), "\\."), tail, 1))

sorted_levels <- unlist(lapply(suffix_order, function(suffix) {
  if (suffix %in% names(grouped_levels)) {
    sorted_group <- sort(grouped_levels[[suffix]])
  }
}))

animal_meta$group <- factor(animal_meta$group, levels = sorted_levels)

levels(animal_meta$group)

celltype_use <- c('B cells', 'T cells', 'NK cells', 'Monocytes', 'pDCs', 'mDCs','Neutrophils','Platelets')

for (i in celltype_use){
  fill_color <- ifelse(grepl(i, levels(animal_meta$group)), cell.colors1[i], "white")
  
  used_cell <- paste0(str_replace_all(i,' ','_'),'_signature_UCell')
  UCell_score <- animal_meta[,c(used_cell,'group','Species')]
  
  colnames(UCell_score)
  colnames(UCell_score)[1] <- 'Score'
  
  ggplot(UCell_score, aes(x = group, y = Score, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "", y = 'UCell score of conserved marker', title = i) +
    theme_classic() +
    scale_fill_manual(values = fill_color) +
    NoLegend() +
    geom_vline(xintercept = c(7,14,18,25,31,36,43,49,57,63,68) + 0.5, 
               linetype = "dashed", color = "red", na.rm = TRUE) +
    scale_y_continuous(limits = c(0, max(UCell_score$Score))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank())
  filename = paste0('./plot/signature/boxplot_', gsub(" ", "_", i), '_UCell.pdf')
  
  ggsave(filename, height = 4)
}

plot.list <- list()
for (i in celltype_use){
  fill_color <- ifelse(grepl(i, levels(animal_meta$group)), cell.colors1[i], "white")
  
  used_cell <- paste0(str_replace_all(i,' ','_'),'_signature_UCell')
  UCell_score <- animal_meta[,c(used_cell,'group','Species')]
  
  colnames(UCell_score)
  colnames(UCell_score)[1] <- 'Score'
  
  y_max <- max(UCell_score$Score)
  y_min <- min(UCell_score$Score)
  
  if (i %in% c(celltype_use[8])) {
    plot.list[[i]] <- 
      ggplot(UCell_score, aes(x = group, y = Score, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      labs(x = "", y = 'UCell score of conserved marker', title = i) +
      theme_classic() +
      scale_fill_manual(values = fill_color) +
      NoLegend() +
      geom_vline(xintercept = c(7,14,18,25,31,36,43,49,57,63,68)+0.5, 
                 linetype = "dashed", color = "red", na.rm = TRUE) +
      scale_y_continuous(limits = c(y_min, y_max), 
                         breaks = seq(y_min, y_max, length.out = 3),
                         labels = function(x) sprintf("%.2f", x)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            axis.title.x = element_blank())
  } else {
    plot.list[[i]] <- 
      ggplot(UCell_score, aes(x = group, y = Score, fill = group)) +
      geom_boxplot(outlier.shape = NA) +
      labs(x = "", y = 'UCell score of conserved marker', title = i) +
      theme_classic() +
      scale_fill_manual(values = fill_color) +
      NoLegend() +
      geom_vline(xintercept = c(7,14,18,25,31,36,43,49,57,63,68)+0.5, 
                 linetype = "dashed", color = "red", na.rm = TRUE) +
      scale_y_continuous(limits = c(y_min, y_max), 
                         breaks = seq(y_min, y_max, length.out = 3),
                         labels = function(x) sprintf("%.2f", x)) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            axis.title.x = element_blank())
  }
}
# First, remove y-axis titles from all plots except the middle one
for (i in 1:8) {
  if (i != 4) {  # Assuming you want to keep the y-axis title on the 4th plot
    plot.list[[i]] <- plot.list[[i]] + 
      theme(axis.title.y = element_blank())
  } else {
    plot.list[[i]] <- plot.list[[i]] + 
      theme(axis.title.y = element_text(size = 14))
  }
}

# Then combine the plots
plots <- plot.list[[1]]/plot.list[[2]]/plot.list[[3]]/plot.list[[4]]/plot.list[[5]]/plot.list[[6]]/plot.list[[7]]/plot.list[[8]]
plots <- plots + plot_layout(ncol = 1, heights = c(1,1,1,1,1,1,1,1))

# Display the plot
plots

ggsave('./plot/signature/boxplot_total_signature_UCell.pdf', width = 11, height = 12)


library(dunn.test)
p_values <- data.frame(comparisons = NA, P.adjusted = NA, signature = NA)

pvalue.list <- list()
for (i in celltype_use){
  fill_color <- ifelse(grepl(i, levels(animal_meta$group)), cell.colors1[i], "white")
  
  used_cell <- paste0(str_replace_all(i,' ','_'),'_signature_UCell')
  UCell_score <- animal_meta[,c(used_cell,'group','Species')]
  
  colnames(UCell_score)
  colnames(UCell_score)[1] <- 'Score'
  for(j in suffix_order){
    tfd_data <- UCell_score[UCell_score$Species == j, ]
    
    kruskal_result <- kruskal.test(Score ~ group, data = tfd_data)
    
    print(kruskal_result)
   
    dunn_result <- dunn.test(tfd_data$Score, tfd_data$group, method = "bonferroni")
    p_values1 <- data.frame(comparisons = dunn_result$comparisons, P.adjusted = dunn_result$P.adjusted)
    p_values1$signature <- i

    p_values <- rbind(p_values, p_values1)
  }
  pvalue.list[[i]] <- p_values
}

pvalue.list <- do.call(rbind, pvalue.list)
pvalue.list <- pvalue.list[-1,]

openxlsx::write.xlsx(pvalue.list, "./file/signature.pvalue.xlsx")
saveRDS(animal_metadata,'./rds/UCell_score.rds')
animal_metadata <- readRDS('./rds/UCell_score.rds')

############### conserved marker heat plot ##########
top_marker <- openxlsx::read.xlsx('./file/conserved_marker_cross_species_orig.xlsx')

heat_marker <- top_marker[,-1]

for (col_name in colnames(heat_marker)[-1]) {
  heat_marker[[col_name]] <- ifelse(!is.na(heat_marker[[col_name]]), TRUE, FALSE)
}

summary_counts <- heat_marker %>%
  group_by(multi_annotation) %>%
  summarize(across(everything(), ~sum(. == TRUE, na.rm = TRUE)))

heatmap_data <- summary_counts[, -1]

row.names(heatmap_data) <- summary_counts$multi_annotation

colnames(heatmap_data) <- colnames(heatmap_data) %>% str_replace('.symbol','') %>% str_replace_all('_', ' ')

heatmap_data <- as.matrix(heatmap_data) %>% t()

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    round(heatmap_data[i, j], 2), 
    x, y,
    gp = gpar(
      fontsize = 10
    ))
}

pdf('./plot/signature/signature_marker_number.pdf', width = 5, height = 10)

Heatmap(heatmap_data, cluster_rows = T,
        heatmap_legend_param = list(title = 'Number',legend_direction = "horizontal"),
        show_heatmap_legend = T,
        cell_fun = cell_fun,
        width = unit(4.5, "cm"), height = unit(14, "cm"),
        row_labels = species)

dev.off()

########## Do heatmap ######
top_marker <- top_marker %>% group_by(multi_annotation) %>% slice_max(score, n = 10)
human <- sp.list[['Human']]
table(human$multi_annotation)
human <- sample_cells(human,group.by = 'multi_annotation',sample.size = 150) # in function.R
levels(as.factor(human$multi_annotation))

human$multi_annotation <- factor(human$multi_annotation,
                                 levels=c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','Platelets'))

top_marker <- top_marker %>%
  arrange(factor(multi_annotation, levels = c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','Platelets')))

heat_human <- top_marker[top_marker$multi_annotation %in%
                           c('B cells','T cells','NK cells','Monocytes','mDCs','pDCs','Platelets'), ]$Human.gene.name

DoHeatmap(human,features = heat_human, slot = 'data', disp.min = 0, disp.max = 1, group.colors = cell.colors1)+
  scale_fill_gradient(low = 'white',high = 'black')

ggsave('./plot/signature/human_heatmap.pdf',height = 10,width = 5)

mouse <- sp.list[['Mouse']]
table(mouse$multi_annotation)
mouse <- sample_cells(mouse,group.by = 'multi_annotation',sample.size = 150)

top_marker <- top_marker %>%
  arrange(factor(multi_annotation, levels = c('B cells','T cells','NK cells','Monocytes','pDCs','Neutrophils')))

heat_mouse <- top_marker[top_marker$multi_annotation %in%
                           c('B cells','T cells','NK cells','Monocytes','pDCs','Neutrophils'), ]$Mouse.gene.name

DoHeatmap(mouse,features = heat_mouse, slot = 'data', disp.min = 0, disp.max = 1, group.colors = cell.colors1)+
  scale_fill_gradient(low = 'white',high = 'black')

ggsave('./plot/signature/mouse_heatmap.pdf',height = 10,width = 5)

ggplot(animal_metadata[[10]], aes(x = UMAP_1, y = UMAP_2, fill = mDCs_signature_UCell)) +
  geom_density_2d() +
  labs(x = "UMAP_1", y = "UMAP_2", title = "Density Plot of UMAP (mDCs Signature)") +
  theme_minimal()


