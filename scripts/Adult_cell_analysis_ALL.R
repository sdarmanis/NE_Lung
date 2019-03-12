source("Adult_cell_analysis_FUNCTIONS.R")
setwd(dir = "/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/170305_reAnalysis")
# Import data  ------------------------------------------------------------
data <- read.table("CIRM_master_expression_count_table.tsv", header=T, row.names=1, sep='\t', check.names = F)
#save(file = "GBM_analysis_workspace.RData", list=ls())
# Subset cells of interest 
data <- data[,grep("zsGreen", colnames(data))]
###########################################################################

# Generate metadata file for cells of interest  ---------------------------
metadata <- as.data.frame(matrix(nrow=ncol(data), ncol=2))
row.names(metadata) <- colnames(data)
row.names(metadata) <- gsub("-zsGreen-1-","-zsGreen-1",row.names(metadata))
# Change cell names for the mis-labelled library
#####
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/Mismatched_submission/")
cor_demux <- read.csv("correct_demux.csv", header=T)
incor_demux  <- read.csv("incorrect_demux.csv", header=T)
# Correct the names of the cells in the metadata file 
metadata[,"initial_name"] <- row.names(metadata)
metadata[,"actual_name"] <- row.names(metadata)
# Given name -> incorrect barcode -> correct name 
for(i in 1:nrow(metadata)){
  # get cell name. That name was assigned to a wrong cell. 
  cell <- metadata$initial_name[i]
  # which index corresponds to that cell in the wrong sample sheet ? 
  index <- as.character(incor_demux[which(incor_demux[,2]==cell),"Index"])
  if(length(index)!=0){
  # Which cell got the actual index instead of the correct cell ? 
  metadata[i,"actual_name"] <-  as.character(cor_demux[which(as.character(cor_demux$Index)==index),"SampleID"])}
}
# Now we need to rename the data fields as well 
# Initial name in metadata matches data column names 
# REplace colnames with actual name 
colnames(data) <- metadata$actual_name
row.names(metadata) <- metadata$actual_name
# 
metadata <- metadata[,-c(3,4)]
colnames(metadata) <- c("plate", "well")
list <- strsplit(row.names(metadata), "-")
for(i in 1:length(list)){
  metadata$plate[i] <- list[[i]][2]
  metadata$well[i] <- list[[i]][1]
}
rm("list")
###########################################################################


# # Import STAR stats to metadata -------------------------------------------
# setwd(dir = "/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/New_analysis/")
# star_meta <- read.csv("../STAR_table_summary_MOUSE.txt", header=F, sep="\t", check.names = F, row.names=1)
# star_meta <- star_meta[grep("zsGreen", row.names(star_meta)),]
# # Subset to columns of interest 
# star_meta <- star_meta[,c(6,7,10,25,27,29,30,31)]
# colnames(star_meta) <- c("Total_reads", "Average_input_length","Mapped_reads","multi-mapping_reads","too-multi-mapping_reads",
#                          "unmapped_reads_mismatch","unmapped_reads_short", "unmapped_reads_other")
# # Add info to metadata 
# metadata <- metadata[order(row.names(metadata)),]
# star_meta <- star_meta[row.names(metadata),]
# metadata <- cbind(metadata, star_meta)
# rm(star_meta)
###########################################################################

# Populate metadata fields ------------------------------------------------
# Remove stat rows that contain the special character "_"
data_stats <- data[grep("no_feature", row.names(data)):nrow(data),]
data_genes <- data[-c(grep("no_feature", row.names(data)):nrow(data)),]
# Remove ERCC rows 
ercc_rows <- grep("ERCC-", row.names(data_genes))
data_ercc <- data_genes[ercc_rows,]
data_genes <- data_genes[-ercc_rows,]
# Calculate fraction of ERCC reads 
metadata[,"ERCC_reads"] <- colSums(data_ercc) 
metadata[,"Non_ERCC_reads"] <- colSums(data_genes) 
metadata[,"ERCC_to_non_ERCC"] <- metadata[,"ERCC_reads"]/metadata[,"Non_ERCC_reads"]
# Number of genes detected
metadata[,"Genes_detected"] <- NA
for(i in 1:ncol(data_genes)){
  a <- which(row.names(metadata)==colnames(data_genes)[i])
  metadata$Genes_detected[a] <- length(which(data_genes[,i] > 1))
}
# Add age field to metadata
metadata[,"Age"] <- NA 
metadata[grep("10005001", metadata$plate),"Age"] <- "E14"
metadata[grep("10005002", metadata$plate),"Age"] <- "E14"
metadata[grep("10005003", metadata$plate),"Age"] <- "E17.5"
metadata[grep("10005004", metadata$plate),"Age"] <- "E19.5"
metadata[grep("10005005", metadata$plate),"Age"] <- "E15"
metadata[grep("10005006", metadata$plate),"Age"] <- "PN21"
metadata[grep("1000500608", metadata$plate),"Age"] <- "PN120"
metadata[grep("1000500609", metadata$plate),"Age"] <- "PN120"
metadata[grep("10005007", metadata$plate),"Age"] <- "PN120"
metadata[grep("10005008", metadata$plate),"Age"] <- "PN90"
metadata$Age <- as.factor(metadata$Age)
# And color 
metadata[,"Age_color"] <- tableau_color_pal("tableau20")(20)[metadata$Age]
###########################################################################


# Subset adult cells ------------------------------------------------------
meta_adult <- metadata[grep("PN", metadata$Age),]
meta_adult[,"plate_color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$plate)]
# Change names 0f data to match metadata
colnames(data_genes) <- gsub("-zsGreen-1-","-zsGreen-1",colnames(data_genes))
data_genes_norm <- apply(X=data_genes, MARGIN=2, function(x) log2((((x/sum(x))*1000000)+1)) )
# Subset adult data and normalized data 
data_adult <- data_genes[,row.names(meta_adult)]
data_norm_adult <- data_genes_norm[,row.names(meta_adult)]
# Number of genes detected
meta_adult[,"Genes_detected"] <- NA
for(i in 1:ncol(data_adult)){
  a <- which(row.names(meta_adult)==colnames(data_adult)[i])
  meta_adult$Genes_detected[a] <- length(which(data_adult[,i] != 0))
}
# Subset fetal cells 
meta_fetal <- metadata[grep("E", metadata$Age),]
data_fetal <- data_genes[,row.names(meta_fetal)]
data_norm_fetal <- data_genes_norm[,row.names(meta_fetal)]
save(list = c("meta_fetal", "data_fetal", "data_norm_fetal"), file="~/Desktop/Fetal_lung_data_Christin.RData")

###########################################################################


# Cell filtering  ---------------------------------------------------------
setwd(dir = "/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/170305_reAnalysis")
pdf("Filtering_and_QC_plots.pdf",15,15)
# Remove cells with NAs and all zero reads 
data_adult <- data_adult[,names(which(is.na(colSums(data_adult))==F))]
data_adult <- data_adult[,which(colSums(data_adult)!=0)]
meta_adult <- meta_adult[colnames(data_adult),]
# Remove cells with too small or too high number of genes 
hist(meta_adult$Genes_detected, 100, main="Genes detected per cell")
abline(v=500, col="red")
abline(v=5000, col="red")
# Remove cells 
cells <- row.names(meta_adult)[which(meta_adult$Genes_detected > 500)]
temp <- meta_adult[cells,]
cells <- row.names(temp)[which(temp$Genes_detected < 5000)]
meta_adult <- meta_adult[cells,]
data_adult <- data_adult[,row.names(meta_adult)]
data_norm_adult <- data_genes_norm[,row.names(meta_adult)]
rm(temp)
# Also remove cells with less than 200000 reads 
hist(meta_adult$Non_ERCC_reads, 100, main="Mapped reads per cell")
abline(v=200000, col="red")
# Remove cells 
cells <- row.names(meta_adult)[which(meta_adult$Non_ERCC_reads > 200000)]
meta_adult <- meta_adult[cells,]
data_adult <- data_adult[,row.names(meta_adult)]
data_norm_adult <- data_genes_norm[,row.names(meta_adult)]
dim(meta_adult)
# # Filter cells for housekeeping genes
# color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
# housekeeping <- read.csv("housekeeping_genes.csv", header=T, sep=";")
# heatmap.2(data_norm_adult[as.character(housekeeping[,1]),], trace="n", main="All cells", cexCol=0.5,  col=color.palette, hclustfun=function(x) hclust(d=x, method='ward.D2'), 
#           ColSideColors=meta_adult$plate)
# # Cluster based on housekeeping genes 
# dist_house <- dist(t(data_norm_adult[as.character(housekeeping[,1]),]), method='euclidean')
# clust_house <- hclust(dist_house, method='ward.D2')
# k <- 2
# tree_house <- cutree(clust_house,k=k)
# meta_adult[,"housekeeping_cluster"] <- tree_house
# meta_adult[,"housekeeping_cluster"] <- as.factor(meta_adult[,"housekeeping_cluster"])
# # Generate color vector for housekeeping genes clustering 
# meta_adult[,"housekeeping_cluster_color"] <- tableau_color_pal("tableau20")(20)[tree_house]
# # Replot heatmap with cluster colors 
# heatmap.2(data_norm_adult[as.character(housekeeping[,1]),], trace="n", main="All cells", cexCol=0.5,  col=color.palette, hclustfun=function(x) hclust(d=x, method='ward.D2'), 
#           ColSideColors=meta_adult$housekeeping_cluster_color)
# Plot QC graphs
plot(meta_adult$Non_ERCC_reads,meta_adult$Genes_detected, xlab="Mapped reads", ylab="Genes detected", 
     pch=19)
# # 
# ggplot(meta_adult, aes(Genes_detected, fill=housekeeping_cluster)) +   geom_density(alpha = 0.2) + theme_bw()
# ggplot(meta_adult, aes(Total_reads, fill=housekeeping_cluster)) +   geom_density(alpha = 0.2) + theme_bw()
# Subset cells belonging to the 'good' cluster and proceed with the analysis
# meta_adult <- meta_adult[which(meta_adult$housekeeping_cluster!=2),]
# data_adult <- data_adult[,row.names(meta_adult)]
# data_norm_adult <- data_norm_adult[,row.names(meta_adult)]
save(file = "Kuo_workspace_LATEST_postQC.RData", list=ls())
dev.off()
####################################



# Overdispersion and TSNE analysis ----------------------------------------
o.adult <- sel.by.cv(data_adult)
data_norm_top_adult <- data_norm_adult[o.adult[1:1000],]
# Remove cells with no expression of any of the overdispersed genes 
cells <- names(which(colSums(data_norm_top_adult)!=0))
meta_adult <- meta_adult[cells,]
data_adult <- data_adult[,row.names(meta_adult)]
data_norm_adult <- data_norm_adult[,row.names(meta_adult)]
# Resubset on overdispersed genes 
# o.adult <- sel.by.cv(data_adult)
# data_norm_top_adult <- data_norm_adult[o.adult[1:1000],]
# Calculate correlation distance 
dist.cor.adult <- as.dist(1-abs(cor(data_norm_top_adult)))
rm(data_norm_top_adult)
# Perform TSNE 
set.seed(123)
my.tsne.2d.adult <- run.tsne(my.dist=dist.cor.adult,k=2, perplexity=32, max_iter=5000) # Last perplexity used = 30
row.names(my.tsne.2d.adult) <- row.names(as.matrix(dist.cor.adult))
# Make sure metadata and counts and TSNE coordinates are ordered
meta_adult <- meta_adult[row.names(my.tsne.2d.adult),]
data_adult <- data_adult[,row.names(meta_adult)]
data_norm_adult <- data_norm_adult[,row.names(meta_adult)]
colnames(my.tsne.2d.adult) <- c("Dim1", "Dim2")
###########################################################################


# Clustering with DBscan --------------------------------------------------
library("fpc")
library("dbscan")
library("factoextra")
# Define the optimal eps value for dbscvan 
dbscan::kNNdistplot(my.tsne.2d.adult, k =  5)
abline(h=c(1,2,3,4,5), lty=2, col="red")
eps <- 3
set.seed(123)
db <- fpc::dbscan(my.tsne.2d.adult, eps = eps, MinPts = 5)
# Plot DBSCAN results
plot(db, my.tsne.2d.adult, main = paste("eps=",eps,"Nclus=",length(unique(db$cluster)),sep=""), frame = FALSE)
# Plot again 
fviz_cluster(db, my.tsne.2d.adult, stand = FALSE, frame = FALSE, geom = "point")
# Add cluster info to metadata
clusters <- db$cluster
meta_adult[,"Cluster_2d"] <- clusters
meta_adult[,"Cluster_2d"] <- as.numeric(meta_adult[,"Cluster_2d"])
meta_adult[,"Cluster_color_2d"]<- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$Cluster_2d)]
# Change cluster color for cluster 0 (unassigned) to black 
meta_adult[which(meta_adult$Cluster_2d==0),"Cluster_color_2d"] <- "grey80"
# And replot TSNE 
plot(my.tsne.2d.adult,  col = meta_adult$Cluster_color_2d, pch=19, cex=0.7)
legend('bottomleft', legend = unique(meta_adult$Cluster_2d), fill=unique(meta_adult$Cluster_color_2d), 
       cex=0.7,ncol=2)
# Remove unnasigned cells from dataset 
unassigned_cells <- row.names(meta_adult)[which(meta_adult$Cluster_2d==0)]
assigned_cells <- row.names(meta_adult)[which(meta_adult$Cluster_2d!=0)]
meta_adult <- meta_adult[assigned_cells,]
my.tsne.2d.adult <- my.tsne.2d.adult[row.names(meta_adult),]
# Also remove non-expressed genes 
genes <- names(which(rowSums(data_adult)!=0))
data_adult <- data_adult[genes,row.names(meta_adult)]
data_norm_adult <- data_norm_adult[genes,row.names(meta_adult)]
###########################################################################

# Crude dif.exp using Wilcoxon to identify cell types  --------------------
mat_wilcox_clus_between <- matrix(nrow=nrow(data_norm_adult), ncol=max(unique(meta_adult$Cluster_2d)))
row.names(mat_wilcox_clus_between) <- row.names(data_norm_adult)
colnames(mat_wilcox_clus_between) <- paste("Cluster", unique(meta_adult$Cluster_2d), sep="_")
for(i in unique(meta_adult$Cluster_2d)){
  print(paste("Cluster is ", i))
  names_1 <- row.names(meta_adult)[which(meta_adult$Cluster_2d==i)]
  names_2 <- row.names(meta_adult)[which(meta_adult$Cluster_2d!=i)]
  for(j in 1:nrow(data_norm_adult)){
    cat("\r",paste("Genes done % is", j/nrow(data_norm_adult)*100))
    test <- wilcox.test(data_norm_adult[j,names_1],data_norm_adult[j,names_2])
    mat_wilcox_clus_between[j,i] <- test$p.value*nrow(mat_wilcox_clus_between)
  }
}
# 
write.table(mat_wilcox_clus_between, "Wilcoxon_between_subclusters_ADULT.csv")
###########################################################################

# Crude dif.exp using ROC analysis to identify cell types  -------------------- 
meta_adult[,"Cluster_2d_character"] <- paste("cluster_", meta_adult$Cluster_2d, sep="")
mat_auc_clus_between <- matrix(nrow=nrow(data_norm_adult), ncol=max(unique(meta_adult$Cluster_2d)))
row.names(mat_auc_clus_between) <- row.names(data_norm_adult)
colnames(mat_auc_clus_between) <- paste("Cluster", unique(meta_adult$Cluster_2d), sep="_")
for(i in 1:length(unique(meta_adult[,"Cluster_2d_character"]))){
  print(paste("Doing", unique(meta_adult[,"Cluster_2d_character"])[i]))
  # Vector of cell type ID for ROC (1:type of interest, 0:all else)
  class_vec <- meta_adult$Cluster_2d_character
  clus <- unique(meta_adult[,"Cluster_2d_character"])[i]
  class_vec[class_vec==clus] <- 1
  class_vec[class_vec!=1] <- 0
  class_vec <- factor(class_vec)
  # Also separate cells into cells_1 and cells_2 (want/dont want)
  cells_1 <- row.names(meta_adult)[which(meta_adult$Cluster_2d_character==clus)]
  cells_2 <- row.names(meta_adult)[which(meta_adult$Cluster_2d_character!=clus)]
 # AUC calculation 
   for(k in 1:nrow(data_norm_adult)){
    pred <- prediction(data_norm_adult[k,], class_vec)
    perf <- performance(prediction.obj = pred, measure = "auc")
    mat_auc_clus_between[k,i] <- slot(perf, "y.values")[[1]]

  }
}
# 
write.table(mat_auc_clus_between, "AUC_between_subclusters_ADULT.csv")
###########################################################################  
  

pdf("TSNE_Adult_cells_with_gene_enrichment.pdf",10,10)
# Plot TSNE and cluster enriched genes ------------------------------------
plot(my.tsne.2d.adult, col=as.character(meta_adult$Cluster_color_2d), pch=19, xlab='', ylab='', cex=0.5,main="TSNE clusters")
legend('topright', legend=unique(meta_adult$Cluster_2d), fill=unique(meta_adult$Cluster_color_2d), 
       cex=0.7, ncol=2)
# Plot tsne with Age
plot(my.tsne.2d.adult, col=as.character(meta_adult$Age_color), pch=19, xlab='', ylab='', cex=0.5,main="Age")
legend('topright', legend=unique(meta_adult$Age), fill=unique(meta_adult$Age_color), 
       cex=0.7, ncol=2)
# Plot tsne with plate
plot(my.tsne.2d.adult, col=tableau_color_pal("tableau20")(20)[as.factor(meta_adult$plate)], pch=19, 
     xlab='', ylab='', cex=0.5,main="Plate")
legend('topright', legend=unique(meta_adult$plate), fill=unique(tableau_color_pal("tableau20")(20)[as.factor(meta_adult$plate)]), 
       cex=0.7, ncol=2)
# Plot the top 5 enriched genes per subcluster 
par(mfcol=c(5,5), mar=c(2,2,2,2))
for(i in sort(unique(meta_adult$Cluster_2d))){
  names <- names(sort(mat_wilcox_clus_between[,i]))[1:25]
  for(j in 1:25){
    boxplot(data_norm_adult[names[j],]~meta_adult$Cluster_2d, outline=F, boxwex=0.5, main=names[j])
    stripchart(data_norm_adult[names[j],]~meta_adult$Cluster_2d,cex=0.5,
               vertical = TRUE, method = "jitter", 
               pch = 21, col = "maroon", bg = "bisque", 
               add = TRUE) 
  }
  names <- names(sort(mat_auc_clus_between[,i],decreasing = T))[1:25]
  for(j in 1:25){
    boxplot(data_norm_adult[names[j],]~meta_adult$Cluster_2d, outline=F, boxwex=0.5, main=names[j])
    stripchart(data_norm_adult[names[j],]~meta_adult$Cluster_2d,cex=0.5,
               vertical = TRUE, method = "jitter", 
               pch = 21, col = "maroon", bg = "bisque", 
               add = TRUE) 
  }
  }
par(mfcol=c(1,1), mar=c(4,4,4,4))
boxplot(meta_adult$Genes_detected~meta_adult$Cluster_2d, main="Number of genes per cluster", xlab="TSNE cluster", ylab="Genes detected")
# Heatmap of top enriched genes (by Wilcoxon)
names <- c()
for(i in sort(unique(meta_adult$Cluster_2d))){
  a <- names(sort(mat_wilcox_clus_between[,i]))[1:10]
  names <- c(names,a)
}
my.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
ord <- order(meta_adult$Cluster_2d)
#pdf("Heatmap_top_10_enriched_genes_by_wilcox.pdf",15,15)
heatmap.2(data_norm_adult[names,ord], trace="n", Colv=F, ColSideColors = meta_adult$Cluster_color_2d[ord], col=my.palette, 
                  hclustfun = function(x) hclust(x, "ward.D2"))
# dev.off()
# And a table with cells per cluster per plate 
table_1 <- table(meta_adult$plate, meta_adult$Cluster_2d)
colnames(table_1) <- paste("Cluster_", colnames(table_1), sep="")
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table_1),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

dev.off()
###########################################################################

# Add cell type info to metadata ------------------------------------------
# This is custom and needs to be adjusted to new data, etc.. 
meta_adult[,"cell.type"] <- meta_adult$Cluster_2d
# These vectors needs to be matched (cluster number to cluster cell type) with the one before it, otherwise shit will happen .. 
vec_1 <- c(1:max(meta_adult$Cluster_2d))
vec_2 <- c("NE", "Club", "Glial", "AT2", "Fibroblasts", "AT1", "Basal","S.Muscle", "Ciliated", "Endothelial") 
for(i in 1:length(vec_1)){
  meta_adult$cell.type[which(meta_adult$cell.type==vec_1[i])] <- vec_2[i]
}
# Convert to colors
meta_adult[,"cell.type.color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$cell.type)]
###########################################################################

# Append FACS data to metadata --------------------------------------------
# Add Index Sort information to metadata 
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/PN FACS data-selected/Index_files/")
# List of all files to read 
files <- list.files()
files <- files[grep(".csv", files)]
# Read all files
list_index <- list()
for (f in 1:length(files)) {
  tempData = read.csv(files[f], header=T, skip=15)
  tempData = apply(tempData, MARGIN = 2, function(x) gsub(",", "", x))
  list_index[[f]] <- tempData
  names(list_index)[[f]]  <- files[f]
} 
# Append plate ID on the row names of each file and standarize number of columns
for(i in 1:length(list_index)){
  temp <- list_index[[i]]
  name <- names(list_index)[i]
  name <- strsplit(name, split = ".", fixed = T)
  barcode <- name[[1]][1]
  row.names(temp) <- paste(temp[,1],"-",barcode,"-zsGreen","-1-",sep="")
  # Keep only columns of interest (containing the string ".Mean")
  temp <- temp[,grep(".Mean",colnames(temp))]
  # Change column names 
  # First remove column Time 
  a <- grep("Time", colnames(temp))
  if(length(a) != 0) {temp <- temp[,-a]}
  colnames(temp) <- c("FSC", "SSC", "APC_Lineage", "DAPI", "PE-Cy7_Epcam", "zsGreen")
  list_index[[i]] <- temp
}
# Unlist and combine in one large file 
mat_index <- do.call(rbind, list_index)
mat_index <- as.data.frame(mat_index)
# Convert columns to numeric
for(i in 1:ncol(mat_index)){
  mat_index[,i] <- as.numeric(as.character(mat_index[,i]))
}
# Remove cells with NA values
mat_index <- mat_index[which(is.na(rowSums(mat_index))==F),]
# Write table for future use
write.table(mat_index, file="Index_data_aggregated.tsv")
# Add index data to metadata and plot 
mat_index_sub <- mat_index[row.names(meta_adult),]
# Append metadata (for adult cells) and index data 
meta_adult <- cbind(meta_adult, mat_index_sub)
# Reset the wd 
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/170305_reAnalysis/")
###########################################################################

pdf("FACS_index_data_plot.pdf",15,15)
# FACS data plot and gating -----------------------------------------------
meta_adult[,"plate.color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$plate)]
# General plots and drawing of cutoffs for EpCam and zsGreen 
lim_1 <- c(0,8)
plate_1 <- grep("10005006", meta_adult$plate)
plate_2 <- grep("10005007", meta_adult$plate)
plate_3 <- grep("10005008", meta_adult$plate)
# Set level for Epcam and zsGreen 
epcam_cut_1 <- 2 # For plates "6"
epcam_cut_2 <- 3.6 # For plates "7"
epcam_cut_3 <- 4 # For plates "8"
green_cut_1 <- 2.3 # For plates "6"
green_cut_2 <- 4 # For plates "7"
green_cut_3 <- 3.5 # For plates "8"

# Plot all cells per "plate"
# Plate 6
plot(log10(meta_adult$zsGreen[plate_1]),log10(meta_adult$`PE-Cy7_Epcam`[plate_1]), pch=19,xlim=lim_1, ylim=lim_1,
     xlab="zsGreen (FACS)", ylab="EpCam (FACS)", main="Plate 6", cex=0.8, col=meta_adult$plate.color[plate_1])
legend("bottomright", legend=unique(meta_adult$plate[plate_1]), fill=unique(meta_adult$plate.color[plate_1]))
abline(h=epcam_cut_1, col="black")
abline(v=green_cut_1, col="black")
# Annotate different groups 
text(1,6,"zsGreen(low)/EpCam(+)", cex=0.8)
text(1,0,"zsGreen(low)/EpCam(-)", cex=0.8)
text(6,6,"zsGreen(high)/EpCam(+)", cex=0.8)
text(6,0,"zsGreen(high)/EpCam(-)", cex=0.8)
# Plate 7
plot(log10(meta_adult$zsGreen[plate_2]),log10(meta_adult$`PE-Cy7_Epcam`[plate_2]), pch=19,xlim=lim_1, ylim=lim_1,
     xlab="zsGreen (FACS)", ylab="EpCam (FACS)", main="Plate 7", cex=0.8, col=meta_adult$plate.color[plate_2])
legend("bottomright", legend=unique(meta_adult$plate[plate_2]), fill=unique(meta_adult$plate.color[plate_2]))
abline(h=epcam_cut_2, col="black")
abline(v=green_cut_2, col="black")
# Annotate different groups 
text(1,6,"zsGreen(low)/EpCam(+)", cex=0.8)
text(1,0,"zsGreen(low)/EpCam(-)", cex=0.8)
text(6,6,"zsGreen(high)/EpCam(+)", cex=0.8)
text(6,0,"zsGreen(high)/EpCam(-)", cex=0.8)
# Plate 8
plot(log10(meta_adult$zsGreen[plate_3]),log10(meta_adult$`PE-Cy7_Epcam`[plate_3]), pch=19,xlim=lim_1, ylim=lim_1,
     xlab="zsGreen (FACS)", ylab="EpCam (FACS)", main="Plate 8", cex=0.8, col=meta_adult$plate.color[plate_3])
legend("bottomright", legend=unique(meta_adult$plate[plate_3]), fill=unique(meta_adult$plate.color[plate_3]))
abline(h=epcam_cut_3, col="black")
abline(v=green_cut_3, col="black")
# Annotate different groups 
text(1,6,"zsGreen(low)/EpCam(+)", cex=0.8)
text(1,0,"zsGreen(low)/EpCam(-)", cex=0.8)
text(6,6,"zsGreen(high)/EpCam(+)", cex=0.8)
text(6,0,"zsGreen(high)/EpCam(-)", cex=0.8)
# Plot again with plate colors
plates <- unique(meta_adult$plate)
par(mfcol=c(3,3),mar=c(5,7,3,7))
for(i in plates){
  cells <- which(meta_adult$plate==i)
  plot(log10(meta_adult$zsGreen[cells]),log10(meta_adult$`PE-Cy7_Epcam`[cells]), pch=19,xlim=lim_1, ylim=lim_1,
       xlab="zsGreen (FACS)", ylab="EpCam (FACS)", main=i, cex=0.6, col=meta_adult$plate.color[cells])
  # Set level for Epcam and zsGreen 
  if(length(grep(10005006, i)) > 0)  {abline(h=epcam_cut_1, col="black") ; abline(v=green_cut_1, col="black")}
  if(length(grep(10005007, i)) > 0)  {abline(h=epcam_cut_2, col="black") ; abline(v=green_cut_2, col="black")}
  if(length(grep(10005008, i)) > 0)  {abline(h=epcam_cut_3, col="black") ; abline(v=green_cut_3, col="black")}
}
# Add annotation to metadata 
# For zsGreen and Epcam using different cutoffs per plate  
meta_adult[,"zsGreen_class"] <- "low"
meta_adult[,"Epcam_class"] <- "low"
for(i in 1:nrow(meta_adult)){
  name <- row.names(meta_adult)[i]
  if(is.na(meta_adult[i,"zsGreen"])==F){
  if(length(grep(10005006, name)) > 0) {
    if(meta_adult$zsGreen[i] > 10^green_cut_1) {meta_adult[i,"zsGreen_class"] <- "high"}
    if(meta_adult$`PE-Cy7_Epcam`[i] > 10^epcam_cut_1) {meta_adult[i,"Epcam_class"] <- "high"}}
  if(length(grep(10005007, name)) > 0) {
    if(meta_adult$zsGreen[i] > 10^green_cut_2) {meta_adult[i,"zsGreen_class"] <- "high"}
    if(meta_adult$`PE-Cy7_Epcam`[i] > 10^epcam_cut_2) {meta_adult[i,"Epcam_class"] <- "high"}}
  if(length(grep(10005008, name)) > 0) {
    if(meta_adult$zsGreen[i] > 10^green_cut_3) {meta_adult[i,"zsGreen_class"] <- "high"}
    if(meta_adult$`PE-Cy7_Epcam`[i] > 10^epcam_cut_3) {meta_adult[i,"Epcam_class"] <- "high"}}}
}
# Add colors 
meta_adult[,"zsGreen_class_color"] <- "black"
meta_adult[which(meta_adult$zsGreen_class=="high"), "zsGreen_class_color"] <- "red" 
meta_adult[,"Epcam_class_color"] <- "black"
meta_adult[which(meta_adult$Epcam_class=="high"), "Epcam_class_color"] <- "red" 
# Create and plot table 
# Add boxplots per TSNE cluster 
# zsGreen
par(mfcol=c(1,1), mar=c(8,4,4,4), las=3)
box.with.strip(log10(meta_adult$zsGreen), group = meta_adult$Cluster_2d, main="zsGreen per TSNE cluster", xlab="", ylab="")
box.with.strip(log10(meta_adult$zsGreen[which(meta_adult$sum_class==1)]), group=meta_adult$final_class[which(meta_adult$sum_class==1)],
        main="zsGreen per final cluster", xlab="", ylab="")
# Plot table TSNE vs zsGreen
table_1 <- table(meta_adult$Cluster_2d, meta_adult$zsGreen_class)
row.names(table_1) <- paste("TSNE", row.names(table_1), sep=".")
colnames(table_1) <-  paste("zsGreen", colnames(table_1), sep=".")
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table_1),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
# Epcam
boxplot(log10(meta_adult$`PE-Cy7_Epcam`)~meta_adult$Cluster_2d, main="Epcam per TSNE cluster")
# Plot table TSNE vs zsGreen
table_1 <- table(meta_adult$Cluster_2d, meta_adult$Epcam_class)
row.names(table_1) <- paste("TSNE", row.names(table_1), sep=".")
colnames(table_1) <-  paste("Epcam", colnames(table_1), sep=".")
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table_1),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
# Also plot TSNE with Epcam and zsGreen colors 
pal=colorRampPalette(c("grey95","blue","red"))
colfunc <- colorRampPalette(c("grey95","blue","yellow"))
my_palette <- colfunc(100)
ncol=100
# ZsGreen 
col <- val2col(log10(meta_adult$zsGreen), col = pal(ncol))
layout(matrix(nrow=1, ncol=2, c(1:2)), widths = c(0.85,0.15))
plot(my.tsne.2d.adult, col=col,pch=19, xlab='', ylab='', cex=0.7, main="zsGreen")
image.scale(log10(meta_adult$zsGreen), col=pal(ncol), axis.pos = 2)
# Epcam
vec <- meta_adult$`PE-Cy7_Epcam`
vec[vec < 1] <- 1
vec <- log10(vec)
col <- val2col(vec, col = pal(ncol))
layout(matrix(nrow=1, ncol=2, c(1:2)), widths = c(0.85,0.15))
plot(my.tsne.2d.adult, col=col,pch=19, xlab='', ylab='', cex=0.7, main="Epcam")
image.scale(vec, col=pal(ncol), axis.pos = 2)
# 
dev.off()
###########################################################################

# Gene markers for each cell type -----------------------------------------
# loop thorugh each individual cell type and compare it to all other cell types 
# find genes that help distinguish each cell type from ALL other cell types 
# All cell types 
# types <- unique(meta_adult$cell.type)
# # Create list to store gene matrix for each cell type 
# list_markers <- list()
# #
# for(i in types){
#   print(paste("Doing", i))
#   # Vector of cell type ID for ROC (1:type of interest, 0:all else)
#   class_vec <- meta_adult$cell.type
#   class_vec[class_vec==i] <- 1
#   class_vec[class_vec!=1] <- 0
#   class_vec <- factor(class_vec)
#   # Also separate cells into cells_1 and cells_2 (want/dont want)
#   cells_1 <- which(meta_adult$cell.type==i)
#   cells_2 <- which(meta_adult$cell.type!=i)
#   # Create and label matrix to store values for each gene
#   mat_marker_temp <- matrix(nrow=nrow(data_norm_adult), ncol=3)  
#   row.names(mat_marker_temp) <- row.names(data_norm_adult)  
#   colnames(mat_marker_temp) <- c("AUC", "Fraction_cell", "Fraction_non_cell")
#   # Go through each gene and recored AUC and fraction of expressing cells 
#   for(k in 1:nrow(data_norm_adult)){
#     pred <- prediction(data_norm_adult[k,], class_vec)  
#     perf <- performance(prediction.obj = pred, measure = "auc")
#     mat_marker_temp[k,1] <- slot(perf, "y.values")[[1]]
#     a <- length(which(data_norm_adult[k,cells_1] > 0))/length(cells_1)
#     if(length(a) > 0) {mat_marker_temp[k,2] <- a} else {mat_marker_temp[k,2] <- 0}
#     b <- length(which(data_norm_adult[k,cells_2] > 0))/length(cells_2)
#     if(length(a) > 0) {mat_marker_temp[k,3] <- b} else {mat_marker_temp[k,3] <- 0}
#   }
#   # Add mat_marker_temp to list_markers
#   list_markers[[which(types==i)]] <- mat_marker_temp
# }
# # Add names to list_markers
# names(list_markers) <- types
# # Go through the list and find genes for each cell type that fulfil the following criteria 
# # Fraction cell > 0.8, Fraction non cell < 0.2, AUC > 0.75
# # Store top genes in list_markers_names
# list_markers_names <- list()
# # 
# for(i in 1:length(list_markers)){
#   mat_temp <- list_markers[[i]]
#   # Select the ones with high AUC values (> 0.75)
#   genes <- names(which(mat_temp[,"AUC"] > 0.75))
#   # Select the genes expressed in more than 80% of cells of interest
#   genes <- which(mat_temp[,"Fraction_cell"] > 0.8)
#   # Select the genes expressed in less than 20% of cells of non interest
#   genes <- names(which(mat_temp[genes,"Fraction_non_cell"] < 0.2))
#   # Import gene list for canonical markers 
#   list_markers_names[[i]] <- genes 
# }
# names(list_markers_names) <- names(list_markers)
###############################################################################################


# Gene markers for each cell type v2 --------------------------------------
# loop thorugh each individual cell type and compare it to each other cell type
# All cell types 
types <- unique(meta_adult$cell.type)
# Matrix to store fraction of expression per gene and per cell type 
mat_markers <- matrix(nrow=nrow(data_norm_adult), ncol=length(types))
row.names(mat_markers) <- row.names(data_norm_adult)
colnames(mat_markers) <- types 
for(i in 1:nrow(data_norm_adult)){
  for(j in 1:length(types)){
    cells <- row.names(meta_adult)[which(meta_adult$cell.type==types[j])]
    mat_markers[i,j] <-  length(which(data_norm_adult[i,cells]!=0))/length(cells)
  }
}
# Go through the list and find genes for each cell type that fulfil the following criteria 
# Fraction cell > 0.85, Max Fraction non cell < 0.35
# Store top genes in list_markers_names
list_markers_names <- list()
for(j in 1:length(types)){
  # Higher than 80% of the correct type 
  pass_1 <- names(which(mat_markers[,types[j]] > 0.85))
  mat_temp <- mat_markers[pass_1,which(colnames(mat_markers)!=types[j])]
  # Find genes with a max fraction of < 0.2 in any other cell type 
  list_markers_names[[j]] <- names(which(apply(mat_temp, MARGIN = 1, mean) < 0.15))
}
names(list_markers_names) <- types
###############################################################################################


# Classify cells based on new specific gene list  -------------------------
# Start with drawing heatmaps for each cell type 
# Plot 
pdf("Extended_marker_list_all_cell_types_UNCUT.pdf",10,10)
for(i in names(list_markers_names)){
# Colors for overlaying on heatmap
col_labels <- cbind(meta_adult$Cluster_color_2d, meta_adult$cell.type.color, meta_adult$cell.type)
col_vec <- col_labels[,3]
col_labels[which(col_vec==i),3] <- "red"
col_labels[which(col_vec!=i),3] <- "black"
colnames(col_labels) <- c("TSNE.color", "Cell.type.color", "Cell.type.shown")
# Get genes and plot 
genes <- list_markers_names[[i]]
heatmap.3(data_norm_adult[genes,], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main=i, 
          labCol = "", cexRow = 0.9, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=7)    
# Add
legend("topright",
       legend=c(unique(meta_adult$Cluster_2d), " ", unique(meta_adult$cell.type)," ", i, paste("Not", i, sep="_")),
       fill=c(unique(meta_adult$Cluster_color_2d), "white", unique(meta_adult$cell.type.color), "white", "red", "black"), 
       border=FALSE, bty="n", y.intersp = 0.7, cex=0.6)
}
dev.off()
# View all heatmaps, decide on how to cut the cell dendrogram 
# Create vector of cuts
types_cuts <- as.data.frame(types)
types_cuts[,"cuts"] <- c(2,2,2,2,2,2,2,3,2) # Vector of tree cuts. Needs to change if trees are different.. 
# Replot all heatmaps with the right number of cuts 
# Look at the heatmaps and decide on the cluster number that contains the cells of interest s
list_cuts <- list()
pdf("Extended_marker_list_all_cell_types_CUT.pdf",10,10)
for(i in names(list_markers_names)){
  genes <- list_markers_names[[i]]
  # Colors for overlaying on heatmap
  col_labels <- cbind(meta_adult$Cluster_color_2d, meta_adult$cell.type.color, meta_adult$cell.type)
  col_vec <- col_labels[,3]
  col_labels[which(col_vec==i),3] <- "red"
  col_labels[which(col_vec!=i),3] <- "black"
  # Generate dendrogram, cut and convert to colors 
  # Store each cut result in a list called list_cuts 
  dist.temp  <- dist(t(data_norm_adult[genes,]), "euclidean")
  clust.temp <- hclust(dist.temp, 'ward.D2')
  cut.temp <- cutree(clust.temp, types_cuts$cuts[which(types_cuts$types==i)])
  list_cuts[[which(types_cuts$types==i)]]  <- cut.temp
  cut.temp.color <- tableau_color_pal("tableau20")(20)[cut.temp]
  col_labels <- cbind(col_labels, cut.temp.color)
  colnames(col_labels) <- c("TSNE.color", "Cell.type.color", "Cell.type.shown", "tree.cut")
  # Get genes and plot 
  genes <- list_markers_names[[i]]
  heatmap.3(data_norm_adult[genes,], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
            hclustfun = function(x) hclust(x, 'ward.D2'), main=i, 
            labCol = "", cexRow = 0.9, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=7)    
  # Add
  legend("topright",
         legend=c(unique(meta_adult$Cluster_2d), " ", unique(meta_adult$cell.type)," ", i, paste("Not", i, sep="_")),
         fill=c(unique(meta_adult$Cluster_color_2d), "white", unique(meta_adult$cell.type.color), "white", "red", "black"), 
         border=FALSE, bty="n", y.intersp = 0.7, cex=0.6)
}
dev.off()
names(list_cuts) <- names(list_markers_names)
# Inspect trees with cluster numbers 
pdf("Classification_of_all_cells_using_specific_gene_lists.pdf",10,10)
# Decide on which cluster to keep for each cell type 
types_cuts[,"cluster_to_keep"] <- c(1,2,2,2,2,2,2,3,2)
# The clustering results are stored in list_cuts
# Go through the list and flag each cell based on the classfication using the new markers 
mat_class <- matrix(nrow=nrow(meta_adult), ncol=nrow(types_cuts))
row.names(mat_class) <- row.names(meta_adult)
colnames(mat_class) <- types_cuts$types
for(i in 1:length(list_cuts)){
  temp <- list_cuts[[i]]
  temp[which(temp==types_cuts$cluster_to_keep[i])] <- names(list_cuts)[i]
  temp[which(temp!=names(list_cuts)[i])] <- "No.clus" 
  mat_class[,i] <-   temp
}
# Convert mat_class to numerical (1: correct classification, 0:no classification)
mat_class_num <- mat_class
mat_class_num[mat_class_num!="No.clus"] <- 1
mat_class_num[mat_class_num=="No.clus"] <- 0
mat_class_num <- apply(mat_class_num, MARGIN = 2, as.numeric)
row.names(mat_class_num) <- row.names(meta_adult)
colnames(mat_class_num) <- types_cuts$types
# Plot each cell's score flag cells and add info to meta_adult 
plot(rowSums(mat_class_num), pch=19, col=meta_adult$Cluster_color_2d, main="TSNE colors")
# Add sum (basically a count of memberships) and final cell.type annotation to meta_adult
meta_adult[,"sum_class"] <- rowSums(mat_class_num)
meta_adult[,"final_class"] <- NA
for(i in 1:nrow(mat_class)){
  vec <- mat_class[i,]
  vec[grep("No.clus",vec)] <- ""
    meta_adult[i,"final_class"] <- paste(vec, collapse = "", sep="_")
}
# Add annotation for rowSums of 0 
meta_adult[which(meta_adult[,"sum_class"]==0),"final_class"] <- "None"
meta_adult$sum_class.color <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$sum_class)]
meta_adult[,"final_class_color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$final_class)]
table <- table( meta_adult[,"final_class"],meta_adult[,"sum_class"])
# Plot table 
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
# Plot classification barplots 
par(las=3, mar=c(15,2,2,2))
for(i in c(1:max(meta_adult$sum_class))){
  barplot(table(meta_adult$final_class[which(meta_adult$sum_class==i)]), main=paste("classifications=",i, sep=""))
}
# Also plot number of multiply classified cells per TSNE cluster 
par(las=3, mar=c(2,2,2,2), mfcol=c(3,4))
for(i in unique(meta_adult$Cluster_2d)){
  barplot(table(meta_adult$sum_class[which(meta_adult$Cluster_2d==i)]), 
          main=paste("TSNE cluster=", i,"/",unique(meta_adult$cell.type[which(meta_adult$Cluster_2d==i)]), sep=""), 
          xlab="Number of classifications")
}

# Finally create a heatmap with all different cell types and cell type specific genes 
all_genes <- unlist(list_markers_names)
class_list <- gsub('[0-9]+', '', names(all_genes))
class_list_col <- t(as.matrix(tableau_color_pal("tableau20")(20)[factor(class_list)]))
row.names(class_list_col) <- "gene.color"
# uniq_class <- unique(gsub('[0-9]+', '', names(all_genes)))
# b <- c()
# for(i in 1:length(uniq_class)){
#   a <- max(which(class_list==uniq_class[i]))
#   b <- c(b,a)
# }
# rowsep.vec <- b
types_ord <- names(list_markers_names)
data_temp <- data_norm_adult[all_genes,]
# Sort cells so they appear nice on the heatmap 
# Set levels 
meta_temp <- meta_adult
meta_temp[,"final_class_fac"] <- factor(meta_temp$final_class, levels=types_ord)
meta_temp <- meta_temp[order(meta_temp$final_class_fac),]
#meta_temp <-  meta_adult[with(meta_adult, order(sum_class, final_class)), ]
data_temp <- data_temp[,row.names(meta_temp)]
# Plot heatmap 
par(mfcol=c(1,1))
# Create bar color vector 
col_vec <- cbind(meta_temp$Cluster_color_2d, meta_temp$cell.type.color, meta_temp$final_class_color, meta_temp$sum_class.color)
colnames(col_vec) <- c("TSNE.cluster", "Cell.type.TSNE", "Cell.type.final", "Number.of.classifications")
heatmap.3(data_temp, trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="All markers all cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,2), ColSideColorsSize=7, Rowv=F, Colv=F, 
          colsep=c(max(which(meta_temp$sum_class==1)), max(which(meta_temp$sum_class==2)), max(which(meta_temp$sum_class==0))),
          sepcolor = "red", RowSideColors = class_list_col)
legend("bottomleft", 
       legend=unique(meta_temp$final_class), 
       fill=unique(meta_temp$final_class_color), cex=0.7)
legend("topright", 
       legend=unique(meta_temp$sum_class), 
       fill=unique(meta_temp$sum_class.color), cex=0.7)
# Also plot heatmap with 10 cells from each type along with double expressors 
meta_temp <- meta_temp[order(meta_temp$sum_class),]
#meta_temp <-  meta_adult[with(meta_adult, order(sum_class, final_class)), ]
data_temp <- data_temp[,row.names(meta_temp)]
# Plot heatmap 
par(mfcol=c(1,1))
# Create bar color vector 
col_vec <- cbind(meta_temp$Cluster_color_2d, meta_temp$cell.type.color, meta_temp$final_class_color, meta_temp$sum_class.color)
colnames(col_vec) <- c("TSNE.cluster", "Cell.type.TSNE", "Cell.type.final", "Number.of.classifications")
heatmap.3(data_temp, trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="All markers all cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,2), ColSideColorsSize=7, Rowv=F, Colv=F, 
          RowSideColors = class_list_col)
# Plot TSNE with overlayed clus sum 
plot(my.tsne.2d.adult, pch=19,  col=meta_adult$sum_class.color, main="Adult cells TSNE (number of classifications)")
legend("bottomleft", 
       legend=unique(meta_adult$sum_class), 
       fill=unique(meta_adult$sum_class.color), cex=0.7)
plot(my.tsne.2d.adult, pch=19, col=meta_adult$final_class_color, main="Adult cells TSNE (cell type final)")
plot(my.tsne.2d.adult, col="black", main="Adult cells TSNE (cell type final)", pch=as.numeric(as.factor(meta_adult$final_class)))
legend("topleft", 
       legend=unique(meta_adult$final_class), 
       pch=unique(as.numeric(as.factor(meta_adult$final_class))), cex=0.5, ncol=2)
# Also plot scatterplot of top markers for each cell type and contrast with canonical markers 
setwd("Gene_lists/Canonical/")
files <- list.files()
files <- files[grep(".txt", files)]
# Read all files of canonical markers 
list_canonical <- list()
for (f in 1:length(files)) {
  tempData = read.csv(files[f], header=F)
  tempData = apply(tempData, MARGIN = 2, function(x) gsub(",", "", x))
  list_canonical[[f]] <- tempData
  names(list_canonical)[[f]]  <- files[f]
} 
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis")
# Set names 
names(list_canonical) <- c("AT1", "AT2","Basal", "Club","Endothelial", "Fibroblast","Ciliated", "Neuroendocrine", "Neuron")
# Produce scatterplots 
for(i in 1:length(list_markers)){
plot(list_markers[[i]][,2], list_markers[[i]][,3], pch=19, cex=0.1, main=names(list_markers)[i], 
     xlab="Fraction of correct cells expressing", ylab="Fraction of other cells expressing")
abline(h=0.2)
abline(v=0.8)
text(list_markers[[i]][list_markers_names[[i]],2], list_markers[[i]][list_markers_names[[i]],3], 
        labels=list_markers_names[[i]], cex=0.6, col="red", adj = 1)

text(list_markers[[i]][list_markers_names[[i]],2], list_markers[[i]][list_markers_names[[i]],3], 
     labels=list_markers_names[[i]], cex=0.6, col="red", adj = 1)
a <- which(names(list_canonical)==names(list_markers)[i])
if(length(a) > 0){
vec <- list_canonical[[a]][,1]
# Find genes in data 
mat_gene_rows <- matrix(nrow=length(vec), ncol=2)
mat_gene_rows[,1] <- vec
for(j in 1:length(vec)){
  a <- which(row.names(list_markers[[i]])==vec[j])
  if(length(a) > 0) {mat_gene_rows[j,2] <- a}
}
mat_gene_rows <- mat_gene_rows[which(is.na(mat_gene_rows[,2])==F),]
text(list_markers[[i]][mat_gene_rows[,1],2], list_markers[[i]][mat_gene_rows[,1],3], 
     labels=vec, cex=0.7, col="blue", adj = 1)
}}

## Plot the double expressors 
# Create a matrix containing the normalized expression of each cell for each of the extended list og genes
uniq_cat <- unique(meta_adult$final_class[which(meta_adult$sum_class==1)])
mat_mean_norm <- matrix(nrow=nrow(meta_adult), ncol=length(uniq_cat))
row.names(mat_mean_norm) <- row.names(meta_adult)
colnames(mat_mean_norm) <- uniq_cat
# Populate matrix 
for(i in 1:nrow(mat_mean_norm)){
  cell <- row.names(mat_mean_norm)[i]
  for(j in 1:ncol(mat_mean_norm)){
    genes_temp <- as.vector(all_genes[grep(colnames(mat_mean_norm)[j],names(all_genes))])
    vec_temp <- data_norm_adult[genes_temp,cell]
      for(k in 1:length(vec_temp)){
        max_k <- max(data_norm_adult[names(vec_temp)[k],])
        vec_temp[k] <- vec_temp[k]/max_k}
    mat_mean_norm[i,j] <-  median(vec_temp)
  }
}
# Which double categories do we get ? 
multi_cat <- unique(meta_adult$final_class[which(meta_adult$sum_class==2)])
# Plot them (only selected categories for now)
# NE and club cells 
cells <- c(which(meta_adult$final_class=="NeuroendocrineClub"),
           which(meta_adult$final_class=="Club"),
           which(meta_adult$final_class=="Neuroendocrine"))
plot(mat_mean_norm[cells,"Club"],mat_mean_norm[cells,"Neuroendocrine"], pch=19,
     col=meta_adult$final_class_color[cells], xlab="Club signature", ylab="NE signature")
legend("bottomleft", legend = unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]))
# Also add heatmap 
genes_temp <- as.vector(all_genes[c(grep("Club", names(all_genes)),
                                    grep("Neuroendocrine", names(all_genes)))])
col_vec <- cbind(meta_adult$final_class_color[cells], meta_adult$sum_class.color[cells])
colnames(col_vec) <- c("Cell.type.final", "Number.of.classifications")

heatmap.3(data_norm_adult[genes_temp,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="NE-Club cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,4), ColSideColorsSize=7, Rowv=F, Colv=F, 
          colsep=c(max(which(meta_adult$sum_class[cells]==1)), max(which(meta_adult$sum_class[cells]==2))),
          sepcolor = "red")
legend("bottomleft", 
       legend=unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]), cex=0.7)

# NE and Endos
cells <- c(which(meta_adult$final_class=="NeuroendocrineEndothelial"),
           which(meta_adult$final_class=="Endothelial"),
           which(meta_adult$final_class=="Neuroendocrine"))
plot(mat_mean_norm[cells,"Endothelial"],mat_mean_norm[cells,"Neuroendocrine"], pch=19,
     col=meta_adult$final_class_color[cells], xlab="Endothelial signature", ylab="NE signature")
legend("bottomleft", legend = unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]))
# Also add heatmap 
genes_temp <- as.vector(all_genes[c(grep("Endothelial", names(all_genes)),
                                    grep("Neuroendocrine", names(all_genes)))])
col_vec <- cbind(meta_adult$final_class_color[cells], meta_adult$sum_class.color[cells])
colnames(col_vec) <- c("Cell.type.final", "Number.of.classifications")

heatmap.3(data_norm_adult[genes_temp,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="NE-Endothelial cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,4), ColSideColorsSize=7, Rowv=F, Colv=F, 
          colsep=c(max(which(meta_adult$sum_class[cells]==1)), max(which(meta_adult$sum_class[cells]==2))),
          sepcolor = "red")
legend("bottomleft", 
       legend=unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]), cex=0.7)

# NE and Glial
cells <- c(which(meta_adult$final_class=="GlialNeuroendocrine"),
           which(meta_adult$final_class=="Glial"),
           which(meta_adult$final_class=="Neuroendocrine"))
plot(mat_mean_norm[cells,"Glial"],mat_mean_norm[cells,"Neuroendocrine"], pch=19,
     col=meta_adult$final_class_color[cells], xlab="Glial signature", ylab="NE signature")
legend("bottomleft", legend = unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]))
# Also add heatmap 
genes_temp <- as.vector(all_genes[c(grep("Glial", names(all_genes)),
                                    grep("Neuroendocrine", names(all_genes)))])
col_vec <- cbind(meta_adult$final_class_color[cells], meta_adult$sum_class.color[cells])
colnames(col_vec) <- c("Cell.type.final", "Number.of.classifications")

heatmap.3(data_norm_adult[genes_temp,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="NE-Glial cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,4), ColSideColorsSize=7, Rowv=F, Colv=F, 
          colsep=c(max(which(meta_adult$sum_class[cells]==1)), max(which(meta_adult$sum_class[cells]==2))),
          sepcolor = "red")
legend("bottomleft", 
       legend=unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]), cex=0.7)

# NE and Glial
cells <- c(which(meta_adult$final_class=="FibroblastsNeuroendocrine"),
           which(meta_adult$final_class=="Fibroblasts"),
           which(meta_adult$final_class=="Neuroendocrine"))
plot(mat_mean_norm[cells,"Fibroblasts"],mat_mean_norm[cells,"Neuroendocrine"], pch=19,
     col=meta_adult$final_class_color[cells], xlab="Fibroblasts signature", ylab="NE signature")
legend("bottomleft", legend = unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]))
# Also add heatmap 
genes_temp <- as.vector(all_genes[c(grep("Fibroblasts", names(all_genes)),
                                    grep("Neuroendocrine", names(all_genes)))])
col_vec <- cbind(meta_adult$final_class_color[cells], meta_adult$sum_class.color[cells])
colnames(col_vec) <- c("Cell.type.final", "Number.of.classifications")

heatmap.3(data_norm_adult[genes_temp,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="NE-Fibroblasts cells", 
          labCol = "", cexRow = 0.9, ColSideColors = col_vec, margins=c(2,4), ColSideColorsSize=7, Rowv=F, Colv=F, 
          colsep=c(max(which(meta_adult$sum_class[cells]==1)), max(which(meta_adult$sum_class[cells]==2))),
          sepcolor = "red")
legend("bottomleft", 
       legend=unique(meta_adult$final_class[cells]), 
       fill=unique(meta_adult$final_class_color[cells]), cex=0.7)

dev.off()
###############################################################################################

# robust-sequential PCA ---------------------------------------------------
# First Subset metadata and data and convert into a data_table compatible with PCA functions 
cells_NE <- row.names(meta_adult)[c(which(meta_adult$final_class=="ClubNeuroendocrine"),
           which(meta_adult$final_class=="Club"),
           which(meta_adult$final_class=="Neuroendocrine"))]
#cells_NE <- row.names(meta_adult)[which(meta_adult$final_class=="Club" | meta_adult$final_class=="ClubNeuroendocrine")]
cells_NE <- row.names(meta_adult)[which(meta_adult$final_class=="Club")]
#common_out <- unique(c(mat_min_max[,1],mat_min_max[,2])[duplicated(c(mat_min_max[,1],mat_min_max[,2]))])
#cells_NE <- cells_NE[!cells_NE %in% common_out]
metadata_pca <- meta_adult[cells_NE,]
metadata_pca$plate.color <- tableau_color_pal("tableau20")(20)[as.factor(metadata_pca$plate)] 
data_pca <- data_adult[,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[,row.names(metadata_pca)]
# Create the data.table
require(data.table)
list_pca <- list()
for(i in 1:ncol(data_pca)){
  print(i)
  df <- as.data.frame(cbind(row.names(data_pca),rep(colnames(data_pca)[i],times=nrow(data_pca)),as.numeric(data_pca[,i]), as.numeric(data_norm_pca[,i])))
  df[,3] <- as.numeric(as.character(df[,3]))
  df[,4] <- as.numeric(as.character(df[,4]))
  list_pca[[i]] <- df
}
data_table_pca <- rbindlist(l = list_pca)
colnames(data_table_pca) <- c("gene", "cell_name", "counts", "log2_cpm")
data_table_pca[,gene:=as.character(gene)]
data_table_pca[,cell_name:=as.character(cell_name)]
# Add additional columns to data_table
data_table_pca[,"plate"] <- NA
data_table_pca[,plate:=as.character(plate)]
data_table_pca[,"plate_color"] <- NA
data_table_pca[,plate_color:=as.character(plate_color)]
data_table_pca[,"Age"] <- NA
data_table_pca[,Age:=as.character(Age)]
data_table_pca[,"Age_color"] <- NA
data_table_pca[,Age_color:=as.character(Age_color)]
data_table_pca[,"num.genes"] <- NA
data_table_pca[,num.genes:=as.numeric(num.genes)]
data_table_pca[,"cluster.tsne"] <- NA
data_table_pca[,cluster.tsne:=as.character(cluster.tsne)]
data_table_pca[,"cluster_color"] <- NA
data_table_pca[,cluster_color:=as.character(cluster_color)]
# And populate them
for(i in 1:nrow(metadata_pca)){
  print(i)
  data_table_pca[cell_name==row.names(metadata_pca)[i], plate := as.character(metadata_pca[i,"plate"])]
  #
  data_table_pca[cell_name==row.names(metadata_pca)[i], plate_color := metadata_pca[i,"plate.color"]]
  ## 
  data_table_pca[cell_name==row.names(metadata_pca)[i], Age := as.character(metadata_pca[i,"final_class"])]
  #
  data_table_pca[cell_name==row.names(metadata_pca)[i], Age_color := metadata_pca[i,"final_class_color"]]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], num.genes := metadata_pca[i,"Genes_detected"]]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], cluster.tsne := as.character(metadata_pca[i,"Cluster_2d"])]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], cluster_color := as.character(metadata_pca[i,"Cluster_color_2d"])]
}
# Start PCA analysis 
# top-level call
cordir="/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/Complete_dataset_analysis/PCA/Club_cells/"
pca.d2 <- dimplots.pca.wExptCols(data_table_pca,cordir,suffix='_Club_50genes',
                                 max.genes = 50,ncp=10,min.cells.detect=10,cexCol = 1)
plot.ggpairs.wExptCols(data_table_pca,dir=cordir,suffix='_Club_50genes',num.pcs = 10,height=22,width=22)
###############################################################################################

# Targeted gene list analysis NE cells ------------------------------------
# Add Index Sort information to metadata 
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/Gene_lists/Final_lists/")
# List of all files to read 
files <- list.files()
files <- files[grep(".txt", files)]
# Read all files
list_genes <- list()
for (f in 1:length(files)) {
  tempData = read.csv(files[f], header=F, sep="\t")
  tempData = apply(tempData, MARGIN = 2, function(x) gsub(",", "", x))
  list_genes[[f]] <- tempData
  names(list_genes)[[f]]  <- files[f]
} 
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/Complete_dataset_analysis/")
# Perform hierarchical clustering on NE cells using each of the gene lists 
cells_NE <- row.names(meta_adult)[which(meta_adult$final_class=="Neuroendocrine")]
#cells_NE <- row.names(meta_adult)[which(meta_adult$sum_class==1)]
#cells_NE <- row.names(meta_adult)[which(meta_adult$final_class=="Neuroendocrine")]
metadata_pca <- meta_adult[cells_NE,]
metadata_pca$plate.color <- tableau_color_pal("tableau20")(20)[as.factor(metadata_pca$plate)] 
data_pca <- data_adult[,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[,row.names(metadata_pca)]
# Find genes in the dataset 
for(i in 1:length(list_genes)){
  genes <- list_genes[[i]][,1]
  gene_f <- as.data.frame(genes)
  gene_f$"Row" <- NA
    for(j in 1:nrow(gene_f)){
      a <- which(row.names(data_norm_pca)==gene_f[j,1])
      if(length(a) > 0) {gene_f[j,"Row"] <- a} 
    }
  list_genes[[i]] <- gene_f
}
# Plot heatmaps 
pdf("ALL_cells_CAMs.pdf",15,15)
for(i in 1:length(list_genes)){
  genes <- unique(list_genes[[i]])
  genes <- as.character(genes[which(is.na(genes$Row)==F),1])
  col_labels <- as.matrix(metadata_pca$final_class_color)
  colnames(col_labels) <- c("cell_type")
  # Get genes and plot 
  # Select genes expressed in more than 5 cells 
  temp <- data_norm_pca[genes,]
  # Generate list of genes and number of expressing cells 
  exp_mat <- matrix(nrow=nrow(temp), ncol=length(unique(metadata_pca$final_class)))
  row.names(exp_mat) <- row.names(temp)
  colnames(exp_mat) <- unique(metadata_pca$final_class)
  for(l in unique(metadata_pca$final_class)){
    cells <- row.names(metadata_pca)[which(metadata_pca$final_class==l)]
    aaa <- which(unique(metadata_pca$final_class)==l)
    for(b in 1:nrow(temp)){
      exp_mat[b,aaa] <- length(which(temp[b,cells] != 0))/length(cells)
    }
  }
  write.table(exp_mat, file=paste(names(list_genes)[i], "cell_counts.txt", sep=""))
  # First remove non-expressed genes 
  genes <- names(which(rowSums(temp)!=0))
  # Plot all genes with no row clusterimg
#   # Heatmap 1 all cell and all expressed genes unclusterd
#   order <- order(metadata_pca$final_class)
#   heatmap.3(data_norm_pca[genes,order], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
#             hclustfun = function(x) hclust(x, 'ward.D2'), main=names(list_genes)[i], 
#             labCol = "", cexRow = 0.9, ColSideColors = col_labels[order,], margins=c(6,12), ColSideColorsSize=1, Rowv=F, Colv=F)    
  # Heatmap 1 all cell and all expressed genes unclusterd
  heatmap.3(data_norm_pca[genes,], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
            hclustfun = function(x) hclust(x, 'ward.D2'), main=names(list_genes)[i], 
            labCol = "", cexRow = 0.9, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=1, Rowv=F)    
  # Remove genes expressed in a few cells 
  temp <- data_norm_pca[genes,]
  count_temp <- apply(temp, 1, function(x) length(which(x!=0)))
  genes <- names(which(count_temp > 5))
  #Plot heatmap 
  # Heatmap 2 all cells and all expressed genes in more than 5 cells clustered
  heatmap.3(data_norm_pca[genes,], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
            hclustfun = function(x) hclust(x, 'ward.D2'), main=names(list_genes)[i], 
            labCol = "", cexRow = 0.9, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=1)    
  # Heatmap 3 all cells and all expressed genes in more than 5 cells unclustered
  heatmap.3(data_norm_pca[genes,], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
            hclustfun = function(x) hclust(x, 'ward.D2'), main=names(list_genes)[i], 
            labCol = "", cexRow = 0.9, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=1, Rowv=F)    
  # Remove cells with no expression of any gene 
  temp <- data_norm_pca[genes,]
  cells <- names(which(colSums(temp)!=0))
  meta_temp <- meta_adult[cells,]
  col_labels <- as.matrix(meta_temp$final_class_color)
  # And replot 
  # Heatmap 4 cells that express at least one gene clustered
  heatmap.3(data_norm_pca[genes,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
            hclustfun = function(x) hclust(x, 'ward.D2'), main=names(list_genes)[i], 
            labCol = "", cexRow = 0.9, margins=c(6,12), ColSideColors = col_labels)
  
}
dev.off()



###############################################################################################




# Various code ------------------------------------------------------------
pdf("plate_cell_assignment.pdf",10,10)
cells <- which(meta_adult$sum_class==1)
table_1 <- table(meta_adult$plate[cells], meta_adult$final_class[cells])
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table_1),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
dev.off()
pdf("Cre_plots.pdf",15,15)
# Plot zsGreen transcript data 
zs_row <- grep("zsGreen", row.names(data_norm_adult))
par(mfcol=c(2,2), mar=c(6,4,2,4), las=3)
box.with.strip(data = data_norm_adult[zs_row,], group = meta_adult$Cluster_2d,
               main=paste("zsGreen mRNA expression","TSNE cluster", sep="\n"), ylab="log2.CPM", xlab="")
box.with.strip(data = data_norm_adult[zs_row,], group = meta_adult$plate,
               main=paste("zsGreen mRNA expression","Plate", sep="\n"), ylab="log2.CPM", xlab="")
box.with.strip(data = data_norm_adult[zs_row,], group = meta_adult$zsGreen_class,
                main=paste("zsGreen mRNA expression","zsGreen class (high/low)", sep="\n"), ylab="log2.CPM", xlab="")
box.with.strip(data = data_norm_adult[zs_row,], group = meta_adult$final_class,
               main=paste("zsGreen mRNA expression","cell.type.final", sep="\n"), ylab="log2.CPM", xlab="")               
dev.off()
# Plot pulse width for plate 604. Doublets 
pdf("Pulse_plots_plate_604_doublets.pdf",15,10)
temp <- read.table("../PN FACS data-selected/Index_files/1000500604_pulse.txt", header=T)
par(mfcol=c(1,2), mar=c(4,4,3,3))
box.with.strip(temp$Pulse.width[which(temp$Cells==1 | temp$Cells==2)], group = temp$Cells[which(temp$Cells==1 | temp$Cells==2)], 
               xlab="Number of classifications", ylab="Pulse width", main="plate 1000500604")

wilcox.test(temp$Pulse.width[which(temp$Cells==1 | temp$Cells==2)]~temp$Cells[which(temp$Cells==1 | temp$Cells==2)])
# ROC analysis
temp_sub <- temp[which(temp$Cells==1 | temp$Cells==2),]
temp_sub$Cells <- as.factor(temp_sub$Cells)
pred <- prediction(temp_sub$Pulse.width, temp_sub$Cells)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))
perf <- performance(pred, measure = "auc")
text(x = 0.4, y=0.5, "AUC=")
text(x = 0.4, y=0.4, slot(perf, "y.values")[[1]])
# Do the same for number of genes 
meta_temp <- meta_adult[which(meta_adult$plate=="1000500604"),]
box.with.strip(meta_temp$Genes_detected, group = meta_temp$sum_class, 
               xlab="Number of classifications", ylab="Number of genes", main="plate 1000500604")
# ROC analysis
meta_temp <- meta_temp[which(meta_temp$sum_class!=3),]
pred <- prediction(meta_temp$Genes_detected, meta_temp$sum_class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col=rainbow(10))
perf <- performance(pred, measure = "auc")
text(x = 0.4, y=0.5, "AUC=")
text(x = 0.4, y=0.4, slot(perf, "y.values")[[1]])
# Combined metric of pulse width and number of genes 
temp <- read.table("../PN FACS data-selected/Index_files/1000500604_pulse.txt", header=T)
temp$Well <- as.character(temp$Well)
meta_temp <- meta_adult[which(meta_adult$plate=="1000500604"),]
# 
mat_temp <- as.data.frame(matrix(nrow=nrow(meta_temp), ncol=5))
colnames(mat_temp) <- c("well", "pulse.width", "genes", "classifications", "class_color")
row.names(mat_temp) <- row.names(meta_temp)
mat_temp[,1] <- meta_temp$well
mat_temp[,3] <- meta_temp$Genes_detected
mat_temp[,4] <- meta_temp$sum_class
mat_temp[,5] <- meta_temp$sum_class.color
for(i in 1:nrow(mat_temp)){
  a <- which(temp$Well==mat_temp[i,1])
  mat_temp[i,2]  <- temp$Pulse.width[a]
}
par(mfcol=c(1,1), mar=c(4,4,3,3))
plot(mat_temp[,2], mat_temp[,3], xlab="Pulse width", ylab="Number of detected genes", pch=19, col=mat_temp[,5])
legend("topleft", legend=unique(mat_temp[,4]), fill=unique(mat_temp[,5]) )
dev.off()


# Plot heatmap DESEq genes 
# Heatmap 2 all cells and all expressed genes in more than 5 cells clustered
deseq_res <- read.csv("DESEq_Club_Upk3a_all_genes.csv", header=T, row.names=1)
cells <- row.names(meta_adult)[which(meta_adult$final_class=="Club" | meta_adult$final_class=="ClubNeuroendocrine")]
#cells <- row.names(meta_adult)[grep("Club", meta_adult$final_class)]
col_labels <- as.matrix(meta_adult[cells,"final_class_color"])
pdf("DESEq_Club_Upk3a(+)_vs_Upk3a(-).pdf",10,10)
for(i in c(0.5,0.1,0.05,0.01)){
genes <- row.names(deseq_res)[which(deseq_res$padj < i)]
heatmap.3(data_norm_adult[genes,cells], trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main=paste("DESEq Club Upk3a(+) vs Upk3a(-)", "p.value=",i), 
          labCol = "", cexRow = 0.7, ColSideColors = col_labels, margins=c(6,12), ColSideColorsSize=1)    
}
dev.off()

###############################################################################################





# Ligand-receptor analysis ------------------------------------------------
# ligands <- read.table("../Gene_lists/neuropeptides_ligands_12-5.txt", sep="\t")
ligands <- as.data.frame(read.table("Classic_and_granin.txt", sep=""))
ligands <- as.data.frame(as.character(unique(ligands$V1)))
# colnames(ligands) <- c("gene", "ligand_number", "class", "description")
colnames(ligands) <- c("gene")
row.names(ligands) <- ligands$gene
receptors <- read.table("../Gene_lists/neuropeptides_receptors_12-5.txt", sep="\t")
colnames(receptors) <- c("gene", "ligand_number", "class", "description")
# Subset only uniqeuly classifed cells 
#cells <- row.names(meta_adult)[which(meta_adult$sum_class==1)]
cells <- row.names(meta_adult)
meta_cells <- meta_adult[cells,]
data_cells <- data_norm_adult[,cells]
# First keep only the genes expressed in the dataset 
genes <- as.character(ligands$gene)
b <- c()
for(i in 1:length(genes)){
  a <- which(row.names(data_cells)==genes[i])
  b <- c(b,a)
}
data_ligand <- data_cells[b,]
# For each list (ligands and receptors count the number of cells expressing each gene)
# temp <- matrix(nrow=nrow(ligands), ncol=length(unique(meta_cells$final_class)))
temp <- matrix(nrow=nrow(ligands), ncol=length(unique(meta_cells$cell.type)))
row.names(temp) <- row.names(ligands)
# colnames(temp) <- paste("Fraction_", unique(meta_cells$final_class), sep="")
colnames(temp) <- paste("Fraction_", unique(meta_cells$cell.type), sep="")
# for(i in 1:length(unique(meta_cells$final_class))){
#   cells <- row.names(meta_cells)[which(meta_cells$final_class==unique(meta_cells$final_class)[i])]  
#     for(j in 1:nrow(data_ligand)){
#     l <- length(which(data_ligand[j,cells] !=0))
#     r <- which(row.names(temp)==row.names(data_ligand)[j])
#     temp[r,i] <- l/length(cells)
#   }
# }
for(i in 1:length(unique(meta_cells$cell.type))){
  cells <- row.names(meta_cells)[which(meta_cells$cell.type==unique(meta_cells$cell.type)[i])]  
  for(j in 1:nrow(data_ligand)){
    l <- length(which(data_ligand[j,cells] !=0))
    r <- which(row.names(temp)==row.names(data_ligand)[j])
    temp[r,i] <- l/length(cells)
  }
}
# ligands <- cbind(ligands,temp)
ligands <- temp
ligands <- ligands[names(which(is.na(rowSums(ligands))==F)),]
write.table(ligands, file="Ligands_with_fraction_all_cell_types.txt")
# Remove non-found genes
pdf("NEuropeptide_ligands_fraction_comparison.pdf",10,10)
# ligands_less <- ligands[which(is.na(ligands$Fraction_Glial)==F),]
ligands_less <- ligands
layout(matrix(nrow=2, ncol=1, c(1,2)), heights = c(0.6,0.4))
par(mar=c(0.5,5,2,2))
boxplot(ligands_less[,grep("Fraction", colnames(ligands_less))], 
        outline = F, ylim=c(0,1), names = F, 
        ylab="Fraction of expressing cells")
stripchart(ligands_less[,grep("Fraction", colnames(ligands_less))],
             vertical = TRUE, method = "jitter", 
             pch = 21, col = "maroon", bg = "bisque", 
             add = TRUE)
par(mar=c(7,5,0.5,2))
names <- gsub("_", "\n",colnames(ligands_less)[grep("Fraction", colnames(ligands_less))])
barplot(apply(ligands_less[,grep("Fraction", colnames(ligands_less))], 2, function(x) length(which(x!=0)))/nrow(ligands_less), 
        ylim=c(0,1), names.arg = names, 
        ylab="Fraction of \n expressed genes \n (total NPs=69)")
for(i in c(1,2,4,5,6,7,8,9)){
  test <- wilcox.test(ligands_less[,grep("Fraction", colnames(ligands_less))][,2], 
              ligands_less[,grep("Fraction", colnames(ligands_less))][,i])
  print(colnames(ligands_less[,grep("Fraction", colnames(ligands_less))])[i])
  print(test$p.value)
}
dev.off()
# 
# Now do the same for the receptor list
# First keep only the genes expressed in the dataset 
temp <- matrix(nrow=nrow(receptors), ncol=length(unique(meta_cells$final_class)))
colnames(temp) <- paste("Fraction_", unique(meta_cells$final_class), sep="")
genes <- as.character(receptors$gene)
for(i in 1:length(unique(meta_cells$final_class))){
  cells <- row.names(meta_cells)[which(meta_cells$final_class==unique(meta_cells$final_class)[i])]  
  for(j in 1:length(genes)){
    a <- which(row.names(data_cells)==genes[j])
    if(length(a) !=0){l <- length(which(data_cells[a,cells] !=0)) ; temp[j,i] <- l/length(cells)}
      }
}
receptors <- cbind(receptors, temp)
write.table(receptors, file="Receptors_with_fraction_all_cell_types.txt")
# Plot ligands and corresponding receptors
pdf("NP_ligands_receptors.pdf", 10,10)
par(las=3, mar=c(5,5,5,2), cex.axis=0.7, mfcol=c(3,3))
for(i in 1:nrow(receptors)){
  lig <- unlist(strsplit(as.character(receptors[i,"ligand_number"]), split = ",", fixed = T))
  b <- c()
  for(j in 1:length(lig)){
    a <- paste(row.names(ligands[which(ligands$ligand_number==lig[j]),]), round(ligands[which(ligands$ligand_number==lig[j]),"Fraction_Neuroendocrine"],2), sep="=")
    b <- c(b,a)
  }
  barplot(as.matrix(receptors[i,c(5:ncol(receptors))]), main=b,ylab=as.character(receptors[i,"gene"]), beside = F, ylim=c(0,1))
}
dev.off()
# Make composite heatmap
# Only expressed genes 
genes <- names(which(rowSums(data_ligand)!=0))
# 
ligands <- ligands[genes, ]
data_ligand <- data_ligand[genes,]
meta_cells <- meta_cells[colnames(data_ligand),]
meta_cells$final_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_cells$final_class)]
# Make heatmap 
pdf("extended_heatmap_elements_Neuropeptide_ligands.pdf",15,15)
# Sort cells based cell type 
# meta_cells <- meta_cells[order(meta_cells$final_class),]
meta_cells <- meta_cells[order(meta_cells$cell.type),]
data_ligand <- data_ligand[,row.names(meta_cells)]
my.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
# test <- heatmap.2(data_ligand, trace="n", Colv=F, ColSideColors = meta_cells$final_class_color, col=my.palette, 
#           hclustfun = function(x) hclust(x, "ward.D2"))
test <- heatmap.2(data_ligand, trace="n", Colv=F, ColSideColors = meta_cells$cell.type.color, col=my.palette, 
                  hclustfun = function(x) hclust(x, "ward.D2"))
# genes in the same order as they appear in heatmap 
genes_on_heat <- row.names(data_ligand[rev(test$rowInd), test$colInd])
# make barplots of fraction of expressing cells 
# ligands_less <- ligands[rev(genes_on_heat),c(5:13)]
ligands_less <- ligands[rev(genes_on_heat),]
par(mfcol=c(1,ncol(ligands_less)), las=1, mar=c(2,1,1,1))

# for(i in 1:ncol(ligands_less)){
#   name <- strsplit(colnames(ligands_less)[i], split = "_")[[1]][2]
#   barplot2(as.matrix(ligands_less[,i]), beside = T, horiz = T, xlim=c(0,1),main=colnames(ligands_less)[i], 
#           names.arg = row.names(ligands_less), col=unique(meta_cells$final_class_color)[which(unique(meta_cells$final_class)==name)])
# }
for(i in 1:ncol(ligands_less)){
  name <- strsplit(colnames(ligands_less)[i], split = "_")[[1]][2]
  barplot2(as.matrix(ligands_less[,i]), beside = T, horiz = T, xlim=c(0,1),main=colnames(ligands_less)[i], 
           names.arg = row.names(ligands_less), col=unique(meta_cells$cell.type.color)[which(unique(meta_cells$cell.type)==name)])
}
dev.off()
# Make paired plots receptors and ligands with a "short" list format 
pdf("Mathced_ligands_receptors_short.pdf",15,15)
# matched_list <- read.csv("../Neuropeptides_matched_ligands_receptors_short.csv", header=T)
matched_list <- read.csv("170314_Ligands_receptors.csv", header=T)
# REmove genes not found in the dataset 
matched_list[,"Row_1"] <- NA
matched_list[,"Row_2"] <- NA
for(i in 1:nrow(matched_list)){
  a <- which(row.names(data_norm_adult)==as.character(matched_list[i,1]))
  if(length(a) > 0) {matched_list[i,"Row_1"] <- a }
  b <- which(row.names(data_norm_adult)==as.character(matched_list[i,2]))
  if(length(b) > 0) {matched_list[i,"Row_2"] <- b }
}
# Remove genes with no expression 
matched_list <- matched_list[which(is.na(rowSums(matched_list[,c(3,4)]))==F),]
# Subset and proceed with analysis 
# meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
meta_temp <- meta_adult
#meta_temp$final_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_temp$final_class)]
data_temp <- data_norm_adult[,row.names(meta_temp)]
# Make two tables
# One for expression of the ligands in NE cells 
mat_ligand <- as.data.frame(matrix(nrow=nrow(matched_list), ncol=2))
mat_ligand[,1] <- as.character(matched_list[,1])
colnames(mat_ligand) <- c("Ligand", "Neuroendocrine")
# cells <- row.names(meta_temp)[which(meta_temp$final_class=="Neuroendocrine")]
cells <- row.names(meta_temp)[which(meta_temp$cell.type=="NE")]
ne_col <- unique(meta_temp[cells,"final_class_color"])
for(i in 1:nrow(mat_ligand)){
  vec <- data_temp[mat_ligand[i,1],cells]
  mat_ligand[i,2] <- round(length(which(vec!=0))/length(vec), 3)
}
# One for expression of the receptors in non_NE cells 
#meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
meta_temp <- meta_adult
#meta_temp$final_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_temp$final_class)]
# cells <- row.names(meta_temp)[which(meta_temp$final_class!="Neuroendocrine")]
data_temp <- data_norm_adult[,row.names(meta_temp)]
classes <- unique(meta_temp$cell.type)
mat_receptor <- as.data.frame(matrix(nrow=nrow(matched_list), ncol=length(classes)+1))
mat_receptor[,1] <- as.character(matched_list[,2])
colnames(mat_receptor) <- c("Receptor", classes)
# for(i in 1:nrow(mat_receptor)){
#   for(j in 1:length(classes)){
#     cells <- row.names(meta_temp)[which(meta_temp$final_class==classes[j])]
#     vec <- data_temp[mat_receptor[i,1],cells]
#     mat_receptor[i,j+1] <- round(length(which(vec!=0))/length(vec), 3)
#   }
# }
for(i in 1:nrow(mat_receptor)){
  for(j in 1:length(classes)){
    cells <- row.names(meta_temp)[which(meta_temp$cell.type==classes[j])]
    vec <- data_temp[mat_receptor[i,1],cells]
    mat_receptor[i,j+1] <- round(length(which(vec!=0))/length(vec), 3)
  }
}
# And plot them 
mat <- matrix(nrow=1, ncol=2, c(1:2))
layout(mat, widths = c(0.3,0.7))
par(las=1)
barplot(mat_ligand[,2], horiz = T, names.arg = mat_ligand[,1], col=ne_col, xlim=c(0,1))
par(las=1)
# barplot(t(as.matrix(mat_receptor[,c(2:ncol(mat_receptor))])), horiz = T, names.arg = mat_receptor[,1], 
#         col=unique(meta_temp$final_class_color), beside = T, xlim=c(0,1))
# legend("topright", legend = unique(meta_temp$final_class), fill=unique(meta_temp$final_class_color), 
#        cex=0.8)
barplot(t(as.matrix(mat_receptor[,c(2:ncol(mat_receptor))])), horiz = T, names.arg = mat_receptor[,1], 
        col=unique(meta_temp$cell.type.color), beside = T, xlim=c(0,1))
legend("topright", legend = unique(meta_temp$cell.type), fill=unique(meta_temp$cell.type.color), 
       cex=0.8)
# Also individual plots per cell type for the repceptors 
# for(i in c(2:ncol(mat_receptor))){
#   mat <- matrix(nrow=1, ncol=2, c(1:2))
#   layout(mat, widths = c(0.3,0.7))
#   par(las=1)
#   barplot(mat_ligand[,2], horiz = T, names.arg = mat_ligand[,1], col=ne_col, xlim=c(0,1))
#   par(las=1)
#   barplot(mat_receptor[,i], horiz = T, names.arg = mat_receptor[,1], 
#           col=unique(meta_temp$final_class_color)[i-1], beside = T, xlim=c(0,1))
#   legend("topright", legend = unique(meta_temp$final_class)[i-1], fill=unique(meta_temp$final_class_color)[i-1], 
#          cex=0.8)
# }
for(i in c(2:ncol(mat_receptor))){
  mat <- matrix(nrow=1, ncol=2, c(1:2))
  layout(mat, widths = c(0.3,0.7))
  par(las=1)
  barplot(mat_ligand[,2], horiz = T, names.arg = mat_ligand[,1], col=ne_col, xlim=c(0,1))
  par(las=1)
  barplot(mat_receptor[,i], horiz = T, names.arg = mat_receptor[,1], 
          col=unique(meta_temp$cell.type.color)[i-1], beside = T, xlim=c(0,1))
  legend("topright", legend = unique(meta_temp$cell.type)[i-1], fill=unique(meta_temp$cell.type.color)[i-1], 
         cex=0.8)
}
dev.off()
# Ligand and receptor complementarity
# Find which ligands are expressed in NE cells 
# Which cells express the corresponding receptors ? 
ligands <- read.table("../Gene_lists/Final_lists/neuropeptides_ligands.txt", header=F, sep="\t")
ligands[,"Row"] <- NA
for(i in 1:nrow(ligands)){
  a <- which(row.names(data_norm_adult)==as.character(ligands$V1)[i]) 
  if(length(a) > 0) {ligands$Row[i] <- a}
}
ligands <- ligands[which(is.na(ligands$Row)==F),]
# 
receptors <- read.table("../Gene_lists/Final_lists/neuropeptides_receptors.txt", header=F, sep="\t")
receptors[,"Row"] <- NA
for(i in 1:nrow(receptors)){
  a <- which(row.names(data_norm_adult)==as.character(receptors$V1)[i]) 
  if(length(a) > 0) {receptors$Row[i] <- a}
}
receptors <- receptors[which(is.na(receptors$Row)==F),]
# Fraction of NE cells expressing each of the ligands 
cells_NE <- row.names(meta_adult)[ which(meta_adult$final_class=="Neuroendocrine")]
data_NE <- data_norm_adult[as.character(ligands$V1),cells_NE]
ligands[,"Fraction_NE"] <- NA
for(i in 1:nrow(ligands)){
  ligands[i,"Fraction_NE"] <- length(which(data_NE[as.character(ligands$V1)[i],] > 0))/ncol(data_NE)
}
# Sort ligand list 
ligands <- ligands[order(ligands$Fraction_NE, decreasing=T),]
# Remove ligands that are not expressed by NE cells 
ligands_ref <- ligands[which(ligands$Fraction_NE!=0),c("V1", "V2", "Fraction_NE")]
# Now find the corresponding receptors 
wmeta_temp <- meta_temp[which(meta_temp$final_class!="Neuroendocrine"),]
receptors[,unique(meta_temp$final_class)] <- 0
data_temp <- data_norm_adult[,row.names(meta_temp)]
for(i in 1:nrow(ligands_ref)){
  a <- which(as.character(receptors$V2)==ligands_ref$V2[i]) 
  if(length(a) > 0) {recs <- as.character(receptors$V1[a])}
  if(length(recs) > 0) {
    for(j in recs){
      for(k in unique(meta_temp$final_class)){
        cells <- row.names(meta_temp)[which(meta_temp$final_class==k)]
        receptors[which(receptors$V1==j),k] <- length(which(data_norm_adult[j,cells] > 0)) / length(cells)
      }
    }
  }
}
write.table(ligands_ref, "Annotated_ligands_NE.csv")
write.table(receptors, "Matching_receptors_NEandnonNE.csv")
############################################################


# Make plots for receptor pairs 
member_1 <- c("Calcrl","Calcr","Calcr","Calcr","Calcrl","Calcrl")
member_2 <- c("Ramp1", "Ramp1", "Ramp2", "Ramp3", "Ramp2", "Ramp3")
mat_coexp <- matrix(nrow=length(member_2), ncol=length(unique(meta_cells$final_class)))
colnames(mat_coexp) <- unique(meta_cells$final_class)
row.names(mat_coexp) <- paste(member_1, member_2, sep="-")

mat_exp_1 <- matrix(nrow=length(member_1), ncol=length(unique(meta_cells$final_class)))
row.names(mat_exp_1) <- member_1
colnames(mat_exp_1) <- unique(meta_cells$final_class)

mat_exp_2 <- matrix(nrow=length(member_2), ncol=length(unique(meta_cells$final_class)))
row.names(mat_exp_2) <- member_2
colnames(mat_exp_2) <- unique(meta_cells$final_class)

# Plot and populate coexpression table 
pdf("Receptor_pair_coexpression.pdf",15,15)
for(i in 1:nrow(mat_coexp)){
par(mfcol=c(3,3), mar=c(2,2,2,2))
  for(j in 1:length(unique(meta_cells$final_class))){
  cells <- row.names(meta_cells)[which(meta_cells$final_class==unique(meta_cells$final_class)[j])]  
  set_1 <- data_cells[member_1[i],cells]  
  set_2 <- data_cells[member_2[i],cells]
  plot(set_1, set_2, pch=19, col=unique(meta_cells$final_class_color)[j], 
       xlab="", ylab="", main=paste(unique(meta_cells$final_class)[j], row.names(mat_coexp)[i], sep="\n"), 
       xlim=c(0,max(set_1)+1), ylim=c(0,max(set_2)+1))
  # Expression 
  mat_exp_1[i,j] <-  length(which(set_1!=0))
  mat_exp_2[i,j] <-  length(which(set_2!=0))
  # Coexpression 
  set_1[which(set_1!=0)] <- 1
  set_2[which(set_2!=0)] <- 1
  set_3 <- set_1+set_2
  mat_coexp[i,j] <- length(which(set_3==2))
}
  }
# And plot tables 
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(mat_exp_1),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(mat_exp_2),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(mat_coexp),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
dev.off()


pdf("Plots_for_figures.pdf",15,15)
# TSNE with gene expression colors 
genes <- c("Calca", "Ascl1", "Syp", "Chga", "Resp18", "Pcsk1", "Snap25", "Scg5")
colfunc <- colorRampPalette(c("grey95","blue","yellow"))
my_palette <- colfunc(100)
ncol=100
for(i in 1:length(genes)){
exp <- data_norm_adult[genes[i], ]
col <- val2col(exp, col = pal(ncol))
layout(matrix(nrow=1, ncol=2, c(1:2)), widths = c(0.85,0.15))
plot(my.tsne.2d.adult, col=col,pch=19, xlab='', ylab='', cex=0.7, main=genes[i])
image.scale(exp, col=pal(ncol), axis.pos = 2)
}
#############
# TSNE colored by metadata 
# By plate
par(mfcol=c(1,1))
plot(my.tsne.2d.adult, col=meta_adult$plate.color,pch=19, xlab='', ylab='', cex=0.7, main="Plate")
legend("topleft", legend=unique(meta_adult$plate), fill=unique(meta_adult$plate.color))
# By injection 
meta_adult[,"Injection"] <- "Prenatal"
meta_adult[grep("80", meta_adult$plate),"Injection"] <- "Postnatal"
meta_adult[,"Injection_color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$Injection)]
plot(my.tsne.2d.adult, col=meta_adult$Injection_color,pch=19, xlab='', ylab='', cex=0.7, main="Injection")
legend("topleft", legend=unique(meta_adult$Injection), fill=unique(meta_adult$Injection_color))
#############
# zsGreen expression summary per plate 
mat_zs <- matrix(nrow=length(unique(meta_adult$plate)), ncol=6)
row.names(mat_zs) <- unique(meta_adult$plate)
colnames(mat_zs) <- c("NE\ncells", "Non_NE\ncells", "NE\nexpressing_ZsGreen","Non_NE\nexpressing_ZsGreen", 
                      "Mean_NE\nexpression","Mean_Non_NE\nexpression")
plates <- row.names(mat_zs)
for(i in 1:length(plates)){
  cells <- row.names(meta_adult)[which(meta_adult$plate==plates[i])]
  meta_temp <- meta_adult[cells,]
  meta_temp <- meta_temp[which(meta_temp$sum_class==1),]
  ne_cells <- row.names(meta_temp)[which(meta_temp$final_class=="Neuroendocrine")]
  non_ne_cells <- row.names(meta_temp)[which(meta_temp$final_class!="Neuroendocrine")]
  zs_ne <- data_norm_adult["zsGreen_transgene", ne_cells]
  zs_non_ne <- data_norm_adult["zsGreen_transgene", non_ne_cells]
  #boxplot(zs_ne, zs_non_ne, labels=c("NE", "Non-NE"))
  mat_zs[i,"NE\ncells"] <- length(ne_cells)
  mat_zs[i,"Non_NE\ncells"] <- length(non_ne_cells)
  mat_zs[i,"NE\nexpressing_ZsGreen"] <- length(which(zs_ne!=0))
  mat_zs[i,"Non_NE\nexpressing_ZsGreen"] <- length(which(zs_non_ne!=0))
  mat_zs[i,"Mean_NE\nexpression"] <- round(mean(zs_ne), 2)
  mat_zs[i,"Mean_Non_NE\nexpression"] <- round(mean(zs_non_ne), 2)
  }
mat_zs <- as.data.frame(mat_zs)
mat_zs[,"Fraction_NE\nexpressing"] <- round(mat_zs[,"NE\nexpressing_ZsGreen"]/mat_zs[,"NE\ncells"], 2)
mat_zs[,"Fraction_non\nNE_expressing"] <- round(mat_zs[,"Non_NE\nexpressing_ZsGreen"]/mat_zs[,"Non_NE\ncells"],2)
##
meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
table_1 <- table(meta_temp$plate, meta_temp$final_class)
# And plot tables 
mat_zs <- mat_zs[row.names(table_1),]
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(mat_zs),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
#
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table(meta_temp$plate, meta_temp$final_class)),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)

# Boxplots of zsgreen per plate and per cell type 
plates <- row.names(mat_zs)
par(mfcol=c(3,4), mar=c(5,3,3,3), las=3)
for(i in 1:length(plates)){
  cells <- row.names(meta_adult)[which(meta_adult$plate==plates[i])]
  meta_temp <- meta_adult[cells,]
  meta_temp <- meta_temp[which(meta_temp$sum_class==1),]
  boxplot(data_norm_adult["zsGreen_transgene", row.names(meta_temp)]~meta_temp$final_class, 
          main=plates[i],outpch = NA)
  stripchart(data_norm_adult["zsGreen_transgene", row.names(meta_temp)]~meta_temp$final_class, 
             vertical = TRUE, method = "jitter", 
             pch = 21, col = "maroon", bg = "bisque", 
             add = TRUE)

}
dev.off()



# FACS plots per plate 
########################
#setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/New_analysis/")
pdf("FACS_plots_per_plate.pdf",15,15)
plates <- unique(meta_adult$plate)
for(i in 1:length(plates)){
par(mfcol=c(2,2), mar=c(4,4,4,4))  
meta_tt <- meta_adult[which(meta_adult$plate==plates[i]),]
plot(meta_tt$FSC, meta_tt$SSC, xlim=c(0, max(meta_adult$FSC,na.rm=T)), ylim=c(0, max(meta_adult$SSC,na.rm=T)), 
     main=plates[i], xlab="FSC-A", ylab="SSC-A", pch=19, col=meta_tt$Cluster_color_2d)
abline(h=seq(from=0,to=3000000,by=50000), col="grey85")
abline(v=seq(from=0,to=3000000,by=50000), col="grey85")
plot(log10(meta_tt$DAPI), log10(meta_tt$APC_Lineage), xlim=c(0, log10(max(meta_adult$DAPI,na.rm=T))), ylim=c(0, log10(max(meta_adult$APC_Lineage,na.rm=T))), 
     main=plates[i], xlab="DAPI", ylab="APC", pch=19, col=meta_tt$Cluster_color_2d)
abline(h=seq(from=0,to=5,by=1), col="grey85")
abline(v=seq(from=0,to=5,by=1), col="grey85")
plot(log10(meta_tt$zsGreen), log10(meta_tt$`PE-Cy7_Epcam`), xlim=c(0, log10(max(meta_adult$zsGreen,na.rm=T))), ylim=c(0, log10(max(meta_adult$`PE-Cy7_Epcam`,na.rm=T))), 
     main=plates[i], xlab="zsGreen", ylab="PE-Cy7", pch=19, col=meta_tt$Cluster_color_2d)
abline(h=seq(from=0,to=5,by=1), col="grey85")
abline(v=seq(from=0,to=5,by=1), col="grey85")
plot(1,1,type="n", axes=F, xlab="", ylab="")
legend('topright', legend=unique(meta_adult$Cluster_2d), fill=unique(meta_adult$Cluster_color_2d),
       cex=0.8, ncol=2)
}
dev.off()



pdf("zsGreen_mRNA_boxplots.pdf",15,15)
par(las=3, mar=c(17,5,4,4))
boxplot(data_norm_adult["zsGreen_transgene",]~meta_adult$final_class, 
        main="zsGreen mRNA",outpch = NA)
stripchart(data_norm_adult["zsGreen_transgene",]~meta_adult$final_class, 
           vertical = TRUE, method = "jitter", 
           pch = 19, col = "maroon", bg = "bisque", 
           add = TRUE, cex=0.5)
dev.off()

save(list=c("meta_adult", "data_norm_adult", "data_adult"), file="Adult_Workspace_170113.RData")


cells <- which(meta_adult$sum_class==1)
table_1 <- table(meta_adult$plate[cells], meta_adult$final_class[cells])
table_2 <- table(meta_adult$plate[cells], meta_adult$cell.type[cells])
pdf("New_analysis/Initial_cell_assignment_by_plate.pdf",10,10)
qplot(1:10, 1:10, geom = "blank") + 
  theme_bw() +
  theme(line = element_blank(),
        text = element_blank(),
        panel.border = element_blank()) +
  annotation_custom(grob = tableGrob(table_2),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
dev.off()



# Correlation of FACS and mRNA values -------------------------------------
#meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
#data_temp <- data_norm_adult[,row.names(meta_temp)]
#meta_temp$final_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_temp$final_class)]
plates <- unique(meta_adult$plate)
pdf("zsGreen_and_Epcam_correlation_per_plate.pdf",15,15)
par(mfrow=c(2,2), mar=c(4,4,4,4))
# Set x and y lim 
meta_adult[,"zsGreen"]
meta_adult[,"zsGreen"]
for(i in plates){
  cells <- row.names(meta_adult)[which(meta_adult$plate==i)]
  plot(log10(meta_adult[cells,"zsGreen"]), data_norm_adult["zsGreen_transgene",cells], 
       col=meta_adult[cells,"Cluster_color_2d"], pch=19, xlab="FACS", ylab="mRNA",
       main=paste(i, "zsGreen", sep="\n"), 
       xlim=c(1,6), ylim=c(0,14))
  plot(log10(meta_adult[cells,"PE-Cy7_Epcam"]), data_norm_adult["Epcam",cells], 
       col=meta_adult[cells,"Cluster_color_2d"], pch=19, xlab="FACS", ylab="mRNA",
       main=paste(i, "EpCam", sep="\n"), 
       xlim=c(1,6), ylim=c(0,14))
  plot(log10(meta_adult[cells,"zsGreen"]), data_norm_adult["zsGreen_transgene",cells], 
       col=meta_adult[cells,"Cluster_color_2d"], pch=19, xlab="FACS", ylab="mRNA",
       main=paste(i, "zsGreen", sep="\n"), type="n", 
       xlim=c(1,6), ylim=c(0,14))
  text(log10(meta_adult[cells,"zsGreen"]), data_norm_adult["zsGreen_transgene",cells],
       meta_adult[cells,"well"], col=meta_adult[cells,"Cluster_color_2d"], 
       xlim=c(1,6), ylim=c(0,14))
  plot(log10(meta_adult[cells,"PE-Cy7_Epcam"]), data_norm_adult["Epcam",cells], 
       col=meta_adult[cells,"Cluster_color_2d"], pch=19, xlab="FACS", ylab="mRNA",
       main=paste(i, "EpCam", sep="\n"), type="n", 
       xlim=c(1,6), ylim=c(0,14))
  text(log10(meta_adult[cells,"PE-Cy7_Epcam"]), data_norm_adult["Epcam",cells],
       meta_adult[cells,"well"], col=meta_adult[cells,"Cluster_color_2d"], 
       xlim=c(1,6), ylim=c(0,14))
}
dev.off()  
################################################################################

# Make paired plots receptors and ligands with a "short" list format 
# For Neurotransmitters
# Since multiple genes participate in the making of a newurotransmitter
# the analysis will be slightly different to that for the neuropeptides 
setwd("New_analysis/")
pdf("Mathced_ligands_receptors_Neurotransmitters_Club.pdf",15,15)
neurot <- read.csv("Neurotransmitters.csv", header=T)
# REmove genes not found in the dataset 
neurot[,"Row"] <- NA
for(i in 1:nrow(neurot)){
  a <- which(row.names(data_norm_adult)==as.character(neurot[i,2]))
  if(length(a) > 0) {neurot[i,"Row"] <- a }
  }
# Remove genes with no expression 
neurot <- neurot[which(is.na(neurot$Row)==F),]
# Subset and proceed with analysis 
meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
meta_temp$final_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_temp$final_class)]
data_temp <- data_norm_adult[,row.names(meta_temp)]
# Make two tables
# One for expression of each of the genes in NE cells 
# Subset NE cells 
#data_temp_ne <- data_temp[,row.names(meta_temp)[which(meta_temp$final_class=="Neuroendocrine")]]
data_temp_ne <- data_temp[,row.names(meta_temp)[which(meta_temp$final_class=="Club")]]
mat_neurot_NE <- as.data.frame(matrix(nrow=nrow(neurot), ncol=ncol(data_temp_ne)))
colnames(mat_neurot_NE) <- colnames(data_temp_ne)
mat_neurot_NE <- cbind(mat_neurot_NE, neurot)
for(i in 1:nrow(mat_neurot_NE)){
  gene <- as.character(mat_neurot_NE[i,"Genes"])
  vec <- data_temp_ne[gene,]
  for(j in 1:length(vec)){
    if(vec[j]!=0){mat_neurot_NE[i,j] <- 1} else {mat_neurot_NE[i,j] <- 0}
  }
}
# Summarize for each gene and for each neurotransmitter 
mat_ligand <- as.data.frame(matrix(nrow=nrow(mat_neurot_NE), ncol=2))
mat_ligand[,1] <- as.character(mat_neurot_NE$Genes)
mat_ligand[,2] <- apply(mat_neurot_NE[,c(1:ncol(data_temp_ne))], MARGIN = 1, function(x) length(which(x==1))/ncol(data_temp_ne))
# For neurotransmitters
transmitters <- unique(as.character(mat_neurot_NE$Neutotransmitter))
mat_ligand_sum <- as.data.frame(matrix(nrow=length(transmitters), ncol=2))
mat_ligand_sum[,1] <- transmitters
for(i in 1:length(transmitters)){
  temp <- mat_neurot_NE[which(mat_neurot_NE$Neutotransmitter==transmitters[i]),c(1:ncol(data_temp_ne))]
  mat_ligand_sum[i,2] <- length(which(colSums(temp)==nrow(temp)))/ncol(data_temp_ne)
}
# Do the same for receptors on all cell types 
neurot_rec <- read.csv("Neurotransmitters_receptors.csv", header=T)
neurot_rec[,"Row"] <- NA
# REmove genes not found in the dataset 
for(i in 1:nrow(neurot_rec)){
  a <- which(row.names(data_norm_adult)==as.character(neurot_rec[i,2]))
  if(length(a) > 0) {neurot_rec[i,"Row"] <- a }
}
# Remove genes with no expression 
neurot_rec <- neurot_rec[which(is.na(neurot_rec$Row)==F),]
# Populate matrix with fraction of expressing cells for each cell type 
types <- unique(meta_temp$final_class)
mat_temp <- as.data.frame(matrix(nrow=nrow(neurot_rec), ncol=length(types)))
colnames(mat_temp) <- types 
for(i in 1:nrow(mat_temp)){
  gene <- as.character(neurot_rec$Receptor)[i]
  for(j in 1:length(types)){
    cells <- row.names(meta_temp)[which(meta_temp$final_class==types[j])]
    mat_temp[i,j] <-length(which(data_temp[gene,cells]!=0))/length(cells)
  }
}
mat_neurot_all <- cbind(neurot_rec, mat_temp)
# Now add the fraction of NE cells expressing each of the neurotransmitters from object : "mat_ligand_sum" 
mat_neurot_all[,"Neurot_fraction_NE(sum)"] <- NA
for(i in 1:nrow(mat_ligand_sum)){
  mat_neurot_all[which(as.character(mat_neurot_all$Neurotransmitter)==mat_ligand_sum[i,1]),"Neurot_fraction_NE(sum)"] <- mat_ligand_sum[i,2]
}
# And plot 
ne_col <- unique(meta_temp$final_class_color[which(meta_temp$final_class=="Neuroendocrine")])
mat <- matrix(nrow=1, ncol=2, c(1:2))
layout(mat, widths = c(0.3,0.7))
par(las=1)
barplot(mat_neurot_all$`Neurot_fraction_NE(sum)`, horiz = T, names.arg = mat_neurot_all$Neurotransmitter, 
        col=ne_col, xlim=c(0,1), cex.names = 0.5)
par(las=1)
barplot(t(as.matrix(mat_neurot_all[,types])), horiz = T, names.arg = mat_neurot_all$Receptor, 
        col=unique(meta_temp$final_class_color), xlim=c(0,1), cex.names = 0.5, beside = T)
legend("topright", legend = unique(meta_temp$final_class), fill=unique(meta_temp$final_class_color), 
       cex=0.8)
# Also one plot per cell type 
for(i in 1:length(types)){
mat <- matrix(nrow=1, ncol=2, c(1:2))
layout(mat, widths = c(0.3,0.7))
par(las=1)
barplot(mat_neurot_all$`Neurot_fraction_NE(sum)`, horiz = T, names.arg = mat_neurot_all$Neurotransmitter, 
        col=ne_col, xlim=c(0,1), cex.names = 0.5)
par(las=1)
barplot(t(as.matrix(mat_neurot_all[,types[i]])), horiz = T, names.arg = mat_neurot_all$Receptor, 
        col=unique(meta_temp$final_class_color)[i], xlim=c(0,1), cex.names = 0.5, beside = T)
legend("topright", legend = unique(meta_temp$final_class)[i], fill=unique(meta_temp$final_class_color)[i], 
       cex=0.8)
}
# Also add one for each gene of the neurotransmitters in NE cells 
mat_ligand[,"combo_name"] <- NA
for(i in 1:nrow(mat_ligand)){
  a <- which(mat_neurot_NE$Genes==mat_ligand[i,1])
  if(length(a) == 1) { mat_ligand[i,"combo_name"] <-   paste(as.character(mat_neurot_NE$Neutotransmitter)[a],mat_ligand[i,1], sep="---")}
  if(length(a) !=1 ) { mat_ligand[i,"combo_name"] <-   mat_ligand[i,1]}

}
par(mfcol=c(1,1), mar=c(3,10,2,2), las=1)
barplot(mat_ligand[,2], horiz = T, names.arg = mat_ligand[,3], xlim=c(0,1), col=ne_col)
dev.off()






# Classic and Granin ------------------------------------------------------
pdf("Number_NPs_expressed_by_each_cell_Granin_and_Classical_ordered.pdf",width = 30,height = 15)
# Classical 
# np <- c(as.character(read.table("Classic_NP_2.txt")[,1]))
#setwd("170305_reAnalysis/")
np <- c(as.character(read.table("Classic_and_granin_ordered.txt")[,1]))
np <- read.table("Classic_and_granin_ordered.txt")
b <- c()
for(i in 1:length(np)){
  a <- which(row.names(data_adult)==np[i])
  if(length(a) > 0 ) {b <- c(b,a)}
}
# 
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
metadata_pca <- meta_adult[cells,]
data_pca <- data_adult[b,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[b,row.names(metadata_pca)]
# Plot barplots 
par(las=3, mar=c(10,5,5,5))
barplot(apply(data_norm_pca, 2, function(x) length(which(x!=0))), cex.names = 0.5, main="Classic_NP_2_no_order")
#ord <- order(apply(data_norm_pca, 2, function(x) length(which(x!=0))))
ord <- order(apply(data_norm_pca[c(7:nrow(data_norm_pca)),], 2, function(x) length(which(x!=0))))
barplot(apply(data_norm_pca[,ord], 2, function(x) length(which(x!=0)))[ord], cex.names = 0.5, main="Classic_NP_2_ordered")
# Binary heatmap 
data_norm_pca[data_norm_pca!=0] <- 1
col=c("white", "black")
heatmap.2(data_norm_pca[,ord], trace="n", Colv=F, Rowv=F, col=col, colsep = c(0:ncol(data_norm_pca)), sepcolor="grey90",
          key=F, rowsep = c(0:nrow(data_norm_pca)), margins = c(10,10))
# Granin family
np <- c(as.character(read.table("Granin_NP.txt")[,1]))
b <- c()
for(i in 1:length(np)){
  a <- which(row.names(data_adult)==np[i])
  if(length(a) > 0 ) {b <- c(b,a)}
}
# 
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
metadata_pca <- meta_adult[cells,]
data_pca <- data_adult[b,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[b,row.names(metadata_pca)]
# Plot barplots 
barplot(sort(apply(data_norm_pca, 2, function(x) length(which(x!=0)))), cex.names = 0.5, main="Granin_NP_ordered")
barplot(apply(data_norm_pca, 2, function(x) length(which(x!=0)))[ord], cex.names = 0.5, main="Granin_NP_ordered_by_classical_NP")
# Binary heatmap 
data_norm_pca[data_norm_pca!=0] <- 1
col=c("white", "black")
heatmap.2(data_norm_pca[,ord], trace="n", Colv=F, Rowv=F, col=col, colsep = c(0:ncol(data_norm_pca)), sepcolor="grey90",
          key=F, rowsep = c(0:nrow(data_norm_pca)))
# Combined binary heatmap 
np <- c(as.character(read.table("Granin_NP.txt")[,1]), as.character(read.table("Classic_NP_2.txt")[,1]))
b <- c()
for(i in 1:length(np)){
  a <- which(row.names(data_adult)==np[i])
  if(length(a) > 0 ) {b <- c(b,a)}
}
# 
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
metadata_pca <- meta_adult[cells,]
data_pca <- data_adult[b,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[b,row.names(metadata_pca)]
data_norm_pca[data_norm_pca!=0] <- 1
col=c("white", "black")
colvec <- c(rep("red", 6), rep("blue", nrow(data_norm_pca)-6))
heatmap.2(data_norm_pca[,ord], trace="n", Colv=F, Rowv=F, col=col, colsep = c(0:ncol(data_norm_pca)), sepcolor="grey90",
          key=F, rowsep = c(0:nrow(data_norm_pca)), RowSideColors = colvec, cexRow = 2, margins = c(10,10))
legend("topleft", legend=c("Classic", "Granin"), fill=c("red", "blue"))
# Remove zero genes 
data_norm_pca <- data_norm_pca[which(rowSums(data_norm_pca)!=0),]
colvec <- c(rep("red", 6), rep("blue", nrow(data_norm_pca)-6))
heatmap.2(data_norm_pca[,ord], trace="n", Colv=F, Rowv=F, col=col, colsep = c(0:ncol(data_norm_pca)), sepcolor="grey90",
          key=F, rowsep = c(0:nrow(data_norm_pca)), RowSideColors = colvec, cexRow = 2, margins = c(10,10),
          main="Non_zero_genes_only")
# 
ord_2 <- order(rowSums(data_norm_pca), decreasing = T)
heatmap.2(data_norm_pca[ord_2,], trace="n", Colv=F, Rowv=F, col=col, colsep = c(0:ncol(data_norm_pca)), sepcolor="grey90",
          key=F, rowsep = c(0:nrow(data_norm_pca)), RowSideColors = colvec, cexRow = 2, margins = c(10,10),
          main="Non_zero_genes_only ordered by gene")
dev.off()
# Also plot as an actual heatmap
np <- read.table("Classic_and_Granin_NP.txt", header=T)
# Subset cells 
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
metadata_pca <- meta_adult[cells,]
data_pca <- data_adult[,row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[,row.names(metadata_pca)]
# Add row info for each gene 
np[,"Row"] <- NA
for(i in 1:nrow(np)){
  a <- which(row.names(data_norm_pca)==as.character(np[i,1]))
  if(length(a) > 0 ) {np[i,"Row"] <- a }
}
# Remove not found genes and Add number of expressing cells 
np <- np[which(is.na(np$Row)==F),]
data_pca <- data_adult[as.character(np$Gene),row.names(metadata_pca)]
data_norm_pca <- data_norm_adult[as.character(np$Gene),row.names(metadata_pca)]
np[,"expressing.cells"] <- apply(data_norm_pca, MARGIN = 1, function(x) length(which(x > 0)))
# Sort Classic and Granin 
np_1 <- np[which(np$Class=="Classic"),]
np_1 <- np_1[order(np_1$expressing.cells, decreasing = T),]
np_2 <- np[which(np$Class=="Granin"),]
np_2 <- np_2[order(np_2$expressing.cells, decreasing = T),]
np_ord <- rbind(np_1, np_2)
# Sort cells 
cell_ord <- apply(data_norm_pca, MARGIN = 2, function(x) length(which(x > 0)))
cells <- names(cell_ord[order(cell_ord, decreasing = T)])
data_norm_pca <- data_norm_adult[as.character(np_ord$Gene),cells]
# Plot heatmap
pdf("Classic_and_Granin_heatmap.pdf", width = 15, height = 10)
#colfunc <- colorRampPalette(c("white","yellow","red"))
colfunc <- colorRampPalette(c("white","black"))
my_palette <- colfunc(10)
ncol=10
heatmap.2(data_norm_pca, trace="n", Colv=F, Rowv=F, col=my_palette, colsep = c(0:ncol(data_norm_pca)), sepcolor="white",
          key=T, rowsep = c(0:nrow(data_norm_pca)), cexRow = 0.7, margins = c(5,5),
          main="All_genes_only ordered by gene")
heatmap.2(data_norm_pca[names(which(rowSums(data_norm_pca)!=0)),], trace="n", Colv=F, Rowv=F, col=my_palette, colsep = c(0:ncol(data_norm_pca)), sepcolor="white",
          key=T, rowsep = c(0:nrow(data_norm_pca)), cexRow = 1, margins = c(5,5),
          main="Non_zero_genes_only ordered by gene")
dev.off()
##################################################


# List of genes expressed by cell and cells by gene 
##################################################
list_genes <- list()
for(i in 1:nrow(data_norm_pca)){
  list_genes[[i]] <- names(which(data_norm_pca[i,] != 0 ))
}
names(list_genes) <- row.names(data_norm_pca)
# 
list_cells <- list()
for(i in 1:ncol(data_norm_pca)){
  list_cells[[i]] <- names(which(data_norm_pca[,i] != 0 ))
}
names(list_cells) <- colnames(data_norm_pca)

save(list = c("data_norm_adult","list_genes", "list_cells", "meta_adult"),file = "NPs_cells_expressing_each_gene.RData" )
##################################################





# Unused code -------------------------------------------------------------
# Analysis of neuropeptides 
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
# 
np <- c(as.character(read.table("Classic_NP_2.txt")[,1]))
b <- c()
for(i in 1:length(np)){
  a <- which(row.names(data_adult)==np[i])
  if(length(a) > 0 ) {b <- c(b,a)}
}
# PCA analysis
metadata_pca <- meta_adult[cells,]
data_pca <- data_adult[b,row.names(metadata_pca)]
data_pca <- data_pca[names(which(rowSums(data_pca)!=0)),]
data_norm_pca <- data_norm_adult[row.names(data_pca),row.names(metadata_pca)]
# Create the data.table
require(data.table)
list_pca <- list()
for(i in 1:ncol(data_pca)){
  print(i)
  df <- as.data.frame(cbind(row.names(data_pca),rep(colnames(data_pca)[i],times=nrow(data_pca)),as.numeric(data_pca[,i]), as.numeric(data_norm_pca[,i])))
  df[,3] <- as.numeric(as.character(df[,3]))
  df[,4] <- as.numeric(as.character(df[,4]))
  list_pca[[i]] <- df
}
data_table_pca <- rbindlist(l = list_pca)
colnames(data_table_pca) <- c("gene", "cell_name", "counts", "log2_cpm")
data_table_pca[,gene:=as.character(gene)]
data_table_pca[,cell_name:=as.character(cell_name)]
# Add additional columns to data_table
data_table_pca[,"plate"] <- NA
data_table_pca[,plate:=as.character(plate)]
data_table_pca[,"plate_color"] <- NA
data_table_pca[,plate_color:=as.character(plate_color)]
data_table_pca[,"Age"] <- NA
data_table_pca[,Age:=as.character(Age)]
data_table_pca[,"Age_color"] <- NA
data_table_pca[,Age_color:=as.character(Age_color)]
data_table_pca[,"num.genes"] <- NA
data_table_pca[,num.genes:=as.numeric(num.genes)]
data_table_pca[,"cluster.tsne"] <- NA
data_table_pca[,cluster.tsne:=as.character(cluster.tsne)]
data_table_pca[,"cluster_color"] <- NA
data_table_pca[,cluster_color:=as.character(cluster_color)]
# And populate them
for(i in 1:nrow(metadata_pca)){
  print(i)
  data_table_pca[cell_name==row.names(metadata_pca)[i], plate := as.character(metadata_pca[i,"plate"])]
  #
  data_table_pca[cell_name==row.names(metadata_pca)[i], plate_color := metadata_pca[i,"plate.color"]]
  ## 
  data_table_pca[cell_name==row.names(metadata_pca)[i], Age := as.character(metadata_pca[i,"final_class"])]
  #
  data_table_pca[cell_name==row.names(metadata_pca)[i], Age_color := metadata_pca[i,"final_class_color"]]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], num.genes := metadata_pca[i,"Genes_detected"]]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], cluster.tsne := as.character(metadata_pca[i,"Cluster_2d"])]
  # #
  data_table_pca[cell_name==row.names(metadata_pca)[i], cluster_color := as.character(metadata_pca[i,"Cluster_color_2d"])]
}
# Start PCA analysis 
# top-level call
source('~/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/New_analysis/Robust_PCA_functions.R')
cordir="/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/New_analysis/NP_analysis/PCA_NPgenes_NEcells/"
pca.d2 <- dimplots.pca.wExptCols(data_table_pca,cordir,suffix='_NP_10_genes',
                                 max.genes = 10,ncp=5,min.cells.detect=10,cexCol = 1)
#plot.ggpairs.wExptCols(data_table_pca,dir=cordir,suffix='_NP_10_genes',num.pcs = 3,height=22,width=22)
###############################################################################################
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/New_analysis/NP_analysis/PCA_NPgenes_NEcells/")
pca.cells <- read.csv("PCscores_NP_10_genes.csv", header=T, row.names = 5)
pca.cells <- pca.cells[,-1]
pca.genes <- read.csv("PCloadings_NP_10_genes.csv", header=T, row.names = 5)
# 
genes_to_plot <- row.names(data_norm_pca)
# Subset 
data_sloan <- data_norm_pca[genes_to_plot,]

pdf("pairwise_PCA_plots_with_gene_expression.pdf",10,10)
for(k in 1:ncol(pca.cells)){
  for(f in 1:ncol(pca.cells)){
    # Make sure that data, metadata and TSNE graph is ordered
    my.tsne.2d.sloan <- pca.cells[cells,c(k,f)]
    # 
    par(mfcol=c(1,1), mar=c(4,4,4,4))
    plot(my.tsne.2d.sloan)
    # Plotting gene expression on TSNE graph ----------------------------------
    #pdf("Immune_markers_in_silico_stains_short.pdf", width = 16  ,height = 10)
    pal=colorRampPalette(c("grey95","dodgerblue2","darkgoldenrod1"))
    ncol=100
    layout(matrix(nrow=3, ncol=6, c(1:18), byrow = T))
    # layout(matrix(nrow=2, ncol=3, c(1:18), byrow = T))
    par(mar=c(2,0,2,0))
    # Create the color matrix 
    col <- val2col(data_sloan, col = pal(ncol))
    col_mat <- matrix(col, nrow=nrow(data_sloan), ncol=ncol(data_sloan))
    colnames(col_mat) <- colnames(data_sloan) 
    min.o <- apply(my.tsne.2d.sloan, MARGIN = 2, min)
    max.o <- apply(my.tsne.2d.sloan, MARGIN = 2, max) 
    for(i in 1:nrow(data_sloan)){
      set_1 <- names(which(data_sloan[i,] == 0))
      set_2 <- names(which(data_sloan[i,] > 0))
      plot(my.tsne.2d.sloan, col="grey55",pch=19, xlab='', ylab='', cex=1.3, 
           main=row.names(data_sloan)[i], xaxt="n", yaxt="n", xlim=c(min.o[1],max.o[1]), 
           ylim=c(min.o[2],max.o[2]))
      par(new=T)
      plot(my.tsne.2d.sloan[set_1,], col=col_mat[i,set_1],pch=19, xlab='', ylab='', cex=0.7, 
           main=row.names(data_sloan)[i], xaxt="n", yaxt="n", xlim=c(min.o[1],max.o[1]), 
           ylim=c(min.o[2],max.o[2]))
      par(new=T)
      plot(my.tsne.2d.sloan[set_2,], col=col_mat[i,set_2],pch=19, xlab='', ylab='', cex=0.7, 
           main="", xaxt="n", yaxt="n", xlim=c(min.o[1],max.o[1]), 
           ylim=c(min.o[2],max.o[2]))
    }
    par(mar=c(2,4.5,2,4.5), las=1)
    image.scale(data_sloan, col = pal(ncol), axis.pos = 2)  
  }}
#
dev.off()



genes_to_plot <- c("Calca", "Dbi", "Cartpt")
# Subset 
order <- order(pca.cells$PC1)
data_sloan <- data_norm_pca[genes_to_plot,order]
# Colors 
data_sloan_norm <- data_sloan
for(i in 1:nrow(data_sloan_norm)){
  data_sloan_norm[1,] <- lowess(data_sloan[1,])$y
}
# Plot 
ncol=100
pal=colorRampPalette(c("grey95","dodgerblue2","darkgoldenrod1"))
for(i in 1:nrow(data_sloan_norm)){
  plot(x=c(1:ncol(data_sloan)), y = rep(i,ncol(data_sloan)), pch=19, cex=0.6, 
       col=val2col(data_sloan_norm[i,], col = pal(ncol)), ylim=c(0,nrow(data_sloan)+1))
  par(new=T)
}
##############################################################################
# Lung development receptors and ligands ----------------------------------
setwd("New_analysis/")
receptors <- read.csv("LungDev_Ligands.txt", header=F)
row.names(receptors) <- receptors[,1]
# find all genes in dataset 
receptors[,"row"] <- NA
for(i in 1:nrow(receptors)){
  a <- which(row.names(data_norm_adult)==row.names(receptors)[i])
  if(length(a) >0) {receptors$row[i] <- a }
}
# Remove not found genes 
receptors <- receptors[-which(is.na(receptors$row)==T),]
# Subset uniquely classifed cells 
meta_temp <- meta_adult[which(meta_adult$sum_class==1),]
# Create expression count matrix 
mat_temp <- matrix(nrow=nrow(receptors), ncol=length(unique(meta_temp$final_class)))
colnames(mat_temp) <- unique(meta_temp$final_class)
# 
for(i in 1:nrow(receptors)){
  data_temp <- data_norm_adult[receptors[i,"row"],]
  for(j in 1:length(unique(meta_temp$final_class))){
    mat_temp[i,j] <- length( which(data_temp[row.names(meta_temp)[which(meta_temp$final_class==unique(meta_temp$final_class)[j])]] != 0))/
      length(which(meta_temp$final_class==unique(meta_temp$final_class)[j]))
  }
}
# 
write.table(cbind(receptors, mat_temp), file="Lung_development_ligands_ALL_fraction_of_expressing_cells.csv")


# # Differential expression Upk3a club cells 
# # Subset cells 
# cells <- row.names(meta_adult)[which(meta_adult$final_class=="Club")]
# cells <- cells[-grep("50060", meta_adult[cells,"plate"])]
# # Separate in two groups based on Upk3a 
# club_zero <- names(which(data_norm_adult["Upk3a",cells]==0))
# club_non_zero <- names(which(data_norm_adult["Upk3a",cells]!=0))
# # Subset data 
# data_club <- data_norm_adult[,cells]
# # Create matrix to store data 
# mat_club <- matrix(nrow=nrow(data_club), ncol=3)
# # 
# row.names(mat_club) <- row.names(data_club)
# colnames(mat_club) <- c("fraction_positive", "fraction_negative", "p.value")
# for(i in 1:nrow(mat_club)){
#   d1 <- data_club[i,club_non_zero]
#   mat_club[i,1]  <- length(which(d1!=0))/length(club_non_zero)
#   d2 <- data_club[i,club_zero]
#   mat_club[i,2]  <- length(which(d2!=0))/length(club_zero)
#   test <- wilcox.test(d1,d2)
#   mat_club[i,3]  <- test$p.value
# }
# # Remove genes not expressed in any of the two populations 
# genes <- names(which(rowSums(mat_club[,c(1,2)])!=0))
# mat_club_small <- mat_club[genes,]
# # 
# write.table(mat_club_small, file="UPK3a_differential_expression.csv")
# 



# Barcode analysis --------------------------------------------------------
genes <- as.character(read.csv("Classic_and_granin.txt", header=F)[,1])
# Find genes 
b <- c()
for(i in 1:length(genes)){
  a <- which(row.names(data_norm_adult)==genes[i])
  b <- c(b,a)
}
cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
data_temp <- data_norm_adult[b, cells]
# Convert to binary
data_temp[data_temp!=0] <- 1
# Sort in order of # of expressing cells 
data_temp <- data_temp[order(rowSums(data_temp), decreasing = T),]
# Generate cell barcodes 
list_cells <- list()
names_t <- c()
temp <- as.data.frame(1:nrow(data_temp))
for(i in 1:ncol(data_temp)){
  list_cells[[i]] <-  paste(as.vector(data_temp[,i]), collapse = "")
  data_t <- as.data.frame(data_temp[,i])
  names_t <- c(names_t,list_cells[[i]])
  temp <- cbind(temp,data_t)
}
genes_barcode <- row.names(data_temp)
temp <- temp[,-1]
colnames(temp) <- names_t
names(list_cells) <- colnames(data_temp)
# Make a list for barcodes 
list_barcodes <- list()
barcodes  <- unlist(unique(list_cells))
for(i in 1:length(barcodes)){
  list_barcodes[[i]] <- names(list_cells)[grep(barcodes[i], list_cells)]
}
names(list_barcodes) <- barcodes
vec <- unlist(lapply(list_barcodes, length))
# Plotting 
pdf("Classic_and_Granin.pdf", width = 20, height = 15)
par(mar=c(5,15,2,2), las=2)
barplot(vec[order(vec, decreasing = T)], horiz = T, cex.names = 0.6, main="Classic_and_Granin", 
        xlab="Number of cells sharing barcode", xlim=c(0,60))
abline(v=seq(from=0,to=60,by=5), col="grey85", lty=3)
table <- do.call(cbind, list_barcodes)
table <- table[,names(vec[order(vec, decreasing = T)])]
colnames(table) <-  paste("X", colnames(table))
write.table(table, file="Barcodes_and_cell_members_G.txt")
write.table(vec, file="Barcodes_and_cell_members_2_G.txt")
write.csv(temp[,names(vec[order(vec, decreasing = T)])], "Barcode_legend_G.csv")
# Add number of barcodes each gene is found in 
mat_genes_bc <- matrix(nrow=length(genes_barcode), ncol=2)
colnames(mat_genes_bc) <- c("number.barcodes", "number.cells")
row.names(mat_genes_bc) <- genes_barcode
for(i in 1:nrow(mat_genes_bc)){
  a <- c()
  b <- c()
  mat_genes_bc[i,2] <- length(which(data_temp[row.names(mat_genes_bc)[i],] !=0))/ncol(data_temp)
  for(j in 1:length(barcodes)){
    a <- strsplit(barcodes[j], split="")[[1]][i]
    b <- c(a,b)
    mat_genes_bc[i,1] <-  sum(as.numeric(b))
    }
}
# Barcode with better legend
layout(matrix(nrow=1, ncol=2, c(1:2)), widths = c(1,2))
# Legend 
par(mar=c(5,15,2,0), las=2)
plot(x=0,y=0,type="n", ylim=c(0,length(barcodes)), xlim=c(0,length(genes_barcode)),
     axes=F, xlab="", ylab="")
for(i in 1:length(barcodes)){
  par(new=T)
  col_bc <- strsplit(barcodes[i], "")[[1]]
  col_bc[col_bc==0] <- "white"
  col_bc[col_bc==1] <- "red"
  plot(x = 1:length(genes_barcode), y= rep(i,length(genes_barcode)), 
       col=col_bc,ylim=c(0,length(barcodes)), xlim=c(0,length(genes_barcode)), 
       axes=F, xlab="", ylab="", pch=15, cex=0.7)
}
abline(v=1:length(genes_barcode), lty=3)
abline(h=1:length(barcodes))
par(cex.axis=0.5)
axis(side = 1, at = 1:length(genes_barcode), labels = genes_barcode)
par(cex.axis=0.5)
axis(side = 2, at = 1:length(barcodes), labels = barcodes)
par(cex.axis=1, las=3)
text(x = 1:length(genes_barcode), y=length(barcodes)+1, labels = mat_genes_bc[,1], col="red", cex=0.5)
# Barplot 
par(mar=c(5.3,0,1.9,2), las=2)
barplot(vec[order(vec, decreasing = T)], horiz = T, cex.names = 0.6, main="", 
        xlab="Number of cells sharing barcode", xlim=c(0,60), names.arg = "")
abline(v=seq(from=0,to=60,by=5), col="grey85", lty=3)
# Number of cells each gene is expressed in 
layout(matrix(nrow=1, ncol=1, 1))
par(mar=c(5,5,5,5), las=2)
barplot(mat_genes_bc[,2], main="Fraction of expressing cells", ylim=c(0,1))
dev.off()
######################



# Number of NPs per cell type plus human data -----------------------------
# Mouse dataset 
genes_m <- as.character(read.csv("All_NP_genes_mouse.txt", header=F)[,1])
# Find genes 
b <- c()
for(i in 1:length(genes_m)){
  a <- which(row.names(data_norm_adult)==genes_m[i])
  b <- c(b,a)
}
# cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
data_m <- data_norm_adult[b,]
mat_np <- matrix(nrow=length(unique(meta_adult$cell.type)), ncol=1)
row.names(mat_np) <- unique(meta_adult$cell.type)
for(i in 1:nrow(mat_np)){
  type <- row.names(mat_np)[i]
  cells <- row.names(meta_adult)[which(meta_adult$cell.type==type)]
  data_temp <- data_m[,cells]
  mat_np[i,1] <- length(which(rowSums(data_temp) !=0))
}
# Human  
human_data <- read.csv("~/Desktop/Datasets_and_metadata/PNAS_dataset/All_cell_counts_brain.csv", header=T, row.names=1, check.names = F)
human_meta <- read.csv("~/Desktop/Datasets_and_metadata/PNAS_dataset/All_cell_info_brain.csv", header=T, row.names=1, sep=";")
# Read human NP gene names 
genes_h <- as.character(read.csv("All_NP_genes_human.txt", header=F)[,1])
# Find genes 
b <- c()
for(i in 1:length(genes_h)){
  a <- which(row.names(human_data)==genes_h[i])
  b <- c(b,a)
}
# cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
data_h <- human_data[b,]
mat_np_h <- matrix(nrow=length(unique(human_meta$Cell_type)), ncol=1)
row.names(mat_np_h) <- unique(human_meta$Cell_type)
for(i in 1:nrow(mat_np_h)){
  type <- row.names(mat_np_h)[i]
  cells <- row.names(human_meta)[which(human_meta$Cell_type==type)]
  data_temp <- data_h[,cells]
  mat_np_h[i,1] <- length(which(rowSums(data_temp) !=0))
}
# Remove hybrids 
mat_np_h <- as.matrix(mat_np_h[-grep('hybrid', row.names(mat_np_h)),])
mat_np_h <- as.matrix(mat_np_h[-grep('Fetal', row.names(mat_np_h)),])
# 
pdf("Number_of_expressed_NPs_per_cell_type_with_human_data.pdf", 10,10)
par(mfcol=c(1,2), mar=c(4,6,4,4), las=2)
barplot(mat_np[order(mat_np[,1]),], beside = T, horiz = T, names.arg = row.names(mat_np)[order(mat_np[,1])], 
        main="Number of NPs (mouse lung)", xlab="total=48")
barplot(mat_np_h[order(mat_np_h[,1]),], beside = T, horiz = T, names.arg = row.names(mat_np_h)[order(mat_np_h[,1])], 
        main="Number of NPs (human brain)", xlab="total=73")
dev.off()
# ERCC plots ------------------------------------------------------------
ercc_table <- read.csv("ERCC_concentrations.csv",header=T, row.names=1)
data_ercc_adult <- data_ercc[,colnames(data_adult)]
data_ercc_adult_norm <- apply(X=data_ercc_adult, MARGIN=2, function(x) log2((((x/sum(x))*1000000)+1)) )
# 
meta_ercc <- matrix(nrow=nrow(data_ercc_adult_norm), ncol=2)
colnames(meta_ercc) <- c("CinA", "CinB")
row.names(meta_ercc) <- row.names(data_ercc_adult_norm)
# 
for(i in 1:nrow(meta_ercc)){
  a <- which(row.names(ercc_table)==row.names(meta_ercc)[i])
  meta_ercc[i,"CinA"] <-  ercc_table[a,"CinA"]
  meta_ercc[i,"CinB"] <- ercc_table[a,"CinB"]
}
meta_ercc <- as.data.frame(meta_ercc)
plot(meta_ercc[,"CinA"],rowMeans(data_ercc_adult_norm), log="x")
# 
meta_ercc[,"non_zero_cells"] <- apply(data_ercc_adult_norm, MARGIN = 1, function(x) length(which(x!=0)))/ncol(data_ercc_adult_norm)

#plot(log10(meta_ercc[,"CinA"]),rowMeans(data_ercc_adult_norm), pch=19, col="red")
lo <- lowess(rowMeans(data_ercc_adult_norm)~log10(meta_ercc[,"CinA"]), f = 0.3)
# plot(lo$x,lo$y, pch=19, col="red", type="n")
par(mar=c(5,5,5,7))
plot(rowMeans(data_ercc_adult_norm)~log10(meta_ercc[,"CinA"]), pch=19, cex=0.2, col="red", 
     xlab="Log10_ERCC_concentration", ylab="Mean ERCC log2 CPM")
abline(h=seq(from=0,to=17,by=1), col="grey85")
lines(lo, col="red", lwd=1) # lowess line (x,y)
par(new=T)
#plot(log10(meta_ercc[,"CinA"]),meta_ercc[,"non_zero_cells"], log="x", pch=19, col="blue", ylim=c(1,0), yaxt="n")
lo <- lowess(meta_ercc[,"non_zero_cells"]~log10(meta_ercc[,"CinA"]), f=0.3)
plot(meta_ercc[,"non_zero_cells"]~log10(meta_ercc[,"CinA"]), pch=19, col="blue", ylim=c(1,0), yaxt="n", cex=0.2, 
     xlab="", ylab="")
abline(h=seq(0,1,0.1), lty=3)
lines(lo, col="blue", lwd=1)
axis(side = 4, at = seq(from=0,to=1,by=0.1))
mtext("Fraction of non zero cells", side=4, line=3, cex.lab=1,las=3, col="black")
# 
plot(rowMeans(data_ercc_adult_norm)~meta_ercc[,"non_zero_cells"], pch=19, col="blue")
lo <- lowess(rowMeans(data_ercc_adult_norm)~meta_ercc[,"non_zero_cells"], f=0.3)
lines(lo)
############################## 


# Annotate metadata for expression of Calca and Scg5 
meta_adult[,"Calca"] <- NA
meta_adult[,"Scg5"] <- NA
# 
for(i in 1:nrow(meta_adult)){
  if(data_norm_adult["Calca",row.names(meta_adult)[i]] !=0 ) {meta_adult[i,"Calca"] <- 1} else {meta_adult[i,"Calca"] <- 0}
  if(data_norm_adult["Scg5",row.names(meta_adult)[i]] !=0 ) {meta_adult[i,"Scg5"] <- 1} else {meta_adult[i,"Scg5"] <- 0}
}
# Create a column for both 
meta_adult[,"Calca_&_Scg5"] <- meta_adult[,"Calca"]+ meta_adult[,"Scg5"]
meta_adult[which(meta_adult[,"Calca_&_Scg5"] !=0),"Calca_&_Scg5"] <- 1
meta_adult[,"Calca_&_Scg5_color"] <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult[,"Calca_&_Scg5"])]
# Store new metadata field 
write.table(meta_adult, file="170305_reAnalysis/updated_metadata_adult.txt")

# Examine significance of othe metadata fields 
wilcox.test(meta_adult$FSC~meta_adult$Calca)
wilcox.test(meta_adult$SSC~meta_adult$Calca)
table_1 <- table(meta_adult$`Calca_&_Scg5`, meta_adult$cell.type)  
table_2 <- rbind((table_1[1,]/colSums(table_1))*100,
                (table_1[2,]/colSums(table_1))*100)
par(las=3, mar=c(10,5,5,5))
barplot(table_2, main="Calca OR Scg5 expressing cells per cell type",
        legend.text = c("Expressing","Not-expressing"))




# Custom heatmaps and barplots  -------------------------------------------
setwd("/Users/spyros/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/170305_reAnalysis/")
# List of all files to read 
files <- list.files("custom_lists/")
files <- files[grep(".txt", files)]
files <- paste("custom_lists/", files, sep="")
classes <- unique(meta_adult$cell.type)
color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
# 
for(i in 1:length(files)){
# Start plotting pdf 
  name <- strsplit(files[i], ".", fixed=T)[[1]][1]
  pdf(paste(name, ".pdf", sep=""),15,15)
# Import list 
table_temp <- read.table(files[i], header=F)  
table_temp[,"Row"] <- NA
# Find genes 
  for(j in 1:nrow(table_temp)){
    a <- which(row.names(data_norm_adult)==as.character(table_temp[j,1]))
    if(length(a) > 0) {table_temp[j,2] <- a}
  }
# Remove genes that were not found 
table_temp <- table_temp[which(is.na(table_temp[,"Row"])==F),]
# Store fraction of expressing cells in the same table for each cell type 
mat_temp <- matrix(nrow=nrow(table_temp), ncol=length(classes))
colnames(mat_temp) <- classes
# Get genes 
genes_temp <- as.character(table_temp[,1])
# Create matrix to store fraction of expressing cells per cell.type 
mat_temp <- matrix(nrow=length(genes_temp), ncol=length(classes))
row.names(mat_temp) <- genes_temp
colnames(mat_temp) <- classes
# Subset data
# Remove dual classified cells 
cells_subset <- row.names(meta_adult)[which(meta_adult$sum_class <= 1)]
data_temp <- data_norm_adult[genes_temp,cells_subset]
meta_temp <- meta_adult[cells_subset, ]
# data_temp <- data_norm_adult[genes_temp,]
# meta_temp <- meta_adult
# Populate matrix for each gene and each cell type 
for(k in 1:length(classes)){
    cells <- row.names(meta_temp)[which(meta_temp$cell.type==classes[k])]
    data_temp_count <- data_temp[,cells]
    mat_temp[,k] <-  apply(data_temp_count, MARGIN = 1, function(x) length(which(x!=0))/length(x))
    if(sum(mat_temp[,k]) !=0){
    heatmap.2(data_temp_count, trace="n", main=classes[k], cexCol=0.5,  col=color.palette, hclustfun=function(x) hclust(d=x, method='ward.D2'), 
                         ColSideColors=meta_temp[cells,"cell.type.color"], Colv=F, Rowv=F)    
    heatmap.2(data_temp_count, trace="n", main=classes[k], cexCol=0.5,  col=color.palette, hclustfun=function(x) hclust(d=x, method='ward.D2'), 
              ColSideColors=meta_temp[cells,"cell.type.color"], Colv=T, Rowv=F)    }
    
  }
# Plot barplots for each gene
par(mfcol=c(4,4),mar=c(4,2,4,2), las=3)
for(k in 1:nrow(mat_temp)){
  xx <- barplot(mat_temp[k,], names.arg = classes, main=table_temp[k,1], ylim=c(0,1))
  text(round(mat_temp[k,],1),x = xx, y =0.3, col="red", cex = 0.5)
}
write.table(cbind(table_temp,mat_temp), file=paste(name, ".csv", sep=""))
dev.off()
}





# Doublet analysis. Latest  -----------------------------------------------
pdf("Doublet_analysis.pdf",15,15)
# Import gene list 
setwd("170305_reAnalysis/")
marker_genes <- read.csv("cell_marker_lists/all_markers.csv", header=T)
# Create matrix to store cell type scores for each cell 
mat_cell_markers <- matrix(nrow=nrow(meta_adult), ncol=length(unique(meta_adult$cell.type)))
colnames(mat_cell_markers) <- unique(meta_adult$cell.type)
row.names(mat_cell_markers) <- row.names(meta_adult)
# Populate the matrix with the cell type score for each cell 
for(i in 1:ncol(mat_cell_markers)){
  # get genes 
  type <- colnames(mat_cell_markers)[i]
  genes <- as.character(marker_genes$Gene[which(marker_genes$cell.type==type)])
  data_temp <- data_norm_adult[genes, row.names(mat_cell_markers)]
  mat_cell_markers[,i] <- colSums(data_temp)/max(colSums(data_temp))
}
# Plot all 
par(mfcol=c(3,4), mar=c(3,3,3,3), las=3)
for(i in 1:ncol(mat_cell_markers)){
  boxplot(mat_cell_markers[,i]~meta_adult$cell.type, main=colnames(mat_cell_markers)[i])
}
# Heatmap plot 
my.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
heatmap.2(mat_cell_markers, trace="n", Colv=F, RowSideColors = meta_adult$cell.type.color, col=my.palette, 
          hclustfun = function(x) hclust(x, "ward.D2"))
# ordered by cell tpe 
names <- row.names(meta_adult)[order(meta_adult$cell.type)]
heatmap.2(mat_cell_markers[names,], trace="n", Colv=F, RowSideColors = meta_adult[names,"cell.type.color"], col=my.palette, 
          hclustfun = function(x) hclust(x, "ward.D2"), Rowv = F)
heatmap.2(mat_cell_markers[names,], trace="n", Colv=F, RowSideColors = meta_adult[names,"cell.type.color"], col=my.palette, 
          hclustfun = function(x) hclust(x, "ward.D2"))
# Plot pair plots 
mat_cell_markers_frame <- as.data.frame(mat_cell_markers)
par(mfcol=c(4,4), mar=c(4,4,4,4))
for(i in 1:ncol(mat_cell_markers)){
  for(j in 1:ncol(mat_cell_markers)){
    plot(mat_cell_markers_frame[,i],mat_cell_markers_frame[,j], 
         xlab=colnames(mat_cell_markers_frame)[i], ylab=colnames(mat_cell_markers_frame)[j],
         pch=19, cex=0.5)
    abline(h=0.5, col="red")
    abline(v=0.5, col="red")
  }
}
# Classify each cell. Minimum cutoff is  0.5 
mat_cell_markers_sum <- apply(mat_cell_markers, 1, function(x) length(which(x > 0.5)))
# Add info to metadata
mat_cell_markers_sum <- mat_cell_markers_sum[row.names(meta_adult)]
meta_adult$sum_class <- mat_cell_markers_sum
# Generate color
meta_adult$sum_class_color <- tableau_color_pal("tableau20")(20)[as.factor(meta_adult$sum_class)]
# Generate and plot table 
table_1 <- table(meta_adult$cell.type, meta_adult$sum_class)
tableplot(table_1)
# 
dev.off()
# Table of dual classified cells 
cells <- row.names(meta_adult)[which(meta_adult$sum_class==2)]
mat_dual <- matrix(nrow=length(cells), ncol=3)
row.names(mat_dual) <- cells
colnames(mat_dual) <- c("class1", "class2", "tsne.class")
for(i in 1:length(cells)){
  mat_dual[i,c(1:2)] <- colnames(mat_cell_markers_frame)[which(mat_cell_markers_frame[cells[i],] > 0.5)]
  mat_dual[i,3] <- meta_adult[cells[i],"cell.type"]
}
write.table(mat_dual, file="Table_of_dual_classified_cells.csv")
# plot of detected genes 
pdf("Number_of_detected_genes_per_classification.pdf",10,10)
boxplot(meta_adult$Genes_detected~meta_adult$sum_class, main="Number of gene detected", 
        xlab="number of classifications", ylab="number of detected genes")
par(las=3)
boxplot(meta_adult$Genes_detected~paste(meta_adult$cell.type, meta_adult$sum_class, sep="."), main="Number of gene detected", 
        xlab="number of classifications", ylab="number of detected genes")
dev.off()
# Extra tables 
mat_dual <- as.data.frame(mat_dual)
mat_dual[,"TSNE_cluster"] <- NA
mat_dual[,"Plate"] <- NA 
for(i in 1:nrow(mat_dual)){
  a <- which(row.names(meta_adult)==row.names(mat_dual)[i])
  if(length(a) > 0)  {mat_dual[i,"TSNE_cluster"] <- meta_adult$Cluster_2d[a]}
  if(length(a) > 0)  {mat_dual[i,"Plate"] <- meta_adult$plate[a]}
}
mat_dual[,"combined_dual"] <- paste(mat_dual$class1, mat_dual$class2, sep = "/")
pdf("Doublet_extra_tables.pdf",15,15)
tableplot(table(mat_dual$combined_dual, mat_dual$TSNE_cluster))
tableplot(table(mat_dual$combined_dual, mat_dual$Plate))
dev.off()
################################################################


# multi-gene ligand receptor analysis -------------------------------------
ligand <- as.character(read.table("NTS_Ach_ligand  .txt", header=F)[,1])
receptor <- as.character(read.table("NTS_Ach_receptors.txt", header=F)[,1])
# Remove dual classified cells 
cells_subset <- row.names(meta_adult)[which(meta_adult$sum_class <= 1)]
data_temp <- data_norm_adult[,cells_subset]
meta_temp <- meta_adult[cells_subset, ]
# Subset ligand and receptor data
# Ligand
ligand_row <- c()
for(i in 1:length(ligand)){
  a <- which(row.names(data_temp)==ligand[i])
  if(length(a) > 0) {ligand_row <- c(ligand_row,a)} else {ligand_row <- c(ligand_row,NA)}
}
ligand <- ligand[which(is.na(ligand_row)==F)]
data_ligand <- data_temp[ligand,]
# Receptor
rec_row <- c()
for(i in 1:length(receptor)){
  a <- which(row.names(data_temp)==receptor[i])
  if(length(a) > 0) {rec_row <- c(rec_row,a)} else {rec_row <- c(rec_row,NA)}
}
receptor <- receptor[which(is.na(rec_row)==F)]
data_receptor <- data_temp[receptor,]
# Caclulate score per cell 
mat_score <- matrix(nrow=nrow(meta_temp), ncol=2)
row.names(mat_score) <- row.names(meta_temp)
colnames(mat_score) <- c("ligand", "receptor")
# Calculate scores 
for(i in 1:nrow(mat_score)){
  mat_score[i,"ligand"]  <- sum(data_ligand[,row.names(mat_score)[i]])
  mat_score[i,"receptor"]  <- sum(data_receptor[,row.names(mat_score)[i]])
}
# Normalize score 
mat_score[,"ligand"] <- mat_score[,"ligand"]/max(mat_score[,"ligand"])
mat_score[,"receptor"] <- mat_score[,"receptor"]/max(mat_score[,"receptor"])
# Plot 
pdf("Ach_ligand_receptor_analysis.pdf",10,10)
boxplot(mat_score[,"ligand"]~meta_temp$cell.type, outline=F, boxwex=0.5, 
        main="Ligand scores")
stripchart(mat_score[,"ligand"]~meta_temp$cell.type,cex=0.5,
           vertical = TRUE, method = "jitter", 
           pch = 21, col = "maroon", bg = "bisque", 
           add = TRUE) 
boxplot(mat_score[,"receptor"]~meta_temp$cell.type, outline=F, boxwex=0.5, 
        main="Receptor scores")
stripchart(mat_score[,"receptor"]~meta_temp$cell.type,cex=0.5,
           vertical = TRUE, method = "jitter", 
           pch = 21, col = "maroon", bg = "bisque", 
           add = TRUE) 
mat_score_temp <- mat_score[names(which(rowSums(mat_score)!=0)),]
heatmap.2(mat_score_temp , trace="n", Colv=F, RowSideColors = meta_temp[row.names(mat_score_temp),"cell.type.color"], col=my.palette, 
          hclustfun = function(x) hclust(x, "ward.D2"), Rowv = T)
plot(mat_score_temp[,"ligand"],mat_score_temp[,"receptor"], 
     pch=19, col=meta_temp[row.names(mat_score_temp),"cell.type.color"], 
     xlab="Ligand score", ylab="Receptor score")
legend("topright", legend = unique(meta_temp[row.names(mat_score_temp),"cell.type"]),
       fill=unique(meta_temp[row.names(mat_score_temp),"cell.type.color"]), cex=0.7)
dev.off()





# Individual cell doublet plots  ------------------------------------------
pdf("Latest_doublet_analysis_plots.pdf",width = 15,height = 10)
# add classification info to metadata 
meta_adult$alt_class <- NA
#meta_adult$class2 <- NA 
for(i in 1:nrow(mat_dual)){
  a <- which(row.names(meta_adult)==row.names(mat_dual)[i])
  if(mat_dual[i,"class1"]!=meta_adult$cell.type[a]) {meta_adult$alt_class[a] <- as.character(mat_dual[i,"class1"])} else {meta_adult$alt_class[a] <- as.character(mat_dual[i,"class2"])}
    }
# Add pulse width to metadata 
temp <- read.table("../PN FACS data-selected/Index_files/1000500604_pulse.txt", header=T)
meta_adult$pulse.width <- NA
for(i in 1:nrow(temp)){
  b <- which(meta_adult$well==as.character(temp$Well[i]))
  meta_adult$pulse.width[b][which(meta_adult$plate[b]=="1000500604")] <- temp$Pulse.width[i]
  }
# Plot pulse width 
meta_temp <- meta_adult[which(meta_adult$plate=="1000500604"),]
box.with.strip(data = meta_temp$pulse.width, group = meta_temp$sum_class, main = "Pulse width for cells of plate 1000500604", xlab = "Number of classifications", 
               ylab = "Pulse width")
# Plot number of genes 
box.with.strip(data = meta_adult$Genes_detected, group = meta_adult$sum_class, main = "Pulse width for all cells", xlab = "Number of classifications", 
               ylab = "Number of detected genes")
# Plot pulse width separately for each cell 
meta_temp <- meta_adult[which(meta_adult$plate=="1000500604"),]
#meta_temp <- meta_temp[c(which(meta_temp$sum_class == 2)),]
meta_temp <- meta_temp[order(meta_temp$pulse.width),]
col_1 <- meta_temp$cell.type
for(i in 1:length(col_1)){col_1[i] <- meta_adult$cell.type.color[which(meta_adult$cell.type==col_1[i])][1]}
col_2 <- meta_temp$alt_class
col_2[which(is.na(col_2)==T)] <- "white"
for(i in 1:length(col_2)){
  a <- meta_temp$cell.type.color[which(meta_temp$cell.type==col_2[i])[1]]
  #a <- colnames(table_col)[which(row.names(table_col)==col_2[i])]
#  if(length(a) > 0) {col_2[i] <- colnames(table_col)[which(row.names(table_col)==col_2[i])]} }
if(length(a) > 0) {col_2[i] <- a}}
col_2[which(is.na(col_2)==T)] <- "white"
# And plot second color 
plot(meta_temp[,"pulse.width"], pch=19, cex=3, col=col_1, xlab="", 
     ylab="Pulse width")
par(new=T)
plot(meta_temp[,"pulse.width"], pch=19, cex=1.5, col=col_2, xlab="", ylab="")
legend("bottomright", legend=unique(meta_adult$cell.type), fill=unique(meta_adult$cell.type.color))
# Same for number of detected genes 
meta_temp <- meta_adult
meta_temp <- meta_temp[order(meta_temp$Genes_detected),]
col_1 <- meta_temp$cell.type
for(i in 1:length(col_1)){col_1[i] <- meta_adult$cell.type.color[which(meta_adult$cell.type==col_1[i])][1]}
col_2 <- meta_temp$alt_class
col_2[which(is.na(col_2)==T)] <- "white"
for(i in 1:length(col_2)){
  a <- colnames(table_col)[which(row.names(table_col)==col_2[i])]
  if(length(a) > 0) {col_2[i] <- colnames(table_col)[which(row.names(table_col)==col_2[i])]} }
plot(meta_temp[,"Genes_detected"], pch=19, cex=1, col=col_1, xlab="", 
     ylab="Genes detected", ylim=c(0, max(meta_temp[,"Genes_detected"] + col_4 + 200, na.rm = T)))
par(new=T)
plot(meta_temp[,"Genes_detected"], pch=19, cex=0.5, col=col_2, xlab="", ylab="",
     ylim=c(0, max(meta_temp[,"Genes_detected"] + col_4 + 200, na.rm = T)))
col_3 <- col_2 
col_3[which(col_3=="white")] <- NA
col_3[which(is.na(col_3)==F)] <- 1
arrows(y0 = meta_temp[,"Genes_detected"]+500+as.numeric(col_3), 
       y1 = meta_temp[,"Genes_detected"]+as.numeric(col_3)+100, 
       x0 = 1:length(col_4), 
       x1 = 1:length(col_4), length = 0.05  )
legend("bottomright", legend=unique(meta_adult$cell.type), fill=unique(meta_adult$cell.type.color), 
       ncol=3, cex=0.5)
# Table of classifications for all cells 
tableplot(table(meta_adult$cell.type, meta_adult$sum_class))

dev.off()



# Receptor ligand breakdown plot  -----------------------------------------
cell.type <- "NE"
ligand <- "Nrtn"
receptor <- c("Ret", "Gfra2")
data <- data_norm_adult
meta <- meta_adult
column.name <- "cell.type"
color.name <- "cell.type.color"
# 
cells <- row.names(meta)[which(meta[,column.name]==cell.type)]
no_cells <- row.names(meta)[which(meta[,column.name]!=cell.type)]
col_cells <- meta[cells,color.name][1]

#
pdf(paste(ligand, "_new_plots_ligand_receptor.pdf", sep=""),15,15)
for(i in 1:length(receptor)){
  mat_l <- matrix(nrow=2, ncol=3, c(1,1,2,3,2,4))
  layout(mat_l)
# Dotchart
    dotchart(sort(data[ligand, cells]), 
           xlab="log2CPM", labels = "", main=ligand, 
           pch=19, col=col_cells)
  
# barplot of receptor expression per cell type 
  temp <- data[receptor[i],]  
  temp[temp!=0] <- 1
  mat_temp <- as.matrix(table(temp, meta[,column.name]))
  par(las=3)
  barplot((mat_temp[2,]/apply(mat_temp,2,sum)), 
          ylim=c(0,1), ylab="fraction of expressing cells", cex.names = 0.7, 
          main=receptor[i])
  # Ligand-receptor scatterplot 
  plot(data[ligand, cells],data[receptor[i], cells], 
       xlab=ligand, ylab=receptor[i], pch=19, col=col_cells)
# table of par and endocrine for cell type 
mat_temp <- cbind(as.data.frame(data[ligand, cells]), as.data.frame(data[receptor[i], cells]))
mat_temp[mat_temp[,1]!=0, 1] <- 2
mat_temp[mat_temp[,2]!=0, 2] <- 1
sum_temp <- apply(mat_temp, 1, sum)
mat_temp <- matrix(nrow=1,ncol=4)
for(i in c(1,2,3,4)){
  mat_temp[1,i] <-  length(which(sum_temp==i-1))
}
colnames(mat_temp) <- c("None", "Receptor only", "Ligand only", "Both")
barplot(mat_temp/length(cells), ylab="Fraction of cells", ylim=c(0,1))
text(x = c(0.7,1.8,3.1,4.3), y = 0.5, labels = round(mat_temp/length(cells), 2), col="red")
}
dev.off()
# 


# TO DO LIST Christin 

# 1
# Heatmaps. Genes pernding. Make for each cell type and only for singlets and doublets 
# Genes are in order and cells are clustered 
# Color 1 is the primary cell type
meta_adult$cell.type
# Color 2 is the number of classifications or 
meta_adult$class1
# Color 2 is the number of classifications or 
meta_adult$class2

# Heatmap of "marker_genes" for all cells that were uniquely classified. Color bar on top is the cell type 
# Color bar on side is cell type as well 
pdf("Marker_genes_uniquely_classified_cells.pdf",15,15)
cells <- row.names(meta_adult)[which(meta_adult$sum_class==1)]
genes <- as.character(marker_genes$Gene)
# Subset and order 
cells <- cells[order(meta_adult[cells,"cell.type"])]
data_to_plot <- data_norm_adult[genes,cells]
col_labels <- as.matrix(meta_adult[cells,"cell.type.color"])
labels <- as.matrix(meta_adult[cells,"cell.type"])
colnames(col_labels) <- c("Cell.type.color")
# Get genes and plot 
heatmap.3(data_to_plot, trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main="Marker genes uniquely classified cells", 
          labCol = "", cexRow = 0.8, ColSideColors = col_labels, margins=c(4,12), ColSideColorsSize=2, 
          Colv=F, Rowv=F)    
# Add legend 
legend("topright",
       legend=unique(labels),
       fill=unique(col_labels), 
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.5)
dev.off()


# NP heatmaps 
pdf("NP_ligands_receptors_3.pdf",15,15)
# Import genes 
genes <- read.table("NP_ligands_receptors_3.txt")
genes$V1 <- as.character(genes$V1)
genes[,"Row"] <- NA 
for(i in 1:nrow(genes)){
  a <- which(row.names(data_norm_adult)==genes$V1[i])
  if(length(a) > 0) {genes$Row[i] <- a}
}
genes <- genes[which(is.na(genes$Row)==F),]
# Subset and order 
cells <- c(row.names(meta_adult)[which(meta_adult$sum_class==1)],
           row.names(meta_adult)[which(meta_adult$sum_class==2)])
data_to_plot <- data_norm_adult[genes$V1,cells]
# Create alt class color labels 
meta_adult$class1.color <- "white"
meta_adult$class2.color <- "white"
# 
for(i in 1:nrow(meta_adult)){
  col <- meta_adult$cell.type.color[which(meta_adult$cell.type==meta_adult$class1[i])[1]]
  if(is.na(col)==F) {meta_adult$class1.color[i] <- col }
  col <- meta_adult$cell.type.color[which(meta_adult$cell.type==meta_adult$class2[i])[1]]
  if(is.na(col)==F) {meta_adult$class2.color[i] <- col }
}
col_labels <- cbind(meta_adult[cells,"cell.type.color"], 
                    meta_adult[cells,"class1.color"],
                    meta_adult[cells,"class2.color"])
row.names(col_labels) <- row.names(meta_adult[cells,])
colnames(col_labels) <- c("Primary_cell_type", "class1", "class2")
# Get genes and plot 
meta_to_plot <- meta_adult[cells,]
types <- unique(meta_to_plot$cell.type)
for(i in 1:length(types)){
  cells <- row.names(meta_to_plot)[which(meta_to_plot$cell.type==types[i])]
  col <- col_labels[cells,]
  data_temp <- data_to_plot[,cells]
# Plot   
heatmap.3(data_temp, trace='n', col=color.palette, distfun = function(x) dist(x, "euclidean"), 
          hclustfun = function(x) hclust(x, 'ward.D2'), main=types[i], 
          labCol = "", cexRow = 0.8, ColSideColors = col, margins=c(4,12), ColSideColorsSize=2, 
          Colv=T, Rowv=F)    
# write table with data used to make heatmap 
write.table(data_temp, file=paste(types[i],"NP_ligands_receptors_3.csv", sep="_"))
# Add legend 
legend("topright",
       legend=unique(meta_adult$cell.type),
       fill=unique(meta_adult$cell.type.color), 
       border=FALSE, bty="n", y.intersp = 0.6, cex=0.5)
}
dev.off()


# Similarity of Nbm expressing NE cells 
cells <- colnames(data_norm_adult)[which(data_norm_adult["Nmb",] != 0)]
cells <- cells[which(meta_adult[cells,"cell.type"]=="NE")]
cells_2 <- row.names(meta_adult)[which(meta_adult$cell.type!="NE")]
dist_less <- as.matrix(dist.cor.adult)[cells,cells_2]
# For cell 1 
neig_1 <- names(dist_less[1,][order(dist_less[1,], decreasing = F)][1:10])
neig_2 <- names(dist_less[2,][order(dist_less[2,], decreasing = F)][1:10])
table(meta_adult[neig_1,"cell.type"])
table(meta_adult[neig_2,"cell.type"])
/-     


# Joy plots 
gene <- "Calca"
cells <- row.names(meta_adult)[which(meta_adult$sum_class==1)]
frame <- cbind(as.data.frame(data_norm_adult[gene,cells]), meta_adult[cells,cell.type])
colnames(frame) <- c(gene, "type")
library(ggplot2)
library(gg1joy)
ggplot(frame, aes(x = Calca, y = type)) + geom_joy() + theme_bw() 


  
# Create table of fractional expression 
clusters <- unique(meta_adult$cell.type)
mat_exp <- matrix(nrow=nrow(data_norm_adult), ncol=length(clusters))
mat_mean <- matrix(nrow=nrow(data_norm_adult), ncol=length(clusters))
for(i in 1:nrow(data_norm_adult)){
  for(j in 1:length(clusters)){
    cells <- row.names(meta_adult)[which(meta_adult$cell.type==clusters[j])]
    mat_mean[i,j] <- median(data_norm_adult[i,cells])
    mat_exp[i,j]  <- length(which(data_norm_adult[i,cells]>0))/length(cells)
  }
}
row.names(mat_exp) <- row.names(data_norm_adult)
colnames(mat_exp) <- clusters
row.names(mat_mean) <- row.names(data_norm_adult)
colnames(mat_mean) <- clusters

write.table(mat_exp, file="~/Desktop/Lung_fraction_expressing.csv") 
write.table(mat_mean, file="~/Desktop/Lung_median_expressing.csv") 


# Fraction of expressing cells 
################
#types <- unique(meta_adult$cell.type)
mat_frac <- matrix(nrow=nrow(data_norm_adult), ncol=2)
row.names(mat_frac) <- row.names(data_norm_adult)
colnames(mat_frac) <- c("NE", "nonNE")
# for(i in 1:ncol(mat_frac)){
  cells <- row.names(meta_adult)[which(meta_adult$cell.type=="NE")]
  data.temp <- data_norm_adult[,cells]
  mat_frac[,1]  <- apply(data.temp, MARGIN = 1, function(x) length(which(x != 0))  )/length(cells)
  noncells <- row.names(meta_adult)[which(meta_adult$cell.type!="NE")]
  data.temp <- data_norm_adult[,noncells]
  mat_frac[,2]  <- apply(data.temp, MARGIN = 1, function(x) length(which(x != 0))  )/length(noncells)
  mat_frac <- as.data.frame(mat_frac)
# Add two columns 
# Cond.1 = genes that are expressed in greater than 90% of NE cells and fewer than 10% of of all other clusters
# Cond.2 =  list of genes expressed in none of other cell types and greater 10 percent of NE cells
genes <- row.names(mat_frac)[which(mat_frac$NE > 0.85)]
mat.temp <- mat_frac[genes,]
mat.temp[which(mat.temp[,"nonNE"] < 0.15),]

# 
genes <- row.names(mat_frac)[which(mat_frac$nonNE == 0)]
mat.temp <- mat_frac[genes,]
mat.temp[which(mat.temp[,"NE"] > 0.1),]


################

save(list = c("data_norm_adult", "meta_adult", "data_adult", "my.tsne.2d.adult", "mat_wilcox_clus_between"), file="~/Documents/Labwork/Projects/Stanford_projects/Single_cells/Kuo_analysis/170305_reAnalysis/Manuscript/NE_adult_Lung_workspace.RData")




# MACA Immune cell analysis 
################
meta.maca$type.color <- tableau_color_pal("tableau10")(4)[meta.maca$type]
pdf("Ramp1_Calcrl_expression_MACA_immune.pdf",10,10)
layout(matrix(nrow=2, ncol=1, c(1,2)), heights = c(3,1))
par(mar=c(4,4,1,3))
plot(data.maca.norm["Calcrl",], data.maca.norm["Ramp1",], pch=19, col=meta.maca$type.color, xlab="Calcrl", ylab="Ramp1")
legend("topright", legend = unique(meta.maca$type), fill = unique(meta.maca$type.color), cex = 0.7)
# Find co-expressing cells 
temp <- cbind(data.maca.norm["Calcrl",], data.maca.norm["Ramp1",])
temp[temp != 0] <- 1
cells <- names(which(rowSums(temp) == 2))
par(mar=c(4,4,1,3))
barplot(table(meta.maca[cells,"type"]), col=unique(meta.maca$type.color), ylab = "Ncells", xpd = F)
layout(matrix(nrow=2, ncol=1, c(1,2)), heights = c(1,1))
par(mar=c(5,4,1,1))
box.with.strip(data = data.maca.norm["Calcrl",], group = meta.maca$type, point.colors = "grey90", main="Calcrl", xlab = "", ylab = "Log2CPM")
box.with.strip(data = data.maca.norm["Ramp1",], group = meta.maca$type, point.colors = "grey90", main="Ramp1", xlab = "", ylab = "Log2CPM")
dev.off()
# Compare gene expression between double expressors and not 
temp <- cbind(data.maca.norm["Calcrl",], data.maca.norm["Ramp1",])
temp[temp != 0] <- 1
cells <- which(meta.maca$type=="macrophages_interstitial")
temp <- temp[cells,]
cells.DE <- names(which(rowSums(temp) == 2))
cells.nDE <- names(which(rowSums(temp) != 2)) 
#
cells <- row.names(meta.maca)[which(meta.maca$type=="macrophages_interstitial")]
data.maca.norm.temp <- data.maca.norm[,cells]
data.maca.norm.temp <- data.maca.norm.temp[names(which(rowSums(data.maca.norm.temp)!=0)),]
mat.DE <- matrix(nrow=nrow(data.maca.norm.temp), ncol=1)
row.names(mat.DE) <- row.names(data.maca.norm.temp)
colnames(mat.DE) <- "macrophages_interstitial.DEvsnonDE"
# 
for(i in 1:nrow(mat.DE)){
  test <- wilcox.test(data.maca.norm.temp[i,cells.DE], data.maca.norm.temp[i,cells.nDE])
  a <- test$p.value*nrow(data.maca.norm.temp)
  if(a > 1) {mat.DE[i,1] <- 1 } else {mat.DE[i,1] <- a }
}

  
require("Seurat")

?CreateSeuratObject()

####### Julien analysis 
genes <- c("Calca", "Sftpc", "Scgb1a1")
cells <- row.names(meta_adult)[which(meta_adult$sum_class==1)]
data.temp <- data_norm_adult[genes,cells]
data.temp[data.temp !=0] <- 1
sums.all <- colSums(data.temp)
table(sums.all,meta_adult[cells,"cell.type"])

# 
require(dplyr)
meta_adult %>% group_by(cell.type) %>% summarize(mean=mean(Genes_detected), median=median(Genes_detected))




length(which(data_norm_adult["Calca",which(meta_adult$cell.type=="NE")]==0))/length(which(meta_adult$cell.type=="NE"))








