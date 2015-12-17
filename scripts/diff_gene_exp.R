##Differential gene expression analysis using DESeq2

#Load libraries
library(DESeq2)
library("RColorBrewer")
library("gplots")
library("ggplot2")
library("cowplot")
library(tikzDevice)

#set working directory
setwd("D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\htseq_counts_union_newnames")

##import all htseq-count files
count_files <- list.files(pattern="*.txt", full.names=FALSE)
list_cf <- lapply(count_files,function(i){read.csv(file=i,header=FALSE, sep="\t")})

num_genes <- nrow(list_cf[[1]])
num_samples <- length(list_cf)
columns <- c("Gene","Count")

list_cf <- lapply(list_cf,function(i){
  colnames(i) <- columns
  return(i)
})
##make count matrix
#make matrix of all counts for each sample and add gene names for the rows
all_counts <- matrix(unlist(lapply(list_cf, function(i){i$Count})),nrow = num_genes, ncol = num_samples)
rownames(all_counts) <- list_cf[[1]]$Gene

#remove last 5 rows (htseq-count special rows such as no_feature etc.)
#could bypass this step by using built-in deseq function to import from htseq but \_0_/
all_counts <- all_counts[1:((dim(all_counts)[1])-5),]

#here, we can generate an alternate all_counts matrix which subsets the data to remove the shandong re-watered samples
#I do this because the SRW samples are in a separate batch and I want to see if they are inflating the batch effects
all_counts_no_srw <- all_counts[,-c(3,4,5,6)]

##set up metadata
parse_names <- strsplit(count_files, split="_")
parse_names <- matrix(unlist(parse_names), nrow=18, ncol=5, byrow=T)
parse_names_no_srw <- matrix(unlist(parse_names), nrow=14, ncol=5, byrow=T)

#add batch information
parse_names <- cbind(parse_names,c("2","2","4","4","4","4","2","3","3","1","3","3","1","3","3","1","3","3"))
parse_names_no_srw <- cbind(parse_names_no_srw,c("2","2","2","3","3","1","3","3","1","3","3","1","3","3"))
col_names_counts <- paste(parse_names[,1], "_", parse_names[,2], "_", parse_names[,3], "_", parse_names[,4], "_", parse_names[,5], "_", parse_names[,6], sep="")
colnames(all_counts) <- col_names_counts

col_names_counts_no_srw <- paste(parse_names_no_srw[,1], "_", parse_names_no_srw[,2], "_", parse_names_no_srw[,3], "_", parse_names_no_srw[,4], "_", parse_names_no_srw[,5], "_", parse_names_no_srw[,6], sep="")
colnames(all_counts_no_srw) <- col_names_counts_no_srw

#remove all genes which had no expression in all of the samples
all_counts <- subset(all_counts,rowSums(all_counts) != 0)
all_counts_no_srw <- subset(all_counts_no_srw,rowSums(all_counts_no_srw) != 0)

#generate metadata data frame
experimental_design = data.frame(
  sample_names = col_names_counts,  # sample name
  ecotype = factor(parse_names[,1]), # shandong or yukon ecotype
  treatment = factor(parse_names[,2]),  # WW, D, or RW for well-watered, drought, or re-watered
  replicate = factor(parse_names[,3]),   # replicate
  lane = factor(parse_names[,4]),      # lane on the Illumina flowcell.
  batch = factor(parse_names[,6])      # which batch of sequencing runs
)

experimental_design_no_srw = data.frame(
  sample_names = col_names_counts_no_srw,  # sample name
  ecotype = factor(parse_names_no_srw[,1]), # shandong or yukon ecotype
  treatment = factor(parse_names_no_srw[,2]),  # WW, D, or RW for well-watered, drought, or re-watered
  replicate = factor(parse_names_no_srw[,3]),   # replicate
  lane = factor(parse_names_no_srw[,4]),      # lane on the Illumina flowcell.
  batch = factor(parse_names_no_srw[,6])      # which batch of sequencing runs
)
##testing batch effects
test_batch_effects <- DESeqDataSetFromMatrix(all_counts, experimental_design, 
                                            design = formula(~ batch))


test_batch_effects <- DESeq(test_batch_effects)

test_batch_effects_results <- results(test_batch_effects, pAdjustMethod="BH")
summary(test_batch_effects_results) #~23% of genes showing a batch effect when SRW is left in

test_batch_effects_no_srw <- DESeqDataSetFromMatrix(all_counts_no_srw, experimental_design_no_srw, 
                                             design = formula(~ batch))


test_batch_effects_no_srw <- DESeq(test_batch_effects_no_srw)

test_batch_effects_results_no_srw <- results(test_batch_effects_no_srw, pAdjustMethod="BH")
summary(test_batch_effects_results_no_srw) #~1.3% genes showing a batch effect when SRW is removed
#this suggests that the batch effect is likely exaggerated because of the lack of mixing of samples in batch 4
#there are multiple ways to approach this, but I will remove the genes showing a batch effect in 
#the absence of SRW, and then check for genes which are interesting later on to see if they are any
#of the genes found in "summary(test_batch_effects_results)"

#taking a look at dispersion estimates and some cluster analysis
plotDispEsts(test_batch_effects) 
plotDispEsts(test_batch_effects_no_srw) 

#rlog transformation in a blind fashion (blind so it doesn't weight the transformation by any of the other varaiables)
for_pca <- rlog(test_batch_effects, blind=TRUE) 
for_pca_no_srw <- rlog(test_batch_effects_no_srw, blind=TRUE) 

#PCA analysis, looking at clustering by batch and by treatment and ecotype
plotPCA(for_pca, intgroup=c("batch"))
plotPCA(for_pca_no_srw, intgroup=c("batch"))
plotPCA(for_pca, intgroup=c("treatment","ecotype"))

#PCA plot to ggplot and to tikzdevice (optional)
#tikz(file = "D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\tex\\pca_batch.tex")
pca_batch <- plotPCA(for_pca, intgroup=c("batch"), returnData=TRUE)
percentVar <- round(100 * attr(pca_batch, "percentVar"))
ggplot(pca_batch, aes(PC1, PC2, color=batch)) +
  geom_point(size=3.5, alpha=0.65) +
  xlab(paste0("PC1: ",percentVar[1]," variance")) +
  ylab(paste0("PC2: ",percentVar[2]," variance")) +
  theme_bw()
#dev.off()

#PCA plot to ggplot and to tikzdevice (optional)
#tikz(file = "D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\tex\\pca_data.tex")
pca_data <- plotPCA(for_pca, intgroup=c("treatment","ecotype"), returnData=TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
ggplot(pca_data, aes(PC1, PC2, color=treatment, shape=ecotype)) +
  geom_point(size=3.5, alpha=0.65) +
  xlab(paste0("PC1: ",percentVar[1]," variance")) +
  ylab(paste0("PC2: ",percentVar[2]," variance")) +
  theme_bw()
#dev.off()

#Making matrix of distances to make a heatmap
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix

rownames(mat) <- colnames(mat) <- with(colData(test_batch_effects), paste(ecotype, treatment, batch, sep=" : "))
hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours

#make heatmap and output to tikzdevice(optional)
#tikz(file = "D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\tex\\heatmap.tex")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
#dev.off()

##Remove genes identified from the batch analysis without SRW samples and regenerate deseq object
##for real analysis

#remove gene
genes_to_remove <- which(test_batch_effects_results_no_srw$padj < 0.1)
all_counts_rg <- all_counts[-genes_to_remove,]
dim(all_counts_rg) #make sure everything is okay...

#re generate deseq object and perform real analysis now
effects <- DESeqDataSetFromMatrix(all_counts_rg, experimental_design, 
                                             design = formula(~ ecotype + treatment + ecotype:treatment))

#re-order factors 
effects$treatment <- factor(effects$treatment, levels=c("WW","D","RW"))
effects$ecotype <- relevel(effects$ecotype,"Y")

effects <- DESeq(effects)

#dispersion estimates and MA plot of effects (plotMA is looking at last term in results)
plotDispEsts(effects)
plotMA(effects)
resultsNames(effects)

##contrast analysis
#overall effect of ecotype
s_vs_y <- results(effects, contrast=list(c("ecotype_S_vs_Y")),pAdjustMethod="BH",alpha = 0.01)
#overall effect of drought and re-water
d_vs_ww<- results(effects, contrast=list(c("treatment_D_vs_WW")),pAdjustMethod="BH",alpha = 0.01)
rw_vs_ww<- results(effects, contrast=list(c("treatment_RW_vs_WW")),pAdjustMethod="BH",alpha = 0.01)
#differential reaction of ecotype to drought and re-water treatment
s_vs_y_d <- results(effects, contrast=list(c("ecotypeS.treatmentD")),pAdjustMethod="BH",alpha = 0.01)
s_vs_y_rw <- results(effects, contrast=list(c("ecotypeS.treatmentRW")),pAdjustMethod="BH",alpha = 0.01)

#change all NA's to 1
s_vs_y$padj[is.na(s_vs_y$padj)] <- 1  
d_vs_ww$padj[is.na(d_vs_ww$padj)] <- 1   
rw_vs_ww$padj[is.na(rw_vs_ww$padj)] <- 1   
s_vs_y_d$padj[is.na(s_vs_y_d$padj)] <- 1  
s_vs_y_rw$padj[is.na(s_vs_y_rw$padj)] <- 1   

summary(s_vs_y)  
summary(d_vs_ww)
summary(rw_vs_ww)
summary(s_vs_y_d)
summary(s_vs_y_rw)

#mean log2foldchanges between ecotype and treatment comparisons
mean(abs(s_vs_y$log2FoldChange[s_vs_y$padj < 0.01]))
mean(abs(d_vs_ww$log2FoldChange[d_vs_ww$padj < 0.01]))

plotCounts(effects, gene=which.min(s_vs_y$padj), intgroup=c("ecotype")) #1002598 is a histidine biosynthesis gene
plotCounts(effects, gene=which.min(d_vs_ww$padj), intgroup=c("treatment"))
plotCounts(effects, gene=which.min(rw_vs_ww$padj), intgroup=c("treatment"))
plotCounts(effects, gene=which.min(s_vs_y_d$padj), intgroup=c("treatment","ecotype"))
plotCounts(effects, gene=which.min(s_vs_y_rw$padj), intgroup=c("treatment","ecotype"))

#lets get the top 20 genes (by p-value) in both interaction contrasts
cool_sy_genes <- rownames(s_vs_y[with(s_vs_y, order(padj)), ])[1:20]
cool_dww_genes <- rownames(d_vs_ww[with(d_vs_ww, order(padj)), ])[1:20]
cool_rwww_genes <- rownames(rw_vs_ww[with(rw_vs_ww, order(padj)), ])[1:20]
cool_drought_genes <- rownames(s_vs_y_d[with(s_vs_y_d, order(padj)), ])[1:20]
cool_rewater_genes <- rownames(s_vs_y_rw[with(s_vs_y_rw, order(padj)), ])[1:20]

#top 20 genes of each contrast. we won't look at all of them though.
write(cool_sy_genes, file="../cool_sy_genes.txt")
write(cool_dww_genes, file="../cool_dww_genes.txt")
write(cool_rwww_genes, file="../cool_rwww_genes.txt")
write(cool_drought_genes, file="../cool_drought_genes.txt")
write(cool_rewater_genes, file="../cool_rewater_genes.txt")
#in system, do perl -pi -e 's/\.v1\.0//' file.txt to get the gene names comaptible with the annotation file
#simple grep -f with each file against the annotation file to get the gene description

#here we can get the lfcs for those genes we just obtained and find out how those genes differ.
#make dataframe of gene id we are interested in and the lfc
gene_lfc_s_vs_y_d <- data.frame(cbind(row.names(data.frame(s_vs_y_d[row.names(s_vs_y_d) %in% cool_drought_genes,])),data.frame(s_vs_y_d[row.names(s_vs_y_d) %in% cool_drought_genes,])$log2FoldChange))

#dehydrin and abscisic acid genes
d_a_genes <- c("Thhalv10000340m.v1.0","Thhalv10001885m.v1.0","Thhalv10001888m.v1.0","Thhalv10003251m.v1.0","Thhalv10004906m.v1.0","Thhalv10008313m.v1.0","Thhalv10008706m.v1.0","Thhalv10010821m.v1.0","Thhalv10011493m.v1.0","Thhalv10011506m.v1.0","Thhalv10012271m.v1.0","Thhalv10015597m.v1.0","Thhalv10019078m.v1.0","Thhalv10019152m.v1.0","Thhalv10023774m.v1.0","Thhalv10023775m.v1.0","Thhalv10024390m.v1.0","Thhalv10024393m.v1.0","Thhalv10025203m.v1.0","Thhalv10026303m.v1.0","Thhalv10026916m.v1.0")
s_y_d_dehydrin_abs_df <- data.frame(s_vs_y_d[row.names(s_vs_y_d) %in% d_a_genes,])
d_ww_dehydrin_abs_df <- data.frame(d_vs_ww[row.names(d_vs_ww) %in% d_a_genes,])

#which of those genes show an lfc in these two comparisons
s_y_d_dehydrin_abs_df[s_y_d_dehydrin_abs_df$padj < 0.01,]
d_ww_dehydrin_abs_df[d_ww_dehydrin_abs_df$padj < 0.01,]
