library(DESeq2)
library("RColorBrewer")
library("gplots")
setwd("D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\htseq_counts_union_newnames")

count_files <- list.files(pattern="*.txt", full.names=FALSE)
list_cf <- lapply(count_files,function(i){read.csv(file=i,header=FALSE, sep="\t")})

num_genes <- nrow(list_cf[[1]])
num_samples <- length(list_cf)
columns <- c("Gene","Count")

list_cf <- lapply(list_cf,function(i){
  colnames(i) <- columns
  return(i)
})

all_counts <- matrix(unlist(lapply(list_cf, function(i){i$Count})),nrow = num_genes, ncol = num_samples)
rownames(all_counts) <- list_cf[[1]]$Gene



parse_names <- strsplit(count_files, split="_")
parse_names <- matrix(unlist(parse_names), nrow=18, ncol=5, byrow=T)
parse_names <- cbind(parse_names,c("2","2","4","4","4","4","2","3","3","1","3","3","1","3","3","1","3","3"))
col_names_counts <- paste(parse_names[,1], "_", parse_names[,2], "_", parse_names[,3], "_", parse_names[,4], "_", parse_names[,5], "_", parse_names[,6], sep="")
colnames(all_counts) = col_names_counts
#remove all genes which had no expression in all of the samples
all_counts <- subset(all_counts,rowSums(all_counts) != 0)
all_counts$condition <- relevel(dds$condition, ref="untreated")


experimental_design = data.frame(
  sample_names = col_names_counts,  # sample name
  ecotype = factor(parse_names[,1]), # each individual beetle
  treatment = factor(parse_names[,2]),  # small or large beetle
  replicate = factor(parse_names[,3]),   # male or females
  lane = factor(parse_names[,4]),      # Whick lane on the Illumina flowcell.
  batch = factor(parse_names[,6])
)


test_batch_effects <- DESeqDataSetFromMatrix(all_counts, experimental_design, 
                                            design = formula(~ batch))

test_batch_effects$treatment <- relevel(test_batch_effects$treatment, ref=1)
test_batch_effects2 <- DESeq(test_batch_effects)

test_batch_effects2_results <- results(test_batch_effects2)
summary(test_batch_effects2_results) # No evidence, but this is a bit incomplete

plotDispEsts(test_batch_effects2) 

for_pca <- rlog(test_batch_effects2, blind=TRUE) 

plotPCA(for_pca, intgroup=c("batch"))
plotPCA(for_pca, intgroup=c("batch","lane"))



rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix

rownames(mat) <- colnames(mat) <- with(colData(test_batch_effects2), paste(ecotype, treatment, batch, sep=" : "))
hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours

heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))



##real analysis
effects <- DESeqDataSetFromMatrix(all_counts, experimental_design, 
                                  design = formula(~ ecotype + treatment + ecotype:treatment))
effects <- DESeq(effects)

plotDispEsts(effects)
plotMA(effects)
res1 <- results(effects, pAdjustMethod="BH")
head(res1)
summary(res1)
resultsNames(effects)
hist(res1$pvalue)


res_contrast_y_s <- results(effects, 
                            contrast=list(c("ecotype_Y_vs_S")), 
                            pAdjustMethod="BH")

idx <- identify(res_contrast_y_s$baseMean, res$log2FoldChange)
rownames(res)[idx]

plotMA(res_contrast_y_s, main="DESeq2", MLE=TRUE)
