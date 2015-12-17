####Gather general read statistics for trimming and mapping for RNA-seq samples to be better initiated with the data
#we are working with. 
###


setwd("D:\\Dropbox\\School\\Graduate\ School\\Masters\\Grad_Courses\\BIO720\\final_project")

##Load libraries
#general
library(reshape2)
library(doBy)
#plotting and export
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(tikzDevice)

##read files in and divide for paired read values
pretrim <- read.csv("pretrim_statistics.txt", header=FALSE)
pretrim <- pretrim / 2  

posttrim <- read.csv("posttrim_statistics.txt", header=FALSE)
posttrim <- posttrim / 2

#sample names are in mapped$V1
mapped <- read.csv("mapped_statistics.txt", header=FALSE)



##Generate dataframe with all values(readcounts) and "name" field for pre-trim, post-trim, or mapped
pretrim <- data.frame(mapped$V1,pretrim)
posttrim <- data.frame(mapped$V1,posttrim)


colnames(pretrim) <- c("Sample","Reads")
colnames(posttrim) <- c("Sample","Reads")
colnames(mapped) <- c("Sample","Reads")

#consolidate rows corresponding to same sample on different lanes
pretrim <- aggregate(pretrim$Reads~pretrim$Sample,FUN=sum)
posttrim <- aggregate(posttrim$Reads~posttrim$Sample,FUN=sum)
mapped <- aggregate(mapped$Reads~mapped$Sample,FUN=sum)


colnames(pretrim) <- c("Sample","Reads")
colnames(posttrim) <- c("Sample","Reads")
colnames(mapped) <- c("Sample","Reads")


pretrim$Name <- "Pre-trimming"
posttrim$Name <- "Post-trimming"
mapped$Name <- "Mapped"

readCounts <- rbind(pretrim,posttrim,mapped)

##Reorder data frame by value of reads so the outputted graph has some logic to it
readCounts <- transform(readCounts, 
                        Sample = reorder(Sample, Reads))

readCounts$Name <- as.factor(readCounts$Name)


##Plot graph (and optionally create .tex file of graph for use with latex)
#
#tikz(file = "D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\tex\\trim_stats_nolane.tex")
p <- ggplot(readCounts, order=-as.numeric(name), aes(Sample, Reads, fill = Name)) + geom_bar(width=0.8, position = "dodge",stat="identity", alpha=0.65)
p + coord_flip() + 
  theme_bw() + 
  scale_y_continuous(labels = comma,expand=c(0,0)) + 
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  labs(fill="", y="Paired Reads") +
  scale_fill_manual(values = rev(brewer.pal(3,"Set1")))
#dev.off()


##Find out the percentage of reads mapped 
#
#Make data frame with percentages of post-trimmed reads which mapped (basically, mapped reads normalized to trimmed reads)
percentages <- as.numeric(mapped$Reads / posttrim$Reads)
percentages.df <- data.frame(as.character(mapped$Sample),percentages)
colnames(percentages.df) <- c("Sample", "Percentage")

#Add ecotype factor
percentages.df$Eco <- NA
percentages.df[percentages.df$Sample == grep("^S", percentages.df$Sample, value=TRUE),][, "Eco"] <- "Shandong"
percentages.df[percentages.df$Sample == grep("^Y", percentages.df$Sample, value=TRUE),][, "Eco"] <- "Yukon"

#Get mean percentage of reads mapped for each accession and do a simple t test to see if they are different 
#(not sure if this is a valid test since it is technically a percentage, but roughly, 
#there doens't seem to be a difference in % of mapped reads between ecotypes)
summaryBy(Percentage~Eco, percentages.df,FUN=c(mean,var))
  t.test(percentages.df$Percentage~percentages.df$Eco)

