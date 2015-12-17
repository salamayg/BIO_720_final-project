library(ggplot2)
library(reshape2)
library(ggthemes)
library(scales)
library(tikzDevice)

#samples <- c("SD1-1","SD1-2","SRW-1_L001","SRW-1_L002","SRW-2_L001","SRW-2_L002","SWW-1","SWW-3_L001","SWW-3_L002","YD1-1","YD1-2_L001","YD1-2_L002","YRW-1","YRW-2_L001","YRW-2_L002","YWW-1","YWW-2_L001","YWW-2_L002")

#samples <- c("SD1-1_L002","SD1-2_L002","SRW-1_L001","SRW-1_L002","SRW-2_L001","SRW-2_L002","SWW-1_L002","SWW-3_L001","SWW-3_L002","YD1-1_L003","YD1-2_L001","YD1-2_L002","YRW-1_L003","YRW-2_L001","YRW-2_L002","YWW-1_L005","YWW-2_L001","YWW-2_L002")
samples <- c("SD11L002","SD12L002","SRW1L001","SRW1L002","SRW2L001","SRW2L002","SWW1L002","SWW3L001","SWW3L002","YD11L003","YD12L001","YD12L002","YRW1L003","YRW2L001","YRW2L002","YWW1L005","YWW2L001","YWW2L002")


pretrim <- c(146275746,166368296,74431708,75070484,69364348,70037384,104683802,
             70074744,71042620,152809518,74713192,75746140,154340426,68067686,
             68970176,157744646,71467056,72467684)
pretrim <- pretrim / 2  
posttrim <- c(141028918,160240026,72129342,72708426,67072990,67681848,101598814,
              68816462,69747432,150633528,73463452,74454464,152215528,66649426,
              67515770,154522436,70053738,71012292)
posttrim <- posttrim / 2

pretrim <- data.frame(samples,pretrim)
posttrim <- data.frame(samples,posttrim)

colnames(pretrim) <- c("Sample","Reads")
colnames(posttrim) <- c("Sample","Reads")

pretrim$name <- "Pre-trimming"
posttrim$name <- "Post-trimming"

readCounts <- rbind(pretrim,posttrim)

readCounts <- transform(readCounts, 
                          Sample = reorder(Sample, Reads))

readCounts$name <- as.factor(readCounts$name)
#readCounts$name <- factor(readCounts$name, levels = rev(levels(readCounts$name)))

tikz(file = "D:\\Dropbox\\School\\Graduate School\\Masters\\Grad_Courses\\BIO720\\final_project\\trim_stats.tex")
p <- ggplot(readCounts, aes(Sample, Reads, fill = name)) + geom_bar(width=0.8, position = "dodge",stat="identity", alpha=0.65)
p + coord_flip() + 
  theme_bw() + 
  scale_y_continuous(labels = comma,expand=c(0,0)) + 
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black")) +
  labs(fill="", y="Paired Reads") 
dev.off()
   #readCounts <- as.data.frame(cbind(samples,pretrim,posttrim))

#colnames(readCounts) <- c("Sample","Pretrim","Posttrim")



