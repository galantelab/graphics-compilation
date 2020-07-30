#Code correlation plot #################################################

#Load libraries
library(corrplot) #function corrplot
library(Hmisc) #function rcorr

#Read file
data = read.table("data_correlation.txt")

#Calculate pairwise correlations of all columns (genes)
res = rcorr(as.matrix(data),type="spearman")

#Make hierarchically clustered plot
pdf("correlation_plot.pdf")
corrplot(res$r, p.mat = res$P, insig = "label_sig", sig.level = c(.001, .01, .05), 
         pch.cex = .9, method = "circle", tl.cex = 0.8, tl.col = "black", tl.srt = 45, type = "lower", order = "hclust")
garbage=dev.off()
