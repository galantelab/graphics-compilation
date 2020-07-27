#MA plot ###################################################################

#Load library
library("ggpubr") #function ggmaplot

#Read file with log2(mean TPM expression + 1) (baseMean), deltaPSI (log2FoldChange) and FDR values (padj) of transcripts
df = read.table("data_MA_plot.txt")

#Make MA plot
pdf("MAplot.pdf", width=3, height=2.5)
ggmaplot(df, main = "MA plot",
         fdr = 0.05, fc = 0.2, size = 0.4,
         palette = c("cadetblue3", "coral2", "darkgray"),
         genenames = as.vector(df$transcript.symbol),
         legend = "bottom", top = 0,
         font.label = c("bold", 10),
         font.main = "bold",
         xlab="log2(mean TPM + 1)",
         ylab="deltaPSI",
      ggtheme = ggplot2::theme_minimal()) + 
      scale_y_continuous(limits=c(-1,1)) + 
      theme(plot.title = element_text(size=9, hjust = 0.5), 
            axis.text.x = element_text(size=8), 
            axis.text.y = element_text(size=8), 
            axis.title.x = element_text(size=9), 
            axis.title.y = element_text(size=9), 
            legend.text = element_text(size=8))
dev.off()
