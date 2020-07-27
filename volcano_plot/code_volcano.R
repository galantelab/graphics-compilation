#Volcano plot ##############################################################

#loading libraries
library(ggplot2)
library(ggrepel)

#Read file (diff expressed genes |log2FC|>2 & FDR<0.05)
df = read.table("data_volcano.txt", header=T)

#Log-transform FDR values (-log10)
df$FDR = log(df$FDR,10)*(-1)

#Add column with direction
df$direction = "down"
df[df$log2FC>0,]$direction <- "up"

#Make plot
pp=ggplot(df, aes(log2FC, FDR)) + geom_point(aes(col = direction), size = 1, alpha=0.5) +
       scale_color_manual(values = c("up" = "dodgerblue3", "down" = "firebrick3"),
       breaks = c("down", "up"), labels = c("down-regulated", "up-regulated")) + 
  xlab("log2(fold change)") +
  ylab("- log10(FDR)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10)) +
  geom_label_repel(data=subset(df, (FDR>100 & abs(log2FC)>4)), aes(label=gene), 
                   size = 2, box.padding = unit(1.0, "lines"), label.padding = unit(0.2, "lines"),
                   fill="gray", color="white", segment.size  = 0.1, segment.color = "black")

ggsave("volcano_plot1.pdf", width=4, height=3)

#Set specific genes to be marked
df$labels = FALSE
df[df$gene=="PRLHR" | df$gene=="SAA2" | df$gene=="TMSB4XP1" | df$gene=="TMSB4XP2",]$labels <- TRUE

#Make plot
pp=ggplot(df, aes(log2FC, FDR)) + geom_point(aes(col = direction), size = 1, alpha=0.5) +
       scale_color_manual(values = c("up" = "dodgerblue3", "down" = "firebrick3"),
       breaks = c("down", "up"), labels = c("down-regulated", "up-regulated")) + 
  xlab("log2(fold change)") +
  ylab("- log10(FDR)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "top", legend.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10)) +
  geom_label_repel(aes(label = ifelse(df$labels, as.character(df$gene),"")), 
                   size = 2, box.padding = unit(1.0, "lines"), label.padding = unit(0.2, "lines"),
                   fill="gray", color="white", segment.size  = 0.1, segment.color = "black")
ggsave("volcano_plot2.pdf", width=4, height=3)
