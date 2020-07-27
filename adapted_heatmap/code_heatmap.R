#Code adapted heatmap ##################################################

#Load libraries
library(pheatmap) #function pheatmap
library(RColorBrewer) #function colorRampPalette

#Read file (matrix with 0s and 1s)
df = read.delim("data_heatmap1.txt")

#Plot heatmap
pdf("heatmap_v1.pdf",width=4,height=1.5)
pheatmap(df,cluster_cols= F, cluster_rows=F, fontsize=5, fontsize_row = 6, fontsize_col=5, angle_col = "90", main="Hallmarks of Cancer", show_rownames=T, show_colnames=T, color = c("white", "cadetblue4"), cellwidth=6, cellheight=6, legend=F)
dev.off()

#Read file (matrix with 0s and FDRs)
df = read.delim("data_heatmap2.txt")

#Define color degrade
colfunc = colorRampPalette(c("white", "cadetblue"))
colors = colfunc(50)
colors = colors[c(1,15:40)]

#Plot heatmap
pdf("heatmap_v2.pdf",width=4.5,height=1.5)
pheatmap(df,cluster_cols= F, cluster_rows=F, fontsize=5, fontsize_row = 6, fontsize_col=5, angle_col = "90", main="Hallmarks of Cancer", show_rownames=T, show_colnames=T, color = colors,cellwidth=6,cellheight=6, legend_breaks = c(min(df),min(df)+max(df)/2,max(df)), legend_labels = c(min(df), "-log10(FDR)\n", round(max(df))))
dev.off()
