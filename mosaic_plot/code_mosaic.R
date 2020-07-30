#Load library
library(vcd)

#Read file
df = read.table("data_mosaic.txt", header=T)

#Set seed
set.seed(4711)

#Make mosaic plot
pdf("mosaic_plot.pdf")
mosaic(~ Treatment + Improved, data=df, gp=shading_max)
garbage=dev.off()

#Make association plot
pdf("association_plot.pdf")
assoc(~ Treatment + Improved, data=df, gp=shading_max)
garbage=dev.off()
