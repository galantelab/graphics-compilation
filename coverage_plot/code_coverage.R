#kpPlotCoverage ############################################################

#Pipeline before: 
 #1. Run bedtools makewindows to create bed file containing adjacent genomic windows of 10000bp for all chromosomes in hg38
 #2. Run samtools depth for BAM file at each region of bed file
 #3. Process result of samtools depth to get "repeated genomic regions" based on coverage of each region

#Load library
library(karyoploteR) #functions plotKaryotype and kpPlotCoverage

#Read file
data = read.table("data_coverage.txt", header=T)

#Convert dataframe to GenomicRanges object
regions = makeGRangesFromDataFrame(data)

#Make plot with all chromosomes
pdf("plot_coverage.pdf")
kp <- plotKaryotype(genome="hg38", main="Plot coverage")
kpPlotCoverage(kp, data=regions)
garbage=dev.off()

#Make plot with some chromosomes
pdf("plot_coverage2.pdf", width=5, height=3)
kp <- plotKaryotype(genome="hg38", main="Plot coverage", chromosomes=c("chr4","chr17","chrX"))
kpPlotCoverage(kp, data=regions)
garbage=dev.off()
