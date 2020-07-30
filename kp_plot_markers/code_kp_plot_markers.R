#kpPlotMarkers ###################################################

#Load libraries
library(biomaRt)
library(regioneR)
library(karyoploteR) 

#Select genes
gene.symbols = c("OR2B2", "OR9A4", "OR2AE1", "OR52K1", "OR2G2", "OR1L6")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes = toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'), filters = 'hgnc_symbol', values = gene.symbols, mart = ensembl))
seqlevelsStyle(genes) <- "UCSC"

#Make plot
pdf("kpPlotMarkers.pdf",width=5, height=5)
kp <- plotKaryotype(genome="hg38", chromosomes = c("chr1", "chr6", "chr7", "chr9", "chr11"), plot.type=1)
kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, text.orientation = "horizontal", r1=0.5, cex=0.8, adjust.label.position = FALSE)
garbage=dev.off()
