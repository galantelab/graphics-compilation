#Downloaded ENCODE data

  #Context info
  #wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz

  #Histone markers H3k4me3 and H3k36me3
  #wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig
  #wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig

  #DNA binding protein EZH2
  #wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562Ezh239875StdSig.bigWig

#Load libraries
library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#Read K562 chromatin state segmentation file
K562.hmm=read.table("wgEncodeBroadHmmK562HMM.bed.gz", comment.char="")
colnames(K562.hmm)[1:3] = c("chr","start","end") #rename columns
K562.hmm = makeGRangesFromDataFrame(K562.hmm,keep.extra.columns=T) #load into GRanges object

#Define region to plot
TP53.region = toGRanges("chr17:7564422-7602719")

#Create karyotype plot (zoom on defined region)
kp = plotKaryotype(zoom = TP53.region, cex=2)

#Get genes
genes.data = suppressMessages(makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene, karyoplot=kp, plot.transcripts=TRUE, plot.transcripts.structure = TRUE))

#Add gene names
genes.data = addGeneNames(genes.data)

#Merge transcripts of each gene into one
genes.data = mergeTranscripts(genes.data)

#Set histone marks and DNA-binding proteins to be plotted
histone.marks = c(H3K4me3="wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig", H3K36me3="wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig")
DNA.binding = c(EZH2="wgEncodeBroadHistoneK562Ezh239875StdSig.bigWig")

#Set plot parameters
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$data1outmargin <- 0

#Make plot
pdf("plot_epigenetics.pdf",width=15)
kp = plotKaryotype(zoom = TP53.region, cex=1.5, plot.params = pp) #plot karyotype (zoom on defined region)
kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000, add.units = TRUE, cex=0.8, tick.len = 3) #add base numbers
kpAddMainTitle(kp, "Epigenetic Regulation in K562", cex=2) #add title
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.1, gene.name.cex = 1.5) #plot genes (r0/1 position genes at the bottom)
kpPlotRegions(kp, K562.hmm, col=K562.hmm$V9, r0=0.15, r1=0.18) #plot context (weak/strong enhancers, promoters, etc) - color at column V9
kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0=0.15, r1=0.18, cex=1.2) #add label to context regions

#Set total number of tracks to display
total.tracks = length(histone.marks)+length(DNA.binding)

#Adjust histone settings
out.at = autotrack(1:length(histone.marks), total.tracks, margin = 0.3, r0=0.23) #automatically set r0/r1 values (vertical position in plot)
kpAddLabels(kp, labels = "Histone marks", r0 = out.at$r0, r1=out.at$r1, cex=1.5, srt=90, pos=1, label.margin = 0.14) #add label to histones

#Plot each histone mark data
for(i in seq_len(length(histone.marks))) {

  at = autotrack(i, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1) #set vertical position
  kp = kpPlotBigWig(kp, data=histone.marks[i], ymax="visible.region", r0=at$r0, r1=at$r1, col = "cadetblue2") #plot data
  computed.ymax = ceiling(kp$latest.plot$computed.values$ymax) #automatically adjust Y axis limits
  kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at$r0, r1=at$r1, cex=0.8) #add Y axis
  kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, cex=1.2, label.margin = 0.035) #add histone name

}

#Adjust DNA binding proteins settings
out.at = autotrack((length(histone.marks)+1):total.tracks, total.tracks, margin = 0.3, r0=0.23)
kpAddLabels(kp, labels = "DNA-binding\nproteins", r0 = out.at$r0, r1=out.at$r1, cex=1.5, srt=90, pos=1, label.margin = 0.14)

#Plot each DNA binding protein data
for(i in seq_len(length(DNA.binding))) {

  at = autotrack(i, length(DNA.binding), r0=out.at$r0, r1=out.at$r1, margin = 0.1) #set vertical position
  kp = kpPlotBigWig(kp, data=DNA.binding[i], ymax="visible.region", r0=at$r0, r1=at$r1, col = "darkolivegreen1") #plot data
  computed.ymax = ceiling(kp$latest.plot$computed.values$ymax) #automatically adjust Y axis limits
  kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at$r0, r1=at$r1, cex=0.8) #add Y axis
  kpAddLabels(kp, labels = names(DNA.binding)[i], r0=at$r0, r1=at$r1, cex=1.2, label.margin = 0.035) #add protein name

}
garbage=dev.off()
