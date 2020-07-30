#kpPlotTranscripts ###########################################################

#Load libraries
library(karyoploteR) #functions plotKaryotype, kpAddBaseNumbers and kpPlotTranscripts
library(GenomicFeatures) #function makeTxDbFromGFF

#Read GENCODE V29 as TxDb object
txdb = makeTxDbFromGFF("toy.gencode.v29.annotation.gtf", format = "gtf")

#Extract transcripts, CDS and 5/3' UTR info from TxDb object
transcripts = transcriptsBy(txdb, by = "gene")
exons = exonsBy(txdb, by = "tx", use.names=T)
cds = cdsBy(txdb, by="tx", use.names=T)
five = fiveUTRsByTranscript(txdb, use.names=T)
three = threeUTRsByTranscript(txdb, use.names=T)

#Select all transcripts only from gene ERN2
selected_gene = transcripts[["ENSG00000134398.14"]]

#Select specific 2 transcripts from gene (ERN2-202 and ERN2-205)
selected_transcripts = selected_gene[(elementMetadata(selected_gene)[,2] %in% c("ENST00000457008.6","ENST00000562562.1"))]
names(selected_transcripts) = c("ERN2-202","ERN2-205")

#Select exons, 5/3' UTRs and CDSs from transcripts
t1_exons = exons[["ENST00000457008.6"]]
t1_five = five[["ENST00000457008.6"]]
t1_three = three[["ENST00000457008.6"]]
t1_cds = cds[["ENST00000457008.6"]] #If any
if(is.null(t1_cds)){t1_cds = GRanges()}

t2_exons = exons[["ENST00000562562.1"]]
t2_five = five[["ENST00000562562.1"]]
t2_three = three[["ENST00000562562.1"]]
t2_cds = cds[["ENST00000562562.1"]] # If any
if(is.null(t2_cds)){t2_cds = GRanges()}

#Create non coding regions for transcripts
t1_noncoding = t1_exons
if(!is.null(t1_three)){ t1_noncoding = union(t1_noncoding,t1_three) }
if(!is.null(t1_five)){ t1_noncoding = union(t1_noncoding, t1_five) }

t2_noncoding = t2_exons
if(!is.null(t2_three)){ t2_noncoding = union(t2_noncoding,t2_three) }
if(!is.null(t2_five)){ t2_noncoding = union(t2_noncoding, t2_five) }

#Merge coding regions from two transcripts
t_coding = list(t1_cds,t2_cds)
names(t_coding) = c("ERN2-202","ERN2-205")

#Merge non coding regions from two transcripts
t_noncoding = list(t1_noncoding,t2_noncoding)
names(t_noncoding) = c("ERN2-202","ERN2-205")

#Create data for plot
mylist <- list(transcripts=selected_transcripts, coding.exons=t_coding, non.coding.exons=t_noncoding)

#Define plot region
chr = as.character(seqnames(selected_transcripts)[1])
min_lim = min(c(start(selected_transcripts)[1], start(selected_transcripts)[2]))-7000
max_lim = max(c(end(selected_transcripts)[1], end(selected_transcripts)[2]))+2000
region = toGRanges(chr, min_lim, max_lim )

#Change plot parameters
pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 0
  pp$bottommargin <- 15
  pp$data1inmargin <- 10
  pp$ideogramheight <- 2
  pp$data1outmargin <- 2

#Make plot
pdf("kpPlotTranscripts_ERN2.pdf",height=3,width=5)
karyoplot <- plotKaryotype(zoom=region, plot.params = pp)
kpAddBaseNumbers(karyoplot, tick.dist = 3000, digits=3)
kpPlotTranscripts(karyoplot, data=mylist, y0=c(0, 0.4), y1=c(0.2, 0.6), r0=0, r1=0.3, 
                  transcript.name.position = "left", transcript.name.cex=0.8)
garbage=dev.off()

#Make colored plot
pdf("kpPlotTranscripts_ERN2_colored.pdf",height=3,width=5)
karyoplot <- plotKaryotype(zoom=region, plot.params = pp)
kpAddBaseNumbers(karyoplot, tick.dist = 3000, digits=3)
kpPlotTranscripts(karyoplot, data=mylist, y0=c(0, 0.4), y1=c(0.2, 0.6), r0=0, r1=0.3, 
                  transcript.name.position = "left", transcript.name.cex=0.8,
                  coding.exons.col="cyan4", coding.exons.border.col="cyan4",
                  non.coding.exons.col="coral3", non.coding.exons.border.col="coral3",
                  introns.col="gray48", marks.col="gray48")
garbage=dev.off()

