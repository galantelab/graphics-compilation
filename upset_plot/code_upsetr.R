#!/usr/bin/R
############################## HEADER ##########################################
#
# upsetr.R: R script to plot an UpSetR plot for inputs of tabular data.
#
#
# Usage:
#
#	$ R --vanilla --slave --args data.tsv [arg2 ...] < upsetr.R
#
# Inputs:
#
#	1. TSV or CSV file with raw data.
#
################################################################################



############################## PREAMBLE ########################################



# Libraries.
#library(tidyverse)
#library(RColorBrewer)
#library(VennDiagram)
#library(pheatmap)
library(UpSetR)


# Import arguments.
args <- commandArgs(trailingOnly = TRUE)



############################## MAIN ############################################

# Import data.
print(paste("Importing files:", args[1], args[2], args[3]))
#rawdata <- read.table(args[1], header = TRUE)
rawdata <- read.table(args[1])



names(rawdata) <- c("ID",
					"CHR", "POS", "POL", "REG", "HOST",
					"ME", "LEN", "DEPTH", "IDENT", "INFO",
					"P_CHR", "P_POS", "P_POL", "FAM", "SUBFAM",
					"TAG", "MULT", "EXCL", "FREQ", "STATS",
					"G1_HM", "G1_HT", "G2_HM", "G2_HT", "G3_HM", "G3_HT",
					"G12", "G13", "G23", "G123",
					"INDIV")

# Remove imprecises and zero counts.
rawdata <- rawdata[!bitwAnd(rawdata$INFO, 1), ]
rawdata <- rawdata[rowSums(rawdata[, c(19, 22:27)]) > 0, ]

rawdata$REG <- as.factor(rawdata$REG)
rawdata$TAG <- as.factor(rawdata$TAG)
str(rawdata)


listInput <- list(
	BLCA = rawdata[rawdata$TAG == "BLCA" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	BLCA_norm = rawdata[rawdata$TAG == "BLCA" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1],
	BRCA = rawdata[rawdata$TAG == "BRCA" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	BRCA_norm = rawdata[rawdata$TAG == "BRCA" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1],
	COAD = rawdata[rawdata$TAG == "COAD" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	COAD_norm = rawdata[rawdata$TAG == "COAD" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1],
	GBM = rawdata[rawdata$TAG == "GBM" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	GBM_norm = rawdata[rawdata$TAG == "GBM" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1],
	LUSC = rawdata[rawdata$TAG == "LUSC" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	LUSC_norm = rawdata[rawdata$TAG == "LUSC" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1],
	LUAD = rawdata[rawdata$TAG == "LUAD" & (rawdata$G1_HM + rawdata$G1_HT) > 1, 1],
	LUAD_norm = rawdata[rawdata$TAG == "LUAD" & (rawdata$G2_HM + rawdata$G2_HT) > 1, 1]
	)



# UpSetR chart.
ups <- upset(fromList(listInput),
	order.by = "freq",
	#sets = c("BLCA", "BRCA", "COAD", "GBM", "LUSC", "LUAD"),
	sets = c("BLCA", "BLCA_norm", "BRCA", "BRCA_norm", "COAD", "COAD_norm", "GBM", "GBM_norm", "LUSC", "LUSC_norm", "LUAD", "LUAD_norm"),
	keep.order = TRUE,
	cutoff = 5,
	#number.angles = 30,
	point.size = 3.5,
	line.size = 2,
	mainbar.y.label = "# Shared Intersections",
	sets.x.label = "Total # of RTC insertions",
	text.scale = c(2, 2, 1.8, 1, 1.8, 1.8)
	)



################################################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	library(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
			ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
	}
	else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
		}
	}
}


################################################################################
# Open device to write.
pdf("upsetr.pdf", width = 29.7*0.4, height = 21*0.4);

# Define plot area.
multiplot(ups, rows = 2)

dev.off();
