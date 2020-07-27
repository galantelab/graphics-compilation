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
library(tidyverse)
#library(RColorBrewer)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)


# Import arguments.
args <- commandArgs(trailingOnly = TRUE)



############################## MAIN ############################################

# Load dataset from github.
data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/13_AdjacencyDirectedWeighted.csv", header=TRUE)
str(data)


# I need a long format.
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_long) <- c("source", "target", "value")
data_long$target <- paste(data_long$target, " ", sep="")


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())


# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1


# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'


# Make the Network
p <- sankeyNetwork(
		Links = data_long, Nodes = nodes,
		Source = "IDsource", Target = "IDtarget",
		Value = "value", NodeID = "name",
		sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20
		)

# save the widget
library(webshot)
library(htmlwidgets)
saveWidget(p, file = "sankey.html")

#install phantom:
#webshot::install_phantomjs()

# Make a webshot in pdf : high quality but can not choose printed zone
#webshot("sankey.html" , "sankey.pdf", delay = 0.2)


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
#pdf("sankey.pdf", width = 29.7*0.4, height = 21*0.4);

# Define plot area.
#multiplot(g1, rows = 2)

#dev.off();
