#Circular bar plot #########################################################

#Load library
library(tidyverse)

#Read file and reorder factor levels
data = read.table("data_circular_barplot.txt", header=T)
data$group = as.factor(data$group)
data$group = factor(data$group,levels(data$group)[c(3,1:2,4)])

#Add 2 empty values (bars) at the end of each group
empty_bar = 2
to_add = data.frame(matrix(NA, empty_bar*nlevels(data$group), ncol(data)))
colnames(to_add) = colnames(data)
to_add$group = rep(levels(data$group), each=empty_bar)
data = rbind(data, to_add)
data = data %>% arrange(group)
data$id = seq(1, nrow(data))
 
#Define position/angle of each gene label
angle = 90-360*(data$id-0.5)/nrow(data)
data$hjust = ifelse(angle < -90,1,0)
data$angle = ifelse(angle < -90,angle+180,angle)

#Define position/angle of each group label
base_data = data %>% group_by(group) %>% summarize(start=min(id), end=max(id) - empty_bar) %>% rowwise() %>% mutate(title=mean(c(start, end)))

#Make the plot
pdf("circular_barplot.pdf", width=5, height=6)
ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +

  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) + ylim(-80,160) + #Create barplot

  theme_minimal() + #Add theme settings
  theme(legend.position = "none", axis.text = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +

  scale_fill_brewer(palette="Set1") + #Set palette of bar colors

  coord_polar() + #Change graph to polar coordinates (circular)

  geom_text(data=data, aes(x=id, y=value+10, label=gene, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= data$angle, inherit.aes = FALSE ) + #Add gene labels

  annotate("text", x=rep(max(data$id),5), y=seq(0,160,by=40), label=as.character(seq(0,160,by=40)), color="grey", size=3 , angle=0, fontface="bold", hjust=1) + #Add y axis

  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), #Add group lines
               colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +

  geom_text(data=base_data, aes(x = title, y = -18, label=group), #Add group labels
            hjust=c(1,1,0,0), colour = "black", alpha=0.8, 
            size=2.5, fontface="bold", inherit.aes = FALSE)
garbage=dev.off()
############################################################################
