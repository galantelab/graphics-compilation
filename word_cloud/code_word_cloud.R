#Word cloud ################################################################

#Load library
library(wordcloud)
 
#Read file
df = read.delim("data_word_cloud.txt", header=T)

#Make plot
pdf("wordcloud.pdf",width=6,height=6)
par(bg="black") #Set background to black
wordcloud(df$word, df$value, col=terrain.colors(length(df$word), alpha=0.9), rot.per=0.3)
dev.off()
