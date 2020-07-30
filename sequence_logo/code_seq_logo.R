#Make sequence logos #######################################################

#Load library
library(ggseqlogo)

#Read position frequency matrix
df = read.table("data_logo1.txt", header=T)

#Make sequence logo
pdf("seq_logo1.pdf", width=4, height=2)
ggseqlogo(as.matrix(df))
garbage=dev.off()

#Read list of DNA sequences (e.g., binding sites of transcription factors)
df = read.table("data_logo2.txt")

#Create custom colour scheme
colors = make_col_scheme(chars=c('A','T','C','G'), groups=c('weak bonds', 'weak bonds', 'strong bonds', 'strong bonds'), cols=c('purple', 'purple', 'blue', 'blue'))

#Make sequence logo
pdf("seq_logo2.pdf", width=4, height=5)
p1 = ggseqlogo(as.character(df$V1), method = 'bits', col_scheme=colors)
p2 = ggseqlogo(as.character(df$V1), method = 'prob', col_scheme=colors)
gridExtra::grid.arrange(p1, p2)
garbage=dev.off()
