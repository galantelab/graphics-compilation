#Code Kaplan-Meier ############################################################

library(survival) #function survfit
library(survminer) #function ggsurvplot

#Read file
df = read.delim("data_survival.txt")

#Stratify patients in low/high based on median TPM expression of LCK
df$LCK = ifelse(df$LCK >= median(df$LCK),"high", "low")

#Fit model with patients stratified by LCK
fit = survfit(Surv(OS.time, OS) ~ LCK, data=df)

#Plot Kaplan-Meier curve
pdf("kaplan_LCK.pdf")
pp = ggsurvplot(fit, risk.table=T, pval=T, censor=T, break.time.by=200, 
           legend.title="", font.legend=c(16,"bold"), font.x=c(18, "bold"), 
           font.y=c(18, "bold"), font.tickslab=16, pval.size=6, pval.coord=c(0,0.05), legend = c(0.8, 0.9)) + 
           xlab("Survival time (days)")
print(pp, newpage = FALSE)
dev.off()

#Fit model with patients stratified by LCK and gender
fit = survfit(Surv(OS.time, OS) ~ LCK + gender, data=df)

#Plot Kaplan-Meier curve
pdf("kaplan_LCK_gender.pdf")
pp = ggsurvplot(fit, risk.table=T, pval=T, censor=T, break.time.by=200, 
           legend.title="", font.legend=c(16,"bold"), font.x=c(18, "bold"), 
           font.y=c(18, "bold"), font.tickslab=16, pval.size=6, pval.coord=c(0,0.05), legend = c(0.8, 0.9)) + 
           xlab("Survival time (days)")
print(pp, newpage = FALSE)
dev.off()
