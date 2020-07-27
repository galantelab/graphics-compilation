#Code Kaplan-Meier ############################################################

library(survival) #functions survfit and coxph
library(survminer) #function ggsurvplot
library(forestmodel) # function forest_model

#Read file
df = read.delim("data_survival.txt")

#Stratify patients in low/high based on median TPM expression of LCK
df$LCK = ifelse(df$LCK >= median(df$LCK),"high", "low")

#Fit multivariate Cox model with all variables together
fit = coxph(Surv(OS.time, OS) ~ LCK + age.diagnosis + gender + MGMT.status, data = df)

#Make forest plot
pdf("forest_plot.pdf")
pp = suppressWarnings(forest_model(model = fit, exponentiate = T, recalculate_width = T, recalculate_height = T))
print(pp, newpage = FALSE)
dev.off()
