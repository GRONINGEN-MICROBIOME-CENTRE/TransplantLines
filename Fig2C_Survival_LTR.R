# ===================================================================
# Mortality analysis (Liver Transplant Recipients [LTR], Figure 2/C)
# ===================================================================
# NOTE: This is implementation on Mock Data published in this github repo,
#       and is intended for demonstration purposes only, the
#       results are not identical to Figures in the manuscript
# ===================================================================
#
# make sure following packages are installed:
#install.packages(c('foreign','tidyverse','ggplot2',
#                   'ggsignif','ggpubr','survival','rms','survminer'))

library(foreign)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(survival)
library(rms)
library(survminer)


# Perform survival analysis for Liver Transplant Recipients
# ============================================================
# load data
CDmeta = read.table("mock_data/MOCK_Meta_Mortality_Cleaned_v2.csv",sep=',',header=T)
CDmeta$ID <- as.character(CDmeta$ID)
#clean whitespaces in IDs
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID
#subset and clean data
CDmetaSurvival_LTR <- subset(CDmeta, Select_samples_mortality == "LTR")
summary(CDmetaSurvival_LTR$Deceased)
CDmetaSurvival_LTR$DeathNumeric[CDmetaSurvival_LTR$Deceased=="Ja"] = 1
CDmetaSurvival_LTR$DeathNumeric[CDmetaSurvival_LTR$Deceased=="Nee"] = 0
print(ftable(CDmetaSurvival_LTR$Deceased_numeric))
# split by median shannon diversity
shannonMedian <- median(CDmetaSurvival_LTR$Shannon_Species)
CDmetaSurvival_LTR$new_diversity_group[CDmetaSurvival_LTR$Shannon_Species<=shannonMedian]="Low diversity"
CDmetaSurvival_LTR$new_diversity_group[CDmetaSurvival_LTR$Shannon_Species>shannonMedian]="High diversity"
print(ftable(CDmetaSurvival_LTR$new_diversity))

#Cox Regression, Unadjusted
# ===============================
Surv<- coxph(formula = Surv(CDmetaSurvival_LTR$Days_since_sampling, CDmetaSurvival_LTR$DeathNumeric) 
             ~ CDmetaSurvival_LTR$new_diversity_group)
summary(Surv)

CDmetaSurvival_LTR$new_diversity_group=as.character(CDmetaSurvival_LTR$new_diversity_group)
CDmetaSurvival_LTR$Gender=as.character(CDmetaSurvival_LTR$Gender)
CDmetaSurvival_LTR$Tx_year=as.character(CDmetaSurvival_LTR$Tx_year)
CDmetaSurvival_LTR$Days_since_sampling=as.numeric(CDmetaSurvival_LTR$Days_since_sampling)
CDmetaSurvival_LTR$Age=as.numeric(CDmetaSurvival_LTR$Age)

#Cox Regression Adjusted for years since tx and age
# =================================================
CDmetaSurvival_LTR$X<- CDmetaSurvival_LTR$new_diversity_group
Surv<- coxph(formula = Surv(CDmetaSurvival_LTR$Days_since_sampling, CDmetaSurvival_LTR$DeathNumeric) 
             ~ CDmetaSurvival_LTR$X + 
               CDmetaSurvival_LTR$Months_since_transplantation + 
               CDmetaSurvival_LTR$Tx_year +
               CDmetaSurvival_LTR$Age +
               CDmetaSurvival_LTR$Gender)
summary(Surv)
a=as.data.frame(summary(Surv)$coefficients)
a$variable=row.names(a)
b=as.data.frame(summary(Surv)$conf.int)
b$variable=row.names(b)
c=merge(b,a)
# save model
write.table(c,file = "mock_data_results/Mockdata_Fig2C_KM_LTR_MortalityModel.csv", sep=',',row.names = F)


# Model for Plot
CDmetaSurvival_LTR$X<- CDmetaSurvival_LTR$new_diversity_group
fit_post <- survfit(Surv(CDmetaSurvival_LTR$Days_since_sampling/30, CDmetaSurvival_LTR$DeathNumeric)
                    ~ CDmetaSurvival_LTR$X, data = CDmetaSurvival_LTR)
print(summary(fit_post))

# Plot
KM_LTR <- ggsurvplot(
  fit_post,                 
  data = CDmetaSurvival_LTR, 
  break.time.by = 6,  xlim = c(0,40), 
  ylim = c(0.7, 1),
  legend.labs=c("High Diversity", "Low Diversity"),
  legend.title= "",
  xlab = "Time in months",
  palette = c("darkgray","#FF7D00"),
  ggtheme = theme_classic(base_size=14, base_family = "Arial"),
  font.family = "Arial",
  break.y.by=0.05)    

# save it
print(KM_LTR)
ggsave(plot = KM_LTR$plot, filename = "mock_data_results/Mockdata_Fig2C_KM_LTR_Months.png", device = "png", type = "cairo", dpi=300)