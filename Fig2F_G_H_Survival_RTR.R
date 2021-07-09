# ======================================================================
# By: Weersma Group, TransplantLines (UMCG)
#
# UMCG Transplantlines Project data analysis
# ======================================================================
#
# Mortality analysis (Renal Transplant Recipients [RTR])
#
# These codes perform survival analysis and build spline plots for RTR 
#  (Fig 2/F, 2/G, 2/H)
#
# NOTE: This is implementation on Mock Data published in this github repo,
#       and is intended for demonstration purposes only, the
#       results are not identical to Figures in the manuscript
#
# NOTE2: Codes are intended to be run from root of this github repo,
#       if running them from different location, make sure to adjust setwd
#       to appropriate location (root of this github repo)
#
# ======================================================================
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


# Survival analysis for RTR,
# Fig 2/F, Kaplan Meier curves for RTR
# ===============================================
setwd(".")

# > load and clean input files, subset RTR
CDmeta = read.table("mock_data/MOCK_Meta_Mortality_Cleaned_v2.csv",sep=',',header=T)
CDmeta$ID <- as.character(CDmeta$ID)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

CDmetaSurvival_RTR <- subset(CDmeta, Select_samples_mortality == "RTR")
summary(CDmetaSurvival_RTR$Deceased)
CDmetaSurvival_RTR$DeathNumeric[CDmetaSurvival_RTR$Deceased=="Ja"] = 1
CDmetaSurvival_RTR$DeathNumeric[CDmetaSurvival_RTR$Deceased=="Nee"] = 0

# > split by median shannon to high & low diversity groups
medianShannon <- median(na.omit(CDmetaSurvival_RTR$Shannon_Species))
CDmetaSurvival_RTR$new_diversity_group[CDmetaSurvival_RTR$Shannon_Species<=medianShannon]="Low diversity"
CDmetaSurvival_RTR$new_diversity_group[CDmetaSurvival_RTR$Shannon_Species>medianShannon]="High diversity"

#Cox Regression, test model without any adjustments (not used for final results)
Surv<- coxph(formula = Surv(CDmetaSurvival_RTR$Days_since_sampling, CDmetaSurvival_RTR$DeathNumeric) 
             ~ CDmetaSurvival_RTR$new_diversity_group)
summary(Surv)

# fix data types for input into coxph() function
CDmetaSurvival_RTR$new_diversity_group=as.character(CDmetaSurvival_RTR$new_diversity_group)
CDmetaSurvival_RTR$Gender=as.character(CDmetaSurvival_RTR$Gender)
CDmetaSurvival_RTR$Tx_year=as.character(CDmetaSurvival_RTR$Tx_year)
CDmetaSurvival_RTR$Days_since_sampling=as.numeric(CDmetaSurvival_RTR$Days_since_sampling)
CDmetaSurvival_RTR$Age=as.numeric(CDmetaSurvival_RTR$Age)

#Cox Regression Adjusted for years since tx and age
CDmetaSurvival_RTR$X<- CDmetaSurvival_RTR$new_diversity_group
Surv<- coxph(formula = Surv(CDmetaSurvival_RTR$Days_since_sampling, CDmetaSurvival_RTR$DeathNumeric) 
             ~ CDmetaSurvival_RTR$X + 
               CDmetaSurvival_RTR$Months_since_transplantation + 
               CDmetaSurvival_RTR$Tx_year +
               CDmetaSurvival_RTR$Age +
               CDmetaSurvival_RTR$Gender)
summary(Surv)

CDmetaSurvival_RTR$X<- CDmetaSurvival_RTR$new_diversity_group
fit_post <- survfit(Surv(CDmetaSurvival_RTR$Days_since_sampling/30, CDmetaSurvival_RTR$DeathNumeric)
                    ~ CDmetaSurvival_RTR$X, data = CDmetaSurvival_RTR)
# Build the plot
KM_RTR <- ggsurvplot(
  fit_post,                 
  data = CDmetaSurvival_RTR, 
  break.time.by = 6,  xlim = c(0,60), 
  ylim = c(0.8, 1),
  legend.labs=c("High Diversity", "Low Diversity"),
  legend.title= "",
  xlab = "Time in months",
  palette = c("darkgray","#45b6fe"),
  ggtheme = theme_classic(base_size=14, base_family = "Arial"),
  font.family = "Arial",
  break.y.by=0.05)    
print(KM_RTR)
ggsave(plot = KM_RTR$plot, filename = "mock_data_results/Mockdata_Fig2F_KM_RTR_Months.png", device = "png", type = "cairo", dpi=300)

# ===================================
# Fig 2/G, Spline plot, Shannon RTR
# ===================================
#summary(CDmetaSurvival_RTR$Shannon_Species)
d4 <- CDmetaSurvival_RTR[c("Shannon_Species")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
options(datadist = "dd")

f1<- with(CDmetaSurvival_RTR, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  Shannon_Species))
pdata1 <- Predict(f1, Shannon_Species, ref.zero = TRUE, fun = exp,np = 10000)
range(pdata1$Shannon_Species)
#range(CDmetaSurvival_RTR$Shannon_Species,na.rm = T)
pdata1 <- Predict(f1, Shannon_Species, ref.zero = F, fun = exp)

#Transform Plot, 1, extract freq per bin, then transform and plot
binFreqTable <- function(x, bins) {
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  return(data.frame(range = ranges, frequency = freq$counts))
}
x <- CDmetaSurvival_RTR$Shannon_Species
FreqDF <- binFreqTable(x,seq(0,3.7,by=0.1))
FreqDF$freqtransf <- FreqDF$frequency/10
summary(FreqDF$rangeNum)
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/10

maxY_Shannon_RTR <- 8
maxYforSpline_Shannon_RTR <- max(pdata1$yhat[pdata1$yhat < maxY_Shannon_RTR])

Spline_Shannon_RTR <- ggplot(CDmetaSurvival_RTR, aes(Shannon_Species))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  geom_ribbon(aes( ymin = lower, ymax = upper), fill="#45b6fe", color="#45b6fe",alpha = .4, data = as.data.frame(pdata1))+
  geom_line(aes(y = yhat), col = "#45b6fe",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Shannon Diversity Index") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  #range(pdata1$Shannon_Species)
  scale_x_continuous(limits=c(1.0,3.6),expand = c(0,0), breaks=c(1.5, 2.0, 2.5, 3.0, 3.5))+
  scale_y_continuous(limits=c(0,maxY_Shannon_RTR),expand = c(0,0), breaks=c(0,1,2,3,4,5,6,7,8),
                     sec.axis = sec_axis(trans=~., name="Frequency", breaks=c(0,1,2,3,4,5,6,7,8)))+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size=14),
        axis.title.x = element_text(color="black", size=14, face="bold", vjust=-0.5),
        axis.title.y.left = element_text(color="black", size=14, face="bold", vjust=4),
        axis.text.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y.right = element_text(color="black", size=14, face="bold", vjust=4),
        plot.title = element_text(face="bold",size=17,hjust = 0.4,vjust=3),plot.margin=unit(c(1.5,0.5,0.5,0.5),"cm"))
print(Spline_Shannon_RTR)
ggsave(plot = Spline_Shannon_RTR, filename = "mock_data_results/Mockdata_Fig2G_Splines_RTR.png", device = "png", type = "cairo", dpi=300)

# ==========================================================================
# Fig 2/H, Spline Plot, Distance to HC RTR
# ==========================================================================
Diss_RCLR <- read.table("mock_data/Dissimilarities_RCLR.csv",sep=',',header = T)
row.names(Diss_RCLR) <- Diss_RCLR$ID

MortalityDissimilarity<- merge(CDmetaSurvival_RTR, Diss_RCLR, by="row.names")

# Build survival model: survival ~ 
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$Deceased_numeric) 
             ~ scale(MortalityDissimilarity$DToHC) + 
               MortalityDissimilarity$Months_since_transplantation + 
               MortalityDissimilarity$Tx_year +
               MortalityDissimilarity$Age +
               MortalityDissimilarity$Gender)
summary(Surv)

summary(MortalityDissimilarity$DToHC)
d4 <- MortalityDissimilarity[c("DToHC")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
options(datadist = "dd")
f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  DToHC))
pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp,np = 10000)

#Transform Plot, 1, extract freq per bin, then transform and plot
binFreqTable <- function(x, bins) {
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  return(data.frame(range = ranges, frequency = freq$counts))
}
x<- MortalityDissimilarity$DToHC
FreqDF <- binFreqTable(x,seq(0,1,by=0.01))
FreqDF$freqtransf <- FreqDF$frequency/10
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/100
maxY_DTHC_RTR  <- 8
maxYforSpline_DTHC_RTR <- max(pdata1$yhat[pdata1$yhat < maxY_DTHC_RTR])

Spline_DTHC_RCLR_RTR <- ggplot(MortalityDissimilarity, aes(DToHC))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  #geom_ribbon(aes( ymin = lower, ymax = upper), fill="#45b6fe", color="#45b6fe",alpha = .4, data = as.data.frame(pdata1))+
  geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline_DTHC_RTR)), fill="#45b6fe", color="#45b6fe",alpha = .4, data = as.data.frame(pdata1))+
  geom_line(aes(y = yhat), col = "#45b6fe",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Dissimilarity to Healthy Control") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  scale_x_continuous(limits=c(0.2,0.4),expand = c(0,0))+
  scale_y_continuous(limits=c(0,maxY_DTHC_RTR),expand = c(0,0),
                     sec.axis = sec_axis(trans=~., name="Frequency", breaks=c(0,1,2,3,4,5,6,7,8)))+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size=14),
        axis.title.x = element_text(color="black", size=14, face="bold", vjust=-0.5),
        axis.title.y.left = element_text(color="black", size=14, face="bold", vjust=4),
        axis.text.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y.right = element_text(color="black", size=14, face="bold", vjust=4),
        plot.title = element_text(face="bold",size=17,hjust = 0.4,vjust=3),plot.margin=unit(c(1.5,0.5,0.5,0.5),"cm"))
print(Spline_DTHC_RCLR_RTR)
ggsave(plot = Spline_DTHC_RCLR_RTR, filename = "mock_data_results/Mockdata_Fig2H_Splines_RTR.png", device = "png", type = "cairo", dpi=300)

