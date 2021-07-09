# ======================================================================
# Mortality analysis (Liver Transplant Recipients [LTR])
#
# These codes build spline plots for LTR (Fig 2/D, 2/E)
#
# NOTE: This is implementation on Mock Data published in this github repo,
#       and is intended for demonstration purposes only, the
#       results are not identical to Figures in the manuscript
#      
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


# Perform survival analysis for LTR
# ============================================================

# 2D Spline Shannon for LTR
# ================================================
# > load and clean input files
CDmeta = read.table("mock_data/MOCK_Meta_Mortality_Cleaned_v2.csv",sep=',',header=T)
CDmeta$ID <- as.character(CDmeta$ID)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

CDmetaSurvival_LTR <- subset(CDmeta, mortality_all_post == "LTR")
summary(CDmetaSurvival_LTR$Deceased)
CDmetaSurvival_LTR$Deceased_numeric
ftable(CDmetaSurvival_LTR$Deceased_numeric)

Diss_RCLR <- read.table("mock_data/Dissimilarities_RCLR.csv",sep=',',header = T)
row.names(Diss_RCLR) <- Diss_RCLR$ID

# LTR
summary(CDmetaSurvival_LTR$Shannon_Species)
d4 <- CDmetaSurvival_LTR[c("Shannon_Species")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
options(datadist = "dd")

f1<- with(CDmetaSurvival_LTR, cph(Surv(Days_since_sampling,Deceased_numeric==1) ~  Shannon_Species))
pdata1 <- Predict(f1, Shannon_Species, ref.zero = TRUE, fun = exp,np = 10000)
range(pdata1$Shannon_Species)
range(CDmetaSurvival_LTR$Shannon_Species)
pdata1 <- Predict(f1, Shannon_Species, ref.zero = F, fun = exp)

#Transform Plot, 1, extract freq per bin, then transform and plot
binFreqTable <- function(x, bins) {
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  return(data.frame(range = ranges, frequency = freq$counts))
  
}
x <- CDmetaSurvival_LTR$Shannon_Species
FreqDF <- binFreqTable(x,seq(0,3.7,by=0.1))
FreqDF$freqtransf <- FreqDF$frequency/10
summary(FreqDF$rangeNum)
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/10

maxY_LTR_Shannon <- 6
maxYforSpline_LTR_Shannon <- max(pdata1$yhat[pdata1$yhat < maxY_LTR_Shannon])

# > generate plot 2D
Spline_Shannon_LTR <- ggplot(CDmetaSurvival_LTR, aes(Shannon_Species))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline_LTR_Shannon)), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  geom_line(aes(y = yhat), col = "#FF7D00",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Shannon Diversity Index") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  #range(pdata1$Shannon_Species)
  scale_x_continuous(limits=c(1.0,3.6),expand = c(0,0), breaks=c(1.5, 2.0, 2.5, 3.0, 3.5))+
  scale_y_continuous(limits=c(0,maxY_LTR_Shannon),expand = c(0,0),
                     sec.axis = sec_axis(trans=~., name="Frequency"))+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size=14),
        axis.title.x = element_text(color="black", size=14, face="bold", vjust=-0.5),
        axis.title.y.left = element_text(color="black", size=14, face="bold", vjust=4),
        axis.text.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y.right = element_text(color="black", size=14, face="bold", vjust=4),
        plot.title = element_text(face="bold",size=17,hjust = 0.4,vjust=3),plot.margin=unit(c(1.5,0.5,0.5,0.5),"cm"))

print(Spline_Shannon_LTR)

ggsave(plot = Spline_Shannon_LTR, filename = "mock_data_results/Mockdata_Fig2D_Splines_LTR.png", device = "png", type = "cairo", dpi=300)


# Fig 2E, Spline Distance to HC LTR
# ======================================
#LTR
# prep data
MortalityDissimilarity<- merge(CDmetaSurvival_LTR, Diss_RCLR, by="row.names")
summary(MortalityDissimilarity$DToHC)
# > test uncorrected model
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$Deceased_numeric) 
             ~ scale(MortalityDissimilarity$DToHC))
summary(Surv)

# > build model corrected for Age, Sex and year of transplantation
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$Deceased_numeric) 
             ~ scale(MortalityDissimilarity$DToHC) + 
               MortalityDissimilarity$Months_since_transplantation + 
               MortalityDissimilarity$Tx_year +
               MortalityDissimilarity$Age +
               MortalityDissimilarity$Gender)
summary(Surv)

d4 <- MortalityDissimilarity[c("DToHC")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
#dd <- datadist(d4)
options(datadist = "dd")
#dd$limits$X <- 1.0
f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,Deceased_numeric==1) ~  DToHC))
pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp,np = 10000)

#Transform Plot, 1, extract freq per bin, then transform and plot
binFreqTable <- function(x, bins) {
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  return(data.frame(range = ranges, frequency = freq$counts))
}
x <- MortalityDissimilarity$DToHC
FreqDF <- binFreqTable(x,seq(0,1,by=0.01))
FreqDF$freqtransf <- FreqDF$frequency/10
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/100
maxY_LTR_DTHC <- 6
maxYforSpline_LTR_DTHC<- max(pdata1$yhat[pdata1$yhat < maxY_LTR_DTHC])

# > build spline plot
Spline_DTHC_RCLR_LTR <- ggplot(MortalityDissimilarity, aes(DToHC))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  #geom_ribbon(aes( ymin = lower, ymax = upper), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline_LTR_DTHC)), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  geom_line(aes(y = yhat), col = "#FF7D00",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Dissimilarity to Healthy Control") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  scale_x_continuous(limits=c(0.2,0.45),expand = c(0,0))+
  scale_y_continuous(limits=c(0,maxY_LTR_DTHC),expand = c(0,0),
                     sec.axis = sec_axis(trans=~., name="Frequency", breaks=c(0,1,2,3,4,5,6,7,8)))+
  theme(axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size=14),
        axis.title.x = element_text(color="black", size=14, face="bold", vjust=-0.5),
        axis.title.y.left = element_text(color="black", size=14, face="bold", vjust=4),
        axis.text.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.y.right = element_text(color="black", size=14, face="bold", vjust=4),
        plot.title = element_text(face="bold",size=17,hjust = 0.4,vjust=3),plot.margin=unit(c(1.5,0.5,0.5,0.5),"cm"))

print(Spline_DTHC_RCLR_LTR)
ggsave(plot = Spline_DTHC_RCLR_LTR, filename = "mock_data_results/Mockdata_Fig2E_Splines_LTR.png", device = "png", type = "cairo", dpi=300)

