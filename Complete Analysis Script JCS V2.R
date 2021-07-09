#0.0      setwd("J:/Metagenomic sequencing/Nature Microbiology/Code and mock data for Github")
#1.0      Microbiome profiling and description of microbiome----
#2.0      Mortality analysis and Figure 2----
library(foreign)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(survival)
library(rms)
library(xlsx)
library(survminer)
#2A Ordination CS
CDmeta = read.table("mock_data/MOCK_Meta_Mortality_Cleaned_v2.csv",sep=',',header=T)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID
CDmeta_All <- subset(CDmeta, Select_samples_mortality == "LTR" | Select_samples_mortality == "RTR")

#Old
centroids_fig2a <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$Cross_sectional_samplesDag3 %in% c("LTR","RTR","Healthy Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, Cross_sectional_samplesDag3), by="ID") %>%
  mutate(Cross_sectional_samplesDag3=factor(Cross_sectional_samplesDag3, levels=c("LTR","RTR","Healthy Control"))) %>% 
  group_by(Cross_sectional_samplesDag3) %>% 
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2))

Fig2A <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDpwys_All/n2/CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDcard_All/n2/CDcard_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDvfdb_All/n2/CDvfdb_All_sample_loadings.csv", row.names=1) %>% 
  
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$Cross_sectional_samplesDag3 %in% c("LTR","RTR","Healthy Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, Cross_sectional_samplesDag3), "ID") %>%
  mutate(Cross_sectional_samplesDag3=factor(Cross_sectional_samplesDag3, levels=c("LTR","RTR","Healthy Control"))) %>% 
  ggplot(aes(x=PC1, y=PC2, color=Cross_sectional_samplesDag3)) +
  geom_point(shape=19, size=3, alpha=0.4) +
  scale_color_manual(values = c("Healthy Control"="darkgray", "RTR"="#45b6fe", "LTR"="#FF7D00"))+
  theme_classic()+ 
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"))+
  theme(legend.text = element_text(face = "bold", size=12),
        legend.title = element_blank(),
        plot.title=element_text(face = "bold",hjust=0.5),
        legend.position="none")+
  #xlab("PC1: 15.4 %") +
  #ylab("PC2: 11.2 %") +
  #ggtitle("Robust Aitchison PCA (DECOIDE)") +
  geom_point(data=centroids_fig2a, aes(PC1_mean, PC2_mean, fill=Cross_sectional_samplesDag3), shape=21, size=6, color="black") +
  scale_fill_manual(values = c("Healthy Control"="darkgray", "RTR"="#45b6fe", "LTR"="#FF7D00"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.06, -0.03, 0.0,0.03, 0.06),limits=c(-0.09, 0.06))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.06,-0.03,0.0,0.03,0.06), limits=c(-0.06,0.07))+
  stat_ellipse(type="norm", level=0.75)
Fig2A

#Distance to healthy PC1 Taxa
centroids_fig2a <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDpwys_All/n2/CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$Cross_sectional_samplesDag3 %in% c("LTR","RTR","Healthy Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, Cross_sectional_samplesDag3), by="ID") %>%
  mutate(Cross_sectional_samplesDag3=factor(Cross_sectional_samplesDag3, levels=c("LTR","RTR","Healthy Control"))) %>% 
  group_by(Cross_sectional_samplesDag3) %>% 
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2)) #%>% 
#mutate(PC2_mean=0)

# Distance to HC along e.g. PC1
# only comparing it in one plane (i.e. along PC1)
euclidean_dist_centroids_fig2a <- as.matrix(dist(centroids_fig2a[,c("PC1_mean","PC2_mean")], method="euclidean"))
#euclidean_dist_centroids_fig2d <- euclidean_dist_centroids_fig2d/max(euclidean_dist_centroids_fig2d)
rownames(euclidean_dist_centroids_fig2a) <- colnames(euclidean_dist_centroids_fig2a) <- c("LTR","RTR","Healthy Control")

TaxaLoadings <- read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) 
TaxaLoadings <- merge(TaxaLoadings, CDmeta, by="row.names")
TaxaLoadings <- subset(TaxaLoadings, Cross_sectional_samplesDag3=="RTR"|Cross_sectional_samplesDag3=="LTR"|Cross_sectional_samplesDag3=="Healthy Control")
summary(TaxaLoadings$Cross_sectional_samplesDag3)

ggdensity(TaxaLoadings$PC1)
ggqqplot(TaxaLoadings$PC1)
ggdensity(TaxaLoadings$PC2)
ggqqplot(TaxaLoadings$PC2)

Fig2A2 <- data.frame(dist_to_hc=euclidean_dist_centroids_fig2a[3,1:3]) %>% 
  rownames_to_column("group") %>% 
  mutate(group=factor(group, levels=c("LTR","RTR","Healthy Control"))) %>% 
  mutate(dist_to_hc_norm=dist_to_hc/max(dist_to_hc),
         X=1) %>% 
  ggplot(aes(y=dist_to_hc_norm, x=X)) +
  coord_flip() +
  geom_vline(xintercept=1) +
  geom_point(size=8, aes(color=group)) +
  #geom_line(aes(group=1)) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black", angle=45, hjust=1),
    axis.text.y = element_blank())+
  theme(legend.text = element_text(face = "bold", size=12),
        plot.title=element_text(size=14, face="bold", colour = "black", hjust=0.5),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"))+
  labs(x=NULL, y="Mean distance to healthy control") +
  scale_colour_manual(values = c("Healthy Control"="darkgray", "RTR"="#45b6fe", "LTR"="#FF7D00"))#+
#scale_y_continuous(limits=c(0,1),expand=c(0.01,0,0,0.049), breaks=seq(from=0, to=1, by=0.25)) #+
#geom_text(aes(label=group), angle=45, hjust=-0.15, vjust=-1, fontface="bold", size=6)
Fig2A2

#setwd("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 2/V2 Januari 2021")
#ggsave("Fig2A2_mean.png", device = "png", type = "cairo", dpi=300)

#2B Boxplot Diversity LTR, RTR, HC
CDmetaOud <- subset(CDmeta, TypeTxN_TxL_TxD_DAG3 == "TxN" |TypeTxN_TxL_TxD_DAG3 == "TxL" | TypeTxN_TxL_TxD_DAG3 == "DAG3 Control")
CDmetaOud <- subset(CDmetaOud, Cross_sectional_samplesDag3 !="NA")
CDmetaOud$Cross_sectional_samplesDag3 <- droplevels(CDmetaOud$Cross_sectional_samplesDag3)
CDmetaOud$Cross_sectional_samplesDag3 <- recode_factor(CDmetaOud$Cross_sectional_samplesDag3, "Healthy Control" = "HC")
summary(CDmetaOud$Cross_sectional_samplesDag3)

CDmetaOud$Cross_sectional_samplesDag3_v2 <- fct_relevel(CDmetaOud$Cross_sectional_samplesDag3, "HC", "LTR", "RTR")
summary(CDmetaOud$Cross_sectional_samplesDag3_v2)

CDmetaOudRTR <- subset(CDmetaOud, Cross_sectional_samplesDag3=="LTR"|Cross_sectional_samplesDag3=="HC")
wilcox <-wilcox.test(Shannon_Species ~ Cross_sectional_samplesDag3, data=CDmetaOudRTR)
wilcox$p.value

Fig2B <- ggplot(CDmetaOud, aes(x=Cross_sectional_samplesDag3_v2, y=Shannon_Species, fill=Cross_sectional_samplesDag3_v2)) +
  geom_jitter(aes(color=Cross_sectional_samplesDag3_v2), alpha=0.0, width=0.1)+
  geom_boxplot(aes(color=Cross_sectional_samplesDag3_v2), alpha=0, width=0.3, size=0.8)+
  theme_classic()+
  theme(
    axis.title.x = element_blank(),    
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"),
    plot.title=element_blank())+
  geom_signif(comparisons = list(c("HC", "RTR"),
                                 c("HC", "LTR"),
                                 c("RTR", "LTR")),
              textsize=4,
              step_increase = 0.09,
              annotations = c("6.3 x 10-21",
                              "3.1 x 10-9",
                              "0.3"))+
  geom_violin(aes(color = Cross_sectional_samplesDag3_v2), trim = FALSE, position = position_dodge(0.9), alpha=0, size=0.8)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")+
  scale_color_manual(values = c("HC"="darkgray", "RTR"="#45b6fe", "LTR"="#FF7D00"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1, decimal.mark = '.'), breaks=c(0,1,2,3,4))+
  ylab("Shannon Diversity Index")
Fig2B

#Combina 2A and 2B
Fig2AB <- ggarrange(Fig2A, Fig2B, nrow=1, align=c("h"), labels = c("A", "B"))
#ggsave("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 2/V2 Januari 2021/Fig2AB__Combined.png", device = "png", type = "cairo", dpi=300, w=12)
Fig2AB2 <- ggarrange(Fig2A2, Fig2A2, nrow=1, align=c("h"))
Fig2ABCombined <- ggarrange(Fig2AB, Fig2AB2, ncol=1, align=c("v"))
Fig2ABCombined
#setwd("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 2/V2 Januari 2021")
ggsave("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 2/V2 Januari 2021/Fig2AB_HC_axis_Combined_MEAN.png", device = "png", type = "cairo", dpi=300,
       h=12, w=12)

#2C Kaplan LTR
setwd("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Mortality/Final For manuscript 22-9-2020")
CDmeta = read.spss("./TxN_metadata_Oud_New_Donor_DAG3_Cleaned_survival.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

CDmetaSurvival_LTR <- subset(CDmeta, mortality_no_duplicate_1year == "LTR")
summary(CDmetaSurvival_LTR$RECOVERL)
CDmetaSurvival_LTR$DeathNumeric[CDmetaSurvival_LTR$RECOVERL=="Ja"] = 1
CDmetaSurvival_LTR$DeathNumeric[CDmetaSurvival_LTR$RECOVERL=="Nee"] = 0
ftable(CDmetaSurvival_LTR$DeathNumeric)

median(CDmetaSurvival_LTR$Shannon_Species)
CDmetaSurvival_LTR$new_diversity_group[CDmetaSurvival_LTR$Shannon_Species<=2.481635]="Low diversity"
CDmetaSurvival_LTR$new_diversity_group[CDmetaSurvival_LTR$Shannon_Species>2.481635]="High diversity"
ftable(CDmetaSurvival_LTR$new_diversity)

#Cox Regression Unadjusted
Surv<- coxph(formula = Surv(CDmetaSurvival_LTR$Days_since_sampling, CDmetaSurvival_LTR$DeathNumeric) 
             ~ CDmetaSurvival_LTR$new_diversity_group)
summary(Surv)

CDmetaSurvival_LTR$new_diversity_group=as.character(CDmetaSurvival_LTR$new_diversity_group)
CDmetaSurvival_LTR$Gender=as.character(CDmetaSurvival_LTR$Gender)
CDmetaSurvival_LTR$Tx_year=as.character(CDmetaSurvival_LTR$Tx_year)
CDmetaSurvival_LTR$Days_since_sampling=as.numeric(CDmetaSurvival_LTR$Days_since_sampling)
CDmetaSurvival_LTR$Age=as.numeric(CDmetaSurvival_LTR$Age)

#Cox Regression Adjusted for years since tx and age
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
#write.csv2(c,file = "LTR.csv", quote = F, sep = ";")

CDmetaSurvival_LTR$X<- CDmetaSurvival_LTR$new_diversity_group
fit_post <- survfit(Surv(CDmetaSurvival_LTR$Months_since_sampling, CDmetaSurvival_LTR$DeathNumeric)
                    ~ CDmetaSurvival_LTR$X, data = CDmetaSurvival_LTR)

# Model for Plot
# ggsurvplot(fit_post, data = CDmetaSurvival_LTR, censor.shape="|", censor.size = 4, pval = TRUE)
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
KM_LTR
KM_LTR<- ggpar(KM_LTR, 
               font.main = c(14, "bold"),
               font.x = c(14, "bold"),
               font.y = c(14, "bold"),
               font.caption = c(14, "bold"), 
               font.legend = c(14, "bold"),
               font.ytickslab = c(14, "bold", "black"),
               font.xtickslab = c(14, "bold", "black"))
KM_LTR
#ggsave(plot = KM_LTR$plot, "KM_LTR_Months.png", device = "png", type = "cairo", dpi=300)
#ggsave("Survival_LTR.png", device = "png", type = "cairo", dpi=300, w=14)

#2D Spline Shannon LTR
setwd("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Mortality/Final For manuscript 22-9-2020")
CDmeta = read.spss("./TxN_metadata_Oud_New_Donor_DAG3_Cleaned_survival.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

CDmetaSurvival_RTR <- subset(CDmeta, mortality_all_post== "RTR")
summary(CDmetaSurvival_RTR$RECOVERL)
CDmetaSurvival_RTR$DeathNumeric
ftable(CDmetaSurvival_RTR$DeathNumeric)
CDmetaSurvival_RTR <- subset(CDmetaSurvival_RTR, Shannon_Species !="NA")

CDmetaSurvival_LTR <- subset(CDmeta, mortality_all_post == "LTR")
summary(CDmetaSurvival_LTR$RECOVERL)
CDmetaSurvival_LTR$DeathNumeric
ftable(CDmetaSurvival_LTR$DeathNumeric)

Diss_RCLR <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Mortality\\Final For manuscript 22-9-2020\\dist_to_hc_by_sampleID_rclr.xlsx",1)
row.names(Diss_RCLR) <- Diss_RCLR$ID

setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 6")
scale_color_manual(values = c("Healthy Control"="darkgray", "RTR"="#45b6fe", "LTR"="#FF7D00"))

#LTR
summary(CDmetaSurvival_LTR$Shannon_Species)
d4 <- CDmetaSurvival_LTR[c("Shannon_Species")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
options(datadist = "dd")

f1<- with(CDmetaSurvival_LTR, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  Shannon_Species))
# RG: Note: increases np in Predict function produces more data points which helps
#     with lil issues with ggplot like tiny cutoff on the top of the plot
#     here it is set to 10k which si extreme, but makes final plot very smooth,
#     1000 is OK for most uses
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
x<- CDmetaSurvival_LTR$Shannon_Species
FreqDF <- binFreqTable(x,seq(0,3.3,by=0.1))
FreqDF$freqtransf <- FreqDF$frequency/10
summary(FreqDF$rangeNum)
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/10
ggplot(FreqDF)+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity")+
  scale_x_continuous()

maxY_LTR_Shannon <- 6
maxYforSpline_LTR_Shannon <- max(pdata1$yhat[pdata1$yhat < maxY_LTR_Shannon])

Spline_Shannon_LTR <- ggplot(CDmetaSurvival_LTR, aes(Shannon_Species))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  # RG: here we introduce limit for ymin for ribbon to force ggplot to render it, otherwise it skips part of spline which goes above (and below) y limits
  geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline_LTR_Shannon)), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  # next line is old code for comparison
  #geom_ribbon(aes( ymin = lower, ymax = upper), fill="#FF7D00", color="#FF7D00",alpha = .4, data = pdataDF)+
  geom_line(aes(y = yhat), col = "#FF7D00",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Shannon Diversity Index") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  #range(pdata1$Shannon_Species)
  scale_x_continuous(limits=c(1.5,3.16),expand = c(0,0), breaks=c(1.5, 2.0, 2.5, 3.0, 3.5))+
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
Spline_Shannon_LTR

#2E Spline Distance to HC LTR
#LTR
#Scale
MortalityDissimilarity<- merge(CDmetaSurvival_LTR, Diss_RCLR, by="row.names")
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$DeathNumeric) 
             ~ scale(MortalityDissimilarity$DToHC))
summary(Surv)

#So if the highter the dissimilarity between LTR and HC the higher the mortality !
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$DeathNumeric) 
             ~ scale(MortalityDissimilarity$DToHC) + 
               MortalityDissimilarity$Months_since_transplantation + 
               MortalityDissimilarity$Tx_year +
               MortalityDissimilarity$Age +
               MortalityDissimilarity$Gender)
summary(Surv)

summary(MortalityDissimilarity$DToHC)
d4 <- MortalityDissimilarity[c("DToHC")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
#dd <- datadist(d4)
options(datadist = "dd")
#dd$limits$X <- 1.0
f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  DToHC))
#f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  DToHC))
pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp,np = 10000)
#pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp)

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
ggplot(MortalityDissimilarity)+
  geom_bar(aes(x=DToHC), stat="bin", binwidth=0.01)

ggplot(FreqDF)+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity")+
  scale_x_continuous()

maxY_LTR_DTHC <- 6
maxYforSpline_LTR_DTHC<- max(pdata1$yhat[pdata1$yhat < maxY_LTR_DTHC])

Spline_DTHC_RCLR_LTR <- ggplot(MortalityDissimilarity, aes(DToHC))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  #geom_ribbon(aes( ymin = lower, ymax = upper), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline_LTR_DTHC)), fill="#FF7D00", color="#FF7D00",alpha = .4, data = as.data.frame(pdata1))+
  geom_line(aes(y = yhat), col = "#FF7D00",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Dissimilarity to Healthy Control") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  scale_x_continuous(limits=c(0.2,0.4),expand = c(0,0))+
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
Spline_DTHC_RCLR_LTR
#,breaks = c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)
#,breaks=c(0,1,2,3,4,5,6,7,8)

#2F Kaplan RTR
setwd("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Mortality/Final For manuscript 22-9-2020")
CDmeta = read.spss("./TxN_metadata_Oud_New_Donor_DAG3_Cleaned_survival.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

CDmetaSurvival_RTR <- subset(CDmeta, mortality_no_duplicate_1year == "RTR")
summary(CDmetaSurvival_RTR$RECOVERL)
CDmetaSurvival_RTR$DeathNumeric[CDmetaSurvival_RTR$RECOVERL=="Ja"] = 1
CDmetaSurvival_RTR$DeathNumeric[CDmetaSurvival_RTR$RECOVERL=="Nee"] = 0
ftable(CDmetaSurvival_RTR$DeathNumeric)

median(na.omit(CDmetaSurvival_RTR$Shannon_Species))
CDmetaSurvival_RTR$new_diversity_group[CDmetaSurvival_RTR$Shannon_Species<=2.479713]="Low diversity"
CDmetaSurvival_RTR$new_diversity_group[CDmetaSurvival_RTR$Shannon_Species>2.479713]="High diversity"
ftable(CDmetaSurvival_RTR$new_diversity)

#Cox Regression Unadjusted
Surv<- coxph(formula = Surv(CDmetaSurvival_RTR$Days_since_sampling, CDmetaSurvival_RTR$DeathNumeric) 
             ~ CDmetaSurvival_RTR$new_diversity_group)
summary(Surv)

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
fit_post <- survfit(Surv(CDmetaSurvival_RTR$Months_since_sampling, CDmetaSurvival_RTR$DeathNumeric)
                    ~ CDmetaSurvival_RTR$X, data = CDmetaSurvival_RTR)

xlim = c(0,1825), # present narrower X axis, but not affect
ylim = c(0.6, 1),
# Model for Plot
# ggsurvplot(fit_post, data = CDmetaSurvival_LTR, censor.shape="|", censor.size = 4, pval = TRUE)
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
KM_RTR
KM_RTR<- ggpar(KM_RTR, 
               font.main = c(14, "bold"),
               font.x = c(14, "bold"),
               font.y = c(14, "bold"),
               font.caption = c(14, "bold"), 
               font.legend = c(14, "bold"),
               font.ytickslab = c(14, "bold", "black"),
               font.xtickslab = c(14, "bold", "black"))
KM_RTR
ggsave(plot = KM_RTR$plot, "KM_RTR_monhts.png", device = "png", type = "cairo", dpi=300)

#ggarrange(KM_RTR$plot, SplineRTR, SplineRTRDis, nrow=1, 
labels = c("A", "B", "C"))

#ggsave("Survival_RTR.png", device = "png", type = "cairo", dpi=300, w=14)

#2G Spline Shannon RTR
#RTR
summary(CDmetaSurvival_RTR$Shannon_Species)
d4 <- CDmetaSurvival_RTR[c("Shannon_Species")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
options(datadist = "dd")

f1<- with(CDmetaSurvival_RTR, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  Shannon_Species))
# RG: Note: increases np in Predict function produces more data points which helps
#     with lil issues with ggplot like tiny cutoff on the top of the plot
#     here it is set to 10k which si extreme, but makes final plot very smooth,
#     1000 is OK for most uses
pdata1 <- Predict(f1, Shannon_Species, ref.zero = TRUE, fun = exp,np = 10000)
range(pdata1$Shannon_Species)
range(CDmetaSurvival_RTR$Shannon_Species)
pdata1 <- Predict(f1, Shannon_Species, ref.zero = F, fun = exp)

#Transform Plot, 1, extract freq per bin, then transform and plot
binFreqTable <- function(x, bins) {
  
  freq = hist(x, breaks=bins, include.lowest=TRUE, plot=FALSE)
  
  ranges = paste(head(freq$breaks,-1), freq$breaks[-1], sep=" - ")
  
  return(data.frame(range = ranges, frequency = freq$counts))
  
}
x<- CDmetaSurvival_RTR$Shannon_Species
FreqDF <- binFreqTable(x,seq(0,3.3,by=0.1))
FreqDF$freqtransf <- FreqDF$frequency/10
summary(FreqDF$rangeNum)
FreqDF$rangeNum <- as.numeric(FreqDF$range)
FreqDF$rangeNum <- FreqDF$rangeNum/10
ggplot(FreqDF)+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity")+
  scale_x_continuous()

maxY_Shannon_RTR <- 8
maxYforSpline_Shannon_RTR <- max(pdata1$yhat[pdata1$yhat < maxY_Shannon_RTR])

Spline_Shannon_RTR <- ggplot(CDmetaSurvival_RTR, aes(Shannon_Species))+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity", data=FreqDF, color="black", fill="white")+
  # RG: here we introduce limit for ymin for ribbon to force ggplot to render it, otherwise it skips part of spline which goes above (and below) y limits
  #geom_ribbon(aes( ymin = pmax(lower,0.00001), ymax = pmin(upper,maxYforSpline)), fill="#45b6fe", color="#45b6fe",alpha = .4, data = as.data.frame(pdata1))+
  geom_ribbon(aes( ymin = lower, ymax = upper), fill="#45b6fe", color="#45b6fe",alpha = .4, data = as.data.frame(pdata1))+
  # next line is old code for comparison
  #geom_ribbon(aes( ymin = lower, ymax = upper), fill="#FF7D00", color="#FF7D00",alpha = .4, data = pdataDF)+
  geom_line(aes(y = yhat), col = "#45b6fe",alpha = 0.8,linetype=1,size=1.3, data = as.data.frame(pdata1))+
  geom_hline(aes(yintercept = 1), linetype = 2,color="grey",size=1.2,alpha=1)+
  xlab("Shannon Diversity Index") + ylab("Hazard ratio Mortality (95% CI)")+
  theme_classic()+
  #range(pdata1$Shannon_Species)
  scale_x_continuous(limits=c(1.5,3.16),expand = c(0,0), breaks=c(1.5, 2.0, 2.5, 3.0, 3.5))+
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
Spline_Shannon_RTR
#2H Spline Distance to HC RTR
#RTR
#Scale
MortalityDissimilarity<- merge(CDmetaSurvival_RTR, Diss_RCLR, by="row.names")
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$DeathNumeric) 
             ~ scale(MortalityDissimilarity$DToHC))
summary(Surv)

#So if the highter the dissimilarity between LTR and HC the higher the mortality !
Surv<- coxph(formula = Surv(MortalityDissimilarity$Days_since_sampling, MortalityDissimilarity$DeathNumeric) 
             ~ scale(MortalityDissimilarity$DToHC) + 
               MortalityDissimilarity$Months_since_transplantation + 
               MortalityDissimilarity$Tx_year +
               MortalityDissimilarity$Age +
               MortalityDissimilarity$Gender)
summary(Surv)

summary(MortalityDissimilarity$DToHC)
d4 <- MortalityDissimilarity[c("DToHC")]
dd <- datadist(d4, q.display=c(0,1),q.effect=c(0,1))
#dd <- datadist(d4)
options(datadist = "dd")
#dd$limits$X <- 1.0
f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  DToHC))
#f1<- with(MortalityDissimilarity, cph(Surv(Days_since_sampling,DeathNumeric==1) ~  DToHC))
pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp,np = 10000)
#pdata1 <- Predict(f1, DToHC, ref.zero = TRUE, fun = exp)

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
ggplot(MortalityDissimilarity)+
  geom_bar(aes(x=DToHC), stat="bin", binwidth=0.01)

ggplot(FreqDF)+
  geom_bar(aes(x=rangeNum, y=freqtransf), stat="identity")+
  scale_x_continuous()

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
Spline_DTHC_RCLR_RTR
#,breaks = c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)
#,breaks=c(0,1,2,3,4,5,6,7,8)
#Arrange All
Fig2AB
Fig2CDE <- ggarrange(KM_LTR$plot, Spline_Shannon_LTR, Spline_DTHC_RCLR_LTR, nrow=1, align=c("h"), labels=c("C","D","E"))
Fig2FGB <- ggarrange(KM_RTR$plot, Spline_Shannon_RTR, Spline_DTHC_RCLR_RTR, nrow=1, align=c("h"), labels=c("F","G","H"))
FinalFig2 <- ggarrange(Fig2ABCombined, Fig2CDE, Fig2FGB, ncol=1, align=c("v"))
FinalFig2
ggsave(plot=FinalFig2, "J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 2/V2 Januari 2021/Figure2_Final_Months_V2.png", device = "png", type = "cairo",
       dpi=300, w=14, h=14)
#3.0      Differential Abundance analysis----
library(foreign)
library(stringr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(scales)
library(gsubfn)
library(xlsx)
library(ggcorrplot)

#3.1      Load Data ----
CDtaxa <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_metaphlan_cleaned_filtered.txt",sep='\t',header=T, fill=T)
CDtaxa$ID <- as.character(CDtaxa$ID)
trimws(CDtaxa$ID)
rownames(CDtaxa) <- CDtaxa$ID

CDpwys <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_humann_pathabundance_filtered.txt",sep='\t',header=T)
CDpwys$ID <- as.character(CDpwys$ID)
trimws(CDpwys$ID)
rownames(CDpwys) <- CDpwys$ID

CDcard <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_CARD_filtered.txt",sep='\t',header=T)
CDcard$ID <- as.character(CDcard$ID)
trimws(CDcard$ID)
CDcard <- CDcard[!duplicated(CDcard$ID),]
rownames(CDcard) <- CDcard$ID 

CDvfdb<- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_VFDB_filtered.txt",sep='\t',header=T)
CDvfdb$ID <- as.character(CDvfdb$ID)
trimws(CDvfdb$ID)
CDvfdb <- CDvfdb[!duplicated(CDvfdb$ID),]
rownames(CDvfdb) <- CDvfdb$ID 

CDmeta = read.spss("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\Metadata_V3_January2021.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

summary(CDmeta$Selection_Long_Cross_S)
CDmeta_All <- subset(CDmeta, Selection_Long_Cross_S == "Yes")

#Taxa
CDtaxa_All <- CDmeta_All["ID"]
CDtaxa_All <- merge(CDtaxa_All, CDtaxa, by="ID")
row.names(CDtaxa_All) <- CDtaxa_All$ID
CDtaxa_All$ID <- NULL
#Pathways
CDpwys_All <- CDmeta_All["ID"]
CDpwys_All <- merge(CDpwys_All, CDpwys, by="ID")
row.names(CDpwys_All) <- CDpwys_All$ID
CDpwys_All$ID <- NULL
#Antibiotic Res. Genes
CDcard_All <- CDmeta_All["ID"]
CDcard_All <- merge(CDcard_All, CDcard, by="ID")
row.names(CDcard_All) <- CDcard_All$ID
CDcard_All$ID <- NULL
CDcard_All <- CDcard_All[!is.na(rowSums(CDcard_All)),]
#Virulence Factors
CDvfdb_All <- CDmeta_All["ID"]
CDvfdb_All <- merge(CDvfdb_All, CDvfdb, by="ID")
row.names(CDvfdb_All) <- CDvfdb_All$ID
CDvfdb_All$ID <- NULL

# this is the RCLR transformed data from DECOIDE
# species
ftbl_rclr_species <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\CDtaxa_All_rclr.csv", row.names=1)
CDtaxa_All_only_s <-  CDtaxa_All[,grep(x=colnames(CDtaxa_All),pattern="s__")]
row.names(CDtaxa_All_only_s)
rownames(ftbl_rclr_species) <- str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7)[,7] #colnames(CDtaxa_All_only_s)
colnames(ftbl_rclr_species) <- rownames(CDtaxa_All_only_s)

ftbl_rclr_species_all_levels <- read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\CDtaxa_All_rclr_all_taxanomic_levels.csv", row.names=1)
ftbl_rclr_species_all_levels <- t.data.frame(ftbl_rclr_species_all_levels)
ftbl_rclr_species_all_levels <- as.data.frame(ftbl_rclr_species_all_levels)

colnames(ftbl_rclr_species_all_levels) <- colnames(CDtaxa_All)  
rownames(ftbl_rclr_species_all_levels) <- rownames(CDtaxa_All)

# pathways
ftbl_rclr_pwys <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\CDpwys_All_rclr.csv", row.names=1)
rownames(ftbl_rclr_pwys) <- colnames(CDpwys_All)
colnames(ftbl_rclr_pwys) <- rownames(CDpwys_All)

#card
ftbl_rclr_card <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\CDcard_All_rclr.csv", row.names=1)
rownames(ftbl_rclr_card) <- colnames(CDcard_All)
colnames(ftbl_rclr_card) <- rownames(CDcard_All)

#CDvfdb
ftbl_rclr_vfdb <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\CDvfdb_All_rclr.csv", row.names=1)
rownames(ftbl_rclr_vfdb) <- colnames(CDvfdb_All)
colnames(ftbl_rclr_vfdb) <- rownames(CDvfdb_All)

#3.2      DA LM for Taxa and Pathways LTR vs HC----
#This is for cross-sectional (ESLD and ESRD) analysis plus covariates
meta_liver <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","LTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)

summary(meta_liver$Group)
sum(is.na(meta_liver$Age))
sum(is.na(meta_liver$BMI))
sum(is.na(meta_liver$Gender))
sum(is.na(meta_liver$Smoking))
summary(meta_liver$Smoking)
sum(is.na(meta_liver$PPI))
sum(is.na(meta_liver$Laxatives))
sum(is.na(meta_liver$Antibiotics))

with(meta_liver, tapply(Smoking, Group, summary))
with(meta_liver, tapply(PPI, Group, summary))
with(meta_liver, tapply(Laxatives, Group, summary))
with(meta_liver, tapply(Antibiotics, Group, summary))

wilcox.test(Age ~ Group, data=meta_liver)
wilcox.test(BMI ~ Group, data=meta_liver)
with(meta_liver, chisq.test(Gender, Group))
with(meta_liver, chisq.test(Smoking, Group))
with(meta_liver, chisq.test(PPI, Group))
with(meta_liver, chisq.test(Laxatives, Group))
with(meta_liver, chisq.test(Antibiotics, Group))

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_species[,colnames(ftbl_rclr_species) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_df <- data.frame(ftbl_rclr)
# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")
tapply(lm_data$s__Faecalibacterium_prausnitzii, lm_data$Group, summary)

lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- CDtaxa_All %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")
tapply(lm_data$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii, lm_data$Group, summary)


lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))#[,c(2,5,7)]
#colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")
colnames(species_taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
species_taxonomy$taxon_id=species_taxonomy$Species
Healthy <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\Healthy_merged.xlsx",1) 
Healthy$Taxonomy <- NULL
dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) 
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat <- left_join(dat, Healthy, by="taxon_id")
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  mutate(Health_Overlap_Gacesa=ifelse(Direction==Unhealthy_Direction_Gacesa,"Yes","No"))%>%
  mutate(Health_Overlap_Gupta=ifelse(Direction==Unhealthy_Direction_Gupta,"Yes","No")) %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`,
         `Butyrare producer` = Butyrate_Prod)
dat <- left_join(dat, species_taxonomy, by="taxon_id")
#view(dat)  

#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\Results_DA_LTR_HC_Full_Model_Species_V2.xlsx")

#Pathways
meta_liver <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","LTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_pwys[,colnames(ftbl_rclr_pwys) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
pathway_classes <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\path_annotated_index_V3.txt", sep='\t',header=T, fill=T)
pathway_classes$taxon_id <- pathway_classes$MetaCyc_complete_name
dat <- left_join(dat, pathway_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_LTR_HC_Full_Model_Pathways.xlsx")


#3.3      DA LTR vs HC CARD VFDB PA----
#CARD
meta_liver <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","LTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)
summary(meta_liver$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_card[,colnames(ftbl_rclr_card) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
#ftbl_rclr_df <- data.frame(ftbl_rclr)
#sum(is.na(ftbl_rclr_df$gb.ACA23185_1.ARO_3000194.tetW__Bifidobacterium_longum_)
#CDcard_All <- data.frame(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0

ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]

glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

#tapply(glm_data$gb.ABA71729_1.ARO_3002959.vanYG1__Enterococcus_faecalis_, glm_data$Group, sum)

#taxon="gb.CAE51745_1.ARO_3000168.tet.D.__Serratia_marcescens_"

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

#lm_ls$gb.AAM09851_1.ARO_3002923.vanRD__Enterococcus_faecium_

#m <- lm(glm_data[,taxon] ~ Group + Age, data=glm_data[,c(taxon,"Group", "Age")])

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
CARD_classes <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\CARD_Annotation.xlsx", 1)
CARD_classes$taxon_id <- CARD_classes$Complete.taxonomy
dat <- left_join(dat, CARD_classes, by="taxon_id")

#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
view(dat)

#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_CARD_Unfiltered.xlsx")

#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_LTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_CARD_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_LTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_CARD_Filtered1%.xlsx")

#VFDB
meta_liver <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","LTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)
summary(meta_liver$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_vfdb[,colnames(ftbl_rclr_vfdb) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]
# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

tapply(glm_data$gb.ABA71729_1.ARO_3002959.vanYG1__Enterococcus_faecalis_, glm_data$Group, sum)

#taxon="gb.CAE51745_1.ARO_3000168.tet.D.__Serratia_marcescens_"

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat$VFID <- dat$taxon_id
dat$VFID <- str_split_fixed(str_split_fixed(dat$VFID, pattern="\\..__", 2)[,2], pattern="\\.", 2)[,1]
VFDB_Anno <- read.csv(file = "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\VFDB_annotation.csv", 
                      stringsAsFactors = FALSE)
dat <- left_join(dat, VFDB_Anno, by="VFID")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_VFDB_Unfiltered.xlsx")

#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_LTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_VFDB_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_LTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_LTR_HC_Full_Model_PA_RCLR_VFDB_Filtered1%.xlsx")


#3.4      DA RTR LM for Taxa and Pathways vs HC----
#This is for cross-sectional (ESLD and ESRD) analysis plus covariates
meta_renal <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","RTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)


summary(meta_renal$Group)
sum(is.na(meta_renal$Age))
sum(is.na(meta_renal$BMI))
sum(is.na(meta_renal$Gender))
sum(is.na(meta_renal$Smoking))
summary(meta_renal$Smoking)
sum(is.na(meta_renal$PPI))
sum(is.na(meta_renal$Laxatives))
sum(is.na(meta_renal$Antibiotics))

with(meta_renal, tapply(Smoking, Group, summary))
with(meta_renal, tapply(PPI, Group, summary))
with(meta_renal, tapply(Laxatives, Group, summary))
with(meta_renal, tapply(Antibiotics, Group, summary))

wilcox.test(Age ~ Group, data=meta_renal)
wilcox.test(BMI ~ Group, data=meta_renal)
with(meta_renal, chisq.test(Gender, Group))
with(meta_renal, chisq.test(Smoking, Group))
with(meta_renal, chisq.test(PPI, Group))
with(meta_renal, chisq.test(Laxatives, Group))
with(meta_renal, chisq.test(Antibiotics, Group))

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_species[,colnames(ftbl_rclr_species) %in% meta_renal$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_df <- data.frame(ftbl_rclr)
# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))#[,c(2,5,7)]
#colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")
colnames(species_taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
species_taxonomy$taxon_id=species_taxonomy$Species
Healthy <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\Healthy_merged.xlsx",1) 
Healthy$Taxonomy <- NULL
dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) 
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat <- left_join(dat, Healthy, by="taxon_id")
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  mutate(Health_Overlap_Gacesa=ifelse(Direction==Unhealthy_Direction_Gacesa,"Yes","No"))%>%
  mutate(Health_Overlap_Gupta=ifelse(Direction==Unhealthy_Direction_Gupta,"Yes","No")) %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`,
         `Butyrare producer` = Butyrate_Prod)
dat <- left_join(dat, species_taxonomy, by="taxon_id")
#view(dat)  
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\Results_DA_RTR_HC_Full_Model_Species_V2.xlsx")

#Pathways
meta_renal <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","RTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)

meta_renal[is.na(meta_renal$Smoking),]$Smoking <- "Nee"


# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_pwys[,colnames(ftbl_rclr_pwys) %in% meta_renal$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
pathway_classes <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\path_annotated_index_V3.txt", sep='\t',header=T, fill=T)
pathway_classes$taxon_id <- pathway_classes$MetaCyc_complete_name
dat <- left_join(dat, pathway_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_RTR_HC_Full_Model_Pathways.xlsx")

#3.5      DA RTR vs HC CARD VFDB PA----
#CARD
meta_renal <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","RTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)
summary(meta_renal$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_card[,colnames(ftbl_rclr_card) %in% meta_renal$ID] # change between meta_renal and meta_renal
ftbl_rclr <- t(ftbl_rclr)
#ftbl_rclr_df <- data.frame(ftbl_rclr)
#sum(is.na(ftbl_rclr_df$gb.ACA23185_1.ARO_3000194.tetW__Bifidobacterium_longum_)
#CDcard_All <- data.frame(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]

glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

#lm_ls$gb.AAM09851_1.ARO_3002923.vanRD__Enterococcus_faecium_

#m <- lm(glm_data[,taxon] ~ Group + Age, data=glm_data[,c(taxon,"Group", "Age")])

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
CARD_classes <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\CARD_Annotation.xlsx", 1)
CARD_classes$taxon_id <- CARD_classes$Complete.taxonomy
dat <- left_join(dat, CARD_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#view(dat)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_CARD_Unfiltered.xlsx")

#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_RTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_CARD_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_RTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_CARD_Filtered1%.xlsx")

#VFDB
meta_renal <- CDmeta_All %>%
  select(ID, Cross_sectional_samplesDag3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(Cross_sectional_samplesDag3 %in% c("Healthy Control","RTR")) %>%
  mutate(Group=relevel(factor(Cross_sectional_samplesDag3), ref="Healthy Control")) %>%
  select(-Cross_sectional_samplesDag3)
summary(meta_renal$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_vfdb[,colnames(ftbl_rclr_vfdb) %in% meta_renal$ID] # change between meta_renal and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]
# Species - Taxon-specific linear models 
# use meta_renal or meta_renal
glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat$VFID <- dat$taxon_id
dat$VFID <- str_split_fixed(str_split_fixed(dat$VFID, pattern="\\..__", 2)[,2], pattern="\\.", 2)[,1]
VFDB_Anno <- read.csv(file = "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\VFDB_annotation.csv", 
                      stringsAsFactors = FALSE)
dat <- left_join(dat, VFDB_Anno, by="VFID")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
view(dat)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_VFDB_Unfiltered.xlsx")

#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_RTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_VFDB_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_RTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA Cross-Sectional\\PA CARD VFDB\\Results_DA_RTR_HC_Full_Model_PA_RCLR_VFDB_Filtered1%.xlsx")






#3.6      DA Immunosuppresive medication analysis Users vs non Users----
#@Shixian codes for Immunosuppresive medication analysis go here
#3.7      DA ESLD vs HC----
#This is for cross-sectional (ESLD and ESRD) analysis plus covariates
meta_liver <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxLTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_liver$Group)
sum(is.na(meta_liver$Age))
sum(is.na(meta_liver$BMI))
sum(is.na(meta_liver$Gender))
sum(is.na(meta_liver$Smoking))
summary(meta_liver$Smoking)
sum(is.na(meta_liver$PPI))
sum(is.na(meta_liver$Laxatives))
sum(is.na(meta_liver$Antibiotics))

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_species[,colnames(ftbl_rclr_species) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

sum(is.na(meta_liver$Age))
sum(is.na(meta_liver$BMI))
sum(is.na(meta_liver$Gender))
sum(is.na(meta_liver$Smoking))
sum(is.na(meta_liver$PPI))
sum(is.na(meta_liver$Laxatives))
sum(is.na(meta_liver$Antibiotics))

summary(meta_liver$Age ~ meta_liver$Group)
sum(is.na(meta_liver$BMI))
sum(is.na(meta_liver$Gender))

with(meta_liver, tapply(Smoking, Group, summary))
with(meta_liver, tapply(PPI, Group, summary))
with(meta_liver, tapply(Laxatives, Group, summary))
with(meta_liver, tapply(Antibiotics, Group, summary))

sum(is.na(meta_liver$Smoking))
sum(is.na(meta_liver$PPI))
sum(is.na(meta_liver$Laxatives))
sum(is.na(meta_liver$Antibiotics))

tapply(meta_liver$Smoking, meta_liver$Group, summary)

wilcox.test(Age ~ Group, data=meta_liver)
wilcox.test(BMI ~ Group, data=meta_liver)
with(meta_liver, chisq.test(Gender, Group))
with(meta_liver, chisq.test(Smoking, Group))
with(meta_liver, chisq.test(PPI, Group))
with(meta_liver, chisq.test(Laxatives, Group))
with(meta_liver, chisq.test(Antibiotics, Group))

#Check intracorrelation covariates
row.names(meta_liver) <- meta_liver$ID
meta_liver$ID <- NULL
meta_liver$Gender <- as.numeric(meta_liver$Gender)
meta_liver$Smoking <- as.numeric(meta_liver$Smoking)
meta_liver$PPI <- as.numeric(meta_liver$PPI)
meta_liver$Laxatives <- as.numeric(meta_liver$Laxatives)
meta_liver$Antibiotics <- as.numeric(meta_liver$Antibiotics)
meta_liver$Group <- as.numeric(meta_liver$Group)

library(ggcorrplot)
corr <- round(cor(meta_liver), 1)
corr
p.mat <- round(cor_pmat(meta_liver), 1)
p.mat

ggcorrplot(corr, hc.order = F,
           type = "full", p.mat = p.mat)+
  ggtitle("ESLD")
ggcorrplot(corr)

#ggsave("J:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/RCLR DA V2/Correlation Heatmap/ESLD_HM.png", device = "png", type = "cairo", dpi=300)

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))#[,c(2,5,7)]
#colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")
colnames(species_taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
species_taxonomy$taxon_id=species_taxonomy$Species

Healthy <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\Healthy_merged.xlsx",1)
Healthy$Taxonomy <- NULL

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) 
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat <- left_join(dat, Healthy, by="taxon_id")
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  mutate(Health_Overlap_Gacesa=ifelse(Direction==Unhealthy_Direction_Gacesa,"Yes","No"))%>%
  mutate(Health_Overlap_Gupta=ifelse(Direction==Unhealthy_Direction_Gupta,"Yes","No")) %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`,
         `Butyrare producer` = Butyrate_Prod)
dat <- left_join(dat, species_taxonomy, by="taxon_id")

view(dat)  
#Write table
#write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_ESLD_HC_Full_Model_Species_V2.xlsx")

#Pathways
meta_liver <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxLTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_pwys[,colnames(ftbl_rclr_pwys) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age+ BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
pathway_classes <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\path_annotated_index_V3.txt", sep='\t',header=T, fill=T)
pathway_classes$taxon_id <- pathway_classes$MetaCyc_complete_name
dat <- left_join(dat, pathway_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_ESLD_HC_Full_Model_Pathways.xlsx")
#3.8      DA ESLD vs HC CARD VFDB PA----
#CARD
meta_liver <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxLTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_liver$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_card[,colnames(ftbl_rclr_card) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
#ftbl_rclr_df <- data.frame(ftbl_rclr)
#sum(is.na(ftbl_rclr_df$gb.ACA23185_1.ARO_3000194.tetW__Bifidobacterium_longum_)
#CDcard_All <- data.frame(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]

glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

tapply(glm_data$gb.ABA71729_1.ARO_3002959.vanYG1__Enterococcus_faecalis_, glm_data$Group, sum)

#taxon="gb.CAE51745_1.ARO_3000168.tet.D.__Serratia_marcescens_"

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

#lm_ls$gb.AAM09851_1.ARO_3002923.vanRD__Enterococcus_faecium_

#m <- lm(glm_data[,taxon] ~ Group + Age, data=glm_data[,c(taxon,"Group", "Age")])

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
CARD_classes <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\CARD_Annotation.xlsx", 1)
CARD_classes$taxon_id <- CARD_classes$Complete.taxonomy
dat <- left_join(dat, CARD_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
view(dat)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESLD_HC_Full_Model_PA_RCLR_CARD_V2_Unfiltered.xlsx")
#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_PreTxLTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESLD_HC_Full_Model_PA_RCLR_CARD_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_PreTxLTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ELSD_HC_Full_Model_PA_RCLR_CARD_Filtered1%.xlsx")

#VFDB
meta_liver <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxLTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_liver$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_vfdb[,colnames(ftbl_rclr_vfdb) %in% meta_liver$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]
# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_liver, by="ID")

tapply(glm_data$gb.ABA71729_1.ARO_3002959.vanYG1__Enterococcus_faecalis_, glm_data$Group, sum)

#taxon="gb.CAE51745_1.ARO_3000168.tet.D.__Serratia_marcescens_"

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat$VFID <- dat$taxon_id
dat$VFID <- str_split_fixed(str_split_fixed(dat$VFID, pattern="\\..__", 2)[,2], pattern="\\.", 2)[,1]
VFDB_Anno <- read.csv(file = "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\VFDB_annotation.csv", 
                      stringsAsFactors = FALSE)
dat <- left_join(dat, VFDB_Anno, by="VFID")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESLD_HC_Full_Model_PA_RCLR_VFDB_Unfiltered.xlsx")
#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_PreTxLTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESLD_HC_Full_Model_PA_RCLR_VFDB_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_PreTxLTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESLD_HC_Full_Model_PA_RCLR_VFDB_Filtered1%.xlsx")

#3.9      DA ESRD vs HC----
#This is for cross-sectional (ESLD and ESRD) analysis plus covariates
meta_renal <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxRTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_renal$Group)
sum(is.na(meta_renal$Age))
sum(is.na(meta_renal$BMI))
sum(is.na(meta_renal$Gender))
sum(is.na(meta_renal$Smoking))
summary(meta_renal$Smoking)
sum(is.na(meta_renal$PPI))
sum(is.na(meta_renal$Laxatives))
sum(is.na(meta_renal$Antibiotics))

with(meta_renal, tapply(Smoking, Group, summary))
with(meta_renal, tapply(PPI, Group, summary))
with(meta_renal, tapply(Laxatives, Group, summary))
with(meta_renal, tapply(Antibiotics, Group, summary))

wilcox.test(Age ~ Group, data=meta_renal)
wilcox.test(BMI ~ Group, data=meta_renal)
with(meta_renal, chisq.test(Gender, Group))
with(meta_renal, chisq.test(Smoking, Group))
with(meta_renal, chisq.test(PPI, Group))
with(meta_renal, chisq.test(Laxatives, Group))
with(meta_renal, chisq.test(Antibiotics, Group))
#meta_renal[is.na(meta_liver$Smoking),]$Smoking <- "Nee"


# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_species[,colnames(ftbl_rclr_species) %in% meta_renal$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_df <- data.frame(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))#[,c(2,5,7)]
#colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")
colnames(species_taxonomy) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
species_taxonomy$taxon_id=species_taxonomy$Species
Healthy <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\Healthy_merged.xlsx",1) 
Healthy$Taxonomy <- NULL

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) 
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat <- left_join(dat, Healthy, by="taxon_id")
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  mutate(Health_Overlap_Gacesa=ifelse(Direction==Unhealthy_Direction_Gacesa,"Yes","No"))%>%
  mutate(Health_Overlap_Gupta=ifelse(Direction==Unhealthy_Direction_Gupta,"Yes","No")) %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`,
         `Butyrare producer` = Butyrate_Prod)
dat <- left_join(dat, species_taxonomy, by="taxon_id")
#view(dat)  
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_ESRD_HC_Full_Model_Species_V2.xlsx")

#Pathways
meta_renal <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxRTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)

# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_pwys[,colnames(ftbl_rclr_pwys) %in% meta_renal$ID] # change between meta_liver and meta_renal
ftbl_rclr <- t(ftbl_rclr)

# Species - Taxon-specific linear models 
# use meta_liver or meta_renal
lm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  lm_data <- data.frame(ftbl_rclr) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lm(lm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=lm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")])
    
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
pathway_classes <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\path_annotated_index_V3.txt", sep='\t',header=T, fill=T)
pathway_classes$taxon_id <- pathway_classes$MetaCyc_complete_name
dat <- left_join(dat, pathway_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\Results_DA_ESRD_HC_Full_Model_Pathways.xlsx")
#3.10     DA ESRD vs HC CARD VFDB PA----
#CARD
meta_renal <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxRTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_renal$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_card[,colnames(ftbl_rclr_card) %in% meta_renal$ID] # change between meta_renal and meta_renal
ftbl_rclr <- t(ftbl_rclr)
#ftbl_rclr_df <- data.frame(ftbl_rclr)
#sum(is.na(ftbl_rclr_df$gb.ACA23185_1.ARO_3000194.tetW__Bifidobacterium_longum_)
#CDcard_All <- data.frame(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]

glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

#lm_ls$gb.AAM09851_1.ARO_3002923.vanRD__Enterococcus_faecium_

#m <- lm(glm_data[,taxon] ~ Group + Age, data=glm_data[,c(taxon,"Group", "Age")])

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
CARD_classes <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\CARD_Annotation.xlsx", 1)
CARD_classes$taxon_id <- CARD_classes$Complete.taxonomy
dat <- left_join(dat, CARD_classes, by="taxon_id")
#view(dat)
dat %>%
  filter(p_adj_sig=="Yes")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#view(dat)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_CARD_Unfiltered.xlsx")
#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_PreTxRTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_CARD_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_PreTxRTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_CARD_Filtered1%.xlsx")

#VFDB
meta_renal <- CDmeta_All %>%
  select(ID, TypePreTx_DAG3, Age, BMI, Gender, Smoking, PPI, Laxatives, Antibiotics) %>%
  filter(TypePreTx_DAG3 %in% c("Healthy Control","PreTxRTR")) %>%
  mutate(Group=relevel(factor(TypePreTx_DAG3), ref="Healthy Control")) %>%
  select(-TypePreTx_DAG3)
summary(meta_renal$Group)
# subset abundance table to focal samples
ftbl_rclr <- ftbl_rclr_vfdb[,colnames(ftbl_rclr_vfdb) %in% meta_renal$ID] # change between meta_renal and meta_renal
ftbl_rclr <- t(ftbl_rclr)
ftbl_rclr_PA <- data.frame(ftbl_rclr)

for (cn in colnames(ftbl_rclr_PA)) {
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] > 0] <- 1
  ftbl_rclr_PA[[cn]][ftbl_rclr_PA[[cn]] <= 0] <- 1
}
ftbl_rclr_PA[is.na(ftbl_rclr_PA)] <- 0
ftbl_rclr_PA <- ftbl_rclr_PA[,colSums(ftbl_rclr_PA)>0]
# Species - Taxon-specific linear models 
# use meta_renal or meta_renal
glm_data <- #data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  glm_data <- data.frame(ftbl_rclr_PA) %>%
  rownames_to_column("ID") %>%
  left_join(meta_renal, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr_PA))
names(lm_ls) <- colnames(ftbl_rclr_PA)
for(taxon in colnames(ftbl_rclr_PA)) {
  tryCatch({
    print(taxon)
    m <- glm(glm_data[,taxon] ~ Group + Age + BMI + Gender + Smoking + PPI + Laxatives + Antibiotics, data=glm_data[,c(taxon,"Group", "Age", "BMI", "Gender", "Smoking", "PPI", "Laxatives", "Antibiotics")], family = binomial)
    mm <- coef(summary(m))[2,c(1,2,4)]
    names(mm)[3] <- "p_val"
    a <- tapply(glm_data[,taxon], glm_data$Group, sum)
    
    mm[4] <- a[1]
    names(mm)[4] <- paste0("n_", names(a)[1])
    mm[5] <- a[1] / summary(glm_data$Group)[1]
    names(mm)[5] <- paste0("%_", names(a)[1])
    mm[6] <- a[2] 
    names(mm)[6] <- paste0("n_", names(a)[2])
    mm[7] <- a[2] / summary(glm_data$Group)[2]
    names(mm)[7] <- paste0("%_", names(a)[2])
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# ERROR : contrasts can be applied only to factors with 2 or more levels
# This is caused by to few samples in e.g. Smoking variable
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>%
  bind_rows(.id="taxon_id") %>%
  mutate(p_adj=p.adjust(p_val, "BH"))

#lm_df[lm_df$taxon_id=="gb.AAG03357_1.ARO_3000596.ErmX__Corynebacterium_striatum_",]

dat <- lm_df %>%
  mutate(sign=sign(Estimate),
         Direction=ifelse(sign=="1","Increased","Decreased"),
         p_adj_sig=ifelse(p_adj<0.1, "Yes", "No"),
         taxon_id=fct_reorder(taxon_id, Estimate, median))
dat$sign <- NULL
dat %>%
  filter(p_adj_sig=="Yes")
dat$VFID <- dat$taxon_id
dat$VFID <- str_split_fixed(str_split_fixed(dat$VFID, pattern="\\..__", 2)[,2], pattern="\\.", 2)[,1]
VFDB_Anno <- read.csv(file = "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Metacyc- VFDB- CARD- Additional Info\\VFDB_annotation.csv", 
                      stringsAsFactors = FALSE)
dat <- left_join(dat, VFDB_Anno, by="VFID")
dat <- dat %>%
  rename(`P-value`=p_val, 
         PFDR = p_adj,
         `Significant (PFDR<0.1)` = p_adj_sig,
         `Standard error` = `Std. Error`)
#view(dat)
#Write table
write.xlsx2(dat, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_VFDB_Unfiltered.xlsx")
#Filter on prevalence here
lm_df_filter01 <- dat %>% 
  filter(`%_PreTxRTR`>0.1 & `%_Healthy Control`>0.1)
summary(as.factor(lm_df_filter01$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter01, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_VFDB_Filtered10%.xlsx")

lm_df_filter001 <- dat %>% 
  filter(`%_PreTxRTR`>0.01 & `%_Healthy Control`>0.01)
summary(as.factor(lm_df_filter001$`Significant (PFDR<0.1)`))
write.xlsx2(lm_df_filter001, "J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\RCLR DA V2\\DA ESD\\PA CARD VFDB ESD\\Results_DA_ESRD_HC_Full_Model_PA_RCLR_VFDB_Filtered1%.xlsx")


#3.11     LONGTIUDINAL DIFFERNETIAL ABUNDANCE ANALYSIS GOES HERE PLEASE----
#4.0      Calculation of microbiome variance explained by phenotypes----
#4.1      ADONIS - Within LTx Oud----
setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/ADONIS TransplantLines")
ADONIS_Within_LTxOld <- subset(CDmeta, Cross_sectional_samples =="LTR")
ADONIS_Within_LTxOld_Taxa <- ADONIS_Within_LTxOld[c("ID")]
inDF <- CDtaxa
spDF <- filterMetaGenomeDF(inDF=subsetMicrobiomeDF(inDF,getTaxa = T,getPWYs = F,getCARDs = F,getVFs = F,getDivs = F,getPhenos = F,verbose = T),
                           minMRelAb = 0,presPerc = 0,minMedRelAb = 0,keepLevels = "S")
ADONIS_Within_LTxOld_Taxa <- merge(ADONIS_Within_LTxOld_Taxa, spDF, by="row.names")
rownames(ADONIS_Within_LTxOld_Taxa) <- ADONIS_Within_LTxOld_Taxa$ID
ADONIS_Within_LTxOld_meta <- ADONIS_Within_LTxOld_Taxa[c("ID")]
rownames(ADONIS_Within_LTxOld_meta ) <- ADONIS_Within_LTxOld_meta$ID
ADONIS_Within_LTxOld_meta <- merge(ADONIS_Within_LTxOld_meta, CDmeta, by="row.names")
rownames(ADONIS_Within_LTxOld_meta) <- ADONIS_Within_LTxOld_meta$ID.x
ADONIS_Within_LTxOld_Taxa$ID <- NULL
ADONIS_Within_LTxOld_Taxa$Row.names <- NULL
ADONIS_Within_LTxOld_meta$ID <- ADONIS_Within_LTxOld_meta$ID.x

adonisVarsTouse <-   select(ADONIS_Within_LTxOld_meta,
                            "Age",
                            "Gender",
                            "LaxativesYears_since_transplantation",
                            "BMI",
                            "BSA",
                            "BIAWAISTHIP",
                            "BIAFATPC",
                            "SBP",
                            "Proteinuria",
                            "Retransplantation",
                            "func_constipation",
                            "IBS",
                            "func_diarrhea",
                            "Hypertension",
                            "Alcohol",
                            "Smoking",
                            "Diabetes",
                            "GGT",
                            "ALT",
                            "AST",
                            "ALP",
                            "Tbil",
                            "CRP",
                            "Leukocytes",
                            "eGFR",
                            "CICLOVBLD",
                            "MYCO_PLAS",
                            "SIROL_VBLD",
                            "Serum_Tacrolimus",
                            "HDL",
                            "LDL",
                            "Total_cholesterol",
                            "NT_proBNP",
                            "Troponin_T",
                            "PPI",
                            "Statins",
                            "Laxatives",
                            "Antibiotics",
                            "Mycophenolic_acid",
                            "Azathiorpine",
                            "Ciclosporin",
                            "Tacrolimus",
                            "Prednisolone",
                            "Sirolimus",
                            "Tac_Pred",
                            "Cyclo_Pred",
                            "MMF_Cyclo",
                            "MMF_Pred",
                            "MMF_Tac",
                            "MMF_Cyclo_Pred",
                            "MMF_Tac_Pred",
                            "Aza_Tac_Pred",)

#ADONIS
rownames(adonisVarsTouse) <- rownames(ADONIS_Within_LTxOld_meta)
# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  varToUse <- colnames(adonisVarsTouse)[[i]]
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  # make new DF with selected phenotype and taxa 
  tmpDF <- merge.data.frame(ADONIS_Within_LTxOld_meta[,colnames(ADONIS_Within_LTxOld_meta) %in% c("ID",varToUse)] ,
                            spDF,by.x="ID",by.y = "row.names")
  # get rid of ID (used only for merging)
  tmpDF$ID <- NULL
  # get rid of rows with NAs
  print('debug, DF size before removing NAs:')
  print(dim(tmpDF))
  tmpDF <- tmpDF[complete.cases(tmpDF),]
  print('debug, DF size after removing NAs:')
  print(dim(tmpDF))
  # calculate BC matrix
  tmpDFforBC <- tmpDF
  tmpDFforBC[[varToUse]] <- NULL
  print('calculating BC matrix')
  bcMat <- vegdist(tmpDFforBC,method = "bray")
  print(paste0('doing adonis for ',varToUse))
  ad<-adonis(bcMat ~ adonisVarsTouse[,i],permutations=9999)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}


adonis_meta= as.data.frame(adonis_meta)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
rownames(adonis_meta) = colnames(adonisVarsTouse)

# View(adonis_meta)
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
plot(g)
ggsave(plot = g,filename = "ADONIS_results_Within_LTxOud_Selected_Variables.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.xlsx(adonis_meta,"ADONIS_results_Within_LTxOud_Selected_Variables.xlsx")




rownames(adonisVarsTouse) <- rownames(ADONIS_Within_LTxOld_meta)
#Drop factors 0 entries in level :)
adonisVarsTouse <- droplevels(adonisVarsTouse)
#Drop factors with only 1 level
adonisVarsTouse <- adonisVarsTouse[, sapply(adonisVarsTouse, nlevels) != 1]
#bcMat <- vegdist(spDF,method = "bray")
#adonis <- adonis(bcurtis ~ .,data=adonisVarsTouse,permutations=100,parallel=F)

# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  varToUse <- colnames(adonisVarsTouse)[[i]]
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  # make new DF with selected phenotype and taxa 
  tmpDF <- merge.data.frame(ADONIS_Within_LTxOld_meta[,colnames(ADONIS_Within_LTxOld_meta) %in% c("ID",varToUse)] ,
                            spDF,by.x="ID",by.y = "row.names")
  # get rid of ID (used only for merging)
  tmpDF$ID <- NULL
  # get rid of rows with NAs
  print('debug, DF size before removing NAs:')
  print(dim(tmpDF))
  tmpDF <- tmpDF[complete.cases(tmpDF),]
  print('debug, DF size after removing NAs:')
  print(dim(tmpDF))
  # calculate BC matrix
  tmpDFforBC <- tmpDF
  tmpDFforBC[[varToUse]] <- NULL
  print('calculating BC matrix')
  bcMat <- vegdist(tmpDFforBC,method = "bray")
  print(paste0('doing adonis for ',varToUse))
  ad<-adonis(bcMat ~ adonisVarsTouse[,i],permutations=2500)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}


adonis_meta= as.data.frame(adonis_meta)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
rownames(adonis_meta) = colnames(adonisVarsTouse)
adonis_meta <- adonis_meta[order(adonis_meta$`FDR(BH)`),]

# View(adonis_meta)
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
plot(g)
ggsave(plot = g,filename = "ADONIS_results_Within_LTxOud_Selected_Variables.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.xlsx(adonis_meta,"ADONIS_results_Within_LTxOud_Selected_Variables.xlsx")

#4.2      ADONIS - Within NTx Oud---- 
setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/ADONIS TransplantLines")
ADONIS_Within_NTxOld <- subset(CDmeta, Cross_sectional_samples =="RTR")
ADONIS_Within_NTxOld_Taxa <- ADONIS_Within_NTxOld[c("ID")]
inDF <- CDtaxa
spDF <- filterMetaGenomeDF(inDF=subsetMicrobiomeDF(inDF,getTaxa = T,getPWYs = F,getCARDs = F,getVFs = F,getDivs = F,getPhenos = F,verbose = T),
                           minMRelAb = 0,presPerc = 0,minMedRelAb = 0,keepLevels = "S")
ADONIS_Within_NTxOld_Taxa <- merge(ADONIS_Within_NTxOld_Taxa, spDF, by="row.names")
rownames(ADONIS_Within_NTxOld_Taxa) <- ADONIS_Within_NTxOld_Taxa$ID
ADONIS_Within_NTxOld_meta <- ADONIS_Within_NTxOld_Taxa[c("ID")]
rownames(ADONIS_Within_NTxOld_meta ) <- ADONIS_Within_NTxOld_meta$ID
ADONIS_Within_NTxOld_meta <- merge(ADONIS_Within_NTxOld_meta, CDmeta, by="row.names")
rownames(ADONIS_Within_NTxOld_meta) <- ADONIS_Within_NTxOld_meta$ID.x
ADONIS_Within_NTxOld_Taxa$ID <- NULL
ADONIS_Within_NTxOld_Taxa$Row.names <- NULL
ADONIS_Within_NTxOld_meta$ID <- ADONIS_Within_NTxOld_meta$ID.x

adonisVarsTouse <-   select(ADONIS_Within_NTxOld_meta,
                            "Age",
                            "Gender",
                            "LaxativesYears_since_transplantation",
                            "BMI",
                            "BSA",
                            "BIAWAISTHIP",
                            "BIAFATPC",
                            "Handgrip_strenght",
                            "SBP",
                            "Proteinuria",
                            "Retransplantation",
                            "func_constipation",
                            "IBS",
                            "func_diarrhea",
                            "Hypertension",
                            "Alcohol",
                            "Smoking",
                            "Diabetes",
                            "GGT",
                            "ALT",
                            "AST",
                            "ALP",
                            "Tbil",
                            "CRP",
                            "Leukocytes",
                            "eGFR",
                            "CICLOVBLD",
                            "EVERO_VBLD",
                            "MYCO_PLAS",
                            "Serum_Tacrolimus",
                            "HDL",
                            "LDL",
                            "Total_cholesterol",
                            "NT_proBNP",
                            "Troponin_T",
                            "PPI",
                            "Statins",
                            "Laxatives",
                            "Antibiotics",
                            "Mycophenolic_acid",
                            "Azathiorpine",
                            "Ciclosporin",
                            "Tacrolimus",
                            "Prednisolone",
                            "Tac_Pred",
                            "Cyclo_Pred",
                            "MMF_Cyclo",
                            "MMF_Pred",
                            "MMF_Tac",
                            "MMF_Cyclo_Pred",
                            "MMF_Tac_Pred",
                            "Aza_Tac_Pred",
                            "Eve_Tac_Pred")

#ADONIS
rownames(adonisVarsTouse) <- rownames(ADONIS_Within_NTxOld_meta)
# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  varToUse <- colnames(adonisVarsTouse)[[i]]
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  # make new DF with selected phenotype and taxa 
  tmpDF <- merge.data.frame(ADONIS_Within_NTxOld_meta[,colnames(ADONIS_Within_NTxOld_meta) %in% c("ID",varToUse)] ,
                            spDF,by.x="ID",by.y = "row.names")
  # get rid of ID (used only for merging)
  tmpDF$ID <- NULL
  # get rid of rows with NAs
  print('debug, DF size before removing NAs:')
  print(dim(tmpDF))
  tmpDF <- tmpDF[complete.cases(tmpDF),]
  print('debug, DF size after removing NAs:')
  print(dim(tmpDF))
  # calculate BC matrix
  tmpDFforBC <- tmpDF
  tmpDFforBC[[varToUse]] <- NULL
  print('calculating BC matrix')
  bcMat <- vegdist(tmpDFforBC,method = "bray")
  print(paste0('doing adonis for ',varToUse))
  ad<-adonis(bcMat ~ adonisVarsTouse[,i],permutations=9999)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}


adonis_meta= as.data.frame(adonis_meta)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
rownames(adonis_meta) = colnames(adonisVarsTouse)

# View(adonis_meta)
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
plot(g)
ggsave(plot = g,filename = "ADONIS_results_Within_NTxOud_Selected_Variables.png")
write.xlsx(adonis_meta,"ADONIS_results_Within_NTxOud_Selected_Variables.xlsx")
#5.0      Longitudinal Analysis----
#5.1      Figure 6----
library(foreign)
library(stringr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(scales)
library(gsubfn)
library(xlsx)

#5.2      Load data----
#Nature Medicine Longitudinal Main Figure
CDtaxa <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_metaphlan_cleaned_filtered.txt",sep='\t',header=T, fill=T)
CDtaxa$ID <- as.character(CDtaxa$ID)
trimws(CDtaxa$ID)
rownames(CDtaxa) <- CDtaxa$ID

CDpwys <- read.table("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_humann_pathabundance_filtered.txt",sep='\t',header=T)
CDpwys$ID <- as.character(CDpwys$ID)
trimws(CDpwys$ID)
rownames(CDpwys) <- CDpwys$ID

CDmeta = read.spss("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\Metadata.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

summary(CDmeta$Selection_Long_Cross_S)
CDmeta_All <- subset(CDmeta, Selection_Long_Cross_S == "Yes")
#Taxa
CDtaxa_All <- CDmeta_All["ID"]
CDtaxa_All <- merge(CDtaxa_All, CDtaxa, by="ID")
row.names(CDtaxa_All) <- CDtaxa_All$ID

#Pathways
CDpwys_All <- CDmeta_All["ID"]
CDpwys_All <- merge(CDpwys_All, CDpwys, by="ID")
row.names(CDpwys_All) <- CDpwys_All$ID
CDpwys_All$ID <- NULL

CDmeta_All <- CDtaxa_All["ID"]
CDmeta_All <- merge(CDmeta_All, CDmeta, by="ID")

CDtaxa_All$ID <- NULL

CDmeta_All <- CDmeta_All %>% 
  mutate(Cross_sectional_samplesDag3=factor(Cross_sectional_samplesDag3, levels=c("Healthy Control","LTR","RTR","")),
         exclude_liver_long=paste0(Longitudinal_Timepoint, "_", TypeTxN_TxL_TxD_DAG3))
#5.3      Ordination Deicode----
#5.4      Fig 6A Ordination Taxa----
centroids_fig2b <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDpwys_All/n2/CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDcard_All/n2/CDcard_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDvfdb_All/n2/CDvfdb_All_sample_loadings.csv", row.names=1) %>% 
  
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), by="ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"))) %>% 
  group_by(exclude_liver_long) %>% 
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2))

Fig4A <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDpwys_All/n2/CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDcard_All/n2/CDcard_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDvfdb_All/n2/CDvfdb_All_sample_loadings.csv", row.names=1) %>% 
  
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), "ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"))) %>% 
  ggplot(aes(x=PC1, y=PC2, color=exclude_liver_long)) +
  geom_point(shape=19, size=3, alpha=0.5) +
  scale_color_manual(values = c("Healthy_Control_DAG3 Control"="darkgray", "PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                                "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  theme_classic()+ 
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"))+
  theme(legend.text = element_text(face = "bold", size=12),
        legend.title = element_blank(),
        plot.title=element_text(face = "bold",hjust=0.5),
        legend.position="none")+
  #xlab("PC1: 15.4 %") +
  #ylab("PC2: 11.2 %") +
  #ggtitle("Robust Aitchison PCA (DECOIDE)") +
  geom_point(data=centroids_fig2b, aes(PC1_mean, PC2_mean, fill=exclude_liver_long), shape=21, size=6, color="black") +
  scale_fill_manual(values = c("Healthy_Control_DAG3 Control"="darkgray", "PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                               "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.08, -0.04, 0.0,0.02, 0.04),limits=c(-0.08, 0.04))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.05,-0.025,0.0,0.025,0.05), limits=c(-0.05,0.05))
#stat_ellipse(type="norm", level=0.95)
Fig4A

#Distance to healthy PC1 Taxa
centroids_fig2d <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_sample_loadings.csv", row.names=1) %>% 
  #read.csv("~/Documents/UMGC/transplant_study/all_data/output/all_samples/CDpwys_All/n2/CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), by="ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"))) %>% 
  group_by(exclude_liver_long) %>% 
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2)) %>% 
  mutate(PC2_mean=0)

# Distance to HC along e.g. PC1
# only comparing it in one plane (i.e. along PC1)
euclidean_dist_centroids_fig2d <- as.matrix(dist(centroids_fig2d[,c("PC1_mean","PC2_mean")], method="euclidean"))
#euclidean_dist_centroids_fig2d <- euclidean_dist_centroids_fig2d/max(euclidean_dist_centroids_fig2d)
rownames(euclidean_dist_centroids_fig2d) <- colnames(euclidean_dist_centroids_fig2d) <- c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL")

Fig4A2 <- data.frame(dist_to_hc=euclidean_dist_centroids_fig2d[1,2:7]) %>% 
  rownames_to_column("group") %>% 
  mutate(group=factor(group, levels=c("LTR_CS_TxL","RTR_CS_TxN","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"))) %>% 
  mutate(dist_to_hc_norm=dist_to_hc/max(dist_to_hc),
         X=1) %>% 
  ggplot(aes(y=dist_to_hc_norm, x=X)) +
  coord_flip() +
  geom_vline(xintercept=1) +
  geom_point(size=8, aes(color=group)) +
  #geom_line(aes(group=1)) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black", angle=45, hjust=1),
    axis.text.y = element_blank())+
  theme(legend.text = element_text(face = "bold", size=12),
        plot.title=element_text(size=14, face="bold", colour = "black", hjust=0.5),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"))+
  labs(x=NULL, y="Average distance to healthy control") +
  scale_colour_manual(values = c("PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                                 "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  scale_y_continuous(limits=c(0,1),expand=c(0.01,0,0,0.049), breaks=seq(from=0, to=1, by=0.25)) +
  geom_text(aes(label=group), angle=45, hjust=-0.15, vjust=-1, fontface="bold", size=6)
Fig4A2
setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 4")
ggsave("Fig4A2_MEAN.png", device = "png", type = "cairo", dpi=300)
#5.5      Fig 6B Ordination Pathways----
centroids_fig2b <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), by="ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"))) %>% 
  group_by(exclude_liver_long) %>% 
  summarise(PC1_mean=mean(PC1), PC2_mean=mean(PC2))

Fig4B <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), "ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL","Healthy_Control_DAG3 Control"))) %>% 
  ggplot(aes(x=PC1, y=PC2, color=exclude_liver_long)) +
  geom_point(shape=19, size=3, alpha=0.5) +
  scale_color_manual(values = c("Healthy_Control_DAG3 Control"="darkgray", "PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                                "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  theme_classic()+ 
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"))+
  theme(legend.text = element_text(face = "bold", size=12),
        legend.title = element_blank(),
        plot.title=element_text(face = "bold",hjust=0.5),
        legend.position="none")+
  #xlab("PC1: 15.4 %") +
  #ylab("PC2: 11.2 %") +
  #ggtitle("Robust Aitchison PCA (DECOIDE)") +
  geom_point(data=centroids_fig2b, aes(PC1_mean, PC2_mean, fill=exclude_liver_long), shape=21, size=6, color="black") +
  scale_fill_manual(values = c("Healthy_Control_DAG3 Control"="darkgray", "PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                               "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.05, -0.025, 0.0,0.025, 0.05),limits=c(-0.05, 0.05))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'), breaks=c(-0.05,-0.025,0.0,0.025,0.05), limits=c(-0.05,0.05))
#stat_ellipse(type="norm", level=0.95)
Fig4B

#Distance to healthy PC1 Taxa
centroids_fig2d <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDpwys_All_sample_loadings.csv", row.names=1) %>% 
  rownames_to_column("ID") %>% 
  select(ID, PC1, PC2) %>% 
  filter(ID %in% CDmeta_All[CDmeta_All$exclude_liver_long %in% c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"),]$ID) %>%
  left_join(select(CDmeta_All, ID, exclude_liver_long), by="ID") %>%
  mutate(exclude_liver_long=factor(exclude_liver_long, levels=c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"))) %>% 
  group_by(exclude_liver_long) %>% 
  summarise(PC1_mean=median(PC1), PC2_mean=median(PC2)) %>% 
  mutate(PC2_mean=0)

# Distance to HC along e.g. PC1
# only comparing it in one plane (i.e. along PC1)
euclidean_dist_centroids_fig2d <- as.matrix(dist(centroids_fig2d[,c("PC1_mean","PC2_mean")], method="euclidean"))
#euclidean_dist_centroids_fig2d <- euclidean_dist_centroids_fig2d/max(euclidean_dist_centroids_fig2d)
rownames(euclidean_dist_centroids_fig2d) <- colnames(euclidean_dist_centroids_fig2d) <- c("Healthy_Control_DAG3 Control","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL")

Fig4B2 <- data.frame(dist_to_hc=euclidean_dist_centroids_fig2d[1,2:7]) %>% 
  rownames_to_column("group") %>% 
  mutate(group=factor(group, levels=c("LTR_CS_TxL","RTR_CS_TxN","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","PreTxLTR_TxL"))) %>% 
  mutate(dist_to_hc_norm=dist_to_hc/max(dist_to_hc),
         X=1) %>% 
  ggplot(aes(y=dist_to_hc_norm, x=X)) +
  coord_flip() +
  geom_vline(xintercept=1) +
  geom_point(size=8, aes(color=group)) +
  #geom_line(aes(group=1)) +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black", angle=45, hjust=1),
    axis.text.y = element_blank())+
  theme(legend.text = element_text(face = "bold", size=12),
        plot.title=element_text(size=14, face="bold", colour = "black", hjust=0.5),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"))+
  labs(x=NULL, y="Average distance to healthy control") +
  scale_colour_manual(values = c("PreTxRTR_TxN"="#FF0000","PreTxLTR_TxL"="#FFA600","3months_TxN"="#eff3ff","6months_TxN"="#c6dbef",
                                 "12months_TxN"="#9ecae1","24months_TxN"="#6bafd6"))+
  scale_y_continuous(limits=c(0,1),expand=c(0.01,0,0,0.049), breaks=seq(from=0, to=1, by=0.25)) +
  geom_text(aes(label=group), angle=45, hjust=-0.15, vjust=-1, fontface="bold", size=6)
Fig4B2
setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 4")
ggsave("Fig4B2.png", device = "png", type = "cairo", dpi=300)
#5.6      Fig 6C Log Ratio Plot Taxa---- 
#Taxa
ftbl_rclr <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_rclr.csv", row.names=1)

CDtaxa_All_only_s <-  CDtaxa_All[,grep(x=colnames(CDtaxa_All),pattern="s__")]
row.names(CDtaxa_All_only_s)
rownames(ftbl_rclr) <- str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7)[,7] #colnames(CDtaxa_All_only_s)
colnames(ftbl_rclr) <- rownames(CDtaxa_All_only_s)

meta <- CDmeta_All %>% 
  select(ID, Timepoint, Patient, exclude_liver_long) %>% 
  filter(exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN")) %>% 
  mutate(TimePoint=recode_factor(Timepoint,"PreTxRTR"="PreTx","3months"="3M","6months"="6M","12months"="12M","24months"="24M")) %>% 
  na.omit()

# subset to renal longi samples
ftbl_rclr <- ftbl_rclr[,colnames(ftbl_rclr) %in% meta$ID]
ftbl_rclr <- t(ftbl_rclr)

tmp_mat <- diag(length(colnames(ftbl_rclr)))
rownames(tmp_mat) <- colnames(tmp_mat) <- colnames(ftbl_rclr)
pairs <- reshape2::melt(tmp_mat) %>% 
  filter(Var1==Var2) %>% 
  select(numerator=Var1, denominator=Var2) %>% 
  arrange(numerator, denominator) %>% 
  mutate(numerator=as.character(numerator), denominator=as.character(denominator))

pairs$TimePoint <- NA
log_ratio_samples <- vector("list", nrow(ftbl_rclr))
names(log_ratio_samples) <- rownames(ftbl_rclr)
log_ratio_samples <- lapply(log_ratio_samples, function(x) pairs)

for(sample in names(log_ratio_samples)) {
  print(sample)
  for(pair in 1:nrow(log_ratio_samples[[sample]])) {
    log_ratio_samples[[sample]][pair,]$TimePoint <- as.character(meta[meta$ID==sample,]$TimePoint)
  }
}

pairs_across_samples <- bind_rows(log_ratio_samples, .id="ID")

log_ratio_self <- vector("list", ncol(ftbl_rclr))
names(log_ratio_self) <- colnames(ftbl_rclr)
log_ratio_self <- lapply(log_ratio_self, function(x) data.frame(clr=rep(NA,(length(unique(meta$TimePoint))-1)),
                                                                clr_std=rep(NA,(length(unique(meta$TimePoint))-1)),
                                                                TimePoint=c("3M","6M","12M","24M")))
for(taxon in colnames(ftbl_rclr)) {
  print(taxon)
  for(v in c("3M","6M","12M","24M")) {
    #log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$log_ratio <- log2( EnvStats::geoMean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon]) / EnvStats::geoMean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon]) )
    #log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$log_ratio_sd <- log( EnvStats::geoSD(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon]) / EnvStats::geoSD(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon]) )
    log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$clr <- mean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon], na.rm=T) - mean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon], na.rm=T)
    log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$clr_std <- sd(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon], na.rm=T) - sd(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon], na.rm=T)
  }
}

## Taxon-specific linear models

#lm_data <- data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
lm_data <- data.frame(ftbl_rclr) %>% 
  rownames_to_column("ID") %>%
  left_join(meta, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lmerTest::lmer(lm_data[,taxon]~TimePoint+(1|Patient), data=lm_data[,c(taxon,"TimePoint","Patient")])
    
    mm <- data.frame(coef(summary(m)))
    colnames(mm)[5] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# remove NULL lists - these are generated by to few samples: ERROR : number of levels of each grouping factor must be < number of observations (problems: Patient) 
# in case of pathways, there are 4
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>% 
  bind_rows(.id="taxon_id") %>% 
  rownames_to_column("group") %>% 
  separate(group, into=c("group",NA)) %>% 
  filter(group!="") %>% 
  mutate(p_adj=p.adjust(p_val, "BH"))

lm_df_sig <- lm_df %>% 
  filter(p_adj<0.1)

lm_sig <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  mutate(taxon_id=fct_reorder(taxon_id, Estimate, median)) %>% 
  select(taxon_id) %>%
  pull() %>% 
  levels() %>% 
  rev()

#filter out taxa with average effect size of below 0.1
lm_sig2 <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  group_by(taxon_id) %>% 
  summarise(mean_es=mean(Estimate)) %>% 
  filter(mean_es >= 0.1 | mean_es <= -0.1) %>% 
  arrange(desc(mean_es)) %>% 
  select(taxon_id) %>%
  pull()

species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))[,c(2,5,7)]
colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")

dat <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  
  #left_join(pathway_classes, by="taxon_id") %>% 
  left_join(species_taxonomy, by="taxon_id") %>% 
  
  
  mutate(sign=sign(Estimate),
         sign=ifelse(sign=="1","+","-"),
         group=str_split(group, "TimePoint", simplify=T)[,2],
         #group=fct_relevel(group, c("3M","6M","12M","24M")),
         group=fct_relevel(group, c("24M","12M","6M","3M")),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) %>%
  filter(taxon_id %in% lm_sig2) 

dat %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile() +
  scale_fill_gradient2(name="Effect size", low = "#AA1111", mid = "white", high = "#1111AA") +
  geom_text(aes(label=sign), size=2) +
  theme(legend.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        #axis.ticks.length=unit(0.25,"cm"), 
        axis.text.x=element_text(color="black",hjust=1),
        axis.text.y=element_text(color="black", size=5)) +
  labs(x=NULL, y=NULL) +
  scale_y_discrete(labels=dat$Category)

## Log-ratio plot using lm_sig/lm_sig2 subset 

bind_rows(log_ratio_self, .id="taxon_id") %>%
  filter(!TimePoint %in% "PreTx") %>% 
  mutate(TimePoint=fct_relevel(TimePoint, c("3M","6M","12M","24M"))) %>%
  filter(taxon_id %in% lm_sig2) %>% 
  mutate(taxon_id=factor(taxon_id, levels=lm_sig2)) %>% 
  ggplot(aes(x=TimePoint, y=clr, group=TimePoint)) + 
  geom_point(color="red") + 
  geom_line(aes(group=1), color="red") +
  geom_errorbar(aes(ymin=clr-clr_std, ymax=clr+clr_std, group=TimePoint), width=0.2) + 
  labs(y="Log fold change", x="Time since transplantation") +
  facet_wrap(.~taxon_id, ncol=6, scales="free_y") +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  labs(title="Log ( Timepoint / PreTx)") +
  theme(plot.title=element_text(size=10),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        #axis.ticks.length=unit(0.25,"cm"), 
        axis.text.x=element_text(size=8, color="black"),
        axis.text.y=element_text(size=8, color="black"),
        strip.text.x=element_text(size=6, colour="black")) 


Fig4C <- bind_rows(log_ratio_self, .id="taxon_id") %>%
  filter(!TimePoint %in% "PreTx") %>% 
  mutate(TimePoint=fct_relevel(TimePoint, c("3M","6M","12M","24M")),
         timepoint=as.numeric(TimePoint)) %>%
  filter(taxon_id %in% lm_sig2) %>% 
  mutate(taxon_id=factor(taxon_id, levels=lm_sig2)) %>% 
  group_by(taxon_id) %>% 
  mutate(majority_sign=median(sign(clr),na.rm=T),
         majority_sign2=recode_factor(majority_sign, `1`="pos",`-1`="neg"),
         taxon_id=fct_reorder(taxon_id, majority_sign, median),
         taxon_id_24M=case_when(TimePoint=="24M" ~ taxon_id)) %>%
  ungroup() %>% 
  
  ggplot(aes(x=timepoint, y=clr, group=taxon_id, color=majority_sign2, label=taxon_id_24M)) + 
  scale_color_manual(values=c("pos"="blue","neg"="red")) +
  geom_point() + 
  geom_line() +
  geom_label_repel(hjust=-1,
                   vjust=0,
                   nudge_x=2,
                   direction="y",
                   segment.size=0.2,
                   size=4,
                   segment.colour = "black",
                   fontface = 'bold', color = 'black') +
  #geom_errorbar(aes(ymin=clr-clr_std, ymax=clr+clr_std, group=TimePoint), width=0.2) + 
  labs(y="Log fold change", x="Time since transplantation") +
  geom_segment(aes(x=1,xend=4,y=0,yend=0), linetype="dashed", color="black") +
  labs(title="Log ( Timepoint / PreTx)") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"))+
  theme(legend.text = element_text(face = "bold", size=12),
        plot.title=element_text(size=14, face="bold", colour = "black", hjust=0.5),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1))+
  scale_x_continuous(expand=c(0.05, 0, 0, 1), breaks=seq(from=1, to=4, by=1), labels=c("3M","6M","12M","24M"))
Fig4C

Fig4C_nolab <- bind_rows(log_ratio_self, .id="taxon_id") %>%
  filter(!TimePoint %in% "PreTx") %>% 
  mutate(TimePoint=fct_relevel(TimePoint, c("3M","6M","12M","24M")),
         timepoint=as.numeric(TimePoint)) %>%
  filter(taxon_id %in% lm_sig2) %>% 
  mutate(taxon_id=factor(taxon_id, levels=lm_sig2)) %>% 
  group_by(taxon_id) %>% 
  mutate(majority_sign=median(sign(clr),na.rm=T),
         majority_sign2=recode_factor(majority_sign, `1`="pos",`-1`="neg"),
         taxon_id=fct_reorder(taxon_id, majority_sign, median),
         taxon_id_24M=case_when(TimePoint=="24M" ~ taxon_id)) %>%
  ungroup() %>% 
  
  ggplot(aes(x=timepoint, y=clr, group=taxon_id, color=majority_sign2, label=taxon_id_24M)) + 
  scale_color_manual(values=c("pos"="blue","neg"="red")) +
  geom_point() + 
  geom_line() +
  #geom_errorbar(aes(ymin=clr-clr_std, ymax=clr+clr_std, group=TimePoint), width=0.2) + 
  labs(y="Log fold change", x="Time since transplantation") +
  labs(title="Log ( Timepoint / PreTx)") +
  theme(
    axis.title.x = element_text(size=14, face="bold", colour = "black"),
    axis.title.y = element_text(size=14, face="bold", colour = "black"),
    axis.text.x = element_text(size=14, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, face="bold", colour = "black"))+
  theme(legend.text = element_text(face = "bold", size=12),
        plot.title=element_text(size=14, face="bold", colour = "black", hjust=0.5),
        legend.position="none",
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1))+
  scale_x_continuous(expand=c(0.05, 0, 0, 1), breaks=seq(from=1, to=4, by=1), labels=c("3M","6M","12M","24M"))
Fig4C_nolab

setwd("I:/Metagenomic sequencing/Analyzes/TxN_microbobiome analysis V1/Tx_N_merged_data/Nature Medicine Plots/Figure 4")
ggsave("Fig4C.png", device = "png", type = "cairo", dpi=300, w=16, h=10)

#5.7      Fig 6D Pathway Heatmap----
#Pathwys
#ftbl_rclr <- 
#  read.csv("I:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDtaxa_All_rclr.csv", row.names=1)

#CDtaxa_All_only_s <-  CDtaxa_All[,grep(x=colnames(CDtaxa_All),pattern="s__")]
#row.names(CDtaxa_All_only_s)
#rownames(ftbl_rclr) <- str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7)[,7] #colnames(CDtaxa_All_only_s)
#colnames(ftbl_rclr) <- rownames(CDtaxa_All_only_s)

ftbl_rclr <- 
  read.csv("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\CDpwys_All_rclr.csv", row.names=1)

# pathways
rownames(ftbl_rclr) <- gsubfn(".", list("'"="","-"=".",":"="."," "=".","("=".",")"=".",";"=".",","="."),colnames(CDpwys_All))

#rownames(ftbl_rclr) <- colnames(CDpwys_All)
colnames(ftbl_rclr) <- rownames(CDpwys_All)

meta <- CDmeta_All %>% 
  select(ID, Timepoint, Patient, exclude_liver_long) %>% 
  filter(exclude_liver_long %in% c("3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN")) %>% 
  mutate(TimePoint=recode_factor(Timepoint,"PreTxRTR"="PreTx","3months"="3M","6months"="6M","12months"="12M","24months"="24M")) %>% 
  na.omit()

# subset to renal longi samples
ftbl_rclr <- ftbl_rclr[,colnames(ftbl_rclr) %in% meta$ID]
ftbl_rclr <- t(ftbl_rclr)

tmp_mat <- diag(length(colnames(ftbl_rclr)))
rownames(tmp_mat) <- colnames(tmp_mat) <- colnames(ftbl_rclr)
pairs <- reshape2::melt(tmp_mat) %>% 
  filter(Var1==Var2) %>% 
  select(numerator=Var1, denominator=Var2) %>% 
  arrange(numerator, denominator) %>% 
  mutate(numerator=as.character(numerator), denominator=as.character(denominator))

pairs$TimePoint <- NA
log_ratio_samples <- vector("list", nrow(ftbl_rclr))
names(log_ratio_samples) <- rownames(ftbl_rclr)
log_ratio_samples <- lapply(log_ratio_samples, function(x) pairs)

for(sample in names(log_ratio_samples)) {
  print(sample)
  for(pair in 1:nrow(log_ratio_samples[[sample]])) {
    log_ratio_samples[[sample]][pair,]$TimePoint <- as.character(meta[meta$ID==sample,]$TimePoint)
  }
}

pairs_across_samples <- bind_rows(log_ratio_samples, .id="ID")

log_ratio_self <- vector("list", ncol(ftbl_rclr))
names(log_ratio_self) <- colnames(ftbl_rclr)
log_ratio_self <- lapply(log_ratio_self, function(x) data.frame(clr=rep(NA,(length(unique(meta$TimePoint))-1)),
                                                                clr_std=rep(NA,(length(unique(meta$TimePoint))-1)),
                                                                TimePoint=c("3M","6M","12M","24M")))
for(taxon in colnames(ftbl_rclr)) {
  print(taxon)
  for(v in c("3M","6M","12M","24M")) {
    #log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$log_ratio <- log2( EnvStats::geoMean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon]) / EnvStats::geoMean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon]) )
    #log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$log_ratio_sd <- log( EnvStats::geoSD(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon]) / EnvStats::geoSD(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon]) )
    log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$clr <- mean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon], na.rm=T) - mean(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon], na.rm=T)
    log_ratio_self[[taxon]][log_ratio_self[[taxon]]$TimePoint==v,]$clr_std <- sd(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint==v,]$ID, taxon], na.rm=T) - sd(ftbl_rclr[pairs_across_samples[pairs_across_samples$numerator==taxon & pairs_across_samples$TimePoint=="PreTx",]$ID, taxon], na.rm=T)
  }
}

## Taxon-specific linear models

lm_data <- data.frame(ftbl_rclr, check.names=F, stringsAsFactors=F) %>%  # for pathways
  #lm_data <- data.frame(ftbl_rclr) %>% 
  rownames_to_column("ID") %>%
  left_join(meta, by="ID")

lm_ls <- vector("list", ncol(ftbl_rclr))
names(lm_ls) <- colnames(ftbl_rclr)
for(taxon in colnames(ftbl_rclr)) {
  tryCatch({
    print(taxon)
    m <- lmerTest::lmer(lm_data[,taxon]~TimePoint+(1|Patient), data=lm_data[,c(taxon,"TimePoint","Patient")])
    
    mm <- data.frame(coef(summary(m)))
    colnames(mm)[5] <- "p_val"
    #mm$p_adj <- p.adjust(mm$p_val, "BH", ncol(ftbl_rclr))
    lm_ls[[taxon]] <- mm
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
lm_ls[sapply(lm_ls, is.null)]
# remove NULL lists - these are generated by to few samples: ERROR : number of levels of each grouping factor must be < number of observations (problems: Patient) 
# in case of pathways, there are 4
lm_ls[sapply(lm_ls, is.null)] <- NULL

lm_df <- lm_ls %>% 
  bind_rows(.id="taxon_id") %>% 
  rownames_to_column("group") %>% 
  separate(group, into=c("group",NA)) %>% 
  filter(group!="") %>% 
  mutate(p_adj=p.adjust(p_val, "BH"))

lm_df_sig <- lm_df %>% 
  filter(p_adj<0.1)

lm_sig <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  mutate(taxon_id=fct_reorder(taxon_id, Estimate, median)) %>% 
  select(taxon_id) %>%
  pull() %>% 
  levels() %>% 
  rev()

#filter out taxa with average effect size of below 0.1
lm_sig2 <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  group_by(taxon_id) %>% 
  summarise(mean_es=mean(Estimate)) %>% 
  filter(mean_es >= 0.1 | mean_es <= -0.1) %>% 
  arrange(desc(mean_es)) %>% 
  select(taxon_id) %>%
  pull()

#species_taxonomy <- data.frame(str_split_fixed(colnames(CDtaxa_All_only_s), "\\.", 7))[,c(2,5,7)]
#colnames(species_taxonomy) <- c("Phylum","Family","taxon_id")

pathway_classes <- data.frame(readxl::read_excel("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\path_annotated_index_edited_to_match_ourdata.xlsx")) %>% 
  select(taxon_id=MetaCyc_complete_name, Class, Category) %>% 
  mutate(taxon_id=gsubfn(".", list("'"="","-"=".",":"="."," "=".","("=".",")"=".",";"=".",","="."), taxon_id))

dat <- lm_df %>% 
  filter(taxon_id %in% lm_df_sig$taxon_id) %>% 
  
  left_join(pathway_classes, by="taxon_id") %>% 
  #left_join(species_taxonomy, by="taxon_id") %>% 
  
  
  mutate(sign=sign(Estimate),
         sign=ifelse(sign=="1","+","-"),
         group=str_split(group, "TimePoint", simplify=T)[,2],
         #group=fct_relevel(group, c("3M","6M","12M","24M")),
         group=fct_relevel(group, c("24M","12M","6M","3M")),
         taxon_id=fct_reorder(taxon_id, Estimate, median)) %>%
  filter(taxon_id %in% lm_sig2) 

#Fill in missing classes due to fail in merge
write.xlsx(dat, "Data_for_pathway_heatmap_Final.xlsx", row.names=F)

dataForPWYSHeatmap  <- read.xlsx("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Longitudinal analysis 13-12-2020\\codedata\\Data_for_pathway_heatmap_Final_filledin_missing_Class.xlsx",1)
dataForPWYSHeatmap$Class <- as.character(dataForPWYSHeatmap$Class)
dataForPWYSHeatmap$Category <- as.character(dataForPWYSHeatmap$Category)
summary(dataForPWYSHeatmap$Class)
unique(dataForPWYSHeatmap$Class)

dataForPWYSHeatmap$Class <- as.factor(dataForPWYSHeatmap$Class)
summary(dataForPWYSHeatmap$Class)

dataForPWYSHeatmap2 <-  dataForPWYSHeatmap %>% 
  mutate(group=fct_relevel(group, c("3M","6M","12M","24M")),
         taxon_id=fct_reorder(taxon_id, Estimate, median))

Fig4D <- dataForPWYSHeatmap2 %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile(colour="black", width=0.5, height=0.5) +
  scale_fill_gradient2(name="Effect size", low = "#1111AA", mid = "white", high = "#AA1111") +
  #geom_text(aes(label=sign), size=3) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust=1, angle=90),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, face="bold", colour = "black"),
        legend.title =element_text(size=14, face="bold", colour = "black")) +
  labs(x=NULL, y=NULL) +
  coord_flip()+
  scale_y_discrete(labels=setNames(dataForPWYSHeatmap2$Class, dataForPWYSHeatmap2$taxon_id))
Fig4D
ggsave("Fig4DFinal.png", device = "png", type = "cairo", dpi=300, w=12, h=8)

Fig4D_noleg <- dataForPWYSHeatmap2 %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile(colour="black", width=0.5, height=0.5) +
  scale_fill_gradient2(name="Effect size", low = "#AA1111", mid = "white", high = "#1111AA") +
  #geom_text(aes(label=sign), size=3) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust=1, angle=90),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, face="bold", colour = "black"),
        legend.title =element_text(size=14, face="bold", colour = "black"), 
        legend.position="none") +
  labs(x=NULL, y=NULL) +
  coord_flip()+
  scale_y_discrete(labels=setNames(dataForPWYSHeatmap2$Class, dataForPWYSHeatmap2$taxon_id))
Fig4D_noleg
ggsave("Fig4DNoLegend.png", device = "png", type = "cairo", dpi=300, w=12, h=8)

Fig4DTaxon_ID <- dataForPWYSHeatmap2 %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile(colour="black", width=0.5) +
  scale_fill_gradient2(name="Effect size", low = "#AA1111", mid = "white", high = "#1111AA") +
  geom_text(aes(label=sign), size=3) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust=1, angle=90),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, face="bold", colour = "black"),
        legend.title =element_text(size=14, face="bold", colour = "black")) +
  labs(x=NULL, y=NULL) +
  coord_flip()
Fig4DTaxon_ID
ggsave("Fig4D_TaxonID.png", device = "png", type = "cairo", dpi=300, w=12, h=8)

Fig4DClass <- dataForPWYSHeatmap2 %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile(colour="black", width=0.5, height=0.5) +
  scale_fill_gradient2(name="Effect size", low = "#1111AA", mid = "white", high = "#AA1111") +
  #geom_text(aes(label=sign), size=3) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust=1, angle=90),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, face="bold", colour = "black"),
        legend.title =element_text(size=14, face="bold", colour = "black")) +
  labs(x=NULL, y=NULL) +
  coord_flip()+
  scale_y_discrete(labels=setNames(dataForPWYSHeatmap2$Class, dataForPWYSHeatmap2$taxon_id))
Fig4DClass
ggsave("Fig4D_Class.png", device = "png", type = "cairo", dpi=300, w=12, h=8)


Fig4DCategory <- dataForPWYSHeatmap2 %>% ggplot(aes(y=taxon_id,x=group,fill=Estimate)) +
  #geom_tile(width=0.3) +
  geom_tile(colour="black", width=0.5) +
  scale_fill_gradient2(name="Effect size", low = "#1111AA", mid = "white", high = "#AA1111")+
  geom_text(aes(label=sign), size=3) +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_rect(fill="white"),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        axis.text.x = element_text(size=6, face="bold", colour = "black", hjust=1, angle=90),
        axis.text.y = element_text(size=12, face="bold", colour = "black"),
        legend.text = element_text(size=12, face="bold", colour = "black"),
        legend.title =element_text(size=14, face="bold", colour = "black"),
        legend.position = "none") +
  labs(x=NULL, y=NULL) +
  coord_flip()+
  scale_y_discrete(labels=setNames(dataForPWYSHeatmap2$Category, dataForPWYSHeatmap2$taxon_id))
Fig4DCategory
ggsave("Fig4D_Category.png", device = "png", type = "cairo", dpi=300, w=12, h=8)

#5.8      Combine Figure 6----
Fig4AB <- ggarrange(Fig4A, Fig4B, nrow=1, align=c("h"))
Fig4AB2 <- ggarrange(Fig4A2, Fig4B2, nrow=1, align=c("h"))
Fig4ABCombined <- ggarrange(Fig4AB, Fig4AB2, ncol=1, align=c("v"))
Fig4ABCombined
ggsave("Fig4AB_Combined_Axis.png", device = "png", type = "cairo", dpi=300, w=12, h=12)
Fig4ABC <- ggarrange(Fig4AB, Fig4C, ncol=1, align=c("v"))
Fig4ABC
ggsave("Fig4ABC_Combined.png", device = "png", type = "cairo", dpi=300, w=12, h=24)
Fig4CD <- ggarrange(Fig4C_nolab, Fig4D_noleg, ncol=1, align=c("v"))
Fig4CD

Fig4ABCD <- ggarrange(Fig4AB, Fig4CD, ncol=1, align=c("v"))
Fig4ABCD
ggsave(plot=Fig4ABCD, "Fig4ABCDCombinedV2.png", device = "png", type = "cairo", dpi=300, w=12, h=28)

ggsave(plot=Fig4ABCD, "Fig4ABCD_forD_CombinedV2.png", device = "png", type = "cairo", dpi=300, w=12, h=12)

Fig4CDCategory <- ggarrange(Fig4C_nolab, Fig4DCategory, ncol=1, align=c("v"))
Fig4ABCDCategory <- ggarrange(Fig4AB, Fig4CDCategory, ncol=1, align=c("v"))
ggsave(plot=Fig4ABCDCategory, "Fig4ABCD_forD_CombinedCategory.png", device = "png", type = "cairo", dpi=300, w=12, h=12)

#5.2      Load data----
#Nature Medicine Longitudinal Main Figure
CDtaxa <- read.table("I:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_metaphlan_cleaned_filtered.txt",sep='\t',header=T, fill=T)
CDtaxa$ID <- as.character(CDtaxa$ID)
trimws(CDtaxa$ID)
rownames(CDtaxa) <- CDtaxa$ID

CDpwys <- read.table("I:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\TxN_merged_humann_pathabundance_filtered.txt",sep='\t',header=T)
CDpwys$ID <- as.character(CDpwys$ID)
trimws(CDpwys$ID)
rownames(CDpwys) <- CDpwys$ID

CDmeta = read.spss("J:\\Metagenomic sequencing\\Analyzes\\TxN_microbobiome analysis V1\\Tx_N_merged_data\\Sending Data To Johannes\\Metadata.sav", to.data.frame=TRUE)
CDmeta$ID <- as.character(CDmeta$ID)
#SPSS String variable for ID puts in stupid white spaces :)
CDmeta$ID <- trimws(CDmeta$ID)
rownames(CDmeta) <- CDmeta$ID

summary(CDmeta$Selection_Long_Cross_S)
CDmeta_All <- subset(CDmeta, Selection_Long_Cross_S == "Yes")
#Taxa
CDtaxa_All <- CDmeta_All["ID"]
CDtaxa_All <- merge(CDtaxa_All, CDtaxa, by="ID")
row.names(CDtaxa_All) <- CDtaxa_All$ID

#Pathways
CDpwys_All <- CDmeta_All["ID"]
CDpwys_All <- merge(CDpwys_All, CDpwys, by="ID")
row.names(CDpwys_All) <- CDpwys_All$ID
CDpwys_All$ID <- NULL

CDmeta_All <- CDtaxa_All["ID"]
CDmeta_All <- merge(CDmeta_All, CDmeta, by="ID")

CDtaxa_All$ID <- NULL

CDmeta_All <- CDmeta_All %>% 
  mutate(Cross_sectional_samplesDag3=factor(Cross_sectional_samplesDag3, levels=c("Healthy Control","LTR","RTR","")),
         exclude_liver_long=paste0(Longitudinal_Timepoint, "_", TypeTxN_TxL_TxD_DAG3))
