# ========================================================================================================
#       Immunosuppressive medication analysis
#
# NOTE: This is implementation on Mock Data published in this github repo,
#       and is intended for demonstration purposes only, the
#       results are not identical to Figures in the manuscript
# ========================================================================================================

library(foreach)

mock_taxa=read.table("mock data/Mocks_CDtaxa_All_rclr.txt",check.names = F,stringsAsFactors = F,header = T)
mock_path=read.table("mock data/Mocks_CDpwys_All_rclr.txt",check.names = F,stringsAsFactors = F,header = T)
tmp_covariate=read.table("mock data/MOCK_Covariate.txt",header = T,stringsAsFactors = F,check.names = F,row.names = 1)
tmp_taxa=as.data.frame(t(mock_taxa))
tmp_pwys=as.data.frame(t(mock_path))
tmp_taxa=tmp_taxa[rownames(tmp_taxa) %in% rownames(mock_covariate),]
tmp_pwys=tmp_pwys[rownames(tmp_pwys) %in% rownames(mock_covariate),]

# test medication individually
taxa_result_pred = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              dum_pred , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_tacrolimus = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              tacrolimus , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_mycofenolzuur = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              mycofenolzuur, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_ciclosporine = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              ciclosporine, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_azathioprine = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              azathioprine, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_case_control = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

pathway_result_pred = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              dum_pred, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_tacrolimus = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              tacrolimus, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_mycofenolzuur = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              mycofenolzuur, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_ciclosporine = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              ciclosporine, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_azathioprine = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              azathioprine, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_case_control = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by="row.names",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics, data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

taxa_result_pred$FDR=p.adjust(taxa_result_pred$Pvalue,method = "BH")
taxa_result_pred=taxa_result_pred[taxa_result_pred$Factor=="dum_predyes" & taxa_result_pred$FDR<0.1,]
taxa_result_tacrolimus$FDR=p.adjust(taxa_result_tacrolimus$Pvalue,method = "BH")
taxa_result_tacrolimus=taxa_result_tacrolimus[taxa_result_tacrolimus$Factor=="tacrolimusyes" & taxa_result_tacrolimus$FDR<0.1,]
taxa_result_mycofenolzuur$FDR=p.adjust(taxa_result_mycofenolzuur$Pvalue,method = "BH")
taxa_result_mycofenolzuur=taxa_result_mycofenolzuur[taxa_result_mycofenolzuur$Factor=="mycofenolzuuryes" & taxa_result_mycofenolzuur$FDR<0.1,]
taxa_result_ciclosporine$FDR=p.adjust(taxa_result_ciclosporine$Pvalue,method = "BH")
taxa_result_ciclosporine=taxa_result_ciclosporine[taxa_result_ciclosporine$Factor=="ciclosporineyes" & taxa_result_ciclosporine$FDR<0.1,]
taxa_result_azathioprine$FDR=p.adjust(taxa_result_azathioprine$Pvalue,method = "BH")
taxa_result_azathioprine=taxa_result_azathioprine[taxa_result_azathioprine$Factor=="azathioprineyes" & taxa_result_azathioprine$FDR<0.1,]
taxa_result_case_control$FDR=p.adjust(taxa_result_case_control$Pvalue,method = "BH")
taxa_result_case_control=taxa_result_case_control[taxa_result_case_control$Factor=="TypeMoreThan12MonthsYes" & taxa_result_case_control$FDR<0.1,]

pathway_result_pred$FDR=p.adjust(pathway_result_pred$Pvalue,method = "BH")
pathway_result_pred=pathway_result_pred[pathway_result_pred$Factor=="dum_predyes" & pathway_result_pred$FDR<0.1,]
pathway_result_tacrolimus$FDR=p.adjust(pathway_result_tacrolimus$Pvalue,method = "BH")
pathway_result_tacrolimus=pathway_result_tacrolimus[pathway_result_tacrolimus$Factor=="tacrolimusyes" & pathway_result_tacrolimus$FDR<0.1,]
pathway_result_mycofenolzuur$FDR=p.adjust(pathway_result_mycofenolzuur$Pvalue,method = "BH")
pathway_result_mycofenolzuur=pathway_result_mycofenolzuur[pathway_result_mycofenolzuur$Factor=="mycofenolzuuryes" & pathway_result_mycofenolzuur$FDR<0.1,]
pathway_result_ciclosporine$FDR=p.adjust(pathway_result_ciclosporine$Pvalue,method = "BH")
pathway_result_ciclosporine=pathway_result_ciclosporine[pathway_result_ciclosporine$Factor=="ciclosporineyes" & pathway_result_ciclosporine$FDR<0.1,]
pathway_result_azathioprine$FDR=p.adjust(pathway_result_azathioprine$Pvalue,method = "BH")
pathway_result_azathioprine=pathway_result_azathioprine[pathway_result_azathioprine$Factor=="azathioprineyes" & pathway_result_azathioprine$FDR<0.1,]
pathway_result_case_control$FDR=p.adjust(pathway_result_case_control$Pvalue,method = "BH")
pathway_result_case_control=pathway_result_case_control[pathway_result_case_control$Factor=="TypeMoreThan12MonthsYes" & pathway_result_case_control$FDR<0.1,]


# test medication combinations
tmp_covariate$group1=NA
tmp_covariate$group1[tmp_covariate$mycofenolzuur=="yes" & tmp_covariate$tacrolimus=="yes" & tmp_covariate$dum_pred=="yes"]="yes"
tmp_covariate$group1[is.na(tmp_covariate$group1)]="no"

tmp_covariate$group2=NA
tmp_covariate$group2[tmp_covariate$mycofenolzuur=="yes" & tmp_covariate$tacrolimus=="yes"]="yes"
tmp_covariate$group2[is.na(tmp_covariate$group2)]="no"

tmp_covariate$group3=NA
tmp_covariate$group3[tmp_covariate$tacrolimus=="yes" & tmp_covariate$dum_pred=="yes"]="yes"
tmp_covariate$group3[is.na(tmp_covariate$group3)]="no"

tmp_covariate$group4=NA
tmp_covariate$group4[tmp_covariate$mycofenolzuur=="yes"  & tmp_covariate$dum_pred=="yes"]="yes"
tmp_covariate$group4[is.na(tmp_covariate$group4)]="no"

tmp_covariate$group5=NA
tmp_covariate$group5[tmp_covariate$ciclosporine =="yes"  & tmp_covariate$dum_pred=="yes"]="yes"
tmp_covariate$group5[is.na(tmp_covariate$group5)]="no"

taxa_result_group1 = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group1 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_group2 = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group2 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_group3 = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group3 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_group4 = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group4 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_group5 = foreach(i=1:ncol(tmp_taxa),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_taxa)[i]
    print(taxon)
    lm_data=merge(tmp_taxa[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group5 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

pathway_result_group1 = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group1 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_group2 = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group2 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_group3 = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group3 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_group4 = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group4 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
pathway_result_group5 = foreach(i=1:ncol(tmp_pwys),.combine = rbind) %do%  {
  tryCatch({
    taxon=colnames(tmp_pwys)[i]
    print(taxon)
    lm_data=merge(tmp_pwys[,taxon,drop=F],tmp_covariate,by.x="row.names",by.y="ID",all=F)
    m <- lm(lm_data[,taxon] ~ TypeMoreThan12Months + ALGLEEFT + BMI + GESLACHT + trans_roke_1 + dum_ppi + dum_laxative + dum_antibiotics +
              group5 , data=lm_data)
    mm <- as.data.frame(coef(summary(m)))
    mm=mm[rownames(mm)!="(Intercept)",]
    
    return.string=data.frame(Factor=rownames(mm),Estimate=mm$Estimate,SD=mm$`Std. Error`,Pvalue=mm$`Pr(>|t|)`,Taxa=taxon)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
taxa_result_group1$FDR=p.adjust(taxa_result_group1$Pvalue,method = "BH")
taxa_result_group2$FDR=p.adjust(taxa_result_group2$Pvalue,method = "BH")
taxa_result_group3$FDR=p.adjust(taxa_result_group3$Pvalue,method = "BH")
taxa_result_group4$FDR=p.adjust(taxa_result_group4$Pvalue,method = "BH")
taxa_result_group5$FDR=p.adjust(taxa_result_group5$Pvalue,method = "BH")

pathway_result_group1$FDR=p.adjust(pathway_result_group1$Pvalue,method = "BH")
pathway_result_group2$FDR=p.adjust(pathway_result_group2$Pvalue,method = "BH")
pathway_result_group3$FDR=p.adjust(pathway_result_group3$Pvalue,method = "BH")
pathway_result_group4$FDR=p.adjust(pathway_result_group4$Pvalue,method = "BH")
pathway_result_group5$FDR=p.adjust(pathway_result_group5$Pvalue,method = "BH")









