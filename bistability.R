

library(microbiome)

# filter data

ps_all <- phyloseq(otu_table(CDtaxa_All_only_s, taxa_are_rows=F),
                   sample_data(CDmeta_All),
                   tax_table(CDtaxatbl_All)) %>% 
  aggregate_rare(level = "species", detection = .1/100, prevalence = 10/100)


ps_hc <- prune_samples(samples=CDmeta_All[CDmeta_All$exclude_liver_long %in% "Healthy_Control_DAG3 Control",]$ID, x=ps_all)

ps_hc_rtr <- prune_samples(samples=CDmeta_All[CDmeta_All$exclude_liver_long %in% c("RTR_CS_TxN","3months_TxN","6months_TxN","12months_TxN","24months_TxN","PreTxRTR_TxN","Healthy_Control_DAG3 Control"),]$ID, x=ps_all)

ps_rtr_long <- prune_samples(samples=CDmeta_All[CDmeta_All$exclude_liver_long %in% c("PreTxRTR_TxN","3months_TxN","6months_TxN","12months_TxN","24months_TxN"),]$ID, x=ps_all)
