By: Weersma Group, Fu Group and Zhernakova Group, UMCG, The Netherlands
# TransplantLines
This github repo describes workflow and codes used in The TransplantLines study:

Contents:
* Microbiome profiling and description of microbiome
* Mortality Analysis
* Differential abundance analysis
* Calculation of microbiome variance explained by phenotypes
* Immunosuppresive medication analysis (Immunosuppresive medication analysis and Figure 5.Rmd)
* Longitudinal analysis (Longitudinal analysis and Figure 6.Rmd)

***Microbiome profiling and description of microbiome***

Illumina adapters and low-quality reads (Phred score <30) were filtered out using KneadData (v0.5.1). Then Bowtie2 (v2.3.4.1) was used to remove reads aligned to the human genome (hg19). The quality of the reads was examined using FastQC toolkit (v0.11.7). Taxonomy alignment was done by MetaPhlAn2 (v2.7.2) against the database of marker genes mpa_v20_m200. Metacyc pathways were profiled by HUMAnN2 (v0.11.1). Bacterial virulence factors and antibiotic resistance genes were identified using shortBRED (shortbred_identify.py (v0.9.5) and shortbred_quantify.py tool (v0.9.5)) against virulence factors of pathogenic bacteria (VFDB) database (http://www.mgc.ac.cn/VFs/main.htm) and comprehensive antibiotic resistance database (CARD) (https://card.mcmaster.ca/) separately. Samples were further excluded by criteria that eukaryotic or viral abundance > 25% of total microbiome content or total read depth < 10 million. In total, we identified 1132 taxa (17 phyla, 27 class, 52 order, 98 family, 231 genera and 705 species), 586 metabolic pathways, 313 virulence factors and 957 antibiotic resistance genes. With a presence of 10% and relative abundance threshold of 0.01%, 384 taxa (8 phyla, 14 class, 20 order, 40 family, 83 genera and 219 species), 351 metabolic pathways, 323 virulence factors and 167 antibiotic resistance genes were left after filtering. Total-sum normalization was applied to all microbiome data after filtering. Analyses were performed using locally installed tools and databases on CentOS (release 6.9) on the high-performance computing infrastructure available at UMCG and University of Groningen (RUG).

***Mortality Analysis***

We assessed the association between gut microbial diversity (Shannon diversity index) and all-cause mortality using multivariable Cox proportional hazard models adjusting for patient age, sex and years since transplantation. To further minimize the effect of the immortal time bias, we used each patient’s first sample collected in 2014 as the start of the follow-up. To reduce surgery related mortality, we removed samples collected <1 year post-transplantation.  We analyzed both the association between all-cause mortality and (1) gut microbial diversity and (2) the gut community dissimilarity (i.e. the Aitchison distance) between transplant recipients and healthy controls. For diversity, patients were stratified into high-diversity and low-diversity groups based on the median of Shannon diversity, and the relationship of each group with patient mortality was determined using Kaplan-Meier curves for LTR and RTR separately.The continuous variables (i.e. Shannon diversity and the Aitchison distance to healthy controls) were scaled to unit variance prior to the analysis. Analyses of all-cause mortality were presented as hazard ratios with 95% CIs. Finally, we checked whether each included variable independently satisfied the assumptions of the Cox model by computing the Schoenfeld residuals (P>0.05). For the survival analysis we used the survival and the rms R package.

***Differential abundance analysis***

We performed differential abundance analysis in LTR vs HC., RTR vs HC, ESLD vs HC, ESRD vs HC and user vs non-users of immunosuppressive drugs. Furthermore, we performed longitudinal analysis in longitudinal data from renal trans plant recipients where we compared post-transplantation samples to pre-transplantation samples. 

For taxa and metabolic pathways we modeled rclr-transformed relative abundances with the lm function in R using the following linear model: 

Feature ~ Cross_sectional_samples + Age  +  BMI + Gender  + Smoking + PPI + Laxatives  + Antitibiotics

Antibiotic resistance genes and virulence factors were extremely sparse in healthy controls. Therefore, we could not apply the same linear models that we used for microbial species and pathways. Instead, we performed logistic regression models on presence-absences of antibiotic resistance genes or virulence factors with a prevalence cutoff at 1%. We used a generalized linear model (family = binomial) from the glm function in R: 

Feature ~ Cross_sectional_samples + Age  +  BMI + Gender  + Smoking + PPI + Laxatives  + Antitibiotics

***Calculation of microbiome variance explained by phenotypes***

The microbiome composition variance explained by phenotypes was calculated by permutational multivariate analysis of variance using distance matrices, implemented in the adonis function for R package vegan (v.2.4-6), using 9999 permutations and a Bray-Curtis distance matrix calculated using relative abundances of microbial species.

***Immunosuppresive medication analysis***

We tested five immunosuppressive drugs individually between users and non-users, correcting for surgery, age, sex, BMI, smoking, and the use of proton pump inhibitor (PPI), laxatives and antibiotics using linear models. 
The microbial data contains the rclr-transformed relative abundances of species and pathways. 
The medication contains prednisolone; tacrolimus; mycophenolic acid; azathioprine; and cyclosporin, and five drugs in combination includes cyclosporin and prednisolone; mycophenolic acid, tacrolimus and prednisolone; mycophenolic acid and tacrolimus; tacrolimus and prednisolone; mycophenolic acid and prednisolone.

***Longitudinal analysis***

Using our longitudinal data, we constructed log-ratios comparing each microbial species’ rclr-transformed relative abundance at a time point post-transplantation to its rclr-transformed relative abundance at pre-transplantation. These log-ratios were constructed by subtracting each focal species’ average rclr-transformed relative abundances across all post transplantation samples from it’s averaged rclr-transformed relative abundances across all pre transplantation samples.
