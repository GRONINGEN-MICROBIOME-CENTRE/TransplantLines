This github repository provides the code to reproduce the analyses presented in "**Gut Microbiome Dysbiosis is Associated with Increased Mortality following Solid Organ Transplantation**"

By: Weersma Group, Fu Group and Zhernakova Group, UMCG, The Netherlands

Last updated 19-07-2021

----------------------------------------------------------------------------------------------------------------------
**File index:**
* mortality.Rmd: R code to reproduce the mortality analysis (Fig. 2C-2H) 
* differenctial_abundance.Rmd: R code to reproduce the differential abundance analysis (Tables S1-S3)
* variance_explained_phenotypes.Rmd: R code to reproduce the PERMANOVA analysis (Fig. 4)
* immunosuppresive_drugs.Rmd: R code to reproduce the immunosuppresive drugs analysis (Fig. 5)
* longitudinal.Rmd: R code reproduce the longitudinal analysis (Fig. 6)
* DEICODE.py: Python code to run the stand-alone version of DEICODE (https://github.com/biocore/DEICODE)
* GraPhlAn.py: Python code to run GraPhlAn (https://github.com/biobakery/graphlan) 

----------------------------------------------------------------------------------------------------------------------
**Dirs:**
* deicode_output: Sample and feature loadings for species and pathways produced by DEICODE 
* mock_data/deicode_input: Input files (mock abundance tables) to run DEICODE
* mock_data/deicode_output: Output files (sample/feature loadings and rclr-transformed abundance tables) produced by DEICODE on mocked data
* r_scripts_library: Collections of various R functions
----------------------------------------------------------------------------------------------------------------------
**Mock Data:**

As the participants' informed consent and privacy regulations for the Lifelines biobank and Transplantlines project (https://umcgresearch.org/en-GB/web/research/w/transplantlines) prevent public sharing of metadata, the majority of analyses in this repository have been implemented with *mock data*. These data include a selection of anonymized and modified phenotypes including microbiome features. Becuase of this, the majority of code will not produce identical results as in the focal manuscript. Exceptions include the robust compositional PCA plots (Figs. 2A; 6A and S7) where we have shared the sample and feature loadings (for species and pathways) produced by DEICODE.     

Data used in the Transplantlines project can be requested from TransplantLines: datarequest.transplantlines@umcg.nl
