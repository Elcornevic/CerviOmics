# BIOINFORMATIC ANALYSIS WITH 'THE CANCER GENOME ATLAS' (TCGA) RNA SEQUENCING DATA IN CERVIX CANCER
## Author: Victor Guillermo Cornejo Villanueva
## First tutorial: https://www.youtube.com/watch?v=UWXv9dUpxNE
setwd("~/Escritorio/R/TCGA_2")

## Install TCGAbiolinks (from GitHub)
### First, uncheck the option 'Tools --> Global options... --> Packages --> Use secure download method for HTTP' https://community.rstudio.com/t/install-packages-unable-to-access-index-for-repository-try-disabling-secure-download-method-for-http/16578
options(timeout = max(1000000, getOption("timeout"))) #https://stackoverflow.com/questions/27583767/warning-downloaded-length-reported-length-in-installing-packages-from-cran
BiocManager::install("TCGAbiolinks")

## Install another packages
BiocManager::install("maftools")
install.packages("pheatmap")
BiocManager::install("SummarizedExperiment")

library(TCGAbiolinks)
library(tidyverse)
library(maftools) ## https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html
library(pheatmap)
library(SummarizedExperiment)
library(writexl)

## Get a list of projects
gdcprojects <- getGDCprojects() ## all projects in TCGA
write_xlsx(gdcprojects, 'GDC_All_Projects_List.xlsx') ## https://datatofish.com/export-dataframe-to-excel-in-r/ 

getProjectSummary('TCGA-BRCA') ## Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (4 types of diseases)
### Reviewed on June 15
### File count = 14 982  //  Case count = 307
### Data categories 
#### Simple nucleotide variation (4 856 files & 306 cases)
#### Transcriptome profiling (1 242 files & 307 cases)
#### DNA Methylation (936 files & 307 cases)
#### Clinical (632 files & 307 cases)

## Building a query (example)
query_TCGA_clinical <- GDCquery(project = 'TCGA-CESC',
                                data.category = 'Clinical',
                                data.type = 'Clinical Supplement')
output_query_TCGA_clinical <- getResults(query_TCGA_clinical)
write_xlsx(output_query_TCGA_clinical, 'Clinical_Data_TCGA-CESC.xlsx')


## Build a query to retrieve gene expression data 
query_TCGA_transcriptome <- GDCquery(project = 'TCGA-CESC',
                                     data.category = 'Transcriptome Profiling',
                                     experimental.strategy = 'RNA-Seq',         ## Or miRNA-Seq
                                     workflow.type = 'STAR - Counts',           ## Or BCGSC miRNA Profiling (.txt), the other is .tsv files
                                     access = 'open')                           ### Includes both data type (Gene Expression Quantification and Splice Junction Quantification)
output_query_TCGA_transcriptome <- getResults(query_TCGA_transcriptome)
#transcriptome_all <- GDCprepare(query_TCGA_transcriptome, summarizedExperiment = TRUE)
write_xlsx(output_query_TCGA_transcriptome, 'TGA_Transcriptome_Profiling_List.xlsx')

### Download data - GDCdownload 
GDCdownload(query_TCGA_transcriptome) ## 309 files
GDCdownload(query_TCGA_clinical)


### Download and visualize mutation data from TCGA-Squamous cell neoplasms
setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms")
query_mutation_squamous <- GDCquery(project = 'TCGA-CESC',
                                    data.category = 'Simple Nucleotide Variation', 
                                    access = 'open',
                                    barcode = c('TCGA-VS-A94Y-01A-11D-A387-09,TCGA-VS-A94Y-10A-01D-A38A-09', 'TCGA-VS-A9V3-01A-11D-A42O-09,TCGA-VS-A9V3-10A-01D-A42R-09', 'TCGA-UC-A7PD-01A-11D-A351-09,TCGA-UC-A7PD-11A-12D-A351-09', 'TCGA-MU-A8JM-01A-11D-A36J-09,TCGA-MU-A8JM-10A-01D-A36M-09', 'TCGA-EA-A5O9-01A-11D-A28B-09,TCGA-EA-A5O9-10A-01D-A28E-09', 'TCGA-VS-A8QF-01A-21D-A37N-09,TCGA-VS-A8QF-10A-01D-A37N-09', 'TCGA-C5-A1M8-01A-21D-A13W-08,TCGA-C5-A1M8-10A-01D-A13W-08', 'TCGA-VS-AA62-01A-11D-A42O-09,TCGA-VS-AA62-10A-01D-A42R-09', 'TCGA-C5-A1BF-01B-11D-A13W-08,TCGA-C5-A1BF-10A-01D-A13W-08', 'TCGA-ZJ-AAX4-01A-11D-A42O-09,TCGA-ZJ-AAX4-10B-01D-A42R-09', 'TCGA-UC-A7PF-01A-11D-A351-09,TCGA-UC-A7PF-11A-31D-A351-09', 'TCGA-C5-A0TN-01A-21D-A14W-08,TCGA-C5-A0TN-10B-01D-A14W-08', 'TCGA-DG-A2KL-01A-11D-A17W-09,TCGA-DG-A2KL-10A-01D-A17W-09', 'TCGA-DS-A1OB-01A-11D-A16Y-08,TCGA-DS-A1OB-10A-01D-A14W-08', 'TCGA-C5-A1BE-01B-11D-A13W-08,TCGA-C5-A1BE-10A-01D-A13W-08', 'TCGA-C5-A7CK-01A-11D-A32I-09,TCGA-C5-A7CK-10A-01D-A32I-09', 'TCGA-DS-A7WI-01A-12D-A351-09,TCGA-DS-A7WI-10A-01D-A351-09',
                                                'TCGA-VS-A94X-01A-11D-A387-09,TCGA-VS-A94X-10A-01D-A38A-09', 'TCGA-C5-A8XH-01A-11D-A37N-09,TCGA-C5-A8XH-10A-01D-A37N-09', 'TCGA-VS-A8QM-01A-11D-A37N-09,TCGA-VS-A8QM-10A-01D-A37N-09', 'TCGA-VS-A8EJ-01A-11D-A37N-09,TCGA-VS-A8EJ-10A-01D-A37N-09', 'TCGA-VS-A8QA-01A-11D-A37N-09,TCGA-VS-A8QA-10A-01D-A37N-09', 'TCGA-C5-A1MN-01A-11D-A14W-08,TCGA-C5-A1MN-10A-01D-A14W-08', 'TCGA-MA-AA3Y-01A-11D-A387-09,TCGA-MA-AA3Y-10A-01D-A38A-09', 'TCGA-C5-A1ML-01A-11D-A14W-08,TCGA-C5-A1ML-10A-01D-A14W-08', 'TCGA-C5-A1MP-01A-11D-A14W-08,TCGA-C5-A1MP-10A-01D-A14W-08', 'TCGA-EA-A5ZF-01A-11D-A28B-09,TCGA-EA-A5ZF-10A-01D-A28E-09', 'TCGA-C5-A7X3-01A-11D-A351-09,TCGA-C5-A7X3-10A-01D-A351-09', 'TCGA-EA-A439-01A-11D-A243-09,TCGA-EA-A439-10A-01D-A243-09', 'TCGA-MU-A5YI-01A-11D-A32I-09,TCGA-MU-A5YI-10A-01D-A32I-09',
                                                'TCGA-VS-A9UD-01A-11D-A42O-09,TCGA-VS-A9UD-10A-01D-A42R-09', 'TCGA-GH-A9DA-01A-21D-A37N-09,TCGA-GH-A9DA-10A-01D-A37N-09', 'TCGA-EX-A1H5-01A-31D-A13W-08,TCGA-EX-A1H5-10A-01D-A13W-08', 'TCGA-EK-A2RC-01A-11D-A18J-09,TCGA-EK-A2RC-10A-01D-A18J-09', 'TCGA-EK-A2RK-01A-11D-A18J-09,TCGA-EK-A2RK-10A-01D-A18J-09', 'TCGA-JX-A5QV-01A-22D-A28B-09,TCGA-JX-A5QV-10A-01D-A28E-09', 'TCGA-C5-A907-01A-11D-A37N-09,TCGA-C5-A907-10A-01D-A37N-09', 'TCGA-EK-A2RO-01A-11D-A18J-09,TCGA-EK-A2RO-10A-01D-A18J-09', 'TCGA-VS-A8EG-01A-11D-A36J-09,TCGA-VS-A8EG-10A-01D-A36M-09', 'TCGA-VS-A9UV-01A-11D-A42O-09,TCGA-VS-A9UV-10A-01D-A42R-09', 'TCGA-DS-A5RQ-01A-11D-A28B-09,TCGA-DS-A5RQ-10A-01D-A28E-09', 'TCGA-VS-A9UY-01A-11D-A42O-09,TCGA-VS-A9UY-10A-01D-A42R-09', 'TCGA-FU-A2QG-01A-11D-A18J-09,TCGA-FU-A2QG-10A-01D-A18J-09', 'TCGA-ZJ-A8QR-01A-11D-A37N-09,TCGA-ZJ-A8QR-10A-01D-A37N-09', 'TCGA-VS-A8EK-01A-12D-A37N-09,TCGA-VS-A8EK-10A-01D-A37N-09', 'TCGA-C5-A8YQ-01A-11D-A37N-09,TCGA-C5-A8YQ-10A-01D-A37N-09', 'TCGA-DS-A0VM-01A-11D-A10S-08,TCGA-DS-A0VM-10A-01D-A10S-08', 'TCGA-C5-A2LX-01A-11D-A18J-09,TCGA-C5-A2LX-10A-01D-A18J-09',
                                                'TCGA-EK-A3GN-01A-11D-A20U-09,TCGA-EK-A3GN-10A-01D-A20U-09', 'TCGA-JW-A5VK-01A-11D-A28B-09,TCGA-JW-A5VK-10A-01D-A28E-09', 'TCGA-MA-AA3Z-01A-11D-A387-09,TCGA-MA-AA3Z-10A-01D-A38A-09', 'TCGA-C5-A7CO-01A-11D-A351-09,TCGA-C5-A7CO-10A-01D-A351-09', 'TCGA-XS-A8TJ-01A-11D-A36J-09,TCGA-XS-A8TJ-10A-01D-A36M-09', 'TCGA-IR-A3LL-01A-11D-A20U-09,TCGA-IR-A3LL-10A-01D-A20U-09', 'TCGA-DG-A2KJ-01A-11D-A18J-09,TCGA-DG-A2KJ-10A-01D-A18J-09', 'TCGA-JW-A852-01A-11D-A351-09,TCGA-JW-A852-10A-01D-A351-09', 'TCGA-C5-A1BK-01B-11D-A13W-08,TCGA-C5-A1BK-10A-01D-A13W-08', 'TCGA-C5-A1BM-01A-11D-A13W-08,TCGA-C5-A1BM-10A-01D-A13W-08', 'TCGA-MA-AA41-01A-11D-A387-09,TCGA-MA-AA41-10A-01D-A38A-09', 'TCGA-C5-A1M6-01A-11D-A13W-08,TCGA-C5-A1M6-10A-01D-A13W-08', 'TCGA-C5-A7UH-01A-11D-A351-09,TCGA-C5-A7UH-10A-01D-A351-09', 'TCGA-EA-A411-01A-11D-A243-09,TCGA-EA-A411-10A-01D-A243-09', 'TCGA-C5-A1BL-01A-11D-A13W-08,TCGA-C5-A1BL-10A-01D-A13W-08', 'TCGA-EK-A2PK-01A-11D-A18J-09,TCGA-EK-A2PK-10A-01D-A18J-09',
                                                'TCGA-VS-A9U5-01A-11D-A42O-09,TCGA-VS-A9U5-10A-01D-A42R-09', 'TCGA-VS-A954-01A-11D-A387-09,TCGA-VS-A954-10A-01D-A38A-09', 'TCGA-C5-A8XI-01A-11D-A37N-09,TCGA-C5-A8XI-10A-01D-A37N-09', 'TCGA-FU-A5XV-01A-11D-A28B-09,TCGA-FU-A5XV-10A-01D-A28E-09', 'TCGA-C5-A7X5-01A-11D-A36J-09,TCGA-C5-A7X5-10A-01D-A36M-09', 'TCGA-EA-A1QT-01A-11D-A14W-08,TCGA-EA-A1QT-10A-01D-A14W-08', 'TCGA-VS-A958-01A-11D-A42O-09,TCGA-VS-A958-10A-01D-A42R-09', 'TCGA-C5-A905-01A-11D-A37N-09,TCGA-C5-A905-10A-01D-A37N-09', 'TCGA-ZJ-A8QQ-01A-11D-A37N-09,TCGA-ZJ-A8QQ-10A-01D-A37N-09', 'TCGA-VS-A9UB-01A-22D-A42O-09,TCGA-VS-A9UB-10A-01D-A42R-09', 'TCGA-VS-A953-01A-11D-A387-09,TCGA-VS-A953-10A-01D-A38A-09', 'TCGA-EA-A5FO-01A-21D-A28B-09,TCGA-EA-A5FO-10A-01D-A28E-09', 'TCGA-ZJ-AAXA-01A-11D-A42O-09,TCGA-ZJ-AAXA-10A-01D-A42R-09', 'TCGA-C5-A7XC-01A-11D-A387-09,TCGA-C5-A7XC-10A-01D-A38A-09', 'TCGA-DS-A0VN-01A-21D-A10S-08,TCGA-DS-A0VN-10A-01D-A10S-08', 'TCGA-C5-A901-01A-11D-A37N-09,TCGA-C5-A901-10A-01D-A37N-09', 'TCGA-VS-A8EL-01A-11D-A37N-09,TCGA-VS-A8EL-10A-01D-A37N-09', 'TCGA-LP-A4AU-01A-32D-A243-09,TCGA-LP-A4AU-10A-01D-A243-09', 'TCGA-EK-A2RB-01A-11D-A18J-09,TCGA-EK-A2RB-10A-01D-A18J-09', 'TCGA-MA-AA3W-01A-11D-A387-09,TCGA-MA-AA3W-10A-01D-A38A-09',
                                                'TCGA-VS-A8QC-01A-11D-A37N-09,TCGA-VS-A8QC-10A-01D-A37N-09', 'TCGA-VS-A9UC-01A-11D-A42O-09,TCGA-VS-A9UC-10A-01D-A42R-09', 'TCGA-C5-A7CH-01A-11D-A33O-09,TCGA-C5-A7CH-10A-01D-A33O-09', 'TCGA-C5-A8XJ-01A-11D-A37N-09,TCGA-C5-A8XJ-10A-01D-A37N-09', 'TCGA-MA-AA43-01A-11D-A42O-09,TCGA-MA-AA43-10A-01D-A45U-09', 'TCGA-JW-A5VJ-01A-11D-A28B-09,TCGA-JW-A5VJ-10A-01D-A28E-09', 'TCGA-VS-A959-01A-11D-A42O-09,TCGA-VS-A959-10A-01D-A42R-09', 'TCGA-Q1-A5R2-01A-11D-A28B-09,TCGA-Q1-A5R2-10A-01D-A28E-09', 'TCGA-BI-A20A-01A-11D-A14W-08,TCGA-BI-A20A-10A-01D-A14W-08', 'TCGA-C5-A7CL-01A-11D-A32I-09,TCGA-C5-A7CL-10A-01D-A32I-09', 'TCGA-MY-A913-01A-11D-A37N-09,TCGA-MY-A913-10A-01D-A37N-09', 'TCGA-EK-A2RM-01A-21D-A18J-09,TCGA-EK-A2RM-10A-01D-A18J-09', 'TCGA-VS-A957-01A-11D-A42O-09,TCGA-VS-A957-10A-01D-A42R-09', 'TCGA-C5-A7UI-01A-11D-A36J-09,TCGA-C5-A7UI-10B-01D-A36M-09', 'TCGA-C5-A8YT-01A-11D-A37N-09,TCGA-C5-A8YT-10A-01D-A37N-09',
                                                'TCGA-FU-A3NI-01A-11D-A21Q-09,TCGA-FU-A3NI-10A-01D-A21Q-09', 'TCGA-EK-A2R7-01A-11D-A18J-09,TCGA-EK-A2R7-10A-01D-A18J-09', 'TCGA-Q1-A73O-01A-11D-A32I-09,TCGA-Q1-A73O-10B-01D-A32I-09', 'TCGA-JW-A5VL-01A-11D-A28B-09,TCGA-JW-A5VL-10A-01D-A28E-09', 'TCGA-DS-A1OD-01A-11D-A14W-08,TCGA-DS-A1OD-10A-01D-A14W-08', 'TCGA-ZJ-AAXU-01A-11D-A42O-09,TCGA-ZJ-AAXU-10A-01D-A42R-09', 'TCGA-C5-A7UE-01A-11D-A33O-09,TCGA-C5-A7UE-10A-01D-A33O-09', 'TCGA-DG-A2KM-01A-11D-A17W-09,TCGA-DG-A2KM-10A-01D-A17W-09', 'TCGA-VS-A9UM-01A-11D-A42O-09,TCGA-VS-A9UM-10A-01D-A42R-09', 'TCGA-C5-A1M7-01A-11D-A13W-08,TCGA-C5-A1M7-10A-01D-A13W-08', 'TCGA-C5-A7CJ-01A-11D-A32I-09,TCGA-C5-A7CJ-10A-01D-A32I-09', 'TCGA-UC-A7PG-01A-11D-A42O-09,TCGA-UC-A7PG-11C-11D-A42R-09', 'TCGA-EA-A3HT-01A-61D-A21Q-09,TCGA-EA-A3HT-10A-01D-A21Q-09', 'TCGA-C5-A7UC-01A-11D-A351-09,TCGA-C5-A7UC-10A-01D-A351-09', 'TCGA-VS-A8Q9-01A-12D-A37N-09,TCGA-VS-A8Q9-10A-01D-A37N-09', 'TCGA-DS-A0VK-01A-21D-A10S-08,TCGA-DS-A0VK-10A-01D-A10S-08', 'TCGA-VS-A94Z-01A-11D-A387-09,TCGA-VS-A94Z-10A-01D-A38A-09',
                                                'TCGA-FU-A3HZ-01A-11D-A20U-09,TCGA-FU-A3HZ-10A-01D-A20U-09', 'TCGA-C5-A902-01A-11D-A37N-09,TCGA-C5-A902-10A-01D-A37N-09', 'TCGA-EK-A2PI-01A-11D-A18J-09,TCGA-EK-A2PI-10A-01D-A18J-09', 'TCGA-HG-A2PA-01A-11D-A20U-09,TCGA-HG-A2PA-10B-01D-A20U-09', 'TCGA-C5-A1MQ-01A-11D-A14W-08,TCGA-C5-A1MQ-10A-01D-A14W-08', 'TCGA-C5-A2LY-01A-31D-A18J-09,TCGA-C5-A2LY-10A-01D-A18J-09', 'TCGA-HM-A3JJ-01A-21D-A21Q-09,TCGA-HM-A3JJ-10A-01D-A21Q-09', 'TCGA-MA-AA42-01A-12D-A387-09,TCGA-MA-AA42-10A-01D-A38A-09', 'TCGA-DS-A0VL-01A-21D-A10S-08,TCGA-DS-A0VL-10A-01D-A10S-08', 'TCGA-VS-A9UU-01A-11D-A42O-09,TCGA-VS-A9UU-10B-01D-A42R-09', 'TCGA-DR-A0ZL-01A-11D-A10S-08,TCGA-DR-A0ZL-10A-01D-A10S-08', 'TCGA-C5-A2LV-01A-11D-A18J-09,TCGA-C5-A2LV-10A-01D-A18J-09', 'TCGA-VS-A9V2-01A-11D-A42O-09,TCGA-VS-A9V2-10A-01D-A42R-09', 'TCGA-BI-A0VS-01A-11D-A10S-08,TCGA-BI-A0VS-10A-01D-A10S-08',
                                                'TCGA-Q1-A6DT-01A-11D-A32I-09,TCGA-Q1-A6DT-10A-01D-A32I-09', 'TCGA-DS-A1OA-01A-11D-A16Y-08,TCGA-DS-A1OA-10A-01D-A16Y-08', 'TCGA-LP-A4AX-01A-12D-A243-09,TCGA-LP-A4AX-10A-01D-A243-09', 'TCGA-EA-A43B-01A-81D-A243-09,TCGA-EA-A43B-10A-01D-A243-09', 'TCGA-EK-A2IP-01A-11D-A17W-09,TCGA-EK-A2IP-10A-01D-A17W-09', 'TCGA-IR-A3L7-01A-21D-A20U-09,TCGA-IR-A3L7-10A-01D-A20U-09', 'TCGA-FU-A3HY-01A-11D-A21Q-09,TCGA-FU-A3HY-10A-01D-A21Q-09', 'TCGA-BI-A0VR-01A-11D-A10S-08,TCGA-BI-A0VR-10A-01D-A10S-08', 'TCGA-ZX-AA5X-01A-11D-A42O-09,TCGA-ZX-AA5X-10A-01D-A42R-09',
                                                'TCGA-EK-A2H1-01A-11D-A17W-09,TCGA-EK-A2H1-10A-01D-A17W-09', 'TCGA-FU-A23K-01A-11D-A16O-08,TCGA-FU-A23K-10A-01D-A16O-08', 'TCGA-VS-A9UJ-01A-11D-A42O-09,TCGA-VS-A9UJ-10A-01D-A42R-09', 'TCGA-C5-A1BI-01B-11D-A13W-08,TCGA-C5-A1BI-10A-01D-A13W-08', 'TCGA-VS-A8EC-01A-11D-A36J-09,TCGA-VS-A8EC-10A-01D-A36M-09', 'TCGA-VS-A8EI-01A-11D-A37N-09,TCGA-VS-A8EI-10A-01D-A37N-09', 'TCGA-EA-A3HU-01A-11D-A20U-09,TCGA-EA-A3HU-10B-01D-A20U-09', 'TCGA-C5-A1MI-01A-11D-A14W-08,TCGA-C5-A1MI-10A-01D-A14W-08', 'TCGA-EK-A2RN-01A-12D-A20U-09,TCGA-EK-A2RN-10A-01D-A20U-09', 'TCGA-RA-A741-01A-11D-A33O-09,TCGA-RA-A741-10B-01D-A33O-09', 'TCGA-LP-A5U3-01A-11D-A28B-09,TCGA-LP-A5U3-10A-01D-A28E-09', 'TCGA-JW-AAVH-01A-11D-A387-09,TCGA-JW-AAVH-10A-01D-A38A-09', 'TCGA-VS-A94W-01A-12D-A37N-09,TCGA-VS-A94W-10A-01D-A37N-09',
                                                'TCGA-C5-A1MH-01A-11D-A14W-08,TCGA-C5-A1MH-10A-01D-A14W-08', 'TCGA-VS-A8EB-01A-11D-A36J-09,TCGA-VS-A8EB-10A-01D-A36M-09', 'TCGA-C5-A8ZZ-01A-11D-A37N-09,TCGA-C5-A8ZZ-10A-01D-A37N-09', 'TCGA-WL-A834-01A-11D-A351-09,TCGA-WL-A834-10A-01D-A351-09', 'TCGA-PN-A8MA-01A-11D-A36J-09,TCGA-PN-A8MA-10A-01D-A36M-09', 'TCGA-VS-A950-01A-11D-A42O-09,TCGA-VS-A950-10B-01D-A42R-09', 'TCGA-Q1-A6DW-01A-11D-A32I-09,TCGA-Q1-A6DW-10B-01D-A32I-09', 'TCGA-VS-A9UH-01A-11D-A42O-09,TCGA-VS-A9UH-10A-01D-A42R-09', 'TCGA-FU-A3YQ-01A-11D-A22X-09,TCGA-FU-A3YQ-10A-01D-A22X-09', 'TCGA-VS-A8EH-01A-11D-A36J-09,TCGA-VS-A8EH-10A-01D-A36M-09', 'TCGA-FU-A3WB-01A-11D-A22X-09,TCGA-FU-A3WB-10A-01D-A22X-09', 'TCGA-VS-A9UL-01A-11D-A42O-09,TCGA-VS-A9UL-10A-01D-A42R-09', 'TCGA-EA-A3HQ-01A-11D-A20U-09,TCGA-EA-A3HQ-10A-01D-A20U-09', 'TCGA-VS-A8Q8-01A-11D-A37N-09,TCGA-VS-A8Q8-10A-01D-A37N-09', 'TCGA-EA-A3QE-01A-21D-A21Q-09,TCGA-EA-A3QE-10A-01D-A21Q-09',
                                                'TCGA-EA-A3QD-01A-32D-A22X-09,TCGA-EA-A3QD-10A-01D-A22X-09', 'TCGA-C5-A8YR-01A-12D-A37N-09,TCGA-C5-A8YR-10A-01D-A37N-09', 'TCGA-FU-A3TQ-01A-11D-A22X-09,TCGA-FU-A3TQ-10A-01D-A22X-09', 'TCGA-Q1-A5R3-01A-11D-A28B-09,TCGA-Q1-A5R3-10B-01D-A28E-09', 'TCGA-C5-A1BJ-01A-11D-A13W-08,TCGA-C5-A1BJ-10A-01D-A13W-08', 'TCGA-LP-A4AW-01A-11D-A243-09,TCGA-LP-A4AW-10A-01D-A243-09', 'TCGA-C5-A8XK-01A-11D-A37N-09,TCGA-C5-A8XK-10A-01D-A37N-09', 'TCGA-DG-A2KK-01A-11D-A17W-09,TCGA-DG-A2KK-10A-01D-A17W-09', 'TCGA-VS-A9UI-01A-11D-A42O-09,TCGA-VS-A9UI-10A-01D-A42R-09', 'TCGA-EX-A69L-01A-11D-A32I-09,TCGA-EX-A69L-10A-01D-A32I-09', 'TCGA-VS-A9U6-01A-11D-A42O-09,TCGA-VS-A9U6-10A-01D-A42R-09', 'TCGA-EK-A2RE-01A-11D-A18J-09,TCGA-EK-A2RE-10A-01D-A18J-09', 'TCGA-MA-AA3X-01A-22D-A42O-09,TCGA-MA-AA3X-10A-01D-A42R-09', 'TCGA-VS-A9U7-01A-11D-A42O-09,TCGA-VS-A9U7-10A-01D-A42R-09', 'TCGA-FU-A3TX-01A-11D-A22X-09,TCGA-FU-A3TX-10A-01D-A22X-09')) 

output_query_mutation_squamous <- getResults(query_mutation_squamous) ## 182 files and 182 cases
write_xlsx(output_query_mutation_squamous, 'TCGA_mutation_squamous_List.xlsx')

setwd("~/Escritorio/R/TCGA_2")
maf_squamous <- GDCprepare(query_mutation_squamous, summarizedExperiment = TRUE)
setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms")
options(max.print=999999)
capture.output(maf_squamous, file = "MAF_List_Squamous.txt")

maftools.input.squamous <- read.maf(maf_squamous)
capture.output(maftools.input.squamous, file = "MAF_List_Squamous_Summary.txt")

maftools.squamous.sample.summary <- getSampleSummary(maftools.input.squamous) ## Sample summary
write_xlsx(maftools.squamous.sample.summary, 'MAF_Squamous_sample_summary.xlsx')

maftools.squamous.gene.summary <- getGeneSummary(maftools.input.squamous)
write_xlsx(maftools.squamous.gene.summary, 'MAF_Squamous_gene_summary.xlsx')

getClinicalData(maftools.input.squamous)

getFields(maftools.input.squamous)

#### Plotting MAF Summary
plotmafSummary(maf = maftools.input.squamous,
               addStat = "median",
               rmOutlier = TRUE,
               dashboard = TRUE,
               titvRaw = FALSE)
##### Save as "Summary_MAF_Squamous.pdf" (Landscape - 7 x 5.50)

#### Oncoplots
col = RColorBrewer::brewer.pal(n = 8, name = 'Set1') ## https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
names(col) = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', ' Multi_Hit', 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
oncoplot(maf = maftools.input.squamous, top = 15, removeNonMutated = TRUE, sortByMutation = TRUE, colors = col, draw_titv = FALSE)
##### Save as "Oncoplot_Squamous.pdf" (Landscape - 7.20 x 4.11)

#### Transition and Transversions
squamous_titv = titv(maf = maftools.input.squamous, plot = TRUE, useSyn = TRUE)
##### Save as "TiTv_Squamous.pdf" (Landscape)

#### Lollipop plots for amino acid changes
TTN_squamous <- lollipopPlot(
  maf = maftools.input.squamous,
  gene = 'TTN',
  AACol = 'HGVSp_Short',
  labelPos = '28251',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(TTN_squamous, file = "TTN_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "1st_TTN_mutations.pdf" (Landscape)

PIK3CA_squamous <- lollipopPlot(
  maf = maftools.input.squamous,
  gene = 'PIK3CA',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  refSeqID = 'NM_006218',
  proteinID = 'P42336',
  labelPos = c(112,542,545,726),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(PIK3CA_squamous, file = "PIK3CA_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "2nd_PIK3CA_mutations.pdf" (Landscape)

KMT2C_squamous <- lollipopPlot(
  maf = maftools.input.squamous,
  gene = 'KMT2C',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  labelPos = c(284,430,2161,2481),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(KMT2C_squamous, file = "KMT2C_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "3rd_KMT2C_mutations.pdf" (Landscape)


#### Rainfall plot (localize hypermutations - kataegis) 
rainfallPlot(maftools.input.squamous, tsb = 'TCGA-JW-A5VL-01A-11D-A28B-09', detectChangePoints = TRUE) ## 2nd - maftools.squamous.sample.summary
##### Save as "Rainfall_plot_Squamous_2nd.pdf" (Landscape)

#### Compare mutation load against TCGA cohorts
laml.mutload.squa = tcgaCompare(maf = maftools.input.squamous, cohortName = 'Squamous', logscale = TRUE, capture_size = 50)
write_xlsx(laml.mutload.squa, 'Mutation_Load_TCGA_Cohorts.xlsx')
##### Save as "Comparison_mutation_load_squamous.pdf" (Landscape)

### Plotting Variant Allele Frequencies
plotVaf(maf = maftools.input.squamous)
##### Save as "VAF_Squamous.pdf" (Landscape)

## Exclusive/co-occurance event analysis on top 20 mutated genes
somatic_interactions_squamous <- somaticInteractions(maf = maftools.input.squamous, top = 20, pvalue = c(0.05, 0.1))
write_xlsx(somatic_interactions_squamous, 'Somatic_interactions_squamous.xlsx')
##### Save as "Co_occurance_Squamous.pdf" (Landscape - 7 x 6)

## Detecting cancer driver genes based on positional clustering
squamous.sig = oncodrive(maf = maftools.input.squamous, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
write_xlsx(squamous.sig, 'Driver_genes_clustering.xlsx')
head(squamous.sig)
plotOncodrive(res = squamous.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.7)
##### Save as "OncoDrive_Squamous.pdf" (Landscape)

## Drug-Gene Interactions
dgi.squamous = drugInteractions(maf = maftools.input.squamous, fontSize = 0.75)
write_xlsx(dgi.squamous, 'Drug_gene_interactions_squamous.xlsx')
##### Save as "Drug_Squamous.pdf" (Landscape)

## Oncogenic Signaling Pathways
OncogenicPathways(maf = maftools.input.squamous)
##### Save as "OncogenicPathways_Squamous.pdf" (Landscape)
PlotOncogenicPathways(maf = maftools.input.squamous, pathways = "RTK-RAS")
##### Save as "OncogenicPathway_RTK-RAS_Squamous.pdf" (Landscape - 8 x 7)

######################################################################################################################################

## Download and visualize mutation data from TCGA-Adenomas and adenocarcinomas
setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas")
query_mutation_adenomas <- GDCquery(project = 'TCGA-CESC',
                                    data.category = 'Simple Nucleotide Variation', 
                                    access = 'open',
                                    barcode = c('TCGA-C5-A1MJ-01A-11D-A14W-08,TCGA-C5-A1MJ-10A-01D-A14W-08', 'TCGA-LP-A5U2-01A-11D-A28B-09,TCGA-LP-A5U2-10A-01D-A28E-09', 
                                                'TCGA-EX-A1H6-01B-11D-A22X-09,TCGA-EX-A1H6-10A-01D-A22X-09', 'TCGA-C5-A2LS-01A-22D-A22X-09,TCGA-C5-A2LS-10A-01D-A22X-09',
                                                'TCGA-FU-A770-01A-11D-A33O-09,TCGA-FU-A770-10A-01D-A33O-09', 'TCGA-FU-A40J-01A-11D-A243-09,TCGA-FU-A40J-10A-01D-A243-09',
                                                'TCGA-C5-A1ME-01A-11D-A13W-08,TCGA-C5-A1ME-10A-01D-A13W-08',
                                                'TCGA-C5-A1M9-01A-11D-A13W-08,TCGA-C5-A1M9-10A-01D-A13W-08', 'TCGA-FU-A3EO-01A-11D-A20U-09,TCGA-FU-A3EO-11A-13D-A20U-09',
                                                'TCGA-C5-A7CM-01A-11D-A33O-09,TCGA-C5-A7CM-10A-01D-A33O-09', 'TCGA-JW-A69B-01A-11D-A32I-09,TCGA-JW-A69B-10A-01D-A32I-09',
                                                'TCGA-EK-A3GK-01A-11D-A20U-09,TCGA-EK-A3GK-10A-01D-A20U-09',
                                                'TCGA-Q1-A73R-01A-11D-A33O-09,TCGA-Q1-A73R-10B-01D-A33O-09', 'TCGA-EK-A2RL-01A-11D-A18J-09,TCGA-EK-A2RL-10A-01D-A18J-09',
                                                'TCGA-UC-A7PI-01A-11D-A42O-09,TCGA-UC-A7PI-10A-01D-A42R-09', 'TCGA-DG-A2KH-01A-21D-A22X-09,TCGA-DG-A2KH-10A-01D-A22X-09',
                                                'TCGA-C5-A2M2-01A-21D-A18J-09,TCGA-C5-A2M2-10A-01D-A18J-09', 'TCGA-Q1-A73P-01A-11D-A32I-09,TCGA-Q1-A73P-10B-01D-A32I-09',
                                                'TCGA-JX-A3Q8-01A-11D-A21Q-09,TCGA-JX-A3Q8-10A-01D-A21Q-09', 'TCGA-IR-A3LF-01A-21D-A22X-09,TCGA-IR-A3LF-10A-01D-A22X-09',
                                                'TCGA-IR-A3LI-01A-11D-A20U-09,TCGA-IR-A3LI-10A-01D-A20U-09',
                                                'TCGA-LP-A7HU-01A-11D-A33O-09,TCGA-LP-A7HU-10A-01D-A33O-09', 'TCGA-C5-A2M1-01A-11D-A18J-09,TCGA-C5-A2M1-10A-01D-A18J-09')) ## From cases column in R

output_query_mutation_adenomas <- getResults(query_mutation_adenomas) ## 23 files and 23 cases
write_xlsx(output_query_mutation_adenomas, 'TCGA_mutation_adenomas_List.xlsx')

setwd("~/Escritorio/R/TCGA_2")
maf_adenomas <- GDCprepare(query_mutation_adenomas, summarizedExperiment = TRUE)
setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas")
options(max.print=999999)
capture.output(maf_adenomas, file = "MAF_List_Adenomas.txt")

maftools.input.adenomas <- read.maf(maf_adenomas)
capture.output(maftools.input.adenomas, file = "MAF_List_Adenomas_Summary.txt")

maftools.adenomas.sample.summary <- getSampleSummary(maftools.input.adenomas) ## Sample summary
write_xlsx(maftools.adenomas.sample.summary, 'MAF_Adenomas_sample_summary.xlsx')

maftools.adenomas.gene.summary <- getGeneSummary(maftools.input.adenomas)
write_xlsx(maftools.adenomas.gene.summary, 'MAF_Adenomas_gene_summary.xlsx')

getClinicalData(maftools.input.adenomas)

getFields(maftools.input.adenomas)

#### Plotting MAF Summary
plotmafSummary(maf = maftools.input.adenomas,
               addStat = "median",
               rmOutlier = TRUE,
               dashboard = TRUE,
               titvRaw = FALSE)
##### Save as "Summary_MAF_Adenomas.pdf" (Landscape - 7 x 5.50)

#### Oncoplots
col = RColorBrewer::brewer.pal(n = 8, name = 'Set1') ## https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
names(col) = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', ' Multi_Hit', 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
oncoplot(maf = maftools.input.adenomas, top = 15, removeNonMutated = TRUE, sortByMutation = TRUE, colors = col, draw_titv = FALSE)
##### Save as "Oncoplot_Adenomas.pdf" (Landscape - 7.20 x 4.11)

#### Transition and Transversions
adenomas_titv = titv(maf = maftools.input.adenomas, plot = TRUE, useSyn = TRUE)
##### Save as "TiTv_Adenomas.pdf" (Landscape)

#### Lollipop plots for amino acid changes
PIK3CA_adenomas <- lollipopPlot(
  maf = maftools.input.adenomas,
  gene = 'PIK3CA',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  refSeqID = 'NM_006218',
  proteinID = 'P42336',
  labelPos = c(542,545),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(PIK3CA_adenomas, file = "PIK3CA_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "1st_PIK3CA_mutations.pdf" (Landscape)

TTN_adenomas <- lollipopPlot(
  maf = maftools.input.adenomas,
  gene = 'TTN',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(TTN_adenomas, file = "TTN_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "2nd_TTN_mutations.pdf" (Landscape)

ERBB2_adenomas <- lollipopPlot(
  maf = maftools.input.adenomas,
  gene = 'ERBB2',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  labelPos = c(310),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(ERBB2_adenomas, file = "ERBB2_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "3rd_ERBB2_mutations.pdf" (Landscape)


#### Rainfall plot (localize hypermutations - kataegis) 
rainfallPlot(maftools.input.adenomas, tsb = 'TCGA-EK-A3GK-01A-11D-A20U-09', detectChangePoints = TRUE) ## 1st - maftools.adenomas.sample.summary
##### Save as "Rainfall_plot_Adenomas_1st.pdf" (Landscape)

#### Compare mutation load against TCGA cohorts
laml.mutload.adenomas = tcgaCompare(maf = maftools.input.adenomas, cohortName = 'Adenomas', logscale = TRUE, capture_size = 50)
write_xlsx(laml.mutload.adenomas, 'Mutation_Load_TCGA_Cohorts.xlsx')
##### Save as "Comparison_mutation_load_adenomas.pdf" (Landscape)

### Plotting Variant Allele Frequencies
plotVaf(maf = maftools.input.adenomas)
##### Save as "VAF_Adenomas.pdf" (Landscape)

# Exclusive/co-occurance event analysis on top 20 mutated genes
somatic_interactions_adenomas <- somaticInteractions(maf = maftools.input.adenomas, top = 20, pvalue = c(0.05, 0.1))
write_xlsx(somatic_interactions_adenomas, 'Somatic_interactions_adenomas.xlsx')
##### Save as "Co_occurance_Adenomas.pdf" (Landscape 7 x 6)

## Detecting cancer driver genes based on positional clustering
adenomas.sig = oncodrive(maf = maftools.input.adenomas, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
write_xlsx(adenomas.sig, 'Driver_genes_clustering.xlsx')
head(adenomas.sig)
plotOncodrive(res = adenomas.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.7)
##### Save as "OncoDrive_Adenomas.pdf" (Landscape)

## Drug-Gene Interactions
dgi.adenomas = drugInteractions(maf = maftools.input.adenomas, fontSize = 0.75)
write_xlsx(dgi.adenomas, 'Drug_gene_interactions_adenomas.xlsx')
##### Save as "Drug_Adenomas.pdf" (Landscape)

## Oncogenic Signaling Pathways
OncogenicPathways(maf = maftools.input.adenomas)
##### Save as "OncogenicPathways_Adenomas.pdf" (Landscape)
PlotOncogenicPathways(maf = maftools.input.adenomas, pathways = "RTK-RAS")
##### Save as "OncogenicPathway_RTK-RAS_Adenomas.pdf" (Landscape)

######################################################################################################################################

## Download and visualize mutation data from TCGA-Cystic, mucinous and serous neoplasms
setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms")
query_mutation_cystic <- GDCquery(project = 'TCGA-CESC',
                                    data.category = 'Simple Nucleotide Variation',
                                    access = 'open',
                                    barcode = c('TCGA-VS-A9UQ-01A-21D-A42O-09,TCGA-VS-A9UQ-10A-01D-A42R-09', 'TCGA-VS-A9V4-01A-12D-A42O-09,TCGA-VS-A9V4-10A-01D-A42R-09',
                                                'TCGA-VS-A9V5-01A-11D-A42O-09,TCGA-VS-A9V5-10A-01D-A42R-09', 'TCGA-Q1-A5R1-01A-11D-A28B-09,TCGA-Q1-A5R1-10A-01D-A28E-09',
                                                'TCGA-VS-A9V1-01A-11D-A42O-09,TCGA-VS-A9V1-10A-01D-A42R-09', 'TCGA-VS-A8QH-01A-11D-A37N-09,TCGA-VS-A8QH-10A-01D-A37N-09',
                                                'TCGA-Q1-A6DV-01A-11D-A32I-09,TCGA-Q1-A6DV-10A-01D-A32I-09', 'TCGA-DS-A7WH-01A-22D-A351-09,TCGA-DS-A7WH-10A-01D-A351-09',
                                                'TCGA-VS-A9UP-01A-11D-A42O-09,TCGA-VS-A9UP-10A-01D-A42R-09', 'TCGA-VS-A9UO-01A-11D-A42O-09,TCGA-VS-A9UO-10A-01D-A42R-09',
                                                'TCGA-C5-A7X8-01A-11D-A36J-09,TCGA-C5-A7X8-10A-01D-A36M-09', 'TCGA-VS-A952-01A-11D-A387-09,TCGA-VS-A952-10A-01D-A38A-09',
                                                'TCGA-VS-A9UT-01A-11D-A42O-09,TCGA-VS-A9UT-10A-01D-A42R-09', 'TCGA-VS-A9UZ-01A-11D-A42O-09,TCGA-VS-A9UZ-10A-01D-A42R-09',
                                                'TCGA-VS-A9UR-01A-11D-A42O-09,TCGA-VS-A9UR-10A-01D-A42R-09', 'TCGA-VS-A9V0-01A-11D-A42O-09,TCGA-VS-A9V0-10A-01D-A42R-09',
                                                'TCGA-C5-A3HF-01A-11D-A20U-09,TCGA-C5-A3HF-10A-01D-A20U-09')) ## From cases column in R

output_query_mutation_cystic <- getResults(query_mutation_cystic) ## 17 files and 17 cases
write_xlsx(output_query_mutation_cystic, 'TCGA_mutation_cystic_List.xlsx')

setwd("~/Escritorio/R/TCGA_2")
maf_cystic <- GDCprepare(query_mutation_cystic, summarizedExperiment = TRUE)
setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms")
options(max.print=999999)
capture.output(maf_cystic, file = "MAF_List_Cystic.txt")

maftools.input.cystic <- read.maf(maf_cystic)
capture.output(maftools.input.cystic, file = "MAF_List_Cystic_Summary.txt")

maftools.cystic.sample.summary <- getSampleSummary(maftools.input.cystic) ## Sample summary
write_xlsx(maftools.cystic.sample.summary, 'MAF_Cystic_sample_summary.xlsx')

maftools.cystic.gene.summary <- getGeneSummary(maftools.input.cystic)
write_xlsx(maftools.cystic.gene.summary, 'MAF_Cystic_gene_summary.xlsx')

getClinicalData(maftools.input.cystic)

getFields(maftools.input.cystic)

#### Plotting MAF Summary
plotmafSummary(maf = maftools.input.cystic,
               addStat = "median",
               rmOutlier = TRUE,
               dashboard = TRUE,
               titvRaw = FALSE)
##### Save as "Summary_MAF_Cystic.pdf" (Landscape - 7 x 5.50)

#### Oncoplots
col = RColorBrewer::brewer.pal(n = 8, name = 'Set1') ## https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
names(col) = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', ' Multi_Hit', 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
oncoplot(maf = maftools.input.cystic, top = 15, removeNonMutated = TRUE, sortByMutation = TRUE, colors = col, draw_titv = FALSE)
##### Save as "Oncoplot_Cystic.pdf" (Landscape - 7.20 x 4.11)

#### Transition and Transversions
cystic_titv = titv(maf = maftools.input.cystic, plot = TRUE, useSyn = TRUE)
##### Save as "TiTv_Cystic.pdf" (Landscape)

#### Lollipop plots for amino acid changes
PIK3CA_cystic <- lollipopPlot(
  maf = maftools.input.cystic,
  gene = 'PIK3CA',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  labelPos = c(542,545),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(PIK3CA_cystic, file = "PIK3CA_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "1st_PIK3CA_mutations.pdf" (Landscape)

TP53_cystic <- lollipopPlot(
  maf = maftools.input.cystic,
  gene = 'TP53',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  labelPos = c(27,331),
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(TP53_cystic, file = "TP53_cystic.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "2nd_TP53_mutations.pdf" (Landscape)

LRP1B_cystic <- lollipopPlot(
  maf = maftools.input.cystic,
  gene = 'LRP1B',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  showDomainLabel = FALSE,
  printCount = TRUE,
  labelOnlyUniqueDoamins = FALSE)
write.table(LRP1B_cystic, file = "LRP1B_mutations.txt", sep = ",", quote = FALSE, row.names = FALSE)
##### Save as "3rd_LRP1B_mutations.pdf" (Landscape)

#### Rainfall plot (localize hypermutations - kataegis) 
## Any of the samples have hypermutations

#### Compare mutation load against TCGA cohorts
laml.mutload.cystic = tcgaCompare(maf = maftools.input.cystic, cohortName = 'Cystic', logscale = TRUE, capture_size = 50)
write_xlsx(laml.mutload.cystic, 'Mutation_Load_TCGA_Cohorts.xlsx')
##### Save as "Comparison_mutation_load_cystic.pdf" (Landscape)

### Plotting Variant Allele Frequencies
plotVaf(maf = maftools.input.cystic)
##### Save as "VAF_Cystic.pdf" (Landscape)

# Exclusive/co-occurance event analysis on top 20 mutated genes
somatic_interactions_cystic <- somaticInteractions(maf = maftools.input.cystic, top = 20, pvalue = c(0.05, 0.1))
write_xlsx(somatic_interactions_cystic, 'Somatic_interactions_cystic.xlsx')
##### Save as "Co_occurance_Cystic.pdf" (Landscape)

## Detecting cancer driver genes based on positional clustering
cystic.sig = oncodrive(maf = maftools.input.cystic, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
write_xlsx(cystic.sig, 'Driver_genes_clustering.xlsx')
head(cystic.sig)
plotOncodrive(res = cystic.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.7)
##### Save as "OncoDrive_Cystic.pdf" (Landscape)

## Drug-Gene Interactions
dgi.cystic = drugInteractions(maf = maftools.input.cystic, fontSize = 0.75)
write_xlsx(dgi.cystic, 'Drug_gene_interactions_cystic.xlsx')
##### Save as "Drug_Cystic.pdf" (Landscape)

## Oncogenic Signaling Pathways
OncogenicPathways(maf = maftools.input.cystic)
##### Save as "OncogenicPathways_Cystic.pdf" (Landscape)
PlotOncogenicPathways(maf = maftools.input.cystic, pathways = "RTK-RAS")
##### Save as "OncogenicPathway_RTK-RAS_Cystic.pdf" (Landscape)


#######################################################################################################
#######################################################################################################

# GISTIC Object

setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms/GISTIC2")
### https://www.youtube.com/watch?v=Ssw7Ryao1x4
tcga.cesc.gistic <- readGistic(gisticAllLesionsFile = "all_lesions.conf_90.txt", 
                               gisticAmpGenesFile = "amp_genes.conf_90.txt", 
                               gisticDelGenesFile = "del_genes.conf_90.txt",
                               gisticScoresFile = "scores.gistic", 
                               isTCGA = TRUE)
tcga.cesc.gistic

data_squamous <- tcga.cesc.gistic@data
write_xlsx(data_squamous, '2_Data_GISTIC_Squamous.xlsx')

cnv_summary_squamous <- tcga.cesc.gistic@cnv.summary
write_xlsx(cnv_summary_squamous, '3_CNV_Summary_Squamous.xlsx')

cytoband_summary_squamous <- tcga.cesc.gistic@cytoband.summary
write_xlsx(cytoband_summary_squamous, '4_Cytoband_Summary_Squamous.xlsx')

gisc_scores_squamous <- tcga.cesc.gistic@gis.scores
write_xlsx(gisc_scores_squamous, '5_GISC_Scores_Squamous.xlsx')

## Genome Plot
gisticChromPlot(gistic = tcga.cesc.gistic, fdrCutOff = 0.05, markBands = "all")
##### Save as "Genome_Plot_Squamous.pdf" (Landscape)

## Bubble Plot
gisticBubblePlot(gistic = tcga.cesc.gistic)
##### Save as "Bubble_Plot_Squamous.pdf" (Landscape)


setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas/GISTIC2")
## GISTIC Object
### https://www.youtube.com/watch?v=Ssw7Ryao1x4
tcga.cesc.gistic <- readGistic(gisticAllLesionsFile = "all_lesions.conf_90.txt", 
                               gisticAmpGenesFile = "amp_genes.conf_90.txt", 
                               gisticDelGenesFile = "del_genes.conf_90.txt",
                               gisticScoresFile = "scores.gistic", 
                               isTCGA = TRUE)
tcga.cesc.gistic

data_adenomas <- tcga.cesc.gistic@data
write_xlsx(data_adenomas, '2_Data_GISTIC_Adenomas.xlsx')

cnv_summary_adenomas <- tcga.cesc.gistic@cnv.summary
write_xlsx(cnv_summary_adenomas, '3_CNV_Summary_Adenomas.xlsx')

cytoband_summary_adenomas <- tcga.cesc.gistic@cytoband.summary
write_xlsx(cytoband_summary_adenomas, '4_Cytoband_Summary_Adenomas.xlsx')

gisc_scores_adenomas <- tcga.cesc.gistic@gis.scores
write_xlsx(gisc_scores_adenomas, '5_GISC_Scores_Adenomas.xlsx')

## Genome Plot
gisticChromPlot(gistic = tcga.cesc.gistic, fdrCutOff = 0.05, markBands = "all")
##### Save as "Genome_Plot_Adenomas.pdf" (Landscape)
gisticBubblePlot(gistic = tcga.cesc.gistic)
##### Save as "Bubble_Plot_Adenomas.pdf" (Landscape)


setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms/GISTIC2")
## GISTIC Object
### https://www.youtube.com/watch?v=Ssw7Ryao1x4
tcga.cesc.gistic <- readGistic(gisticAllLesionsFile = "all_lesions.conf_90.txt", 
                               gisticAmpGenesFile = "amp_genes.conf_90.txt", 
                               gisticDelGenesFile = "del_genes.conf_90.txt",
                               gisticScoresFile = "scores.gistic", 
                               isTCGA = TRUE)
tcga.cesc.gistic

data_cystic <- tcga.cesc.gistic@data
write_xlsx(data_cystic, '2_Data_GISTIC_Cystic.xlsx')

cnv_summary_cystic <- tcga.cesc.gistic@cnv.summary
write_xlsx(cnv_summary_cystic, '3_CNV_Summary_Cystic.xlsx')

cytoband_summary_cystic <- tcga.cesc.gistic@cytoband.summary
write_xlsx(cytoband_summary_cystic, '4_Cytoband_Summary_Cystic.xlsx')

gisc_scores_cystic <- tcga.cesc.gistic@gis.scores
write_xlsx(gisc_scores_cystic, '5_GISC_Scores_Cystic.xlsx')

## Genome Plot
gisticChromPlot(gistic = tcga.cesc.gistic, fdrCutOff = 0.05, markBands = "all")
##### Save as "Genome_Plot_Cystic.pdf" (Landscape)
## Bubble Plot
gisticBubblePlot(gistic = tcga.cesc.gistic)
##### Save as "Bubble_Plot_Cystic.pdf" (Landscape)


#######################################################################################################
#######################################################################################################

# Gene enrichment analysis
## https://rpubs.com/Alexis22/ClusterProfiler
library('org.Hs.eg.db')
library('clusterProfiler')

setwd('~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms/')
library('readxl')
Genes <- read_excel('MAF_Squamous_gene_summary.xlsx') 
genes_total <- Genes[,c(1,11)]
genes_total <- genes_total[order(genes_total$total),]
genes_total$order <- 1:10991

genes_total <- genes_total[genes_total$total>0,]
quantile(genes_total$total, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
# selecting the 10th percentile as low expression genes
genes_top <- genes_total[genes_total$total>1,]

setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms/Gene_Enrichment")
genes_Hugo_squamous <- genes_top$Hugo_Symbol
## Convert gene symbols to Entrez IDs in R
### https://www.gungorbudak.com/blog/2018/08/07/convert-gene-symbols-to-entrez-ids-in-r/
genes_Entrez_squamous <- as.list(mapIds(org.Hs.eg.db, genes_Hugo_squamous, 'ENTREZID', 'SYMBOL'))

## Generate a enrichment with data from KEGG information
KEGG_genes_squamous <- enrichKEGG(gene = genes_Entrez_squamous, organism = 'hsa', pvalueCutoff = 0.05)
head(KEGG_genes_squamous)
KEGG_squamous_results <- KEGG_genes_squamous@result
library('writexl')
write_xlsx(KEGG_squamous_results, 'Gene_enrichment_data_squamous.xlsx')
dotplot(KEGG_genes_squamous, showCategory = 10)
ggsave('KEGG_dotplot_squamous.pdf', width = 6, height = 5, device = 'pdf')

## Relation between enriched pathways
library(enrichplot)
pairwise_relation_squamous <- pairwise_termsim(KEGG_genes_squamous)
enrichplot::emapplot(pairwise_relation_squamous, showCategory = 20)
enrichplot::emapplot(pairwise_relation_squamous, showCategory = c('Human papillomavirus infection', 'Focal adhesion', 'ECM-receptor interaction',
                                          'Regulation of actin cytoskeleton', 'Proteoglycans in cancer', 'cGMP-PKG signaling pathway',
                                          'Platelet activation', 'Growth hormone synthesis, secretion and action', 'Glutamatergic synapse',
                                          'cAMP signaling pathway', 'Calcium signaling pathway', 'Dilated cardiomyopathy', 'Arrhythmogenic right ventricular cardiomyopathy'))
ggsave('KEGG_enrichplot_pathways_squamous.pdf', width = 12, height = 9, device = 'pdf')


setwd('~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas/')
Genes <- read_excel('MAF_Adenomas_gene_summary.xlsx') 
genes_total <- Genes[,c(1,11)]
genes_total <- genes_total[order(genes_total$total),]
genes_total$order <- 1:2279
genes_total <- genes_total[genes_total$total>0,]
quantile(genes_total$total, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
# selecting the 10th percentile as low expression genes
genes_top <- genes_total[genes_total$total>1,]
setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas/Gene_Enrichment")
genes_Hugo_adenomas <- genes_top$Hugo_Symbol
## Convert gene symbols to Entrez IDs in R
genes_Entrez_adenomas <- as.list(mapIds(org.Hs.eg.db, genes_Hugo_adenomas, 'ENTREZID', 'SYMBOL'))
## Generate a enrichment with data from KEGG information
KEGG_genes_adenomas <- enrichKEGG(gene = genes_Entrez_adenomas, organism = 'hsa', pvalueCutoff = 0.5)
head(KEGG_genes_adenomas)
KEGG_adenomas_results <- KEGG_genes_adenomas@result
write_xlsx(KEGG_adenomas_results, 'Gene_enrichment_data_adenomas.xlsx')
dotplot(KEGG_genes_adenomas, showCategory = 10)
## cnetplot(KEGG_genes_adenomas)
ggsave('KEGG_dotplot_adenomas.pdf', width = 6, height = 5, device = 'pdf')
## Relation between enriched pathways
pairwise_relation_adenomas <- pairwise_termsim(KEGG_genes_adenomas)
enrichplot::emapplot(pairwise_relation_adenomas, showCategory = 20)
enrichplot::emapplot(pairwise_relation_adenomas, showCategory = c('Thyroid hormone signaling pathway', 'GnRH secretion',
                                                                  'Aldosterone-regulated sodium reabsorption', 'Endometrial cancer',
                                                                  'mTOR signaling pathway', 'Proteoglycans in cancer',
                                                                  'ErbB signaling pathway', 'EGFR tyrosine kinase inhibitor resistance',
                                                                  'Focal adhesion', 'Long-term potentiation', 'Regulation of actin cytoskeleton',
                                                                  'Growth hormone synthesis, secretion and action',
                                                                  'B cell receptor signaling pathway', 'Neurotrophin signaling pathway',
                                                                  'Carbohydrate digestion and absorption', 'Chemokine signaling pathway',
                                                                  'Inositol phosphate metabolism'))
ggsave('KEGG_enrichplot_pathways_adenomas.pdf', width = 6, height = 5, device = 'pdf')


setwd('~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms/')
Genes <- read_excel('MAF_Cystic_gene_summary.xlsx') 
genes_total <- Genes[,c(1,10)]
genes_total <- genes_total[order(genes_total$total),]
genes_total$order <- 1:1116
genes_total <- genes_total[genes_total$total>0,]
quantile(genes_total$total, probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
# selecting the 10th percentile as low expression genes
genes_top <- genes_total[genes_total$total>1,]
setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms/Gene_Enrichment")
genes_Hugo_cystic <- genes_top$Hugo_Symbol
## Convert gene symbols to Entrez IDs in R
genes_Entrez_cystic <- as.list(mapIds(org.Hs.eg.db, genes_Hugo_cystic, 'ENTREZID', 'SYMBOL'))
## Generate a enrichment with data from KEGG information
KEGG_genes_cystic <- enrichKEGG(gene = genes_Entrez_cystic, organism = 'hsa', pvalueCutoff = 0.05)
head(KEGG_genes_cystic)
KEGG_cystic_results <- KEGG_genes_cystic@result
write_xlsx(KEGG_cystic_results, 'Gene_enrichment_data_cystic.xlsx')
dotplot(KEGG_genes_cystic, showCategory = 10)
ggsave('KEGG_dotplot_cystic.pdf', width = 6, height = 5, device = 'pdf')
## Relation between enriched pathways
pairwise_relation_cystic <- pairwise_termsim(KEGG_genes_cystic)
enrichplot::emapplot(pairwise_relation_cystic, showCategory = 20)
ggsave('KEGG_enrichplot_pathways_cystic.pdf', width = 15, height = 12, device = 'pdf')

#######################################################################################################

# Survival Analysis 
## https://www.youtube.com/watch?v=DnygUTAZFmM
library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

## Getting clinical data for TCGA-CESC cohort
clinical_cesc <- GDCquery_clinic('TCGA-CESC')

any(colnames(clinical_cesc) %in% c('vital_status','days_to_last_follow_up','days_to_death'))
which(colnames(clinical_cesc) %in% c('vital_status','days_to_last_follow_up','days_to_death'))
clinical_cesc[,c(9,38,44)]

### Days_to_death: number of days passed from the initial diagnosis to the patient's death
### Days_to_last_follow_up: number of days passed from the initial diagnosis to the last visit

## Looking at some variables associated with survival
table(clinical_cesc$vital_status) #alive: 235 and death:72
table(clinical_cesc$primary_diagnosis)

## Change certain values the way they are encoded 
clinical_cesc$deceased <- ifelse(clinical_cesc$vital_status == "Alive", FALSE, TRUE)
table(clinical_cesc$deceased) # FALSE (alive): 235 & TRUE (death): 72

## Create an "overall survival" variable that is equal to days_to_death
## for dead patients, and to days_to_last_follow_up for patients who
## are still alive
clinical_cesc$overall_survival <- ifelse(clinical_cesc$vital_status == "Alive",
                                         clinical_cesc$days_to_last_follow_up,
                                         clinical_cesc$days_to_death)
table(clinical_cesc$overall_survival)
setwd("~/Escritorio/R/TCGA_2")
write_xlsx(clinical_cesc, 'Clinical_Data_Information_to_Survival.xlsx')

## Get different cervix cancer types sample barcodes
library(readxl)
setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms")
squamous_list_1 <- read_excel('Cases_ID_List.xlsx')
squamous_vector <- squamous_list_1$`Case ID`
squamous_clinical <- subset(clinical_cesc, clinical_cesc$submitter_id %in% squamous_vector)
write_xlsx(squamous_clinical, 'Clinical_Data_Squamous.xlsx')

setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas")
adenomas_list_2 <- read_excel('Cases_ID_List.xlsx')
adenomas_vector <- adenomas_list_2$`Case ID`
adenomas_clinical <- subset(clinical_cesc, clinical_cesc$submitter_id %in% adenomas_vector)
write_xlsx(adenomas_clinical, 'Clinical_Data_Adenomas.xlsx')

setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms")
cystic_list_3 <- read_excel('Cases_ID_List.xlsx')
cystic_vector <- cystic_list_3$`Case ID`
cystic_clinical <- subset(clinical_cesc, clinical_cesc$submitter_id %in% cystic_vector)
write_xlsx(cystic_clinical, 'Clinical_Data_Cystic.xlsx')


## Get gene expression data from squamous cell neoplasms
setwd("~/Escritorio/R/TCGA_2/Squamous_cell_neoplasms")
query_cesc_squamous = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = squamous_vector)
output_cesc_squamous <- getResults(query_cesc_squamous)
### Download data
GDCdownload(query_cesc_squamous)

### Get counts
tcga_cesc_squamous_data <- GDCprepare(query_cesc_squamous, summarizedExperiment = TRUE)
cesc_matrix_squamous <- assay(tcga_cesc_squamous_data, "unstranded")
cesc_matrix_squamous_all <- as.data.frame(cesc_matrix_squamous)
cesc_matrix_squamous_all <- tibble::rownames_to_column(cesc_matrix_squamous_all, 'Genes.ID')
library(writexl)
write_xlsx(cesc_matrix_squamous_all, 'counts_squamous.xlsx')

## Extract gene and sample metadata from SummarizedExperiment object
gene_squamous_metadata <- as.data.frame(rowData(tcga_cesc_squamous_data))
write_xlsx(gene_squamous_metadata, 'gene_squamous_metadata.xlsx')
coldata <- as.data.frame(colData(tcga_cesc_squamous_data))
write_xlsx(coldata, 'coldata_squamous_metadata.xlsx')

## vst transform counts to be used in survival analysis
### Setting up countData object
dds_squamous <- DESeqDataSetFromMatrix(countData = cesc_matrix_squamous,
                                       colData = coldata,
                                       design = ~ 1)

### Removing genes with sum total of 10 reads across all samples
keep_squamous <- rowSums(counts(dds_squamous)) >= 10
dds_squamous <- dds_squamous[keep_squamous,]

### vst
vsd_squamous <- vst(dds_squamous, blind = FALSE)
cesc_matrix_vst_squamous <- assay(vsd_squamous)


## Get data for TTN gene and add gene metadata information to it 
cesc_squamous_ttn <- cesc_matrix_vst_squamous %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_squamous_metadata, by = "gene_id") %>%
  filter(gene_name == "TTN")

### Get median value
median_value_ttn_cesc <- median(cesc_squamous_ttn$counts)

### Denotate which cases have higher or lower expression than median count
cesc_squamous_ttn$strata <- ifelse(cesc_squamous_ttn$counts >= median_value_ttn_cesc, "HIGH", "LOW")

## Add clinical information to cesc_squamous_ttn
cesc_squamous_ttn$case_id <- gsub('-01.*', '', cesc_squamous_ttn$case_id)
cesc_squamous_ttn <- merge(cesc_squamous_ttn, squamous_clinical, by.x = 'case_id', by.y = 'submitter_id')

# Fitting survival curve
fit_ttn_squamous <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_squamous_ttn)
fit_ttn_squamous
ggsurvplot(fit_ttn_squamous,
           data = cesc_squamous_ttn,
           pval = T,
           risk.table = T)
# Save as 'KM_TTN.pdf' (Landscape 7 x 6)

## Get data for PIK3CA gene and add gene metadata information to it 
cesc_squamous_pik3ca <- cesc_matrix_vst_squamous %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_squamous_metadata, by = "gene_id") %>%
  filter(gene_name == "PIK3CA")
### Get median value
median_value_pik3ca_cesc <- median(cesc_squamous_pik3ca$counts)
### Denotate which cases have higher or lower expression than median count
cesc_squamous_pik3ca$strata <- ifelse(cesc_squamous_pik3ca$counts >= median_value_pik3ca_cesc, "HIGH", "LOW")
## Add clinical information to cesc_squamous_ttn
cesc_squamous_pik3ca$case_id <- gsub('-01.*', '', cesc_squamous_pik3ca$case_id)
cesc_squamous_pik3ca <- merge(cesc_squamous_pik3ca, squamous_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_pik3ca_squamous <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_squamous_pik3ca)
fit_pik3ca_squamous
ggsurvplot(fit_pik3ca_squamous,
           data = cesc_squamous_pik3ca,
           pval = T,
           risk.table = T)
# Save as 'KM_PIK3CA.pdf' (Landscape 7 x 6)

## Get data for KMT2C gene and add gene metadata information to it 
cesc_squamous_kmt2c <- cesc_matrix_vst_squamous %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_squamous_metadata, by = "gene_id") %>%
  filter(gene_name == "KMT2C")
### Get median value
median_value_kmt2c_cesc <- median(cesc_squamous_kmt2c$counts)
### Denotate which cases have higher or lower expression than median count
cesc_squamous_kmt2c$strata <- ifelse(cesc_squamous_kmt2c$counts >= median_value_kmt2c_cesc, "HIGH", "LOW")
## Add clinical information to cesc_squamous_kmt2c
cesc_squamous_kmt2c$case_id <- gsub('-01.*', '', cesc_squamous_kmt2c$case_id)
cesc_squamous_kmt2c <- merge(cesc_squamous_kmt2c, squamous_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_kmt2c_squamous <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_squamous_kmt2c)
fit_kmt2c_squamous
ggsurvplot(fit_kmt2c_squamous,
           data = cesc_squamous_kmt2c,
           pval = T,
           risk.table = T)
# Save as 'KM_KMT2C.pdf'(Landscape 7 x 6)


## Survival analysis with TMB from squamous samples
squamous_TMB <- read_excel('MAF_Squamous_sample_summary.xlsx')
### Get median value
median_value_tmb_squamous <- median(squamous_TMB$total)
### Denotate which cases have higher or lower expression than median count
squamous_TMB$strata <- ifelse(squamous_TMB$total >= median_value_tmb_squamous, "HIGH", "LOW")
## Add clinical information
squamous_TMB$case_id <- gsub('-01.*', '', squamous_TMB$Tumor_Sample_Barcode)
squamous_TMB <- merge(squamous_TMB, squamous_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_squamous <- survfit(Surv(overall_survival, deceased) ~ strata, data = squamous_TMB)
fit_squamous
ggsurvplot(fit_squamous,
           data =  squamous_TMB,
           palette = 'aaas',
           pval = T,
           risk.table = T,
           ggtheme = theme_bw(),
           xlab = 'Time (days)',
           ylab = 'Survival rate',
           legend.labs = c('High TMB','Low TMB'),
           legend.title = '') 
# Save as 'survival_rate_TMB_squamous.pdf'(Landscape 7 x 6)

## Boxplots TMB
### AGE
### Denotate which cases have higher or lower age
squamous_TMB$age_criteria <- ifelse(squamous_TMB$age_at_index > 60, ">60", "<=60")
library(rstatix)
stat.test.age <- squamous_TMB %>% wilcox_test(total ~ age_criteria, p.adjust.method = "none")
stat.test.age
age <- ggboxplot(squamous_TMB, x = "age_criteria", y = "total", fill= "age_criteria", palette = "Set1", xlab = "age", ylab = 'TMB', legend.title = 'age') +
  stat_pvalue_manual(stat.test.age, y.position = c(1810), tip.length = 0.02) + theme_bw()
age
ggsave("age_squamous.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-T
### Denotate which cases have higher or lower age
squamous_TMB$AJCC_T[squamous_TMB$ajcc_pathologic_t == 'T1a1' | squamous_TMB$ajcc_pathologic_t == 'T1b' |squamous_TMB$ajcc_pathologic_t == 'T1b1' |squamous_TMB$ajcc_pathologic_t == 'T1b2' |squamous_TMB$ajcc_pathologic_t == 'T2' |squamous_TMB$ajcc_pathologic_t == 'T2a' |squamous_TMB$ajcc_pathologic_t == 'T2a1' | squamous_TMB$ajcc_pathologic_t == 'T2a2' |squamous_TMB$ajcc_pathologic_t == 'T2b'] <- "T1-2"
squamous_TMB$AJCC_T[squamous_TMB$ajcc_pathologic_t == 'T3a' | squamous_TMB$ajcc_pathologic_t == 'T3b' |squamous_TMB$ajcc_pathologic_t == 'T4'] <- "T3-4"
squamous_TMB <- squamous_TMB[!is.na(squamous_TMB$AJCC_T),]
stat.test.T <- squamous_TMB %>% wilcox_test(total ~ AJCC_T, p.adjust.method = "none")
stat.test.T
AJCC_T <- ggboxplot(squamous_TMB, x = "AJCC_T", y = "total", fill= "AJCC_T", palette = "Set1", xlab = "AJCC-T", ylab = 'TMB', legend.title = 'AJCC-T') +
  stat_pvalue_manual(stat.test.T, y.position = c(1810), tip.length = 0.02) + theme_bw()
AJCC_T
ggsave("AJCC-T_squamous.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-N
### First load again the squamous_TMB matrix
### Denotate which cases have higher or lower age
squamous_TMB$AJCC_N[squamous_TMB$ajcc_pathologic_n == 'N0'] <- "N0"
squamous_TMB$AJCC_N[squamous_TMB$ajcc_pathologic_n == 'N1'] <- "N1"
squamous_TMB <- squamous_TMB[!is.na(squamous_TMB$AJCC_N),]
stat.test.N <- squamous_TMB %>% wilcox_test(total ~ AJCC_N, p.adjust.method = "none")
stat.test.N
AJCC_N <- ggboxplot(squamous_TMB, x = "AJCC_N", y = "total", fill= "AJCC_N", palette = "Set1", xlab = "AJCC-N", ylab = 'TMB', legend.title = 'AJCC-N') +
  stat_pvalue_manual(stat.test.N, y.position = c(1810), tip.length = 0.02) + theme_bw()
AJCC_N
ggsave("AJCC-N_squamous.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-M
### First load again the squamous_TMB matrix
### Denotate which cases have higher or lower age
squamous_TMB$AJCC_M[squamous_TMB$ajcc_pathologic_m == 'M0'] <- "M0"
squamous_TMB$AJCC_M[squamous_TMB$ajcc_pathologic_m == 'M1'] <- "M1"
squamous_TMB <- squamous_TMB[!is.na(squamous_TMB$AJCC_M),]
stat.test.M <- squamous_TMB %>% wilcox_test(total ~ AJCC_M, p.adjust.method = "none")
stat.test.M
AJCC_M <- ggboxplot(squamous_TMB, x = "AJCC_M", y = "total", fill= "AJCC_M", palette = "Set1", xlab = "AJCC-M", ylab = 'TMB', legend.title = 'AJCC-M') +
  stat_pvalue_manual(stat.test.M, y.position = c(1810), tip.length = 0.02) + theme_bw()
AJCC_M
ggsave("AJCC-M_squamous.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")


## Get gene expression data from adenomas and adenocarcinomas
setwd("~/Escritorio/R/TCGA_2/Adenomas_and_adenocarcinomas")
query_cesc_adenomas = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = adenomas_vector)
output_cesc_adenomas <- getResults(query_cesc_adenomas)

### Download data
GDCdownload(query_cesc_adenomas)

### Get counts
tcga_cesc_adenomas_data <- GDCprepare(query_cesc_adenomas, summarizedExperiment = TRUE)
cesc_matrix_adenomas <- assay(tcga_cesc_adenomas_data, "unstranded")

## Extract gene and sample metadata from SummarizedExperiment object
gene_adenomas_metadata <- as.data.frame(rowData(tcga_cesc_adenomas_data))
coldata <- as.data.frame(colData(tcga_cesc_adenomas_data))

## vst transform counts to be used in survival analysis
### Setting up countData object
dds_adenomas <- DESeqDataSetFromMatrix(countData = cesc_matrix_adenomas,
                                       colData = coldata,
                                       design = ~ 1)

### Removing genes with sum total of 10 reads across all samples
keep_adenomas <- rowSums(counts(dds_adenomas)) >= 10
dds_adenomas <- dds_adenomas[keep_adenomas,]

### vst
vsd_adenomas <- vst(dds_adenomas, blind = FALSE)
cesc_matrix_vst_adenomas <- assay(vsd_adenomas)

## Get data for TTN gene and add gene metadata information to it 
cesc_adenomas_ttn <- cesc_matrix_vst_adenomas %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_adenomas_metadata, by = "gene_id") %>%
  filter(gene_name == "TTN")

### Get median value
median_value_ttn_cesc <- median(cesc_adenomas_ttn$counts)

### Denotate which cases have higher or lower expression than median count
cesc_adenomas_ttn$strata <- ifelse(cesc_adenomas_ttn$counts >= median_value_ttn_cesc, "HIGH", "LOW")

## Add clinical information to cesc_adenomas_ttn
cesc_adenomas_ttn$case_id <- gsub('-01.*', '', cesc_adenomas_ttn$case_id)
cesc_adenomas_ttn <- merge(cesc_adenomas_ttn, adenomas_clinical, by.x = 'case_id', by.y = 'submitter_id')

# Fitting survival curve
fit_ttn_adenomas <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_adenomas_ttn)
fit_ttn_adenomas
ggsurvplot(fit_ttn_adenomas,
           data = cesc_adenomas_ttn,
           pval = T,
           risk.table = T)
# Save as 'KM_TTN.pdf' (Landscape 7 x 6)

## Get data for PIK3CA gene and add gene metadata information to it 
cesc_adenomas_pik3ca <- cesc_matrix_vst_adenomas %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_adenomas_metadata, by = "gene_id") %>%
  filter(gene_name == "PIK3CA")
### Get median value
median_value_pik3ca_cesc <- median(cesc_adenomas_pik3ca$counts)
### Denotate which cases have higher or lower expression than median count
cesc_adenomas_pik3ca$strata <- ifelse(cesc_adenomas_pik3ca$counts >= median_value_pik3ca_cesc, "HIGH", "LOW")
## Add clinical information to cesc_adenomas_pik3ca
cesc_adenomas_pik3ca$case_id <- gsub('-01.*', '', cesc_adenomas_pik3ca$case_id)
cesc_adenomas_pik3ca <- merge(cesc_adenomas_pik3ca, adenomas_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_pik3ca_adenomas <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_adenomas_pik3ca)
fit_pik3ca_adenomas
ggsurvplot(fit_pik3ca_adenomas,
           data = cesc_adenomas_pik3ca,
           pval = T,
           risk.table = T)
# Save as 'KM_PIK3CA.pdf' (Landscape 7 x 6)

## Get data for ERBB2 gene and add gene metadata information to it 
cesc_adenomas_erbb2 <- cesc_matrix_vst_adenomas %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_adenomas_metadata, by = "gene_id") %>%
  filter(gene_name == "ERBB2")
### Get median value
median_value_erbb2_cesc <- median(cesc_adenomas_erbb2$counts)
### Denotate which cases have higher or lower expression than median count
cesc_adenomas_erbb2$strata <- ifelse(cesc_adenomas_erbb2$counts >= median_value_erbb2_cesc, "HIGH", "LOW")
## Add clinical information to cesc_adenomas_erbb2
cesc_adenomas_erbb2$case_id <- gsub('-01.*', '', cesc_adenomas_erbb2$case_id)
cesc_adenomas_erbb2 <- merge(cesc_adenomas_erbb2, adenomas_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_erbb2_adenomas <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_adenomas_erbb2)
fit_erbb2_adenomas
ggsurvplot(fit_erbb2_adenomas,
           data = cesc_adenomas_erbb2,
           pval = T,
           risk.table = T)
# Save as 'KM_ERBB2.pdf' (Landscape 7 x 6)

## Survival analysis with TMB from adenomas and adenocarcinomas samples
adenomas_TMB <- read_excel('MAF_Adenomas_sample_summary.xlsx')
### Get median value
median_value_tmb_adenomas <- median(adenomas_TMB$total)
### Denotate which cases have higher or lower expression than median count
adenomas_TMB$strata <- ifelse(adenomas_TMB$total >= median_value_tmb_adenomas, "HIGH", "LOW")
## Add clinical information
adenomas_TMB$case_id <- gsub('-01.*', '', adenomas_TMB$Tumor_Sample_Barcode)
adenomas_TMB <- merge(adenomas_TMB, adenomas_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_adenomas <- survfit(Surv(overall_survival, deceased) ~ strata, data = adenomas_TMB)
fit_adenomas
ggsurvplot(fit_adenomas,
           data =  adenomas_TMB,
           palette = 'aaas',
           pval = T,
           risk.table = T,
           ggtheme = theme_bw(),
           xlab = 'Time (days)',
           ylab = 'Survival rate',
           legend.labs = c('High TMB','Low TMB'),
           legend.title = '')
# Save as 'survival_rate_TMB_adenomas.pdf'(Landscape 7 x 6)

## Boxplots TMB
### AGE
### Denotate which cases have higher or lower age
adenomas_TMB$age_criteria <- ifelse(adenomas_TMB$age_at_index > 60, ">60", "<=60")
library(rstatix)
stat.test.age <- adenomas_TMB %>% wilcox_test(total ~ age_criteria, p.adjust.method = "none")
stat.test.age
age <- ggboxplot(adenomas_TMB, x = "age_criteria", y = "total", fill= "age_criteria", palette = "Set1", xlab = "age", ylab = 'TMB', legend.title = 'age') +
  stat_pvalue_manual(stat.test.age, y.position = c(1050), tip.length = 0.02) + theme_bw()
age
ggsave("age_adenomas.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-T
### Denotate which cases have higher or lower age
adenomas_TMB$AJCC_T[adenomas_TMB$ajcc_pathologic_t == 'T1a1' | adenomas_TMB$ajcc_pathologic_t == 'T1b' | adenomas_TMB$ajcc_pathologic_t == 'T1b1' | adenomas_TMB$ajcc_pathologic_t == 'T1b2' | adenomas_TMB$ajcc_pathologic_t == 'T2' | adenomas_TMB$ajcc_pathologic_t == 'T2a' | adenomas_TMB$ajcc_pathologic_t == 'T2a1' | adenomas_TMB$ajcc_pathologic_t == 'T2a2' | adenomas_TMB$ajcc_pathologic_t == 'T2b'] <- "T1-2"
adenomas_TMB$AJCC_T[adenomas_TMB$ajcc_pathologic_t == 'T3a' | adenomas_TMB$ajcc_pathologic_t == 'T3b' | adenomas_TMB$ajcc_pathologic_t == 'T4'] <- "T3-4"
adenomas_TMB <- adenomas_TMB[!is.na(adenomas_TMB$AJCC_T),]
stat.test.T <- adenomas_TMB %>% wilcox_test(total ~ AJCC_T, p.adjust.method = "none")
stat.test.T
AJCC_T <- ggboxplot(adenomas_TMB, x = "AJCC_T", y = "total", fill= "AJCC_T", palette = "Set1", xlab = "AJCC-T", ylab = 'TMB', legend.title = 'AJCC-T') +
  stat_pvalue_manual(stat.test.T, y.position = c(1050), tip.length = 0.02) + theme_bw()
AJCC_T
ggsave("AJCC-T_adenomas.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-N
### First load again the adenomas_TMB matrix
### Denotate which cases have higher or lower age
adenomas_TMB$AJCC_N[adenomas_TMB$ajcc_pathologic_n == 'N0'] <- "N0"
adenomas_TMB$AJCC_N[adenomas_TMB$ajcc_pathologic_n == 'N1'] <- "N1"
adenomas_TMB <- adenomas_TMB[!is.na(adenomas_TMB$AJCC_N),]
stat.test.N <- adenomas_TMB %>% wilcox_test(total ~ AJCC_N, p.adjust.method = "none")
stat.test.N
AJCC_N <- ggboxplot(adenomas_TMB, x = "AJCC_N", y = "total", fill= "AJCC_N", palette = "Set1", xlab = "AJCC-N", ylab = 'TMB', legend.title = 'AJCC-N') +
  stat_pvalue_manual(stat.test.N, y.position = c(1050), tip.length = 0.02) + theme_bw()
AJCC_N
ggsave("AJCC-N_adenomas.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-M
### First load again the adenomas_TMB matrix
### Denotate which cases have higher or lower age
adenomas_TMB$AJCC_M[adenomas_TMB$ajcc_pathologic_m == 'M0'] <- "M0"
adenomas_TMB$AJCC_M[adenomas_TMB$ajcc_pathologic_m == 'M1'] <- "M1"
adenomas_TMB <- adenomas_TMB[!is.na(adenomas_TMB$AJCC_M),]
stat.test.M <- adenomas_TMB %>% wilcox_test(total ~ AJCC_M, p.adjust.method = "none")
stat.test.M
AJCC_M <- ggboxplot(adenomas_TMB, x = "AJCC_M", y = "total", fill= "AJCC_M", palette = "Set1", xlab = "AJCC-M", ylab = 'TMB', legend.title = 'AJCC-M') +
  stat_pvalue_manual(stat.test.M, y.position = c(1050), tip.length = 0.02) + theme_bw()
AJCC_M
ggsave("AJCC-M_adenomas.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")



## Get gene expression data from cystic, mucinous and serous neoplasms
setwd("~/Escritorio/R/TCGA_2/Cystic_mucinous_and_serous_neoplasms")
query_cesc_cystic = GDCquery(
  project = "TCGA-CESC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open",
  barcode = cystic_vector)
output_cesc_cystic <- getResults(query_cesc_cystic)
### Download data
GDCdownload(query_cesc_cystic)

### Get counts
tcga_cesc_cystic_data <- GDCprepare(query_cesc_cystic, summarizedExperiment = TRUE)
cesc_matrix_cystic <- assay(tcga_cesc_cystic_data, "unstranded")

## Extract gene and sample metadata from SummarizedExperiment object
gene_cystic_metadata <- as.data.frame(rowData(tcga_cesc_cystic_data))
coldata <- as.data.frame(colData(tcga_cesc_cystic_data))

## vst transform counts to be used in survival analysis
### Setting up countData object
dds_cystic <- DESeqDataSetFromMatrix(countData = cesc_matrix_cystic,
                                       colData = coldata,
                                       design = ~ 1)

### Removing genes with sum total of 10 reads across all samples
keep_cystic <- rowSums(counts(dds_cystic)) >= 10
dds_cystic <- dds_cystic[keep_cystic,]

### vst
vsd_cystic <- vst(dds_cystic, blind = FALSE)
cesc_matrix_vst_cystic <- assay(vsd_cystic)

## Get data for PIK3CA gene and add gene metadata information to it 
cesc_cystic_pik3ca <- cesc_matrix_vst_cystic %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_cystic_metadata, by = "gene_id") %>%
  filter(gene_name == "PIK3CA")

### Get median value
median_value_pik3ca_cesc <- median(cesc_cystic_pik3ca$counts)

### Denotate which cases have higher or lower expression than median count
cesc_cystic_pik3ca$strata <- ifelse(cesc_cystic_pik3ca$counts >= median_value_pik3ca_cesc, "HIGH", "LOW")

## Add clinical information to cesc_cystic_pik3ca
cesc_cystic_pik3ca$case_id <- gsub('-01.*', '', cesc_cystic_pik3ca$case_id)
cesc_cystic_pik3ca <- merge(cesc_cystic_pik3ca, cystic_clinical, by.x = 'case_id', by.y = 'submitter_id')

# Fitting survival curve
fit_pik3ca_cystic <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_cystic_pik3ca)
fit_pik3ca_cystic
ggsurvplot(fit_pik3ca_cystic,
           data = cesc_cystic_pik3ca,
           pval = T,
           risk.table = T)
# Save as 'KM_PIK3CA.pdf' (Landscape 7 x 6)

## Get data for TP53 gene and add gene metadata information to it 
cesc_cystic_tp53 <- cesc_matrix_vst_cystic %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_cystic_metadata, by = "gene_id") %>%
  filter(gene_name == "TP53")
### Get median value
median_value_tp53_cesc <- median(cesc_cystic_tp53$counts)
### Denotate which cases have higher or lower expression than median count
cesc_cystic_tp53$strata <- ifelse(cesc_cystic_tp53$counts >= median_value_tp53_cesc, "HIGH", "LOW")
## Add clinical information to cesc_cystic_tp53
cesc_cystic_tp53$case_id <- gsub('-01.*', '', cesc_cystic_tp53$case_id)
cesc_cystic_tp53 <- merge(cesc_cystic_tp53, cystic_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_tp53_cystic <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_cystic_tp53)
fit_tp53_cystic
ggsurvplot(fit_tp53_cystic,
           data = cesc_cystic_tp53,
           pval = T,
           risk.table = T)
# Save as 'KM_TP53.pdf' (Landscape 7 x 6)

## Get data for LRP1B gene and add gene metadata information to it 
cesc_cystic_lrp1b <- cesc_matrix_vst_cystic %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  gather(key ="case_id", value = "counts", -gene_id) %>%
  left_join(., gene_cystic_metadata, by = "gene_id") %>%
  filter(gene_name == "LRP1B")
### Get median value
median_value_lrp1b_cesc <- median(cesc_cystic_lrp1b$counts)
### Denotate which cases have higher or lower expression than median count
cesc_cystic_lrp1b$strata <- ifelse(cesc_cystic_lrp1b$counts >= median_value_lrp1b_cesc, "HIGH", "LOW")
## Add clinical information to cesc_cystic_lrp1b
cesc_cystic_lrp1b$case_id <- gsub('-01.*', '', cesc_cystic_lrp1b$case_id)
cesc_cystic_lrp1b <- merge(cesc_cystic_lrp1b, cystic_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_lrp1b_cystic <- survfit(Surv(overall_survival, deceased) ~ strata, data = cesc_cystic_lrp1b)
fit_lrp1b_cystic
ggsurvplot(fit_lrp1b_cystic,
           data = cesc_cystic_lrp1b,
           pval = T,
           risk.table = T)
# Save as 'KM_LRP1B.pdf' (Landscape 7 x 6)

## Survival analysis with TMB from cystic, mucinous and serous samples
cystic_TMB <- read_excel('MAF_Cystic_sample_summary.xlsx')
### Get median value
median_value_tmb_cystic <- median(cystic_TMB$total)
### Denotate which cases have higher or lower expression than median count
cystic_TMB$strata <- ifelse(cystic_TMB$total >= median_value_tmb_cystic, "HIGH", "LOW")
## Add clinical information
cystic_TMB$case_id <- gsub('-01.*', '', cystic_TMB$Tumor_Sample_Barcode)
cystic_TMB <- merge(cystic_TMB, cystic_clinical, by.x = 'case_id', by.y = 'submitter_id')
# Fitting survival curve
fit_cystic <- survfit(Surv(overall_survival, deceased) ~ strata, data = cystic_TMB)
fit_cystic
ggsurvplot(fit_cystic,
           data =  cystic_TMB,
           palette = 'aaas',
           pval = T,
           risk.table = T,
           ggtheme = theme_bw(),
           xlab = 'Time (days)',
           ylab = 'Survival rate',
           legend.labs = c('High TMB','Low TMB'),
           legend.title = '') 
# Save as 'survival_rate_TMB_cystic.pdf'(Landscape 7 x 6)

## Boxplots TMB
### AGE
### Denotate which cases have higher or lower age
cystic_TMB$age_criteria <- ifelse(cystic_TMB$age_at_index > 60, ">60", "<=60")
library(rstatix)
stat.test.age <- cystic_TMB %>% wilcox_test(total ~ age_criteria, p.adjust.method = "none")
stat.test.age
age <- ggboxplot(cystic_TMB, x = "age_criteria", y = "total", fill= "age_criteria", palette = "Set1", xlab = "age", ylab = 'TMB', legend.title = 'age') +
  stat_pvalue_manual(stat.test.age, y.position = c(160), tip.length = 0.02) + theme_bw()
age
ggsave("age_cystic.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-T
### Denotate which cases have higher or lower age
cystic_TMB$AJCC_T[cystic_TMB$ajcc_pathologic_t == 'T1a1' | cystic_TMB$ajcc_pathologic_t == 'T1b' | cystic_TMB$ajcc_pathologic_t == 'T1b1' | cystic_TMB$ajcc_pathologic_t == 'T1b2' | cystic_TMB$ajcc_pathologic_t == 'T2' | cystic_TMB$ajcc_pathologic_t == 'T2a' | cystic_TMB$ajcc_pathologic_t == 'T2a1' | cystic_TMB$ajcc_pathologic_t == 'T2a2' | cystic_TMB$ajcc_pathologic_t == 'T2b'] <- "T1-2"
cystic_TMB$AJCC_T[cystic_TMB$ajcc_pathologic_t == 'T3a' | cystic_TMB$ajcc_pathologic_t == 'T3b' | cystic_TMB$ajcc_pathologic_t == 'T4'] <- "T3-4"
cystic_TMB <- cystic_TMB[!is.na(cystic_TMB$AJCC_T),]
stat.test.T <- cystic_TMB %>% wilcox_test(total ~ AJCC_T, p.adjust.method = "none")
stat.test.T
AJCC_T <- ggboxplot(cystic_TMB, x = "AJCC_T", y = "total", fill= "AJCC_T", palette = "Set1", xlab = "AJCC-T", ylab = 'TMB', legend.title = 'AJCC-T') +
  stat_pvalue_manual(stat.test.T, y.position = c(160), tip.length = 0.02) + theme_bw()
AJCC_T
ggsave("AJCC-T_cystic.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-N
### First load again the cystic_TMB matrix
### Denotate which cases have higher or lower age
cystic_TMB$AJCC_N[cystic_TMB$ajcc_pathologic_n == 'N0'] <- "N0"
cystic_TMB$AJCC_N[cystic_TMB$ajcc_pathologic_n == 'N1'] <- "N1"
cystic_TMB <- cystic_TMB[!is.na(cystic_TMB$AJCC_N),]
stat.test.N <- cystic_TMB %>% wilcox_test(total ~ AJCC_N, p.adjust.method = "none")
stat.test.N
AJCC_N <- ggboxplot(cystic_TMB, x = "AJCC_N", y = "total", fill= "AJCC_N", palette = "Set1", xlab = "AJCC-N", ylab = 'TMB', legend.title = 'AJCC-N') +
  stat_pvalue_manual(stat.test.N, y.position = c(160), tip.length = 0.02) + theme_bw()
AJCC_N
ggsave("AJCC-N_cystic.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")

### AJCC-M
### First load again the adenomas_TMB matrix
### Denotate which cases have higher or lower age
cystic_TMB$AJCC_M[cystic_TMB$ajcc_pathologic_m == 'M0'] <- "M0"
cystic_TMB$AJCC_M[cystic_TMB$ajcc_pathologic_m == 'M1'] <- "M1"
cystic_TMB <- cystic_TMB[!is.na(cystic_TMB$AJCC_M),]
AJCC_M <- ggboxplot(cystic_TMB, x = "AJCC_M", y = "total", fill= "AJCC_M", palette = "Set1", xlab = "AJCC-M", ylab = 'TMB', legend.title = 'AJCC-M') + theme_bw()
AJCC_M
ggsave("AJCC-M_cystic.pdf", height = 5, width = 4, dpi = 400, bg = "white", device = "pdf")



# Visualization of functional profile comparison
## Reference: https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
library(tidyverse)

setwd("~/Escritorio/R/TCGA_2")
KEGG_genes_ALL <- read_excel("Gene_enrichment_data_ALL_GROUPS_BEST_30.xlsx")
KEGG_genes_ALL <- KEGG_genes_ALL[order(KEGG_genes_ALL$Count, decreasing = TRUE),]

KEGG_genes_ALL$Group <- factor(KEGG_genes_ALL$Group, levels = c("SCC","AA","CMSN")) #to order the Groups
KEGG_genes_ALL$Description <- factor(KEGG_genes_ALL$Description, levels = c("Endocrine resistance", "Non-small cell lung cancer",
                                                                            "Viral myocarditis", "Epstein-Barr virus infection", "Hepatocellular carcinoma", "Cellular senescence", "Kaposi sarcoma-associated herpesvirus infection",
                                                                            "Human cytomegalovirus infection", "Carbohydrate digestion and absorption", "Aldosterone-regulated sodium reabsorption","Inositol phosphate metabolism",
                                                                            "Long-term potentiation", "GnRH secretion", "Lysine degradation", "Nucleocytoplasmic transport", "Regulation of actin cytoskeleton", "ABC transporters",
                                                                            "Arrhythmogenic right ventricular cardiomyopathy", "Protein digestion and absorption", "ECM-receptor interaction", "Thyroid hormone signaling pathway",
                                                                            "Axon guidance", "cAMP signaling pathway", "Motor proteins", "Focal adhesion", "Human papillomavirus infection"))

ALL<- ggplot(KEGG_genes_ALL, aes(x= GeneRatio, y= Description, size=Count, color=p.adjust, group=Group)) + 
  geom_point(alpha = 0.8) + scale_color_gradient(low = "#E34F4F", high="blue")
ALL + facet_grid(~Group, scales = "free_x") + ylab("Description")
ggsave("KEGG_dotplot_ALL.pdf", height = 5, width = 11, dpi = 600, bg = "white", device = "pdf")
ggsave("KEGG_dotplot_ALL.png", height = 5, width = 11, dpi = 600, bg = "white", device = "png")
