library(TCGAbiolinks)
library(tidyverse)
library(pheatmap)
library(SummarizedExperiment)
library(sesameData)


# get a list of projects
gdc_projects <- getGDCprojects()
getProjectSummary('TCGA-MESO')

# Build Query
query_TCGA <- GDCquery('TCGA-MESO','Transcriptome Profiling')
output_query_TCGA <- getResults(query_TCGA)


# build a query to retrieve gene expression data ------------
query_TCGA <- GDCquery(project = 'TCGA-MESO',
                       data.category = 'Transcriptome Profiling',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open',
                       barcode = c('TCGA-MQ-A6BR-01A-11R-A34F-07', 'TCGA-UD-AAC7-01A-11R-A40A-07','TCGA-TS-A8AY-01A-21R-A40A-07'))

# download data - GDCdownload
GDCdownload(query_TCGA)


# prepare data
tcga_MESO_data <- GDCprepare(query_TCGA, summarizedExperiment = TRUE)
MESO_matrix <- assay(tcga_MESO_data, 'fpkm_unstrand')
view(MESO_matrix)

write.table(MESO_matrix, file = "d:/MESO.txt", sep = ",", quote = FALSE, row.names = T)
