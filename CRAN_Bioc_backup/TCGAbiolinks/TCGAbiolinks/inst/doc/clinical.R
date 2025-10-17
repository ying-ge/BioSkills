## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----BCR_Biotab, results='hide', echo=TRUE, message=FALSE, warning=FALSE------
query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
)
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
clinical.BCRtab.all$clinical_drug_acc  %>% 
    head  %>% 
    DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

## ----results = "hide",cache=TRUE, message=FALSE-------------------------------
library(TCGAbiolinks)
query <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
)

GDCdownload(query)
clinical_tab_all <- GDCprepare(query)

## -----------------------------------------------------------------------------
# All available tables
names(clinical_tab_all)

# columns from clinical_patient
dplyr::glimpse(clinical_tab_all$clinical_patient_acc)

## ----results = "hide",cache=TRUE, message=FALSE,warning=FALSE-----------------
# Biospecimen BCR Biotab
query_biospecimen <- GDCquery(
    project = "TCGA-ACC", 
    data.category = "Biospecimen",
    data.type = "Biospecimen Supplement", 
    data.format = "BCR Biotab"
)
GDCdownload(query_biospecimen)
biospecimen_tab_all <- GDCprepare(query_biospecimen)

## -----------------------------------------------------------------------------
# All available tables
names(biospecimen_tab_all)

biospecimen_tab_all$biospecimen_sample_acc  %>% 
    head  %>% 
    DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

## ----results='hide', echo=TRUE, message=FALSE, warning=FALSE------------------
clinical <- GDCquery_clinic(project = "TCGA-ACC", type = "clinical")

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
clinical %>%
    head %>% 
    DT::datatable(
        filter = 'top', 
        options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
        rownames = FALSE
    )

## ----eval=FALSE,results='hide', echo=TRUE, message=FALSE, warning=FALSE-------
#  clinical_beataml <- GDCquery_clinic(
#      project = "BEATAML1.0-COHORT",
#      type = "clinical"
#  )
#  
#  clinical_cptac2 <- GDCquery_clinic(
#      project = "CPTAC-2",
#      type = "clinical"
#  )
#  
#  clinical_genie <- GDCquery_clinic(
#      project = "GENIE-MSK",
#      type = "clinical"
#  )

## ----results = 'hide',echo=TRUE, message=FALSE, warning=FALSE-----------------
query <- GDCquery(
    project = "TCGA-COAD", 
    data.category = "Clinical", 
    data.format = "bcr xml",
    barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")
)
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
clinical %>% 
    datatable(filter = 'top', 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
              rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE----------------
clinical.drug <- GDCprepare_clinic(query, clinical.info = "drug")

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
clinical.drug %>% 
    datatable(filter = 'top', 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
              rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE----------------
clinical.radiation <- GDCprepare_clinic(query, clinical.info = "radiation")

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
clinical.radiation %>% 
    datatable(filter = 'top', 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
              rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE----------------
clinical.admin <- GDCprepare_clinic(query, clinical.info = "admin")

## ----echo = TRUE, message = FALSE, warning = FALSE----------------------------
clinical.admin %>% 
    datatable(filter = 'top', 
              options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),  
              rownames = FALSE)

## ----results = 'hide', echo=TRUE, message=FALSE, warning=FALSE----------------
# Pathology report from harmonized portal 
query_harmonized <- GDCquery(
    project = "TCGA-COAD", 
    data.category = "Biospecimen", 
    data.type = "Slide Image",
    experimental.strategy = "Diagnostic Slide",
    barcode = c("TCGA-RU-A8FL","TCGA-AA-3972")
)  

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
query_harmonized  %>% 
    getResults %>% 
    head  %>% 
    DT::datatable(options = list(scrollX = TRUE, keys = TRUE))

## ----eval = TRUE--------------------------------------------------------------
bar <- c(
    "TCGA-G9-6378-02A-11R-1789-07", "TCGA-CH-5767-04A-11R-1789-07",  
    "TCGA-G9-6332-60A-11R-1789-07", "TCGA-G9-6336-01A-11R-1789-07",
    "TCGA-G9-6336-11A-11R-1789-07", "TCGA-G9-7336-11A-11R-1789-07",
    "TCGA-G9-7336-04A-11R-1789-07", "TCGA-G9-7336-14A-11R-1789-07",
    "TCGA-G9-7036-04A-11R-1789-07", "TCGA-G9-7036-02A-11R-1789-07",
    "TCGA-G9-7036-11A-11R-1789-07", "TCGA-G9-7036-03A-11R-1789-07",
    "TCGA-G9-7036-10A-11R-1789-07", "TCGA-BH-A1ES-10A-11R-1789-07",
    "TCGA-BH-A1F0-10A-11R-1789-07", "TCGA-BH-A0BZ-02A-11R-1789-07",
    "TCGA-B6-A0WY-04A-11R-1789-07", "TCGA-BH-A1FG-04A-11R-1789-08",
    "TCGA-D8-A1JS-04A-11R-2089-08", "TCGA-AN-A0FN-11A-11R-8789-08",
    "TCGA-AR-A2LQ-12A-11R-8799-08", "TCGA-AR-A2LH-03A-11R-1789-07",
    "TCGA-BH-A1F8-04A-11R-5789-07", "TCGA-AR-A24T-04A-55R-1789-07",
    "TCGA-AO-A0J5-05A-11R-1789-07", "TCGA-BH-A0B4-11A-12R-1789-07",
    "TCGA-B6-A1KN-60A-13R-1789-07", "TCGA-AO-A0J5-01A-11R-1789-07",
    "TCGA-AO-A0J5-01A-11R-1789-07", "TCGA-G9-6336-11A-11R-1789-07",
    "TCGA-G9-6380-11A-11R-1789-07", "TCGA-G9-6380-01A-11R-1789-07",
    "TCGA-G9-6340-01A-11R-1789-07", "TCGA-G9-6340-11A-11R-1789-07"
)

S <- TCGAquery_SampleTypes(bar,"TP")
S2 <- TCGAquery_SampleTypes(bar,"NB")

# Retrieve multiple tissue types  NOT FROM THE SAME PATIENTS
SS <- TCGAquery_SampleTypes(bar,c("TP","NB"))

# Retrieve multiple tissue types  FROM THE SAME PATIENTS
SSS <- TCGAquery_MatchedCoupledSampleTypes(bar,c("NT","TP"))

## ----eval = FALSE-------------------------------------------------------------
#  # This code will get all clinical indexed data from TCGA
#  library(data.table)
#  library(dplyr)
#  library(regexPipes)
#  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>%
#      regexPipes::grep("TCGA",value = TRUE) %>%
#      sort %>%
#      plyr::alply(1,GDCquery_clinic, .progress = "text") %>%
#      rbindlist
#  readr::write_csv(clinical,path = paste0("all_clin_indexed.csv"))
#  
#  # This code will get all clinical XML data from TCGA
#  getclinical <- function(proj){
#      message(proj)
#      while(1){
#          result = tryCatch({
#              query <- GDCquery(project = proj, data.category = "Clinical",data.format = "bcr xml")
#              GDCdownload(query)
#              clinical <- GDCprepare_clinic(query, clinical.info = "patient")
#              for(i in c("admin","radiation","follow_up","drug","new_tumor_event")){
#                  message(i)
#                  aux <- GDCprepare_clinic(query, clinical.info = i)
#                  if(is.null(aux) || nrow(aux) == 0) next
#                  # add suffix manually if it already exists
#                  replicated <- which(grep("bcr_patient_barcode",colnames(aux), value = T,invert = T) %in% colnames(clinical))
#                  colnames(aux)[replicated] <- paste0(colnames(aux)[replicated],".",i)
#                  if(!is.null(aux)) clinical <- merge(clinical,aux,by = "bcr_patient_barcode", all = TRUE)
#              }
#              readr::write_csv(clinical,path = paste0(proj,"_clinical_from_XML.csv")) # Save the clinical data into a csv file
#              return(clinical)
#          }, error = function(e) {
#              message(paste0("Error clinical: ", proj))
#          })
#      }
#  }
#  clinical <- TCGAbiolinks:::getGDCprojects()$project_id %>%
#      regexPipes::grep("TCGA",value=T) %>% sort %>%
#      plyr::alply(1,getclinical, .progress = "text") %>%
#      rbindlist(fill = TRUE) %>% setDF %>% subset(!duplicated(clinical))
#  
#  readr::write_csv(clinical,path = "all_clin_XML.csv")
#  # result: https://drive.google.com/open?id=0B0-8N2fjttG-WWxSVE5MSGpva1U
#  # Obs: this table has multiple lines for each patient, as the patient might have several followups, drug treatments,
#  # new tumor events etc...

