"---
title: Midterm Project
author: Kailin Liu
---"
###SETTING UP
#install and load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
library(BiocManager)

if (!require("TCGAbiolinks", quietly = TRUE)) 
  BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

if (!require("SummarizedExperiment", quietly = TRUE)) 
  BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

if (!require("maftools", quietly = TRUE)) 
  BiocManager::install("maftools")
library(maftools)

if (!require(ggplot2)){ 
  install.packages("ggplot2")
}
library(ggplot2)

if (!require(survival)){ 
  install.packages("survival")
}
library(survival)

if (!require(survminer)){ 
  install.packages("survminer")
}
library(survminer)


# make outputs folder
#first create a midsemester_project_lastname folder through terminal, then run the following
dir.create("/Users/kathykliu/Desktop/qbio_490_kailin/midsemester_project_liu/outputs")
setwd("/Users/kathykliu/Desktop/qbio_490_kailin/midsemester_project_liu/outputs")


# set up rna_query and rna_se
data_query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(data_query)
data <- GDCprepare(data_query)


# setting up and saving clinical data frames
clinical_query <- GDCquery(project="TCGA-BRCA",data.category="Clinical",file.type="xml")
GDCdownload(clinical_query)

clinical <- GDCprepare_clinic(clinical_query,clinical.info="patient")

clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")

#--------------------------------------------------
### EDITING CLINICAL DATA FRAMES

#checking for NAs in menopause_status, should return 0 
sum(is.na(clinical$menopause_status))

#masking menopause data for pre and post menopause status
menopause_mask <- ifelse(clinical$menopause_status=="Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)" |
                           clinical$menopause_status== "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)", 
                         TRUE, FALSE)
clinical <- clinical[menopause_mask,]

#renaming the Pre and Post categories to Pre and Post respectively (shorter names, easier to access)
clinical$menopause_status <-ifelse(clinical$menopause_status=="Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)","Pre","Post")
write.csv(clinical, "analysis_clinical_data.csv", row.names = FALSE)

#--------------------------------------------------
###SETTING UP FOR BASIC & SURVIVAL ANALYSIS PLOTS

#create a data frame titled analysis that we will be working with for basic plots and survival analysis
analysis <- data.frame(clinical$bcr_patient_barcode, 
                       clinical$days_to_death, 
                       clinical$days_to_last_followup,
                       clinical$vital_status,
                       clinical$age_at_initial_pathologic_diagnosis,
                       clinical$menopause_status)
#naming the columns of the analysis dataframe
colnames(analysis) <- c("patient_barcode", "days_to_death", "days_to_last_followup","vital_status", "age_at_initial_pathologic_diagnosis", "menopause_status")

#--------------------------------------------------
#CREATING plots to visualize relationship between age and menopause status

#line below will save the boxplot to the outputs folder
jpeg("menopause_age_boxplot.jpg")
#creating the box plot with age as the y axis and menopause status as the x-axis
menopause_age_boxplot <- boxplot(formula=analysis$age_at_initial_pathologic_diagnosis ~ analysis$menopause_status,
                               data=clinical,
                               xlab ="Menopause Status",
                               ylab = "Age",
                               main= "Age vs Menopause Status",
                               cex.axis=0.5)
dev.off()

#--------------------------------------------------
#CREATING plots to analyze survival

#creating a column in analysis titled survival_time
# store either info from the days_to_death column or days_to_last_followup columns in survival_time
#For some patients, we can find survival time using their days_to_death info. 
#But for others, namely those who never had a registered death event (ie days_to_death = NA), 
#our best estimate becomes days_to_last_follow_up, since that is the longest definitive time we can say they survived.
analysis$survival_time <- ifelse(is.na(analysis$days_to_death), analysis$days_to_last_followup, analysis$days_to_death)

#create the death_event column
# based on the vital_status column, store a boolean variable that is TRUE if there was a death event,
#or FALSE if there was no recorded event
analysis$death_event <- ifelse(!is.na(analysis$vital_status),TRUE,FALSE)

#line below will save the survival plot to the outputs folder
jpeg("menopause_survival_plot.jpg")
# Creating the KM Plot for the menopause status
surv_object_menopause <- Surv(time = analysis$survival_time,
                              event = analysis$death_event)

# Create a fit object
menopause_fit <- surv_fit( surv_object_menopause ~ analysis$menopause_status,
                           data = analysis )

# the ggtheme and legend arguments are for formatting. 
survplot_menopause = ggsurvplot(menopause_fit, 
                                pval=TRUE, 
                                ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                                legend = "right")

# create the plot and label it accordingly 
KM_plot_menopause = survplot_menopause$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
KM_plot_menopause
dev.off()

#--------------------------------------------------
#CREATING MAF co-oncoplots to select for genes

#query in the MAF files 
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query)

maf <- GDCprepare(maf_query)

#change the bcr_patient_barcode column name, that way the MAF package can read our clinical file
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
write.csv(clinical, "analysis_clinical_data.csv", row.names = FALSE)

#query-in continued
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

#subsetting maf_object dataframe into two separate data frames( (pre-menopause patients and post-menopause patients)
#creating a boolean mask to differentiate Pre and Post
subset_menopause_mask <- ifelse(maf_object@clinical.data$menopause_status=="Pre", TRUE, FALSE)

#store the patient barcodes of pre-menopause patients in a vector
pre_menopause_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[subset_menopause_mask]
#using the vector in the tsb argument of subsetMaf to create the first data frame
pre_maf <- subsetMaf(maf = maf_object,
                       tsb = pre_menopause_barcodes)

#repeating the process but for post-menopause patients
post_menopause_barcodes <- maf_object@clinical.data$Tumor_Sample_Barcode[!subset_menopause_mask]
post_maf <- subsetMaf(maf = maf_object,
                     tsb = post_menopause_barcodes)

#line below will save the coOncoplot to the outputs folder
jpeg("pre_post_coOncoplot.jpg")
#plotting the pre-menopause and post-menopause patients side by side
coOncoplot(m1 = pre_maf, 
           m2 = post_maf, 
           m1Name = "pre", 
           m2Name = "post")
dev.off()

#notice that the GATA3 gene is twice as mutated in the pre-menopause group vs the post-menopause group
#we will further analyze the GATA3 gene

#--------------------------------------------------
###UTILIZING THE clinical.drug DATAFRAME

#since we are focusing on GATA3 mutations, we will select for Hormone Therapy patients
#Hormone Therapy is typically the therapy type prescribed to GATA3 mutated patients as it is a marker for ER positive BRCA
#to draw conclusions about the efficacy of drugs under Hormone Therapy, we will create a survival plot

#masking the clinical.drug data frame for only Hormone Therapy
hormone_therapy_mask <- ifelse(clinical.drug$therapy_types== "Hormone Therapy", T, F)
clinical.drug_therapy_masked <- clinical.drug[hormone_therapy_mask,]

#many of the rows in clinical.drug are duplicates, we will mask so there is only one of each kind
duplicated(clinical.drug_therapy_masked$bcr_patient_barcode)
clinical.drug_analysis<- clinical.drug_therapy_masked[duplicated(clinical.drug_therapy_masked$bcr_patient_barcode),]
#double check if there are leftover duplicates and mask again as necessaryy
duplicated(clinical.drug_analysis$bcr_patient_barcode)
clinical.drug_final<- clinical.drug_analysis[duplicated(clinical.drug_analysis$bcr_patient_barcode),]

#some of the drug names are written differently (ex. tamoxifen vs Tamoxifen, vs TAMOXIFEN)
unique(clinical.drug_final$drug_name)
#since the survival plot will be case sensitive (and plot a new line for each type of tamoxifen spelling) 
#we need to rename them as the same and save the new names to a column in the dataframe
clinical.drug_final$drug_name_corrected <- ifelse(clinical.drug_final$drug_name == "TAMOXIFEN" |clinical.drug_final$drug_name == "Tamoxifen", "Tamoxifen",
                      ifelse(clinical.drug_final$drug_name == "Letrozol" | clinical.drug_final$drug_name == "Letrozole" | clinical.drug_final$drug_name == "LETROZOLE", "Letrozole",
                             ifelse(clinical.drug_analysis$drug_name == "Arimidex"|clinical.drug_analysis$drug_name == "Arimidex (Anastrozole)"| clinical.drug_analysis$drug_name == "arimidex"|clinical.drug_analysis$drug_name == "ARIMIDEX", "Arimidex",
                                    ifelse(clinical.drug_final$drug_name == "Anastrazole"|clinical.drug_final$drug_name == "Anastrozole"|clinical.drug_final$drug_name == "ANASTROZOLE","Anastrozole",
                                           ifelse(clinical.drug_final$drug_name == "Aromasin"|clinical.drug_final$drug_name == "aromasin","Aromasin",
                                                  ifelse(clinical.drug_final$drug_name=="Femara"|clinical.drug_final$drug_name == "femara","Femara",
                                                         ifelse(clinical.drug_final$drug_name=="EXEMESTANE","Exemestane",
                                                                ifelse(clinical.drug_final$drug_name=="Prednisone","Prednisone","Other"))))))))

# to plot the survival plot, we need to make sure that the dimensions of the drugs (strata) and the dimensions of the time/event match
# thus we need to mask the analysis data frame we created before 
# this way, the patients that were selected for Hormone Therapy are also the same ones we plot for in the survival plot
drug_mask <- ifelse(clinical.drug_analysis$bcr_patient_barcode %in% analysis$patient_barcode , T, F) 
analysis_drug <- analysis[drug_mask,]
#the original analysis dataframe also had duplicates so we will mask them out as well
analysis_drug<- analysis_drug[duplicated(analysis_drug),]

#the line below will save the new survival plot
jpeg("drug_survival_plot.jpg")
# Creating the KM Plot for the type of drug
surv_object_drug <- Surv(time = analysis_drug$survival_time,
                              event = analysis_drug$death_event)

# Create a fit object
drug_fit <- surv_fit( surv_object_drug ~ clinical.drug_final$drug_name_corrected,
                           data = clinical.drug_final )

# the ggtheme and legend arguments are for formatting. 
survplot_drug = ggsurvplot(drug_fit, 
                                pval=TRUE, 
                                ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                                legend = "right")

# create the plot and label it accordingly (drug)
KM_plot_drug = survplot_drug$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
KM_plot_drug
dev.off()
