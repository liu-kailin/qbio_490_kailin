
"---
title: Partner Project
author: Kailin Liu
---"

# set working directory to analysis_data folder
#knitr::opts_knit$set(root.dir = normalizePath("/Users/kathykliu/Desktop/qbio_490_kailin/analysis_data"))

setwd("/Users/kathykliu/Desktop/qbio_490_kailin/analysis_data")

#read in brca_clinical_data.csv file from your local machine
clinical <- read.csv("/Users/kathykliu/Desktop/qbio_490_kailin/analysis_data/brca_clinical_data.csv")

# load in packages
library(TCGAbiolinks)
library(BiocManager)

# query clinical.drug and clinical.rad data frames
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

#--------------------------------------------------

#check for NAs in the lymph_node_examined_count column of clinical
lymph_na <- is.na(clinical$lymph_node_examined_count)
#determine the number of NAs in lymph_node_examined_count
sum(lymph_na)

#check for NAs in the stage_event_pathologic_stage column of clinical
stage_na <- is.na(clinical$stage_event_pathologic_stage)
#determine the number of NAs in stage_event_pathologic_stage
sum(stage_na)

#--------------------------------------------------

#plot the relationship between lymph_node_examined_count (lymph node count) and stage_event_pathologic_stage (stage of cancer)
lymph_stage_boxplot <- boxplot(formula=clinical$lymph_node_examined_count ~ clinical$stage_event_pathologic_stage,
                data=clinical,
                xlab ="Stage",
                ylab = "Lymph Node Count",
                main= "Lymph Node Count vs Stage of Breast Cancer",
                cex.axis=0.5)

#--------------------------------------------------
#CREATING  Kaplan-Meier plots to analyze survival

# load plots necessary to create the KM plots
library(ggplot2)

if (!require(survival)){ 
  install.packages("survival")
}
library(survival)

if (!require(survminer)){ 
  install.packages("survminer")
}
library(survminer)

#--------------------------------------------------

#creating a column in clinical titled survival_time
# store either info from the days_to_death column or days_to_last_followup columns in survival_time
#For some patients, we can find survival time using their days_to_death info. 
#But for others, namely those who never had a registered death event (ie days_to_death = NA), 
#our best estimate becomes days_to_last_follow_up, since that is the longest definitive time we can say they survived.
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)

#create the death_event column
# basted on the vital_status column, store a boolean variable that is TRUE if there was a death event,
#or FALSE if there was no recorded event
clinical$death_event <- ifelse(!is.na(clinical$vital_status),TRUE,FALSE)

#--------------------------------------------------
# Creating the KM Plot for the stage of cancer
surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
age_fit <- surv_fit( surv_object_age ~ clinical$stage_event_pathologic_stage,
                     data = clinical )

# the ggtheme and legend arguments are for formatting. 
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# create the plot and label it accordingly (stage)
KM_plot_stage = survplot_age$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_stage

#--------------------------------------------------
# Creating the KM Plot for the lymph node count

# Since there are too many strata for lymph node count, we will limit it to five strata
# Create an ifelse statement to categorize the data in the clinical$lymph_node_examined_count column into five categories
clinical$lymph_category <- ifelse(clinical$lymph_node_examined_count <= 10,  "minimal", 
                                  ifelse(clinical$lymph_node_examined_count <= 20, "low", 
                                         ifelse(clinical$lymph_node_examined_count <= 30, "medium",
                                                ifelse(clinical$lymph_node_examined_count <= 40, "elevated",
                                                     "high"))))

surv_object_age <- Surv(time = clinical$survival_time,
                        event = clinical$death_event)

# Create a fit object
age_fit <- surv_fit( surv_object_age ~ clinical$lymph_category,
                     data = clinical )

# the ggtheme and legend arguments are for formatting. 
survplot_age = ggsurvplot(age_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

# when you create plots on your own, be sure to name them descriptively
KM_plot_lymph = survplot_age$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_lymph

#--------------------------------------------------
# SAVING WORK

#saving edited clinical data frame
write.csv(clinical, "/Users/kathykliu/Desktop/qbio_490_kailin/week5_kailin/partner_clinical_data.csv", row.names = FALSE)

# saving the box plot
pdf("//Users/kathykliu/Desktop/qbio_490_kailin/week5_kailin/lymph_stage_boxplot.pdf")
lymph_stage_boxplot <- boxplot(formula=clinical$lymph_node_examined_count ~ clinical$stage_event_pathologic_stage,
                               data=clinical,
                               xlab ="Stage",
                               ylab = "Lymph Node Count",
                               main= "Lymph Node Count vs Stage of Breast Cancer",
                               cex.axis=0.5)
dev.off()

#saving the stage KM plot
pdf("/Users/kathykliu/Desktop/qbio_490_kailin/week5_kailin/stage_survival.pdf")
KM_plot_stage
dev.off()

#saving the lymph node count KM plot
pdf("/Users/kathykliu/Desktop/qbio_490_kailin/week5_kailin/lymph_survival.pdf")
KM_plot_lymph
dev.off()

