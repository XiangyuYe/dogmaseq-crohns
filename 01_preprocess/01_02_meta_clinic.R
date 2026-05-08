library(bigreadr)
library(dplyr)
library(stringr)
library(tibble)
library(glue)
##
project_path <- "/ix1/wchen/xiangyu/Projects/03_CD_DOGMA_test/"
setwd(project_path)

##
clinic_df <- read.csv("03_output/02_clean/Final_med_histological_for_submission_with_edits.csv")
# LC     RC     RM     SC     TC     TI 
# Left Colon           Rectum      Right Colon    Sigmoid Colon 
# Terminal Ileum Transverse Colon 
clinic_df$Section_specific <- recode(clinic_df$Anatomic.Location..Specific.,
                                     "Left Colon" = "LC",
                                     "Right Colon" = "RC",
                                     "Rectum" = "RM",
                                     "Sigmoid Colon" = "SC",
                                     "Transverse Colon" = "TC",
                                     "Terminal Ileum" = "TI",)
clinic_df$Section <- ifelse(clinic_df$Section_specific == "TI", "TI", "Colon")
clinic_df$Condition <- recode(clinic_df$Inflammation.State.of.Mucosa,
                              "Inflamed " = "I_I",
                              "Uninflamed in Uninflamed Region" = "N_N",
                              "Uninflamed Adjacent to Ulcers/Erosions" = "N_U")
clinic_df$Sample_exp <- paste(clinic_df$Sample.ID, 
                              clinic_df$Date.Used, 
                              sep = "_")
clinic_df$Sample_ID_exp <- paste(clinic_df$Sample.ID,
                                 "CD",
                                 clinic_df$Condition, 
                                 clinic_df$Section_specific, 
                                 clinic_df$Date.Used, 
                                 sep = "_")
clinic_df$Gender <- as.factor(clinic_df$Gender) %>% as.integer()
clinic_df$Race_specific <- clinic_df$Race
clinic_df$Race <- ifelse(clinic_df$Race_specific == "White", 1, 0)
##
med_list <- c("Corticosteroids", "Aminosalicylates", "Antibiotics", 
              "Immunomodulatory.Drugs", "Methotrexate", "Anti.TNF", 
              "Other.Biologics", "Other.Medications")
clinic_df$med_strategy <- lapply(med_list, function(x){
  ifelse(clinic_df[,x] == "n/a", NA, x)
  }) %>% 
  Reduce("cbind", .) %>%
  apply(., 1, function(xx){
    paste(xx[!is.na(xx)], collapse = "_")
    })
clinic_df$med_strategy[clinic_df$med_strategy == ""] <- "Unrecorded"
tab_med <- table(clinic_df$med_strategy)
clinic_df$med_strategy_comb <- ifelse(clinic_df$med_strategy %in% 
                                        names(tab_med)[tab_med < 10],
                                      "Other_strategies",
                                      clinic_df$med_strategy)
clinic_df$Anti_Inflam <- grepl("Corticosteroids|Aminosalicylates", clinic_df$med_strategy) %>% as.integer()
clinic_df$Immunomodulatory <- grepl("Immunomodulatory|Methotrexate", clinic_df$med_strategy) %>% as.integer()
clinic_df$Biologics <- grepl("Biologics|Anti.TNF", clinic_df$med_strategy) %>% as.integer()
clinic_df$use_Aminosalic <- ifelse(clinic_df$Aminosalicylates == "n/a", 0, 1)
clinic_df$use_Immunomodu <- ifelse(clinic_df$Immunomodulatory.Drugs == "n/a", 0, 1)
clinic_df$use_MTX <- ifelse(clinic_df$Methotrexate == "n/a", 0, 1)
clinic_df$use_AntiTNF <- ifelse(clinic_df$Anti.TNF == "n/a", 0, 1)
clinic_df$use_AntiIL23 <- ifelse(clinic_df$Other.Biologics %in% c("Risankizumab",
                                                                  "Risankizumab-rzaa", 
                                                                  "Ustekinumab"), 
                                 1, 0)
clinic_df$use_AntiIntegrin <- ifelse(clinic_df$Other.Biologics %in% c("Vedolizumab"), 
                                     1, 0)
saveRDS(clinic_df, file = "03_output/02_clean/clinic_df.rds")

