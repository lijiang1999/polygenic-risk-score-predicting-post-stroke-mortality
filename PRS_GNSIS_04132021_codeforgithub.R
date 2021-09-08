#---- This R code was created by Jiang Li (lijiang1999@gmail.com) for the following manuscript entitled 
#"Predicting Mortality among ischmeic stroke patients using pathways-derived polygenic risk scores"

library(readxl)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(broom)
library(splitstackshape)
library(maxstat)
library(glmnet)
library(survC1)
library(forestplot)
library(corrplot)


#####################   Section Title ######################

# 1. command line using PLINK and PRSice-2
# 2. Selection bias check
# 3. input candidate pathways
# 4. Cox regression (univariate)
# 5. calculate Event and Time from EHR (Durgesh)
# 6. LASSO glmnet.CV for modeling (multivariate)
# 7. Comparing models
# 8. Feature QC step
# 9. Subgroup analysis

############################ end ################################

# command line using PLINK and PRSice-2 -----------------------------------


##Command line script using PLINK--------------------------

# ./plink2 --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_merged \
# --keep Phenotype_newselected_remove20_discovery.txt \
# --make-bed \
# --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery \
# --threads 8 \
# --allow-no-sex

#21784 Subjects

# ./plink2 /
# --allow-no-sex /
# --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery /
# --indep-pairwise 50 5 0.2 /
# --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2 /
# --threads 4

# ./plink2 /
# --allow-no-sex /
# --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery /
# --extract GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2.prune.in /
# --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2_pcaapp /
# --pca approx /
# --threads 4

# ./plink2 /
#  --allow-no-sex /
#  --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery /
#  --extract GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2.prune.in /
#  --make-bed /
#  --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2 /
#  --threads 4
# 
# ./plink2 /
#  --allow-no-sex /
#  --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2 /
#  --genome /
#  --min 0.2 /
#  --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_prune50_5_2_genome /
#  --threads 4 

# ./plink2 /
# --bfile GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery /
# --extract extract_AIS_EUR.txt /
# --remove removecases_2nd_discovery.txt /
# --make-bed /
# --out GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_AIS /
# --pheno Phenotype_ALL_STROKE_update_69_PCA_distinct_PRS_2nd_10172020.txt /
# --threads 4 /
# --allow-no-sex


##Command line script using PRSice-2--------------------------
#to create pathway-specific polygenic risk scores

# Rscript PRSice.R --dir . \
# --prsice ./PRSice_linux \
# --bar-levels 0.001,0.005,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,1 \
# --base AIS_EUR_BETA_FINAL1_NEW.txt \
# --beta \
# --binary-target T \
# --clump-kb 250 \
# --clump-p 1.000000 \
# --clump-r2 0.080000 \
# --fastscore  \
# --cov covariate_ALL_STROKE_update_69_PCA_distinct_PRS_2nd_10172020.txt \
# --out bp_withoutgroup_noage_AIS_BETA_r2_3_MAF025_newselected_remove20_discovery_AIS_maf2 \
# --seed 222462130 \
# --target GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_AIS \
# --thread 8

# 2045 are cases and 19739 are controls
#--clump-r2 0.080000 giving the best result

#this is the right one for maf

# Rscript PRSice.R --dir . \ 
# --prsice ./PRSice_linux \
# --bar-levels 0.001,0.005,0.01,0.025,0.05,0.1,0.2,0.3,0.4,0.5,1 \
# --base AS_EUR_FINAL1.txt \
# --binary-target T \
# --clump-kb 250 \
# --clump-p 1.000000 \
# --clump-r2 0.100000 \
# --cov covariate_ALL_STROKE_update_69_PCA_distinct_PRS_2nd_10172020.txt \
# --or \
# --out bp_withoutgroup_noage_AS_BETA_r2_3_MAF025_newselected_remove20_discovery_AS_notmaf \
# --seed 222462130 \
# --stat OR \
# --target GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_AS \
# --thread 4
# 
# Rscript PRSice.R --dir . \ 
# --prsice ./PRSice_linux \
# --bar-levels 0.025,1 \
# --base AIS_EUR_FINAL1.txt \
# --binary-target T \
# --clump-kb 250 \
# --clump-p 1.000000 \
# --clump-r2 0.100000 \
# --fastscore  \
# --feature exon,gene,protein_coding,CDS \
# --gtf Homo_sapiens.GRCh37.87.gtf \
# --msigdb ./msigdb_v7.0_GMTs/c5.bp.v7.0.symbols.gmt \
# --cov covariate_ALL_STROKE_update_69_PCA_distinct_PRS_2nd_10172020.txt \
# --or \
# --out bp_withoutgroup_noage_AIS_BETA_r2_3_MAF025_newselected_remove20_discovery_AIS_notmaf_genesets \
# --seed 222462130 \
# --stat OR \
# --target GHSf6030k_maf01_geno1_mind1_hwe0001_STROKE_newselected_remove20_discovery_AIS \
# --thread 4



# Selection bias check ----------------------------------------------------

##Supplementary Figure 1------------------------- 
##to show the selection bias for MyCode sample when comparing to nonMyCode sample

##Total cumulative incidence of all-cause 3-year or 5-year mortality after ischemic stroke, stratified by age at index (binary), sex, and MyCode identity.
colnames(GNSIS_final_final)
length(unique(GNSIS_final_final$PT_ID))

GNSIS_final_final_complete <- GNSIS_final_final %>%  #15822
  filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3))  #12883 for 3yr; #10582 for 5yr; #15202 for 1yr

GNSIS_final_final$cohort <- ifelse(GNSIS_final_final$PT_ID %in% PRS_cohort_characteristics_update_pathway$FID, "MyCode", "NONMyCode")

GNSIS_final_final_final <- GNSIS_final_final %>%  #15822
  filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3))  #12883 for 3yr; #10582 for 5yr; #15202 for 1yr
colnames(GNSIS_final_final_final)

## KM stratified by cohort
tiff("Result_km_cohort.tiff", units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~as.factor(cohort) + PT_SEX + AGE_AT_INDEX_BINARY, data = GNSIS_final_final_final),                 # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.5),
  pval.size = 4,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(11, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 3.1),         # present narrower X axis, but not affect
  ylim = c(0, 0.5),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("MyCode Female Young", "MyCode Female Old", "MyCode Male Young", "MyCode Male Old", "NONMyCode Female Young", "NONMyCode Female Old", "NONMyCode Male Young", "NONMyCode Male Old"),    # change legend labels.
  font.title    = c(11, "bold.italic", "darkgreen"),
  font.subtitle = c(10, "bold", "darkgreen"),
  font.caption  = c(10, "plain", "darkgreen"),
  font.x        = c(11, "bold.italic", "black"),
  font.y        = c(11, "bold.italic", "black"),
  font.xtickslab = c(11, "bold", "black"),
  font.ytickslab = c(11, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 3yrs mortality",
    subtitle = "MyCode vs Non MyCode patients",
    caption  = "Plotted with R survminer")
dev.off()

# Pairwise comparisons using Log-Rank test to show the result embedded in Supplementary Figure 1
out_pairwise <- pairwise_survdiff(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~cohort + PT_SEX + AGE_AT_INDEX_BINARY, data = GNSIS_final_final_final)                 # survfit object with calculated statistics.
out_pairwise

GNSIS_final_final_final <- GNSIS_final_final %>%  #15822
  filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*5))  #12883 for 3yr; #10582 for 5yr; #15202 for 1yr
colnames(GNSIS_final_final_final)

tiff("Result_km_cohort_5yr.tiff", units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_death_censor_5yr/365.25, event_death_censor_5yr)~as.factor(cohort) + PT_SEX + AGE_AT_INDEX_BINARY, data = GNSIS_final_final_final),                 # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.5),
  pval.size = 4,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(11, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 5.1),         # present narrower X axis, but not affect
  ylim = c(0, 0.65),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("MyCode Female Young", "MyCode Female Old", "MyCode Male Young", "MyCode Male Old", "NONMyCode Female Young", "NONMyCode Female Old", "NONMyCode Male Young", "NONMyCode Male Old"),    # change legend labels.
  font.title    = c(11, "bold.italic", "darkgreen"),
  font.subtitle = c(10, "bold", "darkgreen"),
  font.caption  = c(10, "plain", "darkgreen"),
  font.x        = c(11, "bold.italic", "black"),
  font.y        = c(11, "bold.italic", "black"),
  font.xtickslab = c(11, "bold", "black"),
  font.ytickslab = c(11, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 5yrs mortality",
    subtitle = "MyCode vs Non MyCode patients",
    caption  = "Plotted with R survminer")
dev.off()

out_pairwise <- pairwise_survdiff(Surv(timediff_death_censor_5yr/365.25, event_death_censor_5yr)~cohort + PT_SEX + AGE_AT_INDEX_BINARY, data = GNSIS_final_final_final)                 # survfit object with calculated statistics.
out_pairwise




# input candidate pathways ------------------------------------------------

##identifying candidate pathways-specific PRSs-----------------------
all_pathway <- read.table("Candidate_all_pathway.txt", header = T, stringsAsFactors = F)
head(all_pathway)
list_df_pathway <- merge(Result_GNSIS13_HR_AIS_3yr_PRS_candidate, all_pathway, by = "Category", all.x = T)

list_df_pathway <- read.csv("Result_GNSIS13_HR_AIS_3yr_PRS_candidate_quantitative_70percent_univariate_adjustcovariate.csv", header = T, stringsAsFactors = F)
unique(list_df_pathway$Category)
list_df_pathway <- list_df_pathway %>%
  filter(PHENOTYPE == "3yr_death") %>%
  rename(Set = Header)
head(list_df_pathway)
candidate_pathway <- read.csv("candidate_pathway_maf1_maf25.csv", header = T, stringsAsFactors = F)
head(candidate_pathway)
list_df_pathway_direction <- merge(list_df_pathway, candidate_pathway, by = "Set")
#make sure disease risk and outcome risk in the same direction
list_df_pathway_direction$direction <- (log(list_df_pathway_direction$estimate)/abs(log(list_df_pathway_direction$estimate)))*(list_df_pathway_direction$Effect/abs(list_df_pathway_direction$Effect))
list_df_pathway_direction_diseaserisk <- list_df_pathway_direction[list_df_pathway_direction$direction == 1 & list_df_pathway_direction$p.value <= 0.10, ] #select candidates
list_df_pathway_direction_diseaserisk <- list_df_pathway_direction[list_df_pathway_direction$direction == 1, ]

unique(list_df_pathway_direction_diseaserisk$Category)


# Cox regression (univariate) ----------------------------------------------------------


##This function has been improved to do cox regression on binary variables with 75 percentile fixed cutpoint --- 
### this function not being used in this study
#short format for 3yr, latest version
PRS_TEST <- data.frame()
pathway <- NULL
PRS_GNSIS <- NULL
cox_3yrdeath_prs1 <- NULL
cox_3yrrecur_prs1 <- NULL
coxtidy_3yrdeath_prs1 <- NULL
coxtidy_3yrrecur_prs1 <- NULL
DEATH <- NULL
RECUR <- NULL
dat_DEATH <- data.frame()
BEST_TEST <- data.frame()
dat_RECUR <- NULL
dat <- NULL
PRS_binary <- NULL
PRS_GNSIS_binary <- NULL

Result_HR = function(PRS_file, FAM_file, COV_file, GNSIS_file) {       # for COV_split file
  pathway <- colnames(PRS_file[, 3:(ncol(PRS_file)-3)])    
  for (i in pathway){
    print(i)
    PRS_TEST <- PRS_file %>%
      select(FID, i) %>%
      rename(PRS = i)
    BEST_TEST <- merge(FAM_file[, c("FID", "GROUP")], PRS_TEST[, c("FID", "PRS")], by = "FID")
    BEST_TEST <- merge(BEST_TEST, COV_file, by = "FID")
    #BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP.x==2, ]  #for discovery dataset
    BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP==2, ] #for replication dataset
    #BEST_TEST$PRS_norm <- (BEST_TEST$PRS-min(BEST_TEST$PRS))/(max(BEST_TEST$PRS)-min(BEST_TEST$PRS)) #different normalization method
    colnames(BEST_TEST_CASE)[colnames(BEST_TEST_CASE) == 'FID'] <- 'PT_ID' 
    PRS_GNSIS <- inner_join(BEST_TEST_CASE, GNSIS_file, by = "PT_ID")
    
    #select PT_ID with complete 3_yr observation
    PRS_GNSIS <- PRS_GNSIS %>%
      filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3)) 
    
    #percentile
    PRS_GNSIS$PRS_norm <- perc.rank(PRS_GNSIS$PRS)
    
    # Creating variables PRS_*_cat to denote PRS > 75th percentile and PRS <= 25th percentile from PRS_*. PRS values between these turned to NA and excluded from KM and COX analysis.
    PRS_GNSIS <- PRS_GNSIS %>% 
      mutate(PRS_norm_cat = as.factor(
        case_when(PRS_norm > quantile(PRS_norm, 0.75) ~ "PRS_norm_above75",
                  #PRS_norm <= quantile(PRS_norm, 0.25) ~ "PRS_norm_below25")))
                  PRS_norm <= quantile(PRS_norm, 0.75) ~ "PRS_norm_below75")))
    # Making PRS_*_below25 as reference level for the columns PRS_*_cat
    #PRS_GNSIS$PRS_norm_cat <- relevel(PRS_GNSIS$PRS_norm_cat, ref = "PRS_norm_below25")
    PRS_GNSIS$PRS_norm_cat <- relevel(PRS_GNSIS$PRS_norm_cat, ref = "PRS_norm_below75")
    PRS_GNSIS$Category <- i
    PRS_GNSIS_binary[[i]] <- data.frame(PRS_GNSIS[, c("PT_ID", "PRS_norm_cat", "Category")])
    
    #cox regression model
    tryCatch({
      #print(PRS_GNSIS) #for QC
      # 3yr Motality stroke ~ PRS_norm_cat
      cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX) + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
      coxtidy_3yrdeath_prs1 <- tidy(cox_3yrdeath_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrdeath_prs1$Category <- i
      DEATH[[i]] <- coxtidy_3yrdeath_prs1
      
      #cox.zph(cox_3yrdeath_prs1) # proportional hazards assumptions
      # 3yr Recur stroke ~ PRS_norm_cat
      cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX) + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
      coxtidy_3yrrecur_prs1 <- tidy(cox_3yrrecur_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrrecur_prs1$Category <- i
      RECUR[[i]] <- coxtidy_3yrrecur_prs1
      #cox.zph(cox_3yrrecur_prs1) # proportional hazards assumptions
    }, error=function(e){})
  }
  dat_DEATH <- data.frame(do.call(rbind,DEATH))
  print(dat_DEATH)
  dat_RECUR <- data.frame(do.call(rbind,RECUR))
  #print(dat_RECUR)
  colnames(dat_DEATH) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  colnames(dat_RECUR) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  dat_DEATH$PHENOTYPE <- "3yr_death"
  dat_RECUR$PHENOTYPE <- "3yr_recur"
  dat <- rbind(dat_DEATH, dat_RECUR)
  PRS_binary <- data.frame(do.call(rbind,PRS_GNSIS_binary))
  save(PRS_binary, file="PRS_binary_train_test_70percent_univariate_3yr_adjustcovariate.RData", version = 2)
  return(dat)
}

Result_GNSIS13_HR_AIS_3yr <- Result_HR(BEST_split, FAM_split, COV_split, GNSIS_final_final)  #without filter 3 years.
Result_GNSIS13_HR_AIS_3yr
Result_GNSIS13_HR_AIS_3yr_PRS <- Result_GNSIS13_HR_AIS_3yr[Result_GNSIS13_HR_AIS_3yr$term == "PRS_norm_catPRS_norm_above75", ]


##This function has been improved to do cox regression on binary variables with maxstat cutpoint converted files ---
###short format for 3yr, on training dataset
PRS_TEST <- data.frame()
pathway <- NULL
PRS_GNSIS <- NULL
cox_1yrdeath_prs1 <- NULL
cox_1yrrecur_prs1 <- NULL
coxtidy_1yrdeath_prs1 <- NULL
coxtidy_1yrrecur_prs1 <- NULL
DEATH <- NULL
RECUR <- NULL
dat_DEATH <- data.frame()
BEST_TEST <- data.frame()
dat_RECUR <- NULL
dat <- NULL
PRS_binary <- NULL
PRS_GNSIS_binary <- NULL

Result_HR = function(PRS_file, FAM_file, COV_file, GNSIS_file) {       # for COV_split file
  pathway <- colnames(PRS_file[, 2:(ncol(PRS_file)-1)])    
  #pathway <- c("V6858", "V4451", "V3259", "V3046")
  for (i in pathway){
    print(i)
    PRS_TEST <- PRS_file %>%
      select(FID, i) %>%
      rename(PRS_norm_cat = i) #change into PRS_norm_cat
    BEST_TEST <- merge(FAM_file[, c("FID", "GROUP")], PRS_TEST[, c("FID", "PRS_norm_cat")], by = "FID")
    BEST_TEST <- merge(BEST_TEST, COV_file, by = "FID")
    BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP==2, ] #for replication dataset
    #BEST_TEST$PRS_norm <- (BEST_TEST$PRS-min(BEST_TEST$PRS))/(max(BEST_TEST$PRS)-min(BEST_TEST$PRS)) #different normalization method
    colnames(BEST_TEST_CASE)[colnames(BEST_TEST_CASE) == 'FID'] <- 'PT_ID' 
    #print(dim(BEST_TEST_CASE))
    #PRS_GNSIS <- inner_join(BEST_TEST_CASE[, c("PT_ID", "PRS_norm")], GNSIS_file, by = "PT_ID")
    PRS_GNSIS <- inner_join(BEST_TEST_CASE, GNSIS_file, by = "PT_ID")
    
    #select PT_ID with complete 3_yr observation
    PRS_GNSIS <- PRS_GNSIS %>%
      filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*1)) 
    #print(PRS_GNSIS)
    PRS_GNSIS$Category <- i
    PRS_GNSIS_binary[[i]] <- data.frame(PRS_GNSIS[, c("PT_ID", "PRS_norm_cat", "Category")])
    
    #cox regression model
    tryCatch({
      #print(PRS_GNSIS) #for QC
      # 1yr Motality stroke ~ PRS_norm_cat
      #cox_1yrdeath_prs1 <- coxph(Surv(timediff_death_censor_1yr/365.25, event_death_censor_1yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX_BINARY + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_1yrdeath_prs1 <- coxph(Surv(timediff_death_censor_1yr/365.25, event_death_censor_1yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX), data = PRS_GNSIS)
      cox_1yrdeath_prs1 <- coxph(Surv(timediff_death_censor_1yr/365.25, event_death_censor_1yr)~PRS_norm_cat, data = PRS_GNSIS)
      coxtidy_1yrdeath_prs1 <- tidy(cox_1yrdeath_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_1yrdeath_prs1$Category <- i
      DEATH[[i]] <- coxtidy_1yrdeath_prs1
      
      #cox.zph(cox_1yrdeath_prs1) # proportional hazards assumptions
      # 1yr Recur stroke ~ PRS_norm_cat
      #cox_1yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_1yr/365.25, event_recur_censor_1yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX_BINARY + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_1yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_1yr/365.25, event_recur_censor_1yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX), data = PRS_GNSIS)
      cox_1yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_1yr/365.25, event_recur_censor_1yr)~PRS_norm_cat, data = PRS_GNSIS)
      coxtidy_1yrrecur_prs1 <- tidy(cox_1yrrecur_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_1yrrecur_prs1$Category <- i
      RECUR[[i]] <- coxtidy_1yrrecur_prs1
      #cox.zph(cox_1yrrecur_prs1) # proportional hazards assumptions
    }, error=function(e){})
  }
  dat_DEATH <- data.frame(do.call(rbind,DEATH))
  print(dat_DEATH)
  dat_RECUR <- data.frame(do.call(rbind,RECUR))
  #print(dat_RECUR)
  colnames(dat_DEATH) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  colnames(dat_RECUR) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  dat_DEATH$PHENOTYPE <- "1yr_death"
  dat_RECUR$PHENOTYPE <- "1yr_recur"
  dat <- rbind(dat_DEATH, dat_RECUR)
  save(dat, file="PRS_binary_train_70percent_univariate_1yr.RData", version = 2)
  return(dat)
}

####for all samples (train + test)
colnames(PRS_binary_wideformat)
colnames(PRS_binary_wideformat)[which(names(PRS_binary_wideformat) == "PT_ID")] <- "FID"
Result_GNSIS13_HR_AIS_3yr_train <- Result_HR(PRS_binary_wideformat, FAM_split, COV_split, GNSIS_final_final)

####70percent split
PRS_binary_wideformat_train <- PRS_binary_wideformat[PRS_binary_wideformat$FID %in% BEST_split_train_70percent$FID, ] #1226
colnames(PRS_binary_wideformat_train)
Result_GNSIS13_HR_AIS_3yr_train <- Result_HR(PRS_binary_wideformat_train, FAM_split, COV_split, GNSIS_final_final)

PRS_binary_wideformat_test <- PRS_binary_wideformat[PRS_binary_wideformat$FID %in% BEST_split_test_70percent$FID, ] #530
colnames(PRS_binary_wideformat_test)

##This function has been improved to do cox regression on quantitative variables ---------------------------
###short format for 3yr, on training dataset or testing dataset or both
library(dplyr)
PRS_TEST <- data.frame()
pathway <- NULL
PRS_GNSIS <- NULL
cox_3yrdeath_prs1 <- NULL
cox_3yrrecur_prs1 <- NULL
coxtidy_3yrdeath_prs1 <- NULL
coxtidy_3yrrecur_prs1 <- NULL
DEATH <- NULL
RECUR <- NULL
dat_DEATH <- data.frame()
BEST_TEST <- data.frame()
dat_RECUR <- NULL
dat <- NULL
PRS_binary <- NULL
PRS_GNSIS_binary <- NULL

Result_HR = function(PRS_file, FAM_file, COV_file, GNSIS_file) {       # for COV_split file
  pathway <- colnames(PRS_file[, 3:(ncol(PRS_file)-3)])    
  for (i in pathway){
    print(i)
    PRS_TEST <- PRS_file %>%
      dplyr::select(FID, i) %>%
      dplyr::rename(PRS = i) #change into PRS
    BEST_TEST <- merge(FAM_file[, c("FID", "GROUP")], PRS_TEST[, c("FID", "PRS")], by = "FID")
    BEST_TEST <- merge(BEST_TEST, COV_file, by = "FID")
    BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP==2, ] #for replication dataset
    #BEST_TEST$PRS_norm <- (BEST_TEST$PRS-min(BEST_TEST$PRS))/(max(BEST_TEST$PRS)-min(BEST_TEST$PRS)) #different normalization method
    colnames(BEST_TEST_CASE)[colnames(BEST_TEST_CASE) == 'FID'] <- 'PT_ID' 
    print(dim(BEST_TEST_CASE))
    #PRS_GNSIS <- inner_join(BEST_TEST_CASE[, c("PT_ID", "PRS_norm")], GNSIS_file, by = "PT_ID")
    PRS_GNSIS <- inner_join(BEST_TEST_CASE, GNSIS_file, by = "PT_ID")
    
    #select PT_ID with complete 3_yr observation
    PRS_GNSIS <- PRS_GNSIS %>%
      filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3)) 
    print(dim(PRS_GNSIS))
    
    #cox regression model
    tryCatch({
      #print(PRS_GNSIS) #for QC
      # 3yr Motality stroke ~ PRS
      cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~scale(PRS) + as.factor(PT_SEX.x) + AGE_AT_INDEX_BINARY + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~scale(PRS) + as.factor(PT_SEX.x) + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      coxtidy_3yrdeath_prs1 <- tidy(cox_3yrdeath_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrdeath_prs1$Category <- i
      DEATH[[i]] <- coxtidy_3yrdeath_prs1
      
      #cox.zph(cox_3yrdeath_prs1) # proportional hazards assumptions
      # 3yr Recur stroke ~ PRS
      cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~scale(PRS) + as.factor(PT_SEX.x) + AGE_AT_INDEX_BINARY + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~scale(PRS) + as.factor(PT_SEX.x) + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      coxtidy_3yrrecur_prs1 <- tidy(cox_3yrrecur_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrrecur_prs1$Category <- i
      RECUR[[i]] <- coxtidy_3yrrecur_prs1
      #cox.zph(cox_3yrrecur_prs1) # proportional hazards assumptions
    }, error=function(e){})
  }
  dat_DEATH <- data.frame(do.call(rbind,DEATH))
  print(dat_DEATH)
  dat_RECUR <- data.frame(do.call(rbind,RECUR))
  #print(dat_RECUR)
  colnames(dat_DEATH) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  colnames(dat_RECUR) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  dat_DEATH$PHENOTYPE <- "3yr_death"
  dat_RECUR$PHENOTYPE <- "3yr_recur"
  dat <- rbind(dat_DEATH, dat_RECUR)
  save(dat, file="PRS_quantitative_train_test_70percent_univariate_3yr_adjustcovariate.RData", version = 2)
  return(dat)
}

####for all sample (train + test)
Result_GNSIS13_HR_AIS_3yr_train_test <- Result_HR(BEST_split_train_test_70percent, FAM_split, COV_split, GNSIS_final_final)  #1756*87

####young stroke
GNSIS_final_final_young <- GNSIS_final_final[GNSIS_final_final$AGE_AT_INDEX_BINARY == 0, ]
Result_GNSIS13_HR_AIS_1yr_train_test <- Result_HR(BEST_split_train_test_70percent, FAM_split, COV_split, GNSIS_final_final_young)  #738*87 for 5yr; 877*87 for 3yrs; 909*87 for 1yrs
####old stroke
GNSIS_final_final_old <- GNSIS_final_final[GNSIS_final_final$AGE_AT_INDEX_BINARY == 1, ]
Result_GNSIS13_HR_AIS_5yr_train_test <- Result_HR(BEST_split_train_test_70percent, FAM_split, COV_split, GNSIS_final_final_old)  #642*87 for 5yr; 879*87 for 3yrs; 903*87 for 1yrs



# calculate Event and Time from EHR ----------------------------------------

##This section of codes was originally written by Durgesh Chaudhary

GNSIS <- GNSIS %>% 
  mutate(timediff_death_censor = 
           case_when(!is.na(PT_DEATH_DT) ~ as.numeric(difftime(PT_DEATH_DT, INDEX_DT, units = "days")),
                     is.na(PT_DEATH_DT) ~ as.numeric(difftime(LAST_ACTIVE_DT, INDEX_DT, units = "days"))),
         event_death_censor =
           case_when(PT_STATUS == "ALIVE" ~ 0,
                     PT_STATUS == "DECEASED" ~ 1))


##   ____________________________________________________________________________
##   Limiting to 1/3/5 years for mortality                                   

## Event: Limiting to Death within 5 year (event_death_censor_5yr), 3 year (event_death_censor_3yr), 1 year (event_death_censor_1yr)
## Time to event: timediff_death_censor_5yr, timediff_death_censor_3yr, timediff_death_censor_1yr
GNSIS <- GNSIS %>% 
  mutate(timediff_death_censor_5yr =
           case_when(timediff_death_censor < 5*365.25 ~ timediff_death_censor,
                     timediff_death_censor >= 5*365.25 ~ 5*365.25),
         event_death_censor_5yr =
           case_when(PT_STATUS == "ALIVE" ~ 0,
                     PT_STATUS == "DECEASED" & timediff_death_censor <= 5*365.25 ~ 1,
                     PT_STATUS == "DECEASED" & timediff_death_censor > 5*365.25 ~ 0),
         timediff_death_censor_3yr =
           case_when(timediff_death_censor < 3*365.25 ~ timediff_death_censor,
                     timediff_death_censor >= 3*365.25 ~ 5*365.25),
         event_death_censor_3yr =
           case_when(PT_STATUS == "ALIVE" ~ 0,
                     PT_STATUS == "DECEASED" & timediff_death_censor <= 3*365.25 ~ 1,
                     PT_STATUS == "DECEASED" & timediff_death_censor > 3*365.25 ~ 0),
         timediff_death_censor_1yr =
           case_when(timediff_death_censor < 1*365.25 ~ timediff_death_censor,
                     timediff_death_censor >= 1*365.25 ~ 5*365.25),
         event_death_censor_1yr =
           case_when(PT_STATUS == "ALIVE" ~ 0,
                     PT_STATUS == "DECEASED" & timediff_death_censor <= 1*365.25 ~ 1,
                     PT_STATUS == "DECEASED" & timediff_death_censor > 1*365.25 ~ 0))

GNSIS$timediff_death_censor_5yr[GNSIS$timediff_death_censor_5yr == 0] <- 0.5
GNSIS$timediff_death_censor_3yr[GNSIS$timediff_death_censor_3yr == 0] <- 0.5
GNSIS$timediff_death_censor_1yr[GNSIS$timediff_death_censor_1yr == 0] <- 0.5

## Event: RECUR_STROKE (already available in data)
## Time to event: timediff_recur_censor
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor = 
           case_when(!is.na(RECUR_STROKE_DT) ~ as.numeric(difftime(RECUR_STROKE_DT, INDEX_DT, units = "days")),
                     is.na(RECUR_STROKE_DT) ~ as.numeric(difftime(LAST_ACTIVE_DT, INDEX_DT, units = "days")))) 

##   ____________________________________________________________________________
##   Limiting to 1/3/5 years for recurrence                                 ####

## Event: Limiting to RECUR STROKE in 5 year (event_recur_censor_5yr), 3 year (event_recur_censor_3yr), 1 year (event_recur_censor_1yr)
## Time to event: timediff_recur_censor_5yr, timediff_recur_censor_3yr, timediff_recur_censor_1yr
GNSIS <- GNSIS %>% 
  mutate(timediff_recur_censor_5yr =
           case_when(timediff_recur_censor < 5*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 5*365.25 ~ 5*365.25),
         event_recur_censor_5yr =
           case_when(RECUR_STROKE == 0 ~ 0,
                     RECUR_STROKE == 1 & timediff_recur_censor <= 5*365.25 ~1,
                     RECUR_STROKE == 1 & timediff_recur_censor > 5*365.25 ~0),
         timediff_recur_censor_3yr =
           case_when(timediff_recur_censor < 3*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 3*365.25 ~ 5*365.25),
         event_recur_censor_3yr =
           case_when(RECUR_STROKE == 0 ~ 0,
                     RECUR_STROKE == 1 & timediff_recur_censor <= 3*365.25 ~1,
                     RECUR_STROKE == 1 & timediff_recur_censor > 3*365.25 ~0),
         timediff_recur_censor_1yr =
           case_when(timediff_recur_censor < 1*365.25 ~ timediff_recur_censor,
                     timediff_recur_censor >= 1*365.25 ~ 5*365.25),
         event_recur_censor_1yr =
           case_when(RECUR_STROKE == 0 ~ 0,
                     RECUR_STROKE == 1 & timediff_recur_censor <= 1*365.25 ~1,
                     RECUR_STROKE == 1 & timediff_recur_censor > 1*365.25 ~0))

GNSIS$timediff_recur_censor_5yr[GNSIS$timediff_recur_censor_5yr == 0] <- 0.5
GNSIS$timediff_recur_censor_3yr[GNSIS$timediff_recur_censor_3yr == 0] <- 0.5
GNSIS$timediff_recur_censor_1yr[GNSIS$timediff_recur_censor_1yr == 0] <- 0.5

# LASSO glmnet.CV for modeling --------------------------------------------------


##glmnet modeling including clinical risk factors-----------------------
load("PRS_cohort_characteristics_update_pathway.RData")
colnames(PRS_cohort_characteristics_update_pathway_train)
PRS_GNSIS_train_select_removeNA$PT_ID
PRS_cohort_characteristics_update_pathway_train <- PRS_cohort_characteristics_update_pathway[PRS_cohort_characteristics_update_pathway$FID %in% PRS_GNSIS_train_select_removeNA$PT_ID, ]

y_tr <- Surv(PRS_cohort_characteristics_update_pathway_train$timediff_death_censor_3yr/365.25, PRS_cohort_characteristics_update_pathway_train$event_death_censor_3yr)
results <- lapply(PRS_cohort_characteristics_update_pathway_train[, c(42:72)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48, 71, 51, 56, 66, 61, 44, 69, 52, 47, 59, 38, 53, 65, 67, 72, 49, 60)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48, 71, 51, 56, 66, 61, 44)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48)], function(x) scale(x, center = TRUE, scale = TRUE))

result <- as.data.frame(do.call(cbind, results)) 
colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_train[, c(42:72)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48, 71, 51, 56, 66, 61, 44, 69, 52, 47, 59, 38, 53, 65, 67, 72, 49, 60)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48, 71, 51, 56, 66, 61, 44)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_train[, c(55, 46, 48)])

x_tr <- as.matrix(result)
clinical_tr <- PRS_cohort_characteristics_update_pathway_train[, colnames(PRS_cohort_characteristics_update_pathway_train) %in% c("AGE_AT_INDEX_BINARY","hypertension", "diabetes", "dyslipidemia", "coronary_artery", "atrial_fib", "BMI_Overweight", "smoking")]
x_tr <- cbind(x_tr, clinical_tr)
x_tr[is.na(x_tr)] <- 0
colnames(x_tr) #smoking was  included because of missing value == 0
x_tr <- x_tr[, c(32:39)] #clinical features only

#save(y_tr, x_tr, file = "05172021_MAFALL_training_quantitative_univariate_clinicalgenetic.RData", version = 2)
#save(y_tr, x_tr, file = "05172021_MAFALL_training_quantitative_univariate_clinicalgenetic_p05.RData", version = 2)
#save(y_tr, x_tr, file = "05172021_MAFALL_training_quantitative_univariate_clinicalgenetic_p025.RData", version = 2)
save(y_tr, x_tr, file = "05172021_MAFALL_training_quantitative_univariate_clinicalgenetic_p01.RData", version = 2)

load("05172021_MAFALL_training_quantitative_univariate_clinicalgenetic.RData")

PRS_cohort_characteristics_update_pathway_test <- PRS_cohort_characteristics_update_pathway[PRS_cohort_characteristics_update_pathway$FID %in% PRS_GNSIS_test_select_removeNA$PT_ID, ]


y_te <- Surv(PRS_cohort_characteristics_update_pathway_test$timediff_death_censor_3yr/365.25, PRS_cohort_characteristics_update_pathway_test$event_death_censor_3yr)
results <- lapply(PRS_cohort_characteristics_update_pathway_test[, c(42:72)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48, 71, 51, 56, 66, 61, 44, 69, 52, 47, 59, 38, 53, 65, 67, 72, 49, 60)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48, 71, 51, 56, 66, 61, 44)], function(x) scale(x, center = TRUE, scale = TRUE))
#results <- lapply(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48)], function(x) scale(x, center = TRUE, scale = TRUE))

result <- as.data.frame(do.call(cbind, results)) 
colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_test[, c(42:72)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48, 71, 51, 56, 66, 61, 44, 69, 52, 47, 59, 38, 53, 65, 67, 72, 49, 60)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48, 71, 51, 56, 66, 61, 44)])
#colnames(result) <- colnames(PRS_cohort_characteristics_update_pathway_test[, c(55, 46, 48)])

x_te <- as.matrix(result)
clinical_te <- PRS_cohort_characteristics_update_pathway_test[, colnames(PRS_cohort_characteristics_update_pathway_test) %in% c("AGE_AT_INDEX_BINARY","hypertension", "diabetes", "dyslipidemia", "coronary_artery", "atrial_fib", "BMI_Overweight", "smoking")]
x_te <- cbind(x_te, clinical_te)
x_te[is.na(x_te)] <- 0
x_te
colnames(x_te) #smoking was included because of missing value == 0
x_te <- x_te[, c(32:39)] #clinical features only


library(glmnet)
x_tr <- as.matrix(x_tr)
x_te <- as.matrix(x_te)
colnames(x_tr)


cvfit_tr <- cv.glmnet(x_tr, y_tr, family="cox", maxit = 10000, alpha = 1, lambda = lambdas_to_try, nfolds = 5, type.measure = "C") #Only deviance, C available as type.measure for Cox models; deviance used instead 
cvfit_tr <- cv.glmnet(x_tr, y_tr, family="cox", maxit = 10000, alpha = 1, lambda = lambdas_to_try, nfolds = 5, set.seed = 1234)
fit_tr <- glmnet(x_tr, y_tr, family="cox", maxit = 10000, alpha = 1, lambda = lambdas_to_try, nfolds = 5, type.measure = "C")
fit_tr <- glmnet(x_tr, y_tr, family="cox", maxit = 10000, alpha = 1, lambda = lambdas_to_try, nfolds = 5)

plot(fit_tr, label = T)
plot(fit_tr, xvar="lambda", label = T)
dimnames(fit_tr$beta)[[1]]



##Figure 2-------------------------------
##LASSO process for feature selection

library(plotmo) # for plot_glmnet
tiff("Rplot_plotmo_39_w_age1_clinicalgenetic.tiff", units="in", width=8, height=8, res=300)
plot_glmnet(fit_tr, xvar="lambda", label = 10)                             # default colors and label the 5 biggest final coefs
dev.off()

tiff("Rplot_plot_39_w_age1_clinicalgenetic.tiff", units="in", width=8, height=8, res=300)
plot(cvfit_tr)
dev.off()

cvfit_tr$lambda.min
cvfit_tr$lambda.1se

coef(cvfit_tr, s = "lambda.min")
coef(cvfit_tr, s = "lambda.1se")
#With selected model, make linear predictors on test data, split into two or three groups
co_mat <- as.matrix(coef(cvfit_tr,s="lambda.min"))
co_mat <- as.matrix(coef(cvfit_tr,s="lambda.1se"))
co_mat
# rank(co_mat)
# co_mat[rank(co_mat),]
list <- co_mat[co_mat[,1]!=0,]
list
list_df <- data.frame(stack(unlist(list)))
list_df
colnames(list_df)[which(names(list_df) == "ind")] <- "Category"
all_pathway <- read.table("Candidate_all_pathway.txt", header = T, stringsAsFactors = F)
head(all_pathway)
list_df_pathway <- merge(list_df, all_pathway, by = "Category", all.x = T)
list_df_pathway <- list_df_pathway[order(list_df_pathway$values), ]
list_df_pathway$HR <- exp(list_df_pathway$values)
list_df_pathway
dim(list_df_pathway) 

#glmnet is not able to perform the Efron approximation at the moment. survival's coxph can perform the Breslow approximation by specifying tie = "breslow"
coxph_fit <- coxph(formula = y_tr ~ x_tr[, colnames(x_tr) %in% list_df_pathway$Category], ties = "breslow")
summary(coxph_fit)

coxtidy <- tidy(coxph_fit, exponentiate = T, conf.int = T, conf.level = 0.95)

coxtidy <- cSplit(coxtidy, "term", "]")
coxtidy$estimate_95CI <- paste(round(coxtidy$estimate, 3),"(",round(coxtidy$conf.low, 3),"-",round(coxtidy$conf.high, 3),")")
coxtidy$p.value <- round(coxtidy$p.value, 3)
colnames(coxtidy)

##Figure 3C and Supplementary Figure 2-----------------------
##create a forest plot
tabletext <- coxtidy %>%
  dplyr::select(term_2, p.value, estimate_95CI)

head(tabletext)

header <- c("Risk factor", "P value", "Hazard Ratio (95%CI)")
tabletext <- rbind(header, tabletext)
forestdata <- coxtidy %>%
  dplyr::select(estimate, conf.low, conf.high)
head(forestdata)
header2 <- c(1, 1, 1)
forestdata <- rbind(header2, forestdata)
forestdata_new <- forestdata[c(1:3, 4:29),] 
library(forestplot)
tabletext_new <- tabletext[c(1:3, 4:29),] 

head(tabletext_new)

tiff("Rplot_forestplot_multivariate_clinicalgenetic_new.tiff", units = "in", width = 10, height = 8, res = 300)
p <- forestplot(tabletext_new, 
                graph.pos = 3,
                forestdata_new,new_page = TRUE,
                #is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
                clip=c(-1,4),
                boxsize=0.25,
                zero = 1,
                xlog=FALSE, 
                col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
                xlab = "estimate(95%CI)",
                txt_gp = fpTxtGp(label = list(gpar(fontfamily = "Arial"),
                                              gpar(fontfamily = "",
                                                   col = "#660000")),
                                 ticks = gpar(fontfamily = "", cex = 1.0),
                                 xlab  = gpar(fontfamily = "Arial", cex = 1.0)))

grid::grid.text("3yr Mortality ~ risk factor", .5, 0.99, gp=gpar(cex=1.5))
dev.off()


#Test the Proportional Hazards Assumption of a Cox Regression (Schoenfeld test)
tiff("Rplot_schoenfeld_clinical_genetic.tiff", units="in", width=15, height=8, res=300)
ggcoxdiagnostics(coxph_fit, type = "schoenfeld")
dev.off()

coxph_fit <- coxph(formula = y_tr ~ x_tr, ties = "breslow")
summary(coxph_fit)


colnames(PRS_GNSIS_test_select_removeNA)

importance <- list_df_pathway[abs(list_df_pathway$values) > 0.02, ]$Category
importance
N <- cor(x_tr[, colnames(x_tr) %in% importance])

N <- cor(x_tr_MAF25, x_tr_MAF1) #remove V2007

library(corrplot)
corrplot(N, method="circle")

#identifying correlated predictors
descrCor <- cor(x_tr)
highCorr <- sum(abs(descrCor[upper.tri(descrCor)]) > 0.60)
summary(descrCor[upper.tri(descrCor)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = 0.60)
highlyCorDescr
colnames(x_tr)
filteredDescr <- filteredDescr[, -highlyCorDescr]
descrCor2 <- cor(filteredDescr)
summary(descrCor2[upper.tri(descrCor2)])

## get predictions x^Tbeta
preds <- predict(cvfit_tr,x_tr,s="lambda.1se")
preds <- predict(cvfit_tr,x_te,s="lambda.1se")

preds <- predict(cvfit_tr,x_tr,s="lambda.min")
preds <- predict(cvfit_tr,x_te,s="lambda.min")

as.numeric(preds[,1])

## split into low, medium, high risk groups
library(ggplot2)
levs <- cut_number(preds,3)
head(levs)
unique(levs)

#Make Kaplan-Meier curves and Log Rank p-values for groups
fit <- survfit(y_tr~levs)
out <- survdiff(y_tr~levs)
out

fit <- survfit(y_te~levs)
out <- survdiff(y_te~levs)
out

p.val <- 1 - pchisq(out$chisq, length(out$n) - 1)
p.val

library(ggplot2)

library(ggfortify)
autoplot(fit,xlab="Survival Time (years)", ylab="Survival", conf.int = F, xlim = c(0, 3),
         main=paste0("p-value: ",round(p.val,8)))
autoplot(fit,xlab="Survival Time (years)",ylab="Survival", conf.int = F, xlim = c(0, 3.1), ylim = c(0.65, 1),
         main=paste0("p-value: ",round(p.val,8))) + theme_bw() + theme(axis.line = element_line(colour = "black"),
                                                                       panel.grid.major = element_blank(),
                                                                       panel.grid.minor = element_blank(),
                                                                       #panel.border = element_blank(),
                                                                       panel.background = element_blank())
tr <- data.frame(cbind(y_tr, x_tr))
te <- data.frame(cbind(y_te, x_te))

tiff("Rplot_survival_testing_PRS_MAFALL_05222021_updated_ageonly_multivariate_clinicalgenetic_new.tiff", units="in", width=7, height=8, res=300)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data = te,             # data used to fit survival curves.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 1),
  pval.size = 4,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(11, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 3.1),         # present narrower X axis, but not affect
  ylim = c(0.5, 1),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE, # show bars instead of names in text annotations
  tables.theme = theme_cleantable(),
  # in legend of risk table
  #ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  #ncensor.plot.height = 0.25,
  #conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Low risk", "Intermediate risk", "High risk"),    # change legend labels.
  font.title    = c(11, "bold.italic", "darkgreen"),
  font.subtitle = c(10, "bold", "darkgreen"),
  font.caption  = c(10, "plain", "darkgreen"),
  font.x        = c(11, "bold.italic", "black"),
  font.y        = c(11, "bold.italic", "black"),
  font.xtickslab = c(11, "bold", "black"),
  font.ytickslab = c(11, "bold", "black")
) +
  labs(
    title    = "Survival curves for the testing dataset based on Kaplan-Meier estimates",
    subtitle = "clinical (8) + genetics (21)",
    caption  = "Plotted with R survminer")
dev.off()




# Comparing models --------------------------------------------------------


library(survC1)
mydata <- data.frame(as.matrix(y_te),preds) #change from y_tr to y_te
out <- Est.Cval(mydata, 2000, nofit=TRUE)
cind <- out$Dhat
cind
tau = 365.25*3
C=Inf.Cval(mydata, tau, itr=2000, seed = 1234)
round(c(C$Dhat, C$se, C$low95, C$upp95), digits=3)

## Comaprison of C-statistics
uno.C.delta <- Inf.Cval.Delta(mydata = as.matrix(df_y_te),
                              covs0  = x_te[,38],           # clinical risk factor only model
                              covs1  = x_te[,c(5, 14, 38)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 20, 3, 38)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 18, 17, 12, 26, 8, 19, 2, 27, 1, 4, 38)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 20, 3, 28, 11, 9, 23, 16, 21, 18, 17, 12, 26, 8, 19, 2, 27, 1, 4, 38)],
                              tau    = 3,                       # Trucation time (max time to consider)
                              itr = 2000)                                 # Iteration of perturbation-resampling
uno.C.delta

library(survIDINRI)

res.IDI.INF <- IDI.INF(indata = as.matrix(df_y_te),
                              covs0  = x_te[,32:39],           # clinical risk factor only model
                              covs1  = x_te[,c(5, 14, 32:39)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 32:39)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 18, 12, 26, 8, 19, 32:39)], # integrated model
                              #covs1  = x_te[,c(5, 14, 30, 10, 15, 25, 18, 17, 12, 26, 8, 19, 2, 27, 1, 4, 32:39)], # integrated model
                              t0    = 3,                       # Trucation time (max time to consider)
                              npert = 300, npert.rand = NULL, seed1 = NULL, alpha = 0.05)                                 # Iteration of perturbation-resampling

## M1 IDI; M2 continuous NRI; M3 median improvement
stat <- IDI.INF.OUT(res.IDI.INF)
stat

##Figure 3D and Supplementary Figure 9------------------------
## M1 red area; M2 distance between black points; M3 distance between gray points
tiff("Rplot_IDI_mortality_testing_16.tiff", units="in", width=7, height=8, res=300)
IDI.INF.GRAPH(res.IDI.INF)
title(main = "Improvement in risk prediciton using survIDINRI when comparing clinical(8) + genetic(16) to clinical(8)", sub = "IDI=0.032[0.034-0.134]; NRI=0.275[0.109-0.363]; median improvement = 0.034[0.018-0.124], p < 0.001", 
      # xlab = "X axis", ylab = "Y axis",
      cex.main = 0.75,   font.main= 1, col.main= "red",
      cex.sub = 0.75, font.sub = 1, col.sub = "black",
      col.lab ="darkblue"
)
dev.off()



# Feature QC step ---------------------------------------------------------


##Supplementary Figure 3--------------------------
##density plot for each pathway-specific PRS
tiff("Rplot_densityplot_05172021_update_univariate.tiff", units="in", width=8, height=8, res=300)  ##type = "s" in the following code, can be changed 
PRS_GNSIS_train_select_removeNA[, 2:32] %>%
  #keep(is.numeric) %>%                     # Keep only numeric columns
  tidyr::gather() %>%                             # Convert to key-value pairs
  ggplot(aes(scale(value))) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_density()                         # as density
dev.off()


#This function has been improved to do cox regression on binary variables using maxstat
#short format for 3yr, latest version
PRS_TEST <- data.frame()
pathway <- NULL
PRS_GNSIS <- NULL
cox_3yrdeath_prs1 <- NULL
cox_3yrrecur_prs1 <- NULL
coxtidy_3yrdeath_prs1 <- NULL
coxtidy_3yrrecur_prs1 <- NULL
DEATH <- NULL
RECUR <- NULL
dat_DEATH <- data.frame()
BEST_TEST <- data.frame()
dat_RECUR <- NULL
dat <- NULL
PRS_binary <- NULL
PRS_GNSIS_binary <- NULL

Result_HR = function(PRS_file, FAM_file, COV_file, GNSIS_file) {       # for COV_split file
  pathway <- colnames(PRS_file[, 3:(ncol(PRS_file)-3)])    #select this one for all
  
  for (i in pathway){
    print(i)
    PRS_TEST <- PRS_file %>%
      select(FID, i) %>%
      rename(PRS = i)
    BEST_TEST <- merge(FAM_file[, c("FID", "GROUP")], PRS_TEST[, c("FID", "PRS")], by = "FID")
    BEST_TEST <- merge(BEST_TEST, COV_file, by = "FID")
    #BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP.x==2, ]  #for discovery dataset
    BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP==2, ] #for replication dataset
    #BEST_TEST$PRS_norm <- (BEST_TEST$PRS-min(BEST_TEST$PRS))/(max(BEST_TEST$PRS)-min(BEST_TEST$PRS)) #different normalization method
    colnames(BEST_TEST_CASE)[colnames(BEST_TEST_CASE) == 'FID'] <- 'PT_ID' 
    PRS_GNSIS <- inner_join(BEST_TEST_CASE, GNSIS_file, by = "PT_ID")
    
    #select PT_ID with complete 3_yr observation
    PRS_GNSIS <- PRS_GNSIS %>%
      filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3)) 
    
    #cox regression model
    tryCatch({
      #Creating binary variables
      # conditional Monte-Carlo
      mod <- maxstat.test(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr) ~ PRS, data=PRS_GNSIS,
                          smethod="LogRank", pmethod="condMC", B = 9999)
      PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS <= mod$estimate] <- 0
      PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS > mod$estimate] <- 1
      
      PRS_GNSIS$Category <- i
      #print(PRS_GNSIS)
      PRS_GNSIS_binary[[i]] <- data.frame(PRS_GNSIS[, c("PT_ID", "PRS_norm_cat", "Category")])
      
      #print(PRS_GNSIS) #for QC
      # 3yr Motality stroke ~ PRS_norm_cat
      cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX_BINARY)+ PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrdeath_prs1 <- coxph(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
      coxtidy_3yrdeath_prs1 <- tidy(cox_3yrdeath_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrdeath_prs1$Category <- i
      DEATH[[i]] <- coxtidy_3yrdeath_prs1
      
      #cox.zph(cox_3yrdeath_prs1) # proportional hazards assumptions
      # 3yr Recur stroke ~ PRS_norm_cat
      cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + scale(AGE_AT_INDEX) + PC1 + PC2 + PC3 + PC4 + PC5, data = PRS_GNSIS)
      #cox_3yrrecur_prs1 <- coxph(Surv(timediff_recur_censor_3yr/365.25, event_recur_censor_3yr)~PRS_norm_cat + as.factor(PT_SEX.x) + AGE_AT_INDEX, data = PRS_GNSIS)
      coxtidy_3yrrecur_prs1 <- tidy(cox_3yrrecur_prs1, exponentiate = T, conf.int = T, conf.level = 0.95)
      coxtidy_3yrrecur_prs1$Category <- i
      RECUR[[i]] <- coxtidy_3yrrecur_prs1
      #cox.zph(cox_3yrrecur_prs1) # proportional hazards assumptions
    }, error=function(e){})
  }
  dat_DEATH <- data.frame(do.call(rbind,DEATH))
  dat_RECUR <- data.frame(do.call(rbind,RECUR))
  colnames(dat_DEATH) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  colnames(dat_RECUR) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  dat_DEATH$PHENOTYPE <- "3yr_death"
  dat_RECUR$PHENOTYPE <- "3yr_recur"
  dat <- rbind(dat_DEATH, dat_RECUR)
  PRS_binary <- data.frame(do.call(rbind,PRS_GNSIS_binary))
  save(PRS_binary, file="PRS_binary_train_varypercent_univariate_3yr_adjustcovariate.RData", version = 2)
  #save(PRS_binary, file="PRS_binary_train_test_varypercent_univariate_3yr_adjustcovariate.RData", version = 2)
  return(dat)
}

GNSIS_final_final$AGE_AT_INDEX_BINARY <- ifelse(GNSIS_final_final$AGE_AT_INDEX >= 66.8, 1, 0) #convert AGE_AT_INDEX to binary
head(BEST_split_train_70percent)
Result_GNSIS13_HR_AIS_3yr <- Result_HR(BEST_split_train_test_70percent, FAM_split, COV_split, GNSIS_final_final)  


##Supplementary Figure 4------------------------------
## Identify the cutpoints for each pathway-specific PRS

#short format for 3yr, latest version
dir.create("C:/Users/jli/Desktop/GNSIS/ROCAUC/pathway_PRS_scale/")
PRS_GNSIS <- NULL
cox_3yrdeath_prs1 <- NULL
cox_3yrrecur_prs1 <- NULL
coxtidy_3yrdeath_prs1 <- NULL
coxtidy_3yrrecur_prs1 <- NULL
DEATH <- NULL
RECUR <- NULL
dat_DEATH <- data.frame()
BEST_TEST <- data.frame()
dat_RECUR <- NULL
dat <- NULL
PRS_binary <- NULL
PRS_GNSIS_binary <- NULL

Result_CUTPOINT = function(PRS_file, FAM_file, COV_file, GNSIS_file) {       # for COV_split file
  pathway <- unique(Result_OR_commorbidity_phewas_PRS_remove$PATHWAY)  #select this one for all
  for (i in pathway){
    print(i)
    PRS_TEST <- PRS_file %>%
      select(FID, i) %>%
      rename(PRS = i)
    BEST_TEST <- merge(FAM_file[, c("FID", "GROUP")], PRS_TEST[, c("FID", "PRS")], by = "FID")
    BEST_TEST <- merge(BEST_TEST, COV_file, by = "FID")
    #BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP.x==2, ]  #for discovery dataset
    BEST_TEST_CASE <- BEST_TEST[BEST_TEST$GROUP==2, ] #for replication dataset
    #BEST_TEST$PRS_norm <- (BEST_TEST$PRS-min(BEST_TEST$PRS))/(max(BEST_TEST$PRS)-min(BEST_TEST$PRS)) #different normalization method
    colnames(BEST_TEST_CASE)[colnames(BEST_TEST_CASE) == 'FID'] <- 'PT_ID' 
    PRS_GNSIS <- inner_join(BEST_TEST_CASE, GNSIS_file, by = "PT_ID")
    
    #select PT_ID with complete 3_yr observation
    PRS_GNSIS <- PRS_GNSIS %>%
      filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3)) 
    PRS_GNSIS$PRS <- scale(PRS_GNSIS$PRS)
    #cox regression model
    tryCatch({
      #Creating binary variables
      # conditional Monte-Carlo
      mod <- maxstat.test(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr) ~ PRS, data=PRS_GNSIS,
                          smethod="LogRank", pmethod="condMC", B = 9999)
      PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS <= mod$estimate] <- "PRS_below_threshold"
      PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS > mod$estimate] <- "PRS_above_threshold"
      print(PRS_GNSIS)
      #plotting the PRS screening process
      tiff(paste('TO_MY_PATH/', i, '.tiff', sep=''), units="in", width=10, height=8, res=300)
      #plot(mod)
      plot(mod, main = paste('Pathway = ', unique(Result_OR_commorbidity_phewas_PRS_remove[Result_OR_commorbidity_phewas_PRS_remove$PATHWAY == i, ]$Set), sep=''), col = 'blue') 
      dev.off()
    }, error=function(e){})
  }
}

Result_CUTPOINT(BEST_split_train_test_70percent, FAM_split, COV_split, GNSIS_final_final)  

##Supplementary Figure 5 -----------------
## individual KM curve based on binary pathway-specific PRS converted by the specific cutpoint 

#select PT_ID with complete 3_yr observation
PRS_GNSIS <- PRS_GNSIS %>%
  filter(as.Date(INDEX_DT) < (as.Date("2020-11-15") - 365.25*3)) 
PRS_GNSIS$PRS <- scale(PRS_GNSIS$PRS)
print(PRS_GNSIS)
PRS_GNSIS$PRS
library(maxstat)
#additional parameters to be passed to pmvnorm or B, an integer defining the number of Monte-Carlo replications.
mod <- maxstat.test(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr) ~ PRS, data=PRS_GNSIS,
                    smethod="LogRank", pmethod="condMC", B = 9999)
PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS <= mod$estimate] <- "PRS_below_threshold"
PRS_GNSIS$PRS_norm_cat[PRS_GNSIS$PRS > mod$estimate] <- "PRS_above_threshold"

# KM stratified by cohort
tiff(paste('C:/Users/jli/Desktop/GNSIS/ROCAUC/pathway_PRS_scale_KM_discoverysample/', 'V6899', '.tiff', sep=''), units="in", width=15, height=15, res=600)
ggsurvplot(
  surv_fit(Surv(timediff_death_censor_3yr/365.25, event_death_censor_3yr)~as.factor(PRS_norm_cat), data = PRS_GNSIS), # survfit object with calculated statistics.
  fun = function(x) {1-x},  #plot cumulative probability F(t) = 1 - S(t)
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  pval.coord = c(2, 0.4),
  pval.size = 4,
  conf.int = TRUE,       # show confidence intervals for 
  censor.shape="|", censor.size = 4,
  fontsize = 4, #font size in the table
  font.legend = c(11, "plain", "black"),
  # point estimates of survival curves.
  xlim = c(0, 3.1),         # present narrower X axis, but not affect
  ylim = c(0, 0.4),
  # survival estimates.
  xlab = "Time in years",   # customize X axis label.
  ylab = "Cumulative Incidence",
  break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  n.risk = FALSE,
  legend.labs =
    c("PRS above threshold", "PRS below threshold"),    # change legend labels.
  font.title    = c(11, "bold.italic", "darkgreen"),
  font.subtitle = c(10, "bold", "darkgreen"),
  font.caption  = c(10, "plain", "darkgreen"),
  font.x        = c(11, "bold.italic", "black"),
  font.y        = c(11, "bold.italic", "black"),
  font.xtickslab = c(11, "bold", "black"),
  font.ytickslab = c(11, "bold", "black")
) +
  labs(
    title    = "Cumulative probability for 3yrs mortality",
    subtitle = paste('Pathway = ', unique(Result_OR_commorbidity_phewas_PRS_remove[Result_OR_commorbidity_phewas_PRS_remove$PATHWAY == "V6899", ]$Set)),
    caption  = "Plotted with R survminer")
dev.off()

##Supplementary Figure 6-----------------------
##Correlation among pathway-specific PRS 

library(corrplot)

tiff("Rplot_corrplot_05252021_update_univariate.tiff", units="in", width=8, height=8, res=300)  ##type = "s" in the following code, can be changed 
corrplot(M, method="circle")
dev.off()
pathway_name <- all_pathway[all_pathway$Category %in% colnames(PRS_GNSIS_train_select_removeNA[, 2:32]), ]
save(pathway_name, file = "pathway_name_univariate.RData", version = 2)

pathway_name <- all_pathway[all_pathway$Category %in% colnames(PRS_GNSIS_train_select_removeNA[, 2:26]), ]
save(pathway_name, file = "pathway_name.RData", version = 2)



# Subgroup analysis -------------------------------------------------------


#Figure 4------------------------
## Association of pathway-specific PRS and clinical risk factors
load("PRS_cohort_characteristics_update_pathway.RData")
head(PRS_cohort_characteristics_update_pathway)
PRS_cohort_characteristics_update_pathway_young <- PRS_cohort_characteristics_update_pathway[PRS_cohort_characteristics_update_pathway$AGE_AT_INDEX_BINARY == 0, ]
PRS_cohort_characteristics_update_pathway_old <- PRS_cohort_characteristics_update_pathway[PRS_cohort_characteristics_update_pathway$AGE_AT_INDEX_BINARY == 1, ]

#short format for binary trait
ICD9 <- NULL
mylogit <- NULL
summary_mylogit <- NULL
CI_mylogit <- NULL
summary_mylogit <- data.frame()
Result <- NULL
Result_pathway <- NULL
BEST_commorbidity <- data.frame()
dat_pathway <- NULL
dat <- NULL

Result_commorbidity = function(commorbidity_file) {  #using COV_split      
  ICD9 <- c("hypertension", "diabetes", "dyslipidemia", "smoking", "alcohol", "coronary_artery", "atrial_fib", "alcohol_soc", "BMI_Obesity", "BMI_Overweight", "NIHSS_7above", "NIHSS_10above", "NIHSS_16above")  #using commorbidity file
  PRS <- colnames(commorbidity_file)[46:76]
  for (i in ICD9){
    print(i)
    for (j in PRS) {
      print(j)
      BEST_commorbidity <- commorbidity_file[, c(i, "FID", "PT_SEX", "AGE_AT_INDEX_BINARY", "PC1", "PC2", "PC3", "PC4", "PC5", j)]
      tryCatch({
        #mylogit <- glm(BEST_commorbidity[, 1] ~ as.factor(BEST_commorbidity$PT_SEX) + as.factor(BEST_commorbidity$AGE_AT_INDEX_BINARY) + scale(BEST_commorbidity[[j]]) + PC1 + PC2 + PC3 + PC4 + PC5, data = BEST_commorbidity, family=binomial(link='logit'))
        mylogit <- glm(BEST_commorbidity[, 1] ~ as.factor(BEST_commorbidity$PT_SEX) + scale(BEST_commorbidity[[j]]) + PC1 + PC2 + PC3 + PC4 + PC5, data = BEST_commorbidity, family=binomial(link='logit'))
        
        summary_mylogit <- data.frame(tidy(mylogit, exp = T, conf.int = T, conf.level = 0.95))
        summary_mylogit$PHENOTYPE <- j
        Result[[j]] <- summary_mylogit
      }, error=function(e){})}
    Result_pathway <- as.data.frame(do.call(rbind,Result))
    Result_pathway$PATHWAY <- i
    dat[[i]] <- Result_pathway
  }
  dat_pathway <- as.data.frame(do.call(rbind,dat))
  colnames(dat_pathway) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "PATHWAY", "PHENOTYPE")
  return(dat_pathway)
}

colnames(PRS_cohort_characteristics_update_pathway)
Result_OR_commorbidity_phewas <- Result_commorbidity(PRS_cohort_characteristics_update_pathway)
Result_OR_commorbidity_phewas_PRS <- Result_OR_commorbidity_phewas[Result_OR_commorbidity_phewas$term == "scale(BEST_commorbidity[[j]])", ]
head(Result_OR_commorbidity_phewas_PRS)

#young people
Result_OR_commorbidity_phewas <- Result_commorbidity(PRS_cohort_characteristics_update_pathway_young)
Result_OR_commorbidity_phewas_PRS <- Result_OR_commorbidity_phewas[Result_OR_commorbidity_phewas$term == "scale(BEST_commorbidity[[j]])", ]
head(Result_OR_commorbidity_phewas_PRS)

#old people
Result_OR_commorbidity_phewas <- Result_commorbidity(PRS_cohort_characteristics_update_pathway_old)
Result_OR_commorbidity_phewas_PRS <- Result_OR_commorbidity_phewas[Result_OR_commorbidity_phewas$term == "scale(BEST_commorbidity[[j]])", ]
head(Result_OR_commorbidity_phewas_PRS)

Result_OR_commorbidity_phewas_PRS_remove <- Result_OR_commorbidity_phewas_PRS[Result_OR_commorbidity_phewas_PRS$PHENOTYPE != "alcohol_soc",]
Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select <- Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate[Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate$PHENOTYPE == "3yr_death", ]
Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select$Direction <- log(Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select$estimate)/abs(log(Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select$estimate))
Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select$PATHWAY <- Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select$Category
head(Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select)
Result_OR_commorbidity_phewas_PRS_remove <- merge(Result_GNSIS13_HR_AIS_1yr_3yr_5yr_PRS_candidate_select[, c("PATHWAY", "Direction")], Result_OR_commorbidity_phewas_PRS_remove, by = "PATHWAY")
head(Result_OR_commorbidity_phewas_PRS_remove)
unique(Result_OR_commorbidity_phewas_PRS_remove$PATHWAY)
list_df_pathway_direction_diseaserisk$Header <- list_df_pathway_direction_diseaserisk$Set
list_df_pathway_direction_diseaserisk$PATHWAY <- list_df_pathway_direction_diseaserisk$Category
Result_OR_commorbidity_phewas_PRS_remove <- merge(Result_OR_commorbidity_phewas_PRS_remove, list_df_pathway_direction_diseaserisk[, c("PATHWAY", "Header")], by = "PATHWAY")

tiff("Result_OR_commorbidity_phewas_PRS_remove1_old_NIHSS.tiff", units="in", width=20, height=8, res=600)
ggplot(Result_OR_commorbidity_phewas_PRS_remove, aes(x = as.factor(PHENOTYPE), y = estimate, ymin = conf.low, ymax = conf.high)) + 
  geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  geom_hline(aes(yintercept = 1), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  facet_wrap(. ~ Header.x, scales="free") + 
  theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_blank()) + 
  ylab("Odds ratio & 95% CI") +
  xlab("clinical risk factor for ischemic stroke")
dev.off()

tiff("Result_OR_commorbidity_phewas_PRS_remove1_old_candidate_NIHSS.tiff", units="in", width=20, height=8, res=600)
ggplot(Result_OR_commorbidity_phewas_PRS_remove_candidate, aes(x = as.factor(PHENOTYPE), y = estimate, ymin = conf.low, ymax = conf.high)) + 
  geom_pointrange(aes(col = PHENOTYPE, shape = as.factor(Direction)), position=position_dodge(width=0.8)) + 
  geom_point(aes(shape = as.factor(Direction), color = PHENOTYPE, y = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  #scale_shape_manual(values=c("\u25BA","\u25C4","\u25BC","\u25B2"))
  scale_shape_manual(values=c("\u25BC","\u25B2")) +
  geom_hline(aes(yintercept = 1), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  facet_wrap(. ~ Header.x, scales="free") + 
  theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x=element_blank()) + 
  ylab("Odds ratio & 95% CI") +
  xlab("clinical risk factor for ischemic stroke")
dev.off()

#all 31 pathways
tiff("Result_OR_commorbidity_phewas_PRS_remove_HR_allcases_old.tiff", units="in", width=15, height=10, res=600)
ggplot(Result_OR_commorbidity_phewas_PRS_remove, aes(y = as.factor(PHENOTYPE), x = estimate, xmin = conf.low, xmax = conf.high)) + 
  geom_pointrange(aes(col = PHENOTYPE), position=position_dodge(width=0.8)) + 
  geom_point(aes(color = PHENOTYPE, x = estimate, size = -log10(p.value)), stat="identity", position = position_dodge(0.8))  +
  geom_vline(aes(xintercept = 1), linetype='dashed', col = 'black', width = 1) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme_bw() +
  facet_wrap(. ~ Set, labeller=label_wrap_gen(width=20)) + 
  theme(strip.text.x = element_text(size = 6, colour = "black")) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.y=element_blank()) + 
  xlab("Hazard ratio & 95% CI (n = 1756)") +
  ylab("Observation window")
dev.off()


#Supplementary Figure 10----------------
##logistic regression to test the association between PRSs and clinical risk factors

#short format
commorbidity_select <- data.frame()
ICD9 <- NULL
mylogit <- NULL
summary_mylogit <- NULL
CI_mylogit <- NULL
summary_mylogit <- data.frame()
Result <- NULL
BEST_commorbidity <- data.frame()
BEST_TEST <- data.frame()
commorbidity_select_COV <- data.frame()
dat <- NULL

Result_commorbidity = function(commorbidity_file, PRS_file, COV_file, freq) {  #using COV_split      
  #ICD9 <- colnames(commorbidity_file)[1:(ncol(commorbidity_file)-1)]   #using commorbidity file
  ICD9 <- freq$ICD[400:678]
  ICD9 <- ICD9[ICD9 != "X195"]
  ICD9 <- ICD9[ICD9 != "X318"]
  ICD9 <- ICD9[ICD9 != "X822"]
  ICD9 <- ICD9[ICD9 != "X886"]
  ICD9 <- ICD9[ICD9 != "X933"]
  for (i in ICD9){
    print(i)
    commorbidity_select <- commorbidity_file[, c(i, "FID")]
    commorbidity_select_COV <- merge(commorbidity_select, COV_file[, c("FID", "PT_SEX", "PC1", "PC2", "PC3", "PC4", "PC5")], by = "FID")
    BEST_TEST <- PRS_file %>%  #using BEST_split file
      select(FID, V5425) %>%
      rename(PRS = V5425)
    BEST_commorbidity <- merge(commorbidity_select_COV, BEST_TEST[, c("FID", "PRS")], by = "FID") #1329*9
    BEST_commorbidity$PRS_norm <- scale(BEST_commorbidity$PRS)
    #BEST_commorbidity[, 2] <- as.factor(BEST_commorbidity[, 2])
    #BEST_commorbidity <- BEST_commorbidity[complete.cases(BEST_commorbidity),]
    # print(BEST_commorbidity)
    # print(length(BEST_commorbidity$PRS_norm))
    #BEST_commorbidity[, 2][BEST_commorbidity[, 2] == 0] <- " no disease"
    #BEST_commorbidity[, 2][BEST_commorbidity[, 2] == 1] <- "disease"
    mylogit <- glm(BEST_commorbidity[, 2] ~ BEST_commorbidity$PRS_norm + BEST_commorbidity$PT_SEX + BEST_commorbidity$PC1 + BEST_commorbidity$PC2 + BEST_commorbidity$PC3 + BEST_commorbidity$PC4 + BEST_commorbidity$PC5, data = BEST_commorbidity, family = "binomial")
    summary_mylogit <- data.frame(tidy(mylogit))
    CI_mylogit <- data.frame(confint(mylogit))
    summary_mylogit <- cbind(summary_mylogit, CI_mylogit) 
    summary_mylogit$Category <- i
    Result[[i]] <- summary_mylogit
  }
  dat <- as.data.frame(do.call(rbind,Result))
  print(dat)
  colnames(dat) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  return(dat)
}

Result_OR_LAS_phewas <- Result_commorbidity(commorbidity, BEST_split, COV_split, freq_max)
head(Result_OR_LAS_phewas)

Result_OR_LAS_phewas_PRS <- Result_OR_LAS_phewas[Result_OR_LAS_phewas$term == "BEST_commorbidity$PRS_norm", ]


#short format
commorbidity_select <- data.frame()
ICD9 <- NULL
mylogit <- NULL
summary_mylogit <- NULL
CI_mylogit <- NULL
summary_mylogit <- data.frame()
Result <- NULL
BEST_commorbidity <- data.frame()
BEST_TEST <- data.frame()
commorbidity_select_COV <- data.frame()
dat <- NULL

Result_commorbidity = function(commorbidity_file, PRS_file, COV_file, freq) {  #using COV_split      
  #ICD9 <- colnames(commorbidity_file)[1:(ncol(commorbidity_file)-1)]   #using commorbidity file
  ICD9 <- freq$ICD
  for (i in ICD9){
    print(i)
    commorbidity_select <- commorbidity_file[, c(i, "FID")]
    commorbidity_select_COV <- merge(commorbidity_select, COV_file[, c("FID", "PT_SEX", "PC1", "PC2", "PC3", "PC4", "PC5")], by = "FID")
    BEST_TEST <- PRS_file %>%  #using BEST_split file
      select(FID, V5425) %>%
      rename(PRS = V5425)
    BEST_commorbidity <- merge(commorbidity_select_COV, BEST_TEST[, c("FID", "PRS")], by = "FID") #1329*9
    BEST_commorbidity$PRS_norm <- scale(BEST_commorbidity$PRS)
    #BEST_commorbidity[, 2] <- as.factor(BEST_commorbidity[, 2])
    #BEST_commorbidity <- BEST_commorbidity[complete.cases(BEST_commorbidity),]
    # print(BEST_commorbidity)
    # print(length(BEST_commorbidity$PRS_norm))
    #BEST_commorbidity[, 2][BEST_commorbidity[, 2] == 0] <- " no disease"
    #BEST_commorbidity[, 2][BEST_commorbidity[, 2] == 1] <- "disease"
    tryCatch({
      mylogit <- glm(BEST_commorbidity[, 2] ~ BEST_commorbidity$PRS_norm + BEST_commorbidity$PT_SEX + BEST_commorbidity$PC1 + BEST_commorbidity$PC2 + BEST_commorbidity$PC3 + BEST_commorbidity$PC4 + BEST_commorbidity$PC5, data = BEST_commorbidity, family = "binomial")
      summary_mylogit <- data.frame(tidy(mylogit))
      CI_mylogit <- data.frame(confint(mylogit))
      summary_mylogit <- cbind(summary_mylogit, CI_mylogit) 
      summary_mylogit$Category <- i
      Result[[i]] <- summary_mylogit
    }, error=function(e){})
  }
  dat <- as.data.frame(do.call(rbind,Result))
  print(dat)
  colnames(dat) <- c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high", "Category")
  return(dat)
}

Result_OR_LAS_phewas <- Result_commorbidity(commorbidity, BEST_split, COV_split, freq_max)
head(Result_OR_LAS_phewas)

Result_OR_LAS_phewas_PRS <- Result_OR_LAS_phewas[Result_OR_LAS_phewas$term == "BEST_commorbidity$PRS_norm", ]
save(Result_OR_LAS_phewas_PRS, file = "Result_OR_SVS_phewas_PRS_ALL.RData")

