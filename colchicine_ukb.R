gc()
rm(list=ls())

library(dplyr)
library(tidyverse)
library(data.table)
#library(ukbtools) 
library(survival)
library(survminer)
library(ggrepel)
library(gridExtra)
#library(plinkQC)
library(tableone)
library(stringr)

#########################################
###############UKBB############################
##############################################
###read data
#colchicine
colchicine_1=fread("/medpop/esp2/btruong/Projects/tentative/colchicine/data/colchicine_UKB_230525.txt")
colchicine_1=colchicine_1%>%select(eid, colchicine=colchicine_all)
table(colchicine_1$colchicine, useNA="always")

colchicine_2=fread("/medpop/esp2/SarahUrbut/colchicine.txt")
colchicine_2=colchicine_2%>%mutate(colchicine=1)%>%select(eid=V1, colchicine, date=V3)
# Convert date column to Date type
colchicine_2$date <- as.Date(colchicine_2$date, format = "%d/%m/%Y")

# Filter for the earliest date per participant
colchicine_2 <- colchicine_2 %>% arrange(eid, date) %>% group_by(eid) %>% slice(1)

check=intersect(colchicine_2$eid, colchicine_1[which(colchicine_1$colchicine==1),]$eid)
length(colchicine_1[which(colchicine_1$colchicine==1),]$eid)
colchicine_2 <- colchicine_2[!colchicine_2$eid %in% check,]

colchicine_1=colchicine_1%>%mutate(date=NA, baseline="yes")%>%filter(colchicine==1)
colchicine_2=colchicine_2%>%mutate(baseline="no")

# Convert the date columns to Date objects
colchicine_1$date <- as.Date(colchicine_1$date, format = "%d/%m/%Y")
colchicine_2$date <- as.Date(colchicine_2$date, format = "%d/%m/%Y")

# Then try the rbind operation
colchicine <- rbind(colchicine_1, colchicine_2)

# Read in Maryam's master pheno file
pheno = fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.REAL.txt")
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) |  pheno$Non_Consented== 1),1,0) ## N:
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

pheno = pheno%>%select(id, in_white_British_ancestry_subset, age, Sex, contains("smok"), contains("Coronary_Artery_Disease"), PC1:PC10, 
genotyping_array, contains("composite_mi_cad_stroke"), contains("Myocardial_Infarction"), contains("heme"), BMI, contains("Diabetes_Type_2"), 
contains("Gout"), contains("Ischemic_stroke"), contains("leukemia"), contains("Leukemia"), contains("anycancer"), contains("Hypertension"), 
contains("Chronic_kidney_disease"), father_heart_disease, mother_heart_disease, sibling_heart_disease, family_disease_history, Haematocrit.percentage, Mean.corpuscular.volume,
Red.blood.cell..erythrocyte..distribution.width, Mean.corpuscular.haemoglobin.concentration, Platelet.count,Red.blood.cell..erythrocyte..count,
White.blood.cell..leukocyte..count, Haemoglobin.concentration, Total.Cholesterol, HDL.Cholesterol,Total.Triglycerides,LDL.Cholesterol)

pheno_2023=fread("/medpop/esp2/zyu/others/ukb_pheno_curation/CAD_Phenotypes_All_Participants_use.csv")
pheno_2023_more=fread("/medpop/esp2/zyu/others/ukb_pheno_curation/CAD_Phenotypes_All_Participants_hardcomposite_use.csv")
pheno_2023_more=pheno_2023_more%>%select(f.eid, contains("Composite"))
pheno_2023=merge(pheno_2023, pheno_2023_more, by="f.eid")

pheno_2023=pheno_2023%>%mutate(Prev_all_cause_mortality=0)%>%select(id=f.eid, Prev_all_cause_mortality, Inc_all_cause_mortality=all_cause_mortality, 
Years_To_all_cause_mortality, Years_To_CAD_HARD, Prev_CAD_HARD, Inc_CAD_HARD, Any_CAD_HARD, Years_To_CAD_INTERMEDIATE, 
Prev_CAD_INTERMEDIATE, Inc_CAD_INTERMEDIATE, Any_CAD_INTERMEDIATE, Years_To_CAD_SOFT, Prev_CAD_SOFT, Inc_CAD_SOFT, 
Any_CAD_SOFT, Years_To_Heart_Failure, Prev_Heart_Failure, Inc_Heart_Failure, 
Any_Heart_Failure, Years_To_Ischemic_Stroke, Prev_Ischemic_Stroke, Inc_Ischemic_Stroke, Any_Ischemic_Stroke, Years_To_MI, Prev_MI, 
Inc_MI, Any_MI, contains("Composite"), time_to_follow_up)
pheno_2023=pheno_2023%>%mutate(Any_all_cause_mortality=Inc_all_cause_mortality)

# Rename the columns, adding "_2023" to all columns except the first one
pheno_2023 <-pheno_2023 %>% 
  rename_with(~if_else(. == names(pheno_2023)[1], ., paste0(., "_2023")), -all_of(names(pheno_2023)[1]))

pheno_2023 <- pheno_2023 %>% 
  rename_with(~str_replace_all(., c("Years_To_" = "FollowUp_", "Inc_" = "Incd_")))

# Get column names that start with "FollowUp_"
followup_cols <- grep("^FollowUp_", names(pheno_2023), value = TRUE)

# Get column names that start with "Prev_"
prev_cols <- grep("^Prev_", names(pheno_2023), value = TRUE)

# Get column names that start with "Incd_"
incd_cols <- grep("^Incd_", names(pheno_2023), value = TRUE)

# Loop through columns
for(i in seq_along(followup_cols)) {
  # Replace NA values with 'time_to_follow_up_2023' in 'FollowUp_' columns
  na_indices <- which(is.na(pheno_2023[[followup_cols[i]]]))
  pheno_2023[na_indices, followup_cols[i]] <- pheno_2023[na_indices, "time_to_follow_up_2023"]

  # Find indices where 'FollowUp_' columns have negative values and 'Prev_' columns have 1
  ind <- which(pheno_2023[[prev_cols[i]]] == 1)

  # Make the negative values in 'FollowUp_' as NA and the corresponding 'Incd_' variables as NA
  pheno_2023[ind, followup_cols[i]] <- NA
  pheno_2023[ind, incd_cols[i]] <- NA
}

pheno=merge(pheno, pheno_2023, by="id")

pheno=merge(pheno, colchicine, by.x="id", by.y="eid", all.x=T)

firstdate=fread("/medpop/esp2/projects/UK_Biobank/baskets/ukb672714/ukb672714.tab.gz")
firstdate=firstdate%>%select(eid=f.eid, firstdate=f.53.0.0)
pheno=merge(pheno, firstdate, by.x="id", by.y="eid", all.x=T)

load("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda")
newchip_overall=ukb450k_ch
newchip_overall[which(newchip_overall$CHIP==1),][is.na(newchip_overall[which(newchip_overall$CHIP==1),])] <- 0

pheno_chip=merge(pheno, newchip_overall, by.x="id", by.y="eid_7089")
pheno_chip$colchicine[is.na(pheno_chip$colchicine)] <- 0

pheno_chip$firstdate <- as.Date(pheno_chip$firstdate, format = "%Y-%m-%d")
pheno_chip$colchicine[is.na(pheno_chip$colchicine)] <- 0

#################################################CHIP caitlyn###############################
chip_caitlyn=fread("/medpop/esp2/tnakao/data/UKBB/phenotype/CHIP/450k/2022Oct/CHIP_calls_Oct16_2022.txt")

CHIP_5= unique(chip_caitlyn[chip_caitlyn$ minAD>=5,]$Broad_ID)
DNMT3A_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "DNMT3A"& chip_caitlyn$ minAD>=5,]$Broad_ID)
non_DNMT3A_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene != "DNMT3A"& chip_caitlyn$ minAD>=5,]$Broad_ID)
TET2_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "TET2" & chip_caitlyn$ minAD>=5,]$Broad_ID)
ASXL1_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "ASXL1"& chip_caitlyn$ minAD>=5,]$Broad_ID)

CHIP_10_5= unique(chip_caitlyn[chip_caitlyn$ minAD>=5 &chip_caitlyn$AF>=0.1,]$Broad_ID)
DNMT3A_10_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "DNMT3A"& chip_caitlyn$ minAD>=5 &chip_caitlyn$AF>=0.1,]$Broad_ID)
non_DNMT3A_10_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene != "DNMT3A"& chip_caitlyn$ minAD>=5 &chip_caitlyn$AF>=0.1,]$Broad_ID)
TET2_10_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "TET2" & chip_caitlyn$ minAD>=5 &chip_caitlyn$AF>=0.1,]$Broad_ID)
ASXL1_10_5 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "ASXL1"& chip_caitlyn$ minAD>=5 &chip_caitlyn$AF>=0.1,]$Broad_ID)

CHIP_3= unique(chip_caitlyn[chip_caitlyn$ minAD>=3,]$Broad_ID)
DNMT3A_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "DNMT3A"& chip_caitlyn$ minAD>=3,]$Broad_ID)
non_DNMT3A_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene != "DNMT3A"& chip_caitlyn$ minAD>=3,]$Broad_ID)
TET2_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "TET2" & chip_caitlyn$ minAD>=3,]$Broad_ID)
ASXL1_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "ASXL1"& chip_caitlyn$ minAD>=3,]$Broad_ID)

CHIP_10_3= unique(chip_caitlyn[chip_caitlyn$ minAD>=3 &chip_caitlyn$AF>=0.1,]$Broad_ID)
DNMT3A_10_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "DNMT3A"& chip_caitlyn$ minAD>=3 &chip_caitlyn$AF>=0.1,]$Broad_ID)
non_DNMT3A_10_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene != "DNMT3A"& chip_caitlyn$ minAD>=3 &chip_caitlyn$AF>=0.1,]$Broad_ID)
TET2_10_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "TET2" & chip_caitlyn$ minAD>=3 &chip_caitlyn$AF>=0.1,]$Broad_ID)
ASXL1_10_3 = unique(chip_caitlyn[chip_caitlyn$Gene.refGene == "ASXL1"& chip_caitlyn$ minAD>=3 &chip_caitlyn$AF>=0.1,]$Broad_ID)

###################read in Caityln
pheno_chip$CHIP_5 <- ifelse(pheno_chip$id %in% CHIP_5,1,0)
pheno_chip$DNMT3A_5 <- ifelse(pheno_chip$id %in% DNMT3A_5,1,0)
pheno_chip$non_DNMT3A_5 <- ifelse(pheno_chip$id %in% non_DNMT3A_5,1,0)
pheno_chip$TET2_5 <- ifelse(pheno_chip$id %in% TET2_5,1,0)
pheno_chip$ASXL1_5 <- ifelse(pheno_chip$id %in% ASXL1_5,1,0)

pheno_chip$CHIP_3 <- ifelse(pheno_chip$id %in% CHIP_3,1,0)
pheno_chip$DNMT3A_3 <- ifelse(pheno_chip$id %in% DNMT3A_3,1,0)
pheno_chip$non_DNMT3A_3 <- ifelse(pheno_chip$id %in% non_DNMT3A_3,1,0)
pheno_chip$TET2_3 <- ifelse(pheno_chip$id %in% TET2_3,1,0)
pheno_chip$ASXL1_3 <- ifelse(pheno_chip$id %in% ASXL1_3,1,0)

pheno_chip$CHIP_10_5 <- ifelse(pheno_chip$id %in% CHIP_10_5,1,0)
pheno_chip$DNMT3A_10_5 <- ifelse(pheno_chip$id %in% DNMT3A_10_5,1,0)
pheno_chip$non_DNMT3A_10_5 <- ifelse(pheno_chip$id %in% non_DNMT3A_10_5,1,0)
pheno_chip$TET2_10_5 <- ifelse(pheno_chip$id %in% TET2_10_5,1,0)
pheno_chip$ASXL1_10_5 <- ifelse(pheno_chip$id %in% ASXL1_10_5,1,0)

pheno_chip$CHIP_10_3 <- ifelse(pheno_chip$id %in% CHIP_10_3,1,0)
pheno_chip$DNMT3A_10_3 <- ifelse(pheno_chip$id %in% DNMT3A_10_3,1,0)
pheno_chip$non_DNMT3A_10_3 <- ifelse(pheno_chip$id %in% non_DNMT3A_10_3,1,0)
pheno_chip$TET2_10_3 <- ifelse(pheno_chip$id %in% TET2_10_3,1,0)
pheno_chip$ASXL1_10_3 <- ifelse(pheno_chip$id %in% ASXL1_10_3,1,0)

#names(chip_caitlyn)[2]="Gene"
#chip_caitlyn_DNMT3A=chip_caitlyn%>%filter(Gene=="DNMT3A")
#chip_caitlyn_DNMT3A=chip_caitlyn_DNMT3A%>%select(Broad_ID, AF)%>% group_by(Broad_ID) %>% summarise(DNMT3A_VAF=sum(AF))
#pheno_chip=merge(pheno_chip, chip_caitlyn_DNMT3A, by.x="id", by.y="Broad_ID", all.x=T)
#pheno_chip=pheno_chip%>%mutate(DNMT3A_VAF=ifelse(DNMT3A_3==0, 0, ifelse(DNMT3A_3==1, DNMT3A_VAF, NA)))
#chip_caitlyn_TET2=chip_caitlyn%>%filter(Gene=="TET2")
#chip_caitlyn_TET2=chip_caitlyn_TET2%>%select(Broad_ID, AF)%>% group_by(Broad_ID) %>% summarise(TET2_VAF=sum(AF))
#pheno_chip=merge(pheno_chip, chip_caitlyn_TET2, by.x="id", by.y="Broad_ID", all.x=T)
#pheno_chip=pheno_chip%>%mutate(TET2_VAF=ifelse(TET2_3==0, 0, ifelse(TET2_3==1, TET2_VAF, NA)))
#chip_caitlyn_ASXL1=chip_caitlyn%>%filter(Gene=="ASXL1")
#chip_caitlyn_ASXL1=chip_caitlyn_ASXL1%>%select(Broad_ID, AF)%>% group_by(Broad_ID) %>% summarise(ASXL1_VAF=sum(AF))
#pheno_chip=merge(pheno_chip, chip_caitlyn_ASXL1, by.x="id", by.y="Broad_ID", all.x=T)
#pheno_chip=pheno_chip%>%mutate(ASXL1_VAF=ifelse(ASXL1_3==0, 0, ifelse(ASXL1_3==1, ASXL1_VAF, NA)))

###################Exclude relatedness######################################
bd_related_excluded=fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt")
pheno_chip_unrelated=pheno_chip[!pheno_chip$id %in% bd_related_excluded$IID,]
pheno_chip=pheno_chip_unrelated

########################Table 1
dim(pheno_chip[which(pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_MI_2023==0 ),])

factorVars <- c("CHIP_3", "DNMT3A_3","TET2_3","ASXL1_3", "in_white_British_ancestry_subset", "Sex", "ever_smoked", 
"Prev_Diabetes_Type_2", "Prev_Gout", "Incd_MI_2023",  "Prev_anycancer", "Prev_Hypertension", "Prev_Chronic_kidney_disease", 
 "family_disease_history")

vars <- c("CHIP_3", "DNMT3A_3","TET2_3","ASXL1_3", "in_white_British_ancestry_subset", "Sex", "ever_smoked", 
"Prev_Diabetes_Type_2", "Prev_Gout", "Incd_MI_2023", "age", "Prev_anycancer", "Prev_Hypertension", "Prev_Chronic_kidney_disease", 
 "family_disease_history", "Haematocrit.percentage", "Mean.corpuscular.volume",
"Red.blood.cell..erythrocyte..distribution.width", "Mean.corpuscular.haemoglobin.concentration", "Platelet.count","Red.blood.cell..erythrocyte..count",
"White.blood.cell..leukocyte..count", "Haemoglobin.concentration", "Total.Cholesterol", "HDL.Cholesterol","Total.Triglycerides","LDL.Cholesterol", "BMI")

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),],factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables_nostrata.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$CHIP_3 ==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables_CHIP.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$DNMT3A_3 ==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables_DNMT3A.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$TET2_3 ==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables_TET2.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$ASXL1_3 ==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_ukb_unrealted_noprevdisease_2023_removeleukemia_morevariables_ASXL1.txt", sep="\t", quote=FALSE, row.names=TRUE)

vars <- c("TET2_VAF")
tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$TET2_3==1 & pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),], strata="colchicine")
print(tableOne)
print(tableOne,nonnormal = vars)

vars <- c("DNMT3A_VAF")
tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$DNMT3A_3==1 & pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),], strata="colchicine")
print(tableOne)
print(tableOne,nonnormal = vars)

vars <- c("ASXL1_VAF")
tableOne <- CreateTableOne(vars = vars, data = pheno_chip[which(pheno_chip$ASXL1_3==1 & pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),], strata="colchicine")
print(tableOne)
print(tableOne,nonnormal = vars)

###check follow-up time
summary(pheno_chip[which(pheno_chip$Prev_MI_2023==0 & pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),]$FollowUp_MI_2023)

###################KM plot###########################
pheno_chip$tet2colchicine=ifelse(pheno_chip$TET2_3==0&pheno_chip$colchicine==0,0,
			ifelse(pheno_chip$TET2_3==1&pheno_chip$colchicine==0,1,
				ifelse(pheno_chip$TET2_3==0&pheno_chip$colchicine==1,2, 
					ifelse(pheno_chip$TET2_3==1&pheno_chip$colchicine==1,3, NA))))

pdf(file = '/medpop/esp2/zyu/others/tet2colchicine/output/UKB_KPcurve.pdf',width=7, height=7) 
fit <- survfit(Surv(FollowUp_MI_2023, Incd_MI_2023) ~ tet2colchicine, data = pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0),])
gg=ggsurvplot(fit, pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0),], conf.int = F, ggtheme = theme_classic(), palette = c("#BABABA","#404040","#F4A582","#CA0020"),  xlab = "Years", fun = "event", 
ylab="Cumulative event rate", ylim=c(0,0.12),surv.scale="percent", censor=F, legend=c(0.3,0.9), legend.title=element_blank(),
legend.labs=c("No TET2, Colchicine", "TET2, No Colchicine", "No TET2, No Colchicine", "TET2, Colchicine"), font.legend = list(size = 16, color = "black"),font.main = c(16,"plain", "black"),
   font.x = c(16, "plain", "black"),
   font.y = c(16, "plain", "black"),
   font.tickslab = c(16, "plain", "black"))
gg$plot
dev.off()

pheno_chip$tet2colchicine <- ifelse(pheno_chip$TET2_3 == 0 & pheno_chip$colchicine == 1, 0,
                            ifelse(pheno_chip$TET2_3 == 1 & pheno_chip$colchicine == 0, 1,
                            ifelse(pheno_chip$TET2_3 == 0 & pheno_chip$colchicine == 0, 2,
                            ifelse(pheno_chip$TET2_3 == 1 & pheno_chip$colchicine == 1, 3, NA))))

pdf(file = '/medpop/esp2/zyu/others/tet2colchicine/output/UKB_KPcurve.pdf', width=7, height=7) 
fit <- survfit(Surv(FollowUp_MI_2023, Incd_MI_2023) ~ tet2colchicine, data = pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$Prev_Myeloid_leukemia==0),])
gg <- ggsurvplot(fit, data = pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$Prev_Myeloid_leukemia==0),], conf.int = FALSE, ggtheme = theme_classic(),
                 palette = c("#404040", "#F4A582", "#BABABA", "#CA0020"), xlab = "Years", fun = "event", ylab = "Cumulative event rate", ylim = c(0, 0.12), surv.scale = "percent",
                 censor = FALSE, legend = c(0.3, 0.9), legend.title = element_blank(),
                 legend.labs = c("No TET2, Colchicine", "TET2, No Colchicine", "No TET2, No Colchicine", "TET2, Colchicine"),
                 font.legend = list(size = 16, color = "black"), font.main = c(16, "plain", "black"),
                 font.x = c(16, "plain", "black"), font.y = c(16, "plain", "black"), font.tickslab = c(16, "plain", "black"))
print(gg$plot)
dev.off()

##############################################################
###############no interaction##################################
#pheno_chip=pheno_chip[which(pheno_chip$age>50),]

chip_new=c("CHIP_3", "DNMT3A_3","TET2_3","ASXL1_3")

chip_fulllist=c(chip_new)

outcome_list=c("MI_2023")
#"Coronary_Artery_Disease_INTERMEDIATE_2023", "Coronary_Artery_Disease_HARD_2023", "Ischemic_stroke_2023"


h<-data.frame(chip=character(), outcome=character(), n=double(), hr=double(), se=double(), z=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

for(m in 1:length(outcome_list)) {
	for(i in 1:length(chip_fulllist)) {
      		h[i+(m-1)*length(chip_fulllist),"chip"]<-chip_fulllist[i]
		h[i+(m-1)*length(chip_fulllist),"outcome"]<-outcome_list[m]
		Incd_outcome<-paste("Incd_", outcome_list[m], sep="")
      		fmla<-as.formula(paste("Surv(as.numeric(FollowUp_", outcome_list[m], "),Incd_", outcome_list[m], ") ~ ", chip_fulllist[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout+ever_smoked",sep=""))
      		if(!is.na(Incd_outcome)){
			f<-coxph(fmla,data=pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0),])
      			h[i+(m-1)*length(chip_fulllist),"n"]<-summary(f)$n
      			h[i+(m-1)*length(chip_fulllist),"hr"]<-summary(f)$coefficient[1,2]
      			h[i+(m-1)*length(chip_fulllist),"se"]<-summary(f)$coefficient[1,3]
			h[i+(m-1)*length(chip_fulllist),"z"]<-summary(f)$coefficient[1,4]
      			h[i+(m-1)*length(chip_fulllist),"pval"]<-summary(f)$coefficient[1,5]
      			h[i+(m-1)*length(chip_fulllist),"lci"]<-exp(summary(f)$coefficient[1,1]-1.96*summary(f)$coefficient[1,3])
      			h[i+(m-1)*length(chip_fulllist),"uci"]<-exp(summary(f)$coefficient[1,1]+1.96*summary(f)$coefficient[1,3])
		} 	
	}
}
write.table(h, "/medpop/esp2/zyu/others/checkforFuster/output/final_inci_ukb_chip_CVDs_v4_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	

##################################interaction#############################
#by CHIP
h_chip = data.frame(chip=character(), outcome=character(), n_nochip=double(), hr_nochip=double(), se_nochip=double(), z_nochip=double(), pval_nochip=double(), lci_nochip=double(), uci_nochip=double(), 
               n_chip=double(), hr_chip=double(), se_chip=double(), z_chip=double(), pval_chip=double(),lci_chip=double(), uci_chip=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_chip[index,"chip"] <- chip
	      h_chip[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(FollowUp_", outcome, ",Incd_", outcome, ") ~ ", "colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout+ever_smoked",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==0& pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0 ),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==1& pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0 ),])
        fmla_inter = as.formula(paste("Surv(FollowUp_", outcome, ",Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout+ever_smoked",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 & pheno_chip$Prev_Myeloid_leukemia==0 ),])
        h_chip[index,"n_nochip"] <- summary(f1)$n
        h_chip[index,"hr_nochip"] <- summary(f1)$coefficient[1,2]
        h_chip[index,"se_nochip"] <- summary(f1)$coefficient[1,3]
        h_chip[index,"z_nochip"] <- summary(f1)$coefficient[1,4]
        h_chip[index,"pval_nochip"] <- summary(f1)$coefficient[1,5]
        h_chip[index,"lci_nochip"] <- exp(confint(f1)[1,1])
        h_chip[index,"uci_nochip"] <- exp(confint(f1)[1,2])
        h_chip[index,"n_chip"] <- summary(f2)$n
        h_chip[index,"hr_chip"] <- summary(f2)$coefficient[1,2]
        h_chip[index,"se_chip"] <- summary(f2)$coefficient[1,3]
        h_chip[index,"z_chip"] <- summary(f2)$coefficient[1,4]
        h_chip[index,"pval_chip"] <- summary(f2)$coefficient[1,5]
        h_chip[index,"lci_chip"] <- exp(confint(f2)[1,1])
        h_chip[index,"uci_chip"] <- exp(confint(f2)[1,2])
        h_chip[index,"beta_inter"] <- summary(f3)$coefficient[19,2]
        h_chip[index,"se_inter"] <- summary(f3)$coefficient[19,3]
        h_chip[index,"z_inter"] <- summary(f3)$coefficient[19,4]
        h_chip[index,"pval_inter"] <- summary(f3)$coefficient[19,5]
        index <- index + 1
    }
}

write.table(h_chip, "/medpop/esp2/zyu/others/checkforFuster/output/final_inci_ukb_chipcolchicine_stratify_withmorecolchicine_v4_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	

#by colchicine
h_colchicine = data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), hr_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), hr_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_colchicine[index,"chip"] <- chip
	      h_colchicine[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout+ever_smoked",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==0& pheno_chip$Prev_AnyHemeCa==0),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==1& pheno_chip$Prev_AnyHemeCa==0),])
        fmla_inter = as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout+ever_smoked",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_Lymphoid_leukemia==0 &pheno_chip$Prev_Myeloid_leukemia==0),])
        h_colchicine[index,"n_nocolchicine"] <- summary(f1)$n
        h_colchicine[index,"hr_nocolchicine"] <- summary(f1)$coefficient[1,2]
        h_colchicine[index,"se_nocolchicine"] <- summary(f1)$coefficient[1,3]
        h_colchicine[index,"z_nocolchicine"] <- summary(f1)$coefficient[1,4]
        h_colchicine[index,"pval_nocolchicine"] <- summary(f1)$coefficient[1,5]
        h_colchicine[index,"lci_nocolchicine"] <- exp(confint(f1)[1,1])
        h_colchicine[index,"uci_nocolchicine"] <- exp(confint(f1)[1,2])
        h_colchicine[index,"n_colchicine"] <- summary(f2)$n
        h_colchicine[index,"hr_colchicine"] <- summary(f2)$coefficient[1,2]
        h_colchicine[index,"se_colchicine"] <- summary(f2)$coefficient[1,3]
        h_colchicine[index,"z_colchicine"] <- summary(f2)$coefficient[1,4]
        h_colchicine[index,"pval_colchicine"] <- summary(f2)$coefficient[1,5]
        h_colchicine[index,"lci_colchicine"] <- exp(confint(f2)[1,1])
        h_colchicine[index,"uci_colchicine"] <- exp(confint(f2)[1,2])
        h_colchicine[index,"beta_inter"] <- summary(f3)$coefficient[19,2]
        h_colchicine[index,"se_inter"] <- summary(f3)$coefficient[19,3]
        h_colchicine[index,"z_inter"] <- summary(f3)$coefficient[19,4]
        h_colchicine[index,"pval_inter"] <- summary(f3)$coefficient[19,5]
        index <- index + 1
  }
}

write.table(h_colchicine, "/medpop/esp2/zyu/others/checkforFuster/output/final_inci_ukb_chipcolchicine_stratifycolchicine_withmorecolchicine_v4_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	

########################Scan################################

##############################################################
###############no interaction##################################
chip_new=c("CHIP_3", "DNMT3A_3","TET2_3","ASXL1_3")
#, "CHIP_5", "DNMT3A_5","TET2_5","ASXL1_5","CHIP_10_5", "DNMT3A_10_5",
#"TET2_10_5","ASXL1_10_5", "CHIP_10_3", "DNMT3A_10_3","TET2_10_3","ASXL1_10_3")
#chip_list=c("CHIP", "DNMT3A", "TET2", "ASXL1")
#chip_large_list=c("expandedCHIP", "expandedDNMT3A", "expandedTET2", "expandedASXL1")

#chip_fulllist=c(chip_new, chip_list, chip_large_list)
chip_fulllist=chip_new

#outcome_list=c("composite_mi_cad_stroke_2023", "Myocardial_Infarction", "Coronary_Artery_Disease_SOFT", "Coronary_Artery_Disease_INTERMEDIATE", "Coronary_Artery_Disease_HARD", 
#"Ischemic_stroke", "composite_mi_cad_stroke_2023", "Myocardial_Infarction_2023", "Coronary_Artery_Disease_SOFT_2023", "Coronary_Artery_Disease_INTERMEDIATE_2023", "Coronary_Artery_Disease_HARD_2023", 
#"Ischemic_stroke_2023")

outcome_list=c("all_cause_mortality_2023", "CAD_HARD_2023", "Heart_Failure_2023", "Ischemic_Stroke_2023", "MI_2023", "Composite_MI_CAD_HARD_Ischemic_Stroke_all_cause_mortality_2023")

#"CAD_INTERMEDIATE_2023", "CAD_SOFT_2023", 
#"Composite_MI_CAD_INTERMEDIATE_2023", "Composite_MI_CAD_INTERMEDIATE_Ischemic_Stroke_2023", 
#"Composite_MI_CAD_INTERMEDIATE_all_cause_mortality_2023", "Composite_MI_CAD_INTERMEDIATE_Ischemic_Stroke_all_cause_mortality_2023", 
#"Composite_MI_CAD_INTERMEDIATE_Ischemic_Stroke_all_cause_mortality_2023", "Composite_MI_CAD_HARD_2023", "Composite_MI_CAD_HARD_Ischemic_Stroke_2023", 
#"Composite_MI_CAD_HARD_all_cause_mortality_2023", "Composite_MI_CAD_HARD_Ischemic_Stroke_all_cause_mortality_2023")

h<-data.frame(chip=character(), outcome=character(), n=double(), hr=double(), se=double(), z=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

for(m in 1:length(outcome_list)) {
	for(i in 1:length(chip_fulllist)) {
      		h[i+(m-1)*length(chip_fulllist),"chip"]<-chip_fulllist[i]
		h[i+(m-1)*length(chip_fulllist),"outcome"]<-outcome_list[m]
		Incd_outcome<-paste("Incd_", outcome_list[m], sep="")
      		fmla<-as.formula(paste("Surv(as.numeric(FollowUp_", outcome_list[m], "),Incd_", outcome_list[m], ") ~ ", chip_fulllist[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
      		if(!is.na(Incd_outcome)){
			f<-coxph(fmla,data=pheno_chip[which(pheno_chip$Prev_AnyHemeCa==0),])
      			h[i+(m-1)*length(chip_fulllist),"n"]<-summary(f)$n
      			h[i+(m-1)*length(chip_fulllist),"hr"]<-summary(f)$coefficient[1,2]
      			h[i+(m-1)*length(chip_fulllist),"se"]<-summary(f)$coefficient[1,3]
			h[i+(m-1)*length(chip_fulllist),"z"]<-summary(f)$coefficient[1,4]
      			h[i+(m-1)*length(chip_fulllist),"pval"]<-summary(f)$coefficient[1,5]
      			h[i+(m-1)*length(chip_fulllist),"lci"]<-exp(summary(f)$coefficient[1,1]-1.96*summary(f)$coefficient[1,3])
      			h[i+(m-1)*length(chip_fulllist),"uci"]<-exp(summary(f)$coefficient[1,1]+1.96*summary(f)$coefficient[1,3])
		} 	
	}
}
write.table(h, "/medpop/esp2/zyu/others/checkforFuster/inci_ukb_chip_CVDs_v5_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	

###remove leukemia
h<-data.frame(chip=character(), outcome=character(), n=double(), hr=double(), se=double(), z=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

for(m in 1:length(outcome_list)) {
	for(i in 1:length(chip_fulllist)) {
      		h[i+(m-1)*length(chip_fulllist),"chip"]<-chip_fulllist[i]
		h[i+(m-1)*length(chip_fulllist),"outcome"]<-outcome_list[m]
		Incd_outcome<-paste("Incd_", outcome_list[m], sep="")
      		fmla<-as.formula(paste("Surv(as.numeric(FollowUp_", outcome_list[m], "),Incd_", outcome_list[m], ") ~ ", chip_fulllist[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
      		if(!is.na(Incd_outcome)){
			f<-coxph(fmla,data=pheno_chip[which(pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
      			h[i+(m-1)*length(chip_fulllist),"n"]<-summary(f)$n
      			h[i+(m-1)*length(chip_fulllist),"hr"]<-summary(f)$coefficient[1,2]
      			h[i+(m-1)*length(chip_fulllist),"se"]<-summary(f)$coefficient[1,3]
			h[i+(m-1)*length(chip_fulllist),"z"]<-summary(f)$coefficient[1,4]
      			h[i+(m-1)*length(chip_fulllist),"pval"]<-summary(f)$coefficient[1,5]
      			h[i+(m-1)*length(chip_fulllist),"lci"]<-exp(summary(f)$coefficient[1,1]-1.96*summary(f)$coefficient[1,3])
      			h[i+(m-1)*length(chip_fulllist),"uci"]<-exp(summary(f)$coefficient[1,1]+1.96*summary(f)$coefficient[1,3])
		} 	
	}
}
write.table(h, "/medpop/esp2/zyu/others/checkforFuster/inci_ukb_chip_CVDs_v5_2023_removeleukemia.txt", sep="\t", quote=FALSE, row.names=FALSE)	

##################################interaction#############################
#by CHIP
h_chip = data.frame(chip=character(), outcome=character(), n_nochip=double(), hr_nochip=double(), se_nochip=double(), z_nochip=double(), pval_nochip=double(), lci_nochip=double(), uci_nochip=double(), 
               n_chip=double(), hr_chip=double(), se_chip=double(), z_chip=double(), pval_chip=double(),lci_chip=double(), uci_chip=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_chip[index,"chip"] <- chip
	h_chip[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", "colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==0& pheno_chip$Prev_AnyHemeCa==0),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==1& pheno_chip$Prev_AnyHemeCa==0),])
        fmla_inter = as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_AnyHemeCa==0),])
        h_chip[index,"n_nochip"] <- summary(f1)$n
        h_chip[index,"hr_nochip"] <- summary(f1)$coefficient[1,2]
        h_chip[index,"se_nochip"] <- summary(f1)$coefficient[1,3]
        h_chip[index,"z_nochip"] <- summary(f1)$coefficient[1,4]
        h_chip[index,"pval_nochip"] <- summary(f1)$coefficient[1,5]
        h_chip[index,"lci_nochip"] <- exp(confint(f1)[1,1])
        h_chip[index,"uci_nochip"] <- exp(confint(f1)[1,2])
        h_chip[index,"n_chip"] <- summary(f2)$n
        h_chip[index,"hr_chip"] <- summary(f2)$coefficient[1,2]
        h_chip[index,"se_chip"] <- summary(f2)$coefficient[1,3]
        h_chip[index,"z_chip"] <- summary(f2)$coefficient[1,4]
        h_chip[index,"pval_chip"] <- summary(f2)$coefficient[1,5]
        h_chip[index,"lci_chip"] <- exp(confint(f2)[1,1])
        h_chip[index,"uci_chip"] <- exp(confint(f2)[1,2])
        h_chip[index,"beta_inter"] <- summary(f3)$coefficient[18,2]
        h_chip[index,"se_inter"] <- summary(f3)$coefficient[18,3]
        h_chip[index,"z_inter"] <- summary(f3)$coefficient[18,4]
        h_chip[index,"pval_inter"] <- summary(f3)$coefficient[18,5]
        index <- index + 1
    }
}

write.table(h_chip, "/medpop/esp2/zyu/others/checkforFuster/inci_ukb_chipcolchicine_stratify_withmorecolchicine_v5_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	

#by colchicine

outcome_list=("MI_2023")

h_colchicine = data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), hr_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), hr_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_colchicine[index,"chip"] <- chip
	      h_colchicine[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==0& pheno_chip$Prev_AnyHemeCa==0),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==1& pheno_chip$Prev_AnyHemeCa==0),])
        fmla_inter = as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_AnyHemeCa==0),])
        h_colchicine[index,"n_nocolchicine"] <- summary(f1)$n
        h_colchicine[index,"hr_nocolchicine"] <- summary(f1)$coefficient[1,2]
        h_colchicine[index,"se_nocolchicine"] <- summary(f1)$coefficient[1,3]
        h_colchicine[index,"z_nocolchicine"] <- summary(f1)$coefficient[1,4]
        h_colchicine[index,"pval_nocolchicine"] <- summary(f1)$coefficient[1,5]
        h_colchicine[index,"lci_nocolchicine"] <- exp(confint(f1)[1,1])
        h_colchicine[index,"uci_nocolchicine"] <- exp(confint(f1)[1,2])
        h_colchicine[index,"n_colchicine"] <- summary(f2)$n
        h_colchicine[index,"hr_colchicine"] <- summary(f2)$coefficient[1,2]
        h_colchicine[index,"se_colchicine"] <- summary(f2)$coefficient[1,3]
        h_colchicine[index,"z_colchicine"] <- summary(f2)$coefficient[1,4]
        h_colchicine[index,"pval_colchicine"] <- summary(f2)$coefficient[1,5]
        h_colchicine[index,"lci_colchicine"] <- exp(confint(f2)[1,1])
        h_colchicine[index,"uci_colchicine"] <- exp(confint(f2)[1,2])
        h_colchicine[index,"beta_inter"] <- summary(f3)$coefficient[18,2]
        h_colchicine[index,"se_inter"] <- summary(f3)$coefficient[18,3]
        h_colchicine[index,"z_inter"] <- summary(f3)$coefficient[18,4]
        h_colchicine[index,"pval_inter"] <- summary(f3)$coefficient[18,5]
        index <- index + 1
  }
}

write.table(h_colchicine, "/medpop/esp2/zyu/others/tet2colchicine/output/inci_ukb_chipcolchicine_stratifycolchicine_withmorecolchicine_v5_2023.txt", sep="\t", quote=FALSE, row.names=FALSE)	


#####removing leukemia
#by CHIP
h_chip = data.frame(chip=character(), outcome=character(), n_nochip=double(), hr_nochip=double(), se_nochip=double(), z_nochip=double(), pval_nochip=double(), lci_nochip=double(), uci_nochip=double(), 
               n_chip=double(), hr_chip=double(), se_chip=double(), z_chip=double(), pval_chip=double(),lci_chip=double(), uci_chip=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_chip[index,"chip"] <- chip
	h_chip[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", "colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==0& pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip[[chip]]==1& pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        fmla_inter = as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        h_chip[index,"n_nochip"] <- summary(f1)$n
        h_chip[index,"hr_nochip"] <- summary(f1)$coefficient[1,2]
        h_chip[index,"se_nochip"] <- summary(f1)$coefficient[1,3]
        h_chip[index,"z_nochip"] <- summary(f1)$coefficient[1,4]
        h_chip[index,"pval_nochip"] <- summary(f1)$coefficient[1,5]
        h_chip[index,"lci_nochip"] <- exp(confint(f1)[1,1])
        h_chip[index,"uci_nochip"] <- exp(confint(f1)[1,2])
        h_chip[index,"n_chip"] <- summary(f2)$n
        h_chip[index,"hr_chip"] <- summary(f2)$coefficient[1,2]
        h_chip[index,"se_chip"] <- summary(f2)$coefficient[1,3]
        h_chip[index,"z_chip"] <- summary(f2)$coefficient[1,4]
        h_chip[index,"pval_chip"] <- summary(f2)$coefficient[1,5]
        h_chip[index,"lci_chip"] <- exp(confint(f2)[1,1])
        h_chip[index,"uci_chip"] <- exp(confint(f2)[1,2])
        h_chip[index,"beta_inter"] <- summary(f3)$coefficient[18,2]
        h_chip[index,"se_inter"] <- summary(f3)$coefficient[18,3]
        h_chip[index,"z_inter"] <- summary(f3)$coefficient[18,4]
        h_chip[index,"pval_inter"] <- summary(f3)$coefficient[18,5]
        index <- index + 1
    }
}

write.table(h_chip, "/medpop/esp2/zyu/others/tet2colchicine/output/inci_ukb_chipcolchicine_stratifychip_withmorecolchicine_v5_2023_removeleukemia.txt", sep="\t", quote=FALSE, row.names=FALSE)	


#by colchicine
h_colchicine = data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), hr_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), hr_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

index = 1

for(chip in chip_fulllist) {
  for(outcome in outcome_list) {
        h_colchicine[index,"chip"] <- chip
	      h_colchicine[index, "outcome"] <- outcome
        fmla <- as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f1 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==0& pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        f2 <- coxph(fmla,data=pheno_chip[which(pheno_chip$colchicine==1& pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        fmla_inter = as.formula(paste("Surv(as.numeric(FollowUp_", outcome, "),Incd_", outcome, ") ~ ", chip, "*colchicine+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+Prev_Diabetes_Type_2+Prev_Gout",sep=""))
        f3 <- coxph(fmla_inter,data=pheno_chip[which(pheno_chip$Prev_Myeloid_leukemia==0 & pheno_chip$Prev_Lymphoid_leukemia==0),])
        h_colchicine[index,"n_nocolchicine"] <- summary(f1)$n
        h_colchicine[index,"hr_nocolchicine"] <- summary(f1)$coefficient[1,2]
        h_colchicine[index,"se_nocolchicine"] <- summary(f1)$coefficient[1,3]
        h_colchicine[index,"z_nocolchicine"] <- summary(f1)$coefficient[1,4]
        h_colchicine[index,"pval_nocolchicine"] <- summary(f1)$coefficient[1,5]
        h_colchicine[index,"lci_nocolchicine"] <- exp(confint(f1)[1,1])
        h_colchicine[index,"uci_nocolchicine"] <- exp(confint(f1)[1,2])
        h_colchicine[index,"n_colchicine"] <- summary(f2)$n
        h_colchicine[index,"hr_colchicine"] <- summary(f2)$coefficient[1,2]
        h_colchicine[index,"se_colchicine"] <- summary(f2)$coefficient[1,3]
        h_colchicine[index,"z_colchicine"] <- summary(f2)$coefficient[1,4]
        h_colchicine[index,"pval_colchicine"] <- summary(f2)$coefficient[1,5]
        h_colchicine[index,"lci_colchicine"] <- exp(confint(f2)[1,1])
        h_colchicine[index,"uci_colchicine"] <- exp(confint(f2)[1,2])
        h_colchicine[index,"beta_inter"] <- summary(f3)$coefficient[18,2]
        h_colchicine[index,"se_inter"] <- summary(f3)$coefficient[18,3]
        h_colchicine[index,"z_inter"] <- summary(f3)$coefficient[18,4]
        h_colchicine[index,"pval_inter"] <- summary(f3)$coefficient[18,5]
        index <- index + 1
  }
}

write.table(h_colchicine, "/medpop/esp2/zyu/others/tet2colchicine/output/inci_ukb_chipcolchicine_stratifycolchicine_withmorecolchicine_v5_2023_removeleukemia.txt", sep="\t", quote=FALSE, row.names=FALSE)	























