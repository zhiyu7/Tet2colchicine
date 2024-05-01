gc()
rm(list=ls())

library(dplyr)
library(tidyverse)
library(data.table)
library(ukbtools) 
library(survival)
library(survminer)
library(ggrepel)
library(gridExtra)
library(plinkQC)
library(tableone)
####################################################################
#####################MGBB###########################################
####################################################################
##########################read in data############################
colchicine_mgbb=fread("/medpop/esp2/btruong/Projects/tentative/colchicine/data/colchicine_cvd_chip_MGBB_221017.txt")
colchicine_mgbb_2=fread("/medpop/esp2/zyu/others/tet2colchicine/cvd_chip_230909.txt")
colchicine_mgbb_2=colchicine_mgbb_2%>%select(Biobank_Subject_ID, Deceased, age_death, Has_hf, Incd_hf, age_Incd_hf, Has_CAD_MI, Has_CAD_MI_allcauseMortality,
age_CAD_MI, age_CAD_MI_allcauseMortality, Inc_CAD_MI, Inc_CAD_MI_allcauseMortality, Has_CAD_MI_IsStroke_allcauseMortality, age_CAD_MI_IsStroke_allcauseMortality,
Inc_CAD_MI_IsStroke_allcauseMortality)

colchicine_mgbb=merge(colchicine_mgbb, colchicine_mgbb_2, by="Biobank_Subject_ID", all.x=T)

#CHIP data
colchicine_mgbb$hasTET2[is.na(colchicine_mgbb$hasTET2)] <- 0
colchicine_mgbb$hasDNMT3A[is.na(colchicine_mgbb$hasDNMT3A)] <- 0
colchicine_mgbb$hasASXL1[is.na(colchicine_mgbb$hasASXL1)] <- 0
colchicine_mgbb$hasJAK2[is.na(colchicine_mgbb$hasJAK2)] <- 0

mgbb_chip_vaf=fread("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/mgbb53k.chip_var_N5910_nchip7145.csv.gz")
colnames(mgbb_chip_vaf)[1]="Biobank_Subject_ID"

mgbb_chip_vaf_tet2=mgbb_chip_vaf %>% filter(Gene=="TET2") %>%select(Biobank_Subject_ID, VAF) %>% group_by(Biobank_Subject_ID) %>%summarise(TET2_VAF=sum(VAF))
colchicine_mgbb=merge(colchicine_mgbb, mgbb_chip_vaf_tet2, by="Biobank_Subject_ID", all.x=T)
colchicine_mgbb=colchicine_mgbb%>%mutate(TET2_VAF=ifelse(hasTET2==0, 0, ifelse(hasTET2==1, TET2_VAF, NA)))
mgbb_chip_vaf_dnmt3a=mgbb_chip_vaf %>% filter(Gene=="DNMT3A") %>%select(Biobank_Subject_ID, VAF) %>% group_by(Biobank_Subject_ID) %>%summarise(DNMT3A_VAF=sum(VAF))
colchicine_mgbb=merge(colchicine_mgbb, mgbb_chip_vaf_dnmt3a, by="Biobank_Subject_ID", all.x=T)
colchicine_mgbb=colchicine_mgbb%>%mutate(DNMT3A_VAF=ifelse(hasDNMT3A==0, 0, ifelse(hasDNMT3A==1, DNMT3A_VAF, NA)))
mgbb_chip_vaf_asxl1=mgbb_chip_vaf %>% filter(Gene=="ASXL1") %>%select(Biobank_Subject_ID, VAF) %>% group_by(Biobank_Subject_ID) %>%summarise(ASXL1_VAF=sum(VAF))
colchicine_mgbb=merge(colchicine_mgbb, mgbb_chip_vaf_asxl1, by="Biobank_Subject_ID", all.x=T)
colchicine_mgbb=colchicine_mgbb%>%mutate(ASXL1_VAF=ifelse(hasASXL1==0, 0, ifelse(hasASXL1==1, ASXL1_VAF, NA)))
mgbb_chip_vaf_jak2=mgbb_chip_vaf %>% filter(Gene=="JAK2") %>%select(Biobank_Subject_ID, VAF) %>% group_by(Biobank_Subject_ID) %>%summarise(JAK2_VAF=sum(VAF))
colchicine_mgbb=merge(colchicine_mgbb, mgbb_chip_vaf_jak2, by="Biobank_Subject_ID", all.x=T)
colchicine_mgbb=colchicine_mgbb%>%mutate(JAK2_VAF=ifelse(hasJAK2==0, 0, ifelse(hasJAK2==1, JAK2_VAF, NA)))

#biomarker
colchicine_mgbb_biomarker=fread("/medpop/esp2/btruong/Projects/tentative/colchicine/data/colchicine_curated_221025.txt")
colchicine_mgbb_biomarker = as.data.frame(apply(colchicine_mgbb_biomarker, 2, as.numeric))
colnames(colchicine_mgbb_biomarker)[1]="Biobank_Subject_ID"

colchicine_mgbb=merge(colchicine_mgbb, colchicine_mgbb_biomarker, by="Biobank_Subject_ID")

#additional biomarker required in revision
colchicine_mgbb_biomarker_2=fread("/medpop/esp2/btruong/Projects/zhi/general/zyu_28102023.txt")
colchicine_mgbb_biomarker_2 = as.data.frame(apply(colchicine_mgbb_biomarker_2, 2, as.numeric))
colnames(colchicine_mgbb_biomarker_2)[1]="Biobank_Subject_ID"

colchicine_mgbb=merge(colchicine_mgbb, colchicine_mgbb_biomarker_2, by="Biobank_Subject_ID")

#further processing
colchicine_mgbb$in_white_British_ancestry_subset=ifelse(colchicine_mgbb$ancestry_pred=="EUR", 1, 0)
colchicine_mgbb$bmi=ifelse(colchicine_mgbb$bmi>100, NA, colchicine_mgbb$bmi)

#Satoshi curated relatedness
related_remove=fread("/medpop/esp2/projects/MGB_Biobank/genotype/53K_GSA/release/GSA_53K.related.tsv", header=F)
related_remove=related_remove%>%separate(V1, c(NA, NA, "IID"))

colchicine_mgbb=colchicine_mgbb[!colchicine_mgbb$Biobank_Subject_ID %in% related_remove$IID,]

#add new gout
gout=fread("/medpop/esp2/btruong/Projects/zhi/gout/gout.txt")
colchicine_mgbb$Biobank_Subject_ID=as.character(colchicine_mgbb$Biobank_Subject_ID)
#gout$eid=as.integer(gout$eid)
colchicine_mgbb=merge(colchicine_mgbb, gout, by.x="Biobank_Subject_ID",by.y="eid")

#add NLP MI
mi=fread("/medpop/esp2/zyu/others/tet2colchicine/MI_fullMGBB_BiobankPortal_mch24_2022-10-14-175126.csv")
mi=mi[,c(1,7)]
colnames(mi)=c("Biobank_Subject_ID", "Has_NLP_MI")
mi$Biobank_Subject_ID=as.character(mi$Biobank_Subject_ID)
colchicine_mgbb=merge(colchicine_mgbb, mi, by="Biobank_Subject_ID",all.x=T)
colchicine_mgbb$Has_NLP_MI[colchicine_mgbb$Has_NLP_MI=="Yes"]=1
colchicine_mgbb$Has_NLP_MI[is.na(colchicine_mgbb$Has_NLP_MI)]=0
colchicine_mgbb$Has_NLP_MI=as.numeric(colchicine_mgbb$Has_NLP_MI)

#add T2D
t2d=fread("/medpop/esp2/btruong/Projects/zhi/t2d/t2d.txt")
#t2d$eid=as.integer(t2d$eid)
colchicine_mgbb=merge(colchicine_mgbb, t2d, by.x="Biobank_Subject_ID",by.y="eid", all.x=T)

#add cancer
cancer=fread("/medpop/esp2/btruong/Projects/zhi/anycancer/any_cancer_raw.txt")
colchicine_mgbb=merge(colchicine_mgbb, cancer, by.x="Biobank_Subject_ID", by.y="eid")
leukemia=fread("/medpop/esp2/btruong/Projects/zhi/anycancer/leukemia.txt")
colchicine_mgbb=merge(colchicine_mgbb, leukemia, by.x="Biobank_Subject_ID", by.y="eid")

###reformat Buu's dataset
colchicine_mgbb$Time_Censor_CAD_MI=colchicine_mgbb$age_CAD_MI-colchicine_mgbb$age_enroll
colchicine_mgbb$Time_Censor_death=colchicine_mgbb$age_death-colchicine_mgbb$age_enroll
colchicine_mgbb$Time_Censor_CAD_MI_allcauseMortality=colchicine_mgbb$age_CAD_MI_allcauseMortality-colchicine_mgbb$age_enroll
colchicine_mgbb$Time_Censor_hf=colchicine_mgbb$age_Incd_hf-colchicine_mgbb$age_enroll
colchicine_mgbb$Time_Censor_CAD_MI_IsStroke_allcauseMortality=colchicine_mgbb$age_CAD_MI_IsStroke_allcauseMortality-colchicine_mgbb$age_enroll
colchicine_mgbb$Incd_death=colchicine_mgbb$Deceased
colchicine_mgbb=colchicine_mgbb%>%rename(Incd_CAD_MI=Inc_CAD_MI, Incd_CAD_MI_allcauseMortality=Inc_CAD_MI_allcauseMortality, 
				              Incd_CAD_MI_IsStroke_allcauseMortality=Inc_CAD_MI_IsStroke_allcauseMortality, Has_death=Deceased)

#list of event types
outcome=c("cvd", "cad", "ischemic_stroke", "mi", "CAD_MI", "hf", "CAD_MI_allcauseMortality", "CAD_MI_IsStroke_allcauseMortality", "death")

#replace non-positive values with NA in columns that start with "Time_Censor_"
colchicine_mgbb <- colchicine_mgbb %>% mutate(across(starts_with("Time_Censor_"), ~ ifelse(. <= 0, NA, .)))

#replace prevalent cases follow up time to NA
for (outcome_use in outcome) {
  	has_col <- paste0("Has_", outcome_use)
  	incd_col <- paste0("Incd_", outcome_use)
  	time_censor_col <- paste0("Time_Censor_", outcome_use)
  
  	colchicine_mgbb <- colchicine_mgbb %>% mutate(!!time_censor_col := ifelse((!!as.name(has_col) == 1) & (!!as.name(incd_col) == 0), NA, !!as.name(time_censor_col)))
}

#add medication
ramed=fread("/medpop/esp2/btruong/Projects/zhi/RA/data/zhi_RA_220620.csv")
colnames(ramed)[1]="IID"
colnames(ramed)[7]="Polymyalgia_rheumatica"
colnames(ramed)[8]="Giant_cell_arteritis"
ramed$steroid_RA=ifelse(ramed$Prednisolone == 1 | ramed$Prednisone == 1 | ramed$Dexamethasone == 1 | ramed$Methylprednisolone ==1, 1, 0)
ramed$csDMARD_RA=ifelse(ramed$Hydroxychloroquine == 1 | ramed$Sulfasalazine == 1 | ramed$Methotrexate == 1 | ramed$Leflunomide ==1, 1, 0)
ramed$bDMARD_RA=ifelse(ramed$Etanercept == 1 | ramed$Adalimumab == 1 | ramed$Infliximab == 1 | ramed$Golimumab ==1 | ramed$Certolizumab==1 | ramed$Certolizumab_pegol_syringe==1 | ramed$Abatacept==1 |  ramed$Tocilizumab==1 | ramed$Rituximab==1, 1, 0)
ramed$tsDMARD_RA=ifelse(ramed$Tofacitinib ==1, 1, 0)
ramed$anyRAmed=ifelse(ramed$csDMARD_RA==1 | ramed$steroid_RA==1 | ramed$bDMARD_RA==1 | ramed$tsDMARD_RA==1, 1, 0)
ramed$anyDMARD=ifelse(ramed$csDMARD_RA==1 | ramed$bDMARD_RA==1 | ramed$tsDMARD_RA==1, 1, 0)
ramed$antiCCP[is.na(ramed$antiCCP)] <- 0
ramed$antiCCP_pos=ifelse(ramed$antiCCP>0, 1, 0)
ramed=ramed%>%select(Biobank_Subject_ID=IID, steroid_RA, csDMARD_RA, bDMARD_RA, tsDMARD_RA, anyRAmed, anyDMARD)
colchicine_mgbb=merge(colchicine_mgbb, ramed, by="Biobank_Subject_ID", all.x=T)

aspirin=fread("/medpop/esp2/zyu/others/tet2colchicine/data/aspirin_biobankid.txt")
names(aspirin)[1]="Biobank_Subject_ID"
insulin=fread("/medpop/esp2/zyu/others/tet2colchicine/data/insulin_biobankid.txt")
names(insulin)[1]="Biobank_Subject_ID"

aspirin$Biobank_Subject_ID=as.character(aspirin$Biobank_Subject_ID)
insulin$Biobank_Subject_ID=as.character(insulin$Biobank_Subject_ID)

colchicine_mgbb=merge(colchicine_mgbb, aspirin, by="Biobank_Subject_ID", all.x=T)
colchicine_mgbb=merge(colchicine_mgbb, insulin, by="Biobank_Subject_ID", all.x=T)

######please note here#####
#in revision we subset to >40 in MGBB to have consisten age range as UKB#
#subset by age#
colchicine_mgbb <- colchicine_mgbb%>%filter(Age_Genotyping>40)

###output study id for merging with enrollment data on laptop
id=colchicine_mgbb[which(colchicine_mgbb$leu==0),]$Biobank_Subject_ID
write.table(id,"/medpop/esp2/zyu/others/tet2colchicine/output/id.txt", sep="\t", quote=FALSE, row.names=TRUE)

########################Table 1################################
factorVars <- c("sex","current_smoker", "gout", "t2d", "in_white_British_ancestry_subset",  "hasCHIP", "hasDNMT3A", "hasTET2", "hasASXL1","hasJAK2",
"Has_cvd", "Has_cad", "Has_mi", "Has_ischemic_stroke", "Has_NLP_MI", "any_cancer", "any_cancer_RPDR", "leu", "HTN", "FH_CVD", "CKD", "renal_dialysis", 
"anyDMARD", "aspirin", "insulin" )

vars <- c("Age_Genotyping", "sex","current_smoker", "gout", "t2d", "in_white_British_ancestry_subset", 
"HCT", "MCV", "RDW", "MCHC", "PLT", "RBC", "WBC", "HGB", "Cholesterol", "Triglycerides", "HDL", "LDL",
"hasCHIP", "hasDNMT3A", "hasTET2", "hasASXL1","hasJAK2","Has_cvd", "Has_cad", "Has_mi", "Has_ischemic_stroke",
"Has_NLP_MI", "any_cancer", "any_cancer_RPDR", "leu", "HTN", "FH_CVD", "CKD", "renal_dialysis", "CRP", "Creatinine", 
"bmi", "eGFR", "anyDMARD","aspirin", "insulin","BMI")

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0),], factorVars = factorVars)
write.table(print(tableOne, contDigits = 1, pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_nostratification_older40_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_v2witht2d_older40_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0 & colchicine_mgbb$hasCHIP==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_v2witht2d_older40_CHIP_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0 & colchicine_mgbb$hasTET2==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_v2witht2d_older40_TET2_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0 & colchicine_mgbb$hasDNMT3A==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_v2witht2d_older40_DNMT3A_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$leu==0 & colchicine_mgbb$hasASXL1==1),], strata="colchicine", factorVars = factorVars)
write.table(print(tableOne, contDigits = 1,  pDigits=2),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_MGBB_withgout_v2witht2d_older40_ASXL1_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)

vars <- c("TET2_VAF")
tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$hasTET2==1 & colchicine_mgbb$leu==0),], strata="colchicine")
write.table(print(tableOne),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_tet2vaf_withgout_older40_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)
print(tableOne,nonnormal = vars)

vars <- c("DNMT3A_VAF")
tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$hasDNMT3A==1 & colchicine_mgbb$leu==0),], strata="colchicine")
write.table(print(tableOne),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_dnmt3avaf_withgout_older40_v1.txt", sep="\t", quote=FALSE, row.names=TRUE)
print(tableOne,nonnormal = vars)

vars <- c("ASXL1_VAF")
tableOne <- CreateTableOne(vars = vars, data = colchicine_mgbb[which(colchicine_mgbb$hasASXL1==1 & colchicine_mgbb$leu==0),], strata="colchicine")
write.table(print(tableOne),"/medpop/esp2/zyu/others/tet2colchicine/output/table1_asxl1vaf_withgout.txt_older40_v1", sep="\t", quote=FALSE, row.names=TRUE)
print(tableOne,nonnormal = vars)

#############################################################################
########################Cross-sectional analyis###############################
###########################################################################
####check
chip=c("hasCHIP", "hasTET2", "hasDNMT3A", "hasASXL1", "colchicine")#colchicine is here is just for addressing reviewer's comments
outcome=c("cvd", "cad", "ischemic_stroke", "mi", "CAD_MI", "hf", "CAD_MI_allcauseMortality", "CAD_MI_IsStroke_allcauseMortality", "death")

##########this is final analysis###############
#with removing leukemia#
h_test<-data.frame(outcome=character(), chip_exposure=character(), n=double(), or=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      
for(i in 1:length(outcome)) {
	for(j in 1:length(chip)) {
      		h_test[j+(i-1)*length(chip),"outcome"]<-outcome[i]
		h_test[j+(i-1)*length(chip),"chip_exposure"]<-chip[j]
     		fmla<-as.formula(paste("Has_", outcome[i], "~", chip[j], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
      		f<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$leu==0),], family = "binomial")
		h_test[j+(i-1)*length(chip),"n"]<-summary(f)$df.null+1
     	 	h_test[j+(i-1)*length(chip),"or"]<-exp(summary(f)$coefficient[2,1])
   	   	h_test[j+(i-1)*length(chip),"pval"]<-summary(f)$coefficient[2,4]
   	   	h_test[j+(i-1)*length(chip),"lci"]<-exp(summary(f)$coefficient[2,1]-1.96*summary(f)$coefficient[2,2])
      		h_test[j+(i-1)*length(chip),"uci"]<-exp(summary(f)$coefficient[2,1]+1.96*summary(f)$coefficient[2,2])
	}
}

#write.table(h_test, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_revmoingleukemia_tet2_chip_witht2d_remove18_v4moreoutcomes.txt", sep="\t", quote=FALSE, row.names=FALSE)	
write.table(h_test, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_revmoingleukemia_tet2_chip_witht2d_remove40_v4moreoutcomes_notadjustforcancer.txt", sep="\t", quote=FALSE, row.names=FALSE)		

#sensitivity with removing RPDR cancer
h_test<-data.frame(chip=character(), exposure=character(), n=double(), or=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

h_test<-data.frame(outcome=character(), chip_exposure=character(), n=double(), or=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      
for(i in 1:length(outcome)) {
	for(j in 1:length(chip)) {
      		h_test[j+(i-1)*length(chip),"outcome"]<-outcome[i]
		h_test[j+(i-1)*length(chip),"chip_exposure"]<-chip[j]
     		fmla<-as.formula(paste("Has_", outcome[i], "~", chip[j], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
      		f<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$any_cancer_RPDR==0),], family = "binomial")
		h_test[j+(i-1)*length(chip),"n"]<-summary(f)$df.null+1
     	 	h_test[j+(i-1)*length(chip),"or"]<-exp(summary(f)$coefficient[2,1])
   	   	h_test[j+(i-1)*length(chip),"pval"]<-summary(f)$coefficient[2,4]
   	   	h_test[j+(i-1)*length(chip),"lci"]<-exp(summary(f)$coefficient[2,1]-1.96*summary(f)$coefficient[2,2])
      		h_test[j+(i-1)*length(chip),"uci"]<-exp(summary(f)$coefficient[2,1]+1.96*summary(f)$coefficient[2,2])
	}
}

write.table(h_test, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_revmoingrpdrcancer_tet2_chip_witht2d_v3withmoreoutcomes.txt", sep="\t", quote=FALSE, row.names=FALSE)	

####by colchicine
#with removing leukemia

h_colchicine_noremovinggout<-data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), or_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), or_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

for(i in 1:length(chip)) {
	for(j in 1:length(outcome)) {
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"chip"]<-chip[i]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome), "outcome"]<-outcome[j]
     		fmla<-as.formula(paste("Has_", outcome[j], "~", chip[i], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
      		f1<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==0 & colchicine_mgbb$leu==0),], family = "binomial")
      		f2<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==1 & colchicine_mgbb$leu==0),], family = "binomial")
		fmla_inter=as.formula(paste("Has_", outcome[j], "~", chip[i], "*colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
		f3<-glm(fmla_inter,data=colchicine_mgbb[which(colchicine_mgbb$leu==0),], family = "binomial")
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_nocolchicine"]<-summary(f1)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"or_nocolchicine"]<-exp(summary(f1)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_nocolchicine"]<-summary(f1)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_nocolchicine"]<-summary(f1)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_nocolchicine"]<-summary(f1)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]-1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]+1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_colchicine"]<-summary(f2)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"or_colchicine"]<-exp(summary(f2)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_colchicine"]<-summary(f2)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_colchicine"]<-summary(f2)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_colchicine"]<-summary(f2)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_colchicine"]<-exp(summary(f2)$coefficient[2,1]-1.96*summary(f2)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_colchicine"]<-exp(summary(f2)$coefficient[2,1]+1.96*summary(f2)$coefficient[2,2])
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"beta_inter"]<-exp(summary(f3)$coefficient[23,1])
     		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_inter"]<-summary(f3)$coefficient[23,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_inter"]<-summary(f3)$coefficient[23,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_inter"]<-summary(f3)$coefficient[23,4]
	}
}
#write.table(h_colchicine_noremovinggout, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_excludingleukemia_tet2_chipcolchicine_noremovinggout_stratify_witht2d_remove18_v4moreoutcomes.txt", sep="\t", quote=FALSE, row.names=FALSE)	
write.table(h_colchicine_noremovinggout, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_excludingleukemia_tet2_chipcolchicine_noremovinggout_stratify_witht2d_remove40_v4moreoutcomes_notadjustforcancer.txt", sep="\t", quote=FALSE, row.names=FALSE)	

#sensitivity adjusting for DMARD
h_colchicine_noremovinggout<-data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), or_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), or_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

for(i in 1:length(chip)) {
	for(j in 1:length(outcome)) {
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"chip"]<-chip[i]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome), "outcome"]<-outcome[j]
     		fmla<-as.formula(paste("Has_", outcome[j], "~", chip[i], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d+anyDMARD",sep=""))
      		f1<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==0 & colchicine_mgbb$leu==0),], family = "binomial")
      		f2<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==1 & colchicine_mgbb$leu==0),], family = "binomial")
		fmla_inter=as.formula(paste("Has_", outcome[j], "~", chip[i], "*colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d+anyDMARD",sep=""))
		f3<-glm(fmla_inter,data=colchicine_mgbb[which(colchicine_mgbb$leu==0),], family = "binomial")
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_nocolchicine"]<-summary(f1)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"or_nocolchicine"]<-exp(summary(f1)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_nocolchicine"]<-summary(f1)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_nocolchicine"]<-summary(f1)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_nocolchicine"]<-summary(f1)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]-1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]+1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_colchicine"]<-summary(f2)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"or_colchicine"]<-exp(summary(f2)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_colchicine"]<-summary(f2)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_colchicine"]<-summary(f2)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_colchicine"]<-summary(f2)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_colchicine"]<-exp(summary(f2)$coefficient[2,1]-1.96*summary(f2)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_colchicine"]<-exp(summary(f2)$coefficient[2,1]+1.96*summary(f2)$coefficient[2,2])
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"beta_inter"]<-exp(summary(f3)$coefficient[24,1])
     		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_inter"]<-summary(f3)$coefficient[24,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_inter"]<-summary(f3)$coefficient[24,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_inter"]<-summary(f3)$coefficient[24,4]
	}
}
write.table(h_colchicine_noremovinggout, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_excludingleukemia_tet2_chipcolchicine_noremovinggout_stratify_witht2d_remove40_v4moreoutcomes_notadjustforcancer_adjustingDMARD.txt", sep="\t", quote=FALSE, row.names=FALSE)	

#sensitivity with removing RPDR cancer

h_colchicine_noremovinggout<-data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), hr_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), hr_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

for(i in 1:length(chip)) {
	for(j in 1:length(outcome)) {
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"chip"]<-chip[i]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome), "outcome"]<-outcome[j]
     		fmla<-as.formula(paste("Has_", outcome[j], "~", chip[i], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
      		f1<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==0 & colchicine_mgbb$any_cancer_RPDR==0),], family = "binomial")
      		f2<-glm(fmla,data=colchicine_mgbb[which(colchicine_mgbb$colchicine==1 & colchicine_mgbb$any_cancer_RPDR==0),], family = "binomial")
		fmla_inter=as.formula(paste("Has_", outcome[j], "~", chip[i], "*colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
		f3<-glm(fmla_inter,data=colchicine_mgbb[which(colchicine_mgbb$any_cancer_RPDR==0),], family = "binomial")
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_nocolchicine"]<-summary(f1)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"hr_nocolchicine"]<-exp(summary(f1)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_nocolchicine"]<-summary(f1)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_nocolchicine"]<-summary(f1)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_nocolchicine"]<-summary(f1)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]-1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_nocolchicine"]<-exp(summary(f1)$coefficient[2,1]+1.96*summary(f1)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"n_colchicine"]<-summary(f2)$df.null+1
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"hr_colchicine"]<-exp(summary(f2)$coefficient[2,1])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_colchicine"]<-summary(f2)$coefficient[2,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_colchicine"]<-summary(f2)$coefficient[2,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_colchicine"]<-summary(f2)$coefficient[2,4]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"lci_colchicine"]<-exp(summary(f2)$coefficient[2,1]-1.96*summary(f2)$coefficient[2,2])
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"uci_colchicine"]<-exp(summary(f2)$coefficient[2,1]+1.96*summary(f2)$coefficient[2,2])
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"beta_inter"]<-exp(summary(f3)$coefficient[23,1])
     		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"se_inter"]<-summary(f3)$coefficient[23,2]
		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"z_inter"]<-summary(f3)$coefficient[23,3]
      		h_colchicine_noremovinggout[j+(i-1)*length(outcome),"pval_inter"]<-summary(f3)$coefficient[23,4]
	}
}
write.table(h_colchicine_noremovinggout, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_excludingcancer_tet2_chipcolchicine_noremovinggout_stratify_v2witht2d.txt", sep="\t", quote=FALSE, row.names=FALSE)	


###########biomarkers
chip=c("hasCHIP", "hasTET2", "hasDNMT3A", "hasASXL1")
outcome_list=c("HCT", "MCV","RDW","MCHC","PLT","RBC", "WBC", "HGB", "Cholesterol", "Triglycerides", "HDL", "LDL")

temp=colchicine_mgbb 
temp$HCT=scale(log(temp$HCT))
temp$MCV=scale(log(temp$MCV))
temp$RDW=scale(log(temp$RDW))
temp$MCHC=scale(log(temp$MCHC))
temp$PLT=scale(log(temp$PLT))
temp$RBC=scale(log(temp$RBC))
temp$WBC=scale(log(temp$WBC))
temp$HGB=scale(log(temp$HGB))
temp$Cholesterol=scale(log(temp$Cholesterol))
temp$Triglycerides=scale(log(temp$Triglycerides))
temp$HDL=scale(log(temp$HDL))
temp$LDL=scale(log(temp$LDL))

#with removing leukemia
h_colchicine_noremovinggout_markers<-data.frame(chip=character(), outcome=character(), n_nocolchicine=double(), hr_nocolchicine=double(), se_nocolchicine=double(), z_nocolchicine=double(), pval_nocolchicine=double(), lci_nocolchicine=double(), uci_nocolchicine=double(), 
               n_colchicine=double(), hr_colchicine=double(), se_colchicine=double(), z_colchicine=double(), pval_colchicine=double(),lci_colchicine=double(), uci_colchicine=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

for(j in 1:length(outcome_list)) {
	for(i in 1:length(chip)) {
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"chip"]<-chip[i]
		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list), "outcome"]<-outcome_list[j]
		fmla=as.formula(paste(outcome_list[j]," ~ ", chip[i], "+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
      		f1<-glm(fmla,data=temp[which(temp$colchicine==0 & temp$leu==0),], family = "gaussian")
      		f2<-glm(fmla,data=temp[which(temp$colchicine==1 & temp$leu==0),], family = "gaussian")
		fmla_inter=as.formula(paste(outcome_list[j]," ~ ", chip[i], "*colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d",sep=""))
		f3<-glm(fmla_inter,data=temp[which(temp$leu==0),], family = "gaussian")
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"n_nocolchicine"]<-summary(f1)$df.null+1
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"hr_nocolchicine"]<-summary(f1)$coefficient[2,1]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"se_nocolchicine"]<-summary(f1)$coefficient[2,2]
		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"z_nocolchicine"]<-summary(f1)$coefficient[2,3]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"pval_nocolchicine"]<-summary(f1)$coefficient[2,4]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"lci_nocolchicine"]<-summary(f1)$coefficient[2,1]-1.96*summary(f1)$coefficient[2,2]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"uci_nocolchicine"]<-summary(f1)$coefficient[2,1]+1.96*summary(f1)$coefficient[2,2]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"n_colchicine"]<-summary(f2)$df.null+1
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"hr_colchicine"]<-summary(f2)$coefficient[2,1]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"se_colchicine"]<-summary(f2)$coefficient[2,2]
		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"z_colchicine"]<-summary(f2)$coefficient[2,3]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"pval_colchicine"]<-summary(f2)$coefficient[2,4]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"lci_colchicine"]<-summary(f2)$coefficient[2,1]-1.96*summary(f2)$coefficient[2,2]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"uci_colchicine"]<-summary(f2)$coefficient[2,1]+1.96*summary(f2)$coefficient[2,2]
		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"beta_inter"]<-summary(f3)$coefficient[23,1]
     		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"se_inter"]<-summary(f3)$coefficient[23,2]
		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"z_inter"]<-summary(f3)$coefficient[23,3]
      		h_colchicine_noremovinggout_markers[j+(i-1)*length(outcome_list),"pval_inter"]<-summary(f3)$coefficient[23,4]
	}
}
write.table(h_colchicine_noremovinggout_markers, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_tet2_chipcolchicine_removeleukemia_noremovinggout_markers_stratify_v2witht2d_older40_noadjustcancer.txt", sep="\t", quote=FALSE, row.names=FALSE)	


####by CHIP
chip=c("hasCHIP", "hasTET2", "hasDNMT3A", "hasASXL1")
outcome_list=c("cvd", "cad", "ischemic_stroke", "mi", "CAD_MI", "hf", "CAD_MI_allcauseMortality", "CAD_MI_IsStroke_allcauseMortality", "death")

#with removing leukemia

h_chip_noremovinggout<-data.frame(chip=character(), outcome=character(), n_nochip=double(), hr_nochip=double(), se_nochip=double(), z_nochip=double(), pval_nochip=double(), lci_nochip=double(), uci_nochip=double(), 
               n_chip=double(), hr_chip=double(), se_chip=double(), z_chip=double(), pval_chip=double(),lci_chip=double(), uci_chip=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

for(i in 1:length(chip)) {
	for(j in 1:length(outcome_list)) {
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"chip"]<-chip[i]
		h_chip_noremovinggout[j+(i-1)*length(outcome_list), "outcome"]<-outcome_list[j]
            	temp=colchicine_mgbb
     		fmla<-as.formula(paste("Has_", outcome_list[j], "~colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d+any_cancer_RPDR",sep=""))
      		f1<-glm(fmla,data=temp[temp[[chip[i]]] == 0 & temp$leu==0], family = "binomial")
      		f2<-glm(fmla,data=temp[temp[[chip[i]]] == 1 & temp$leu==0], family = "binomial")
		fmla_inter=as.formula(paste("Has_", outcome_list[j], "~", chip[i], "*colchicine+in_white_British_ancestry_subset+Age_Genotyping+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+current_smoker+as.factor(batch)+gout+t2d+any_cancer_RPDR",sep=""))
		f3<-glm(fmla_inter,data=temp[which(temp$leu==0),], family = "binomial")
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"n_nochip"]<-summary(f1)$df.null+1
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"hr_nochip"]<-exp(summary(f1)$coefficient[2,1])
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"se_nochip"]<-summary(f1)$coefficient[2,2]
		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"z_nochip"]<-summary(f1)$coefficient[2,3]
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"pval_nochip"]<-summary(f1)$coefficient[2,4]
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"lci_nochip"]<-exp(summary(f1)$coefficient[2,1]-1.96*summary(f1)$coefficient[2,2])
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"uci_nochip"]<-exp(summary(f1)$coefficient[2,1]+1.96*summary(f1)$coefficient[2,2])
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"n_chip"]<-summary(f2)$df.null+1
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"hr_chip"]<-exp(summary(f2)$coefficient[2,1])
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"se_chip"]<-summary(f2)$coefficient[2,2]
		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"z_chip"]<-summary(f2)$coefficient[2,3]
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"pval_chip"]<-summary(f2)$coefficient[2,4]
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"lci_chip"]<-exp(summary(f2)$coefficient[2,1]-1.96*summary(f2)$coefficient[2,2])
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"uci_chip"]<-exp(summary(f2)$coefficient[2,1]+1.96*summary(f2)$coefficient[2,2])
		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"beta_inter"]<-exp(summary(f3)$coefficient[24,1])
     		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"se_inter"]<-summary(f3)$coefficient[24,2]
		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"z_inter"]<-summary(f3)$coefficient[24,3]
      		h_chip_noremovinggout[j+(i-1)*length(outcome_list),"pval_inter"]<-summary(f3)$coefficient[24,4]
	}
}
write.table(h_chip_noremovinggout, "/medpop/esp2/zyu/others/tet2colchicine/output/crosssectional_mgbb_tet2_chipcolchicine_removeleukemia_noremovinggout_stratifybychip_v2witht2d_older40.txt", sep="\t", quote=FALSE, row.names=FALSE)	


