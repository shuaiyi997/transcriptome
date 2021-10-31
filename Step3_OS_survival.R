#exstract the differential expressed genes (up regulation in tumor)

#Three methods can be used to download the TCGA-LIVE data
##Methods one: TCGAbiolinks, offered by TCGAbiolinks

### install the TCGAbiolinks
library(survival)
library(survminer)
library(biomaRt)
library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
devtools::install_github("BioinformaticsFMRP/TCGAbiolinks") #??????TCGAbiolinks
library("TCGAbiolinks")

?GDCquery 
TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
TCGAbiolinks:::getGDCprojects()$project_id

exp <- GDCquery(  project = "TCGA-LIHC", 
                  legacy = TRUE, 
                  experimental.strategy = "RNA-Seq", 
                  data.category = "Gene expression", 
                  data.type = "Gene expression quantification", 
                  file.type  = "normalized_results",
                  sample.type = c("Primary Tumor"))

###download the expression data
GDCdownload(exp, method = "api", files.per.chunk = 10)

###arrange the information
liver_exp<-GDCprepare(exp)
###import the summarizedExperiment

library(SummarizedExperiment)
count_matrix<-assay(liver_exp)

###extract the patients' names
library(tidyverse)
colnames(count_matrix)<-str_sub(colnames(count_matrix),1,12)

###get the clinical information 
if(FALSE){cli <- GDCquery(  project = "TCGA-LIHC",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  )
GDCdownload(cli)
liver_cli <- GDCprepare(cli)
names(liver_cli)}

clinical <- GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")

### match the expression information and clinical info
####delete the patients without either expresion information and clinical info
#length(intersect(colnames(count_matrix),clinical$submitter_id)) #371 patients 

patients<-intersect(colnames(count_matrix),clinical$submitter_id)

###code provided by TCGAlinks manual  
clinical_patient_Cancer <- GDCquery_clinic("TCGA-LIHC","clinical")

clinical_patient_Cancer<-clinical_patient_Cancer[clinical_patient_Cancer$submitter_id %in% patients,]
dim(clinical_patient_Cancer) #371 patients, 69 information colums

count_matrix<-count_matrix[,colnames(count_matrix) %in% patients] #19947   371

dataLIHCcomplete <- log2(count_matrix)
##head(dataLIHCcomplete)

###gene should intersected among expression profile in TCGA-LIHC,differential expressed genes and secreted protein
gene<-Reduce(intersect, list(rownames(dataLIHCcomplete),rownames(DGE[DGE$change=="UP",]),secretome$Gene)) #remember to change

##VEGFD not in the dataLIHCcomplete
save(DGE_BT_BTN_RPA,DGE_BT_BTN,overlapped,count_matrix,expr,clinical_patient_Cancer,sig,TPM,unfavo_sig,file="./variables.Rdata")
load('./variables.Rdata')
union<-c(gene_DGE,gene_RPA)

#length(BT_BTN_RPA) 42  
#length(BT_BTN_DGE) 34 
dataLIHCcomplete.1<-read.csv("./dataLIHCcomplete.1.csv",sep=",",row.names = 1)

dataLIHCcomplete.1<-dataLIHCcomplete[row.names(dataLIHCcomplete) %in% gene,]

##head(dataLIHCcomplete.1,2)
# TCGA-FV-A3R3 TCGA-DD-A4NV TCGA-ZS-A9CE TCGA-DD-AADP TCGA-CC-5262 TCGA-DD-A4NA TCGA-DD-AAC9 TCGA-BC-A112 TCGA-2Y-A9HA TCGA-BD-A3EP TCGA-DD-AAEI
#COL15A1    11.577855     9.603015     6.368138     9.818107    10.089811     9.254359    10.490446     10.13262     7.505405    10.057185    10.501051
#EFNA4       7.309054     6.995461     7.294641     7.877346     9.001412     7.642365     8.416119 
#data form:normalzed_data

#dim(dataLIHCcomplete.1)
#[1] 371 49,row respents gene symbol and col respents samples

###because of lots of potential genes need to be analysis,so it needs to break down part

tabSurvKM<-TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
                                    dataLIHCcomplete.1,
                                    Genelist = rownames(dataLIHCcomplete.1),
                                    Survresult = T,
                                    ThreshTop=0.5,
                                    ThreshDown=0.5) ##should order by yourself or system
#head(tabSurvKM)
#pvalue Group2 Deaths Group2 Deaths with Top Group2 Deaths with Down Mean Group2 Top Mean Group2 Down Mean Group1
#COL15A1 0.01543882           130                     92                      38        9.362809         6.548728    8.657392
#EFNA4   0.04203990           130                    104                      26        8.357312         6.802964    7.967678
  
tabSurvKMcomplete <- rbind(tabSurvKMcomplete,tabSurvKM)


#View(tabSurvKMcomplete) (p-value<0.01)
#2 significant gene
tabSurvKMcomplete <- tabSurvKMcomplete[order(tabSurvKMcomplete$pvalue, decreasing=F),]
dataLIHCcomplete.1 %>% select(.,rownames(tabSurvKMcomplete))


####the figures produced above are very very very ugly!!!!
dataLIHCcomplete.1<-as.data.frame(t(dataLIHCcomplete.1))
###keep the dataLIHCcomplete.1
write.csv(dataLIHCcomplete.1,"./dataLIHCcomplete.1.csv",quote = F)
dataLIHCcomplete.1<-read.csv("C:/Users/37140/Documents/figures2/dataLIHCcomplete.1.csv",sep=",",row.names = 1)
round(ncol(dataLIHCcomplete.1))
setwd(("C:/Users/37140/Documents/figures2/"))

sig<-c()
for( i in 1:ncol(dataLIHCcomplete.1)){
  message(paste( i, "of ", round(nrow(dataLIHCcomplete.1))))
  name<-names(dataLIHCcomplete.1)[i] 
  print(name)
  sub<-dataLIHCcomplete.1 %>% dplyr::select(.,name)
  sub<-arrange(sub,desc(sub[,name])) #"arrange" is better than "order" 
  #sometimes,for specific gene,  
  high<-mutate(clinical_patient_Cancer[clinical_patient_Cancer$submitter_id %in% rownames(head(sub,round(nrow(sub)*0.5))),],Group=paste0(name,"_High"))
  low<-mutate(clinical_patient_Cancer[clinical_patient_Cancer$submitter_id %in% rownames((tail(sub,ceiling(nrow(sub)*0.5)))),],Group=paste0(name,"_Low"))
  cli<-dplyr::bind_rows(high,low)
  cli$vital_status<-gsub("Alive",0,cli$vital_status)
  cli$vital_status<-gsub("Dead",1,cli$vital_status)
  cli<-cli[cli$vital_status=="0" | cli$vital_status=="1",]
  cli$vital_status<-as.numeric(cli$vital_status)
  fit=survdiff(Surv(cli$days_to_death,as.numeric(cli$vital_status)) ~ as.factor(cli$Group),data=cli)
  #if many factors needs to be considered or the factor is not categorical variables, with different types of variables
  model=coxph(Surv(cli$days_to_death,as.numeric(cli$vital_status)) ~ as.factor(cli$Group),data=cli)
  p_value=pchisq(fit$chisq, length(fit$n)-1, lower.tail = FALSE)
  print(p_value) 
  if(p_value<0.05){
    sur_plot(cli,name)
    sig<-c(sig,name)
  }
}

sur_plot<-function(data,i){
  print(i)
  #print(fit2)
  legend=c(paste0(i,"_High"),paste0(i,"_Low"))
  fit2<- survfit(Surv(cli$days_to_death,as.numeric(cli$vital_status)) ~ as.factor(cli$Group), data = cli)
  #summary(fit2)
  g<-ggsurvplot(fit2, data = data,surv.median.line = "hv",legend.labs=legend,risk.table = F,xlab = "Time",ylab="Overall Survival", pval =TRUE,pval.method=T,plot.title = element_text(size = 1,hjust = 0.5),title=paste0("Overall survival analysis of ",i),risk.table.col="black", risk.table.height = 0.3,font.title=c(18,"bold","black"),font.tickslab=13,conf.int = F) 
  #conf.int=TRUE (show the confident interval or not) break.x.by=100 #set the break of x-axis,legend.labs(set the labels of legend)
  png(paste0("C:/Users/37140/Documents/os_analysis/",i,".png"),width =500,height = 500,res=100)
  print(g)
  dev.off()
  }
#legend.labs = c(paste0(i,"_MUT"),paste0(i,"_WILD")),legend=c("right"),legend.title=""legend = "none"

###survminer plots forest figures by ggforest
ggforest(model,date=NULL,main="",fontsize=0.7,refLabel="",anoDigits="")






