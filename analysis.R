#load the required libraries
library("xlsx")
library("preprocessCore")
library("msImpute")
library("survival")
library("survminer")
library("RegParallel")
library("stringr")
library("org.Hs.eg.db")
library("ReactomePA")
library("pathfindR")

#set directory

setwd("C://Users/kptl632/OneDrive - AZCollaboration/Documents/personal/task_dfkz/")
list.files("Bioinformatic_task/")


#load the input data
prot_abun=read.csv("Bioinformatic_task/Protein_abundance.tsv",
                   sep="\t", 
                   stringsAsFactors = F, 
                   check.names = F,
                   row.names = 1)

metadata=read.xlsx(file = "Bioinformatic_task/sampleAnnotation.xls",
                   sheetIndex = 1)


#basic check

##% of missing data

print(paste0("Percentage missing=",round(100*(sum(is.na(prot_abun))/(nrow(prot_abun)*ncol(prot_abun))),2), " %"))

###dimention
dim(prot_abun)
dim(metadata)

print("Metadtaa has 1 sample extra")

missing_sample=metadata$sample.ID[!metadata$sample.ID %in% colnames(prot_abun)]

metadata[metadata$sample.ID==missing_sample,]

print("most of metadata from missing sample is not available")

metadata_filter=metadata[metadata$sample.ID %in% colnames(prot_abun),]
row.names(metadata_filter)=metadata_filter$sample.ID
metadata_filter$sample.ID=NULL

prot_abun_order=prot_abun[,row.names(metadata_filter)]



#Data processing
qn_pa=normalize.quantiles(as.matrix(prot_abun_order))

row.names(qn_pa)=row.names(prot_abun_order)
colnames(qn_pa)=colnames(prot_abun_order)

log_qn_pa=log2(qn_pa)

max(rowSums( is.na(log_qn_pa)))
table(rowSums( is.na(log_qn_pa)))

#removed features
row.names(log_qn_pa[rowSums( is.na(log_qn_pa))>=46,])


filtered= log_qn_pa[ rowSums( is.na(log_qn_pa)) <= 45, ]

library(msImpute)
imputed=msImpute(filtered,method = "v1")

imputed_mnar=msImpute(filtered,group = as.character(metadata_filter$died))

#lowest value is also quite high suggesting strong signals and no need for MNAR
min(na.omit(prot_abun_order))
min(na.omit(filtered))
min(imputed_mnar)
min(imputed)



###
#Quality Assessment


##prepare cox ph data

#prepare cox data
coxdata=cbind(t(imputed),metadata_filter)
#join1=joined[,1:3942]
coxdata$time=as.numeric(coxdata$last.known.alive-coxdata$date.of.diagnosis)
class(coxdata$total.protein.concentration)
coxdata$total.protein.concentration=as.numeric(coxdata$total.protein.concentration)

class(coxdata$freeThawCycle)
coxdata$freeThawCycle=as.factor(coxdata$freeThawCycle)
table(coxdata$freeThawCycle)

class(coxdata$batch)
table(coxdata$batch)
coxdata[coxdata$batch=="0","batch"]="7"
coxdata$batch=as.factor(coxdata$batch)

class(coxdata$operator)
table(coxdata$operator)
coxdata[coxdata$operator=="MG/FMA","operator"]="MG"
coxdata$operator=as.factor(coxdata$operator)


class(coxdata$died)
coxdata$event=1
coxdata[coxdata$died==TRUE,"event"]=2


#Univriate coxph

batch_cox <- coxph(Surv(time,event) ~ batch, data = coxdata)
op_cox <- coxph(Surv(time,event) ~ operator, data = coxdata)
thaw_cox <- coxph(Surv(time,event) ~ freeThawCycle, data = coxdata)
conc_cox <- coxph(Surv(time,event) ~ total.protein.concentration, data = coxdata)



cox_summary=function(x){ 
  x <- summary(x)
  p.value<-signif(x$wald["pvalue"], digits=2)
  wald.test<-signif(x$wald["test"], digits=2)
  beta<-signif(x$coef[1], digits=2);#coeficient beta
  HR <-signif(x$coef[2], digits=2);#exp(beta)
  
  res<-c(beta, HR, wald.test, p.value)
  names(res)<-c("beta", "HR ", "wald.test", 
                "p.value")
  return(res)
  
}


univariate_summary=data.frame(rbind(cox_summary(op_cox),cox_summary(thaw_cox),cox_summary(conc_cox),cox_summary(batch_cox)),row.names = c("operator","thaw","protein concentration","batch"))

ggsurvplot(survfit(op_cox),data =coxdata,conf.int = FALSE)
new_df=data.frame("operator"=c("CB","MG"),"freeThawCycle"=c("0"),"total.protein.concentration"=mean(coxdata$total.protein.concentration),"batch"="3")
fit=survfit(op_cox,newdata = new_df)
op_survplot=ggsurvplot(fit, legend.labs=c("CB","MG"),data = coxdata,ggtheme = theme_minimal(),conf.int = FALSE,pval = TRUE)

op_survplot

ggsurvplot(survfit(thaw_cox),data =coxdata,conf.int = FALSE)
new_df=data.frame("operator"=c("CB"),"freeThawCycle"=c("0","1"),"total.protein.concentration"=mean(coxdata$total.protein.concentration),"batch"="3")
fit=survfit(thaw_cox,newdata = new_df)
thaw_survplot=ggsurvplot(fit, legend.labs=c("0","1"),data = coxdata,ggtheme = theme_minimal(),conf.int = FALSE,pval = TRUE)

thaw_survplot


##conclusion: technical factor are not having effect on survival modelling


##PCA analysis to check effect of technical factors on abundance based clusters

t_imp=t(imputed)
pca_imp=prcomp(t_imp,scale. = T,center = T)
autoplot(pca_imp)

autoplot(pca_imp, data=coxdata,colour='batch')
autoplot(pca_imp, data=coxdata,colour='died')
autoplot(pca_imp, data=coxdata,colour='freeThawCycle')
autoplot(pca_imp, data=coxdata,colour='total.protein.concentration')
autoplot(pca_imp, data=coxdata,colour='operator')

##No obvious cluster separation observed in samples based on technical factors


### Protein Markers ###

BiocManager::install("RegParallel")

library(RegParallel)
#just checking which columns to remove
colnames(coxdata)[dim(imputed)[1]:dim(coxdata)[2]]

protein_coxdata=coxdata[,! colnames(coxdata) %in% c("batch","operator","total.protein.concentration","conc","died","last.known.alive","date.of.diagnosis","freeThawCycle")]

colnames(protein_coxdata)=make.names(colnames(protein_coxdata))


prognosis <- RegParallel(
  data = protein_coxdata,
  formula = 'Surv(time, event) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(protein_coxdata)[1:(ncol(protein_coxdata)-2)],
  blocksize = 200,
  cores = 3,
  nestedParallel = FALSE,
  conflevel = 95)

prognosis2=prognosis[!is.na(prognosis$P),]
prognosis2

prognosis2 <- prognosis2[order(prognosis2$LogRank, decreasing = FALSE),]
signif_feature <- subset(prognosis2, LogRank < 0.01)
signif_feature2 <- subset(prognosis2, LogRank < 0.05)

library(stringr)

prognostic_proteins <- str_split_fixed(signif_feature2$Variable,pattern = "\\.",3)[,2]
prognostic_gene     <- str_remove(str_split_fixed(signif_feature2$Variable,pattern = "\\.",3)[,3],"_HUMAN")



###Pathway Enrichment

library(org.Hs.eg.db)

hs <- org.Hs.eg.db
my.symbols <- c(prognostic_gene)
entrez_map=select(hs, 
                  keys = my.symbols,
                  columns = c("ENTREZID", "SYMBOL"),
                  keytype = "SYMBOL")

prognostic_gene_entrez=entrez_map[!is.na(entrez_map$ENTREZID),2]


#reactomePA
x <- enrichPathway(gene=prognostic_gene_entrez,
                   pvalueCutoff=0.2,
                   qvalueCutoff=0.2,
                   readable=T)


reactome_pa_res=x@result

print(reactome_pa_res[1:15,c(2,5,6,8)])


##pathfindr
library(pathfindR)

pfr_input=signif_feature2[,c(1,3,9)]
pfr_input$Variable=str_remove(str_split_fixed(pfr_input$Variable,pattern = "\\.",3)[,3],"_HUMAN")
colnames(pfr_input)=c("Gene.symbol","logFC","adj.P.Val")

pfr_out <- run_pathfindR(pfr_input)

