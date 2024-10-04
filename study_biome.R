# rm(list = ls())

require('readr')
set.seed(12345)

setwd('C:/Users/DELL/Desktop/blastim/project/')
biomeData <- read.table('otu_table_L7.txt', header = T)
biomeData <- as.data.frame(t(biomeData))

### 2) ПОЛУЧЕНИЕ ПРОЦЕНТОВ ПРЕДСТАВЛЕННОСТИ

#count value*100%/10 000 rarefied depth of counts = count value/100
biomeData <- biomeData / 100
###

library(dplyr)
library(tidyverse)

biomeData <- biomeData %>%
  add_column(run_accession = rownames(biomeData), .before = 1)

rownames(biomeData) <- NULL

biomeData$run_accession <- sub("\\..*", "", biomeData$run_accession) 

metaData <- read_delim('metadata.tsv.txt',  "\t", col_names = TRUE, escape_double = FALSE, trim_ws = TRUE)
metaData$control_case <- sub(".*(normal|tumor)$", "\\1", metaData$sample_alias) 
alphaDivData <- read_delim('alpha_chao.tsv',  "\t", col_names = TRUE,escape_double = FALSE, trim_ws = TRUE)
colnames(alphaDivData)[1]<-colnames(metaData)[1]
alphaDivData$run_accession <- sub("\\..*", "", alphaDivData$run_accession) 


commonSamples<-Reduce(intersect, list(biomeData$run_accession,alphaDivData$run_accession,metaData$run_accession))
commonSamples

rownames(biomeData) <-biomeData$run_accession
rownames(alphaDivData)<- alphaDivData$run_accession
rownames(metaData)<-metaData$run_accession


alphaDivDataS<-alphaDivData[commonSamples,]

biomeDataS<-biomeData[commonSamples,]

metaDataS<-metaData[commonSamples,]


biomeDataS<-biomeDataS[,-1]


#######################Stat tests

#check normality
str(alphaDivDataS)
shapiro.test(alphaDivDataS$chao1) #p<0.05 means that it is not a normal distribution of data

#U-test
tumor_run_accession <- metaDataS[grepl("tumor", metaDataS$sample_alias), "run_accession"]
normal_run_accession <- metaDataS[grepl("normal", metaDataS$sample_alias), "run_accession"]
# Print the result
tumor_run_accession
normal_run_accession

cases<-tumor_run_accession[[1]]
cases

controls<-normal_run_accession[[1]]

controls

#### 3) ДИФФЕРЕНЦИАЛЬНЫЙ АНАЛИЗ, выявлены различающиеся таксоны в группах: - больные и здоровые 
wilcoxRes<-matrix(ncol=3,nrow=0)


for (i in colnames(biomeDataS)){ #для каждого таксона я делаю тест
  wt<-wilcox.test(biomeDataS[cases,i],biomeDataS[controls,i])
  wilcoxRes<-rbind(wilcoxRes, c(i,wt$statistic, wt$p.value))
}


wicoxPvalAdj<-p.adjust(wilcoxRes[,3], method = 'fdr')
wilcoxRes<-cbind(wilcoxRes, wicoxPvalAdj)
colnames(wilcoxRes)<-c('tax','stat','pval','pval_adj')
wilcoxRes<-as.data.frame(wilcoxRes)
wilcoxRes$pval<- round(as.numeric(as.character(wilcoxRes$pval)), 3)
wilcoxRes$pval_adj<- round(as.numeric(as.character(wilcoxRes$pval_adj)), 3)
wilcoxRes
wilcoxResSign<-wilcoxRes[which(wilcoxRes$pval_adj <0.5),]
wilcoxResSign
wt<-wilcox.test(alphaDivDataS[which(alphaDivDataS$run_accession %in% cases),'chao1'][[1]], alphaDivDataS[which(alphaDivDataS$run_accession %in% controls),'chao1'][[1]])
wt
#p-value = 1, Since 0.1299 > 0.05, there is not enough evidence to say that there is the  difference in alpha diversity.

### МОДЕЛЬ-КЛАССИФИКАТОР для определения здоровой и патологической микробиоты 
#Generalized linear model
require(dplyr)
glmDF<-inner_join(metaDataS[,c('run_accession','control_case')], biomeData, by = 'run_accession')

resGLM<-matrix(nrow=0, ncol=6)

for (i in colnames(biomeDataS)) {
  model<- glm(glmDF[,i][[1]] ~ glmDF[,'control_case'][[1]])
  tr<-summary(model)
  tr<-tr$coefficients
  tr[,'Pr(>|t|)']<-round(as.numeric(tr[,'Pr(>|t|)']),2)
  tr<-cbind(rep(i,nrow(tr)),tr)
  tr<-cbind(c('control_case'),tr)
  rownames(tr)<-c()
  resGLM<-rbind(resGLM, tr)
}

colnames(resGLM)<-c('factor','tax','estimate','std_error','t_val','pval')
resGLM<-as.data.frame(resGLM)
resGLM<-cbind(resGLM, p.adjust(resGLM$pval, method = 'fdr'))
colnames(resGLM)[length(colnames(resGLM))]<-'pval_adj'

# resGLM<-resGLM[which(resGLM$factor == 'control_vs_case'),]
resGLMFilt<-resGLM[which(resGLM$pval_adj<0.05),]
resGLMFilt$tax[which(resGLMFilt$tax %in% wilcoxResSign$tax)]
View(resGLMFilt)


#PERMANOVA
# install.packages("PERMANOVA")
require(PERMANOVA)

biomeDist<- DistContinuous(biomeDataS, coef = 'Bray_Curtis')
biomePERM=PERMANOVA(biomeDist, as.factor(metaDataS$control_case),  CoordPrinc = TRUE, PostHoc = 'fdr')
biomePERM$pvalue
summary(biomePERM)


#Boxplot
require(data.table)
dfGG<-biomeDataS[,wilcoxResSign$tax]
dfGG<-cbind(rownames(biomeDataS),dfGG)
colnames(dfGG)[1]<-'run_accession'
dfGGm<-melt(dfGG, id.vars = 1)


plotBplot <- ggplot(dfGGm, aes(x=variable, y=value )) + geom_boxplot(outlier.colour="red", outlier.shape=16,outlier.alpha=0.5, outlier.size=0, notch=FALSE)+geom_jitter(alpha=0.5)+scale_x_discrete(guide = guide_axis(angle = 90))+ theme_bw()
plotBplot

dfGGm<-join(dfGGm,metaDataS[,c('run_accession','control_case')], by="run_accession")

plotBplot <- ggplot(dfGGm, aes(x=variable, y=value, fill=control_case) ) + geom_boxplot(outlier.colour="red", outlier.shape=16,outlier.alpha=0.5, outlier.size=0, notch=FALSE)+geom_jitter(alpha=0.1, color='green')+scale_x_discrete(guide = guide_axis(angle = 90))+ theme_bw()
plotBplot

#MDS
require(MASS)


myMDS<-isoMDS(biomeDistBC)
dfMDS <- data.frame(X = as.vector(myMDS$points[,1]), Y = as.vector(myMDS$points[,2]), run_accession = rownames(myMDS[[1]]))
plotMDS <- ggplot(dfMDS, aes(x=X, y=Y,label = run_accession)) + geom_point(size=1) +geom_text(hjust = 0, nudge_x = 0.05) + theme_bw()
plotMDS
library(dplyr)
dfMDS<-full_join(dfMDS, metaDataS, by = 'run_accession')
plotMDS <- ggplot(dfMDS, aes(x=X, y=Y, color=control_case)) + geom_point(size=1) + theme_bw()
plotMDS

