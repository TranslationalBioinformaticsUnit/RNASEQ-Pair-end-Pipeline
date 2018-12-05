## screen -r InCure

setwd("~/Dropbox/WORKING_DOCUMENTS/PROJECT_SUPERVISION/2017_InnCure/CODE")
setwd("~/Dropbox/WORKING_DOCUMENTS/PROJECT_SUPERVISION/2017_InnCure/InVivo_2018_08")



###############
### END STEP 1: READ DATA

M.union <- read.csv("DataCountTables/CountTable_InCure_D2_Vivo.txt",sep="\t",row.names=1)
SAMPLE_DESCRIPTION <- read.csv("DataCountTables/DESIGN_InCure_D2_Vivo.txt",sep="\t")

###############
### STEP 3: SELECT THE GENES WITH ENOUGH READS
# ACTION: DEFINE THE LIMIT OF COMPUTATION.
# By now: AT LEAST 10 reads in at least 20 samples
# Only using Intersect

library('biomaRt')
library("NOISeq")

martmmusculus <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")  #
type.mart<-"ensembl_gene_id"
biotypes<-getBM(attributes=c("ensembl_gene_id","chromosome_name","strand",
                             "start_position","end_position","gene_biotype"),
                filters="ensembl_gene_id",values=rownames(M.union),mart=martmmusculus)
save(biotypes,file="biotypesVIVO_2018_08.Rda")
rownames(biotypes) <- biotypes[,1]
#gene.INFO<-read.csv("../GENE_DESCRIPTION_red.txt",sep="\t")

gene.INFO<-read.csv("DATA/Gene_details.txt",sep="\t")
sum(rownames(M.union) %in% as.character(gene.INFO[,1])) / nrow(M.union)


############## 
### CHECK WITH CPM > 1 in at least 3 samples.
require(edgeR)
M.go.cpm<-cpm(M.union)
keep = rowSums(M.go.cpm > 1) >= 4 
M.goF2 <- M.union[keep,]
dim(M.goF2)
sum(rownames(M.goF2) %in% as.character(gene.INFO[,1])) / nrow(M.goF2)
M.goF3 <- M.goF2[rownames(M.goF2) %in% as.character(gene.INFO[,1]),]
sum(rownames(M.goF3) %in% as.character(gene.INFO[,1])) / nrow(M.goF3)
gene.INFOF3 <- gene.INFO[as.character(gene.INFO[,1]) %in% rownames(M.goF3),]

a1 <- table(as.character(gene.INFOF3[,1]))
ENSMUSG00000111375
gene.INFOF3[as.character(gene.INFOF3[,1]) %in% "ENSMUSG00000111375",]
gene.INFOF3<- gene.INFOF3[-which(as.character(gene.INFOF3[,7])=="A830010M20Rik"),]

nrow(gene.INFOF3)
nrow(M.goF3)
rownames(gene.INFOF3) <- as.character(gene.INFOF3[,1])
gene.INFOF3 <- gene.INFOF3[rownames(M.goF3),]


save(M.goF3,gene.INFOF3,SAMPLE_DESCRIPTION,file="RNASeq_PROCESS_STEP2_Incure_Prenorm_inVivo_2018_08.Rda")

##############################
##############################
#####  STEP 5: NORMALIZATION
##############################
##############################

load("RNASeq_PROCESS_STEP2_Incure_Prenorm_inVivo_2018_08.Rda")
  # M.goF3,gene.INFOF3,SAMPLE_DESCRIPTION

###### 
## MDS

dist_go <- 1-cor(M.goF3,method="spearman")

fit <- cmdscale(dist_go, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y,
     xlim=c(-1,1),
     ylim=c(-1,1))#,
     #col=TYPE_col,
     #pch=HOUR_col)

names_samples <- colnames(M.goF3)
text(x, y, pos = 4, labels = names_samples)

###### 
## DISCARD SAMPLES

SAMPLE_DESCRIPTION_FILTERED <- SAMPLE_DESCRIPTION
rownames(SAMPLE_DESCRIPTION_FILTERED) <- SAMPLE_DESCRIPTION_FILTERED$NAME
rownames(SAMPLE_DESCRIPTION_FILTERED) == colnames(M.goF3)

S_18_0863_21-A16217-6-m-APPPS1_N702-S517
APPPS1_6m__6_1_21

SAMPLE_DESCRIPTION_FILTERED <- SAMPLE_DESCRIPTION_FILTERED[-which(rownames(SAMPLE_DESCRIPTION_FILTERED)=="APPPS1_6m__6_1_21"),]
M.goF3 <- M.goF3[,-which(colnames(M.goF3)=="APPPS1_6m__6_1_21"),]

###### 
## MDS

dist_go <- 1-cor(M.goF3,method="spearman")

fit <- cmdscale(dist_go, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y,
     xlim=c(-0.1,0.1),
     ylim=c(-0.1,0.1))#,
#col=TYPE_col,
#pch=HOUR_col)

names_samples <- colnames(M.goF3)
text(x, y, pos = 4, labels = names_samples)

save(M.goF3,SAMPLE_DESCRIPTION_FILTERED,gene.INFOF3,file="RNASeq_PROCESS_STEP3_Incure_Filtered_inVivo_2018_08.Rda")


#######
## ORDER SAMPLES

SAMPLE_DESCRIPTION_go <- SAMPLE_DESCRIPTION_FILTERED[order(SAMPLE_DESCRIPTION_FILTERED$MUT,SAMPLE_DESCRIPTION_FILTERED$Mxo4,SAMPLE_DESCRIPTION_FILTERED$TIME),]
M.goF4                <- M.goF3[,rownames(SAMPLE_DESCRIPTION_go)]

require(edgeR)
library(cqn)
library(scales)
library(limma)

sizeFactors.GO<-apply(M.goF4,2,sum)

BATCH<-SAMPLE_DESCRIPTION_go$LIB_PREP
LANE<-SAMPLE_DESCRIPTION_go$LANE
CONDITION<-paste(SAMPLE_DESCRIPTION_go$MUT,
                 SAMPLE_DESCRIPTION_go$Mxo4,
                 SAMPLE_DESCRIPTION_go$TIME,
                 sep="_")
TYPE       <- as.character(SAMPLE_DESCRIPTION_go$MUT)
VAR  <- as.character(SAMPLE_DESCRIPTION_go$Mxo4)
MONTH       <- as.character(SAMPLE_DESCRIPTION_go$TIME) 

TYPE_col <- TYPE
TYPE_col[TYPE_col=="WT"] <- "blue"
TYPE_col[TYPE_col=="NLRP3"] <- "red"
TYPE_col[TYPE_col=="APPPS1"] <- "orange"
TYPE_col[TYPE_col=="APPPS1_NLRP3"] <- "grey"


BATCH_col <- as.character(BATCH)
BATCH_col[BATCH_col=="1"] <- "green"
BATCH_col[BATCH_col=="2"] <- "orange"

VAR_col <- VAR
VAR_col[VAR_col==""] <- "black"
VAR_col[VAR_col=="Mxo4"] <- "grey"

MONTH_col <- MONTH
MONTH_col[MONTH_col=="4m"]  <- 15
MONTH_col[MONTH_col=="6m"] <- 16
MONTH_col[MONTH_col=="12m"] <- 17
MONTH_col <- as.numeric(MONTH_col)

#### CQN
lengthgo<- as.numeric(gene.INFOF3$Gene.end..bp.)-as.numeric(gene.INFOF3$Gene.start..bp.)
cqn.subset <- cqn(M.goF4, lengths = lengthgo,
                  x = as.numeric(as.character(gene.INFOF3$Gene...GC.content)), sizeFactors = sizeFactors.GO,
                  verbose = TRUE)
M.go.ORD.cqn<- cqn.subset$y + cqn.subset$offset
head(M.go.ORD.cqn)
min(M.go.ORD.cqn)
M.go.ORD.cqn <- M.go.ORD.cqn + abs(min(M.go.ORD.cqn)) + 1
min(M.go.ORD.cqn)

#### COMBAT
require(sva)
design.cb<-model.matrix(~ CONDITION )

M.go.ORD.cqn.ADJ<-ComBat(M.go.ORD.cqn, batch=BATCH, mod=design.cb) 
n.sv.CQN = num.sv(M.go.ORD.cqn,design.cb,method="leek")
n.sv.CORRECTED.CQN = num.sv(M.go.ORD.cqn.ADJ,design.cb,method="leek")

save(M.go.ORD.cqn.ADJ,M.go.ORD.cqn,M.goF4,
     SAMPLE_DESCRIPTION_go,gene.INFOF3,
     CONDITION,
     TYPE,TYPE_col,
     VAR,VAR_col,
     MONTH,MONTH_col,
     BATCH,BATCH_col,
     file="RNASeq_PROCESS_STEP4_Incure_NORM_inVivo_2018_08.Rda")
colnames(M.go.ORD.cqn.ADJ) == rownames(SAMPLE_DESCRIPTION_go)


#### 

## CREATE MDS
#name    <- "cqn.combat"

name    <- "cqn.combat_discarded_inVivo_2018_08"
dist_go <- 1-cor(M.go.ORD.cqn.ADJ,method="spearman")

name    <- "cqn"
dist_go <- 1-cor(M.go.ORD.cqn,method="spearman")

fit <- cmdscale(dist_go, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y,
     xlim=c(-0.1,0.1),
     ylim=c(-0.1,0.1),
     col=TYPE_col,
     pch=MONTH_col)

names_samples <- CONDITION
text(x, y, pos = 4, labels = CONDITION)

system("mkdir PLOT_2018_09")
png(paste("PLOT_2018_09/MDS_analysis_",name,"_2018_09.png",sep=""),
      width = 9, height = 9, units = "in",res=300)
par(mfrow = c(2, 2))

plot(x, y,
      xlim=c(-0.075,0.075),
      ylim=c(-0.075,0.075),
      col=TYPE_col,
      pch=MONTH_col)
legend("topleft",
       c("WT","NLRP3","APPPS1","APPPS1_NLRP3"),
       col=c("blue","red","orange","grey"),
       pch=19)
legend("topright",
       c("4m","6m","12m"),
       pch=c(15,16,17),
       col="black")

plot(x, y,
     xlim=c(-0.075,0.075),
     ylim=c(-0.075,0.075),
     col=VAR_col,
     pch=MONTH_col)
legend("topleft",
       c("","Mxo4"),
       col=c("black","grey"),pch=19)
legend("topright",
       c("4m","6m","12m"),
       pch=c(15,16,17),
       col="black")

plot(x, y,
     xlim=c(-0.075,0.075),
     ylim=c(-0.075,0.075),
     col=BATCH_col,
     pch=MONTH_col)
legend("topleft",
       c("Batch1","Batch2"),
       col=c("green","orange"),pch=19)

legend("topright",
       c("4m","6m","12m"),
       pch=c(15,16,17),
       col="black")

dev.off()




png(paste("PLOT_2018_09/HCLUST_analysis_",name,"_2018_09.png",sep=""),
    width = 9, height = 9, units = "in",res=300)

plot(hclust(as.dist(dist_go)))

dev.off()



#######################################
###  QC
#######################################
gene.INFOF3[grep("Pse",as.character(gene.INFOF3$MGI.symbol)),]

Appps1 -> 
  APP    ENSMUSG00000022892
  PSEN1  ENSMUSG00000019969

  DEFINITION:APPPS1 mice contain human transgenes for both APP bearing the Swedish 
             mutation and PSEN1 containing an L166P mutation, both under the control 
             of the Thy1 promoter. In these mice, expression of the human APP transgene 
             is approximately 3-fold higher than endogenous murine APP. Human Aβ42 is 
             preferentially generated over Aβ40, but levels of both increase with age. 
             In the brain, the Aβ42/Aβ40 decreases with the onset of amyloid deposition 
             (Radde et al., 2006; Maia et al., 2013). 
  
library(ggplot2)
library(gridExtra)            

genes_plot <- c("ENSMUSG00000022892","ENSMUSG00000019969","ENSMUSG00000032691")             
genes_plot_mgi <- c("App","Psen1","Nlrp3")             

forplot <- t(as.data.frame(M.go.ORD.cqn.ADJ))
forplot <- cbind(forplot,data.frame(CONDITION,TYPE,BATCH,VAR))

for(i in 1:length(genes_plot))
{#i<-3
 png(paste("PLOT_2018_09/QC_Genes_",genes_plot[i],"_",genes_plot_mgi[i],"_2018_09.png",sep=""),
                 width = 9, height = 9, units = "in",res=300)
   par(mfrow = c(2, 2))
   
   forplot2 <- forplot[,c("CONDITION",genes_plot[i])]
   forplot2$CONDITION <- as.factor(forplot2$CONDITION)
   colnames(forplot2)[2]<- "GENE"
   forplot3 <- forplot[,c("TYPE",genes_plot[i])]
   forplot3$TYPE <- as.factor(forplot3$TYPE)
   colnames(forplot3)[2]<- "GENE"
   forplot4 <- forplot[,c("BATCH",genes_plot[i])]
   forplot4$BATCH <- as.factor(forplot4$BATCH)
   colnames(forplot4)[2]<- "GENE"
   forplot5 <- forplot[,c("VAR",genes_plot[i])]
   forplot5$VAR <- as.factor(forplot5$VAR)
   colnames(forplot5)[2]<- "GENE"
   
   
   p <- ggplot(forplot2, aes(factor(CONDITION), GENE),trim=FALSE)
   p1 <- p +  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  
 
   p <- ggplot(forplot3, aes(factor(TYPE), GENE))
   p2 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)

   p <- ggplot(forplot4, aes(factor(BATCH), GENE))
   p3 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)

   p <- ggplot(forplot5, aes(factor(VAR), GENE))
   p4 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)
   
   grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
   
 dev.off()   
                        
}            

#######################################
###  LOAD THE INFORMATION
#######################################

load("RNASeq_PROCESS_STEP4_Incure_NORM_inVivo_2018_08.Rda")
  #M.go.ORD.cqn.ADJ,M.go.ORD.cqn,M.goF4,
  #SAMPLE_DESCRIPTION_go,gene.INFOF3,
  #CONDITION,
  #TYPE,TYPE_col,
  #VAR,VAR_col,
  #MONTH,MONTH_col,
  #BATCH,BATCH_col,


#######################################
###  DIFFERENTIAL EXPRESSION
#######################################

design<-model.matrix(~0 + CONDITION)
colnames(design) <- gsub("CONDITION","",colnames(design))
fit <- lmFit(M.go.ORD.cqn.ADJ,design)#,weigths=v$weights)
cont.matrix.TIME <- makeContrasts(
  WT_4_6=WT__6m-WT__4m,
  WT_6_12=WT__12m-WT__6m,
  WT_4_12=WT__12m-WT__4m,
  NLRP3_4_6=NLRP3__6m-NLRP3__4m,
  NLRP3_6_12=NLRP3__12m-NLRP3__6m,
  NLRP3_4_12=NLRP3__12m-NLRP3__4m,
  APPPS1_4_6=APPPS1__6m-APPPS1__4m,
  APPPS1_6_12=APPPS1__12m-APPPS1__6m,
  APPPS1_4_12=APPPS1__12m-APPPS1__4m,
  APPPS1_NLRP3_4_6=APPPS1_NLRP3__6m-APPPS1_NLRP3__4m,
  APPPS1_NLRP3_6_12=APPPS1_NLRP3__12m-APPPS1_NLRP3__6m,
  APPPS1_NLRP3_4_12=APPPS1_NLRP3__12m-APPPS1_NLRP3__4m,
  WT_Mox4=WT_Mxo4_12m-WT__12m,
  NLRP3_Mox4=NLRP3_Mxo4_12m-NLRP3__12m,
  APPPS1_Mox4=APPPS1_Mxo4_12m-APPPS1__12m,
  APPPS1_NLRP34_Mox4=APPPS1_NLRP3_Mxo4_12m-APPPS1_NLRP3__12m,
  WT_NLRP3_4=NLRP3__4m-WT__4m,
  WT_APPPS1_4=APPPS1__4m-WT__4m,
  WT_APPPS1_NLRP3_4=APPPS1_NLRP3__4m-WT__4m,
  APPPS1_NLRP3_4=NLRP3__4m-APPPS1__4m,
  APPPS1_APPPS1_NLRP3_4=APPPS1_NLRP3__4m-APPPS1__4m,
  NLRP3_APPPS1_NLRP3_4=APPPS1_NLRP3__4m-NLRP3__4m,
  #add 6m
  WT_NLRP3_6=NLRP3__6m-WT__6m,
  WT_APPPS1_6=APPPS1__6m-WT__6m,
  WT_APPPS1_NLRP3_6=APPPS1_NLRP3__6m-WT__6m,
  APPPS1_NLRP3_6=NLRP3__6m-APPPS1__6m,
  APPPS1_APPPS1_NLRP3_6=APPPS1_NLRP3__6m-APPPS1__6m,
  NLRP3_APPPS1_NLRP3_6=APPPS1_NLRP3__6m-NLRP3__6m,
  #add 12m
  WT_NLRP3_12=NLRP3__12m-WT__12m,
  WT_APPPS1_12=APPPS1__12m-WT__12m,
  WT_APPPS1_NLRP3_12=APPPS1_NLRP3__12m-WT__12m,
  APPPS1_NLRP3_12=NLRP3__12m-APPPS1__12m,
  APPPS1_APPPS1_NLRP3_12=APPPS1_NLRP3__12m-APPPS1__12m,
  NLRP3_APPPS1_NLRP3_12=APPPS1_NLRP3__12m-NLRP3__12m,
  Incr_WT_NLRP3_4_6=(NLRP3__6m-NLRP3__4m)-(WT__6m-WT__4m),
  Incr_WT_NLRP3_6_12=(NLRP3__12m-NLRP3__6m)-(WT__12m-WT__6m),
  Incr_WT_NLRP3_4_12=(NLRP3__12m-NLRP3__4m)-(WT__12m-WT__4m),
  Incr_WT_APPPS1_4_6=(APPPS1__6m-APPPS1__4m)-(WT__6m-WT__4m),
  Incr_WT_APPPS1_6_12=(APPPS1__12m-APPPS1__6m)-(WT__12m-WT__6m),
  Incr_WT_APPPS1_4_12=(APPPS1__12m-APPPS1__4m)-(WT__12m-WT__4m),
  Incr_WT_APPPS1_NLRP3_4_6=(APPPS1_NLRP3__6m-APPPS1_NLRP3__4m)-(WT__6m-WT__4m),
  Incr_WT_APPPS1_NLRP3_6_12=(APPPS1_NLRP3__12m-APPPS1_NLRP3__6m)-(WT__12m-WT__6m),
  Incr_WT_APPPS1_NLRP3_4_12=(APPPS1_NLRP3__12m-APPPS1_NLRP3__4m)-(WT__12m-WT__4m),
  levels=design)
fit.TIME <- contrasts.fit(fit, cont.matrix.TIME)
fit.TIME <- eBayes(fit.TIME)

Diff_TIME.ALL<-topTable(fit.TIME,coef=1:43,number=nrow(M.go.ORD.cqn.ADJ),sort.by="none")

Diff_InVivoAll <- Diff_TIME.ALL[,c("P.Value","adj.P.Val")]

for(i in 1:43)
{#i<-1
  rm(Diff_TIMEgo)
  Diff_TIMEgo<-topTable(fit.TIME,coef=i,
                        number=nrow(M.go.ORD.cqn.ADJ),
                        sort.by="none")
  Diff_TIMEgo<-Diff_TIMEgo[,c("logFC","t","P.Value","adj.P.Val")]
  colnames(Diff_TIMEgo) <- paste(colnames(cont.matrix.TIME)[i],colnames(Diff_TIMEgo),sep="_")
  Diff_InVivoAll <- cbind(Diff_InVivoAll,Diff_TIMEgo)
}

Diff_InVivoAll <- cbind(gene.INFOF3[rownames(Diff_InVivoAll),"MGI.symbol"],Diff_InVivoAll)
colnames(Diff_InVivoAll)[1]<-"MGI"

save(Diff_InVivoAll,file="RNASeq__Dif_InVivoAll_2018_10.Rda")

write.table(Diff_InVivoAll,"2Diff_InVivoAll_summaryA_InVivo_2018_08.txt",
            sep="\t",quote=F,row.names = T)
######################################################
############ TOP 20 GENES FOR EACH CONTRAST
library(ggplot2)
library(gridExtra)
forplot <- t(as.data.frame(M.go.ORD.cqn.ADJ))
forplot <- cbind(forplot,data.frame(CONDITION,TYPE,BATCH,VAR))

Diff_cols_adjp<-grep("_adj.P",colnames(Diff_InVivoAll))

col<-1
i<-1
rm(Diff_top20, Diff_top20_1)
for(i in 1:43)
{
  #i=1
  condition_name<-strsplit(colnames(Diff_InVivoAll)[Diff_cols_adjp[i]], "_adj.P.Val")
  contrast<-t(t(rep(condition_name[1], 20)))
   if(i>1){
     Diff_top20_1<-cbind(rownames(Diff_InVivoAll[order(Diff_InVivoAll[,Diff_cols_adjp[i]]),])[1:20], 
                       Diff_InVivoAll[order(Diff_InVivoAll[,Diff_cols_adjp[i]])[1:20],c(1, Diff_cols_adjp[i]-3,Diff_cols_adjp[i]-2, Diff_cols_adjp[i]-1, Diff_cols_adjp[i])],
                       contrast)
     colnames(Diff_top20_1)<-c("Ensembl ID", "MGI","logFC","t.Value", "P.Value", "adj.P.Val", "contrast")
     Diff_top20<-rbind(Diff_top20, Diff_top20_1) 
   }
   else{
     Diff_top20<-cbind(rownames(Diff_InVivoAll[order(Diff_InVivoAll[,Diff_cols_adjp[i]]),])[1:20],
                       Diff_InVivoAll[order(Diff_InVivoAll[,Diff_cols_adjp[i]])[1:20],c(1, Diff_cols_adjp[i]-3,Diff_cols_adjp[i]-2, Diff_cols_adjp[i]-1, Diff_cols_adjp[i])],
                       contrast)
     colnames(Diff_top20)<-c("Ensembl ID", "MGI","logFC","t.Value", "P.Value", "adj.P.Val", "contrast")         
   
  }
  # for (gene in 1:20){
  #   if(gene==1){
  #     if(!dir.exists(paste("PLOT_TOP20_GENES/", condition_name[1],  sep=""))){
  #       dir.create(paste("PLOT_TOP20_GENES/", condition_name[1],  sep=""), mode="0777")
  #     }
  #   }
  #   gene_to_plot<-rownames(Diff_InVivoAll[order(Diff_InVivoAll[,Diff_cols_adjp[i]]),])[gene]
  #   png(paste("PLOT_DIFF_TOP20_GENES/", condition_name[1],"/","GENE_", gene_to_plot
  #             ,"_",Diff_InVivoAll[gene_to_plot,1],"_2018_10.png",sep=""),
  #       width = 9, height = 9, units = "in",res=300)
  #     par(mfrow = c(2, 2))
  #     
  #     forplot2 <- forplot[,c("CONDITION",gene_to_plot)]
  #     forplot2$CONDITION <- as.factor(forplot2$CONDITION)
  #     colnames(forplot2)[2]<- "GENE"
  #     forplot3 <- forplot[,c("TYPE",gene_to_plot)]
  #     forplot3$TYPE <- as.factor(forplot3$TYPE)
  #     colnames(forplot3)[2]<- "GENE"
  #     forplot4 <- forplot[,c("BATCH",gene_to_plot)]
  #     forplot4$BATCH <- as.factor(forplot4$BATCH)
  #     colnames(forplot4)[2]<- "GENE"
  #     forplot5 <- forplot[,c("VAR",gene_to_plot)]
  #     forplot5$VAR <- as.factor(forplot5$VAR)
  #     colnames(forplot5)[2]<- "GENE"
  #     
  #     
  #     p <- ggplot(forplot2, aes(factor(CONDITION), GENE),trim=FALSE)
  #     p1 <- p +  geom_violin(trim=FALSE) + geom_boxplot(width=0.1)  
  #     
  #     p <- ggplot(forplot3, aes(factor(TYPE), GENE))
  #     p2 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)
  #     
  #     p <- ggplot(forplot4, aes(factor(BATCH), GENE))
  #     p3 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)
  #     
  #     p <- ggplot(forplot5, aes(factor(VAR), GENE))
  #     p4 <- p + geom_violin(trim=FALSE)  + geom_boxplot(width=0.1)
  #     
  #     grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)
  #   
  #   dev.off() 
  # }
  
}
rm(condition_name, col, i)
rownames(Diff_top20)<-NULL
#dim(20, 258)
dim(Diff_top20)
save(Diff_top20,file="Diff_top20_genes2.Rda")

Diff_top20$contrast<-unlist(Diff_top20$contrast, use.names= FALSE)
write.table(Diff_top20,"Diff_InVivoAll_top20genes.txt", sep="\t",quote=F,row.names = F, col.names = T)


###########################################################
############ ANALYSIS
Diff_InVivoAlla <- Diff_InVivoAll[,grep("adj",colnames(Diff_InVivoAll)),]

num_apval <- matrix(NA,31,5)
rownames(num_apval) <- colnames(cont.matrix.TIME)
numbersgo <- c(0.1,0.05,0.01,0.005,0.001)
colnames(num_apval)  <- as.character(numbersgo)

for(i in 1:nrow(num_apval))
{
  for(j in 1:length(numbersgo))
  {
    num_apval[i,j] <- sum(Diff_InVivoAlla[,i+1]<numbersgo[j])
  }
}



############ PLOT SPECIFICS

system("mkdir PLOTS_logFC_2018_06")

gene_nams <- c("Tnfsf18","Tnfsf4","Il27","Spint1","Arg1",
               "Pcdhgb2","Pcdhgb5","Pcdhgb6","Pcdhga4","Pcdhga8","Pcdhgc4",
               "Cdh24","Cdhr4","Cdh22","Cdh11","Cdh13",
               "Csf1","Csf2","Csf3",
               "Mybl2","Plk1","Foxm1","Ccne1","Ccnd1","Mcm10","Mcm9",
               "Trem2","Ptk2b","Apoc1","Plxna4",
               "Abca7","Atxn1","Bin1","Cd33","Glod4","Inpp5d","Picalm","Zcwpw1")

              # "","","","","","","","","","",)


Diff_TIME.A <- Diff_TIME.summary2[,grep("logFC",colnames(Diff_TIME.summary2))]

Diff_TIME.A <- Diff_TIME.A[!(Diff_TIME.summary2[,1]=="") & 
                             !(is.na(Diff_TIME.summary2[,1])) &
                             !(as.character(Diff_TIME.summary2[,1]) %in% c("Crybg3","Itgam","Lilrb4a","U2af1l4")),]

colnames(Diff_TIME.A) <- paste("Contrast_A_",colnames(Diff_TIME.A),sep="")

Diff_TIME.B <- Diff_TIME.summaryT2[,grep("logFC",colnames(Diff_TIME.summaryT2))]

Diff_TIME.B <- Diff_TIME.B[!(Diff_TIME.summaryT2[,1]=="") & 
                             !(is.na(Diff_TIME.summaryT2[,1])) &
                             !(as.character(Diff_TIME.summaryT2[,1]) %in% c("Crybg3","Itgam","Lilrb4a","U2af1l4")),]

colnames(Diff_TIME.B) <- paste("Contrast_B_",colnames(Diff_TIME.B),sep="")


Diff_TIME.C <- Diff_TIME.summaryNLRP3_DMSO_minus_WT_DMSO[,grep("logFC",colnames(Diff_TIME.summaryT2))]

Diff_TIME.C <- Diff_TIME.C[!(Diff_TIME.summaryNLRP3_DMSO_minus_WT_DMSO[,1]=="") & 
                             !(is.na(Diff_TIME.summaryNLRP3_DMSO_minus_WT_DMSO[,1])) &
                             !(as.character(Diff_TIME.summaryNLRP3_DMSO_minus_WT_DMSO[,1]) %in% c("Crybg3","Itgam","Lilrb4a","U2af1l4")),]
colnames(Diff_TIME.C) <- paste("Contrast_C_",colnames(Diff_TIME.C),sep="")

sum(rownames(Diff_TIME.A)==rownames(Diff_TIME.B))/nrow(Diff_TIME.B)
sum(rownames(Diff_TIME.B)==rownames(Diff_TIME.C))/nrow(Diff_TIME.C)

Dif_LogFC <- cbind(Diff_TIME.A,Diff_TIME.B,Diff_TIME.C)
rownames(Dif_LogFC) <- Diff_TIME.summary2[rownames(Dif_LogFC),1]

colnames(Dif_LogFC)<-gsub("Contrast_","",colnames(Dif_LogFC))
colnames(Dif_LogFC)<-gsub("_logFC","",colnames(Dif_LogFC))


for(i in 1:length(gene_nams))
{#i <- 1
  png(paste("PLOTS_logFC_2018_06/GENE",gene_nams[i],"_2018_06.png",sep=""),
      width = 12, height = 6, units = "in",res=200)
  par(mar=c(14, 5, 5, 5))
  t0 <- as.numeric(Dif_LogFC[gene_nams[i],])
  valuego <- abs(max(c(-min(t0),max(t0))))*1.1
  barplot(t0,las=3,
          ylim=c(-valuego,valuego),names.arg =colnames(Dif_LogFC),
          ylab="logFC",main=gene_nams[i])
  
  dev.off()
  
}


for(i in 1:nrow(Dif_LogFC))
{#i <- 1
  png(paste("PLOTS_logFC_2018_06_ALL/GENE",rownames(Dif_LogFC)[i],"_2018_06.png",sep=""),
      width = 12, height = 6, units = "in",res=200)
  par(mar=c(14, 5, 5, 5))
  t0 <- as.numeric(Dif_LogFC[i,])
  valuego <- abs(max(c(-min(t0),max(t0))))*1.1
  barplot(t0,las=3,
          ylim=c(-valuego,valuego),names.arg =colnames(Dif_LogFC),
          ylab="logFC",main=rownames(Dif_LogFC)[i])
  
  dev.off()
  
}




########### POWER ANALYSSI
  
require(pwr)
power_info <- matrix(NA,200000,2)
go <- 0
for(i in 1:nrow(Dif_LogFC))
{
  for(h in 1:ncol(Dif_LogFC))
  {
    go <- go+ 1 
    power_info[go,1] <- pwr.t.test(d=Dif_LogFC[i,h],
                                   n=4,sig.level=0.01,type="one.sample",
                                   alternative="two.sided")$power   
    power_info[go,2] <- colnames(Dif_LogFC)[h]   
  }
  
}
power_info <- power_info[1:go,]
power_info[,2] <- gsub("logFC_","",power_info[,2])

png(paste("PLOT_2018_06/POWER_ANALYSIS","_replot_2018_06.png",sep=""),
    width = 12, height = 6, units = "in",res=300)
par(mar=c(14, 5, 5, 5))
boxplot(as.numeric(power_info[,1]) ~ power_info[,2],las=2)
dev.off()


require("SSPA")

Diff_TIME1<-topTable(fit.TIME,coef=1,number=nrow(M.go.ORD.cqn.ADJ),sort.by="none")
Diff_TIME2<-topTable(fit.TIME,coef=2,number=nrow(M.go.ORD.cqn.ADJ),sort.by="none")
Diff_TIME3<-topTable(fit.TIME,coef=3,number=nrow(M.go.ORD.cqn.ADJ),sort.by="none")


pdD <-  pilotData(statistics = Diff_TIME3$t,
                  samplesize = 4,#sqrt(1/(1/4 +1/4)),
                  distribution="t",df=3)
ssD <- sampleSize(pdD,method='congrad',
                  control=list(from=-6, to=6))#, resolution=250))#, control=list(from=-6, to=6))
Jpred <- seq(10, 20, by=2)

pwrD <- predictpower(ssD, samplesizes=sqrt(Jpred/2), alpha=0.05)
matplot(Jpred, pwrD, type="b", pch=16, ylim=c(0, 1),
        ylab="predicted power", xlab="sample size (per group)")
grid()



