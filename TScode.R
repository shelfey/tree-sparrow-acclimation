
###Determine the candidate genes 

### WGCNA analysis 
### !!!WGCNA analysis using linear model is according to the regular code of WGCNA
### WGCNA analysis using non-linear model

setwd("F:\\PAPER\\bulk_rna_seq\\NH_treated\\Fgenes\\flight\\WGCNA")
library(BiocGenerics)
library(org.Gg.eg.db)
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

tpm<-read.csv("flight.tpm.csv",header = T,sep = ",")

##filter genes with low level of expression 
data <- tpm[rowSums(as.matrix(tpm[,2:5])>0) >=3 & 
              rowSums(as.matrix(tpm[,6:9])>0)>=3 & 
              rowSums(as.matrix(tpm[,10:13])>0)>=3,]

dat<- data[apply(data[,-1], 1, mean) > 1,] 

###### 
titData <- dat
colnames(titData)[1]<-"Genes"
head(titData)
dim(titData)
names(titData) 
datExpr0 = as.data.frame(t(titData[,-c(1)]))
names(datExpr0) =titData$Genes;
rownames(datExpr0) = names(titData)[-c(1)]
str(datExpr0)
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
str(gsg)
write.table(names(datExpr0)[!gsg$goodGenes], file="removeGene.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)  
write.table(names(datExpr0)[!gsg$goodSamples], file="removeSample.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)

sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(12,9)
pdf("SampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 15, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
keepSamples = (clust==0) 
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = read.csv("stage.txt",header = T,sep = "\t");
head(traitData)
dim(traitData)
names(traitData)
allTraits = traitData;
dim(allTraits)
names(allTraits)
titSamples = rownames(datExpr);
traitRows = match(titSamples, allTraits$Sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(12,9)
pdf("dendrogram and trait heatmap.pdf", width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, file = "tit-01-dataInput.RData")
lnames = load(file = "tit-01-dataInput.RData");
lnames

powers = c(c(1:10), seq(from = 11, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

softPower <- sft$powerEstimate
softPower

sizeGrWindow(9, 5)
pdf("Soft-thresholding powers.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(datExpr, power = softPower, maxBlockSize = 9000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "titTOM",
                       verbose = 3)
save(net, file = "tit-02-networkConstruction-net.RData")

# open a graphics window
sizeGrWindow(12, 9)
pdf("Dendrogram and the modules.pdf", width = 12, height = 9);
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "tit-02-networkConstruction-auto.RData")
####################################################################################################################################
lnames = load(file = "tit-02-networkConstruction-auto.RData");
lnames
#####################Relating modules to external clinical traits and identifying important genes
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
write.csv(moduleColors,"moduleColors.csv")
a<-table(moduleColors)
write.csv(a,"modulecolors_table.csv")


# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
write.csv(MEs,"MEs.csv")
moduleTraitCor = cor(MEs, datTraits, use = "p");
write.csv(moduleTraitCor,"cor-moduletrait.csv")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
write.csv(moduleTraitPvalue,"P-moduletrait.csv")

###Non-linear regression
me<-read.csv("MEs.csv",header = T,sep = ",")
stage<-read.csv("stage.txt",header = T,sep = "\t")

Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)}

pall<-c()
rall<-c()
RMSEall<-c()
for (i in 2:ncol(stage)) {
  x=stage[,i]
  p<-c()
  R2<-c()
  RMSE<-c()
  for (j in 2:ncol(me)) {
    y=me[,j]
    data<-data.frame(y=y,x=x)
    M<- lm(y ~ poly(x, 5, raw = TRUE),data=data)
    p[j]<-Regressionp(M)
    R2[j]<-summary(M)$adj.r.squared
    RMSE[j]<-sqrt(mean(M$residuals^2))
  }
  pall<-rbind(pall,p)
  rall<-rbind(rall,R2)
  RMSEall<-rbind(RMSEall,RMSE)
}
pall1<-t(as.data.frame(pall)[,-1])
rownames(pall1)<-names(me)[-1]
colnames(pall1)<-names(stage)[-c(1)]
rall1<-t(as.data.frame(rall)[,-1])
rownames(rall1)<-names(me)[-1]
colnames(rall1)<-names(stage)[-c(1)]
RMSEall1<-t(as.data.frame(RMSEall)[,-1])
rownames(RMSEall1)<-names(me)[-1]
colnames(RMSEall1)<-names(stage)[-c(1)]
write.table(pall1,"non-linear_pvalue.txt",quote = F)
write.table(rall1,"non-linear_cor.txt",quote = F)
write.table(RMSEall1,"non-linear_RMSEalla.txt",quote = F)

###Plot non-linear regression
RMSEall1<-round(RMSEall1, digits = 2)
rall1<-round(rall1, digits = 2)
pall1<-round(pall1, digits = 2)
textMatrix =  paste(rall1, "\n(",pall1, ")", sep = "");
dim(textMatrix) = dim(rall1)
par(mar = c(4, 6, 1, 1));
labeledHeatmap(Matrix = rall1,
               xLabels = names(stage)[-c(1)],
               yLabels = names(me)[-1],
               ySymbols = names(me)[-1],
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### Get the genes from target modules   
gene<-c()
for (i in c("turquoise","brown","red","grey","yellow","green")){
  which.module=i; 
  ME=MEs[, paste("ME",which.module, sep="")]
  par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
  plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
          nrgcols=30,rlabels=F,rcols=which.module,
          main=which.module, cex.main=2)
  par(mar=c(2,2.3,0.5,0.8))
  barplot(ME, col=which.module, main="", cex.main=2,
          ylab="eigengene expression",xlab="array sample")
  module = i;
  # Select module probes
  probes = colnames(datExpr) 
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  name<-paste0(i,".txt")
  write.table(modProbes,name,row.names = F,col.names = F,quote = F)
  gene<-c(gene,modProbes)
}
gene<-gene[!duplicated(gene)]
write.table(gene,"gene.txt",row.names = F,col.names = F,quote = F)

####
####Output the GS-value  and P-values of stage associated genes 
for(i in names(datTraits)){
  da = as.data.frame(datTraits[,i]);
  names(da) = i
  geneTraitSignificance = as.data.frame(cor(datExpr, da, use = "p"));
  GS<-abs(geneTraitSignificance)
  GSP = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  all<-cbind(GS,GSP)
  colnames(all)<-c("R","P")
  na<-paste0("GS_",i,"GSP.txt")
  write.table(all,na,col.names = T,quote = F,sep="\t")
}

####
####
#### Using pearson correlation to identify phenotypic associated genes

Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)}

  data<-data.frame(tpm=gene_tpm,trait=ind_trait)##matrix of expression value of single gene and corresponding trait value
  colnames(data)[c(1,2)]<-c("TPM","TRAIT")
  M<-cor.test(data$TPM,data$TRAIT,method = c("pearson"))
  p<-M$p.value ## P value of pearson correlation 
  r<-M$estimate ## Pearson correlation coefficient 


####
#### Candidate genes were filtered according to the thresholds used in the paper. 
#### 

#### Classify the candidate genes with different threshold 
setwd("F:\\PAPER\\bulk_rna_seq\\NH_treated\\Fgenes\\flight")
tpm<-read.csv("flight.tpm.csv",header = T,sep = ",")
adp<-read.table("candidate_genes.txt",header = T,sep = "\t")
dat<-merge(tpm,adp,by=c("Genes"))
dim(dat)
a<-c()
for (j in c(0.5,1,1.5,2)) {
  for (i in 1:nrow(dat)) {
    low<-median(t(dat[i,2:5]))
    hy<-median(t(dat[i,6:9]))
    high<-median(t(dat[i,10:13]))
    cut<-j*low
    {if (abs(hy-low) > cut & abs(high-hy) > cut){
      a[i]<-as.character(dat[i,1]) 
    } else
    {a[i]<-"NA"}}
    
  }
  name<-paste0(j,"_flight.txt")
  gene<-as.data.frame(a)
  gene<- gene[- grep("NA", gene$a),]
  gene<-as.data.frame(gene)
  colnames(gene)<-"Genes"
  write.table(gene,name,row.names = F,quote = F)}

####
####
####

#### Classify the genes with different threshold into groups of reinforcement or reversion plasticity 

tpm<-read.csv("flight.tpm.csv",header = T,sep = ",")
colnames(tpm)[1]<-"Genes"
dim(tpm)
file<-list.files(pattern = "*_flight.txt$",recursive = F) ## the files of genes with different threshold 
file
r<-c()
f<-c()
all<-c()
for (j in 1:length(file)) {
  gene<-read.table(file[j],header = T,sep = " ")
  dat<-merge(tpm,gene,by=c("Genes"))
  dat$Genes<-as.character(dat$Genes)
  dim(dat)
  all[j]<-nrow(gene)
  a<-c()
  a1<-c()
  b<-c()
  b1<-c()
  for (i in 1:nrow(dat)) {
    low<-median(t(dat[i,2:5]))
    hy<-median(t(dat[i,6:9]))
    high<-median(t(dat[i,10:13]))
    {if (hy > low & hy > high){
      b[i]<-as.character(dat[i,1]) 
    } else
    {b[i]<-"NA"}}
    
    {if (hy < low & hy < high){
      b1[i]<-as.character(dat[i,1])
    } else
    {b1[i]<-"NA"}}
    
    {if ( low < hy & hy < high){
      a[i]<-as.character(dat[i,1])
    } else
    {a[i]<-"NA"}}
    
    {if (low > hy & hy > high){
      a1[i]<-as.character(dat[i,1])
    } else
    {a1[i]<-"NA"}}
  }
  reinforce<-as.data.frame(a)
  reinforce1<-as.data.frame(a1)
  colnames(reinforce1)<-"a"
  reinforce<-rbind(reinforce,reinforce1)
  reinforce<- reinforce[- grep("NA", reinforce$a),]
  reinforce<-reinforce[!duplicated(reinforce)]
  reinforce<-as.data.frame(reinforce)
  f[j]<-nrow(reinforce)
  colnames(reinforce)<-"Genes"
  name<-paste0(file[j],"_reinforce.txt")
  write.table(reinforce,name,row.names = F,quote = F)
  
  reverse<-as.data.frame(b)
  reverse1<-as.data.frame(b1)
  colnames(reverse1)<-"b"
  reverse<-rbind(reverse,reverse1)
  reverse<- reverse[- grep("NA", reverse$b),]
  reverse<- reverse[!duplicated(reverse)]
  reverse<-as.data.frame(reverse)
  colnames(reverse)<-"Genes"
  r[j]<-nrow(reverse)
  name1<-paste0(file[j],"_reverse.txt")
  write.table(reverse,name1,row.names = F,quote = F)
}
rf<-data.frame(reinforce=f,reverse=r,all=all)
write.table(rf,"reinforce and reversion.txt",row.names = F,quote = F,sep = "\t")

####Binomial test
data<-read.csv("reinforce and reversion.txt",header = T,sep = "\t")
head(data)
p<-c()
for (i in 1:nrow(data)) {
  test<-binom.test(x=data[i,2],n=data[i,3],p=data[i,1]/data[i,3],alternative = "two.side")
  p[i]<-test$p.value
}
p1<-cbind(data,p)
write.table(p1,"reinforce and reversion_Pvalues.txt",row.names = F,quote = F,sep = "\t")

### Bootstrap method

se<- function(x){
  sd(x)/sqrt(length(x))
}

file<-list.files(pattern = "*_flight.dependent.txt$",recursive = F)  ## reinforcement and reversion geneset 
file
all<-c()
for (i in 1:length(file)) {
  data<-read.table(file[i],header = T,sep = "\t")
  data$group<- rep(file[i],nrow(data))
  all<-rbind(all,data)
}

data<-read.csv("F:\\PAPER\\bulk_rna_seq\\NH_treated\\1118\\flight.tpm.csv",header = T,sep = ",")
colnames(data)[1]<- "Genes"
dat<- left_join(all,data,by=c("Genes"))

set.seed(123)
num_reinforce<-c()
num_reversion<-c()
for (i in 1:nrow(dat)) {
  me_low<- rowMeans(dat[i,3:6])
  me_hy<- rowMeans(dat[i,7:10])
  me_high<- rowMeans(dat[i,11:14])
  
  sd_low<- se(dat[i,3:6])
  sd_hy<- se (dat[i,7:10])
  sd_high<- se (dat[i,11:14])
  
  simu_low<- rnorm(1000, mean = me_low, sd = sd_low) 
  simu_hy<- rnorm(1000, mean = me_hy, sd = sd_hy) 
  simu_high<- rnorm(1000, mean = me_high, sd = sd_high) 
  
  simu_new<- as.data.frame(cbind(simu_low,simu_hy,simu_high))  
  simu_new$group_reinforce <- ifelse((simu_new$simu_high - simu_new$simu_hy > 0 & simu_new$simu_hy - simu_new$simu_low >0) | 
                             (simu_new$simu_high - simu_new$simu_hy < 0 & simu_new$simu_hy - simu_new$simu_low < 0),"YES","NO")
  num_reinforce[i]<- sum(simu_new$group_reinforce=="YES")
  
  simu_new$group_reversion <- ifelse((simu_new$simu_high - simu_new$simu_hy > 0 & simu_new$simu_hy - simu_new$simu_low <0) | 
                             (simu_new$simu_high - simu_new$simu_hy < 0 & simu_new$simu_hy - simu_new$simu_low > 0),"YES","NO")
  num_reversion[i]<- sum(simu_new$group_reversion=="YES")
}
res_reinforce<-data.frame(Genes=dat$Genes,group=dat$group,numYES=num_reinforce)
res_reversion<-data.frame(Genes=dat$Genes,group=dat$group,numYES=num_reversion)
write.table(res,"flight_allgroups_reinforce_gennes.txt",row.names = F,sep = "\t",quote = F)
write.table(res,"flight_allgroups_reversion_gennes.txt",row.names = F,sep = "\t",quote = F)

reinforce_geneNum<-nrow(subset(res_reinforce,res$numYES>=950))
reversion_geneNum<-nrow(subset(res_reversion,res$numYES>=950))

#### Random sampling for permuted FST distributions 
# FST distributions for the genes within each category by random sampling the same number of genes, the same number of SNPs(allowing for a fluctuation of 5% of total SNPs)

gff<-read.csv("Shuma.genome.gff_gb2k.txt",header=F,sep="\t") ## annotation file of genic region, the columns includes "chrom","startpos","endpos","GeneID","SNP_number"
head(gff)
list<-read.csv("genelist_all.txt",header = F,sep = "\t") # all annotated genes with specific GeneID

f<-list.files(pattern = "*_flight.dependent.txt$",recursive = F)  ## reinforcement and reversion geneset 
f

for (j in 1:length(f)) {
  g<-paste0("F:\\PAPER\\bulk_rna_seq\\NH_treated\\1118\\D\\genic_region\\",f[j])
  gn<-read.table(g,header = T,sep = "\t")
  colnames(gn)<-"V4"
  num<-nrow(gn)
  gnm<-merge(gn,gff,by=c("V4"))
  gf<-paste0("F:\\PAPER\\bulk_rna_seq\\NH_treated\\1118\\D\\genic_region\\",substr(f[j],1,19))
  dir.create(gf)
  setwd(gf)
  while (length(list.files())<=99) {
    for (i in 1) {
      gene <- sample(list$V1,num,replace = F)
      gene<-as.data.frame(gene)
      colnames(gene)[1]<-"V4"
      gene1<-merge(gene,gff,by=c("V4"))
      if (sum(gene1$V5) < sum(gnm$V5)+sum(gnm$V5)*0.05 & sum(gene1$V5) > sum(gnm$V5)-sum(gnm$V5)*0.05){
        sample<-gene1[,c(2,3,4)]
        colnames(sample)<-c("CHROM","start","end")
        name<-paste0(sample(c(1:1000000), 1, replace = TRUE),"_",i,"_g.pos")
        write.table(sample,name,row.names=F,quote=F,sep = "\t")}}
  }
}


####
####

#### Calculate average FST using Vcftools

#### Calculate observed fst 
#! bin bash
setwd("/home/shehs/sparrow_accli_reseq/genic_region")
cd ../genic_region
for line in `ls *g.txt`;
do
dos2unix ${line}
vcftools --vcf tree_sparrow_analysis_ready.vcf --out ${line}.vcf --recode --bed ${line}
vcftools --vcf ${line}.vcf.recode.vcf --weir-fst-pop highland.txt \
--weir-fst-pop lowland.txt \
--out ${line}.fst 
rm ${line}.vcf.recode.vcf 
done

#### Calculate permuted fst 
#! bin bash
cat ../genic_region/file_list.txt | while read a
do
cd ../genic_region/${a}
for line in `ls *.pos`;
do
dos2unix ${line}
vcftools --vcf tree_sparrow_analysis_ready.vcf --out ${line}.vcf --recode --bed ${line}
vcftools --vcf ${line}.vcf.recode.vcf --weir-fst-pop highland.txt \
--weir-fst-pop lowland.txt \
--out ${line}.fst 
rm ${line}.vcf.recode.vcf 
done
done
####
