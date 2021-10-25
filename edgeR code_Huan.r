#####read in counts and traits data
library(edgeR)
alldesign=read.csv("colData_final.csv")
head(alldesign)
allcount=read.csv("20180208peter_batch1-2-3_allcounts.csv")
counts=allcount[,-c(1:4)]
row.names(allcount)=allcount[,3]
min(colSums(counts))  #460708   #0.460708   #5069308 in batch2------1142903
cpms = cpm(counts)
keep = rowSums(cpms >1) >=4   ####it requires 5~10 counts in one library to indicate the gene is expressed. the smallest total counts sample contain A counts, A/1000000 means ONE cpm has A/1000000 counts, and "1" indicates 5/(A/1000000), and "3" means the least samples in one group, like here there are at least 3 samples in one condition. Hence we keep genes with at least two counts per million (CPM) in at least four samples:
counts=allcount[keep,-c(1:4)]
dim(counts)   #[1] 16810  
names(counts)
#colnames(counts)[1:ncol(counts)]=alldesign$DiseaseActivity
head(counts)
##############################**********************************************edger begin

#Create a design matrix (see 'Experimental design' for further details) to specify the factors that are expected to affect expression levels:
d = DGEList(counts=counts,genes=row.names(counts),group=alldesign$DiseaseActivity)

###########delete duplicate gene, while keep the genes left with the largest rowSums of counts 
nrow(d)
o <- order(rowSums(d$counts), decreasing=TRUE)
d <- d[o,]
duplicate <- duplicated(d$genes)
d <- d[!duplicate,]
nrow(d)

# d$samples$lib.size <- colSums(d$counts)
rownames(d$counts) <- rownames(d$genes) <- d$genes$genes

d = calcNormFactors(d)

#design = model.matrix(~DiseaseActivity,data=samples)
design = model.matrix(~0+DiseaseActivity,data=alldesign)
colnames(design)
#[1] "DiseaseActivityCD_Ac" "DiseaseActivityCD_In" "DiseaseActivityHC_In" "DiseaseActivityUC_Ac" "DiseaseActivityUC_In"

#plotMDS(d)
#Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood7, 53, as follows:
d2 = estimateGLMTrendedDisp(d, design)
d2 = estimateGLMTagwiseDisp(d2, design)

#Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)

#Perform a likelihood ratio test, specifying the difference of interest (here, knockdown versus control, which corresponds to the third column of the above design matrix):
#"CD_Ac" "CD_In" "HC_In" "UC_Ac" "UC_In"
de1 = glmLRT(f, contrast=c(1,-1,0,0,0))   #####CD active / CD inactive    
de2 = glmLRT(f, contrast=c(0,0,0,1,-1))   #####UC active / UC inactive    
de3 = glmLRT(f, contrast=c(1,0,0,-1,0))   #####active CD vs active UC     
de4 = glmLRT(f, contrast=c(0,-1,0,0,1))   #####inactive UC vs inactive CD     
de5 = glmLRT(f, contrast=c(1,0,-1,0,0))   #####active CD vs controls    
de6 = glmLRT(f, contrast=c(0,0,-1,1,0))   #####active UC versus controls    
de7 = glmLRT(f, contrast=c(0,0,-1,0,1))   #####inactive UC versus controls 
de8 = glmLRT(f, contrast=c(0,1,-1,0,0))   #####inactive CD versus controls

sigGenes=unique(c(subset(topTags(de1, n=nrow(de1))$table,FDR<=0.1)[,1],
                  subset(topTags(de2, n=nrow(de2))$table,FDR<=0.1)[,1],
                  subset(topTags(de3, n=nrow(de3))$table,FDR<=0.1)[,1],
                  subset(topTags(de4, n=nrow(de4))$table,FDR<=0.1)[,1],
                  subset(topTags(de5, n=nrow(de5))$table,FDR<=0.1)[,1],
                  subset(topTags(de6, n=nrow(de6))$table,FDR<=0.1)[,1],
                  subset(topTags(de7, n=nrow(de7))$table,FDR<=0.1)[,1],
                  subset(topTags(de8, n=nrow(de8))$table,FDR<=0.1)[,1]))
length(sigGenes)  ###11420

NameTrans=function(de,NAME){
  tt = topTags(de, n=nrow(de))
  tt=tt$table
  tt=tt[c(1,2,5,6)]
  print(c("siggene number",nrow(subset(tt,FDR<=0.1))))
  colnames(tt)[2:4]=c(paste(NAME,"logFC",sep=":"),paste(NAME,"PValue",sep=":"),paste(NAME,"FDR",sep=":"))
  return(tt)
}
de1=NameTrans(de1,"CD_ac|CD_in")  ##"69----372
de2=NameTrans(de2,"UC_ac|UC_in")  ##"2608"-----7734
de3=NameTrans(de3,"CD_ac|UC_ac")   ###"1356" ------3472
de4=NameTrans(de4,"UC_in|CD_in")  ##"51"------137
de5=NameTrans(de5,"CD_ac|control")  ###"346"-----2510
de6=NameTrans(de6,"UC_ac|control")  ### "5485" -------9723
de7=NameTrans(de7,"UC_in|control")   ###"176" -------1345
de8=NameTrans(de8,"CD_in|control")  ##"57" -----608

results=cbind(de1,de2[,-1],de3[,-1],de4[,-1],de5[,-1],de6[,-1],de7[,-1],de8[,-1])
results=cbind(results,counts[,order(alldesign$DiseaseActivity)])


############## GO analysis 

  ######GO
  ego2 <- enrichGO(gene         = as.character(inputGENElist),
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'ENSEMBL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   #pvalueCutoff  = 0.01,
                   #qvalueCutoff  = 0.05)
                   pvalueCutoff  = 1,
                   qvalueCutoff  = 0.1)
  png(filename=paste(OutputName,".GO_enrichTERM_top10.png",sep=""),width = 8,height = 4,units = "in",res=600 )
  print(barplot(ego2 , drop=TRUE, showCategory=10))
  dev.off()
  