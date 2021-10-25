################### main function

#### do not use the function below in my own computer, for some multi-thread options are on!!!!!
WGCNA=function(DATATraits,DATARpkm,OUTNAME){
  
  #=====================================================================================
  #
  #  Code chunk 1
  #
  #=====================================================================================
  
  # Load the WGCNA package
  library(WGCNA);
  # The following setting is important, do not omit.
  options(stringsAsFactors = FALSE);
  #Read in the female liver data set
  femData=DATARpkm
  
  #=====================================================================================
  #
  #  Code chunk 2
  #
  #=====================================================================================
  
  
  datExpr0 = as.data.frame(t(femData));
  
  
  
  #=====================================================================================
  #
  #  Code chunk 3
  #
  #=====================================================================================
  
  
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  
  #=====================================================================================
  #
  #  Code chunk 4
  #
  #=====================================================================================
  
  
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
  
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  
  
  sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  #sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  png(paste(OUTNAME,".sampleClustering.png",sep=""),width = 15,heigh=8,units="in",res=300)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
       cex.axis = 1.5, cex.main = 2)
  
  
  #=====================================================================================
  #
  #  Code chunk 6
  #
  #=====================================================================================
  
  
  # Plot a line to show the cut
  abline(h = 15, col = "red");
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 1)
  dev.off()
  table(clust)
  # clust 1 contains the samples we want to keep.
  #keepSamples = (clust==1)
  #datExpr = datExpr0[keepSamples, ]
  #nGenes = ncol(datExpr)
  #nSamples = nrow(datExpr)
  datExpr=datExpr0
  
  #=====================================================================================
  #
  #  Code chunk 7
  #
  #=====================================================================================
  
  
  # Form a data frame analogous to expression data that will hold the clinical traits.
  
  femaleSamples = rownames(datExpr);
  #traitRows = match(femaleSamples, row.names(allTraits));
  DATATraits =DATATraits 
  indx <- sapply(DATATraits, is.factor)
  DATATraits[indx] <- lapply(DATATraits[indx], function(x) as.numeric(x))
  collectGarbage();
  
  
  #=====================================================================================
  #
  #  Code chunk 8
  #
  #=====================================================================================
  
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  #traitColors = numbers2colors(DATATraits, signed = FALSE);
  traitColors = numbers2colors(DATATraits, signed = TRUE);
  # Plot the sample dendrogram and the colors underneath.
  png(paste(OUTNAME,".dendrogram.png",sep=""),width = 15,heigh=8,units="in",res=300)
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = colnames(DATATraits), 
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  
  #=====================================================================================
  #
  #  Code chunk 9
  #
  #=====================================================================================
  
  
  # Choose a set of soft-thresholding powers
  # powers = c(c(1:10), seq(from = 12, to=20, by=2))
  powers =  c(c(1:10), seq(from = 12, to=30, by=2))
  
  #dim(datExpr[,colSums(datExpr)>100])
  #datExpr.sub=datExpr[,colSums(datExpr)>100]
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  png(paste(OUTNAME,".power.",sft$powerEstimate,".png",sep=""),width = 15,heigh=8,units="in",res=300)
  #sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  power = sft$powerEstimate
  dev.off()
  print(power)
  #power=9
  #power=18
  #=====================================================================================
  #
  #  Code chunk 3
  #
  #=====================================================================================
  cor = WGCNA::cor
  #####-----power is important, and type is unsigned for default. But we need to judge the correlation is positive or negative, so we should use unsigned.
  net = blockwiseModules(datExpr, power = power,
                         #TOMType = "unsigned", 
                         TOMType = "signed", 
                         minModuleSize = 30,
                         #reassignThreshold = 0, 
                         reassignThreshold = 0.7, 
                         mergeCutHeight = 0.35, #Maximum dissimilarity (i.e., 1-correlation) that qualifies modules for merging.
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         #saveTOMs = TRUE,
                         #maxBlockSize = 20000, 
                         #saveTOMFileBase = "femaleMouseTOM", 
                         verbose = 3)
  cor = stats::cor
  
  #=====================================================================================
  #
  #  Code chunk 4
  #
  #=====================================================================================
  
  
  # open a graphics window
  #sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Plot the dendrogram and the module colors underneath
  png(paste(OUTNAME,".clustereddendrogram.png",sep=""),width = 15,height=8,units="in",res=300)
  print(plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05))
  dev.off()
  
  #=====================================================================================
  #
  #  Code chunk 5
  #
  #=====================================================================================
  
  
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  #save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-auto.RData")
  
  
  
  
  
  #=====================================================================================
  #
  #  Code chunk 1
  #
  #=====================================================================================
  
  
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, DATATraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  
  #=====================================================================================
  #
  #  Code chunk 3
  #
  #=====================================================================================
  
  
  #sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  #png("test3.png",width=16,heigh=28,units="in",res=600)
  #png(paste(OUTNAME,".WGCNA.png",sep=""),width=30,height=30,units="in",res=600)
  png(paste(OUTNAME,".WGCNA.png",sep=""),width=10,height=15,units="in",res=600)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(DATATraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors=blueWhiteRed(50),
                 #colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  
  dev.off()
  
  save(list = ls(all.names = TRUE),file="my_work_WGCNA.RData")
  hub=chooseTopHubInEachModule(datExpr = datExpr,colorh =moduleColors,
                               #omitColors = "grey",
                               power=power,
                               type="signed")
  #For modular membership:
  
  MM = as.data.frame(cor(datExpr, MEs, use = "p"));
  modNames = substring(names(MEs), 3)
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(MM), nSamples));
  #names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  #intramodular connectivity
  connet=abs(cor(datExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  
  
  cluster.colors=function(ll,OUTNAME){
    tmmp=NULL
    hubGene=NULL
    membership=NULL
    for (i in 1:length(ll)){
      tmmp=unique(c(tmmp,names(datExpr)[moduleColors==ll[i]]))
      
      hubGene=c(hubGene,hub[names(hub)==ll[i]])
      
      coll=paste("ME",ll[i],sep="")
      
      #DIR=paste(coll,"_",length(tmmp),sep="")
      #dir.create(DIR)
      membership=unique(c(membership,row.names(MM[order(MM[[coll]],decreasing = T),][1:10,])))
      
      GeneID=sapply(tmmp,function(x){unlist(strsplit(x,split = "_"))[1]})
      module=ll[1]
      moduleGenes = moduleColors==module
      membership_file=MM[moduleGenes,]
      MMPvalue_file=MMPvalue[moduleGenes,coll]
      geneTraitSignificance = as.data.frame(cor(datExpr, DATATraits, use = "p"));
      GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
      connectivity=Alldegrees1$kWithin[moduleGenes]
      names(GSPvalue) = paste("p.GTS", names( GSPvalue), sep="");
      outfile=data.frame(Gene=row.names(membership_file),geneModuleMembership=membership_file[,coll],
                         MMPvalue=MMPvalue_file,geneTraitSignificance[moduleGenes,],GSPvalue[moduleGenes,],
                         connectivity=connectivity)
      write.csv(outfile,paste(OUTNAME,"_",paste(ll,collapse = "-"),".allResults.csv",sep=""),row.names = F,quote=F)
      write.csv(outfile[outfile$Gene%in%membership,],paste(OUTNAME,"_",paste(ll,collapse = "-"),".topResults.csv",sep=""),row.names = F,quote=F)
      
      
      ###Finding genes with high gene significance and high intramodular connectivity in
      # interesting modules
      # abs(GS1)> .9 
      # abs(datKME$MM.black)>.8 
      #FilterGenes= abs(GS1)> .9 & abs(datKME$MM.black)>.8
      
      #（3） Diagnostics: displaying module heatmap and the eigengene
      png(paste(OUTNAME,".",".ModuleHeatmap_eigengene.png",sep=""),width=10,height=10,units="in",res=600)
      #sizeGrWindow(8,7);
      ME=MEs[, paste("ME",module, sep="")]
      par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
      plotMat(t(scale(datExpr[,moduleColors==module ]) ),
              nrgcols=30,rlabels=F,rcols=module,
              main=module, cex.main=2)
      par(mar=c(5, 4.2, 0, 0.7))
      barplot(ME, col=module, main="", cex.main=2,
              ylab="eigengene expression",xlab="Samples")
      dev.off()
      
      
      
      
      for ( j in 1:length(DATATraits)){
        CB = as.data.frame(DATATraits[,j]);
        names(CB)=colnames(DATATraits)[j]
        geneTraitSignificance = as.data.frame(cor(datExpr, CB, use = "p"));
        GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
        names(geneTraitSignificance) = paste("GS.", names(CB), sep="");
        names(GSPvalue) = paste("p.GS.", names(CB), sep="")
        #module = "turquoise"
        column = match(module, modNames);
        moduleGenes = moduleColors==module;
        #sizeGrWindow(7, 7);
        par(mfrow = c(1,1));
        png(paste(OUTNAME,".",names(DATATraits)[j],".",module,".module.png",sep=""),width=10,height=10,units="in",res=600)
        verboseScatterplot(abs(MM[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for",names(DATATraits)[i]),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        dev.off()
        
        # Add the weight to existing module eigengenes
        MET = orderMEs(cbind(MEs, CB))
        # Plot the relationships among the eigengenes and the trait
        png(paste(OUTNAME,".",names(DATATraits)[j],".",module,".eigengenes_trait.png",sep=""),width=10,height=10,units="in",res=600)
        #sizeGrWindow(5,7.5);
        par(cex = 0.9)
        plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                              = 90)
        
        dev.off()
        
      }
      
    }
    print(length(tmmp))
    write.csv(tmmp,paste(OUTNAME,"_",paste(ll,collapse = "-"),".allGene.csv",sep=""),row.names = F,quote=F)
    GeneID=sapply(tmmp,function(x){unlist(strsplit(x,split = "_"))[1]})
    write.csv(hubGene,paste(OUTNAME,"_",paste(ll,collapse = "-"),".hub.csv",sep=""),row.names = F,quote=F)
    write.csv(membership,paste(OUTNAME,"_",paste(ll,collapse = "-"),".topMembership.csv",sep=""),row.names = F,quote=F)
  }
  
  
  for (i in 1:length(colnames(MEs))){
    MEs.color=colnames(MEs)[i]
    dir.create(MEs.color)
    cluster.colors(unlist(strsplit(MEs.color,split="ME"))[2],
                   paste(MEs.color,"/",unlist(strsplit(MEs.color,split="ME"))[2],sep=""))
  }
  
  
}

###################read in traits and counts data
#setwd("~/Desktop/uab/dkcross2/")
#library(edgeR)
alldesign=read.csv("colData_final.csv")
head(alldesign)
allcount=read.csv("20180208peter_batch1-2-3_allcounts.csv")
counts=allcount[,-c(1:4)]
row.names(allcount)=allcount[,1]
min(colSums(counts))  #460708   #0.460708   #5069308 in batch2------1142903
#cpms = cpm(counts)
#keep = rowSums(cpms >1) >=4   ####it requires 5~10 counts in one library to indicate the gene is expressed. the smallest total counts sample contain A counts, A/1000000 means ONE cpm has A/1000000 counts, and "1" indicates 5/(A/1000000), and "3" means the least samples in one group, like here there are at least 3 samples in one condition. Hence we keep genes with at least two counts per million (CPM) in at least four samples:
allcount=allcount[,-c(1,3:4)]
tpm3 <- function(counts,len) {
  x <- counts/len
  x=x[rowSums(x)!="NaN",]
  return(t(t(x)*1e6/colSums(x,na.rm = T)))
}
tpm=tpm3(allcount[,-1],allcount[,1])

dim(tpm)   #[1] 54135  89
#head(tpm,1)
colnames(tpm)[1:ncol(tpm)]=as.character(alldesign$RNASeq.code.)

tpm=tpm[,alldesign$batch!="batch1"]

alldesign=alldesign[alldesign$batch!="batch1",]


######### read in comparison file
test=read.delim("Huan_Trait_data_for_AcCD_for_batch_2_3.txt")
head(test,1)
############
tpm=tpm[,colnames(tpm)%in%as.character(test$Sample.code)]
test=test[test$Sample.code%in%colnames(tpm),]
head(test)
head(tpm)
dim(test)
dim(tpm)
tpm=tpm[,order(colnames(tpm))]
test=test[order(test$Sample.code),]
colnames(tpm)
test$Sample.code
names(test)
Data.num=test[,-c(1:5,7)]
row.names(Data.num)=test[,2]

rpkm=as.matrix(tpm)
traitsOrg=as.matrix(Data.num)
traitsOrg=apply(traitsOrg[,1:ncol(traitsOrg)],2,function(x){as.numeric(x)})

row.names(Data.num)
colnames(rpkm)
###select traits
colnames(traitsOrg)
##IL-13, SOCS1, IL-17, ITAC, endo score, IFNg, IL-8, MIP1a, IL-23, CRP, IL-12

#cbind(alldesign[alldesign$Sample.code.%in%test$RNA.Seq.code,c(2,6)],test[,c(1,20)])
dir.create("2020March.addIL10.addBio.v1.details")
setwd("2020March.addIL10.addBio.v1.details")

WGCNA(as.data.frame(traitsOrg[,c(10,8,11,17,4,9,34,35,30,5,23,21,14)]),rpkm,"addIL10.addBio")####10



