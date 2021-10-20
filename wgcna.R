## https://www.jianshu.com/p/f0409a045dab

setwd('/users/leump/Desktop/')

lncRNAexpr <- read.csv("All_sample_LncRNA_exp_RPKM.csv",sep=",",header = T)

head(lncRNAexpr)

dim(lncRNAexpr)

##去掉Annotation这列
fpkm <- lncRNAexpr[,-66]
head(fpkm)
##重命名行名和列名
rownames(fpkm)=fpkm[,1]
fpkm=fpkm[,-1]
fpkm[1:4,1:4]

##Sample Info
subname=sapply(colnames(fpkm),function(x) strsplit(x,"_")[[1]][1])
datTraits = data.frame(gsm=names(fpkm),
                       subtype=subname)
rownames(datTraits)=datTraits[,1]
head(datTraits)

library(WGCNA)

RNAseq_voom <- fpkm 
WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
datExpr <- WGCNA_matrix  ## top 5000 mad genes
datExpr[1:4,1:4]

sampleNames = rownames(datExpr);
traitRows = match(sampleNames, datTraits$gsm)  
rownames(datTraits) = datTraits[traitRows, 1]
datTraits


net = blockwiseModules(datExpr, power = 3,
                       maxBlockSize = 6000,TOMType = "unsigned", 
                       minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)

table(net$colors)

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.

#明确样本数和基因数
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#首先针对样本做个系统聚类树
datExpr_tree<-hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1, cex.main = 1,cex.lab=1)
## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
#针对前面构造的样品矩阵添加对应颜色
sample_colors <- numbers2colors(as.numeric(factor(datTraits$subtype)), 
                                colors = c("grey","blue","red","green"),signed = FALSE)
## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目。
#  sample_colors <- numbers2colors( datTraits ,signed = FALSE)
## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，
##当然，这样给的颜色不明显，意义不大。
#构造10个样品的系统聚类树及性状热图
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors,
                    groupLabels = colnames(sample),
                    cex.dendroLabels = 0.8,
                    marAll = c(1, 4, 3, 1),
                    cex.rowText = 0.01,
                    main = "Sample dendrogram and trait heatmap")


design=model.matrix(~0+ datTraits$subtype)

colnames(design)=levels(factor(datTraits$subtype))
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));



# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

head(design)


# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="")

## 只有连续型性状才能只有计算
## 这里把是否属于 CB 表型这个变量用0,1进行数值化。
CB = as.data.frame(design[,2]);
names(CB) = "CB"
geneTraitSignificance = as.data.frame(cor(datExpr, CB, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CB), sep="");
names(GSPvalue) = paste("p.GS.", names(CB), sep="")

module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for CB",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# STEP7:网络的可视化
#首先针对所有基因画热图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 3); # 设置power=sft$powerEstimate，最佳beta值，此处是3
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#然后随机选取部分基因作图
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")



#最后画模块和性状的关系
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属于 Luminal 表型这个变量用0,1进行数值化。
CB = as.data.frame(design[,2]);
names(CB) = "CB"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, CB))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
## 模块的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



# STEP8:提取指定模块的基因名
# Select module
module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]

































