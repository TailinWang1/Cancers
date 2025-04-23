library(readxl)
library(dplyr)
library(limma)
library(GEOquery)
library(impute)
library(readr)
library(preprocessCore)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(WGCNA)
#DGE analysis
PTC_top_table<- read.delim("GSE.top.table.tsv", header = TRUE, sep = "\t")
PTC_top_table$diffexpressed <- "NS"
PTC_top_table$diffexpressed[PTC_top_table$adj.P.Val< 0.05 & PTC_top_table$logFC > 1] <- "Up"
PTC_top_table$diffexpressed[PTC_top_table$adj.P.Val< 0.05 & PTC_top_table$logFC< -1] <- "Down"
top_genes <- rbind(
  PTC_top_table %>% 
    filter(diffexpressed == "Up") %>% 
    top_n(10, wt = logFC),
  PTC_top_table %>% 
    filter(diffexpressed == "Down") %>% 
    top_n(10, wt = -logFC)
)
volcano_plot <- ggplot(data=PTC_top_table, 
                       aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) +
  geom_point(size=1.5, alpha=0.7) +
  scale_color_manual(values=c("Down"="#4575b4", "NS"="#cccccc", "Up"="#d73027")) +
  geom_vline(xintercept=c(-1, 1), col="grey50", linetype="dashed", size=0.5) +
  geom_hline(yintercept=-log10(0.05), col="grey50", linetype="dashed", size=0.5) +
  geom_text_repel(data = top_genes, aes(label = Gene.symbol), 
                  size = 3, box.padding = 0.5, point.padding = 0.3, 
                  segment.color = "grey50", show.legend = FALSE) +
  labs(title="Differential Gene Expression in PTC",
       x="log2 Fold Change",
       y="-log10 Adjusted P-value",
       color="Expression") +
  theme_minimal() +
  theme(legend.position="right",
        plot.title = element_text(size=16, face="bold", hjust=0.5),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        panel.grid.major = element_line(colour = "grey90"),
        panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 12))  # 根据数据范围调整
print(volcano_plot)
PTC_DEGs <- PTC_top_table[which(abs(PTC_top_table$logFC) > 1 & PTC_top_table$adj.P.Val < 0.05), c(2, 7)]
str(PTC_DEGs)
DEG31 <- unique(PTC_DEGs$Gene.symbol[!is.na(PTC_DEGs$Gene.symbol) & PTC_DEGs$Gene.symbol != ""])
write.csv(DEG31, 
          file = "PTC_DEGs_155个.csv", 
          row.names = FALSE,
          quote = FALSE)
#WGCNA
T1D_gset <- getGEO('GSE', destdir=".",AnnotGPL = T,getGPL = T)
T1D_exp<-exprs(T1D_gset[[1]])
T1D_GPL<-fData(T1D_gset[[1]])
T1D_gpl<- T1D_GPL[, c(1, 3)]
T1D_exp<-as.data.frame(T1D_exp)
T1D_exp$ID<-rownames(T1D_exp)
T1D_exp_symbol<-merge(T1D_exp,T1D_gpl,by="ID")
T1D_exp_symbol<-na.omit(T1D_exp_symbol)
table(duplicated(T1D_exp_symbol$`Gene symbol`))
T1D_datExpr02<-avereps(T1D_exp_symbol[,-c(1,ncol(T1D_exp_symbol))],ID=T1D_exp_symbol$`Gene symbol`)
datExpr0<- t(T1D_datExpr02)
m.vars = apply(datExpr0, 2, var)
cutoff = quantile(m.vars, prob = 0.8)
expro.upper = datExpr0[, which(m.vars > cutoff)]
datExpr1 <- data.matrix(expro.upper)
gsg = goodSamplesGenes(datExpr1, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr1), method = "average")
par(cex = 0.7);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
datExpr = as.data.frame(datExpr1)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("1PTC02Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft$powerEstimate
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       #saveTOMFileBase = "MyTOM",
                       verbose = 3)
table(net$colors)
mergedColors = labels2colors(net$colors)
pdf("2PTC02module.pdf",width = 12, height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
text <- unique(moduleColors)
for (i  in 1:length(text)) {
  y=t(assign(paste(text[i],"expr",sep = "."),datExpr[moduleColors==text[i]]))
  write.csv(y,paste(text[i],"csv",sep = "."),quote = F)
}
T1Dsamples <- read_excel('PTC临床信息全部18.xlsx')
row_names <- T1Dsamples[[1]]
T1Dsamples <- T1Dsamples %>% select(-1)
T1Dsamples_df <- as.data.frame(T1Dsamples)
row.names(T1Dsamples_df) <- row_names
print(T1Dsamples_df)
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
dim(MEsWW)  
dim(T1Dsamples_df) 
modTraitCor = cor(MEsWW, T1Dsamples_df, use = "p")
colnames(MEsWW)
modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf("3PTC02Module-trait.pdf",width = 6, height = 12)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(T1Dsamples_df), yLabels = names(MEsWW), cex.lab = 0.8,  yColorWidth=0.02, 
               xColorWidth = 0.04,
               ySymbols = colnames(modlues), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1)
               , main = paste("Module-trait relationships"))
dev.off()
modTraitCor_df <- as.data.frame(modTraitCor)
modTraitP_df <- as.data.frame(modTraitP)
max_cor <- apply(abs(modTraitCor_df), 1, max)
significant_modules <- data.frame(
  Module = names(max_cor),
  MaxAbsCor = max_cor,
  Trait = colnames(modTraitCor_df)[apply(abs(modTraitCor_df), 1, which.max)]
)
significant_modules <- significant_modules[order(-significant_modules$MaxAbsCor), ]
min_p <- apply(modTraitP_df, 1, min)
significant_modules_p <- data.frame(
  Module = names(min_p),
  MinPvalue = min_p,
  Trait = colnames(modTraitP_df)[apply(modTraitP_df, 1, which.min)]
)
significant_modules_p <- significant_modules_p[order(significant_modules_p$MinPvalue), ]
print("按相关系数排序的模块：")
print(head(significant_modules))
print("\n按P值排序的模块：")
print(head(significant_modules_p))
magentaGenes <- colnames(magenta.expr)
write.csv(magentaGenes, "PTC_WGCNA最相关基因.csv", row.names = FALSE, quote = FALSE)
#Combine DEG and WGCNA
kGenes31<- intersect(DEG31,magentaGenes)
key_genes_df <- data.frame(Gene_Symbol = kGenes31)
write.csv(key_genes_df, "PTC_key_genes_52个.csv", row.names = FALSE)
