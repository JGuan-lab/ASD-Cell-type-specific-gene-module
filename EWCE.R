install.packages("devtools")
library(devtools)
install_github("nathanskene/ewce",force = TRUE)

options(repos="https://CRAN.R-project.org")
install.packages("usethis")
library("usethis")
install.packages("ggdendro")

library("scater")
library("SingleCellExperiment")
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
getwd()


setwd("D:/Yiping/allen_MTG_sce/R code")
sce <- readRDS("allen_MTG_sce_protein_coding.rds")  #17120 12506
y<-read.csv("SFARI-Genes.csv")# asd candidate genes 1054
example_genelist=y1$gene.symbol

exp_combat1=sce@assays$data$combat_data
symbol <- sce@rowRanges@elementMetadata$gene
rownames(exp_combat1)<-symbol
ENTREZID<-rownames(sce@assays$data$combat_data)
trans <- data.frame(symbol,ENTREZID)
hvg <- readRDS("hvg_protein_coding.rds")     #7011
hvg_genes <-rownames(hvg)
ix = match(hvg_genes,trans[,2])
exp_combat <- exp_combat1[ix,]
cell_cluster <- sce@colData$cluster
unique_cluster=unique(cell_cluster)
cell_type = colData(sce)$Cell_type
unique_cell = unique(cell_type)
#### Loading datasets
expData = as.data.frame(exp_combat) 
l1 = colData(sce)$Cell_type   #cell type
l2 = cell_cluster     #sub cell type 
annotLevels = list(l1=l1,l2=l2)

setwd("D:/Yiping/allen_MTG_sce/allen_MTG_sce(2)/review/EWCE")
fNames_ALLCELLS = generate.celltype.data(exp=expData,annotLevels,"allKImouse",no_cores=1)
load(fNames_ALLCELLS)

#############minfc
setwd("D:/Yiping/allen_MTG_sce/R code/find_specificity/cluster_minfc")
cluster_minfc = matrix(nrow=nrow(expData),ncol=length(unique_cluster))
genes = rownames(expData)
for(i in 1:length(unique_cluster)){
  data = read.csv(paste(unique_cluster[i],"_cluster_minfc.csv",sep=""))
  fc = data[order(data[,1]),"cluster_min_FC"]
  cluster_minfc[,i] = fc
}
rownames(cluster_minfc) = genes
colnames(cluster_minfc) = unique_cluster

cell_minfc = matrix(nrow=nrow(expData),ncol=length(unique_cell))
for(i in 1:length(unique_cell)){
  data = read.csv(paste(unique_cell[i],"_minfc.csv",sep=""))
  fc = data[order(data[,1]),"min_FC"]
  cell_minfc[,i] = fc
}
rownames(cell_minfc) = genes
colnames(cell_minfc) = unique_cell

cluster_minfc[is.na(cluster_minfc)] <- 0
cell_minfc[is.na(cell_minfc)] <- 0
ctd[[1]][["specificity"]] = cell_minfc
ctd[[2]][["specificity"]] = cluster_minfc

m2h = genes
human.hits = m2h[m2h %in% example_genelist]
human.bg = m2h
# Set the parameters for the analysis
reps=20000 # <- Use 100 bootstrap lists so it runs quickly, for publishable analysis use >10000
# Bootstrap significance testing controlling for transcript length and GC content

full_results = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,
                                         bg=human.bg,reps=reps,annotLevel=2,
                                         geneSizeControl=TRUE,
                                         sctSpecies="human",genelistSpecies="human")
p.adj = p.adjust(full_results$results$p,method = "fdr")
merged_results = cbind(full_results$results,p.adj)

setwd("D:/Yiping/allen_MTG_sce/allen_MTG_sce(2)/review/EWCE")
print(full_results$results[order(full_results$results$p),3:5][1:6,])
write.csv(full_results$results[order(full_results$results$p),],
          paste(reps,"_full_results_cluster.csv",sep=""))
write.csv(merged_results[order(merged_results$p.adj),],
          paste(reps,"_merged_results_cluster.csv",sep=""))

full_results_cell = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,
                                              bg=human.bg,reps=reps,annotLevel=1,
                                              geneSizeControl=TRUE,
                                              sctSpecies="human",genelistSpecies="human")
p.adj_cell = p.adjust(full_results_cell$results$p,method = "fdr")
merged_results_cell = cbind(full_results_cell$results,p.adj_cell)
write.csv(full_results_cell$results[order(full_results_cell$results$p),],
          paste(reps,"_full_results_class.csv",sep=""))
write.csv(merged_results_cell[order(merged_results_cell$p.adj_cell),],
          paste(reps,"_merged_results_class.csv",sep=""))

print(ewce.plot(full_results$results,mtc_method="BH"))
print(ewce.plot(full_results_cell$results,mtc_method="BH"))

generate.bootstrap.plots(sct_data=ctd,hits=human.hits,bg=human.bg,reps=100,annotLevel=1,full_results=full_results,listFileName="VignetteGraphs")