
library("R.matlab")
library("scater")
library("scran")
library("SingleCellExperiment")
library("preprocessCore")
library("edgeR")
options(stringsAsFactors = FALSE)

#######
#specificity
#####  normalization
cpm_tmm <- function(counts){
  d = DGEList(counts = counts)
  d = calcNormFactors(d, method = "TMM")
  return(cpm(d, prior.count = 10, normalized.lib.sizes = TRUE))
}
setwd("D:/Yiping/allen_MTG_sce/R code")
sce <- readRDS("allen_MTG_sce_protein_coding.rds")  #17120 12506
hvg <- readRDS("hvg_protein_coding.rds")     #7011
hvg_genes <-rownames(hvg)
ix = match(hvg_genes,rownames(sce@assays$data$combat_data))
sce0 <- sce[ix,]
exp_combat1=sce@assays$data$combat_data
symbol <- sce@rowRanges@elementMetadata$gene
rownames(exp_combat1)<-symbol
ENTREZID<-rownames(sce@assays$data$combat_data)
trans <- data.frame(symbol,ENTREZID)
genes=symbol[match(hvg_genes,ENTREZID)]

count_minfc <- function(cell_type,list_other_cells)
{ #compute fold-changes for CTOI vs others
  fc_df = data.frame(rep(1:nrow(filter_count)))
  for(i in 1:length(list_other_cells)){
    fc = rowMeans(filter_count[ , which(sample_types %in% cell_type)]) /
      rowMeans(filter_count[ , which(sample_types %in% list_other_cells[[i]])])
    if(!exists("fc_df")){
      fc_df = data.frame(fc)
    } else {
      fc_df[,i] = fc
    }
  }
  fc_min = apply(fc_df, 1, min)
  fc_min_names = data.frame(genes, fc_min)
  names(fc_min_names) = c("symbol", "min_FC")
  toptable = fc_min_names[order(fc_min_names$min_FC,decreasing = TRUE ), ]
  return(toptable)
}
count_matrix <- assay(sce0) 
filter_count = count_matrix
cell_type = colData(sce)$Cell_type
unique_cell = unique(cell_type)
filter_count <- cpm_tmm(filter_count)
sample_types= colData(sce)$Cell_type
for (i in 1:length(unique_cell)){
  cell = unique_cell[i]
  cell_other = unique_cell[which(unique_cell!=cell)]
  specificity = count_minfc(cell_type = cell,list_other_cells = cell_other)
  write.csv(specificity,file = paste0("cell class/",cell,"_minfc.csv"))
}

count_matrix <- assay(sce0) 
filter_count = count_matrix
cell_cluster = sce0@colData@listData[["cluster"]]
unique_cluster = unique(cell_cluster)
filter_count <- cpm_tmm(filter_count)
sample_types= sce0@colData@listData[["cluster"]]
for (i in 1:length(unique_cluster)){
  cell = unique_cluster[i]
  cell_other = unique_cluster[which(unique_cluster!=cell)]
  specificity = count_minfc(cell_type = cell,list_other_cells = cell_other)
  write.csv(specificity,file = paste0(cell,".csv"))
}
#######

### hypergeometric test
setwd("D:/Yiping/allen_MTG_sce/R code")
sce <- readRDS("allen_MTG_sce_protein_coding.rds")  #17120 12506
exp_combat1=sce@assays$data$combat_data
symbol <- sce@rowRanges@elementMetadata$gene
rownames(exp_combat1)<-symbol
ENTREZID<-rownames(sce@assays$data$combat_data)
trans <- data.frame(symbol,ENTREZID)

hvg <- readRDS("hvg_protein_coding.rds")     #7011
hvg_genes <-rownames(hvg)
ix = match(hvg_genes,trans[,2])
exp_combat <- exp_combat1[ix,]
exp0=as.data.frame(exp_combat) 
cell_type=colData(sce)$Cell_type
unique_cell=unique(cell_type)

unique_cell<-c("Glutamatergic neuron", "GABAergic interneuron","astrocyte","Oligodendrocyte precursor cell", "Oligodendrocyte" ,"Microglial")
markfile = c("GABA","Gluta","Ast","OPC","Oli","Mic")
mlength=c("13","8","28","71","70","20")
x<-rownames(exp0)
n=length(x)

l = length(markfile)
enriched = matrix(nrow=l,ncol=l)
rownames(enriched) = markfile
colnames(enriched) = unique_cell
library(stringr)
for (i_mark in 1:length(markfile)) {
  filename<-paste("d15_",markfile[i_mark],"_hvg.csv",sep = "")
  y<-read.csv(file = filename)# cell type  markgenes 
  y=y$gene_symbol[1:500]  ##top500 specificity
  y = symbol[match(y,ENTREZID)]
  
  for (i in 1:length(unique_cell)){
    m_data=intersect(x,y)
    m=length(m_data)
    p<-c(1);nq_data<-c(1)
    nmodule = mlength[i]
    Module_array <- as.character(rep(1:nmodule))
    for (t in 1:nmodule) {
      mgene <- read.csv(paste("M", t, "_gene.csv")) #
      k_data = mgene[,2]
      k=length(k_data)
      q_data=intersect(y,k_data)
      q=length(q_data)
      nq_data[t] =q
      p_value=1-phyper(q,m,n-m,k)
      p_value
      p[t]<-p_value
    }
    p_adj=p.adjust(p,method = "fdr")
    result_data = rbind(p,p_adj,nq_data)
    rowname<-c("p_value","p_corrected","Intersect(spec_hvg,module)")
    dimnames(result_data)=list(rowname,Module_array)
    result_data = t(result_data)
    result = data.frame(ix = rep(1:nrow(result_data)), p_adj = result_data[,2])
    ix = result[which(result$p_adj<0.1),1]
    
    if (length(ix)!=0) {enriched[i_mark,i] = str_c(ix,collapse=',') }
    if (length(ix)==0) {enriched[i_mark,i] = "na" }
    filename <- paste0(unique_cell[i],"-",markfile[i_mark],"_spechvg.csv")
    write.csv(result_data,file = filename)
  }
}
write.csv(enriched,"enriched module.csv")
getwd()
