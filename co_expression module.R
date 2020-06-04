####co_expression module
setwd("E:/Yiping/allen_MTG_sce/R code")
sce <- readRDS("allen_MTG_sce_protein_coding.rds")  #17120 12506
exp_combat1=sce@assays$data$combat_data
symbol <- sce@rowRanges@elementMetadata$gene
rownames(exp_combat1)<-symbol
ENTREZID<-rownames(sce@assays$data$combat_data)
trans <- data.frame(symbol,ENTREZID)
#7011 high variable genes
hvg <- readRDS("hvg_protein_coding.rds")     
hvg_genes <-rownames(hvg)
ix = match(hvg_genes,ENTREZID)

exp_combat <- exp_combat1[ix,]
exp0=as.data.frame(exp_combat) 
cell_type=colData(sce)$Cell_type
unique_cell=unique(cell_type)

y<-read.csv("SFARI-Genes.csv")# asd candidate genes 1054
y=y$gene.symbol

power_list <- vector()
modulenum <- vector()
type = "unsigned"
corFnc = "cor"

for (i in 1:length(unique_cell)){  
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  #Choose a set of soft-thrsholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(dataExpr, powerVector = powers, corFnc = corFnc,
                          networkType = type, verbose = 5)
  sizeGrWindow(9, 5)
  cex1 = 0.9
  figurename <- paste(unique_cell[i],".pdf")
  pdf(file = figurename)
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = unique_cell[i])
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,col="red")
  abline(h=0.85,col="red")
  dev.off()
  power = sft$powerEstimate
  
  if (is.na(power)){
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
    power_list[i]<-paste(unique_cell[i],"NA",power)
  }else { 
    power_list[i]<-paste(unique_cell[i],power,sep = "")
  }
  power 
  
  net = blockwiseModules(dataExpr,  maxBlockSize = nGenes,
                         TOMType="signed", corType= "pearson", #corFnc ="cor"
                         power = power,    #type = "unsigned"
                         networkType = type, minModuleSize = 30, deepSplit = 2,
                         mergeCutHeight = 0.2, numericLabels = TRUE,
                         minKMEtoStay = 0.2,verbose = 5,
                         saveTOMs  = TRUE )###for further analyze
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  table(net$colors)
  table(moduleColors)
  
  write.csv(moduleLabels,file = paste(unique_cell[i],"_moduleLabels.csv",sep = ""))
  write.csv(moduleColors,file = paste(unique_cell[i],"_moduleColors.csv",sep = ""))
  
  mlength = length(unique(net$colors))-1
  modulenum[i] <- paste(unique_cell[i],"_",mlength)
  Module_array <- as.character(rep(1:mlength))
  
  x<-rownames(exp0)
  m_data=intersect(x,y)
  m=length(m_data)
  n=length(x)
  
  p<-c()
  nq_data<-c()
  for(t in 1:mlength) {
    #save moduel infor
    module_name<-Module_array[t]
    m_index<-which(moduleLabels==module_name)
    module_genes<-row.names(exp)
    m_gene<- module_genes[m_index]
    m_gene=data.frame(m_gene)
    filename = paste("M", t, "_gene.csv")
    write.csv(m_gene,file = filename)
    ##hypergeometric analyze
    k_data = m_gene
    k=length(k_data[,1])
    q_data=intersect(y,k_data[,1])
    q=length(q_data)
    nq_data[t] =q
    p_value=1-phyper(q,m,n-m,k)
    p_value
    p[t]<-p_value
  }
  p_adj=p.adjust(p,method = "fdr")
  result_data = rbind(p,p_adj,nq_data)
  rowname<-c("p_value","p_adj","Intersect(ASD,module)"  )
  dimnames(result_data)=list(rowname,Module_array)
  result_data = t(result_data)
  filename <- paste(unique_cell[i],"_ASD_associated.csv")
  write.csv(result_data,file = filename)
}