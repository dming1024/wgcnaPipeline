#' @title Pipeline of WGCNA
#'
#' @description This package used to generate data whose distribution is normal
#'
#' @param GeneMatrix,probeAnnotation,traits,method,sftEstimate
#'
#' @return NULL
#'
#' @examples  wgcnaPipeline(GeneMatrix,probeAnnotation,traits,method=mean,sftEstimate=20)
#'
#' @export wgcnaPipeline
#' @import GEOquery,tidyr,clusterProfiler,biomaRt,enrichplot,pathview,WGCNA,stringr

options(device.ask.default = F)
cwd=getwd()
#GSE24759_series_matrix.txt.gz
#probe_annotation.tsv
#mean,median
#如果数据不收敛，直接使用sft=20进行构建基因共表达网络
#临床数据格式这块，需要使用者事前进行整理：多分类变量一定要做出二分类变量；连续变量就是连续变量
#其他应该就没问题了,目前只能分析GEO中的数据

wgcnaPipeline <- function(GeneMatrix,probeAnnotation,traits,method=mean,sftEstimate=20) {
  #读入表达矩阵
      ex_matrix = getGEO(
      filename = GeneMatrix,
      getGPL = F,
      parseCharacteristics = F
    ) %>% exprs(.)


  #输出表达矩阵
  if(!dir.exists("data")){dir.create("data")}
  write.csv(ex_matrix,paste0(cwd,"/data/expression_matrix.csv"))

  #所有样本的基因表达作图
  #保持样本质控图1
  if(!dir.exists("plots")){dir.create("plots")}
  path01=paste(cwd,"plots",sep="/")
  pdf(file = paste(path01,"expressionDistribution.pdf",sep="/"),
      width = ncol(ex_matrix)/10,height = (max(ex_matrix)-min(ex_matrix))/2)
  par(mar=c(7,2,2,2))
  boxplot(ex_matrix[,c(1:ncol(ex_matrix))],
          outline = F,col="lightblue",
          las=2,ylim=c(min(ex_matrix),max(ex_matrix)))
  dev.off()

  #对探针进行注释
  gpl=read.table(probeAnnotation,fill = T,stringsAsFactors = F,sep="\t",quote = "\"",header = T)
  colnames(gpl)=c("ID_REF","GeneSymbol")
  #head(gpl)
  ex_matrix=data.frame(ex_matrix)
  ex_matrix$ID_REF=rownames(ex_matrix)
  ex_matrix_genesymbol=merge(ex_matrix,gpl,by.x = "ID_REF",by.y="ID_REF",all.x=T)


  #对重复的gene求均值，或者中位数
  temp_group=ex_matrix_genesymbol$GeneSymbol
  temp_all=ncol(ex_matrix_genesymbol)-1
  temp_data=ex_matrix_genesymbol[2:temp_all]
  temp_list=apply(temp_data,
                  2,
                  function(x)aggregate(x~temp_group,data=temp_data,method)
  )
  temp_merge=c()
  for(i in 1:length(temp_list)){
    temp_merge=cbind(temp_merge,temp_list[[i]]$x)
  }
  temp_df=data.frame(temp_merge)
  colnames(temp_df)=names(temp_list)
  temp_df$Genesymbol=temp_list[[1]]$temp_group
  filter_genes_matrix=temp_df
  rm(list=ls(pattern = "temp*"))

  #去除空值genes
  filter_genes_matrix=filter_genes_matrix[filter_genes_matrix$Genesymbol !="",]
  rownames(filter_genes_matrix)=filter_genes_matrix$Genesymbol
  filter_genes_matrix=filter_genes_matrix[,-ncol(filter_genes_matrix)]
  write.csv(filter_genes_matrix,
            paste0(cwd,"/data/SymbolGeneMatrix.csv"))

  #读入临床信息
  samples_info=read.csv(traits,stringsAsFactors = F,
                        row.names = 1)

  #WGCNA分析开始
  datExpr0=t(filter_genes_matrix)
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
  sampleTree = hclust(dist(datExpr0), method = "average")
  pdf(file = "plots/01sampleClustering.pdf", width = nrow(datExpr0)/10, height = 9);
  par(mar = c(2,6,4,2),cex=0.6)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  #采用聚类树99%分位数，进行cutoff的设置
  abline(h = quantile(sampleTree$height,0.99), col = "red");
  dev.off()

  clust = cutreeStatic(sampleTree, cutHeight = quantile(sampleTree$height,0.99), minSize = 10)
  #table(clust)
  keepSamples = (clust!=0)
  datExpr = datExpr0[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

#临床信息与表达矩阵进行matched
  Samples = rownames(datExpr)
  traitRows = na.omit(match(Samples, rownames(samples_info)))
  datTraits = samples_info[traitRows,]
  #rownames(datTraits) = datTraits$Sample_geo_accession
  #datTraits=datTraits[,-c(1)]
  collectGarbage()

  #再聚类，同时加上临床信息
  sampleTree2 = hclust(dist(datExpr), method = "average")
  #traitColors = labels2colors(datTraits,
  #                            commonColorCode = F);
  #如果数字居多，则采用numbers2colors()
  traitColors=numbers2colors(datTraits)


  pdf(file="plots/02traits_samples.pdf", nrow(datExpr0)/10, height = 10)
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),cex.dendroLabels = 0.8,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()


  #选择阈值
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  #sizeGrWindow(9, 5)
  #par(mfrow = c(1,2));
  pdf(file="plots/03sft_selected.pdf",width = 4, height = 6)
  cex1 = 0.8;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  #构建网络
  sft$powerEstimate=ifelse(sft$powerEstimate<20,sft$powerEstimate,sftEstimate)
  net = blockwiseModules(datExpr, power = sft$powerEstimate,
                         TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "data/blockwise",
                         verbose = 3,nThreads = 4)
  mergedColors = labels2colors(net$colors)
  pdf(file = "plots/04Net_dendrograms.pdf",width = 15, height = 10)
  for(i in 1:length(table(net$blocks))){
    plotDendroAndColors(net$dendrograms[[i]], mergedColors[net$blockGenes[[i]]],
                        "Module colors",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
  }
  dev.off()


  MEs = net$MEs
  geneTree = net$dendrograms[[1]]
  moduleColors = labels2colors(net$colors)
  save(datExpr,datExpr0,datTraits,sft,net,mergedColors,moduleColors,sft,
       file = "data/step01_networkConstruction-auto.RData")


  #colnames(datTraits)=c('title','ch1','batch','type_')
  #datTraits_factor=as.data.frame(model.matrix(~0+type_,datTraits))
  #当traits为binary时，需要设置robustY=F
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  cor_results=bicorAndPvalue(MEs, datTraits, use = "p",robustY=F)
  moduleTraitCor=cor_results$bicor
  moduleTraitPvalue=cor_results$p


  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  pdf(file = "plots/05Multi_trait_correlation.pdf",width =10, height = 10)
  par(mar = c(10, 10, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  save(datExpr,datTraits,traitColors,sft,net,mergedColors,moduleColors,MEs,moduleTraitCor,moduleTraitPvalue,sft,
       file = "data/step02_networkConstruction-auto.RData")

  #导出所有modules下的genes list
  #grey module不要进行富集分析了
  modNames = substring(names(MEs), 3)
  all_modules=names(table(moduleColors))
  all_modules=all_modules[all_modules!="grey"]
  for (module in all_modules) {
    #确定interesting modules所在的列
    column = match(module, modNames);
    #选择这个modules下的基因
    moduleGenes = moduleColors==module;
    #输出某一module下所有的genes
    interesting_modules_genes=colnames(datExpr)[moduleGenes]
    if(!dir.exists(paste("data",module,sep="/"))){dir.create(paste("data",module,sep="/"))}
    write.csv(interesting_modules_genes,paste0("data/",paste(module,paste0(toupper(module),"_Interesting_Genes.csv"),sep="/")))
    #进行富集分析
    a=bitr(interesting_modules_genes,fromType = "SYMBOL",toType = c("ENSEMBL","ENTREZID"),
           OrgDb = "org.Hs.eg.db")
    #570genes成功注释获得479个
    #clusterprofile似乎经历过大的更新
    go_mf=enrichGO(a$ENTREZID,
                 ont="MF",
                 pvalueCutoff = 0.05,
                 OrgDb = org.Hs.eg.db)
    go_cc=enrichGO(a$ENTREZID,
                 ont="CC",
                 pvalueCutoff = 0.05,
                 OrgDb = org.Hs.eg.db)
    go_bp=enrichGO(a$ENTREZID,
                 ont="BP",
                 pvalueCutoff = 0.05,
                 OrgDb = org.Hs.eg.db)
    go_kegg=enrichKEGG(a$ENTREZID,
                     pvalueCutoff = 0.05,
                     organism="hsa")
    #如果不用adjustPvalue，需要将pAdjustMethod="none"添加到参数中
    setwd(paste(paste(cwd,"data",sep="/"),module,sep="/"))
    #气泡图

    pdf(file = paste0(paste(paste(cwd,"data",sep="/"),module,sep="/"),"/Enrichment_bubble.pdf"),
        width = 7,height = 5)
    #plot.new()
    #CairoPDF(file = paste0(paste(paste(cwd,"data",sep="/"),module,sep="/"),"/enrichment_bubble.pdf"),width = 5,height = 5)
    temp01=dotplot(go_bp)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Biological Process")+theme(plot.title = element_text(hjust = 0.5))
    print(temp01)
    temp02=dotplot(go_cc)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Cell Component")+theme(plot.title = element_text(hjust = 0.5))
    print(temp02)
    temp03=dotplot(go_mf)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Molecular Function")+theme(plot.title = element_text(hjust = 0.5))
    print(temp03)
    temp04=dotplot(go_kegg)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="KEGG Pathway")+theme(plot.title = element_text(hjust = 0.5))
    print(temp04)
    dev.off()

    #输出富集分析结果
    write.csv(go_bp@result,"ENRICHMENT_BIOLOGICALPROGRESS.csv")
    write.csv(go_cc@result,"ENRICHMENT_CELLCOMPONENT.csv")
    write.csv(go_mf@result,"ENRICHMENT_MOLECULARFUNCTION.csv")
    write.csv(go_kegg@result,"ENRICHMENT_KEGG.csv")

    #柱状图
    #pdf(file = paste0(paste(paste(cwd,"data",sep="/"),module,sep="/"),"/Enrichment_barplot.pdf"),
    #    width = 7,height = 5)
    #temp01=barplot(go_bp)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Biological Process")+theme(plot.title = element_text(hjust = 0.5))
    #print(temp01)
    #temp02=barplot(go_cc)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Cell Component")+theme(plot.title = element_text(hjust = 0.5))
    #print(temp02)
    #temp03=barplot(go_mf)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="Molecular Function")+theme(plot.title = element_text(hjust = 0.5))
    #print(temp03)
    #temp04=barplot(go_kegg)+scale_y_discrete(labels=function(y)str_wrap(y,width = 40))+labs(title="KEGG Pathway")+theme(plot.title = element_text(hjust = 0.5))
    #print(temp04)
    #dev.off()

    result_kegg=go_kegg@result
    result_kegg_id=result_kegg[result_kegg$p.adjust<0.05,"ID"]
    #KEGG pathway 可视化
    if(!dir.exists("kegg")){dir.create("kegg")}
    setwd(paste(paste(paste(cwd,"data",sep="/"),module,sep="/"),"kegg",sep="/"))
    source(paste(cwd,"R/pathview_modified.R",sep="/"))
    for (keggid in result_kegg_id) {
      pathview(
        gene.data  = a$ENTREZID,
        pathway.id = keggid,
        species    = "hsa",
        kegg.dir="."
      )
    }
    setwd(cwd)
  }


}

