#!/usr/bin/env Rscript 
source("/data/wqihao/Project.sc/scrna_analysis_prepare.R")

if(!interactive()){
  option_list <- list(make_option(c("-f", "--fastqDir"),type="character", help=".*fastq.gz file location"), 
                      make_option(c("-i", "--inputDir"),type="character", help="Cellranger Results Directory"),
                      make_option(c("-o", "--outDir"),type="character", help="output Directory"),
                      make_option(c("-s", "--species"),type="character", help="species of single cell [hsa, mmu]"),
                      make_option(c("-a", "--annoDB"),type="character", help="DB use to anno celltype"),
                      make_option(c("-d", "--DBdir"),type="character", default = "//data/wqihao/database/refForScrna/singleR.DB/",help="DB directory use to anno celltype"))
  opts <- parse_args(OptionParser(option_list=option_list))
}else{
  opts= list()
  opts$fastqDir = "/data/wqihao/Project.sc/cellrangeRun/inputFile/GP1_fastqs/"
  opts$inputDir = "/data/wqihao/Project.sc/cellrangeRun/runDir/scRNA_20210511_GP1/"
  opts$outDir   = "/data/wqihao/Project.sc/testpipline/"
  opts$species  = "mmu"
  opts$annoDB   = "ImmGenData"
  opts$DBdir    = "/data/wqihao/database/refForScrna/singleR.DB/"
}

BIN = "/data/wqihao/Project.sc/workflow/bin/"
CATLOG("0.Parse args")

pipeline.log = file.path(opts$outDir,"pipeline.log")
DIRCREATE(pipeline.log)

###### Check options 
if(length(opts) < 4){
  cat ("Use: Rscript  %prog -h see more Useage \n")
  cat ("Version: v1.0\n")
  cat ("Date: 2021-06-03\n")
  break
}

INPUTFASTQ.DIR = opts$fastqDir
INPUT.DIR = opts$inputDir
INPUT.OUTS.DIR = file.path(INPUT.DIR,"outs")
sampleName = basename(INPUT.DIR)
species = opts$species

# 整理后结果文件夹 -------
OUT.DIR = opts$outDir
DIRCREATE(OUT.DIR)
OUT.RAWDATA.DIR = file.path(OUT.DIR,"raw",c("bamfile","fastqfile"))
OUT.RESTULE.DIR = file.path(OUT.DIR,"out")
OUT.RAWDATA.SAMPLEN.DIR = file.path(OUT.RAWDATA.DIR,sampleName)
OUT.RESTULE.SAMPLEN.DIR = file.path(OUT.RESTULE.DIR, sampleName) 

# Setp1 整理输出原始bam,fastq-------
setwd(OUT.DIR)
DIRCREATE(file.path(OUT.RAWDATA.SAMPLEN.DIR[1],c("count","vdj_b")))

bamFiles = list.files(INPUT.DIR,pattern = "bam",recursive = T) %>% grep(pattern = "^outs",value = T)
bamFiles.raw = bamFiles %>% grep(pattern = "multi",value = T)
bamFiles.filter = bamFiles %>% grep(pattern = "per_sample_outs",value = T)
OUT.bamFiles.raw = file.path(FDIR(bamFiles.raw),paste0("raw_",basename(bamFiles.raw))) %>% file.path(OUT.RAWDATA.SAMPLEN.DIR[1],.)
OUT.bamFiles.filter = file.path(FDIR(bamFiles.filter),paste0("filtered_",basename(bamFiles.filter))) %>% file.path(OUT.RAWDATA.SAMPLEN.DIR[1],.)
COPYFILE(from.name = file.path(INPUT.DIR,bamFiles),to.name = c(OUT.bamFiles.raw,OUT.bamFiles.filter),step = "Copy bam files")


fastqFiles = list.files(INPUT.DIR,pattern = "fastq",recursive = T) %>% grep(pattern = "^outs",value = T)
fastqFiles.raw = fastqFiles %>% grep(pattern = "multi",value = T)
fastqFiles.filter = fastqFiles %>% grep(pattern = "per_sample_outs",value = T)

DIRCREATE(file.path(OUT.RAWDATA.SAMPLEN.DIR[2],c("vdj_b")))
OUT.fastqFiles.raw = file.path(FDIR(fastqFiles.raw),paste0("raw_",basename(fastqFiles.raw))) %>% file.path(OUT.RAWDATA.SAMPLEN.DIR[2],.)
OUT.fastqFiles.filter = file.path(FDIR(fastqFiles.filter),paste0("filtered_",basename(fastqFiles.filter))) %>% file.path(OUT.RAWDATA.SAMPLEN.DIR[2],.)
COPYFILE(from.name = file.path(INPUT.DIR,fastqFiles),to.name = c(OUT.fastqFiles.raw,OUT.fastqFiles.filter),step = "Copy fastq files")

# Setp2 整理输出文件-------
OUT.TREE.DIR = file.path(OUT.RESTULE.SAMPLEN.DIR,c("1.Data_assess","2.Basic_analysis","3.Cellranger_advanced","4.Seurat","5.Enrichment","6.GSEA","src"))
# 1.Data_assess --------
CATLOG("1.Data_assess")
OUT.Data_assess.DIR = OUT.TREE.DIR[1]
DIRCREATE(OUT.Data_assess.DIR)
if(!file.exists(GetoptLong::qq("@{pipeline.log}/1.Data_assess.flage"))){
  rawFastqFiles = list.files(INPUTFASTQ.DIR,pattern = "_R.*.fastq.gz",full.names = T,recursive = T)
  rawFastqFiles.5GEX = rawFastqFiles %>% grep(pattern = "5gex",value = T) 
  rawFastqFiles.BCR = rawFastqFiles %>% grep(pattern = "BCR",value = T) 
  rawFastqFiles.5GEX.tmp = paste(rawFastqFiles.5GEX,collapse = " ")
  rawFastqFiles.5GEX.QC = file.path(OUT.Data_assess.DIR,gsub(pattern = ".fastq.gz",replacement = "_fastqc.zip",basename(rawFastqFiles.5GEX)))
  cmd = GetoptLong::qq("fastqc -o @{OUT.Data_assess.DIR} -t 8 -f fastq @{rawFastqFiles.5GEX.tmp}")
  Cmd(cmd = cmd,cl = 10,step = "fastq QC",out = rawFastqFiles.5GEX.QC[1])
  flage(pipeline.log,"1.Data_assess")
}

# 2.Basic_analysis -------- 
CATLOG("2.Basic_analysis")
OUT.Basic_analysis.DIR = file.path(OUT.TREE.DIR[2],c("2.1.raw_feature_bc_matrix","2.2.filtered_feature_bc_matrix","2.3.h5_files"))

rawFeatureMatrix = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% 
  grep(pattern = "raw_feature_bc_matrix",value = T) 
rawFeatureMatrix.files = list.files(rawFeatureMatrix,full.names = T)
rawH5 = list.files(INPUT.OUTS.DIR,pattern = "raw.*\\.h5$",recursive = T,full.names = T)
filteredFeatureMatrix = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% 
  grep(pattern = "sample_feature_bc_matrix",value = T)
filteredFeatureMatrix.files = list.files(filteredFeatureMatrix,full.names = T)
filteredH5 = list.files(INPUT.OUTS.DIR,pattern = "sample.*\\.h5$",recursive = T,full.names = T)
DIRCREATE(OUT.Basic_analysis.DIR[1])
DIRCREATE(OUT.Basic_analysis.DIR[2])
DIRCREATE(OUT.Basic_analysis.DIR[3])

metricSummaryJson = GETFILE(input.dir = INPUT.DIR,pattern = "metrics_summary_json") %>% grep(pattern = "SUMMARIZE_REPORTS",value = T) %>% '['(3) # -----------------------
metricSummaryCsv = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "metrics_summary.csv")
web_summary.html = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "web_summary.html")

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Basic_analysis.flage"))){
  COPYFILE(from.name = rawFeatureMatrix.files[1],to.name = OUT.Basic_analysis.DIR[1],step = "Copy raw barcodes.tsv.gz")
  COPYFILE(from.name = rawFeatureMatrix.files[2],to.name = OUT.Basic_analysis.DIR[1],step = "Copy raw features.tsv.gz")
  COPYFILE(from.name = rawFeatureMatrix.files[3],to.name = OUT.Basic_analysis.DIR[1],step = "Copy raw matrix.mtx.gz")
  COPYFILE(from.name = filteredFeatureMatrix.files[1],to.name = OUT.Basic_analysis.DIR[2],step = "Copy filtered barcodes.tsv.gz")
  COPYFILE(from.name = filteredFeatureMatrix.files[2],to.name = OUT.Basic_analysis.DIR[2],step = "Copy filtered features.tsv.gz")
  COPYFILE(from.name = filteredFeatureMatrix.files[3],to.name = OUT.Basic_analysis.DIR[2],step = "Copy filtered matrix.mtx.gz")
  COPYFILE(from.name = rawH5,to.name = OUT.Basic_analysis.DIR[3],step = "Copy rawH5")
  COPYFILE(from.name = filteredH5,
           to.name = file.path(OUT.Basic_analysis.DIR[3],gsub(pattern = "sample_",replacement = "filtered_",basename(filteredH5))),
           step = "Copy filteredH5")
  COPYFILE(from.name = metricSummaryCsv,to.name = file.path(OUT.TREE.DIR[2], basename(metricSummaryCsv)),step = "Copy metricSummaryCsv")
  COPYFILE(from.name = metricSummaryJson,to.name = file.path(OUT.TREE.DIR[2], basename(metricSummaryJson)),step = "Copy metricSummaryJson")
  COPYFILE(from.name = web_summary.html,to.name = file.path(OUT.TREE.DIR[2], basename(web_summary.html)),step = "Copy web_summary.html")
  
  seurat_object = CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),project = "cellranger")
  seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
  colnames(seurat_object@meta.data)[2:3] = c("nUMI","nGene")
  PLOTHISTOGRAM <- function(x,xlab,ylab,title){
    data <- data.frame(x=x)
    r <- hist(x,plot = F,breaks = 30)
    x.aes = max(x) * 0.9
    y.aes = max(r$counts) * 0.9
    summary.x = round(summary(x),3)
    ggplot2::ggplot( data=data, ggplot2::aes(x=x)) + 
      ggplot2::geom_histogram(fill="white",col = "black",alpha=1,bins=30) +
      ggplot2::ggtitle(title) +
      ggplot2::theme_classic() + 
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) +
      ggplot2::annotate("text",x = x.aes,y = y.aes,label = paste0(names(summary.x)[1],": ",summary.x[1])) +
      ggplot2::annotate("text",x = x.aes,y = y.aes-10,label = paste0(names(summary.x)[2],": ",summary.x[2])) +
      ggplot2::annotate("text",x = x.aes,y = y.aes-20,label = paste0(names(summary.x)[3],": ",summary.x[3])) +
      ggplot2::annotate("text",x = x.aes,y = y.aes-30,label = paste0(names(summary.x)[4],": ",summary.x[4])) +
      ggplot2::annotate("text",x = x.aes,y = y.aes-40,label = paste0(names(summary.x)[5],": ",summary.x[5])) +
      ggplot2::annotate("text",x = x.aes,y = y.aes-50,label = paste0(names(summary.x)[6],": ",summary.x[6]))
  }
  p1 = PLOTHISTOGRAM(x = seurat_object$nGene,xlab = "Genes per cell",ylab = "cell counts",title = "Gene Distribution")
  p2 = PLOTHISTOGRAM(x = seurat_object$nUMI,xlab = "UMIs per cell",ylab = "cell counts",title = "UMI Distribution")
  p3 = VlnPlot(seurat_object,features = c("nGene", "nUMI"),group.by = "orig.ident", ncol = 2,pt.size = 0.4) 
  p4 = FeatureScatter(seurat_object, feature1 = "nUMI", feature2 = "nGene",cols = "black") + ggplot2::theme(legend.title=ggplot2::element_blank(),legend.text = ggplot2::element_blank())
  nGene.histogram = file.path(OUT.TREE.DIR[2],"nGene.histogram.pdf")
  nUMI.histogram = file.path(OUT.TREE.DIR[2],"nUMI.histogram.pdf")
  nGene_and_nUMI_correlation = file.path(OUT.TREE.DIR[2],"nGene_and_nUMI_correlation.pdf")
  nGene_and_nUMI_distribution = file.path(OUT.TREE.DIR[2],"nGene_and_nUMI_distribution.pdf")
  ggplot2::ggsave(filename = nGene.histogram,plot =p1 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nUMI.histogram,plot =p2 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nGene_and_nUMI_correlation,plot =p4 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nGene_and_nUMI_distribution,plot =p3 ,width = 6,height = 6,dpi = "retina")
  
  flage(pipeline.log,"2.Basic_analysis")
}


# 3.Cellranger_advanced --------
CATLOG("3.Cellranger_advanced")
OUT.Cellranger_advanced.DIR = file.path(OUT.TREE.DIR[3],c("cloupe","Cluster_specific_genes","correlation","pca_cluster","tsne_cluster","umap_cluser"))

# cloupe
cloupeFiles = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "cloupe") %>% grep(pattern = "per_sample_outs",value = T)
vloupeFiles = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "vloupe")
DIRCREATE(file.path(OUT.Cellranger_advanced.DIR[1],c("count","vdj_b")))
OUT.cloupeFiles = file.path(FDIR(cloupeFiles),paste0(basename(cloupeFiles))) %>% file.path(OUT.Cellranger_advanced.DIR[1],.)
OUT.vloupeFiles = file.path(FDIR(vloupeFiles),paste0(basename(vloupeFiles))) %>% file.path(OUT.Cellranger_advanced.DIR[1],.)
COPYFILE(from.name = cloupeFiles,to.name = OUT.cloupeFiles,step = "Copy cloupeFiles")
COPYFILE(from.name = vloupeFiles,to.name = OUT.vloupeFiles,step = "Copy vloupeFiles")
# Cluster_specific_genes
analysisDIR = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% grep(pattern = "analysis",value = T) 
analysis.TSNE = analysisDIR %>% grep(pattern = "tsne/",value = T) %>% list.files(full.names = T)
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
analysis.diffexp = analysisDIR %>% grep(pattern = "diffexp/",value = T) %>% list.files(full.names = T)
seurat_object = CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),project = "cellranger")
seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
log2.gene.exp = log2(seurat_object@assays$RNA@counts+1) 
if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Cellranger_advanced.Cluster_specific_genes.flage"))){
  CATLOG("  Run Cluster_specific_genes")
  for (n in 1:10) {
    cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    cat(clusterType)
    tmp.DIR = file.path(OUT.Cellranger_advanced.DIR[2],clusterType)
    DIRCREATE(tmp.DIR)
    tmp = read.table(analysis.TSNE,sep = ",",header = T,stringsAsFactors = F)
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    tmp.diffexp = read.table(analysis.diffexp[n],sep = ",",header = T,stringsAsFactors = F)
    clusters = colnames(tmp.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()
    COPYFILE(from.name = analysis.diffexp[n],to.name = file.path(tmp.DIR,paste0(clusterType,".",basename(analysis.diffexp[n]))),step = "")
    tmp.matrix.pvalue =  tmp.diffexp  %>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) %>% 
      dplyr::select(matches("(ID)|(Adjusted.p.value)"))
    sig.gene.ids = apply(tmp.matrix.pvalue[,-1] < 0.05, 1, sum)!=0
    tmp.matrix.log2FC =  tmp.diffexp[sig.gene.ids,]  %>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) %>% 
      dplyr::select(matches("(ID)|(Log2.fold)")) %>% data.frame(stringsAsFactors = F,check.names = F)
    rownames(tmp.matrix.log2FC) = tmp.matrix.log2FC$ID
    tmp.matrix.log2FC$ID = NULL
    tmp.matrix.log2FC = t(tmp.matrix.log2FC)
    tmp.gene.group.log2FC = subset(tmp.diffexp[sig.gene.ids,]%>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) ,
                                   ID%in%colnames(tmp.matrix.log2FC))%>% 
      dplyr::select(matches("(ID)|(Log2.fold)")) %>% data.frame(stringsAsFactors = F,check.names = F)
    tmp.gene.group.pvalue = subset(tmp.diffexp[sig.gene.ids,]%>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) ,
                                   ID%in%colnames(tmp.matrix.log2FC))%>% 
      dplyr::select(matches("(ID)|(Adjusted.p.value)")) %>% data.frame(stringsAsFactors = F,check.names = F)
    tmp.gene.group.log2FC = tmp.matrix.log2FC
    set.seed(1233)
    ht = ComplexHeatmap::Heatmap(matrix = tmp.gene.group.log2FC,show_column_names = F)
    p = as.ggplot(~draw(ht))
    outname = paste0(clusterType,".pheatmap.pdf")
    ggplot2::ggsave(filename = file.path(tmp.DIR,outname),plot =  p ,width = 20,height =6 )
    
    plot.list2 = list()
    for (cluster in clusters) {
      cat(cluster)
      # filter top 20 genes
      p.cols = grep(pattern = paste0(cluster,".Adjusted.p.value"),colnames(tmp.diffexp))
      exp.cols = grep(pattern = paste0(cluster,".Mean.Counts"),colnames(tmp.diffexp))
      fc.cols = grep(pattern = paste0(cluster,".Log2.fold.change"),colnames(tmp.diffexp))
      top20genes = tmp.diffexp[order(tmp.diffexp[,fc.cols],decreasing = T)&tmp.diffexp[,p.cols]<0.05,] %>% '['(,2) %>% head(n=20)
      if(length(top20genes)==0){
        write.table(x = outname,file = file.path(OUT.RESTULE.DIR,"err.txt"),append = T,quote = F,col.names = F)
        next
      }
      plot.list = list()
      length(tmp.clusters$Cluster)
      plot.list[[1]] = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp.clusters$Cluster,type = "Cluster",title ="t-SNE")
      
      plot.list2[[cluster]] = volcanoPlot(logFC = tmp.diffexp[,fc.cols],pValue =tmp.diffexp[,p.cols] ,FC = 3,main = paste0(clusterType,":",cluster),size = 0.3)
      outname = paste0(clusterType,".",gsub(pattern = "\\.",replacement = "_",cluster),".pdf")
      
      for (gene in top20genes) {
        plot.list[[gene]] = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = as.numeric(log2.gene.exp[gene,]),
                                      title = paste(c(as.character(tmp.diffexp[tmp.diffexp$Feature.Name==gene,1:2]),"Log2 Total UMI Count Per Cell"),collapse = ","),type = "exp")
      }
      ncols = 3
      pdf(file = file.path(tmp.DIR,outname), width = 25,height =40,onefile = F)
      multiplot(plotlist = plot.list,cols = ncols)
      dev.off()
    }
    outname = paste0(clusterType,".diffexp.volvano.pdf")
    if(length(plot.list2)==0){
      write.table(x = outname,file = file.path(OUT.RESTULE.DIR,"err.txt"),append = T,quote = F,col.names = F)
      next
    }
    ncols = 2
    pdf(file = file.path(tmp.DIR,outname),width = 7,height = 16,onefile = F)
    multiplot(plotlist = plot.list2,cols = ncols)
    dev.off()
  }
  flage(pipeline.log,"3.Cellranger_advanced.Cluster_specific_genes")
  
}else{
  CATLOG("  Skip Cluster_specific_genes")
}

# correlation
DIRCREATE(OUT.Cellranger_advanced.DIR[3])

if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Cellranger_advanced.correlation.flage"))){
  CATLOG("  Run correlation")
  for (n in 1:10) {
    cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    # cat(clusterType)
    tmp.diffexp = read.table(analysis.diffexp[n],sep = ",",header = T,stringsAsFactors = F)
    
    tmp.Mean.Counts = tmp.diffexp %>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) %>% 
      dplyr::select(matches("(ID)|(Mean.Counts)"))
    
    corr <- cor(apply(t(tmp.Mean.Counts[,-1]),1,scale), method = 'pearson')
    p <-  pheatmap::pheatmap(corr,cluster_cols = F,cluster_rows = F,display_numbers = TRUE,
                             color = colorRampPalette(c("yellow","red"))(100))
    outname = paste0(clusterType,".","subpopulation.correlation",".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[3],outname),plot =  p ,width = 6,height =6 )
  }
  flage(pipeline.log,"3.Cellranger_advanced.correlation")
}else{
  CATLOG("  Skip correlation")
}

#  pca
DIRCREATE(OUT.Cellranger_advanced.DIR[4])
analysis.PCA = analysisDIR %>% grep(pattern = "pca/",value = T) 
analysis.PCA.clusters =  GETFILE(input.dir = analysis.PCA,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.PCA.clusters,sep = ",",header = T,stringsAsFactors = F)
if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Cellranger_advanced.pca.flage"))){
  CATLOG("  Run pca")
  for (n in 1:10) {
    cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$PC.1 ,tsneY = tmp$PC.2,Clusters = tmp.clusters$Cluster,type = "Cluster",
                  xlab = "PC1",ylab = "PC2",
                  title =paste0(clusterType,"\n","PCA projection of Cells Colored by Autometed Clustering"))
    outname = paste0("pca.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[4],outname),plot = p,width = 6,height =6 )
  }
  flage(pipeline.log,"3.Cellranger_advanced.pca")
}else{
  CATLOG("  Skip pca")
}

# tsne
DIRCREATE(OUT.Cellranger_advanced.DIR[5])
analysis.TSNE = analysisDIR %>% grep(pattern = "tsne/",value = T) 
analysis.TSNE.clusters =  GETFILE(input.dir = analysis.TSNE,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.TSNE.clusters,sep = ",",header = T,stringsAsFactors = F)
if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Cellranger_advanced.tsne.flage"))){
  CATLOG("  Run tsne")
  for (n in 1:10) {
    cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    # cat(clusterType)
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp.clusters$Cluster,type = "Cluster",xlab = "TSNE.1",ylab = "TSNE.2",title =paste0(clusterType,"\n","t-SNE projection of Cells Colored by Autometed Clustering"))
    outname = paste0("tsne.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[5],outname),plot = p,width = 6,height =6 )
  }
  p = PLOT.TSNE2(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = seurat_object$nFeature_RNA,type = "exp",
                 xlab = "TSNE.1",ylab = "TSNE.2",title =paste0( "t-SNE Projection of Cells Colored by Gene Counts"),
                 labs = "Gene Counts")
  outname = paste0("tsne.nGene",".pdf")
  ggplot2::ggsave(filename = file.path(OUT.TREE.DIR[3],outname),plot = p,width = 6,height =6 )
  p = PLOT.TSNE2(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = seurat_object$nCount_RNA,type = "exp",
                 xlab = "TSNE.1",ylab = "TSNE.2",title =paste0( "t-SNE Projection of Cells Colored by UMI Counts"),
                 labs = "UMI Counts")
  outname = paste0("tsne.nUMI",".pdf")
  ggplot2::ggsave(filename = file.path(OUT.TREE.DIR[3],outname),plot = p,width = 6,height =6 )
  flage(pipeline.log,"3.Cellranger_advanced.tsne")
}else{
  CATLOG("  Skip tsne")
}

# umap
DIRCREATE(OUT.Cellranger_advanced.DIR[6])
analysis.UMAP = analysisDIR %>% grep(pattern = "umap/",value = T) 
analysis.UMAP.clusters =  GETFILE(input.dir = analysis.UMAP,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.UMAP.clusters,sep = ",",header = T,stringsAsFactors = F)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Cellranger_advanced.umap.flage"))){
  CATLOG("  Run umap")
  for (n in 1:10) {
    cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    # cat(clusterType)
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$UMAP.1 ,tsneY = tmp$UMAP.2,Clusters = tmp.clusters$Cluster,type = "Cluster",xlab = "UMAP.1",ylab = "UMAP.2",title =paste0(clusterType,"\n","UMAP projection of Cells Colored by Autometed Clustering"))
    outname = paste0("umap.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[6],outname),plot = p,width = 6,height =6 )
  }
  # 
  analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
  stat.cell = function(n) {
    # cat(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    # cat(clusterType)
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    res = data.frame(Cluster.Method = clusterType,table(tmp.clusters$Cluster),stringsAsFactors = F,check.names = F)
    colnames(res)[2:3] = c("Cluster.ID","Cell_Counts")
    return(res)
  }
  nCell_per_cluster = do.call(rbind,lapply(1:10, stat.cell))
  write.table(nCell_per_cluster,file = file.path(OUT.TREE.DIR[3],"nCell_per_cluster.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  # n=analysis.clusters[1]
  f.read <- function(n){
    clusterType = basename(dirname(n))
    res = data.table::fread(n)
    colnames(res)[2] = clusterType
    return(res)
  }
  tmp.TSNE = read.table(analysis.TSNE.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.UMP = read.table(analysis.UMAP.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.PCA = read.table(analysis.PCA.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.seurat = data.frame(Barcode=paste0(rownames(seurat_object@meta.data),"-1"),nGene=seurat_object@meta.data$nFeature_RNA,nUMI=seurat_object@meta.data$nCount_RNA,mt_percent=seurat_object@meta.data$percent.mt)
  
  res.list = lapply(analysis.clusters, f.read)
  res.list = append(res.list,list(tmp.TSNE))
  res.list = append(res.list,list(tmp.UMP))
  res.list = append(res.list,list(tmp.PCA))
  res.list = append(res.list,list(tmp.seurat))
  summary_cell.pca_tsne_clustering = purrr::reduce(res.list,merge)
  write.table(summary_cell.pca_tsne_clustering,file = file.path(OUT.TREE.DIR[3],"summary_cell.pca_tsne_clustering.csv"),quote = F,sep = ",",col.names = T,row.names = F)
  flage(pipeline.log,"3.Cellranger_advanced.umap")
}else{
  CATLOG("  Skip umap")
}

# 4.Seurat--------
CATLOG("4.Seurat")
OUT.Seurat.DIR = file.path(OUT.TREE.DIR[4],c("1.Expression","2.PCA_analysis","3.Cluster_and_diff","4.celltype_annotation"))
DIRCREATE(OUT.Seurat.DIR)

if(!file.exists(file.path(OUT.DIR,"seurat_object.filted.RDS"))){
  seurat_object <- CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),min.cells = 3,min.features = 200,project = "cellranger")
  seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
  Max.nUMI  = stats::quantile(seurat_object$nCount_RNA,probs=c(0.99))  %>%  as.numeric()
  seurat_object <- subset(seurat_object, percent.mt <= 25&nFeature_RNA>200& nCount_RNA < Max.nUMI)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(object = seurat_object,selection.method = "vst")
  seurat_object <- ScaleData(seurat_object,features = rownames(seurat_object))     
  seurat_object <- RunPCA(object= seurat_object,pc.genes=VariableFeatures(object = seurat_object),npcs = 50,nfeatures.print = 20,verbose=F,ndims.print = 1:50)  
  seurat_object <- JackStraw(seurat_object, num.replicate = 100,dims = 50)
  seurat_object <- ScoreJackStraw(seurat_object, dims = 1:50)
  pcsNumber     <- which(seurat_object@reductions$pca@jackstraw@overall.p.values[,2] < 0.01)
  seurat_object <- FindNeighbors(seurat_object, dims = pcsNumber)
  seurat_object <- FindClusters(seurat_object, resolution = 1)
  seurat_object <- RunUMAP(seurat_object, dims = pcsNumber)
  seurat_object <- RunTSNE(seurat_object, dims = pcsNumber)
  saveRDS(seurat_object,file = file.path(OUT.DIR,"seurat_object.filted.RDS"))
}else{
  seurat_object <- readRDS(file.path(OUT.DIR,"seurat_object.filted.RDS"))  
}


#1.Expression
if(!file.exists(GetoptLong::qq("@{pipeline.log}/4.Seurat.Expression.flage"))){
  CATLOG("  Run Expression")
  p1  = VlnPlot(seurat_object,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident", ncol = 3,pt.size = 0.1)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[1],"1.nGene_nUMI_percent.mito.pdf"),plot = p1,width = 6,height =5)
  p2. = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by =  "orig.ident")
  p3. = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by =  "orig.ident")
  P   = p2.+p3.
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[1],"2.Percent.mito_and_nGene_associate_with_nUMI.pdf"),plot = P,width = 10,height =5)
  top10 = head(VariableFeatures(seurat_object), 10)
  plot1 = VariableFeaturePlot(seurat_object,pt.size = 0.05)
  p3    = LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[1],"3.Detect_variable_genes.pdf"),plot = p3,width = 6,height =5)
  flage(pipeline.log,"4.Seurat.Expression")
}else{
  CATLOG("  Skip Expression")
}

# 2.pca_analysis
if(!file.exists(GetoptLong::qq("@{pipeline.log}/4.Seurat.pca_analysis.flage"))){
  CATLOG("  Run pca_analysis")
  pca50_vargene = seurat_object@reductions$pca@feature.loadings
  pca50_vargene = data.frame(Gene = rownames(pca50_vargene),pca50_vargene)
  write.table(pca50_vargene,file = file.path(OUT.Seurat.DIR[2],"1.pca50_vargene.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  sink(file.path(OUT.Seurat.DIR[2],"2.PC50_top20_sep.txt"))
  print(seurat_object@reductions$pca, dims = 1:50, nfeatures = 20)
  sink()
  p4 = VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca",ncol = 2)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],"3.Top30_gene_with_PC1_PC2.pdf"),plot = p4,width = 6,height =5)
  p5 = ggplotify::as.ggplot(~DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE))
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],"4.Top30_gene_PC1_heatmap.pdf"),plot = p5,width = 6,height =5)
  pdf(file = file.path(OUT.Seurat.DIR[2],"5.PC50_heatmap.pdf"),width = 10,height = 25,onefile = F)
  DimHeatmap(seurat_object, dims = 1:50,cells = 500,ncol = 5)
  dev.off()
  pcsNumber = which(seurat_object@reductions$pca@jackstraw@overall.p.values[,2] < 0.01)
  p7 = JackStrawPlot(seurat_object,reduction = "pca",dims = pcsNumber)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],"6.PC_significant.pdf"),plot = p7,width = 9,height =8)
  p8 = ElbowPlot(seurat_object,ndims = 50)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],"7.PC50_sd.pdf"),plot = p8,width = 8,height =8)
  flage(pipeline.log,"4.Seurat.pca_analysis")
}else{
  CATLOG("  Skip pca_analysis")
}


# 3.Cluster_and_diff
OUT.ClusterAndDiff.DIR = file.path(OUT.Seurat.DIR[3],c("2.Cluster_diff","3.Cluster_marker_gene_plot"))
DIRCREATE(OUT.ClusterAndDiff.DIR)
if(!file.exists(GetoptLong::qq("@{pipeline.log}/4.Seurat.Cluster_and_diff.flage"))){
  CATLOG("  Run Cluster_and_diff")
  p9  = DimPlot(seurat_object, reduction = "tsne", label = T,pt.size = 0.5)
  p10 = DimPlot(seurat_object, reduction = "umap",label = T, pt.size = 0.5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"1.Cell_cluster_byTSNE.pdf"),plot = p9,width = 6,height =5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"1.Cell_cluster_byUMAP.pdf"),plot = p10,width = 6,height =5)
  seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  cluster.averages <- AverageExpression(seurat_object)
  
  cluster.averages <- data.frame(gene=rownames(cluster.averages$RNA),cluster.averages$RNA,check.names = F)
  All_cluster_up_marker = data.frame(gene=rownames(seurat_object.markers),seurat_object.markers[,c(1:4,6)],check.names = F)
  All_cluster_up_marker = merge(All_cluster_up_marker,cluster.averages)
  All_cluster_up_marker = All_cluster_up_marker %>% dplyr::arrange(cluster,desc(cluster))
  write.table(All_cluster_up_marker,file = file.path(OUT.Seurat.DIR[3],"1.All_cluster_up_marker.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  top20 = seurat_object.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  write.table(top20[,-5],file = file.path(OUT.Seurat.DIR[3],"4.Top20_cluster_pos_marker_forheatmap.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  # pdf(file = file.path(OUT.Seurat.DIR[3],"5.Cluster_top_express_heatmap.pdf"),width = 10,height = 10,onefile = F)
  
  # dev.off()
  p = as.ggplot(Seurat::DoHeatmap(seurat_object, features = top20$gene,label=FALSE))
  ggplot2::ggsave(filename =file.path(OUT.Seurat.DIR[3],"5.Cluster_top_express_heatmap.pdf"),width = 10,height = 10 )
  
  Cluster_info = data.frame(barcode = colnames(seurat_object),seurat_object@reductions$umap@cell.embeddings,seurat_object@reductions$tsne@cell.embeddings,cluster=seurat_object@meta.data$seurat_clusters)
  write.table(Cluster_info,file = file.path(OUT.Seurat.DIR[3],"Cluster_info.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  Cluster_num_stat = data.frame(table(seurat_object@meta.data$seurat_clusters),check.names = F)
  colnames(Cluster_num_stat) = c("cluster","number")
  write.table(Cluster_num_stat,file = file.path(OUT.Seurat.DIR[3],"Cluster_num_stat.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  for (i in unique(seurat_object$seurat_clusters)) {
    cat(i)
    cluster.markers <- FindMarkers(seurat_object, ident.1 = i, min.pct = 0.25,logfc.threshold = 0.25)
    cluster.markers = data.frame(gene=rownames(cluster.markers),cluster.markers,check.names = F)
    cluster.markers = merge(cluster.markers,cluster.averages)
    outname = paste0("Cluster",i,".differential_and_annoation.xls")
    write.table(cluster.markers,file = file.path(OUT.ClusterAndDiff.DIR[1],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    top20.tmp = subset(top20,cluster==i)
    cluster.markers.top20 = subset(cluster.markers,gene%in%top20.tmp$gene)
    outname = paste0("Cluster",i,"_top20_gene_diffInf.xls")
    write.table(cluster.markers.top20,file = file.path(OUT.ClusterAndDiff.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    Cluster_info.tmp  = subset(Cluster_info,cluster==i)
    seurat_object.tmp = seurat_object[top20.tmp$gene,Cluster_info.tmp$barcode]
    seurat_object.tmp = data.frame(barcode=colnames(seurat_object.tmp),t(seurat_object.tmp@assays$RNA@scale.data))
    Cluster_info.tmp  = merge(Cluster_info.tmp[,-6],seurat_object.tmp)
    outname = paste0("Cluster",i,"_top20_gene_featureExp.xls")
    write.table(Cluster_info.tmp,file = file.path(OUT.ClusterAndDiff.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    outname = paste0("Cluster",i,"_top_gene_exp.pdf")
    p1 = VlnPlot(seurat_object, features = top20.tmp$gene,ncol = 4) + NoLegend()
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p1,width = 9,height =20)
    outname = paste0("Cluster",i,"_top_gene_feature_byTSNE.pdf")
    p2 = FeaturePlot(seurat_object, features = top20.tmp$gene,reduction = "tsne")  
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p2,width = 9,height =20)
    outname = paste0("Cluster",i,"_top_gene_feature_byUMAP.pdf")
    p3 = FeaturePlot(seurat_object, features = top20.tmp$gene,reduction = "umap")  
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p3,width = 9,height =20)
  }
  flage(pipeline.log,"4.Seurat.Cluster_and_diff")
}else{
  CATLOG("  Skip Cluster_and_diff")
}

# 4.celltype_annotation
OUT.celltypeAnnotation.DIR = file.path(OUT.Seurat.DIR[4],c("2.celltype_diff","3.celltype_marker_gene_plot"))
DIRCREATE(OUT.celltypeAnnotation.DIR)

# ref.rds = "//data/wqihao/database/refForScrna/singleR.DB/ImmGenData.rds"
if(!file.exists(file.path(OUT.DIR,"seurat_object.filted.celltype.RDS"))){
  ref.rds    = file.path(opts$DBdir,paste0(opts$annoDB,".rds"))
  ref_object = readRDS(ref.rds)
  sc_object  = SingleCellExperiment(assays=list(counts=seurat_object@assays$RNA@counts))
  sc_object  = scater::logNormCounts(sc_object)
  common     = intersect(rownames(sc_object), rownames(ref_object))
  ref_object = ref_object[common,]
  sc_object  = sc_object[common,]
  pred <- SingleR(test = sc_object, ref = ref_object, labels = ref_object$label.main)
  anno = data.frame(barcode=colnames(sc_object),
                    labels=pred$labels,
                    first.tuning.scores=pred$tuning.scores$first,
                    second.tuning.scores=pred$tuning.scores$second)
  seurat_object <- AddMetaData(seurat_object,metadata =pred$labels,col.name = "celltype" )
  Idents(object = seurat_object ) <- 'celltype'
  saveRDS(seurat_object,file = file.path(OUT.DIR,"seurat_object.filted.celltype.RDS"))
}else{
  seurat_object = readRDS(file = file.path(OUT.DIR,"seurat_object.filted.celltype.RDS"))
}

if(!file.exists(GetoptLong::qq("@{pipeline.log}/4.Seurat.celltype_annotation.flage"))){
  CATLOG("  Run celltype_annotation")
  p9  = DimPlot(seurat_object, reduction = "tsne",group.by = "celltype", label = F,pt.size = 0.5)
  p10 = DimPlot(seurat_object, reduction = "umap",group.by = "celltype", label = F, pt.size = 0.5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[4],"1.Celltype_byTSNE.pdf"),plot = p9,width = 6,height =5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[4],"1.Celltype_byUMAP.pdf"),plot = p10,width = 6,height =5)
  seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,min.cells.group = 1)
  cluster.averages <- AverageExpression(seurat_object)
  cluster.averages <- data.frame(gene=rownames(cluster.averages$RNA),cluster.averages$RNA,check.names = F)
  All_cluster_up_marker = data.frame(gene=rownames(seurat_object.markers),seurat_object.markers[,c(1:4,6)],check.names = F)
  All_cluster_up_marker = merge(All_cluster_up_marker,cluster.averages)
  write.table(All_cluster_up_marker,file = file.path(OUT.Seurat.DIR[4],"1.All_celltype_up_marker.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  top20 = seurat_object.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_log2FC)
  write.table(top20[,-5],file = file.path(OUT.Seurat.DIR[4],"4.Top20_celltype_pos_marker_forheatmap.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  # pdf(file = file.path(OUT.Seurat.DIR[4],"5.celltype_top_express_heatmap.pdf"),width = 10,height = 10,onefile = F)
  # dev.off()
  p = as.ggplot(Seurat::DoHeatmap(seurat_object, features = top20$gene,label=FALSE))
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[4],"5.celltype_top_express_heatmap.pdf"),plot = p,width = 10,height = 10)
  
  celltype_info = data.frame(barcode = colnames(seurat_object),seurat_object@reductions$umap@cell.embeddings,seurat_object@reductions$tsne@cell.embeddings,celltype=seurat_object@meta.data$celltype)
  write.table(celltype_info,file = file.path(OUT.Seurat.DIR[4],"celltype_info.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  celltype_num_stat = data.frame(table(seurat_object@meta.data$celltype),check.names = F)
  colnames(celltype_num_stat) = c("celltype","number")
  write.table(celltype_num_stat,file = file.path(OUT.Seurat.DIR[4],"celltype_num_stat.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  for (i in unique(seurat_object$celltype)) {
    cat(i)
    cluster.markers <- FindMarkers(seurat_object, ident.1 = i, min.pct = 0.25,logfc.threshold = 0.25,min.cells.group = 1)
    cluster.markers = data.frame(gene=rownames(cluster.markers),cluster.markers,check.names = F)
    cluster.markers = merge(cluster.markers,cluster.averages)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),".differential_and_annoation.xls")
    write.table(cluster.markers,file = file.path(OUT.celltypeAnnotation.DIR[1],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    top20.tmp = subset(top20,cluster==i)
    cluster.markers.top20 = subset(cluster.markers,gene%in%top20.tmp$gene)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top20_gene_diffInf.xls")
    write.table(cluster.markers.top20,file = file.path(OUT.celltypeAnnotation.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    celltype_info.tmp = subset(celltype_info,celltype==i)
    seurat_object.tmp = seurat_object[top20.tmp$gene,celltype_info.tmp$barcode]
    seurat_object.tmp = data.frame(barcode=colnames(seurat_object.tmp),t(seurat_object.tmp@assays$RNA@scale.data))
    celltype_info.tmp = merge(celltype_info.tmp[,-6],seurat_object.tmp)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top20_gene_featureExp.xls")
    write.table(celltype_info.tmp,file = file.path(OUT.celltypeAnnotation.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top_gene_exp.pdf")
    p1 = VlnPlot(seurat_object, features = top20.tmp$gene,ncol = 4) + NoLegend()
    ggplot2::ggsave(filename = file.path(OUT.celltypeAnnotation.DIR[2],outname),plot = p1,width = 9,height =20)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top_gene_feature_byTSNE.pdf")
    p2 = FeaturePlot(seurat_object, features = top20.tmp$gene,reduction = "tsne")  
    ggplot2::ggsave(filename = file.path(OUT.celltypeAnnotation.DIR[2],outname),plot = p2,width = 9,height =20)
    outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top_gene_feature_byUMAP.pdf")
    p3 = FeaturePlot(seurat_object, features = top20.tmp$gene,reduction = "umap")  
    ggplot2::ggsave(filename = file.path(OUT.celltypeAnnotation.DIR[2],outname),plot = p3,width = 9,height =20)
  }
  flage(pipeline.log,"4.Seurat.celltype_annotation")
}else{
  CATLOG("  Skip celltype_annotation")
}

# 5.Enrichment --------
CATLOG("5.Enrichment")
OUT.Enrichment.DIR = file.path(OUT.TREE.DIR[5],c("graphclust","seurat","celltype"))
DIRCREATE(OUT.Enrichment.DIR)
graphclust.diffexp = read.table(analysis.diffexp[1],sep = ",",header = T,stringsAsFactors = F)
clusters = colnames(graphclust.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()

# CATLOG(clusters)
if(!file.exists(GetoptLong::qq("@{pipeline.log}/5.Enrichment.graphclust.flage"))){
  CATLOG("  Run graphclust")
  for (cluster in clusters) {
    # filter top 20 genes
    cat(cluster)
    tmp.DIR = file.path(OUT.Enrichment.DIR[1],cluster)
    DIRCREATE(tmp.DIR)
    
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
      next
    }
    
    p.cols = grep(pattern = paste0(cluster,".Adjusted.p.value"),colnames(graphclust.diffexp))
    exp.cols = grep(pattern = paste0(cluster,".Mean.Counts"),colnames(graphclust.diffexp))
    fc.cols = grep(pattern = paste0(cluster,".Log2.fold.change"),colnames(graphclust.diffexp))
    de.genes = graphclust.diffexp[graphclust.diffexp[,fc.cols]>1 &graphclust.diffexp[,p.cols]<0.05,]
    de.genes.tmp = de.genes[,c(1,2,exp.cols,fc.cols,p.cols)]
    colnames(de.genes.tmp)  = c("ID","symbol","mean_counts","logFC","adj_pval")
    write.table(de.genes.tmp,file = file.path(tmp.DIR,paste0(cluster,".gene.txt")),quote = F,sep = "\t",row.names = F,col.names = T)
    
    if(species=="hsa"){
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      gene.df <- bitr(de.genes$Feature.ID, fromType = "ENSEMBL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Hs.eg.db)
      gene.df<- merge(gene.df,de.genes[,c(1,fc.cols)],by.x="ENSEMBL",by.y="Feature.ID")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      # GO
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        ego.BP <- enrichGO(gene          = gene.df$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           ont           = n,
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = TRUE) 
        
        # head(ego.BP.tab)
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(ego.BP, showCategory=30)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(ego.BP, foldChange=geneList)
        p5 <- cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        p8 = as.ggplot(~plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
        
        #KEGG 
        CATLOG("    Run KEGG")
        OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
        pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
        DIRCREATE(OUT.KEGG.DIT)
        DIRCREATE(pathviewDIR)
        kk <- enrichKEGG(gene         = names(geneList),
                         pAdjustMethod = "none",qvalueCutoff = 1,
                         organism     = 'hsa',
                         pvalueCutoff = 0.05)
        kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        kk.tab = kk@result
        kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
        write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(kk, showCategory=30)
        p2 <- clusterProfiler::dotplot(kk, showCategory=30)
        x2 <- pairwise_termsim(kk)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(kk, foldChange=geneList)
        p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
        
        for (ID in kk.tab$ID[1:5] ) {
          setwd(pathviewDIR)
          tryCatch(
            { pathview(gene.data  = geneList,
                       pathway.id = ID,
                       species    = "hsa",
                       kegg.dir = pathviewDIR,
                       limit      = list(gene=max(abs(geneList)), cpd=1))
            },
            warning = function(w) { cat("warning") },
            error = function(e) { "next" }, 
            finally = {cat("finished")}
          )
        }
        
        # Reactome
        CATLOG("    Run Reactome")
        OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
        DIRCREATE(OUT.Reactome.DIT)
        de <- names(geneList)
        x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "human", readable=TRUE)
        p1 <- barplot(x, showCategory=30)
        p2 <- clusterProfiler::dotplot(x, showCategory=30)
        x2 <- pairwise_termsim(x)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(x, foldChange=geneList)
        p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
        x.tab = x@result
        x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
        write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      }
    }
    
    if(species=="mmu"){
      suppressPackageStartupMessages(library(org.Mm.eg.db))
      
      gene.df <- bitr(de.genes$Feature.ID, fromType = "ENSEMBL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Mm.eg.db)
      gene.df<- merge(gene.df,de.genes[,c(1,fc.cols)],by.x="ENSEMBL",by.y="Feature.ID")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      # GO
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        cat(n)
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        if(file.exists(file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")))){
          next
        }
        
        ego.BP <- enrichGO(gene          = gene.df$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = n,
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = TRUE) 
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(ego.BP, showCategory=30)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        p4 <- clusterProfiler::heatplot(ego.BP, foldChange=geneList)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        p5 <- clusterProfiler::cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- clusterProfiler::cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        p8 = as.ggplot(~clusterProfiler::plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
        
      }
      
      
      #KEGG
      CATLOG("    Run KEGG")
      OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
      pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
      DIRCREATE(OUT.KEGG.DIT)
      DIRCREATE(pathviewDIR)
      kk <- enrichKEGG(gene         = names(geneList),
                       pAdjustMethod = "none",qvalueCutoff = 1,
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kk.tab = kk@result
      kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
      write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      p1 <- barplot(kk, showCategory=30)
      p2 <- clusterProfiler::dotplot(kk, showCategory=30)
      x2 <- pairwise_termsim(kk)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")

      p4 <- heatplot(kk, foldChange=geneList)
      p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
      
      
      for (ID in kk.tab$ID[1:5] ) {
        setwd(pathviewDIR)
        tryCatch(
          { pathview(gene.data  = geneList,
                     pathway.id = ID,
                     species    = "mmu",
                     kegg.dir = pathviewDIR,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
          },
          warning = function(w) { cat("warning") },
          error = function(e) { "next" }, 
          finally = {cat("finished")}
        )
      }
      
      # Reactome
      CATLOG("    Run Reactome")
      OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
      DIRCREATE(OUT.Reactome.DIT)
      de <- names(geneList)
      x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "mouse", readable=TRUE)
      p1 <- barplot(x, showCategory=30)
      p2 <- clusterProfiler::dotplot(x, showCategory=30)
      x2 <- pairwise_termsim(x)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(x, foldChange=geneList)
      p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
      x.tab = x@result
      x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
      write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
    }
    # flag for all finished -----------------
    system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
  }
  flage(pipeline.log,"5.Enrichment.graphclust")
  
}else{
  CATLOG("  Skip graphclust")
}

# SEURAT
if(!file.exists(GetoptLong::qq("@{pipeline.log}/5.Enrichment.SEURAT.flage"))){
  CATLOG("  Run SEURAT")
  for (i in unique(seurat_object$seurat_clusters)) {
    cat(i)
    tmp.DIR = file.path(OUT.Enrichment.DIR[2],paste0("Cluster.",i))
    DIRCREATE(tmp.DIR)
    
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/Cluster.@{i}.finished.flag"))){
      next
    }
    
    
    inputname = paste0("Cluster",i,".differential_and_annoation.xls")
    cluster.markers = read.table(file.path(OUT.ClusterAndDiff.DIR[1],inputname),header = T,check.names = F)
    cluster.markers.sig = subset(cluster.markers,p_val<0.05)
    
    if(species=="hsa"){
      library(org.Hs.eg.db)
      gene.df <- bitr(cluster.markers.sig$gene, fromType = "SYMBOL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Hs.eg.db)
      gene.df<- merge(gene.df,cluster.markers.sig[,c(1,3)],by.x="SYMBOL",by.y="gene")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      # GO
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        ego.BP <- enrichGO(gene          = gene.df$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           ont           = n,
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = TRUE) 
        
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(ego.BP, showCategory=30)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(ego.BP, foldChange=geneList)
        p5 <- cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        p8 = as.ggplot(~plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
      }
      
      # KEGG
      OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
      pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
      DIRCREATE(OUT.KEGG.DIT)
      DIRCREATE(pathviewDIR)
      kk <- enrichKEGG(gene         = names(geneList),
                       pAdjustMethod = "none",qvalueCutoff = 1,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kk.tab = kk@result
      kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
      write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      p1 <- barplot(kk, showCategory=30)
      p2 <- clusterProfiler::dotplot(kk, showCategory=30)
      x2 <- pairwise_termsim(kk)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(kk, foldChange=geneList)
      p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
      
      for (ID in kk.tab$ID[1:5]) {
        setwd(pathviewDIR)
        tryCatch(
          { pathview(gene.data  = geneList,
                     pathway.id = ID,
                     species    = "hsa",
                     kegg.dir = pathviewDIR,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
          },
          warning = function(w) { cat("warning") },
          error = function(e) { "next" }, 
          finally = {cat("finished")}
        )
      }
      
      OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
      DIRCREATE(OUT.Reactome.DIT)
      de <- names(geneList)
      x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "human", readable=TRUE)
      p1 <- barplot(x, showCategory=30)
      p2 <- clusterProfiler::dotplot(x, showCategory=30)
      x2 <- pairwise_termsim(x)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(x, foldChange=geneList)
      p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 12,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
      x.tab = x@result
      x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
      write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/Cluster.@{i}.finished.flag"))
    }
    
    if(species=="mmu"){
      library(org.Mm.eg.db)
      gene.df <- bitr(cluster.markers.sig$gene, fromType = "SYMBOL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Mm.eg.db)
      gene.df<- merge(gene.df,cluster.markers.sig[,c(1,3)],by.x="SYMBOL",by.y="gene")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      # GO
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        ego.BP <- enrichGO(gene          = gene.df$ENTREZID,
                           OrgDb         = org.Mm.eg.db,
                           ont           = n,
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = TRUE) 
        
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(ego.BP, showCategory=30)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(ego.BP, foldChange=geneList)
        p5 <- cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        p8 = as.ggplot(~plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
      }
      
      # KEGG
      OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
      pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
      DIRCREATE(OUT.KEGG.DIT)
      DIRCREATE(pathviewDIR)
      kk <- enrichKEGG(gene         = names(geneList),
                       pAdjustMethod = "none",qvalueCutoff = 1,
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kk.tab = kk@result
      kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
      write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      p1 <- barplot(kk, showCategory=30)
      p2 <- clusterProfiler::dotplot(kk, showCategory=30)
      x2 <- pairwise_termsim(kk)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(kk, foldChange=geneList)
      p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
      
      for (ID in kk.tab$ID[1:5]) {
        setwd(pathviewDIR)
        tryCatch(
          { pathview(gene.data  = geneList,
                     pathway.id = ID,
                     species    = "mmu",
                     kegg.dir = pathviewDIR,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
          },
          warning = function(w) { cat("warning") },
          error = function(e) { "next" }, 
          finally = {cat("finished")}
        )
      }
      
      OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
      DIRCREATE(OUT.Reactome.DIT)
      de <- names(geneList)
      x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "mouse", readable=TRUE)
      p1 <- barplot(x, showCategory=30)
      p2 <- clusterProfiler::dotplot(x, showCategory=30)
      x2 <- pairwise_termsim(x)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(x, foldChange=geneList)
      p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 12,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
      x.tab = x@result
      x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
      write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/Cluster.@{i}.finished.flag"))  
    }
    

  }
  
  flage(pipeline.log,"5.Enrichment.SEURAT")
  
}else{
  CATLOG("  Skip SEURAT")
}

#celltype
if(!file.exists(GetoptLong::qq("@{pipeline.log}/5.Enrichment.celltype.flage"))){
  CATLOG("  Run celltype")
  for (i in make.names(unique(seurat_object$celltype))) {
    print(i)
    tmp.DIR = file.path(OUT.Enrichment.DIR[3],make.names(i))
    DIRCREATE(tmp.DIR)
    
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/celltype.@{i}.finished.flag"))){
      next
    }
    
    inputname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),".differential_and_annoation.xls")
    cluster.markers = read.table(file.path(OUT.celltypeAnnotation.DIR[1],inputname),header = T,sep = "\t",check.names = F)
    cluster.markers.sig = subset(cluster.markers,p_val_adj<0.05)
    colnames(cluster.markers.sig) %<>% make.names()
    de.genes.tmp = cluster.markers.sig[,c(1,which(colnames(cluster.markers.sig)==i),3,6)]
    colnames(de.genes.tmp)  = c("symbol","mean_counts","logFC","adj_pval")
    write.table(de.genes.tmp,file = file.path(tmp.DIR,paste0("celltype.",i,".gene.txt")),quote = F,sep = "\t",row.names = F,col.names = T)
    
    
    if(species=="hsa"){
      library(org.Hs.eg.db)
      gene.df <- bitr(cluster.markers.sig$gene, fromType = "SYMBOL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Hs.eg.db) %>% invisible()
      # head(gene.df)
      gene.df<- merge(gene.df,cluster.markers.sig[,c(1,3)],by.x="SYMBOL",by.y="gene")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      
      # GO --------------
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        ego.BP <- enrichGO(gene          = gene.df$ENTREZID,
                           OrgDb         = org.Hs.eg.db,
                           ont           = n,
                           pAdjustMethod = "none",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 1,
                           readable      = TRUE) 
        
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        p1 <- barplot(ego.BP, showCategory=30)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(ego.BP, foldChange=geneList)
        p5 <- cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        p8 = ggplotify::as.ggplot(~plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
      }
      
      # KEGG -------------
      CATLOG("    Run KEGG")
      OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
      pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
      DIRCREATE(OUT.KEGG.DIT)
      DIRCREATE(pathviewDIR)
      kk <- enrichKEGG(gene         = names(geneList),
                       pAdjustMethod = "none",qvalueCutoff = 1,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kk.tab = kk@result
      kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
      write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      p1 <- barplot(kk, showCategory=30)
      p2 <- clusterProfiler::dotplot(kk, showCategory=30)
      x2 <- pairwise_termsim(kk)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(kk, foldChange=geneList)
      p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
      
      library(pathview)
      for (ID in kk.tab$ID[1:5]) {
        setwd(pathviewDIR)
        tryCatch(
          { pathview(gene.data  = geneList,
                     pathway.id = ID,
                     species    = "hsa",
                     kegg.dir = pathviewDIR,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
          },
          warning = function(w) { cat("warning") },
          error = function(e) { "next" }, 
          finally = {cat("finished\n")}
        )
      }
      CATLOG("    Run Reactome")
      OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
      DIRCREATE(OUT.Reactome.DIT)
      de <- names(geneList)
      x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "human", readable=TRUE)
      p1 <- barplot(x, showCategory=30)
      p2 <- clusterProfiler::dotplot(x, showCategory=30)
      x2 <- pairwise_termsim(x)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(x, foldChange=geneList)
      p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 12,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
      x.tab = x@result
      x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
      write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/celltype.@{i}.finished.flag"))
    }
    
    if(species=="mmu"){
      library(org.Mm.eg.db)
      gene.df <- bitr(cluster.markers.sig$gene, fromType = "SYMBOL",
                      toType = c("ENTREZID"),drop = F,
                      OrgDb = org.Mm.eg.db) 
      # head(gene.df)
      gene.df<- merge(gene.df,cluster.markers.sig[,c(1,3)],by.x="SYMBOL",by.y="gene")
      colnames(gene.df)[3] = "FC"
      gene.df = gene.df[!is.na(gene.df$ENTREZID),]
      geneList = structure(gene.df$FC,names=gene.df$ENTREZID)
      
      # GO --------------
      CATLOG("    Run GO")
      for (n in c("BP","CC","MF")) {
        DIRCREATE(file.path(tmp.DIR,"GO",n))
        ego.BP <- clusterProfiler::enrichGO(
          gene          = gene.df$ENTREZID,
          OrgDb         = org.Mm.eg.db,
          ont           = n,
          pAdjustMethod = "none",
          pvalueCutoff  = 0.05,
          qvalueCutoff  = 1,
          readable      = TRUE) 
        ego.BP.tab = ego.BP@result
        ego.BP.tab$Hyperlink = paste0("http://amigo.geneontology.org/amigo/term/",ego.BP.tab$ID)
        write.table(ego.BP.tab,file = file.path(tmp.DIR,"GO",n,"Enrichment_reuslt.xlsx"),quote = F,sep = "\t",row.names = F,col.names = T)
        
        p1 <- barplot(ego.BP, showCategory=30)
        p2 <- clusterProfiler::dotplot(ego.BP, showCategory=30)
        x2 <- pairwise_termsim(ego.BP)
        p3 <- clusterProfiler::emapplot(x2,layout="kk")
        p4 <- heatplot(ego.BP, foldChange=geneList)
        p5 <- cnetplot(ego.BP, categorySize="pvalue", foldChange=geneList)
        p6 <- cnetplot(ego.BP, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
        p7 <- cowplot::plot_grid(p5, p6, ncol=2)
        p8 <- ggplotify::as.ggplot(~plotGOgraph(ego.BP))
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".heatmap.pdf")),plot = p4,width = 20,height =8)
        ggplot2::ggsave(filename = file.path(tmp.DIR,"GO",n,paste0(n,".DAG.pdf")),plot = p8,width = 10,height =10)
      }
      
      # KEGG -------------
      CATLOG("    Run KEGG")
      OUT.KEGG.DIT =  file.path(tmp.DIR,"KEGG")
      pathviewDIR = file.path(OUT.KEGG.DIT,"pathview")
      DIRCREATE(OUT.KEGG.DIT)
      DIRCREATE(pathviewDIR)
      kk <- enrichKEGG(gene         = names(geneList),
                       pAdjustMethod = "none",qvalueCutoff = 1,
                       organism     = 'mmu',
                       pvalueCutoff = 0.05)
      kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kk.tab = kk@result
      kk.tab$Hyperlink = paste0("https://www.genome.jp/kegg-bin/show_pathway?",kk.tab$ID)
      write.table(kk.tab,file = file.path(tmp.DIR,"KEGG",paste0("KEGG","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      p1 <- barplot(kk, showCategory=30)
      p2 <- clusterProfiler::dotplot(kk, showCategory=30)
      x2 <- pairwise_termsim(kk)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(kk, foldChange=geneList)
      p5 <- cnetplot(kk, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Enrichment_Map.pdf")),plot = p3,width = 8,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"KEGG",paste0("KEGG",".heatmap.pdf")),plot = p4,width = 20,height =8)
      
      for (ID in kk.tab$ID[1:5]) {
        setwd(pathviewDIR)
        tryCatch(
          { pathview(gene.data  = geneList,
                     pathway.id = ID,
                     species    = "mmu",
                     kegg.dir = pathviewDIR,
                     limit      = list(gene=max(abs(geneList)), cpd=1))
          },
          warning = function(w) { cat("warning") },
          error = function(e) { "next" }, 
          finally = {cat("finished\n")}
        )
      }
      CATLOG("    Run Reactome")
      OUT.Reactome.DIT =  file.path(tmp.DIR,"Reactome")
      DIRCREATE(OUT.Reactome.DIT)
      de <- names(geneList)
      x <- enrichPathway(gene=de, pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 1,organism = "mouse", readable=TRUE)
      p1 <- barplot(x, showCategory=30)
      p2 <- clusterProfiler::dotplot(x, showCategory=30)
      x2 <- pairwise_termsim(x)
      p3 <- clusterProfiler::emapplot(x2,layout="kk")
      p4 <- heatplot(x, foldChange=geneList)
      p5 <- cnetplot(x, categorySize="pvalue", foldChange=geneList)
      p6 <- cnetplot(x, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
      p7 <- cowplot::plot_grid(p5, p6, ncol=2)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".barplot.p.adjust.pdf")),plot = p1,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".dotplot.p.adjust.pdf")),plot = p2,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Enrichment_Map.pdf")),plot = p3,width = 12,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".Gene-Term_Network.pdf")),plot = p7,width = 20,height =8)
      ggplot2::ggsave(filename = file.path(tmp.DIR,"Reactome",paste0("Reactome",".heatmap.pdf")),plot = p4,width = 20,height =8)
      x.tab = x@result
      x.tab$Hyperlink = paste0("http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=",x.tab$ID)
      write.table(x.tab,file = file.path(tmp.DIR,"Reactome",paste0("Reactome","_Enrichment_reuslt.xlsx")),quote = F,sep = "\t",row.names = F,col.names = T)
      
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/celltype.@{i}.finished.flag"))
    }
  }
  
  flage(pipeline.log,"5.Enrichment.celltype")
  
}else{
  CATLOG("  Skip celltype")
}

# 6.Convert --------
CATLOG("6.Convert pdf2png")
all.pdf = list.files(OUT.DIR,pattern = ".pdf$",recursive = T,full.names = T)
all.pdf.rm = grep(pattern = "Rplots.pdf",all.pdf,value = T)
file.remove(all.pdf.rm) %>% invisible()
all.pdf = grep(pattern = "Rplots.pdf",all.pdf,value = T,invert = T)
all.png = gsub(pattern = "pdf",replacement = "png",all.pdf)
# file.remove(all.png) %>% invisible()
# all.png.old =  gsub(pattern = "pdf",replacement = "ppi.png",all.pdf)
# file.rename(from = all.png,to = all.png.old) %>% invisible()

cmd = sprintf("/usr/bin/convert -density 300 %s %s",all.pdf,all.png)
Cmd(cmd = cmd,cl =16 ,out = all.png,step = "Convert pdf 2 png",f = system)


