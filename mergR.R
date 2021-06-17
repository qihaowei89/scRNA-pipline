#!/usr/bin/env Rscript 
source("/data/wqihao/Project.sc/scrna_analysis_prepare.R")
if(!interactive()){
  option_list <- list(
                      make_option(c("-i", "--inputDir"),type="character", help="Cellranger Results Directory"),
                      make_option(c("-o", "--outDir"),type="character", help="output Directory"),
                      make_option(c("-s", "--species"),type="character", help="species of single cell [hsa, mmu]"),
                      make_option(c("-a", "--annoDB"),type="character", help="DB use to anno celltype"),
                      make_option(c("-S", "--subdivide"),type="character", help="subdivide sample names"),
                      make_option(c("-d", "--DBdir"),type="character", default = "//data/wqihao/database/refForScrna/singleR.DB/",help="DB directory use to anno celltype"))
  opts <- parse_args(OptionParser(option_list=option_list))
}else{
  opts= list()
  opts$inputDir = "/data/wqihao/Project.sc/cellrangeRun/runDir/scRNA_20210511_aggr_GP1and2/"
  opts$outDir   = "/data/wqihao/Project.sc/testpipline/"
  opts$species  = "mmu"
  opts$annoDB   = "ImmGenData"
  opts$DBdir    = "/data/wqihao/database/refForScrna/singleR.DB/"
  opts$subdivide= "GP1,GP2" 
}
CATLOG("0.Parse args")
BIN = "/data/wqihao/Project.sc/workflow/bin/"

###### Check options 
if(length(opts) < 5){
  cat ("Use: Rscript  %prog -h see more Usage \n")
  cat ("Version: v1.0\n")
  cat ("Date: 2021-06-07\n")
  break
}

# input 
INPUT.DIR = opts$inputDir 
INPUT.OUTS.DIR = file.path(INPUT.DIR,"outs")

# output
OUT.DIR = opts$outDir
DIRCREATE(OUT.DIR)
OUT.RESTULE.DIR = file.path(OUT.DIR,"out","merged_samples")

pipeline.log = file.path(OUT.RESTULE.DIR,"pipeline.log")
DIRCREATE(pipeline.log)

subdivide.map = structure(strsplit(opts$subdivide,",")[[1]],names=c(1,2))
species = opts$species

setwd(OUT.DIR)
OUT.TREE.DIR = file.path(OUT.RESTULE.DIR,c("1.Basic_analysis","2.Cellranger_advanced","3.Seurat","4.CellAnnotation","5.Enrichment","6.GSEA","7.PPI","8.TFBS"))
DIRCREATE(OUT.TREE.DIR)

##### ------------------------------- 1.Basic_analysis --------------------------------
CATLOG("1.Basic_analysis")
OUT.Basic_analysis.DIR = file.path(OUT.TREE.DIR[1],c("1.1.raw_feature_bc_matrix","1.2.filtered_feature_bc_matrix","1.3.h5_files","1.4.BCR"))
DIRCREATE(OUT.Basic_analysis.DIR)
rawFeatureMatrix = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% grep(pattern = "raw_feature_bc_matrix",value = T) 
rawFeatureMatrix.files = list.files(rawFeatureMatrix,full.names = T)
rawH5 = list.files(INPUT.OUTS.DIR,pattern = "raw.*\\.h5$",recursive = T,full.names = T)

filteredFeatureMatrix = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% grep(pattern = "filtered_feature_bc_matrix",value = T)
filteredFeatureMatrix.files = list.files(filteredFeatureMatrix,full.names = T)

filteredH5 = list.files(INPUT.OUTS.DIR,pattern = "filtered.*\\.h5$",recursive = T,full.names = T)
SummaryJson = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "summary.json") 
aggregation.csv = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "aggregation.csv")
web_summary.html = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "web_summary.html")
bcfFiles  = GETFILE(input.dir = INPUT.OUTS.DIR, pattern = "(clonotypes.csv)|(consensus_annotations.csv)|(consensus.fasta)|(filtered_contig_annotations.csv)")

DIRCREATE(OUT.Basic_analysis.DIR[1])
DIRCREATE(OUT.Basic_analysis.DIR[2])
DIRCREATE(OUT.Basic_analysis.DIR[3])
DIRCREATE(OUT.Basic_analysis.DIR[4])

if(!file.exists(GetoptLong::qq("@{pipeline.log}/1.Basic_analysis.flage"))){
  COPYFILE(from.name = rawFeatureMatrix.files,to.name = OUT.Basic_analysis.DIR[1],step = "Copy rawFeatureMatrix")
  COPYFILE(from.name = filteredFeatureMatrix.files,to.name = OUT.Basic_analysis.DIR[2],step = "Copy filteredFeatureMatrix")

  COPYFILE(from.name = rawH5,to.name = OUT.Basic_analysis.DIR[3],step = "Copy rawH5")
  COPYFILE(from.name = filteredH5,to.name = file.path(OUT.Basic_analysis.DIR[3]),step = "Copy filteredH5")
  
  COPYFILE(from.name = bcfFiles,to.name = OUT.Basic_analysis.DIR[4],step = "Copy bcfFiles")
  COPYFILE(from.name = aggregation.csv,to.name =  OUT.TREE.DIR[1],step = "Copy aggregation.csv")
  COPYFILE(from.name = SummaryJson,to.name = OUT.TREE.DIR[1],step = "Copy SummaryJson")
  COPYFILE(from.name = web_summary.html,to.name = OUT.TREE.DIR[1],step = "Copy web_summary.html")
  
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
  p3 = VlnPlot(seurat_object,features = c("nGene", "nUMI"), ncol = 2,pt.size = 0.4) 
  p4 = FeatureScatter(seurat_object, feature1 = "nUMI", feature2 = "nGene",cols = "black") + ggplot2::theme(legend.title=ggplot2::element_blank(),legend.text = ggplot2::element_blank())
  nGene.histogram = file.path(OUT.TREE.DIR[1],"nGene.histogram.pdf")
  nUMI.histogram = file.path(OUT.TREE.DIR[1],"nUMI.histogram.pdf")
  nGene_and_nUMI_correlation = file.path(OUT.TREE.DIR[1],"nGene_and_nUMI_correlation.pdf")
  nGene_and_nUMI_distribution = file.path(OUT.TREE.DIR[1],"nGene_and_nUMI_distribution.pdf")
  ggplot2::ggsave(filename = nGene.histogram,plot =p1 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nUMI.histogram,plot =p2 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nGene_and_nUMI_correlation,plot =p4 ,width = 6,height = 6,dpi = "retina")
  ggplot2::ggsave(filename = nGene_and_nUMI_distribution,plot =p3 ,width = 6,height = 6,dpi = "retina")
  flage(pipeline.log,"1.Basic_analysis")
}


##### ------------------------------- 2.Cellranger_advanced --------------------------------
CATLOG("2.Cellranger_advanced")
OUT.Cellranger_advanced.DIR = file.path(OUT.TREE.DIR[2],c("cloupe","Cluster_specific_genes","correlation","pca_cluster","tsne_cluster","tsne_subdivide","umap_cluser","umap_subdivide"))

# cloupe
cloupeFiles = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "cloupe")
vloupeFiles = GETFILE(input.dir = INPUT.OUTS.DIR,pattern = "vloupe")
DIRCREATE(file.path(OUT.Cellranger_advanced.DIR[1],c("count","vdj_b")))
# OUT.cloupeFiles = file.path(FDIR(cloupeFiles),paste0(basename(cloupeFiles))) %>% file.path(OUT.Cellranger_advanced.DIR[1],.)
# OUT.vloupeFiles = file.path(FDIR(vloupeFiles),paste0(basename(vloupeFiles))) %>% file.path(OUT.Cellranger_advanced.DIR[1],.)
COPYFILE(from.name = cloupeFiles,to.name = file.path(OUT.Cellranger_advanced.DIR[1],"count"),step = "Copy cloupeFiles")
COPYFILE(from.name = vloupeFiles,to.name = file.path(OUT.Cellranger_advanced.DIR[1],"vdj_b"),step = "Copy vloupeFiles")

# Cluster_specific_genes
analysisDIR = list.dirs(INPUT.OUTS.DIR,recursive = T,full.names = T) %>% grep(pattern = "analysis",value = T) 
analysis.TSNE = analysisDIR %>% grep(pattern = "tsne/",value = T) %>% list.files(full.names = T)
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
analysis.diffexp = analysisDIR %>% grep(pattern = "diffexp/",value = T) %>% list.files(full.names = T)
seurat_object = CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),project = "cellranger")
seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
log2.gene.exp = log2(seurat_object@assays$RNA@counts+1)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.Cluster_specific_genes.flage"))){
  CATLOG("  Run Cluster_specific_genes")
  for (n in 1:10) {
    clusterType = basename(dirname(analysis.clusters[n]))
    CATLOG(clusterType)
    tmp.DIR = file.path(OUT.Cellranger_advanced.DIR[2],clusterType)
    DIRCREATE(tmp.DIR)
    
    if(file.exists(file.path(tmp.DIR,paste0(clusterType,".diffexp.volvano.pdf")))){
      next
    }
    
    tmp = read.table(analysis.TSNE,sep = ",",header = T,stringsAsFactors = F)
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    tmp.diffexp = read.table(analysis.diffexp[n],sep = ",",header = T,stringsAsFactors = F)
    
    clusters = colnames(tmp.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()
    COPYFILE(from.name = analysis.diffexp[n],to.name = tmp.DIR,step = "Copy diffexp")
    
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
    p = ggplotify::as.ggplot(~ComplexHeatmap::draw(ht))
    outname = paste0(clusterType,".pheatmap.pdf")
    ggplot2::ggsave(filename = file.path(tmp.DIR,outname),plot =  p ,width = 20,height =6 )
    plot.list2 = list()
    for (cluster in clusters) {
      CATLOG(cluster)
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
      plot.list[[1]] = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp.clusters$Cluster,type = "Cluster",title ="t-SNE")
      plot.list2[[cluster]] = volcanoPlot(logFC = tmp.diffexp[,fc.cols],pValue =tmp.diffexp[,p.cols] ,FC = 3,main = paste0(clusterType,":",cluster),size = 0.3)
      outname = paste0(clusterType,".",gsub(pattern = "\\.",replacement = "_",cluster),".pdf")
      for (gene in top20genes) {
        plot.list[[gene]] = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = as.numeric(log2.gene.exp[gene,]),
                                      title = paste(c(as.character(tmp.diffexp[tmp.diffexp$Feature.Name==gene,1:2]),"Log2 Total UMI Count Per Cell"),collapse = ","),type = "exp")
      }
      ncols = 3
      pdf(file = file.path(tmp.DIR,outname),width = 25,height =40,onefile = F)
      scater::multiplot(plotlist = plot.list,cols = ncols)
      dev.off()
    }
    
    outname = paste0(clusterType,".diffexp.volvano.pdf")
    if(length(plot.list2)==0){
      write.table(x = outname,file = file.path(OUT.RESTULE.DIR,"err.txt"),append = T,quote = F,col.names = F)
      next
    }
    ncols = 2
    pdf(file = file.path(tmp.DIR,outname),width = 7,height = 16,onefile = F)
    scater::multiplot(plotlist = plot.list2,cols = ncols)
    dev.off()
    
  }
  flage(pipeline.log,"2.Cellranger_advanced.Cluster_specific_genes")
  
}else{
  CATLOG("  Skip Cluster_specific_genes")
}

# correlation 
DIRCREATE(OUT.Cellranger_advanced.DIR[3])
if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.correlation.flage"))){
  CATLOG("  Run correlation")
  for (n in 1:10) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    if(file.exists(file.path(OUT.Cellranger_advanced.DIR[3],paste0(clusterType,".subpopulation.correlation.pdf")))){
      next
    }
    
    tmp.diffexp = read.table(analysis.diffexp[n],sep = ",",header = T,stringsAsFactors = F)
    
    tmp.Mean.Counts = tmp.diffexp %>% tidyr::unite("ID", Feature.ID:Feature.Name,sep = ":",remove = T) %>% 
      dplyr::select(matches("(ID)|(Mean.Counts)"))
    
    corr <- cor(apply(t(tmp.Mean.Counts[,-1]),1,scale), method = 'pearson')
    p <-  pheatmap::pheatmap(corr,cluster_cols = F,cluster_rows = F,display_numbers = TRUE,
                             color = colorRampPalette(c("yellow","red"))(100))
    outname = paste0(clusterType,".","subpopulation.correlation",".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[3],outname),plot =  p ,width = 6,height =6 )
  }
  flage(pipeline.log,"2.Cellranger_advanced.correlation")
}else{
  CATLOG("  Skip correlation")
}

# pca 
DIRCREATE(OUT.Cellranger_advanced.DIR[4])
analysis.PCA = analysisDIR %>% grep(pattern = "pca/",value = T) 
analysis.PCA.clusters =  GETFILE(input.dir = analysis.PCA,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.PCA.clusters,sep = ",",header = T,stringsAsFactors = F)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.pca.flage"))){
  CATLOG("  Run pca")
  for (n in 1:10) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    if(file.exists(file.path(OUT.Cellranger_advanced.DIR[4],paste0("pca.",clusterType,".pdf")))){
      next
    }
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$PC.1 ,tsneY = tmp$PC.2,Clusters = tmp.clusters$Cluster,type = "Cluster",xlab = "PC1",ylab = "PC2",title =paste0(clusterType,"\n","PCA projection of Cells Colored by Autometed Clustering"))
    outname = paste0("pca.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[4],outname),plot = p,width = 6,height =6 )
    
  }
  flage(pipeline.log,"2.Cellranger_advanced.pca")
}else{
  CATLOG("  Skip pca")
}


# tsne -----
DIRCREATE(OUT.Cellranger_advanced.DIR[5])
analysis.TSNE = analysisDIR %>% grep(pattern = "tsne/",value = T) 
analysis.TSNE.clusters =  GETFILE(input.dir = analysis.TSNE,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.TSNE.clusters,sep = ",",header = T,stringsAsFactors = F)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.tsne.flage"))){
  CATLOG("  Run tsne")
  for (n in 1:10) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp.clusters$Cluster,type = "Cluster",xlab = "TSNE.1",ylab = "TSNE.2",title =paste0(clusterType,"\n","t-SNE projection of Cells Colored by Autometed Clustering"))
    outname = paste0("tsne.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[5],outname),plot = p,width = 6,height =6 )
  }
  p = PLOT.TSNE2(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = seurat_object$nFeature_RNA,type = "exp",
                 xlab = "TSNE.1",ylab = "TSNE.2",title =paste0( "t-SNE Projection of Cells Colored by Gene Counts"),labs = "Gene Counts")
  outname = paste0("tsne.nGene",".pdf")
  ggplot2::ggsave(filename = file.path(OUT.TREE.DIR[2],outname),plot = p,width = 6,height =6 )
  p = PLOT.TSNE2(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = seurat_object$nCount_RNA,type = "exp",
                 xlab = "TSNE.1",ylab = "TSNE.2",title =paste0( "t-SNE Projection of Cells Colored by UMI Counts"),labs = "UMI Counts")
  outname = paste0("tsne.nUMI",".pdf")
  ggplot2::ggsave(filename = file.path(OUT.TREE.DIR[2],outname),plot = p,width = 6,height =6 )
  
  flage(pipeline.log,"2.Cellranger_advanced.tsne")
}else{
  CATLOG("  Skip tsne")
}


# tsne subdivide -----
DIRCREATE(OUT.Cellranger_advanced.DIR[6])
analysis.TSNE = analysisDIR %>% grep(pattern = "tsne/",value = T) 
analysis.TSNE.clusters =  GETFILE(input.dir = analysis.TSNE,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.TSNE.clusters,sep = ",",header = T,stringsAsFactors = F)
tmp$subdivide.ID = gsub(pattern = ".*-(\\d{1,2})",replacement = "\\1",tmp$Barcode)  
tmp$subdivide = subdivide.map[tmp$subdivide.ID]
p = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp$subdivide,size = 0.5,
              type = "Cluster",xlab = "TSNE.1",ylab = "TSNE.2",
              title = "colored by sample")
outname =  "subdivide.sample.pdf"
ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[6],outname),plot = p,width = 6,height =6 )

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.tsne.subdivide.flage"))){
  CATLOG("  Run tsne subdivide")
  for (n in 1:10) {
    # n=1
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    tmp.clusters = merge(tmp.clusters,tmp)
    tmp.clusters$ClusterAndSample = paste(tmp.clusters$Cluster,tmp.clusters$subdivide,sep = "+")
    
    
    tmp.subdivided = tmp.clusters[,c(1,7)]
    colnames(tmp.subdivided)[2] = paste0(clusterType,"-","subdivided")
    # head(tmp.subdivided)
    outname = paste0(clusterType,".subdivided.csv")
    write.table(tmp.clusters,file = file.path(OUT.Cellranger_advanced.DIR[6],outname),quote = F,sep = ",",col.names = T,row.names = F)
    
    tmp.stat = tmp.clusters[,c(2,6)] %>% 
      dplyr::group_by(subdivide,Cluster) %>% 
      dplyr::mutate(nCells=n()) %>% 
      dplyr::arrange(Cluster) %>% unique() %>%
      dplyr::select(subdivide,nCells,Cluster) %>%
      dplyr::rename(Sample=subdivide)
    outname = paste0(clusterType,".stat.xls")
    write.table(tmp.stat,file = file.path(OUT.Cellranger_advanced.DIR[6],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    
    p = PLOT.TSNE(tsneX =tmp$TSNE.1 ,tsneY = tmp$TSNE.2,Clusters = tmp.clusters$ClusterAndSample,
                  type = "Cluster",xlab = "TSNE.1",ylab = "TSNE.2",
                  title =paste0(clusterType, ": colored by cluster and sample"))
    outname = paste0(clusterType,".subdivide.pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[6],outname),plot = p,width = 6,height =6 )
    
    tmp.clusters$type = clusterType
    tmp.list = list(tmp.clusters)
    # dim(tmp.clusters)
    for(i in subdivide.map){
      # i = subdivide.map[2]
      sub.tmp = subset(tmp.clusters,subdivide==i)
      dim(sub.tmp)
      sub.tmp$type = i
      tmp.list = append(tmp.list,list(sub.tmp))
    }
    length.tmp.list = length(tmp.list)
    tmp.list = do.call(rbind,tmp.list)
    tmp.list$type = factor(tmp.list$type,levels = c(clusterType,subdivide.map))
    outname = paste0(clusterType,".subdivide.samples.pdf")
    l = 4
    if(length.tmp.list%%2==0){
      p = PLOT.TSNE3(tsneX =tmp.list$TSNE.1 ,tsneY = tmp.list$TSNE.2,Clusters = tmp.list$Cluster,
                     type = tmp.list$type,xlab = "TSNE.1",ylab = "TSNE.2",size = 0.2,
                     title ="")  
      p = p +  ggplot2::facet_wrap(vars(type),ncol = 2)
      ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[6],outname),plot = p,height = l,width =l*3*ceiling(length.tmp.list/2))
    }else{
      p = PLOT.TSNE3(tsneX =tmp.list$TSNE.1 ,tsneY = tmp.list$TSNE.2,Clusters = tmp.list$Cluster,
                     type = tmp.list$type,xlab = "TSNE.1",ylab = "TSNE.2",size = 0.2,
                     title ="")  
      p = p + ggplot2::facet_wrap(vars(type),ncol = 3)
      ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[6],outname),plot = p,height = l,width =l*3*ceiling(length.tmp.list/3))
    }
    
  }
  
  flage(pipeline.log,"2.Cellranger_advanced.tsne.subdivide")
}else{
  CATLOG("  Skip tsne subdivide")
}


# umap -----
DIRCREATE(OUT.Cellranger_advanced.DIR[7])
analysis.UMAP = analysisDIR %>% grep(pattern = "umap/",value = T) 
analysis.UMAP.clusters =  GETFILE(input.dir = analysis.UMAP,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.UMAP.clusters,sep = ",",header = T,stringsAsFactors = F)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.umap.flage"))){
  CATLOG("  Run umap")
  
  for (n in 1:10) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    p = PLOT.TSNE(tsneX =tmp$UMAP.1 ,tsneY = tmp$UMAP.2,Clusters = tmp.clusters$Cluster,type = "Cluster",xlab = "UMAP.1",ylab = "UMAP.2",title =paste0(clusterType,"\n","UMAP projection of Cells Colored by Autometed Clustering"))
    outname = paste0("umap.",clusterType,".pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[7],outname),plot = p,width = 6,height =6 )
  }
  
  flage(pipeline.log,"2.Cellranger_advanced.umap")
}else{
  CATLOG("  Skip umap")
}

# umap subdivide ----
DIRCREATE(OUT.Cellranger_advanced.DIR[8])
subdivide.map = structure(c("GP1","GP2"),names=c(1,2))
analysis.UMAP = analysisDIR %>% grep(pattern = "umap/",value = T) 
analysis.UMAP.clusters =  GETFILE(input.dir = analysis.UMAP,pattern = "projection.csv")
analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
tmp = read.table(analysis.UMAP.clusters,sep = ",",header = T,stringsAsFactors = F)
tmp$subdivide.ID = gsub(pattern = ".*-(\\d{1,2})",replacement = "\\1",tmp$Barcode)  
tmp$subdivide = subdivide.map[tmp$subdivide.ID]
p = PLOT.TSNE(tsneX =tmp$UMAP.1 ,tsneY = tmp$UMAP.2,Clusters = tmp$subdivide,size = 0.5,
              type = "Cluster",xlab = "UMAP.1",ylab = "UMAP.2",
              title = "colored by sample")
outname =  "subdivide.sample.pdf"
ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[8],outname),plot = p,width = 6,height =6 )


if(!file.exists(GetoptLong::qq("@{pipeline.log}/2.Cellranger_advanced.umap.subdivide.flage"))){
  CATLOG("  Run umap subdivide")
  for (n in 1:10) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    tmp.clusters = merge(tmp.clusters,tmp)
    tmp.clusters$ClusterAndSample = paste(tmp.clusters$Cluster,tmp.clusters$subdivide,sep = "+")
    # head(tmp.clusters)
    
    tmp.subdivided = tmp.clusters[,c(1,7)]
    colnames(tmp.subdivided)[2] = paste0(clusterType,"-","subdivided")
    # head(tmp.subdivided)
    outname = paste0(clusterType,".subdivided.csv")
    write.table(tmp.clusters,file = file.path(OUT.Cellranger_advanced.DIR[8],outname),quote = F,sep = ",",col.names = T,row.names = F)
    
    tmp.stat = tmp.clusters[,c(2,6)] %>% 
      dplyr::group_by(subdivide,Cluster) %>% 
      dplyr::mutate(nCells=n()) %>% 
      dplyr::arrange(Cluster) %>% unique() %>%
      dplyr::select(subdivide,nCells,Cluster) %>%
      dplyr::rename(Sample=subdivide)
    outname = paste0(clusterType,".stat.xls")
    write.table(tmp.stat,file = file.path(OUT.Cellranger_advanced.DIR[8],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    
    
    
    p = PLOT.TSNE(tsneX =tmp$UMAP.1 ,tsneY = tmp$UMAP.2,Clusters = tmp.clusters$ClusterAndSample,
                  type = "Cluster",xlab = "UMAP.1",ylab = "UMAP.2",
                  title =paste0(clusterType, ": colored by cluster and sample"))
    outname = paste0(clusterType,".subdivide.pdf")
    ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[8],outname),plot = p,width = 6,height =6 )
    tmp.clusters$type = clusterType
    tmp.list = list(tmp.clusters)
    for(i in subdivide.map){
      sub.tmp = subset(tmp.clusters,subdivide==i)
      dim(sub.tmp)
      sub.tmp$type = i
      tmp.list = append(tmp.list,list(sub.tmp))
    }
    length.tmp.list = length(tmp.list)
    tmp.list = do.call(rbind,tmp.list)
    tmp.list$type = factor(tmp.list$type,levels = c(clusterType,subdivide.map))
    outname = paste0(clusterType,".subdivide.samples.pdf")
    l = 4
    if(length.tmp.list%%2==0){
      p = PLOT.TSNE3(tsneX =tmp.list$UMAP.1 ,tsneY = tmp.list$UMAP.2,Clusters = tmp.list$Cluster,
                     type = tmp.list$type,xlab = "UMAP.1",ylab = "UMAP.2",size = 0.2,
                     title ="")  
      p = p +  ggplot2::facet_wrap(vars(type),ncol = 2)
      ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[8],outname),plot = p,height = l,width =l*3*ceiling(length.tmp.list/2))
    }else{
      p = PLOT.TSNE3(tsneX =tmp.list$UMAP.1 ,tsneY = tmp.list$UMAP.2,Clusters = tmp.list$Cluster,
                     type = tmp.list$type,xlab = "UMAP.1",ylab = "UMAP.2",size = 0.2,
                     title ="")  
      p = p + ggplot2::facet_wrap(vars(type),ncol = 3)
      
      ggplot2::ggsave(filename = file.path(OUT.Cellranger_advanced.DIR[8],outname),plot = p,height = l,width =l*3*ceiling(length.tmp.list/3))
    }
  }
  
  flage(pipeline.log,"2.Cellranger_advanced.umap.subdivide")
}else{
  CATLOG("  Skip umap subdivide")
}

if(!file.exists(file.path(OUT.TREE.DIR[2],"summary_cell.pca_tsne_clustering.csv"))){
  CATLOG("  Run summary ")
  analysis.clusters = analysisDIR %>% grep(pattern = "clustering/",value = T) %>% list.files(full.names = T)
  stat.cell = function(n) {
    CATLOG(n)
    clusterType = basename(dirname(analysis.clusters[n]))
    
    tmp.clusters = read.table(analysis.clusters[n],sep = ",",header = T,stringsAsFactors = F)
    res = data.frame(Cluster.Method = clusterType,table(tmp.clusters$Cluster),stringsAsFactors = F,check.names = F)
    colnames(res)[2:3] = c("Cluster.ID","Cell_Counts")
    return(res)
  }
  nCell_per_cluster = do.call(rbind,lapply(1:10, stat.cell))
  write.table(nCell_per_cluster,file = file.path(OUT.TREE.DIR[2],"nCell_per_cluster.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  
  f.read <- function(n){
    clusterType = basename(dirname(n))
    res = data.table::fread(n)
    colnames(res)[2] = clusterType
    return(res)
  }
  
  tmp.TSNE = read.table(analysis.TSNE.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.UMP = read.table(analysis.UMAP.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.PCA = read.table(analysis.PCA.clusters,sep = ",",header = T,stringsAsFactors = F)
  tmp.seurat = data.frame(Barcode=paste0(rownames(seurat_object@meta.data)),nGene=seurat_object@meta.data$nFeature_RNA,nUMI=seurat_object@meta.data$nCount_RNA,mt_percent=seurat_object@meta.data$percent.mt)
  
  res.list = lapply(analysis.clusters, f.read)
  res.list = append(res.list,list(tmp.TSNE))
  res.list = append(res.list,list(tmp.UMP))
  res.list = append(res.list,list(tmp.PCA))
  res.list = append(res.list,list(tmp.seurat))
  
  summary_cell.pca_tsne_clustering = purrr::reduce(res.list,merge)
  write.table(summary_cell.pca_tsne_clustering,file = file.path(OUT.TREE.DIR[2],"summary_cell.pca_tsne_clustering.csv"),quote = F,sep = ",",col.names = T,row.names = F)
  
}else{
  CATLOG("  Skip summary")
}


##### ------------------------------- 3.seurat --------------------------------
CATLOG("3.seurat")
# OUT.Seurat.DIR = file.path(OUT.TREE.DIR[3],c("1.Expression","2.PCA_analysis","3.Cluster_and_diff","4.Conserved_analysis","5.Sample_and_diff","6.ClusterSample_and_diff"))
OUT.Seurat.DIR = file.path(OUT.TREE.DIR[3],c("1.Expression","2.PCA_analysis","3.Cluster_and_diff"))
DIRCREATE(OUT.Seurat.DIR)

npcs=30

if(!file.exists(file.path(OUT.TREE.DIR[3],"seurat_object.filted.RDS"))){
  seurat_object = CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),min.cells = 3,min.features = 200,project = "cellranger")
  seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
  Max.nUMI  = stats::quantile(seurat_object$nCount_RNA,probs=c(0.99))  %>%  as.numeric()
  seurat_object <- subset(seurat_object, percent.mt <= 25&nFeature_RNA>200& nCount_RNA < Max.nUMI)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(object = seurat_object,selection.method = "vst")
  seurat_object <- ScaleData(seurat_object,features = rownames(seurat_object))     
  seurat_object <- RunPCA(object= seurat_object,pc.genes=VariableFeatures(object = seurat_object),npcs = npcs,verbose = F)  
  seurat_object <- JackStraw(seurat_object, num.replicate = 100,dims = npcs)
  seurat_object <- ScoreJackStraw(seurat_object, dims = 1:npcs)
  pcsNumber = which(seurat_object@reductions$pca@jackstraw@overall.p.values[,2] < 0.01)
  seurat_object <- FindNeighbors(seurat_object, dims = pcsNumber)
  seurat_object <- FindClusters(seurat_object, resolution = 1)
  seurat_object <- RunUMAP(seurat_object, dims = pcsNumber)
  seurat_object <- RunTSNE(seurat_object, dims = pcsNumber)
  saveRDS(seurat_object,file = file.path(OUT.TREE.DIR[3],"seurat_object.filted.RDS"))
}else{
  seurat_object <- readRDS(file.path(OUT.TREE.DIR[3],"seurat_object.filted.RDS"))  
}


# 1.Expression-----------
Intergrated_matrix.DIR = file.path(OUT.Seurat.DIR[1],"1.Intergrated_matrix")
DIRCREATE(Intergrated_matrix.DIR)

if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Seurat.Expression.flage"))){
  CATLOG("  Run Expression")
  Intergrated_matrix = seurat_object@assays$RNA@counts %>% as.matrix()
  write.csv(Intergrated_matrix,file = file.path(Intergrated_matrix.DIR,"gene_expression.csv"),quote = F,sep = ",",row.names = T,col.names = T)
  p1 = VlnPlot(seurat_object,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "orig.ident",ncol = 3,pt.size = 0.1)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[1],"2.nGene_nUMI_percent.mito.pdf"),plot = p1,width = 6,height =5)
  p2. <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3. <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  P = p2.+p3.
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[1],"3.Percent.mito_and_nGene_associate_with_nUMI.pdf"),plot = P,width = 10,height =5)
  
  flage(pipeline.log,"3.Seurat.Expression")
}else{
  CATLOG("  Skip Expression")
}

# 2.pca_analysis ----------
if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Seurat.pca_analysis.flage"))){
  CATLOG("  Run pca_analysis")
  pca50_vargene = seurat_object@reductions$pca@feature.loadings
  pca50_vargene = data.frame(Gene = rownames(pca50_vargene),pca50_vargene)
  write.table(pca50_vargene,file = file.path(OUT.Seurat.DIR[2],sprintf("1.pca%s_vargene.xls",npcs)),quote = F,sep = "\t",col.names = T,row.names = F)
  sink(file.path(OUT.Seurat.DIR[2],sprintf("2.PC%s_top20_sep.txt",npcs)))
  print(seurat_object@reductions$pca, dims = 1:npcs, nfeatures = 20)
  sink()
  p4 = VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca",ncol = 2)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],sprintf("3.Top%s_gene_with_PC1_PC2.pdf",npcs)),plot = p4,width = 6,height =5)
  p5 = ggplotify::as.ggplot(~DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE))
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],sprintf("4.Top%s_gene_PC1_heatmap.pdf",npcs)),plot = p5,width = 6,height =5)
  pdf(file = file.path(OUT.Seurat.DIR[2],sprintf("5.PC%s_heatmap.pdf",npcs)),width = 10,height = 25,onefile = F)
  DimHeatmap(seurat_object, dims = 1:npcs,cells = 500,ncol = 5)
  dev.off()
  pcsNumber = which(seurat_object@reductions$pca@jackstraw@overall.p.values[,2] < 0.01)
  p7 = JackStrawPlot(seurat_object,reduction = "pca",dims = pcsNumber)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],"6.PC_significant.pdf"),plot = p7,width = 9,height =8)
  p8 = ElbowPlot(seurat_object,ndims = npcs)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[2],sprintf("7.PC%s_sd.pdf",npcs)),plot = p8,width = 8,height =8)
  
  flage(pipeline.log,"3.Seurat.pca_analysis")
}else{
  CATLOG("  Skip pca_analysis")
}

# 3.Cluster_and_diff ----------
OUT.ClusterAndDiff.DIR = file.path(OUT.Seurat.DIR[3],c("2.Cluster_diff","3.Cluster_marker_gene_plot"))
DIRCREATE(OUT.ClusterAndDiff.DIR)

subdivide = colnames(seurat_object) %>% gsub(pattern = ".*-(\\d{1,2})",replacement = "\\1")
subdivide = subdivide.map[subdivide]
# group = subdivide 

if(!file.exists(GetoptLong::qq("@{pipeline.log}/3.Seurat.Cluster_and_diff.flage"))){
  CATLOG("  Run Cluster_and_diff")
  seurat_object = AddMetaData(seurat_object,metadata = subdivide,col.name = "subdivide")
  # seurat_object = AddMetaData(seurat_object,metadata = group,col.name = "group")
  
  pt.size=0.3
  p9  = DimPlot(seurat_object, reduction = "tsne", label = T,pt.size = pt.size)
  p10 = DimPlot(seurat_object, reduction = "umap",label = T, pt.size = pt.size)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_byTSNE.pdf"),plot = p9,width = 6,height =5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_byUMAP.pdf"),plot = p10,width = 6,height =5)
  
  
  # subdivide 
  n.subdivide = length(subdivide.map)
  p11  = DimPlot(seurat_object, reduction = "tsne", label = T,pt.size = 0.3,split.by ="subdivide")
  p11 = p9+p11 + plot_layout(widths = c(1, n.subdivide))
  p12 = DimPlot(seurat_object, reduction = "umap",label = T, pt.size = 0.3,split.by ="subdivide")
  p12 = p10+p12  + plot_layout(widths = c(1, n.subdivide))
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_with_sample_byTSNE.pdf"),plot = p11,width = 5+5*n.subdivide,height =5)
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_with_sample_byUMAP.pdf"),plot = p12,width = 5+5*n.subdivide,height =5)
  
  # # group
  # n.group = length(group.map)
  # p11  = DimPlot(seurat_object, reduction = "tsne", label = T,pt.size = 0.3,split.by ="subdivide")
  # p11 = p9+p11 + plot_layout(widths = c(1, n.subdivide))
  # p12 = DimPlot(seurat_object, reduction = "umap",label = T, pt.size = 0.3,split.by ="subdivide")
  # p12 = p10+p12  + plot_layout(widths = c(1, n.subdivide))
  # ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_with_group_byTSNE.pdf"),plot = p11,width = 5+5*n.subdivide,height =5)
  # ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR[3],"Cell_cluster_with_group_byUMAP.pdf"),plot = p12,width = 5+5*n.subdivide,height =5)
  
  if(!file.exists(file.path(OUT.TREE.DIR[3],"seurat_object.markers.RDS"))){
    seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    saveRDS(seurat_object.markers,file = file.path(OUT.TREE.DIR[3],"seurat_object.markers.RDS"))
  }else{
    seurat_object.markers = readRDS(file.path(OUT.TREE.DIR[3],"seurat_object.markers.RDS"))
  }
  
  cluster.averages <- AverageExpression(seurat_object)
  cluster.averages <- data.frame(gene=rownames(cluster.averages$RNA),cluster.averages$RNA,check.names = F)
  
  All_cluster_up_marker = data.frame(gene=rownames(seurat_object.markers),seurat_object.markers[,c(1:4,6)],check.names = F)
  All_cluster_up_marker = merge(All_cluster_up_marker,cluster.averages)
  All_cluster_up_marker = All_cluster_up_marker %>% dplyr::arrange(cluster,desc(cluster))
  write.table(All_cluster_up_marker,file = file.path(OUT.Seurat.DIR[3],"1.All_cluster_up_marker.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  
  topN=10
  
  top20 = seurat_object.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = topN, wt = avg_log2FC)
  write.table(top20[,-5],file = file.path(OUT.Seurat.DIR[3],sprintf("4.Top%s_cluster_pos_marker_forheatmap.xls",topN)),quote = F,sep = "\t",col.names = T,row.names = F)
  
  Cluster_info = data.frame(barcode = colnames(seurat_object),
                            seurat_object@reductions$umap@cell.embeddings,
                            seurat_object@reductions$tsne@cell.embeddings,
                            cluster=seurat_object@meta.data$seurat_clusters)
  write.table(Cluster_info,file = file.path(OUT.Seurat.DIR[3],"Cluster_info.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  
  Cluster_num_stat = data.frame(table(seurat_object@meta.data$seurat_clusters),check.names = F)
  colnames(Cluster_num_stat) = c("cluster","number")
  write.table(Cluster_num_stat,file = file.path(OUT.Seurat.DIR[3],"Cluster_num_stat.xls"),quote = F,sep = "\t",col.names = T,row.names = F)
  
  for (i in unique(seurat_object$seurat_clusters)) {
    CATLOG(i)
    
    if(file.exists(file.path(OUT.ClusterAndDiff.DIR[2],paste0("Cluster",i,"_top_gene_feature_byUMAP.pdf")))){
      next
    }
    
    cluster.markers <- FindMarkers(seurat_object, ident.1 = i, min.pct = 0.25,logfc.threshold = 0.25)
    cluster.markers = data.frame(gene=rownames(cluster.markers),cluster.markers,check.names = F)
    cluster.markers = merge(cluster.markers,cluster.averages)
    head(cluster.markers)
    outname = paste0("Cluster",i,".differential_and_annoation.xls")
    write.table(cluster.markers,file = file.path(OUT.ClusterAndDiff.DIR[1],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    top20.tmp = subset(top20,cluster==i)
    cluster.markers.top20 = subset(cluster.markers,gene%in%top20.tmp$gene)
    outname = paste0("Cluster",i,"_top20_gene_diffInf.xls")
    write.table(cluster.markers.top20,file = file.path(OUT.ClusterAndDiff.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    
    Cluster_info.tmp = subset(Cluster_info,cluster==i)
    seurat_object.tmp = seurat_object[top20.tmp$gene,Cluster_info.tmp$barcode]
    seurat_object.tmp = data.frame(barcode=colnames(seurat_object.tmp),t(seurat_object.tmp@assays$RNA@scale.data))
    Cluster_info.tmp = merge(Cluster_info.tmp[,-6],seurat_object.tmp)
    outname = paste0("Cluster",i,"_top20_gene_featureExp.xls")
    write.table(Cluster_info.tmp,file = file.path(OUT.ClusterAndDiff.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
    
    outname = paste0("Cluster",i,"_top_gene_exp.pdf")
    p1 = VlnPlot(seurat_object, features = top20.tmp$gene,ncol = 5,pt.size = 0.1) + NoLegend()
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p1,width =19,height =8)
    
    outname = paste0("Cluster",i,"_top_gene_feature_byTSNE.pdf")
    p2 = FeaturePlot(seurat_object, features = top20.tmp$gene,pt.size = 0.3,reduction = "tsne",ncol = 5)  
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p2,width = 19,height =8)
    
    outname = paste0("Cluster",i,"_top_gene_feature_byUMAP.pdf")
    p3 = FeaturePlot(seurat_object, features = top20.tmp$gene,pt.size = 0.3,reduction = "umap",ncol = 5)  
    ggplot2::ggsave(filename = file.path(OUT.ClusterAndDiff.DIR[2],outname),plot = p3,width = 19,height =8)
  }
  flage(pipeline.log,"3.Seurat.Cluster_and_diff")
}else{
  CATLOG("  Skip Cluster_and_diff")
}

# stat.xls  为提供格式
# 4.Conserved_analysis 目前忽略不做
# 5.Sample_and_diff
# 6.ClusterSample_and_diff

##### ------------------------------- 4.CellAnnotation --------------------------------
CATLOG("4.CellAnnotation")

OUT.celltypeAnnotation.DIR = file.path(OUT.TREE.DIR[4] ,c("2.celltype_diff","3.celltype_marker_gene_plot"))
DIRCREATE(OUT.celltypeAnnotation.DIR)

if(!file.exists(file.path(OUT.TREE.DIR[4],"seurat_object.filted.RDS"))){
  CATLOG("  Run TSNE")
  seurat_object = CreateSeuratObject(counts = Read10X(data.dir = filteredFeatureMatrix),min.cells = 3,min.features = 200,project = "cellranger")
  seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[mM][tT]-") # mm10 mt-  hg19 Mt
  Max.nUMI  = stats::quantile(seurat_object$nCount_RNA,probs=c(0.99))  %>%  as.numeric()
  seurat_object <- subset(seurat_object, percent.mt <= 25&nFeature_RNA>200& nCount_RNA < Max.nUMI)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_object <- FindVariableFeatures(object = seurat_object,selection.method = "vst")
  seurat_object <- ScaleData(seurat_object,features = rownames(seurat_object))     
  seurat_object <- RunPCA(object= seurat_object,pc.genes=VariableFeatures(object = seurat_object),npcs = 50,verbose = F)  
  seurat_object <- JackStraw(seurat_object, num.replicate = 100,dims = 50)
  seurat_object <- ScoreJackStraw(seurat_object, dims = 1:50)
  pcsNumber = which(seurat_object@reductions$pca@jackstraw@overall.p.values[,2] < 0.01)
  seurat_object <- FindNeighbors(seurat_object, dims = pcsNumber)
  seurat_object <- FindClusters(seurat_object, resolution = 1)
  seurat_object <- RunUMAP(seurat_object, dims = pcsNumber)
  seurat_object <- RunTSNE(seurat_object, dims = pcsNumber)
  saveRDS(seurat_object,file = file.path(OUT.TREE.DIR[4] ,"seurat_object.filted.RDS"))
  # saveRDS(seurat_object,file = file.path(OUT.DIR,"seurat_object.filted.R3.RDS"))
}else{
  CATLOG("  Skip TSNE")
  seurat_object <- readRDS(file.path(OUT.TREE.DIR[4] ,"seurat_object.filted.RDS"))  
}

if(!file.exists(file.path(OUT.TREE.DIR[4] ,"seurat_object.filted.celltype.RDS"))){
  ref.rds    = file.path(opts$DBdir,paste0(opts$annoDB,".rds"))
  ref_object = readRDS(ref.rds)
  sc_object = SingleCellExperiment(assays=list(counts=seurat_object@assays$RNA@counts))
  sc_object <- scater::logNormCounts(sc_object)
  common <- intersect(rownames(sc_object), rownames(ref_object))
  ref_object <- ref_object[common,]
  sc_object <- sc_object[common,]
  
  pred <- SingleR(test = sc_object, ref = ref_object, labels = ref_object$label.main)
  anno = data.frame(barcode=colnames(sc_object),
                    labels=pred$labels,
                    first.tuning.scores=pred$tuning.scores$first,
                    second.tuning.scores=pred$tuning.scores$second)
  # write.csv(anno,file = file.path(OUT.Seurat.DIR[4],"anno.csv"),quote = F,row.names = F)
  
  seurat_object <- AddMetaData(seurat_object,metadata =pred$labels,col.name = "celltype" )
 
  CATLOG("111")
  
  p9  = as.ggplot( DimPlot(seurat_object, reduction = "tsne",group.by = "celltype", label = F,pt.size = 0.5))
  print(p9)
  dev.off()
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR,"1.Celltype_byTSNE.pdf"),plot = p9,width = 6,height =5)
  
  CATLOG("222")
  p10 = as.ggplot(DimPlot(seurat_object, reduction = "umap",group.by = "celltype", label = F, pt.size = 0.5))
  dev.off()
  ggplot2::ggsave(filename = file.path(OUT.Seurat.DIR,"1.Celltype_byUMAP.pdf"),plot = p10,width = 6,height =5)
  
  Idents(object = seurat_object ) <- 'celltype'
  saveRDS(seurat_object,file = file.path(OUT.TREE.DIR[4] ,"seurat_object.filted.celltype.RDS"))
}else{
  seurat_object = readRDS(file = file.path(OUT.TREE.DIR[4] ,"seurat_object.filted.celltype.RDS"))
  # subdivide = ifelse(colnames(seurat_object) %>% grepl(pattern = "-1"),"GP1","GP2")  # ------------------------ 修改
  # seurat_object <- Seurat::AddMetaData(seurat_object,metadata =subdivide,col.name = "subdivide" )
}

if(!file.exists(file.path(OUT.TREE.DIR[4] ,"seurat_object.markers.RDS"))){
  seurat_object.markers <- Seurat::FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,min.cells.group = 1)
  saveRDS(seurat_object.markers,file = file.path(OUT.TREE.DIR[4] ,"seurat_object.markers.RDS"))
}else{
  seurat_object.markers = readRDS(file.path(OUT.TREE.DIR[4] ,"seurat_object.markers.RDS"))
}


cluster.averages <- AverageExpression(seurat_object)
cluster.averages <- data.frame(gene=rownames(cluster.averages$RNA),cluster.averages$RNA,check.names = F)

All_cluster_up_marker = data.frame(gene=rownames(seurat_object.markers),seurat_object.markers[,c(1:4,6)],check.names = F)
All_cluster_up_marker = merge(All_cluster_up_marker,cluster.averages)

All_cluster_up_marker = All_cluster_up_marker %>% dplyr::arrange(cluster,desc(cluster))
All_cluster_up_marker = write.table(All_cluster_up_marker,file = file.path(OUT.TREE.DIR[4],"1.All_cluster_up_marker.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

write.table(All_cluster_up_marker,file = file.path(OUT.TREE.DIR,"1.All_celltype_up_marker.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

top20 = seurat_object.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_logFC)
write.table(top20[,-5],file = file.path(OUT.TREE.DIR,"4.Top20_celltype_pos_marker_forheatmap.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

pdf(file = file.path(OUT.TREE.DIR,"5.celltype_top_express_heatmap.pdf"),width = 10,height = 10,onefile = F)
Seurat::DoHeatmap(seurat_object, features = top20$gene,label=FALSE)
dev.off()

celltype_info = data.frame(barcode = colnames(seurat_object),seurat_object@reductions$umap@cell.embeddings,seurat_object@reductions$tsne@cell.embeddings,celltype=seurat_object@meta.data$celltype)
write.table(celltype_info,file = file.path(OUT.TREE.DIR,"celltype_info.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

celltype_num_stat = data.frame(table(seurat_object@meta.data$celltype),check.names = F)
colnames(celltype_num_stat) = c("celltype","number")
write.table(celltype_num_stat,file = file.path(OUT.TREE.DIR,"celltype_num_stat.xls"),quote = F,sep = "\t",col.names = T,row.names = F)

for (i in unique(seurat_object$celltype)) {
  CATLOG(i)
  cluster.markers <- FindMarkers(seurat_object, ident.1 = i, min.pct = 0.25,logfc.threshold = 0.25,min.cells.group = 1)
  cluster.markers = data.frame(gene=rownames(cluster.markers),cluster.markers,check.names = F)
  cluster.markers = merge(cluster.markers,cluster.averages)
  # head(cluster.markers)
  outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),".differential_and_annoation.xls")
  write.table(cluster.markers,file = file.path(OUT.celltypeAnnotation.DIR[1],outname),quote = F,sep = "\t",col.names = T,row.names = F)
  top20.tmp = subset(top20,cluster==i)
  
  cluster.markers.top20 = subset(cluster.markers,gene%in%top20.tmp$gene)
  # head(cluster.markers.top20)
  outname = paste0("celltype.",gsub(pattern = " ",replacement = ".",i),"_top20_gene_diffInf.xls")
  write.table(cluster.markers.top20,file = file.path(OUT.celltypeAnnotation.DIR[2],outname),quote = F,sep = "\t",col.names = T,row.names = F)
  
  celltype_info.tmp = subset(celltype_info,celltype==i)
  # head(celltype_info)
  
  seurat_object.tmp = seurat_object[top20.tmp$gene,celltype_info.tmp$barcode]
  seurat_object.tmp = data.frame(barcode=colnames(seurat_object.tmp),t(seurat_object.tmp@assays$RNA@scale.data))
  # head(seurat_object.tmp)
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

##### ------------------------------- 5.Enrichment --------------------------------
CATLOG("5.Enrichment")
OUT.Enrichment.DIR = file.path(OUT.TREE.DIR[5],c("graphclust","seurat","celltype"))
DIRCREATE(OUT.Enrichment.DIR)
graphclust.diffexp = read.table(analysis.diffexp[1],sep = ",",header = T,stringsAsFactors = F)
clusters = colnames(graphclust.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()

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
        kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
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
      kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
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
  # CATLOG(unique(seurat_object$celltype))
  # unique(seurat_object$celltype)
  
  for (i in unique(seurat_object$celltype)) {
    i=makeNames(i)
    print(i)
    # CATLOG(i)
    tmp.DIR = file.path(OUT.Enrichment.DIR[3],makeNames(i))
    DIRCREATE(tmp.DIR)
    
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/celltype.@{i}.finished.flag"))){
      next
    }
    
    inputname = paste0("celltype.",makeNames(i),".differential_and_annoation.xls")
    cluster.markers = read.table(file.path(OUT.celltypeAnnotation.DIR[1],inputname),header = T,sep = "\t",check.names = F)
    cluster.markers.sig = subset(cluster.markers,p_val_adj<0.05)
    colnames(cluster.markers.sig) %<>% makeNames()
    de.genes.tmp = cluster.markers.sig[,c(1,which(colnames(cluster.markers.sig)==i),3,6)]
    colnames(de.genes.tmp)  = c("symbol","mean_counts","logFC","adj_pval")
    write.table(de.genes.tmp,file = file.path(tmp.DIR,paste0("celltype.",i,".gene.txt")),quote = F,sep = "\t",row.names = F,col.names = T)
    
    
    if(species=="hsa"){
      suppressPackageStartupMessages(library(org.Hs.eg.db))
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
      kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
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
      
      suppressPackageStartupMessages(library(pathview))
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



##### ------------------------------- 6.GSEA --------------------------------
CATLOG("6.GSEA")
OUT.GSEA.DIR = file.path(OUT.TREE.DIR[6],c("graphclust","seurat","celltype"))
DIRCREATE(OUT.GSEA.DIR)

GSEA.sh = "/data/wqihao/soft/GSEA_Linux_4.1.0/gsea-cli.sh" 

# list2gmt(list = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c2.cp.reactome.v7.1.entrez.rds")

gmx.kegg = switch (species,
  "mmu" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Mm.c2.cp.kegg.v7.1.entrez.gmt",
  "hsa" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c2.cp.kegg.v7.1.entrez.gmt"
)
gmx.bp = switch (species,
                   "mmu" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Mm.c5.bp.v7.1.entrez.gmt",
                   "hsa" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c5.bp.v7.1.entrez.gmt"
)
gmx.cc = switch (species,
                 "mmu" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Mm.c5.cc.v7.1.entrez.gmt",
                 "hsa" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c5.cc.v7.1.entrez.gmt"
)
gmx.mf = switch (species,
                 "mmu" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Mm.c5.mf.v7.1.entrez.gmt",
                 "hsa" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c5.mf.v7.1.entrez.gmt"
)
gmx.reactome = switch (species,
                 "mmu" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Mm.c2.cp.reactome.v7.1.entrez.gmt",
                 "hsa" = "/data/wqihao/soft/GSEA_Linux_4.1.0/GSEAdb/MSigDB/v7.1/Hs.c2.cp.reactome.v7.1.entrez.gmt"
)


graphclust.diffexp = read.table(analysis.diffexp[1],sep = ",",header = T,stringsAsFactors = F)
clusters = colnames(graphclust.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()
# graphclust ------
for (cluster in clusters) {
  CATLOG(cluster)
  tmp.DIR = file.path(OUT.GSEA.DIR[1],paste0(cluster,"_VS_clusterRest"))
  DIRCREATE(tmp.DIR)
  if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
    next
  }
  
  p.cols = grep(pattern = paste0(cluster,".Adjusted.p.value"),colnames(graphclust.diffexp))
  exp.cols = grep(pattern = paste0(cluster,".Mean.Counts"),colnames(graphclust.diffexp))
  fc.cols = grep(pattern = paste0(cluster,".Log2.fold.change"),colnames(graphclust.diffexp))
  
  counts.ID =  colnames(graphclust.diffexp) %>% grep(pattern = ".Mean.Counts",value = T) 
  
  de.genes = graphclust.diffexp[abs(graphclust.diffexp[,fc.cols])>1&graphclust.diffexp[,p.cols]<0.05,]
  xls = data.frame(gene =graphclust.diffexp$Feature.Name,V2="Na",V3="Na",V4="Na",V5="Na",V6="Na", graphclust.diffexp[,counts.ID])
  colnames(xls) = colnames(xls) %>%  gsub(pattern = ".Mean.Counts",replacement = "") 
  xls.tmp = tempfile()
  write.table(xls,file = xls.tmp,quote = F,sep = "\t",row.names = F,col.names = T)
  # out = file.path(tmp.DIR,cluster)
  
  xls2gct(xls = xls.tmp,species =species ,out = tmp.DIR,celltype =cluster )
  
  res = file.path(tmp.DIR,paste0(cluster,".gct"))
  cls = file.path(tmp.DIR,paste0(cluster,".cls"))
  
  GSEA.list = switch (species,
    "mmu" = c("KEGG","GO_BP","GO_CC","GO_MF"),
    "hsa" = c("KEGG","GO_BP","GO_CC","GO_MF","REACTOME")
  )
  i=GSEA.list[1]
  for (i in GSEA.list) {
    rpt_label = switch (i,
    "KEGG"    = "c2.KEGG",
    "GO_BP"   = "c5.GO_BP",
    "GO_CC"   = "c5.GO_CC",
    "GO_MF"   = "c5.GO_MF",
    "REACTOME" = "c2.REACTOME"
    )
    
    gmx = switch (i,
    "KEGG"    = gmx.kegg,
    "GO_BP"   = gmx.bp,
    "GO_CC"   = gmx.cc,
    "GO_MF"   = gmx.mf,
    "REACTOME" = gmx.reactome
    )
    
    cmd = GetoptLong::qq("@{GSEA.sh} GSEA -res @{res} -gmx @{gmx} -cls @{cls} -collapse false -mode Max_probe -norm meandiv -nperm 1000  -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label @{rpt_label}  -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false  -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 0 -zip_report false -out @{tmp.DIR}")
    system(cmd)
    system(GetoptLong::qq("cd @{tmp.DIR};rename @{rpt_label}.Gsea.* @{rpt_label}.Gsea *"))
    
  }
  system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
  
}

# SEURAT ---------
Cluster_diff.file = list.files(OUT.ClusterAndDiff.DIR[1],pattern = ".xls",full.names = T)
clusters = Cluster_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
cluster = clusters[1]
for (cluster in clusters) {
  # filter top 20 genes
  CATLOG(cluster)
  tmp.DIR = file.path(OUT.GSEA.DIR[2],paste0(cluster,"_VS_clusterRest"))
  DIRCREATE(tmp.DIR)
  if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
    next
  }
  
  xls.tmp = grep(pattern = paste0(cluster,".differential_and_annoation.xls"),Cluster_diff.file,value = T)
  
  xls2gct(xls = xls.tmp,species =species ,out = tmp.DIR,celltype =cluster )
  
  res = file.path(tmp.DIR,paste0(cluster,".gct"))
  cls = file.path(tmp.DIR,paste0(cluster,".cls"))
  
  GSEA.list = switch (species,
                      "mmu" = c("KEGG","GO_BP","GO_CC","GO_MF"),
                      "hsa" = c("KEGG","GO_BP","GO_CC","GO_MF","REACTOME")
  )
  i=GSEA.list[1]
  for (i in GSEA.list) {
    rpt_label = switch (i,
                        "KEGG"    = "c2.KEGG",
                        "GO_BP"   = "c5.GO_BP",
                        "GO_CC"   = "c5.GO_CC",
                        "GO_MF"   = "c5.GO_MF",
                        "REACTOME" = "c2.REACTOME"
    )
    
    gmx = switch (i,
                  "KEGG"    = gmx.kegg,
                  "GO_BP"   = gmx.bp,
                  "GO_CC"   = gmx.cc,
                  "GO_MF"   = gmx.mf,
                  "REACTOME" = gmx.reactome
    )
    
    cmd = GetoptLong::qq("@{GSEA.sh} GSEA -res @{res} -gmx @{gmx} -cls @{cls} -collapse false -mode Max_probe -norm meandiv -nperm 1000  -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label @{rpt_label}  -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false  -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 0 -zip_report false -out @{tmp.DIR}")
    system(cmd)
    system(GetoptLong::qq("cd @{tmp.DIR};rename @{rpt_label}.Gsea.* @{rpt_label}.Gsea *"))
  }
  system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
}


# celltype --------
celltype_diff.file = list.files(OUT.celltypeAnnotation.DIR[1],pattern = ".xls",full.names = T)
celltypes = celltype_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
# celltype = celltypes[1]
for (celltype in celltypes) {
  # filter top 20 genes
  CATLOG(celltype)
  tmp.DIR = file.path(OUT.GSEA.DIR[3],paste0(celltype,"_VS_celltypeRest"))
  DIRCREATE(tmp.DIR)
  if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{celltype}.finished.flag"))){
    next
  }
  
  xls.tmp = grep(pattern = paste0(celltype,".differential_and_annoation.xls"),celltype_diff.file,value = T)
  
  xls2gct(xls = xls.tmp,species =species ,out = tmp.DIR,celltype =celltype )
  
  res = file.path(tmp.DIR,paste0(celltype,".gct"))
  cls = file.path(tmp.DIR,paste0(celltype,".cls"))
  
  GSEA.list = switch (species,
                      "mmu" = c("KEGG","GO_BP","GO_CC","GO_MF"),
                      "hsa" = c("KEGG","GO_BP","GO_CC","GO_MF","REACTOME")
  )
  i=GSEA.list[1]
  for (i in GSEA.list) {
    rpt_label = switch (i,
                        "KEGG"    = "c2.KEGG",
                        "GO_BP"   = "c5.GO_BP",
                        "GO_CC"   = "c5.GO_CC",
                        "GO_MF"   = "c5.GO_MF",
                        "REACTOME" = "c2.REACTOME"
    )
    
    gmx = switch (i,
                  "KEGG"    = gmx.kegg,
                  "GO_BP"   = gmx.bp,
                  "GO_CC"   = gmx.cc,
                  "GO_MF"   = gmx.mf,
                  "REACTOME" = gmx.reactome
    )
    
    cmd = GetoptLong::qq("@{GSEA.sh} GSEA -res @{res} -gmx @{gmx} -cls @{cls} -collapse false -mode Max_probe -norm meandiv -nperm 1000  -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label @{rpt_label}  -metric Signal2Noise -sort real -order descending -include_only_symbols true -make_sets true -median false  -num 100 -plot_top_x 20 -save_rnd_lists false -set_max 500 -set_min 0 -zip_report false -out @{tmp.DIR}")
    system(cmd)
    system(GetoptLong::qq("cd @{tmp.DIR};rename @{rpt_label}.Gsea.* @{rpt_label}.Gsea *"))
  }
  system(GetoptLong::qq("touch @{tmp.DIR}/@{celltype}.finished.flag"))
}

##### ------------------------------- 7.PPI --------------------------------
CATLOG("7.PPI")
library(STRINGdb)
species = "mmu"
proteins.DB = switch (species,
  "hsa" = "/data/wqihao/database/PPI/9606__proteins.tsv.gz",
  "mmu" = "/data/wqihao/database/PPI/10090__proteins.tsv.gz"
)
protein_links.DB = switch (species,
  "hsa" = "/data/wqihao/database/PPI/9606__protein_links.tsv.gz",
  "mmu" = "/data/wqihao/database/PPI/10090__protein_links.tsv.gz"
)

OUT.PPI.DIR = file.path(OUT.TREE.DIR[7],c("seurat","celltype","graphclust"))
DIRCREATE(OUT.PPI.DIR)

# graphclust

graphclust.diffexp = read.table(analysis.diffexp[1],sep = ",",header = T,stringsAsFactors = F)
clusters = colnames(graphclust.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()
for (cluster in clusters) {
  CATLOG(cluster)
  if(file.exists(file.path(OUT.PPI.DIR[3],paste0(cluster,".propro.node.xls")))|file.exists(file.path(OUT.PPI.DIR[3],paste0(cluster,".err.txt")))){
    next
  }
  
  p.cols = grep(pattern = paste0(cluster,".Adjusted.p.value"),colnames(graphclust.diffexp))
  exp.cols = grep(pattern = paste0(cluster,".Mean.Counts"),colnames(graphclust.diffexp))
  fc.cols = grep(pattern = paste0(cluster,".Log2.fold.change"),colnames(graphclust.diffexp))
 
  diff_exp = graphclust.diffexp[abs(graphclust.diffexp[,fc.cols])>1&graphclust.diffexp[,p.cols]<0.05,c(p.cols,fc.cols,2)]
  colnames(diff_exp) =  c("p_val_adj","avg_logFC","gene")
  string_db <- STRINGdb::STRINGdb$new( version="10", 
                                       species=ifelse(species=="hsa",9606,10090),
                                       score_threshold=0, 
                                       input_directory= "/data/wqihao/database/PPI/") # default the score threshold is set to 400 
  if(nrow(diff_exp)>=3){
    diff_exp_mapped <- string_db$map( diff_exp, "gene", removeUnmappedRows = TRUE ) %>% invisible()
    mapped = structure(diff_exp_mapped$gene,names=diff_exp_mapped$STRING_id)
    proteins.tsv.gz = data.table::fread(proteins.DB)
    protein_links.tsv.gz = data.table::fread(protein_links.DB)
    protein_links.tsv.filter = subset(protein_links.tsv.gz,protein1%in%diff_exp_mapped$STRING_id&protein2%in%diff_exp_mapped$STRING_id)
    protein_links.tsv.filter$protein1 = mapped[protein_links.tsv.filter$protein1]
    protein_links.tsv.filter$protein2 = mapped[protein_links.tsv.filter$protein2]
    n = nrow(protein_links.tsv.filter)
    if(n>1500){
      probs = 1500/n
    }else{
      probs = 0
    }
    FC = quantile(protein_links.tsv.filter$combined_score,probs = probs)
    protein_links.tsv.filter.FC = subset(protein_links.tsv.filter,combined_score >= FC)
    plot.ppi(PPI = protein_links.tsv.filter.FC,outdir = OUT.PPI.DIR[3],cluster = cluster) 
  }else{
    logFile=file.path(OUT.PPI.DIR[3],paste0(cluster,".err.txt"))
    write.table(paste0("[error: PPI] No PPI target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
  }
  
}

# seurat
Cluster_diff.file = list.files(OUT.ClusterAndDiff.DIR[1],pattern = ".xls",full.names = T)
clusters = Cluster_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
for (cluster in clusters) {
  # filter top 20 genes
  CATLOG(cluster)
  if(file.exists(file.path(OUT.PPI.DIR[1],paste0(cluster,".propro.node.xls")))|file.exists(file.path(OUT.PPI.DIR[1],paste0(cluster,".err.txt")))){
    next
  }
  string_db <- STRINGdb::STRINGdb$new( version="10", 
                                       species=ifelse(species=="hsa",9606,10090),
                                       score_threshold=0, 
                                       input_directory= "/data/wqihao/database/PPI/") # default the score threshold is set to 400 

  Cluster_diff.tmp = data.table::fread(grep(pattern = paste0(cluster,".differential_and_annoation.xls"),Cluster_diff.file,value = T))
  diff_exp = subset(Cluster_diff.tmp,p_val_adj<0.05,select = c("p_val_adj","avg_logFC","gene"))
  if(nrow(diff_exp)!=0){
    diff_exp_mapped <- string_db$map( diff_exp, "gene", removeUnmappedRows = TRUE ) %>% invisible()
    mapped = structure(diff_exp_mapped$gene,names=diff_exp_mapped$STRING_id)
    proteins.tsv.gz = data.table::fread(proteins.DB)
    protein_links.tsv.gz = data.table::fread(protein_links.DB)
    protein_links.tsv.filter = subset(protein_links.tsv.gz,protein1%in%diff_exp_mapped$STRING_id&protein2%in%diff_exp_mapped$STRING_id)
    protein_links.tsv.filter$protein1 = mapped[protein_links.tsv.filter$protein1]
    protein_links.tsv.filter$protein2 = mapped[protein_links.tsv.filter$protein2]
    n = nrow(protein_links.tsv.filter)
    if(n>1500){
      probs = 1500/n
    }else{
      probs = 0
    }
    FC = quantile(protein_links.tsv.filter$combined_score,probs = probs)
    protein_links.tsv.filter.FC = subset(protein_links.tsv.filter,combined_score >= FC)
    plot.ppi(PPI = protein_links.tsv.filter.FC,outdir = OUT.PPI.DIR[1],cluster = cluster) 
  }else{
    logFile=file.path(OUT.PPI.DIR[1],paste0(cluster,".err.txt"))
    write.table(paste0("[error: PPI] No PPI target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
  }
  # string_db=NULL
}

# celltype
celltype_diff.file = list.files(OUT.celltypeAnnotation.DIR[1],pattern = ".xls",full.names = T)
celltypes = celltype_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
for (celltype in celltypes) {
  # filter top 20 genes
  CATLOG(celltype)
  if(file.exists(file.path(OUT.PPI.DIR[2],paste0(celltype,".propro.node.xls")))|file.exists(file.path(OUT.PPI.DIR[2],paste0(celltype,".err.txt")))){
    next
  }
  string_db <- STRINGdb::STRINGdb$new( version="10", 
                                       species=ifelse(species=="hsa",9606,10090),
                                       score_threshold=0, 
                                       input_directory= "/data/wqihao/database/PPI/") # default the score threshold is set to 400 

  celltype_diff.tmp = data.table::fread(grep(pattern = paste0(celltype,".differential_and_annoation.xls"),celltype_diff.file,value = T))
  diff_exp = subset(celltype_diff.tmp,p_val_adj<0.05,select = c("p_val_adj","avg_logFC","gene"))
  if(nrow(diff_exp)!=0){
    diff_exp_mapped <- string_db$map( diff_exp, "gene", removeUnmappedRows = TRUE ) %>% invisible()
    mapped = structure(diff_exp_mapped$gene,names=diff_exp_mapped$STRING_id)
    proteins.tsv.gz = data.table::fread(proteins.DB)
    protein_links.tsv.gz = data.table::fread(protein_links.DB)
    protein_links.tsv.filter = subset(protein_links.tsv.gz,protein1%in%diff_exp_mapped$STRING_id&protein2%in%diff_exp_mapped$STRING_id)
    protein_links.tsv.filter$protein1 = mapped[protein_links.tsv.filter$protein1]
    protein_links.tsv.filter$protein2 = mapped[protein_links.tsv.filter$protein2]
    n = nrow(protein_links.tsv.filter)
    if(n>1500){
      probs = 1500/n
    }else{
      probs = 0
    }
    FC = quantile(protein_links.tsv.filter$combined_score,probs = probs)
    protein_links.tsv.filter.FC = subset(protein_links.tsv.filter,combined_score >= FC)
    plot.ppi(PPI = protein_links.tsv.filter.FC,outdir = OUT.PPI.DIR[2],cluster = celltype) 
  }else{
    logFile=file.path(OUT.PPI.DIR[2],paste0(celltype,".err.txt"))
    write.table(paste0("[error: PPI] No PPI target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
  }
  # string_db=NULL
}


##### ------------------------------- 8.TFBS --------------------------------
CATLOG("8.TFBS")

#1 Cluster (graphclust seurat)
OUT.TFBS.DIR = file.path(OUT.TREE.DIR[8],c("seurat","celltype","graphclust"))
DIRCREATE(OUT.TFBS.DIR)

# seurat
Cluster_diff.file = list.files(OUT.ClusterAndDiff.DIR[1],pattern = ".xls",full.names = T)
clusters = Cluster_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
# cluster = clusters[4]
for (cluster in clusters) {
  # filter top 20 genes
  CATLOG(cluster)
  tmp.DIR = file.path(OUT.TFBS.DIR[1],cluster)
  DIRCREATE(tmp.DIR)
  
  if(file.exists(file.path(tmp.DIR,"SCENIC.tmp/output","Step2_regulonTargetsInfo.tsv"))){
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
      next  
    }else{
      if(!file.exists(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))){
        system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
        next
      }
      tab = fread(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))
      plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=cluster)
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
    }
  }else{
    
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
      next  
    }
    Cluster_diff.tmp = fread(grep(pattern = paste0(cluster,".differential_and_annoation.xls"),Cluster_diff.file,value = T))
    call.SCENIC(tmp.DIR = tmp.DIR,org = "mgi",Cluster_diff.tmp = Cluster_diff.tmp,cluster = cluster,seurat_object = seurat_object)
    if(!file.exists(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))){
      system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
      next
    }
    tab = fread(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))
    plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=cluster)
    # flag for all finished -----------------
    system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
  }
}


#2 celltype 
celltype_diff.file = list.files(OUT.celltypeAnnotation.DIR[1],pattern = ".xls",full.names = T)
celltypes = celltype_diff.file %>% basename() %>% gsub(pattern = ".differential_and_annoation.xls",replacement = "") %>% unique()
# celltype=celltypes[6]
for (celltype in celltypes) {
  # filter top 20 genes
  CATLOG(celltype)
  tmp.DIR = file.path(OUT.TFBS.DIR[2],celltype)
  DIRCREATE(tmp.DIR)
  
  if(file.exists(file.path(tmp.DIR,"SCENIC.tmp/output","Step2_regulonTargetsInfo.tsv"))){
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{celltype}.finished.flag"))){
      next  
    }else{
      if(!file.exists(file.path(tmp.DIR,paste0(celltype,".predicted_TF.xls")))){
        system(GetoptLong::qq("touch @{tmp.DIR}/@{celltype}.finished.flag"))
        next
      }
      tab = fread(file.path(tmp.DIR,paste0(celltype,".predicted_TF.xls")))
      plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=celltype)
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/@{celltype}.finished.flag"))
    }
  }else{
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{celltype}.finished.flag"))){
      next  
    }
    
    celltype_diff.tmp = data.table::fread(grep(pattern = paste0(celltype,".differential_and_annoation.xls"),celltype_diff.file,value = T))
    call.SCENIC(tmp.DIR = tmp.DIR,org = "mgi",Cluster_diff.tmp = celltype_diff.tmp,cluster = celltype,seurat_object = seurat_object)
    
    if(!file.exists(file.path(tmp.DIR,paste0(celltype,".predicted_TF.xls")))){
      system(GetoptLong::qq("touch @{tmp.DIR}/@{celltype}.finished.flag"))
      next
    }
    tab = fread(file.path(tmp.DIR,paste0(celltype,".predicted_TF.xls")))
    plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=celltype)
    # flag for all finished -----------------
    system(GetoptLong::qq("touch @{tmp.DIR}/@{celltype}.finished.flag"))
  }
  
}

# graphclust 
graphclust.diffexp = read.table(analysis.diffexp[1],sep = ",",header = T,stringsAsFactors = F)
clusters = colnames(graphclust.diffexp)[-c(1,2)] %>% grep(pattern = ".Mean.Counts",value = T) %>% gsub(pattern = ".Mean.Counts",replacement = "") %>% unique()
for (cluster in clusters) {
  CATLOG(cluster)
  tmp.DIR = file.path(OUT.TFBS.DIR[3],cluster)
  DIRCREATE(tmp.DIR)
  
  if(file.exists(file.path(tmp.DIR,"SCENIC.tmp/output","Step2_regulonTargetsInfo.tsv"))){
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
      next  
    }else{
      if(!file.exists(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))){
        system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
        next
      }
      tab = fread(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))
      plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=cluster)
      # flag for all finished -----------------
      system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
    }
    
    
  }else{
    if(file.exists(GetoptLong::qq("@{tmp.DIR}/@{cluster}.finished.flag"))){
      next  
    }
    
    p.cols = grep(pattern = paste0(cluster,".Adjusted.p.value"),colnames(graphclust.diffexp))
    exp.cols = grep(pattern = paste0(cluster,".Mean.Counts"),colnames(graphclust.diffexp))
    fc.cols = grep(pattern = paste0(cluster,".Log2.fold.change"),colnames(graphclust.diffexp))
    
    Cluster_diff.tmp = graphclust.diffexp[abs(graphclust.diffexp[,fc.cols])>1&graphclust.diffexp[,p.cols]<0.05,c(p.cols,fc.cols,2)]
    colnames(Cluster_diff.tmp) =  c("p_val_adj","avg_logFC","gene")

    call.SCENIC(tmp.DIR = tmp.DIR,org = "mgi",Cluster_diff.tmp = Cluster_diff.tmp,cluster = cluster,seurat_object = seurat_object)
    
    if(!file.exists(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))){
      system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
      next
    }
    
    tab = fread(file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")))
    plot.TFBS(TFBS=tab,outdir=tmp.DIR,cluster=cluster)
    # flag for all finished -----------------
    system(GetoptLong::qq("touch @{tmp.DIR}/@{cluster}.finished.flag"))
    
  }
}


# 9.Convert --------
CATLOG("9.Convert pdf2png")
all.pdf = list.files(OUT.DIR,pattern = ".pdf$",recursive = T,full.names = T)
all.pdf.rm = grep(pattern = "Rplots.pdf",all.pdf,value = T)
file.remove(all.pdf.rm) %>% invisible()
all.pdf = grep(pattern = "Rplots.pdf",all.pdf,value = T,invert = T)
all.png = gsub(pattern = "pdf",replacement = "png",all.pdf)
cmd = sprintf("/usr/bin/convert -density 300 %s %s",all.pdf,all.png)
Cmd(cmd = cmd,cl =16 ,out = all.png,step = "Convert pdf 2 png",f = system)
flage(pipeline.log,"9.Convert.pdf2png")


# 10.html report
srcDIR = file.path(OUT.RESTULE.SAMPLEN.DIR,"src")
DIRCREATE(dirs = srcDIR)

srcfiles = list.files("/data/wqihao/Project.sc/workflow/src/",full.names = T,recursive = T)
COPYFILE(from.name = srcfiles,to.dir = srcDIR,step = "Copy src files")

if(any(grepl(pattern = "vdj",fastqFiles))){
  if(check.flage(pipeline.log,"10.html")){
    CATLOG("10.html Rmd2html")
    file.copy(from = "/data/wqihao/Project.sc/workflow/singlecell_vdj_report.Rmd",to = file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_vdj_report.Rmd"))
    rmarkdown::render(file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_vdj_report.Rmd"))
    file.remove(file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_vdj_report.Rmd"))
    file.remove(file.path(OUT.TREE.DIR[1],"per_base_quality.png"))
    file.remove(file.path(OUT.TREE.DIR[1],"per_base_sequence_content.png"))
    flage(pipeline.log,"10.html")
  }else{
    CATLOG("10.html Rmd2html")
  }
}else{
  if(check.flage(pipeline.log,"10.html")){
    CATLOG("10.html Rmd2html")
    file.copy(from = "/data/wqihao/Project.sc/workflow/singlecell_report.Rmd",to = file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_report.Rmd"))
    rmarkdown::render(file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_report.Rmd"))
    file.remove(file.path(OUT.RESTULE.SAMPLEN.DIR,"singlecell_report.Rmd"))
    file.remove(file.path(OUT.TREE.DIR[1],"per_base_sequence_content.png"))
    flage(pipeline.log,"10.html")
  }else{
    CATLOG("10.html Rmd2html")
  }
  
}


##### ------------------------------- 9.monocle(待完善)
##### ------------------------------- 10.细胞间相互作用(待完善)


