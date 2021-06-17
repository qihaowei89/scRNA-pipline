#!/usr/bin/env Rscript
# library 
options(warn = -1)
suppressPackageStartupMessages({
  if(interactive()){
    .libPaths(new = "~/miniconda3/lib/R/library")
  }
  packages <- c("doMC","iterators","clusterProfiler","enrichplot","ComplexHeatmap",
                "GetoptLong","data.table","ReactomePA","patchwork","ggplotify",
                "optparse","tximport","magrittr","pathview","biomaRt",
                "scRNAseq","SingleR","Seurat","scater","dplyr",
                "dplyr","hdf5r","Rmisc","plyr","sva","DOSE","crayon")
  for (p in packages) {
    if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
      if(!"BiocManager"%in%rownames(installed.packages())) install.packages("BiocManager")
      # p = "dplyr"
      BiocManager::install(p,BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor",update = F)
      suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
    }
  }
})

# function
CATLOG <- function(text){
  suppressPackageStartupMessages({
    require(GetoptLong)
    library(crayon)
  })
  error <- red $ bold
  warn <- magenta $ underline
  note <- cyan $ bold
  # cat(error("Error: subscript out of bounds!\n"))
  # cat(warn("Warning: shorter argument was recycled.\n"))
  # cat(note("Note: no such directory.\n"))
  qq.options("cat_prefix" = function(x) format(Sys.time(), paste0("[",note("INFO ") ,"%Y-%m-%d %H:%M:%S]: ")))
  # qqcat("@{text}\n")
  qqcat(blue("@{text}\n"))
  qq.options(RESET = TRUE)
}

CATERR <- function(text){
  suppressPackageStartupMessages({
    require(GetoptLong)
    library(crayon)
  })
  error <- red $ bold
  warn <- magenta $ underline
  note <- cyan $ bold
  # cat(error("Error: subscript out of bounds!\n"))
  # cat(warn("Warning: shorter argument was recycled.\n"))
  # cat(note("Note: no such directory.\n"))
  qq.options("cat_prefix" = function(x) format(Sys.time(), paste0("[",error("Error ") ,"%Y-%m-%d %H:%M:%S]: ")))
  # qqcat("@{text}\n")
  qqcat(warn("@{text}\n"))
  qq.options(RESET = TRUE)
}

save.plots <- function(p,file,width,height){
  require(ggplot2)
  dir.create(dirname(file),recursive = T)
  file = paste0(file,c(".pdf",".jpg")) 
  ggsave(file[1],plot = print(p),width = width,height = height,onefile=F)
  ggsave(file[2],plot = print(p),width = width,height = height,units = "in",dpi = "print")
}
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))
  
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))
  head(clono)
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3)]
  rownames(tcr) <- tcr[,1]
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
add_clonotype2 <- function(tcr_prefix, seurat_obj, type="b"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep="/"))
  # Remove the -1 at the end of each barcode.
  # Subsets so only the first line of each barcode is kept,
  # as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcr <- tcr[,c("barcode", "raw_clonotype_id","chain","v_gene","d_gene","j_gene","c_gene")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep="/"))
  head(clono)
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"
  head(tcr)
  
  # Reorder so barcodes are first column and set them as rownames.
  tcr <- tcr[, c(2,1,3:8)]
  
  rownames(tcr) <- tcr[,1] %>% gsub(pattern = "-1",replacement ="" )
  tcr[,1] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

Cmd <- function(cmd,cl=8,step,f=System,out) 
{ 
  System <- function(cmd)
    {
       message(cmd)
       system(cmd)
    }
  library(doParallel)
  # check if out file exists
  cmd_need_run <- cmd[!file.exists(out)]
  out_need_get <- out[!file.exists(out)]
  doMC::registerDoMC(cl)
  if (length(cmd_need_run)!=0){
    CATLOG(sprintf("Start %s total(%s), run(%s), skip(%s)",step,length(cmd),length(cmd_need_run),length(cmd)-length(cmd_need_run)))
    stat = foreach(a=iter(cmd_need_run),.combine = c) %dopar% f(a)
    err_cmd <- which(stat!=0)
    if(length(err_cmd!=0)) {
      system(sprintf("rm -rf %s",out_need_get[err_cmd]))
      CATERR(step)
      break
    }
    
    CATLOG(sprintf("End %s total(%s), run(%s), skip(%s)",step,length(cmd), length(cmd_need_run)-length(err_cmd),length(err_cmd)))
  }else{
    CATLOG(sprintf("  Skip %s",step))
  }
}
PLOT.TSNE  <- function(tsneX,tsneY,Clusters,type,title,size=0.2,xlab="TSNE.1",ylab="TSNE.2"){
  if(type=="exp") {
    summarys =  round(summary(Clusters),digits = 3)
    col <- colorRampPalette(c("gray","red"))(100)
    tsne_plot <- data.frame(x = tsneX, y = tsneY, col = Clusters)
    p1 = ggplot2::ggplot(tsne_plot) + 
      ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=col),size=size) + 
      ggplot2::theme_test()+
      ggplot2::theme(legend.title = ggplot2::element_blank()) + 
      ggplot2::ggtitle(title)+
      ggplot2::scale_colour_gradientn(colours = col) +
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) + 
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+6*max(tsneY)/10,
                        label=paste0(names(summarys)[1],":",summarys[1])) +
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+5*max(tsneY)/10,
                        label=paste0(names(summarys)[2],":",summarys[2])) + 
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+4*max(tsneY)/10,
                        label=paste0(names(summarys)[3],":",summarys[3])) + 
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+3*max(tsneY)/10,
                        label=paste0(names(summarys)[4],":",summarys[4])) + 
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+2*max(tsneY)/10,
                        label=paste0(names(summarys)[5],":",summarys[5])) +
      ggplot2::annotate("text",x = min(tsneX)+max(tsneX)/10,y = min(tsneY)+max(tsneY)/10,
                        label=paste0(names(summarys)[6],":",summarys[6]))
    
    
  }
  if(type=="Cluster") {
    # col <- colorRampPalette(c("gray","red"))(100)
    tsne_plot <- data.frame(x = tsneX, y = tsneY, col = as.character(Clusters))
    p1 = ggplot2::ggplot(tsne_plot) + 
      ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=col),size=size) + 
      ggplot2::theme_test()+
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(title)+
      # ggplot2::scale_colour_gradientn(colours = col) +
      ggplot2::xlab(xlab)  + 
      ggplot2::ylab(ylab)
  }
  return(p1)
}
PLOT.TSNE2 <- function(tsneX,tsneY,Clusters,type,title,size=0.5,labs,xlab="TSNE.1",ylab="TSNE.2"){
  if(type=="exp") {
    summarys =  round(summary(Clusters),digits = 3)
    col <- colorRampPalette(c("blue","green","orange","red"))(round(max(Clusters),0))
    
    tsne_plot <- data.frame(x = tsneX, y = tsneY, col = Clusters)
    
    p1 = ggplot2::ggplot(tsne_plot) + 
      ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=col),size=size) + 
      ggplot2::theme_test()+
      # ggplot2::theme(legend.title = ggplot2::element_blank()) + 
      ggplot2::ggtitle(title)+
      ggplot2::scale_colour_gradientn(colours = col) +
      ggplot2::labs(col=labs)+
      ggplot2::xlab(xlab) + 
      ggplot2::ylab(ylab) 

  }
  return(p1)
}
PLOT.TSNE3 <- function(tsneX,tsneY,Clusters,type,title,size=0.2,xlab="TSNE.1",ylab="TSNE.2"){
  # col <- colorRampPalette(c("gray","red"))(100)
  tsne_plot <- data.frame(x = tsneX, y = tsneY,type = type,col = as.character(Clusters))
  p1 = ggplot2::ggplot(tsne_plot) + 
    ggplot2::geom_point(ggplot2::aes(x=x, y=y, color=col),size=size) + 
    ggplot2::theme_test()+
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::ggtitle(title)+
    # ggplot2::facet_wrap() +
    ggplot2::xlab(xlab)  + 
    ggplot2::ylab(ylab)
  return(p1)
}
volcanoPlot <- function(logFC,pValue,FC,main,size){
  options(na.rm = T)
  dat.mRNA <- data.frame(logFC,pValue)
  dat.mRNA$significant <- "normal"
  dat.mRNA$significant[dat.mRNA$logFC >  FC & dat.mRNA$pValue<0.05] <- "up"
  dat.mRNA$significant[dat.mRNA$logFC < -FC & dat.mRNA$pValue<0.05] <- "down"
  dat.mRNA$significant <- factor(dat.mRNA$significant,levels = c("up","down","normal"))
  dat.mRNA$significant %>% table()
  dat.mRNA$NeglogP <- -log10(dat.mRNA$pValue+0.0000000001)
  significant <- dat.mRNA$significant
  p <- ggplot2::qplot(x=dat.mRNA$logFC,y=dat.mRNA$NeglogP,
                      xlab = expression(paste(log[2],"(FC)")),
                      ylab = expression(paste(-log[10],"(p-value)")),
                      size=I(size),colour =significant,
                      xlim = c(-max(abs(dat.mRNA$logFC),na.rm = T),max(abs(dat.mRNA$logFC),na.rm = T)))
  p <- p+ ggplot2::scale_color_manual(values = c("up"="red","down"="#0DAB00","normal"="grey"))
  xline <- c(-FC,FC)
  yline <- -log10(0.05)
  p <- p + ggplot2::geom_vline(xintercept=xline,lty=2,size=I(0.5),colour="grey11")
  p <- p + ggplot2::geom_hline(yintercept=yline,lty=2,size=I(0.5),colour="grey11")
  p <- p + ggplot2::theme_bw() + 
    ggplot2::labs(x=expression(paste(log[2],"(FC)")),
                  y=expression(paste(-log[10],"(p-value)")),
                  title=main)+
    ggplot2::theme(axis.text.y = ggplot2::element_text(angle=90, hjust=0.5,size = 15),
                   axis.text = ggplot2::element_text(colour = "black",size = 15),
                   title =  ggplot2::element_text(colour = "black",size = 13,face = "bold"),
                   panel.background = ggplot2::element_rect(colour = "black",size = 1,fill = "white"),
                   panel.grid=ggplot2::element_blank(),
                   plot.margin=unit(rep(1.5,4),'lines'),
                   legend.position="none",
                   legend.title = ggplot2::element_text(colour = "black",size = 14),
                   legend.text = ggplot2::element_text(colour = "black",size = 12),
                   axis.title.x=ggplot2::element_text(colour = "black",size = 15,vjust = -1),
                   axis.title.y=ggplot2::element_text(colour = "black",size = 15,vjust = 1))
  p <- p + ggplot2::annotate("text", 
                             x = c(-max(abs(dat.mRNA$logFC),na.rm = T)+max(abs(dat.mRNA$logFC),na.rm = T)/3,
                                   -max(abs(dat.mRNA$logFC),na.rm = T)+max(abs(dat.mRNA$logFC),na.rm = T)/3), 
                             y =  c(max(dat.mRNA$NeglogP,na.rm = T)-max(dat.mRNA$NeglogP,na.rm = T)/16, 
                                    max(dat.mRNA$NeglogP,na.rm = T)-max(dat.mRNA$NeglogP,na.rm = T)*2/16 ), 
                             size = 4,
                             label = c(paste0("bold(down): ","bold(",sum(dat.mRNA$significant=="down"),")"), 
                                       paste0("bold(up): ","bold(",sum(dat.mRNA$significant=="up"),")")),parse = TRUE,
                             colour = c("#0DAB00","red"))
  p
}
DIRCREATE = function(dirs){
  invisible(sapply(dirs, dir.create,recursive = T))
}
COPYFILE = function(from.name,to.name="",to.dir="",step=""){
  
  if(to.name[1]==""&to.dir[1]!=""){
    if(!dir.exists(to.dir)) dir.create(to.dir,recursive = T)
    to.name = file.path(to.dir,basename(from.name))
    cmd = GetoptLong::qq("cp -r @{from.name} @{to.name}\n") %>% strsplit("\n") %>% '[['(1) 
    Cmd(cmd = cmd,cl = 10,step = step,out = to.name)
  }
  
  if(to.dir[1]==""&to.name[1]!=""){
    to.dir = base::dirname(to.name)
    if(!dir.exists(to.dir)) dir.create(to.dir,recursive = T)
    cmd = GetoptLong::qq("cp -r @{from.name} @{to.name}\n") %>% strsplit("\n") %>% '[['(1) 
    Cmd(cmd = cmd,cl = 10,step = step,out = to.name)
  }
  
}

MOVEFILE = function(from.name,to.name="",to.dir="",step=""){
  
  if(to.name[1]==""&to.dir[1]!=""){
    if(!dir.exists(to.dir)) dir.create(to.dir,recursive = T)
    to.name = file.path(to.dir,basename(from.name))
    cmd = GetoptLong::qq("mv @{from.name} @{to.name}\n") %>% strsplit("\n") %>% '[['(1) 
    Cmd(cmd = cmd,cl = 10,step = step,out = to.name)
  }
  
  if(to.dir[1]==""&to.name[1]!=""){
    to.dir = base::dirname(to.name)
    if(!dir.exists(to.dir)) dir.create(to.dir,recursive = T)
    cmd = GetoptLong::qq("mv @{from.name} @{to.name}\n") %>% strsplit("\n") %>% '[['(1) 
    Cmd(cmd = cmd,cl = 10,step = step,out = to.name)
  }
  
}



GETFILE = function(input.dir,pattern){
  list.files(input.dir,pattern = pattern,recursive = T,full.names = T)
}

FDIR = function(file){
  basename(dirname(file))
}
geneFiltering. = function (exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * 
                             ncol(exprMat), minSamples = ncol(exprMat) * 0.01) 
{
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if (class(scenicOptions) == "ScenicOptions") {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  }
  else {
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if (is.null(dbFilePath)) 
    stop("dbFilePath")
  if (is.data.frame(exprMat)) {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", 
                                   "", methods("AUCell_buildRankings")), collapse = ", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    stop("'exprMat' should be one of the following classes: ", 
         supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if (any(table(rownames(exprMat)) > 1)) 
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  nCountsPerGene <- rowSums(exprMat, na.rm = T)
  nCellsPerGene <- rowSums(exprMat > 0, na.rm = T)
  message("Maximum value in the expression matrix: ", max(exprMat, 
                                                          na.rm = T))
  message("Ratio of detected vs non-detected: ", signif(sum(exprMat > 
                                                              0, na.rm = T)/sum(exprMat == 0, na.rm = T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  message("\nNumber of genes left after applying the following filters (sequential):")
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > 
                                                      minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", 
          minCountsPerGene)
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > 
                                                      minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", 
          minSamples, " cells")
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath)
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% 
                                                               genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  genesKept <- genesLeft_minCells_inDatabases
  if (!is.null(outFile_genesKept)) {
    saveRDS(genesKept, file = outFile_genesKept)
    if (getSettings(scenicOptions, "verbose")) 
      message("Gene list saved in ", outFile_genesKept)
  }
  return(genesKept)
}

viewMotifs. = function (tableSubset, motifCol = c("motif", "bestMotif", "MotifID"), 
                        dbVersion = "v9", nSignif = 3, 
                        colsToShow = c(motifEnrichment = c("motifDb", "logo", "NES", "geneSet", "TF_highConf"), 
                                       regulonTargets = c("TF", "gene", "nMotifs", "logo", "bestMotif","NES", "highConfAnnot", "Genie3Weight")),
                        options = list(pageLength = 50), ...) 
{
  if (!is.null(motifCol)) {
    motifCol <- motifCol[which(motifCol %in% colnames(tableSubset))]
    if (length(motifCol) == 1) {
      tableSubset <- RcisTarget::addLogo(tableSubset, 
                                         motifCol = motifCol, dbVersion = dbVersion, 
                                         addHTML = TRUE)
      if (!is.null(colsToShow)) 
        colsToShow <- c("logo", colsToShow)
    }
    else {
      stop("Please indicate the column containing the motif id (argument 'motifCol') or set it to NULL.")
    }
    tableSubset <- tableSubset[grep("transfac_pro__", tableSubset[[motifCol]], 
                                    invert = T), ]
  }
  for (i in which(sapply(tableSubset, is.numeric))) {
    tableSubset[[i]] <- signif(tableSubset[[i]], nSignif)
  }
  if (!is.null(colsToShow)) {
    colsToShow <- unique(unname(unlist(colsToShow)))
    colsToShow <- colsToShow[which(colsToShow %in% colnames(tableSubset))]
    tableSubset <- tableSubset[, colsToShow, with = F]
  }
  
  # DT::datatable(tableSubset, escape = FALSE, filter = "top", options = options)
  
  tableSubset$logo = tableSubset$logo %>% gsub(pattern = ".*(http.*png).*",replacement = "\\1")
  
  
  library(kableExtra)
  tableSubset = tableSubset[,c("gene","TF","nMotifs","NES","bestMotif","logo","highConfAnnot","Genie3Weight")]
  tableSubset = tableSubset[!is.na(Genie3Weight)&order(gene),]
  logo = tableSubset$logo
  tableSubset$logo = ""
  tableSubset %>% kbl() %>%
    kable_styling("striped", full_width = T, font_size = 15,
                  html_font = "Glyphicons Halflings") %>%
    row_spec(0, bold = TRUE, background = '#E2EDFA',hline_after = T)  %>%
    column_spec(6, image = spec_image(path =logo ,res = 25,width = 40,height = 18))
}
call.SCENIC <- function(tmp.DIR,org,Cluster_diff.tmp,cluster,seurat_object,dbDir="/data/wqihao/Project.sc/seuratRun/tmp/cisTarget_databases"){
  logFile=file.path(tmp.DIR,"err.txt")
  write.table("",quote = F,row.names = F,col.names = F,file = logFile,append = F) 
  pwd = getwd()
  setwd(tmp.DIR)
  # 初始化SCENIC设置 
  # org="mgi" # or hgnc, or dmel
  # dbDir="/data/wqihao/Project.sc/seuratRun/tmp/cisTarget_databases" # RcisTarget databases location
  myDatasetTitle="SCENIC" # choose a name for your analysis
  checkNeedRun = !file.exists(list.files(tmp.DIR,pattern = "Step2_regulonTargetsInfo.tsv",full.names = T,recursive = T))
  if(length(checkNeedRun)==0){
    dir.create("SCENIC.tmp")
    setwd("SCENIC.tmp") # Or `knitr::opts_knit$set(root.dir = 'example_results/SCENIC_MouseBrain')` in the first chunk if running a notebook
    
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20) 
    # scenicOptions <- readRDS("int/scenicOptions.Rds")
    scenicOptions@settings$verbose <- TRUE
    scenicOptions@settings$nCores <- 20
    scenicOptions@settings$seed <- 123
    
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    aa= Cluster_diff.tmp
    
    if(nrow(aa)<=6){
      write.table(paste0("[error: diff number < 6]   No TF target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
      
    }else{
      # aa = subset(aa,p_val<0.05&(pct.1>0.5&pct.2<0.2))
      aa = subset(aa,p_val_adj<0.05)
      # seurat_object <- readRDS(file.path("/data/wqihao/Project.sc/resultsDir/mergeGP1.GP2/","seurat_object.filted.celltype.RDS"))
      exprMat. = seurat_object@assays$RNA@counts
      both = intersect(aa$gene,rownames(exprMat.))
      
      exprMat. = exprMat.[both,]
      
      
      genesKept <- geneFiltering.(exprMat = exprMat., scenicOptions=scenicOptions,
                                  minCountsPerGene=3*.01*ncol(exprMat.),
                                  minSamples=ncol(exprMat.)*.01)
      
      exprMat_filtered <- as.matrix(exprMat.[genesKept, ])
      
      
      # step1: 计算相关性---------------------#区分潜在的激活和抑制
      runCorrelation(as.matrix(exprMat_filtered), scenicOptions) 
      
      # step2: Run GENIE3---------------------
      exprMat_filtered <- log2(exprMat_filtered+1) 
      nexts = T
      tryCatch(
        { 
          runGenie3(exprMat_filtered, scenicOptions)
        },
        # warning = function(w) { cat("warning") },
        error = function(e) { 
          write.table(paste0("[error: runGenie3]   No TF target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
          nexts <<-  F
        }, 
        finally = {cat("finished")}
      )
      
      if(nexts){
        tryCatch(
          { 
            runSCENIC_1_coexNetwork2modules(scenicOptions)
          },
          # warning = function(w) { cat("warning") },
          error = function(e) { 
            write.table(paste0("[error: runSCENIC_1]  No TF target DE gene in :",cluster ),quote = F,row.names = F,col.names = F,file = logFile,append = T)
            nexts <<-  F
          }, 
          finally = {cat("finished")}
        )
        
        if(nexts){
          tryCatch(
            { 
              runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
            },
            # warning = function(w) { cat("warning") },
            error = function(e) { 
              write.table(paste0("[error: runSCENIC_2]  No TF target DE gene in :",cluster ),quote = F,row.names = F,col.names = F,file = logFile,append = T) 
              nexts <<-  F
            }, 
            finally = {cat("finished")}
          )
        }
      }
      
      if(nexts){
        # regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
        regulonTargetsInfo = data.table::fread("output/Step2_regulonTargetsInfo.tsv")
        tableSubset <- regulonTargetsInfo[highConfAnnot==TRUE]
        fwrite(tableSubset,file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")),quote = F,sep = "\t")
      }
      setwd(pwd)
    }

  }else{
    if(checkNeedRun){
      dir.create("SCENIC.tmp")
      setwd("SCENIC.tmp") # Or `knitr::opts_knit$set(root.dir = 'example_results/SCENIC_MouseBrain')` in the first chunk if running a notebook
      
      data(defaultDbNames)
      dbs <- defaultDbNames[[org]]
      scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20) 
      # scenicOptions <- readRDS("int/scenicOptions.Rds")
      scenicOptions@settings$verbose <- TRUE
      scenicOptions@settings$nCores <- 20
      scenicOptions@settings$seed <- 123
      
      saveRDS(scenicOptions, file="int/scenicOptions.Rds")
      aa= Cluster_diff.tmp
      # aa = subset(aa,p_val<0.05&(pct.1>0.5&pct.2<0.2))
      aa = subset(aa,p_val_adj<0.05)
      # seurat_object <- readRDS(file.path("/data/wqihao/Project.sc/resultsDir/mergeGP1.GP2/","seurat_object.filted.celltype.RDS"))
      exprMat. = seurat_object@assays$RNA@counts
      exprMat. = exprMat.[aa$gene,]
      
      genesKept <- geneFiltering.(exprMat = exprMat., scenicOptions=scenicOptions,
                                  minCountsPerGene=3*.01*ncol(exprMat.),
                                  minSamples=ncol(exprMat.)*.01)
      
      exprMat_filtered <- as.matrix(exprMat.[genesKept, ])
      
      
      # step1: 计算相关性---------------------#区分潜在的激活和抑制
      runCorrelation(as.matrix(exprMat_filtered), scenicOptions) 
      
      # step2: Run GENIE3---------------------
      exprMat_filtered <- log2(exprMat_filtered+1) 
      nexts = T
      tryCatch(
        { 
          runGenie3(exprMat_filtered, scenicOptions)
        },
        # warning = function(w) { cat("warning") },
        error = function(e) { 
          write.table(paste0("[error: runGenie3]   No TF target DE gene in :",cluster ),file = logFile,quote = F,row.names = F,col.names = F,append = T) 
          nexts <<-  F
        }, 
        finally = {cat("finished")}
      )
      
      if(nexts){
        tryCatch(
          { 
            runSCENIC_1_coexNetwork2modules(scenicOptions)
          },
          # warning = function(w) { cat("warning") },
          error = function(e) { 
            write.table(paste0("[error: runSCENIC_1]  No TF target DE gene in :",cluster ),quote = F,row.names = F,col.names = F,file = logFile,append = T)
            nexts <<-  F
          }, 
          finally = {cat("finished")}
        )
        
        if(nexts){
          tryCatch(
            { 
              runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
            },
            # warning = function(w) { cat("warning") },
            error = function(e) { 
              write.table(paste0("[error: runSCENIC_2]  No TF target DE gene in :",cluster ),quote = F,row.names = F,col.names = F,file = logFile,append = T) 
              nexts <<-  F
            }, 
            finally = {cat("finished")}
          )
        }
      }
      
      if(nexts){
        # regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
        regulonTargetsInfo = data.table::fread("output/Step2_regulonTargetsInfo.tsv")
        tableSubset <- regulonTargetsInfo[highConfAnnot==TRUE]
        fwrite(tableSubset,file.path(tmp.DIR,paste0(cluster,".predicted_TF.xls")),quote = F,sep = "\t")
      }
      setwd(pwd)
    }else{
      message("SCENIC already finished")
    }
  }
  
}
plot.ppi <- function(PPI,outdir,cluster){
  library(igraph)
  links <- data.frame(
    source=PPI$protein1,
    target=PPI$protein2,
    importance = PPI$combined_score/1000
  )
  nodes <- data.frame(node_name = unique(c(PPI$protein1,PPI$protein2)),res="protein")
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  deg <- igraph::degree(network, mode="all")
  deg.sort = sort(deg,decreasing = T)
  deg.top100 = head(deg.sort,100)
  # head(deg.sort)
  nodes <- subset(nodes,node_name%in%names(deg.top100))
  # dim(nodes)
  links <- subset(links,source%in%nodes$node_name & target%in%nodes$node_name)
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  deg <- igraph::degree(network, mode="all")
  nodes$nodes_degree=deg
  
  if(mean(deg)<5) {
    k = (12-2)/(max(deg)-min(deg))
    deg = 1+k*(deg-min(deg))
  }
  
  if(mean(deg)>5) {
    # [1,6]
    k = (12-2)/(max(deg)-min(deg))
    deg = 1+k*(deg-min(deg))
  }
  # plot it
  coul  <- RColorBrewer::brewer.pal(3, "Set1")[1:2]
  coul = structure(coul,names=c("protein"))
  my_color = coul[V(network)$res]
  set.seed(124)
  pdf(file = file.path(outdir,paste0(cluster,".propro.network.pdf")),width = 8,height =8 )
  plot.igraph(network,
              # xlim = 1:5, 
              # vertex.size=ifelse(V(network)$carac>5,8,2),
              # vertex.size=V(network)$carac,
              vertex.size=deg,
              edge.width=ifelse(E(network)$importance>0.8,E(network)$importance,0),
              # edge.lty="solid",
              edge.color="gray",
              edge.lty = 1, 
              # edge.curved=0.2  ,
              # edge.arrow.width=1000,
              vertex.color=my_color,
              vertex.frame.color = "gray", 
              vertex.frame.size = 0.5, 
              vertex.label.cex=0.45,
              vertex.label.dist=0.2,
              vertex.label.font=2,  
              vertex.label.family="Times",
              vertex.label.color="black",
              layout=layout.sphere
  )
  dev.off()
  data.table::fwrite(links,file = file.path(outdir,paste0(cluster,".propro.network.xls")),quote = F,sep = "\t",col.names = T)
  
  data.table::fwrite(nodes,file = file.path(outdir,paste0(cluster,".propro.node.xls")),quote = F,sep = "\t",col.names = T)
}
plot.TFBS <- function(TFBS,outdir,cluster){
  TFBS = subset(TFBS,Genie3Weight!="^$")
  library(igraph)
  links <- data.frame(
    source=TFBS$TF,
    target=TFBS$gene
  )
  nodes <- data.frame(node_name =c( unique(TFBS$TF), unique(TFBS$gene) ),res=c(rep("TF",length(unique(TFBS$TF))),rep("gene",length(unique(TFBS$gene)))))
  nodes = nodes[order(nodes$node_name,nodes$res),]
  nodes = nodes[!duplicated(nodes$node_name),]
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  deg <- igraph::degree(network, mode="all")
  deg.sort = sort(deg,decreasing = T)
  deg.top100 = head(deg.sort,100)
  # head(deg.sort)
  nodes <- subset(nodes,node_name%in%names(deg.top100))
  # dim(nodes)
  links <- subset(links,source%in%nodes$node_name & target%in%nodes$node_name)
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  deg <- igraph::degree(network, mode="all")
  nodes$degree = deg
  nodes = nodes[order(nodes$res,nodes$node_name),]
  if(mean(deg)<5) {
    k = (12-4)/(max(deg)-min(deg))
    deg = 4+k*(deg-min(deg))
  }
  
  if(mean(deg)>5) {
    # [1,6]
    k = (12-4)/(max(deg)-min(deg))
    deg = 4+k*(deg-min(deg))
  }
  # plot it
  coul  <- RColorBrewer::brewer.pal(3, "Set2")[1:2]
  coul = structure(coul,names=c("TF","gene"))
  my_color = coul[V(network)$res]
  set.seed(124)
  pdf(file = file.path(outdir,paste0(cluster,".TF.network.pdf")),width = 8,height =8 )
  plot.igraph(network,
              # xlim = 1:5, 
              # vertex.size=ifelse(V(network)$carac>5,8,2),
              # vertex.size=V(network)$carac,
              vertex.size=deg,
              # edge.width=ifelse(E(network)$importance>0.8,E(network)$importance,0),
              # edge.lty="solid",
              edge.color="gray",
              edge.lty = 1, 
              # edge.curved=0.2  ,
              # edge.arrow.width=1000,
              vertex.color=my_color,
              vertex.frame.color = "gray", 
              vertex.frame.size = 0.5, 
              vertex.label.cex=0.45,
              vertex.label.dist=0.2,
              vertex.label.font=2,  
              vertex.label.family="Times",
              vertex.label.color="black",
              layout=layout.sphere
  )
  legend(x=0.70, y=-0.75, 
         legend=c(names(coul)), 
         col = c(coul) , 
         bty = "n", pch=20 , pt.cex = 1.5, cex = 0.8,
         text.col="black" , horiz = F)
  
  dev.off()
  data.table::fwrite(links,file = file.path(outdir,paste0(cluster,".TF.network.xls")),quote = F,sep = "\t",col.names = T)
  data.table::fwrite(nodes,file = file.path(outdir,paste0(cluster,".TF.node.xls")),quote = F,sep = "\t",col.names = T)
}
list2gmt <- function(list){
  library(magrittr)
  rds = readRDS(list)
  list %<>% gsub(pattern = "rds",replacement = "gmt")
  tmp= data.frame(V1=names(rds),V2="na",V3=do.call(c,lapply(rds, paste,collapse = "\t")),stringsAsFactors = F)
  library(tidyr)
  aa = sprintf("%s\t%s\t%s",tmp$V1,tmp$V2,tmp$V3)
  writeLines(aa,con =list,sep = "\n")
}
xls2gct = function(xls,species,out,celltype){
  xls.table = data.table::fread(xls,stringsAsFactors = F,check.names = F)
  colnames(xls.table) %<>% gsub(pattern = " ",replacement = ".")
  if(species=="mmu"){
    library(org.Mm.eg.db)    
    xls.table$gene = mapIds(org.Mm.eg.db,keys = xls.table$gene,column = "ENTREZID",keytype = "SYMBOL")
    xls.table = xls.table[,c(-3,-4,-5,-6)]
    colnames(xls.table)[2] = "V2"
    xls.table$V2 = "NA"
    colnames(xls.table)[1:2] = c("NAME","Description")
    
    heads = c("#1.2",paste(nrow(xls.table),ncol(xls.table)-2,sep = "\t"))
    writeLines(heads,con = file.path(out,paste0(celltype,".gct")),sep = "\n")
    write.table(xls.table,file =file.path(out,paste0(celltype,".gct")),quote = F,append = T,sep = "\t",row.names = F,col.names = T)
  }
  if(species=="hsa"){
    library(org.Hs.eg.db)    
    xls.table$gene = mapIds(org.Hs.eg.db,keys = xls.table$gene,column = "ENTREZID",keytype = "SYMBOL")
    xls.table = xls.table[,c(-3,-4,-5,-6)]
    xls.table$p_val = "NA"
    colnames(xls.table)[1:2] = c("NAME","Description")
    heads = c("#1.2",paste(nrow(xls.table),ncol(xls.table)-2,sep = "\t"))
    writeLines(heads,con =file.path(out,paste0(celltype,".gct")),sep = "\n")
    write.table(xls.table,file =file.path(out,paste0(celltype,".gct")),quote = F,append = T,sep = "\t",row.names = F,col.names = T)
  }
  
  celltype.tmp = celltype %>% gsub(pattern = "(cluster)|(celltype\\.)",replacement = "",ignore.case = T)
  heads = c(paste(ncol(xls.table)-2,2,1,sep = "\t"),
            paste("#",celltype,"others",sep = "\t"),
            paste(ifelse(colnames(xls.table)[c(-1,-2)]!=celltype.tmp,0,1),collapse = "\t")
  )
  writeLines(heads,con = file.path(out,paste0(celltype,".cls")),sep = "\n")
}


flage <- function(pipeline.log,step){
  invisible(file.create(GetoptLong::qq("@{pipeline.log}/@{step}.flage"),showWarnings=F))
}
check.flage <- function(pipeline.log,step){
  !file.exists(GetoptLong::qq("@{pipeline.log}/@{step}.flage"))
}


makeNames <- function(n){
  n = gsub(pattern ="[ ]+",replacement = ".",n)
  
  return(n)
}


 
