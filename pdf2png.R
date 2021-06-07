Cmd <- function(cmd,cl=8,step,f=System,out) 
{ 
  packages <- c("doMC","iterators")
  invisible(purrr::map(packages,require,character.only = T))
  
  # check if out file exists
  cmd_need_run <- cmd[!file.exists(out)]
  out_need_get <- out[!file.exists(out)]
  doMC::registerDoMC(cl)
  if (length(cmd_need_run)!=0){
    cat(sprintf("[ %s ] Start %10s total(%s), run(%s), skip(%s)\n",Sys.time(),step,length(cmd),length(cmd_need_run),length(cmd)-length(cmd_need_run)))
    stat = foreach(a=iter(cmd_need_run),.combine = c) %dopar% f(a)
    
    err_cmd <- which(stat!=0)
    if(length(err_cmd!=0)) system(sprintf("rm -rf %s",out_need_get[err_cmd]))
    cat(sprintf("[ %s ] End   %10s finished(%s), failed(%s)\n",Sys.time(),step,length(cmd_need_run)-length(err_cmd),length(err_cmd)))
  }else{
    cat(sprintf("[ %s ] Skip  %10s\n",Sys.time(),step))
  }
}


OUT.DIR = "/mnt/d/百度云同步盘/国药/单细胞相关/210001_01_单样本分析结果_02210517/分析结果/mergeGP1.GP2/"

all.pdf = list.files(OUT.DIR,pattern = ".pdf$",recursive = T,full.names = T)
all.png = gsub(pattern = "pdf",replacement = "png",all.pdf)
cmd = sprintf("/usr/bin/convert %s %s",all.pdf,all.png)
Cmd(cmd = cmd,cl =16 ,out = all.png,step = "aaa",f = system)

library(magrittr)
OUT.DIR = "/mnt/d/百度云同步盘/国药/单细胞相关/210001_01_单样本分析结果_02210517/分析结果/mergeGP1.GP2/6.GSEA/"
f.name = list.dirs(OUT.DIR,recursive = T) %>% grep(pattern = "Gsea",value = T) %>% grep(pattern = "edb",invert = T,value = T)
t.name = f.name %>% gsub(pattern ="(.*Gsea)(\\.\\d*)$" ,replacement = "\\1")
file.rename(from = f.name,to = t.name)

OUT.DIR = "/mnt/d/百度云同步盘/国药/单细胞相关/210001_01_单样本分析结果_02210517/分析结果/mergeGP1.GP2/"
f.name = list.dirs(OUT.DIR,recursive = T) %>% grep(pattern = "SCENIC.tmp$",value = T) 

cmd = sprintf("rm -rf  %s",f.name)
Cmd(cmd = cmd,cl =16 ,out = as.character(1:30),step = "aaa",f = system)

