#!/usr/bin/env Rscript  
args = commandArgs(T)

h5Files = list.files(path = args[1],pattern = "sample_molecule_info.h5",recursive = T,full.names = T)
samples = gsub(pattern = '.*per_sample_outs/(.*)/count.*',replacement = '\\1',h5Files)
csvTab = data.frame(sample_id=samples,molecule_h5=h5Files)
write.csv(csvTab,file = args[2],quote = F,row.names = F)

