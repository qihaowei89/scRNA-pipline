#! /usr/bin/env nextflow

/*
 * scPipe: A nextflow-based single cell analysis pipeline
 * by    : Qihao Wei
 * data  : 2021-06-07
*/

//pre-defined functions for render command

ANSI_RESET = "\u001B[0m";
ANSI_BLACK = "\u001B[30m";
ANSI_RED = "\u001B[31m";
ANSI_GREEN = "\u001B[32m";
ANSI_YELLOW = "\u001B[33m";
ANSI_BLUE = "\u001B[34m";
ANSI_PURPLE = "\u001B[35m";
ANSI_CYAN = "\u001B[36m";
ANSI_WHITE = "\u001B[37m";

def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }

//=======================================================================================
// Nextflow  version
version="1.0"

// Nextflow Version check
if( !nextflow.version.matches('21.04.0') ) {
    println print_yellow("This workflow requires Nextflow version 19.0.0 or greater -- You are running version ") + 
    print_red(nextflow.version)
}


//help information
params.help = null
if (params.help) {
    println ''
    println print_purple('-------------------------------------------------------------------------' )
    println print_purple("   scPipe: A nextflow-based single cell analysis pipeline v$version      " )
    println print_purple('-------------------------------------------------------------------------' )
    println print_yellow('Usage: ')
    println print_yellow('    The typical command for running the pipeline is as follows (we do not\n' ) + 
    print_yellow        ('    recommend users passing configuration parameters through command line,\n') +
    print_yellow        ('    please modify the config.file instead):\n'                               ) +
    print_purple        ('    Nextflow run main.nf '                                                   ) +
    print_yellow        ('    General arguments:Input and output setting\n'                            ) +
    print_cyan          ('      --run   <run1,run2>\n'                                                 ) + 
    print_green         ('        input sample(s) OR run(s) name\n'                                    ) +
    
    print_cyan          ('      --csv   <csv1,csv2>\n'                                                 ) + 
    print_green         ('        Path to input cellranger multi config.csv(s)\n'                      ) +
    print_cyan          ('      --fastqdirs   <path1,path2>\n'                                         ) + 
    print_green         ('        Path to input fastq file(s)\n'                                       ) +
    print_cyan          ('      --reads   <path1,path2>\n'                                             ) + 
    print_green         ('        Path to input data(s)\n'                                             ) +
    print_cyan          ('      --species   <species>\n'                                               ) + 
    print_green         ('        mmu or hsa           \n'                                             ) +
    print_cyan          ('      --DBdir   <DBdir>\n'                                                   ) + 
    print_green         ('        database dir for singleR annotation celltypes\n'                     ) +
    print_cyan          ('      --outdir   <path>\n'                                                   ) + 
    print_green         ('        Path to output files\n'                                              ) +
    print_cyan          ('      --skipAggr  <Boolean>\n'                                               ) + 
    print_green         ('        true:skip cellranger aggr; false: run cellranger aggr\n'             )
    println print_purple('-------------------------------------------------------------------------'   )
    exit 0
}


/*
 * Step1: Cell Ranger 
 *        1. task Run_cellranger_mkfastq
 *        2. task Run_cellranger_count
 *        3. task Run_cellranger_vdj
 *        4. Or Run_cellranger_multi
 */

runs       = Channel.from( params.run.split(','))
csvs       = Channel.from( params.csv.split(','))
channelObj = runs.merge( csvs )

runs2       = Channel.from( params.run.split(','))
fastqs     = Channel.from( params.fastqdirs.split(','))
singlerObj = runs2.merge (fastqs)  


process RunCellrangerMulti {
  publishDir pattern: "*", path: "${params.outdir}", mode: 'copy'
  tag "$run"
  
  input:
  set (run, csv) from channelObj 
  
  output:
  file ("*") into ValidateOutput1s,GeneratAggrCsvObj

  shell:
  """
  cellranger multi --id $run --csv $csv 
  """
}

process GeneratAggrCsv {
  publishDir pattern: "*", path: "${params.outdir}", mode: 'copy'
  
  when: !params.skipAggr
  
  input:
  file h5s from GeneratAggrCsvObj.collect()
  output:
  file 'aggr.csv' into CellrangerAggrObj
  
  shell:
  """
  ${params.bin}/aggr.R ${params.outdir} aggr.csv
  """
}

process RunCellrangerAggr{
  publishDir pattern: "*", path: "${params.outdir}", mode: 'copy'
  tag "merged"
  
  when: !params.skipAggr
  
  input:
  file aggrCsv from CellrangerAggrObj 
  
  output:
  file("*") into ValidateOutput2
  
  shell:
  """
  cellranger aggr --id merged --csv $aggrCsv 
  """
}

process RunSingleR {
  publishDir pattern: "*", path: "${params.outdir}/report", mode: 'copy'
  
  input:
  file aggrCsv  from ValidateOutput1s 
  set (run, fastqdir) from  singlerObj
  
  output:
  set val(csv), file("*") into ValidateOutput3
  
  shell:
    
  """
  Rscript ${params.bin}/singleSample.R  --fastqDir $fastqdir --inputDir ${params.outdir}/$run --outDir ./ --species hsa --annoDB BlueprintEncodeData --DBdir ${params.DBdir}
  """
}

process RunMergR {
  publishDir pattern: "*", path: "${params.outdir}/report", mode: 'copy'
  
  input:
  file aggrCsv  from ValidateOutput2 
  
  output:
  set val(csv), file("*") into ValidateOutput4
  
  shell:
    
  """
  Rscript ${params.bin}/mergR.R  --inputDir ${params.outdir}/merged --outDir ./ --species hsa  --subdivide ${params.run} --annoDB BlueprintEncodeData
  """
}


workflow.onComplete {                                                           
  log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/pipeline.report.html\n $params.outdir/report/singlecell_report.html\n" : "Oops .. something went wrong" )
} 

