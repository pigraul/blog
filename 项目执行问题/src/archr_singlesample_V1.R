#conda activate ArchR
#options(synchronous = NULL)
#library(DBI)
#library(RSQLite)
#on.exit({
#  lapply(dbListConnections(SQLite()), dbDisconnect)
#})
library(readr)
library(argparser)
library(ArchR)
library(pheatmap,lib.loc="/SGRNJ01/Public/Software/conda_env/ArchR/lib/R/library/")
library(AnnotationDbi, lib.loc = "/Personal/wangxiang/R_lib/")
library(Seurat)
library(tidyr)
library(dplyr)
library(chromVARmotifs)
library(yaml)  
library(purrr)   

## parameters
argv <- arg_parser('')
argv <- add_argument (argv,"--fragment", help="the cell gene express file list")
argv <- add_argument (argv,"--spname",help="the samples name list")
argv <- add_argument (argv,"--threads",help="threads used",default = 1)
argv <- add_argument (argv,"--species",help="hg19 or hg38 or mm10 or mm9",default = "hg38")
argv <- add_argument (argv,"--genomeSize",help="effective genome size", default = NULL)
argv <- add_argument (argv,"--outdir",help="path of outdir")
argv <- add_argument (argv,"--archrdir",help="path of archr project dir")
argv <- add_argument (argv,"--prefix",help="output file prefix")
argv <- add_argument (argv,"--resolution", help="the cluster resolution,default=1.2; if ncell < 1500, the default resolution would be 0.8",default=1.2)
argv <- add_argument (argv,"--dim", help="use cluster number 1:20",default=20)
argv <- add_argument (argv,"--minTSS", help="min signalvsbackground ratio,default=1",default=0.1)
argv <- add_argument (argv,"--minFrags", help="min unique fragment number,default=100",default=200)
argv <- add_argument (argv,"--scRNA", help ="path of scRNA rds or no",default = "no")
argv <- add_argument (argv,"--testMethod",help = "wilcoxon, ttest or binomial",default = "wilcoxon")
argv <- add_argument (argv,"--peakcells", help="call peaks for all cluster with more than 50 cells",default=50)
argv <- add_argument (argv,"--extraLib", help="extra R lib dir", default="/Personal/wangxiang/R_lib/")
argv <- add_argument (argv,"--fdr",help="fdr from step4, required if start from step5", default = NULL)
argv <- add_argument (argv,"--startPoint",default=0,help="0:create,1:basic analysis,2:Dimensionality Reduction,3:scRNA-seq based annotation,4:Identify peaks and plot marker peak,5:Motif analysis,6:Co-accessibility and Peak2GeneLinkage,7:Trajectory Analysis(TBC)")

argv <- parse_args(argv)


#//                          _ooOoo_                               //
#//                         o8888888o                              //
#//                         88" . "88                              //
#//                         (| ^_^ |)                              //
#//                         O\  =  /O                              //
#//                      ____/`---'\____                           //
#//                    .'  \\|     |//  `.                         //
#//                   /  \\|||  :  |||//  \                        //
#//                  /  _||||| -:- |||||-  \                       //
#//                  |   | \\\  -  /// |   |                       //
#//                  | \_|  ''\---/''  |   |                       //
#//                  \  .-\__  `-`  ___/-. /                       //
#//                ___`. .'  /--.--\  `. . ___                     //
#//              ."" '<  `.___\_<|>_/___.'  >'"".                  //
#//            | | :  `- \`.;`\ _ /`;.`/ - ` : | |                 //
#//            \  \ `-.   \_ __\ /__ _/   .-` /  /                 //
#//      ========`-.____`-.___\_____/___.-`____.-'========         //
#//                           `=---='                              //
#//      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        //
#//             佛祖保佑       永不宕机     永无BUG                 //

outdir <- argv$outdir
dir.create (outdir,showWarnings=T)
archrdir <- argv$archrdir
resol <- as.numeric(argv$resolution)
cnum <- as.numeric(argv$dim)
start_point <- as.numeric(argv$startPoint)
scRNA <- argv$scRNA
minTSS<-argv$minTSS
minFrags<-argv$minFrags
testmethod <- argv$testMethod
prefix <- argv$prefix
peakfilter <- as.numeric(argv$peakcells)
elib <- argv$extraLib
#filterdoublet <- argv$filterdoublet
#obj <- as.numeric(argv$object)
# wd <- getwd()

# set seed for randomized opperations for consistency
set.seed(1)
#addArchRLocking(locking = TRUE)
# Setting default number of Parallel threads.
addArchRThreads(threads = as.numeric(argv$threads)) 


dirnameScript <- function(){
    # get full directory path of current script located
    cmd = commandArgs(trailingOnly = FALSE)
    scriptName = sub('--file=', "", cmd[grep('^--file=', cmd)])
    if (length(scriptName) > 0) {
        path = normalizePath(scriptName)
        dirname = dirname(path)
    } else {
        print('Warning: not a runtime environment, using current directory instead.')
        dirname = getwd()
    }
    return(dirname)
}
scriptDirname = dirnameScript()



# source 
yaml_file <- paste0(scriptDirname,'/species_config.yml') 
config <- yaml.load_file(yaml_file)
funcR<-paste0(scriptDirname,'/archr_func.R') 
source(funcR)


# Check species
# Setting a Genome and GeneAnnotation
species_Builtin <- c("hg19", "hg38", "mm10", "mm9")
genome_dict <- list(
  "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
  "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
  "mm9" = "BSgenome.Mmusculus.UCSC.mm9",
  "mm10" = "BSgenome.Mmusculus.UCSC.mm10"
)
species_cfg <- names(config)

if (argv$species %in% species_Builtin) {
    # Built-in species
    message("[INFO] Species '", argv$species, "' is a built-in species. Using default annotations.") 
    addArchRGenome(argv$species)
} else if (argv$species %in% names(config)) {
    # Config species
    message("[INFO] Species '", argv$species, "' found in config file. Loading custom annotations.")
    species_config <- config[[argv$species]]
    packages_to_load <- unlist(species_config, use.names = FALSE)
    for (pkg in packages_to_load) {
        library(pkg, character.only = TRUE, lib.loc=elib)
    }
    # Creating a Custom ArchRGenome
    genomeAnnotation <- createGenomeAnnotation(genome = get(species_config$BSgenome))
    geneAnnotation <- createGeneAnnotation(
        TxDb = get(species_config$TxDb), 
        OrgDb = get(species_config$OrgDb)
    )
} else {
    # species not found
    stop("[ERROR] Species '", argv$species, "' is not supported. ",
        "Check spelling or add it to the config file. ",
        "Available options: ", 
        paste(c(species_Builtin, names(config)), collapse = ", "))
}

## create output dirs
dir.create (file.path(outdir,"/Plots"), showWarnings=T) 
outrds <- file.path (outdir,"0.rds") 
dir.create (outrds,showWarnings = T)


## step 0 or load
if(start_point == 0){
    message("[INFO] Start from step ", start_point, " Creating Arrow Files and ArchRProject.")
    inputFiles <- c(argv$fragment)
    samplename <- argv$spname
    names (inputFiles) <- samplename

    if (argv$species %in% species_Builtin){
        ArrowFiles <- createArrow(argv$species, inputFiles, minTSS, minFrags)
        projHeme0 <- ArchRProject(
            ArrowFiles = ArrowFiles, 
            outputDirectory = file.path(outdir),
            copyArrows = TRUE  #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
        )
    } else {
        ArrowFiles <- createArrow(argv$species, inputFiles, minTSS, minFrags,
                   genomeAnnotation = genomeAnnotation, geneAnnotation = geneAnnotation) 
        projHeme0 <- ArchRProject(
            ArrowFiles = ArrowFiles,
            outputDirectory = file.path(outdir),
            copyArrows = TRUE,
            geneAnnotation = geneAnnotation,
            genomeAnnotation = genomeAnnotation,
        )
    }
    saveArchRProject(ArchRProj = projHeme0, load = FALSE)
    message("[INFO] Create Arrow Files and ArchRProject DONE!")
} else{
    message("[INFO] Start from step ", start_point, " Loading ArchRProject from ", archrdir)
    pre_no = start_point - 1
    var_name <- paste0("projHeme", pre_no)
    assign(var_name, loadArchRProject(path = archrdir))  
    if ( ! identical(normalizePath(outdir), normalizePath(archrdir)) ){
        message("[INFO] Save ArchRProject from ", archrdir ," to ", outdir, ".")
        assign(var_name, saveArchRProject(ArchRProj = get(var_name), outputDirectory = outdir, load = TRUE))
    }
    message("[INFO] Loading ArchRProject DONE!")
} 


x <- c(0, 1, 2, 3, 4, 5, 6, 7)
run_steps <- x[x >= start_point]  

###################################-------------------------------------------projHeme1
### Per-cell Quality Control (TSS enrichment score greater than 4 and more than 1000 unique nuclear fragments)
if( 1 %in% run_steps){
    message("[INFO] Starting execution of 1 Per-cell Quality Control.")
    projHeme1 <- step1_Percell_Quality_Control(projHeme0, outdir)
    message("[INFO] Per-cell Quality Control execution finished.")
}

###################################-------------------------------------------projHeme2
### Dimensionality Reduction
if( 2 %in% run_steps){
    message("[INFO] Starting execution of 2 Dimensionality Reduction.")
    projHeme2 <- step2_Dimensionality_Reduction(projHeme1, outdir, prefix, resol=resol, dims=cnum)
    message("[INFO] Dimensionality Reduction execution finished.")
}

###################################-------------------------------------------projHeme3
### Cell Annotation
if( 3 %in% run_steps){
    message("[INFO] Starting execution of 3 Cell Annotation.")
    projHeme3 <- step3_Cell_Annotation(projHeme2, outdir, prefix, scRNA=scRNA)
    message("[INFO] Cell Annotation execution finished.")
}

###################################-------------------------------------------projHeme4
### Calling Peaks with ArchR
if( 4 %in% run_steps){
    message("[INFO] Starting execution of 4 Peak Calling.")
    ## reload genome libraries
    if(argv$species %in% species_Builtin) {
        pkg <- genome_dict[[argv$species]]
        library(pkg, character.only = TRUE)
    }else if (argv$species %in% names(config)) {
        for (pkg in packages_to_load) {
            library(pkg, character.only = TRUE, lib.loc=elib)
        }
        # load pacakage in extra lib dir
        original_libpaths <- .libPaths()
        .libPaths(c(elib, original_libpaths))
    }

    if (!is.null(argv$genomeSize) && !is.na(argv$genomeSize) && argv$genomeSize != "") {
        genome_size <- as.numeric(argv$genomeSize)
        message("Using provided genome size: ", genome_size)
    } else {
        message("--genomeSize not provided. Using default calculation.")
        if (argv$species %in% species_Builtin){
            genome_size <- NULL
        } else {
            genome_size <- get_genome_size(projHeme3)
            message("Calculated genome size: ", genome_size)
        }
    }

    result <- step4_Call_Peaks(projHeme3, outdir, prefix, genomeSize=genome_size, minCells=peakfilter, testmethod=testmethod)
    projHeme4 <- result$archr_proj        
    FDR <- result$fdr
    message("[INFO] Peak Calling execution finished.")
}

###################################-------------------------------------------projHeme5
#### Motif and Feature Enrichment with ArchR
if( 5 %in% run_steps){
    message("[INFO] Starting execution of 5 Motif and Feature Enrichment.")
    ## check fdr
    if (exists("FDR")) {
        message("[INFO] Using FDR from step 4: ", FDR)
    } else if (!is.null(argv$fdr) && !is.na(argv$fdr) && argv$fdr != "") {
        FDR <- as.numeric(argv$fdr)
        message("[INFO] Using FDR from provided: ", FDR)
    } else {
        stop("[ERROR] --fdr not provided !")
    }
    ## reload genome libraries
    if(argv$species %in% species_Builtin) {
        pkg <- genome_dict[[argv$species]]
        library(pkg, character.only = TRUE)
    }else if (argv$species %in% names(config)) {
        for (pkg in packages_to_load) {
            library(pkg, character.only = TRUE, lib.loc=elib)
        }
        # load pacakage in extra lib dir
        original_libpaths <- .libPaths()
        .libPaths(c(elib, original_libpaths))
    }

    projHeme5 <- step5_Motif(projHeme4, outdir, prefix, FDR, testmethod=testmethod)
    message("[INFO] Motif and Feature Enrichment execution finished.")
}

###################################-------------------------------------------projHeme6
### Co-accessibility with ArchR
if( 6 %in% run_steps){
    message("[INFO] Starting execution of 6 Co-accessibility.")
    projHeme6 <- step6_Co_Accessibility(projHeme5, outdir, prefix)
    message("[INFO] Co-accessibility execution finished.")
}


















