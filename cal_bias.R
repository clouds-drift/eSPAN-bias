required.pack = c("getopt","optparse","openxlsx","reshape2","zoo","ggplot2","plyr","pheatmap","RColorBrewer","officer","rvg")
for(p in required.pack){
  sig=suppressMessages(require(p,character.only = TRUE,quietly = T))
  #if(!sig){stop(paste("package",p,"needed!"))}
  if(!sig){
    cat("install",p,"...\n")
    install.packages(p)
  }
}
bioc.pack=c("GenomicRanges","Rsamtools","rtracklayer")
for(p in bioc.pack){
  sig=suppressMessages(require(p,character.only = TRUE,quietly = T))
  if(!sig){
    cat("install",p,"...\n")
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(p)
  }
}
opt.p=OptionParser(description="Calculate read coverage bias between plus and minus strands of bigwig files.eSPAN experiments will generate strand-specific DNA libraries.Separated read coverage on plus and minus strands can be used for bias calculation.")
opt.p=add_option(opt.p,c("--watson","-w"),type="character",
                 help="Bigwig file of plus strand. (Required) 
                 Multiple files can be joined by ',', like 'sample1_plus.bw, sample2_plus.bw'.")
opt.p=add_option(opt.p,c("--crick","-c"),type="character",
                 help="Bigwig file of minus strand. (Required)
                 Multiple files can be joined by ',', like 'sample1_minus.bw, sample2_minus.bw'. Should have same length as '--watson'.")
opt.p=add_option(opt.p,c("--reference_point","-r"),type="character",
                 help="A BED file of interested annotations and the middle point will be aligned. (Required)
                 Default is the yeast G1 replication origins.")
opt.p=add_option(opt.p,c("--method","-m"),type="character",
                 help="Method for bias calculation, should be 'partition' or 'log'. (Optional, default='partition')
                 'partition' = (w-c) / (w+c)
                 'log' = log2 (w / c)", default="partition")
opt.p=add_option(opt.p,c("--bw"),type="character",
                 help="Output folder for bigwig signals of genome-wide bias.It can be loaded by IGV to check individual locus. (Optional, default='bw_bias')",
                 default="bw_bias")
opt.p=add_option(opt.p,c("--matrix","-o"),type="character",
                 help="Output folder for bias matrix. Will be used for bias normalization or bias plot. (Optional, default='bias_matrix')",
                 default="bias_matrix")
opt.p=add_option(opt.p,c("--sumMatFile", "-s"),type="character",
                 help="Summary of output matrix file. (Optional, default='my_bias_mat.xlsx')",
                 default="my_bias_mat.xlsx")
opt.p=add_option(opt.p,c("--flank"),type="character",
                 help="Flanking regions around the reference point. It should be two numeric values joined by ','. (Optional, default='-2000,2000')
                 Positive value means downstream and negative value means upstream to the reference point.",default="-2000,2000")
opt.p=add_option(opt.p,c("--bin","-b"),type="integer",help="Bin size for matrix calculation. (Optional, default=100)",default=100)
opt.p=add_option(opt.p,c("--smooth"),type="integer",help="Number of bins on each side will be used for smooth. (Optional, default=0)",default=0)
opt.p=add_option(opt.p,c("--threshold","-t"),type="integer",help="Regions with coverage less than this value will be sikpped. (Optional, default=0)",default=0)
opt.p=add_option(opt.p,c("--N_thread","-N"),type="integer",help="Number of parallel processing. (Optional, default=2)",default=2)
opt.p=add_option(opt.p,c("--force","-f"),type="logical", help="If not forced, existing files will be skipped for calculation again. (Optional, default=T)",default=T)
opt=parse_args(opt.p)
if(is.null(opt$watson)|is.null(opt$crick)){print_help(opt.p);stop()}
print(opt)

where=function(){
  spath <-parent.frame(2)$ofile
  
  if (is.null(spath)) {
    args <- commandArgs()
    filearg <- args[grep("^--file=", args)]
    fname <- strsplit(filearg, "=")[[1]][2]
  } else {
    fname <- spath
  }
  
  dirname(normalizePath(fname))
}
source(file.path(where(),"source/bias_function.R"))
##initialize value
w=unlist(strsplit(opt$watson,","));w=trimws(w)
c=unlist(strsplit(opt$crick,","));c=trimws(c)
if(is.null(opt$reference_point)){
  reference.point=file.path(where(),"annotation/saccer3_G1_origin.bed")
}else{reference.point=opt$reference_point}
method=opt$method
bw.bias.dir=opt$bw
bias.raw.dir=opt$biasRaw
bias.dir=opt$matrix
sum.mat.file=opt$sumMatFile
flank=as.numeric(unlist(strsplit(opt$flank,",")))
bin=opt$bin
sm=opt$smooth
threshold=opt$threshold
N_thread=opt$N_thread
force=opt$force

cat("Cal bias in smooth line...\n")
bias.smooth.df=Bias_in_smooth( w=w,c=c,reference.point=reference.point,method=method,
            bw.bias.dir=bw.bias.dir,bias.dir=bias.dir,sum.mat.file=sum.mat.file,
            flank=c(flank[1],flank[2]+bin),bin=bin,sm=sm,N_thread=N_thread,threshold=threshold,force=force)


