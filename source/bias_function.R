required.pack = c("openxlsx","GenomicRanges","Rsamtools","rtracklayer","reshape2","zoo"
)
for(p in required.pack){
  sig=suppressMessages(require(p,character.only = TRUE))
  if(!sig){stop(paste("package",p,"needed!"))}
}


Bias_in_smooth=function(w,c,reference.point,method="partition",
                        bw.bias.dir="bw_bias",bias.dir="bias_matrix",sum.mat.file="my_bias_mat.xlsx",
                        flank=c(-2000,2000),bin=100,sm=0,N_thread=2,threshold=0,force=F){
  if(!dir.exists(bw.bias.dir)){dir.create(bw.bias.dir,recursive = T)}
  if(!dir.exists(bias.dir)){dir.create(bias.dir,recursive = T)}
  
  cat("\ncal bias bw...\n")
  bias.df=cal_partition_bw(w=w, c=c, outdir=bw.bias.dir,threshold=threshold,force=force,method=method)
  cat("\ncal bias matrix...\n")
  infile=bias.df$File
  names(infile)=bias.df$Sample
  bias.mat.file=profile_point(infile, reference.point,flanking=flank, bs=bin, outdir=bias.dir, N_thread=N_thread,force=force) ##get profile with n bp
  bias.mat.df=read_deeptools_mat(sample=names(infile) ,infile=bias.mat.file, outdir=bias.dir, flanking=flank, bs=bin,force=force)
  mat.file=bias.mat.df$Mat
  names(mat.file)=bias.mat.df$Sample
  cat("\nsmooth bias matrix...\n")
  bias.smooth.df=cal_nuc_profile(mat.file, outdir=bias.dir, bs=sm,fill=NA)##smooth profile with neighbour n bins
  write.xlsx(bias.smooth.df, sum.mat.file,overwrite=T)
  sapply(bias.mat.file,function(x){if(file.exists(x)) file.remove(x)})
  
  return(bias.smooth.df)
}

##sum reads of two strands need > threshold
cal_partition_bw=function(w, c, outdir,threshold=10, force=F,method="partition"){
  library("GenomicRanges")
  library("Rsamtools")
  library("rtracklayer")
  outdir=sub("/$", "", outdir)
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  if(length(w) != length(c)){
    cat("watson and crick files are not paired!")
    stop()
  }
  
  out.df=c()
  out.num=0
  for (i in 1:length(w)) {
    out.name=string_overlap(basename(w[i]), basename(c[i]))
    out.name=sub("_$", "", out.name)
    bw.file=paste(outdir, "/", out.name, "_partion.bw", sep="")
    one.df=data.frame(Sample=out.name, File=bw.file, stringsAsFactors = F)
    out.num=out.num+1
    cat(out.num, "cal bias of", out.name, "\n")
    if(file.exists(bw.file) & !force){
      cat(bw.file, "exists\n")
    }else{
      w.cov=import(w[i], as="RleList")
      c.cov=import(c[i], as="RleList")
      if(grepl("partition",method,ignore.case = T)){
        bias.cov=(abs(w.cov)-abs(c.cov)) / (abs(w.cov)+abs(c.cov))
      }else if(grepl("log",method,ignore.case = T)){
        bias.cov=log2( abs(w.cov) / (abs(c.cov)+0.1) )
      }else{
        stop("unknow method:",method,"!\n")
      }
      
      bias.cov[(abs(c.cov)+abs(w.cov)) < threshold]=NA
      export.bw(bias.cov, bw.file, format="bw")
    }
    out.df=rbind(out.df, one.df, stringsAsFactors=F)
  }
  return(out.df)
}


##deeptools required
profile_point=function(bw, bed, flanking=c(-100,100), bs=10, outdir, N_thread=2, yMin="", yMax="", zMin="", zMax="", N_color=5, kmeans="",sort=T, force=FALSE){
  if(!dir.exists(outdir)){dir.create(outdir, recursive = T)}
  ##check for input file
  if(is.null(names(bw))){names(bw)=sub("\\..*$","",basename(bw))}
  bw=sapply(bw,function(x){if(grepl("\\.bw$",x)){return(x)}else{return(paste(x,"/",basename(x),".bw",sep=""))}})
  file.ind=file.exists(bw); if(!all(file.ind)){stop(paste(bw[!file.ind],collapse = " "),"not found!\n")}
  if(nchar(yMin)==0){arg.ymin=""}else{arg.ymin=paste("--yMin ",yMin)}
  if(nchar(yMax)==0){arg.ymax=""}else{arg.ymax=paste("--yMax ",yMax)}
  if(nchar(zMin)==0){arg.zmin=""}else{arg.zmin=paste("--zMin ",zMin)}
  if(nchar(zMax)==0){arg.zmax=""}else{arg.zmax=paste("--zMax ",zMax)}
  if(nchar(kmeans)==0){arg.kmeans=""}else{arg.kmeans=paste("--kmeans ", kmeans)}
  if(sort==T){arg.sort=""}else{arg.sort="--sortRegions no"}
  if(length(bw)!=length(bed)){bed=rep(bed[1],length(bw));cat(bed[1],"is aligned...\n")}
  flanking=as.numeric(flanking)
  
  result=c()
  for(i in 1:length(bw)){
    cat(i, bw[i], "profile\n")
    out.name=names(bw)[i]
    if(is.null(out.name)){out.name=sub("(\\.bw)|(\\.bigWig)$", "",  basename(bw[i]), ignore.case = T)}
    outfile=paste(outdir, "/",out.name, "_on_", sub("\\.bed$","",basename(bed[i])), ".mat.gz", sep="")
    final.file=paste(outdir,"/",out.name,".xlsx",sep="")
    out.region=paste(outdir, "/",out.name, "_on_", sub("\\.bed$","",basename(bed[i])), "_region.txt", sep="")
    graph.file=paste(outdir, "/", out.name, "_on_", sub("\\.bed$","",basename(bed[i])), ".pdf", sep="")
    graph.file2=paste(outdir, "/", out.name, "_on_", sub("\\.bed$","",basename(bed[i])), "_heatmap.pdf", sep="")
    if(!file.exists(final.file) | force==TRUE){
      cmd=paste("computeMatrix  reference-point", "-S", paste("'",bw[i],"'",sep=""), "-R", paste(paste("'",bed[i],"'",sep=""), collapse = " "),
                "-b", sprintf("%d", abs(flanking[1])),
                "--referencePoint center",
                "-a", sprintf("%d", abs(flanking[2])),
                "-bs", bs,
                "-o", paste("'",outfile,"'",sep=""),
                #"--outFileNameMatrix", outfile2,
                "-p", N_thread,
                #"--missingDataAsZero",
                sep=" ")
      system(cmd)
    }else{
      cat(outfile, "exists\n")
    }
    if(!file.exists(graph.file)|!file.exists(graph.file2)|force==TRUE){
      cmd=paste("plotProfile", "-m", paste("'",outfile,"'",sep=""),
                "-out", paste("'",graph.file,"'",sep=""),
                "--plotType=lines",
                "--plotHeight", 8,
                "--plotWidth", 10,
                arg.ymin,
                arg.ymax,
                arg.kmeans,
                sep=" ")
      #system(cmd)
      cmd=paste("plotHeatmap", "-m", paste("'",outfile,"'",sep=""),
                "-out", paste("'",graph.file2,"'",sep=""),
                "--whatToShow", "\'heatmap and colorbar\'",
                "--heatmapWidth", 8, 
                "--colorList blue,white,red",
                "--colorNumber", N_color,
                arg.zmin,
                arg.zmax,
                arg.kmeans,
                arg.sort,
                "--outFileSortedRegions", paste("'",out.region,"'",sep=""),
                sep=" ")
      #system(cmd)
    }
    
    result[i]=outfile
    names(result)[i]=out.name
  }
  return(result)
}

read_deeptools_mat=function(sample, infile, outdir, flanking, bs, force=F){
  library("openxlsx")
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  bin.num=ceiling((flanking[2]-flanking[1]) / bs)
  position=seq(flanking[1], flanking[1]+(bin.num-1)*bs, by = bs)
  
  mat.file=paste(outdir,"/", sample, ".mat", sep="")
  #paste(outdir, "/", sub("\\.gz$", ".xlsx", basename(infile)), sep="")
  #avg.file=paste(outdir, "/", sub("\\.mat\\.gz$", ".xlsx", basename(infile)), sep="")
  result=c()
  for(i in 1:length(sample)){
    if(!file.exists(mat.file[i])|force){
    mat=read.table(gzfile(infile[i]), header = F, sep="\t", quote = "\"", skip=1, stringsAsFactors = F,check.names = F)
    #mat[is.na(mat)]=0
    rownames(mat)=mat[,4]
    mat=mat[, 7:ncol(mat)]
    colnames(mat)=position
    #avg=apply(mat, 2, mean)
    #avg.df=data.frame(aligned_position=names(avg),average_signal=avg, stringsAsFactors = F)
    #write.xlsx(mat, mat.file[i], col.names=T, row.names=T,keepNA=T,overwrite=T)
    write.table(mat,mat.file[i],quote = F,sep="\t",row.names = T,col.names = T)
    #write.xlsx(avg.df, avg.file[i])
    }
    temp.df=data.frame(Sample=sample[i], Mat=mat.file[i], stringsAsFactors = F)
    result=rbind(result, temp.df)
  }
  return(result)
}

cal_nuc_profile=function(mat.file, outdir, bs=5,fill=NA,sort=F){
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  library("reshape2")
  library("zoo")
  
  result=c()
  for(i in 1:length(mat.file)){
    cat(i,mat.file[i],"...\n")
    if(is.matrix(mat.file[[i]])){
      mat=mat.file[[i]]
      sample=names(mat.file)[[i]]
    }else if(is.data.frame(mat.file[[i]])){
      mat=mat.file[[i]]
      sample=names(mat.file)[[i]]
    }else if(file.exists(mat.file[[i]])){
      if(grepl("\\.xlsx$", mat.file[[i]])){
        mat=read.xlsx(mat.file[[i]],colNames = T, rowNames = T)
        sample=names(mat.file)[[i]]
        #sample=sub("(\\.mat)*\\.xlsx$", "", basename(mat.file[[i]]))
      }else{
        mat=read.table(mat.file[[i]], header = T, sep="\t", quote = "\"", stringsAsFactors = F,check.names = F)
        sample=names(mat.file)[[i]]
        #sample=sub("\\..*$", "", basename(mat.file[[i]]))
      }
    }else{
      stop("Need one matrix!\n")
    }
    if(sort==T){
      w1=apply(abs(mat), 1, function(x){mean(x,na.rm=T)})
      w2=apply(abs(mat), 1, function(x){sd(x,na.rm = T)})
      w2[w2==0]=NA
      w2[is.na(w2)]=min(w2,na.rm = T)
      w=w1-w2
      mat=mat[order(w, decreasing = T),]
    }
    #oi=c(20)
    #draw_mat_profiles(mat[oi,], paste(outdir,"/",sample,"_", rownames(mat)[oi][1], "_origins.pdf", sep=""),bs=bs, yMin=0)
    #mat=t(scale(t(abs(mat))))
    #mat=mat / weight
    #mat[!complete.cases(mat)]=0
    if(!is.na(fill)){mat[is.na(mat)]=fill} ##fill NA
    if(bs>0){
      mat.smooth=apply(mat, 1, function(x){rollapply(x, width=bs*2+1,FUN=function(y){mean(y,na.rm=T)},fill=0,partial=T)})
      mat.smooth=t(mat.smooth)
      colnames(mat.smooth)=colnames(mat)
      mat=mat.smooth
    }
    mat.smooth.file=paste(outdir,"/", sample, ".mat", sep="")
      #paste(outdir,"/", sample, ".xlsx", sep="")
    #write.xlsx(mat, mat.smooth.file, colNames=T, rowNames=T,keepNA=T)
    write.table(mat,mat.smooth.file,quote = F,sep="\t",row.names = T,col.names = T)
    
    avg=apply(mat, 2, function(x){mean(x,na.rm=T)})
    avg.df=data.frame(aligned_position=names(avg),average_signal=avg, stringsAsFactors = F)
    avg.file=paste(outdir, "/", sample, ".txt", sep="")
    write.table(avg.df, avg.file, quote = F, sep="\t", row.names = F, col.names = T)
    
    temp.df=data.frame(Sample=sample, Avg=avg.file, Mat=mat.smooth.file, stringsAsFactors = F)
    result=rbind(result, temp.df)
  }
  return(result)
}

##return the overlap between two strings
string_overlap=function(s1, s2){
  s1.char=unlist(strsplit(s1, ""))
  s2.char=unlist(strsplit(s2, ""))
  length(s1.char)=length(s2.char)=max(length(s1.char), length(s2.char))
  diff.index=which(s1.char != s2.char) ##split the string and find out the position with difference
  same=substr(s1, 1,diff.index[1]-1)
  return(paste(same, collapse = ""))
}


Bias_normalize=function(mat.file,sample.name,
                        sample.class,sample.condition,
                        target=c("H3K4me3_eSPAN","H3K56ac_eSPAN","PCNA_eSPAN","RNApolII_eSPAN","cac1_eSPAN","cac2_eSPAN","sir2_eSPAN"),
                        control=c("BrdU","BrdU","BrdU","BrdU","BrdU","BrdU","BrdU"),
                        adjbias.dir,sum.adjmat.file
                        ){
  if(!dir.exists(adjbias.dir)){dir.create(adjbias.dir,recursive = T)}
  if(is.null(sample.name)){sample.name=sub("\\..*$","",basename(mat.file))}
  if(length(sample.condition)!=length(mat.file)){sample.condition=rep('',length(mat.file));cat("empty sample condition\n")}
  
  my.treat=paste(sample.class,":",sample.condition,sep="")
  sample.pair=pair_sample(sample.name, my.treat, target=target, control=control)
  #write.xlsx(sample.pair, file=paste(adj.ll.dir,"/",project.name,"_profile_adjust_sample_pair.xlsx", sep=""))
  cat("Samples are paired as:\n")
  print(sample.pair)
  
  names(mat.file)=sample.name
  adj.df=adjust_signals(sample.pair, mat.file, adjbias.dir)
  ##merge info
  adj.df=merge(adj.df, sample.pair)
  write.xlsx(adj.df, sum.adjmat.file,overwrite=T)
  
  return(sample.pair)
}

##Pair up the samples based on treatments. Return a dataframe (Target_file, Control_file, Target_treat, Control_treat)
pair_sample=function(sample, treat, target, control){
  result=c()
  target=toupper(target)
  control=toupper(control)
  treat=toupper(treat)
  for (i in 1:length(treat)) {
    treat.prefix=sub("^([^:]+)\\:.*$", "\\1", treat[i])
    treat.suffix=sub("^[^:]+\\:(.*)$", "\\1", treat[i])
    t.flag=treat.prefix %in% target  ##if the sample is a target
    if(t.flag){
      con.prefix=control[target==treat.prefix]
      con.treat=paste(con.prefix,":", treat.suffix, sep="")
      con.index=which(treat==con.treat)
      if(length(con.index)!=0){ ##if have the paired control
        one.df=data.frame(Target=sample[i], Control=sample[con.index], Target_treat=treat[i], Control_treat=treat[con.index], stringsAsFactors = F)
        result=rbind(result, one.df)
      }
    }
  }
  return(result)
}

##MODE="SUB" or "FC"
adjust_signals=function(sample.pair, mat.file, outdir,MODE="SUB",fill=0,ps.count=0.5){
  if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
  t.index=sapply(sample.pair[, 1], function(x){which(names(mat.file)==x)})
  target.mat=mat.file[t.index]
  #names(target.mat)=bias.df$Sample[t.index]
  c.index=sapply(sample.pair[, 2], function(x){which(names(mat.file)==x)})
  control.mat=mat.file[c.index]
  adj.mat=mat_operation(target.mat, control.mat, outdir,MODE,fill,ps.count)
  adj.avg=sapply(adj.mat, function(x){
    if(grepl("\\.xlsx$", x)){
      mat=read.xlsx(x,colNames = T, rowNames = T)
    }else{
      mat=read.table(x, header = T, sep="\t", quote = "\"", stringsAsFactors = F, check.names = F)
    }
    sig=apply(mat, 2, function(x){mean(x,na.rm=T)})
    result=data.frame(aligned_position=colnames(mat), average_signal=sig, stringsAsFactors = F);
    outfile=sub("(\\.xlsx)*$", ".txt", x);
    names(outfile)=names(x);
    write.table(result, outfile, quote = F, sep="\t", row.names = F, col.names = T);
    return(outfile)})
  
  result=data.frame(Sample=names(target.mat), Mat=adj.mat, Avg=adj.avg, 
                    Target=names(target.mat), Control=names(control.mat),
                    stringsAsFactors = F)
  return(result)
}


##subtract two matrix file
mat_operation=function(target.mat, control.mat, outdir,MODE="SUB",fill=0,ps.count=0.5){
  if(!dir.exists(outdir)){dir.create(outdir)}
  adj.mat=c()
  for(i in 1:length(target.mat)){
    cat(i, "operating MODE=",MODE, target.mat[i], "with", control.mat[i], "\n")
    outfile=paste(outdir, "/", sub("(\\.[^.]*)*$", "_adjusted.mat", basename(target.mat[i])), sep="")
      #paste(outdir, "/", sub("(\\.xlsx)*$", "_adjusted.mat.xlsx", basename(target.mat[i])), sep="")
    if(file.exists(target.mat[i]) & file.exists(control.mat[i])){
      if(grepl("\\.xlsx$", target.mat[i])){
        t.df=read.xlsx(target.mat[i], colNames = T, rowNames = T)
      }else{
        t.df=read.table(target.mat[i], header = T, sep="\t", quote = "\"", stringsAsFactors = F, check.names = F)
      }
      if(grepl("\\.xlsx$", control.mat[i])){
        c.df=read.xlsx(control.mat[i], colNames = T, rowNames = T)
      }else{
        c.df=read.table(control.mat[i], header = T, sep="\t", quote = "\"", stringsAsFactors = F, check.names = F)
      }
      t.df[is.na(t.df)]=fill
      c.df[is.na(c.df)]=fill
      
      if(MODE=="SUB"){
        adj.df=as.matrix(t.df)-as.matrix(c.df)
      }else if(MODE=="FC"){
        adj.df=as.matrix(t.df)/(as.matrix(c.df)+ps.count)
      }
      #write.xlsx(adj.df, file=outfile, colNames=T, rowNames=T,keepNA=T)
      write.table(adj.df, file=outfile,quote = F,sep="\t",row.names = T,col.names = T)
      names(outfile)=names(target.mat[i])
      adj.mat=c(adj.mat, outfile)
    }
  }
  return(adj.mat)
}