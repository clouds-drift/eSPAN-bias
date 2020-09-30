required.pack = c("openxlsx","ggplot2","plyr","pheatmap","RColorBrewer")
for(p in required.pack){
  sig=suppressMessages(require(p,character.only = TRUE))
  if(!sig){stop(paste("package",p,"needed!"))}
}


##need improve:sig.class when file names are null; ppt and pdf
draw_profiles=function(sig.file, graph.file, sig.class=c(),group=c(),interval=F, 
                       xMin=NA,xMax=NA,yMin=NA,yMax=NA, xlab="Relative position", ylab="Log2 ratio of Watson/Crick",
                       x.breaks=NA,y.breaks=NA,width=NULL,height=NULL,fontsize=12,legend.ncol=2){
  library("ggplot2")
  
  if(!dir.exists(dirname(graph.file))){dir.create(dirname(graph.file), recursive = T)}
  if(length(sig.class)!=length(sig.file)){sig.class=sub("(\\.txt)|(\\.xlsx)$", "", basename(sig.file))}##valid for files
  if(is.vector(sig.file)){ ##input files
    if(is.null(names(sig.file))){names(sig.file)=sub("(\\.txt)|(\\.xlsx)$", "", basename(sig.file))}
    graph.list=tapply(sig.file, sig.class, function(x){x})
  }else{
    stop("Please input a file list!\n")
  }
  cat(length(sig.file), "files are drawn on", length(graph.list), "graphs:", names(graph.list), "\n")
  if(length(group)!=length(sig.file)){group=names(sig.file)}
  names(group)=names(sig.file)
  read_filelist=function(filelist){
    page.data=c()
    for (j in 1:length(filelist)) {
      one.file=filelist[j]
      if(is.null(names(filelist)[j])){
        names(one.file)=sub("(\\.txt)|(\\.xlsx)$", "", basename(one.file))
      }else{
        names(one.file)=names(filelist)[j]
      }
      if(is.list(one.file)){
        one.df=data.frame(position=names(one.file[[1]]),signal=one.file[[1]],stringsAsFactors = F)
      }else if(grepl("\\.xlsx$", one.file)){
        one.df=read.xlsx(one.file)
      }else{
        one.df=read.table(one.file, header = T, sep="\t", quote = "\"", stringsAsFactors = F)
      }
      ##group the profiles
      one.df=cbind(one.df,File=names(one.file),stringsAsFactors=F)
      page.data=rbind(page.data, one.df)
    }
    return(page.data)
  }
  plot_profile=function(sig.df,param){
    x.var=param[1];y.var=param[2]
    graph.title=param[3];fontsize=as.numeric(param[4])
    xmin=as.numeric(param[5]);xmax=as.numeric(param[6]);ymin=as.numeric(param[7]);ymax=as.numeric(param[8])
    xlab=param[9];ylab=param[10]
    x.breaks=param[11];y.breaks=param[12]
    interval=param[13]
    legend.ncol=as.numeric(param[14])
    
    sig.plot= ggplot(sig.df,aes(x=sig.df[,x.var], y=sig.df[,y.var])) +
      geom_line(aes(color=sig.df[,"Group"], linetype=sig.df[,"Group"]),size=1)+
      scale_color_manual(name="type", breaks=levels(sig.df[,"Group"]), 
                         #labels=sub("^([^_]+)_.+$", "\\1", unique(group.df$type)), 
                         values = c("black", "blue", "red","#009E73","orange", "yellow","purple", "grey"))+
      scale_linetype_manual(name="type", breaks=levels(sig.df[,"Group"]), values = rep("solid",8))+
      theme_bw()+
      theme(
        axis.text=element_text(size=fontsize, face="bold", colour = "black"),
        axis.title = element_text(size=16,face="bold",colour="black"),
        #axis.text.x = element_text(angle=90),
        axis.ticks = element_line(size = 1, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid=element_blank(), panel.border=element_blank(), panel.background = element_blank(), #remove grid
        axis.line.x = element_blank(),
        axis.line.y=element_line(colour="black",size=1),
        legend.position = "top",
        legend.text = element_text(size = fontsize)
      ) +
      guides(col=guide_legend(ncol=legend.ncol))+
      geom_vline(xintercept = c(0), linetype="dashed", color = "black", size=1)+
      geom_hline(yintercept = c(0), color="black", size=1)+
      scale_y_continuous(limits = c(ymin,ymax))+
      scale_x_continuous(limits = c(xmin,xmax),expand=expansion(mult=0))+
      labs(x=xlab, y=ylab,title=graph.title)
    if(interval=="TRUE"){
      sig.plot=sig.plot+geom_ribbon(aes(ymin=sig.df[,y.var]-sig.df[,"se"],ymax=sig.df[,y.var]+sig.df[,"se"],fill=sig.df[,"Group"]),alpha=0.5)+
        scale_fill_manual(name="type", breaks=levels(sig.df[,"Group"]), values = c("#0073C2FF","#EFC000FF","red","#009E73","orange","purple", "grey","black"))
    }
    if(!is.na(y.breaks)){
      y.breaks=as.numeric(y.breaks)
      sig.plot=sig.plot+scale_y_continuous(limits = c(ymin,ymax),breaks = seq(ymin*10e4,ymax*10e4,y.breaks*10e4)/10e4)
    }
    if(!is.na(x.breaks)){
      x.breaks=as.numeric(x.breaks)
      sig.plot=sig.plot+scale_x_continuous(limits = c(xmin,xmax),breaks = seq(xmin*10e4,xmax*10e4,x.breaks*10e4)/10e4)
    }
    
    return(sig.plot)
  }
  
  data.list=list()
  param.list=list()
  q=1
  for(i in 1:length(graph.list)){
    graph.title=names(graph.list)[i]
    file.list=graph.list[[i]]
    #names(file.list)=names(graph.list)[i]
    page.data=read_filelist(file.list)
    page.data$Group=group[page.data$File]
    if(is.na(xMin)){xmin=range(as.numeric(page.data[,1]))[1]}else{xmin=xMin}
    if(is.na(xMax)){xmax=range(as.numeric(page.data[,1]))[2]}else{xmax=xMax}
    if(is.na(yMin)){ymin=range(as.numeric(page.data[,2]))[1]}else{ymin=yMin}
    if(is.na(yMax)){ymax=range(as.numeric(page.data[,2]))[2]}else{ymax=yMax}
    page.data[,3]=factor(page.data[,3], levels = unique(page.data[,3]), ordered = T)
    page.data[,4]=factor(page.data[,4], levels = unique(page.data[,4]), ordered = T)
    #######seperate lines to groups of maxium 8 line
    line.type=unique(as.character(page.data[,"Group"]))
    split.group=split(line.type,rep(1:length(line.type),each=8,length.out=length(line.type)))##each plot maximum 8 lines
    page.list=lapply(1:length(split.group), function(x){
      return(page.data[as.character(page.data$Group) %in% split.group[[x]], ])})
    ################################################
    for(n in 1:length(page.list)){
      one.page=page.list[[n]]
      x.var=colnames(one.page)[1]
      y.var=colnames(one.page)[2]
      sig.df=summarySE(one.page, measurevar=y.var, groupvars=c(x.var,"Group"))
      data.list[[q]]=sig.df
      param.list[[q]]=c(x.var,y.var,graph.title,fontsize,xmin,xmax,ymin,ymax,xlab,ylab,x.breaks,y.breaks,interval,legend.ncol)
      q=q+1
    }
  }
  g2plot.list=lapply(1:length(data.list),
                     function(x){plot_profile(data.list[[x]],param.list[[x]])})
  if(grepl("\\.ppt$|\\.pptx$",graph.file)){
    graph2pptx(g2plot.list,graph.file,width=width,height=height)
  }else if(grepl("\\.pdf",graph.file)){
    pdf(graph.file,width=width,height=height)
    #par(mai=c(2,2,2,2))
    lapply(g2plot.list,print)
    dev.off()
  }
}

##list of ggplot objects to pptx
graph2pptx=function(g2plot.list,graph.file,width=NULL,height=NULL,picture=F){
  library("officer")
  library("rvg")
  if(!dir.exists(dirname(graph.file))){dir.create(dirname(graph.file),recursive = T)}
  if(grepl("\\.ppt$|\\.pptx$",graph.file)){
    graph.file=sub("\\.ppt$",".pptx",graph.file)
  }else{
    stop("Only support .pptx file!\n")
  }
  
  my.pres=read_pptx()
  for(i in  1:length(g2plot.list)){
    my.pres=add_slide(my.pres)
    if(picture){##export as picture of 300dpi
      my.pres=ph_with(my.pres,value=g2plot.list[[i]],location = ph_location(width=width,height=height),res=300)
    }else{
      my.pres=ph_with(my.pres,value=dml(ggobj=g2plot.list[[i]]),location = ph_location(width=width,height=height))
      #my.pres=ph_with(my.pres,value=g2plot.list[[i]],location = ph_location(width=width,height=height)) ##not editable
    }
  }
  print(my.pres, target = graph.file)
}

##Description:
##  Draw heatmap for aligned reads density.
##Args:
##  mat.file, files of reads density matrix.
##  graph.dir, output dir
##  my.clust, a list of clustering object
draw_heatmap=function(mat.file, graph.file, 
                      zMin=NA, zMax=NA,zRange=c(0,1),smooth=NA,fill=NA, my.color=NA, ncolor=100, 
                      sort=FALSE,show_rownum=1,show_colnum=1,width=NA,height=NA,...){
  library("openxlsx")
  library("pheatmap")
  library("RColorBrewer")
  
  obj.list=c()
  for(i in 1:length(mat.file)){
    cat(i,"loading data of ")
    if(is.matrix(mat.file[[i]])){
      #mat.list[[i]]=mat.file[[i]]
      plot.mat=mat.file[[i]]
      #names(mat.list)[i]=names(mat.file)[[i]]
      plot.name=names(mat.file)[[i]]
    }else if(is.data.frame(mat.file[[i]])){
      plot.mat=mat.file[[i]]
      plot.name=names(mat.file)[[i]]
    }else if(is.vector(mat.file[[i]]) & length(mat.file[[i]])==1){
      if(file.exists(mat.file[[i]])){
        if(grepl("\\.xlsx$", mat.file[[i]])){
          plot.mat=read.xlsx(mat.file[[i]],colNames = T, rowNames = T)
          plot.name=sub("(\\.mat)*\\.xlsx$", "", basename(mat.file[[i]]))
        }else{
          plot.mat=read.table(mat.file[[i]], header = T, sep="\t", quote = "\"", stringsAsFactors = F)
          plot.name=sub("\\..*$", "", basename(mat.file[[i]]))
        }
      }else{
        stop(paste(i,mat.file[[i]],"not exists!\n"))
      }
    }else{
      stop(paste(i,names(mat.file)[i],"is not a matrix!\n"))
    }
    cat(plot.name,"\n")
    
    ##scale in row and remove NAN rows
    #plot.mat=t(scale(t(plot.mat))) ##scale in row
    remove.ind=which(apply(plot.mat, 1, function(x){all(is.nan(x))}))
    if(length(remove.ind)>0){
      cat(paste(rownames(plot.mat)[remove.ind], collapse = ","), "removed because of NaN\n")
      plot.mat=plot.mat[-c(remove.ind), ]
    }
    ##set zMin and zMax
    if(is.na(zMin) | is.na(zMax)){
      #row.min=apply(plot.mat, 1, function(x){quantile(x, zRange[1], na.rm=TRUE)})
      #row.max=apply(plot.mat, 1, function(x){quantile(x, zRange[2], na.rm=TRUE)})
      #my.range=c(quantile(row.min,0.1,na.rm=T), quantile(row.max,0.9,na.rm=T))
      vmin=max(min(unlist(plot.mat),na.rm=T), mean(unlist(plot.mat),na.rm=T)-2*sd(unlist(plot.mat),na.rm=T))
      vmax=min(max(unlist(plot.mat),na.rm=T), mean(unlist(plot.mat),na.rm=T)+2*sd(unlist(plot.mat),na.rm=T))
      my.range=c(vmin,vmax)
    }
    if(is.na(zMin)){zMin=my.range[1]}
    if(is.na(zMax)){zMax=my.range[2]}
    ##sort the mat by row average
    if(!is.na(smooth)){ plot.mat=apply(plot.mat, 1, function(x){rollapply(x, width=2*smooth, function(y){mean(y,na.rm=T)})})}
    if(!is.na(fill)){plot.mat[is.na(plot.mat)]=fill} ##fill NA
    if(sort==T){
      cat("sorting the matrix...\n")
      w1=apply(abs(plot.mat), 1, function(x){mean(x,na.rm=T)})
      w2=apply(abs(plot.mat), 1, function(x){sd(x,na.rm = T)})
      w2[w2==0]=NA
      w2[is.na(w2)]=min(w2,na.rm = T)
      w=w1
      #w1-w2
      #w1*(w2)
      plot.mat=plot.mat[order(w, decreasing = T),]
    }
    plot.mat[plot.mat > zMax]=zMax
    plot.mat[plot.mat < zMin]=zMin
    ##set color. my.color=c("green", "black", "red")
    if(is.na(my.color[1])){
      my.color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(ncolor)
    }else{
      if(all(grepl("^#",my.color))){
        ##keep color
      }else{
        my.color=colorRampPalette(my.color)(ncolor)
      }
    }
    if((show_rownum >=0)&(show_rownum <= 1)){
      show_rownum=nrow(plot.mat)*show_rownum
    }
    if(show_rownum>1){
      cut=ceiling(nrow(plot.mat)/show_rownum)
      show.ind=seq(cut,nrow(plot.mat),by=cut)
      tmp.name=rownames(plot.mat)[show.ind]
      plot.mat=as.matrix(plot.mat)
      rownames(plot.mat)=rep("",nrow(plot.mat))
      rownames(plot.mat)[show.ind]=tmp.name
    }
    if((show_colnum >=0)&(show_colnum <= 1)){
      show_colnum=ncol(plot.mat)*show_colnum
    }
    if(show_colnum>1){
      cut=ceiling(ncol(plot.mat)/show_colnum)
      show.ind=seq(cut,ncol(plot.mat),by=cut)
      tmp.name=colnames(plot.mat)[show.ind]
      plot.mat=as.matrix(plot.mat)
      colnames(plot.mat)=rep("",ncol(plot.mat))
      colnames(plot.mat)[show.ind]=tmp.name
    }
    my.break=seq(zMin, zMax, length.out = ncolor+1)
    #cat(i,"draw",plot.name,"...\n")
      p.obj=pheatmap(plot.mat, 
                     color=my.color,
                     breaks=my.break,
                     border_color = "NA",
                     main = paste(nrow(plot.mat),"rows for",plot.name, sep=" "),
                     width=width,height=height,
                     ...)
      #gg=as.ggplot(p.obj)
      #gg.list[[i]]=gg
      #names(gg.list)[i]=plot.name
      #if(i==1){graph2ppt(x=p.ggplot,file=graph.file,append=F,vector.graphic=F,width=width,height=height)}else{graph2ppt(x=p.ggplot,file=graph.file,append=T,vector.graphic=F,width=width,height=height)}
      obj.list[[i]]=p.obj
      names(obj.list)[i]=plot.name
  }
  
  ##export to file
  if(grepl("\\.ppt",basename(graph.file))){
    if(!dir.exists(dirname(graph.file))){dir.create(dirname(graph.file),recursive = T)}
    library('ggplotify')
    gg.list=lapply(obj.list, as.ggplot)
    names(gg.list)=names(obj.list)
    graph2pptx(gg.list,graph.file, width=width, height=height, picture=T)
  }else if(grepl("\\.pdf",basename(graph.file))){
    if(!dir.exists(dirname(graph.file))){dir.create(dirname(graph.file),recursive = T)}
    pdf(graph.file,width=width,height=height,onefile = T,bg="white")
    #par(mai=c(10,10,10,10))
    lapply(obj.list,function(x){plot(x[[4]])}) ##plot gtable object
    dev.off()
  }else{
    if(!dir.exists(graph.file)){dir.create(graph.file,recursive = T)}
    for(i in 1:length(obj.list)){
      jpeg(filename =file.path(graph.file,paste(names(obj.list)[i],".jpg",sep="")),res=300)
      print(obj.list[[i]])
      dev.off()
    }
  }
    
  result=list(obj.list)
  return(result)
}


summarySE=function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=0.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac=ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac=rename(datac, c("mean" = measurevar))
  
  datac$se=datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is 0.95, use 0.975 (above/below), and use df=N-1
  ciMult = qt(conf.interval/2 + .5, datac$N-1)
  datac$ci = datac$se * ciMult
  
  return(datac)
}