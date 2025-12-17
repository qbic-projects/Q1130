##
rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path
setwd(path)
#setwd("/home/heumos/QBiC/Q1130/Results_Proteus_LFQ_CBOLT/Results_Proteus_LFQ_CBOLT/R/")

#load packages and functions, installation of some of them is quite time consuming
#e.g. devtools, pca3d and proteus
#To install proteus you might need rmarkdown and devtools if not yet installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("rmarkdown")){
  BiocManager::install("rmarkdown")
  library("rmarkdown")
}
if (!require("devtools")){
  BiocManager::install("devtools")
  library("devtools")
}

if (!require("limma")){
  BiocManager::install("limma")
  library("limma")
}

if (!require("proteus")){
  #In order to run examples or vignette code, 
  #additional packages with example data need to be installed:
  devtools::install_github("bartongroup/proteusLabelFree")
  devtools::install_github("bartongroup/proteusTMT")
  devtools::install_github("bartongroup/proteusSILAC")
  devtools::install_github("bartongroup/proteus", 
                           build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)
  library("proteus")
  library("proteusLabelFree")
  library("proteusSILAC")
  library("proteusTMT")
}

if (!require("dendextend")){
  BiocManager::install("dendextend")
  library("dendextend")
}
if (!require("ggplot2")){
  BiocManager::install("ggplot2")
  library("ggplot2")
}
if (!require("reshape")){
  BiocManager::install("reshape")
  library("reshape")
}

if (!require("rgl")){
  BiocManager::install("rgl")
  library("rgl")
}
if (!require("biomaRt")){
  BiocManager::install("biomaRt")
  library("biomaRt")
}
if (!require("gplots")){
  BiocManager::install("gplots")
  library("gplots")
}
if (!require("genefilter")){
  BiocManager::install("genefilter")
  library("genefilter")
}
if (!require("RColorBrewer")){
  BiocManager::install("RColorBrewer")
  library("RColorBrewer")
}
if (!require("EnhancedVolcano")){
  BiocManager::install("EnhancedVolcano")
  library("EnhancedVolcano")
}
if (!require("xtable")){
  BiocManager::install("xtable")
  library("xtable")
}
if (!require("UpSetR")){
  BiocManager::install("UpSetR")
  library("UpSetR")
}
if (!require("dbplyr")){
  BiocManager::install("dbplyr")
  library("dbplyr")
}


#begin main analysis
mt <- data.frame(norms=c("p1n <- normalizeData(p1,norm.fun = normalizeMedian)",
                         "p1n <- normalizeData(p1,norm.fun = normalizeQuantiles)"),
                 filt = c(0,0,2,2,6,6))

mt$norms <- as.character(mt$norms)
mt$filt <- as.character(mt$filt)
k <- dim(mt)[1]
k=2  #run on quantile data only


#loop was removed
  ki=k
  print(ki)
  #prepare folders
  wd = "Results_Proteus_LFQ_CBOLT"
  tmp = "tmp"
  system(paste("rm -rf ",wd,sep = ""))
  dir.create(wd)
  dir.create(paste(wd,tmp,sep = "/"))
  file.copy("CBOLT_sample_preparations.txt",to = wd)
  file.copy("proteinGroups_CBOLT.txt",to = wd)
  file.copy("function.R",to = wd)
  file.copy("Q1130_Metaproteomics_SH.R", to = wd)
  
  source(paste(wd,"/function.R",sep = ""))
  #
  # Create experimentalDesign
  m <- read.delim(paste(wd,"CBOLT_sample_preparations.txt",sep = "/"),header = T)
  m$sample = m$QBiC.Code
  m$condition = m$Condition..culture
  m$condition
  m$culture = m$Condition..culture
  m$Secondary.Name <- gsub(" ","",m$Secondary.Name)
  #create Identifier readable by Proteus
  measure.cols.LFQ <- setNames(paste("LFQ intensity", m$QBiC.Code, sep = " "), m$QBiC.Code)
  measure.cols.LFQ
  
  str(proteinColumns)
  
  p1 <- readProteinGroups(paste(wd,"proteinGroups_CBOLT.txt", sep = "/"),
                          m, measure.cols = measure.cols.LFQ, data.cols = proteinColumns)     
  #adapt colnames and sample attribute in metadata
  colnames(p1$tab) = paste(p1$metadata$Secondary.Name,"-",p1$metadata$Lab.ID,"_",colnames(p1$tab),sep = "")
  colnames(p1$tab)
  p1$metadata$sample = colnames(p1$tab) 

  #loop through mt:
  mt_sub <- mt[ki,]
  #define some vars:
  id_for_plots <- gsub(".* = normalize","",mt_sub[,1])
  id_for_plots <- gsub("\\)","",id_for_plots)
  filt <- paste("log2_expression > 0 in ",mt_sub[,2]," samples",sep = "")
  
  ### (3) Protein intensities in this file are summed over protein groups, so skip peptide and protein aggregation steps.
  # Normalize intensities (use standard = median)
  eval(parse(text=mt_sub[,1]))
  if (class(p1n)[2] != "proteusData") 
    stop("p1n object has to be a proteusData type")
  #
  png(paste(wd,tmp,'plot_SampleDistributions_NotNormalized.png',sep = "/"), width = 5*300, height = 5*300, res = 300, pointsize = 8)
  print(plotSampleDistributions(p1, title="Not normalized", fill="condition", method="violin"))
  dev.off()
  #
  png(paste(wd,tmp,'plot_SampleDistributions_Normalized.png',sep = "/"), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
  print(plotSampleDistributions(p1n, title=paste(id_for_plots," normalization",sep = ""), fill="condition", method="violin"))
  dev.off()
  #
  png(paste(wd,tmp,'plot_Samples_Mean-variance relationship.png', sep = "/"),width = 5*300, height = 5*300, res = 300, pointsize = 8) 
  print(plotMV(p1n, with.loess=TRUE))             
  dev.off()
  #
  png(paste(wd,tmp,'plot_Samples_ClusteringDendrogram.png',  sep = "/"),width = 5*300, height = 5*300, res = 300, pointsize = 8)
  print(plotClustering(p1n))                            
  dev.off()
  
  #own dendrogram
  #clustering
  #for arrays, problem: arrays are not row names but at column position, thus transpose is needed
  corr.mat <- cor(p1n$tab, use="complete.obs")
  dis <- as.dist(1 - corr.mat)  # dissimilarity matrix
  hc <- hclust(dis)
  dend <- as.dendrogram(hc)
  #remember groups, assign a new color code
  pdf(paste(wd,tmp,"Samples_ClusteringDendrogram_own.pdf",sep = "/"),height = 12,width = 10)
  par(oma=c(2,2,2,2))
  dend %>% set("labels_cex",0.4) %>% plot()
  dev.off()
  

  #own ggplot
  x = p1n
  dat <- x$tab
  condition = x$metadata$condition
  condition = relevel(condition, ref="Biculture")
  culture = x$metadata$culture
  sample = colnames(dat)
  #sample = gsub("-.*","",sample)
  tx=which(apply(t(na.omit(log10(dat))), 2, var)==0)
  tx=row.names(dat) %in% names(tx)
  tx
  dat=dat[!tx,]
  pca <- prcomp(t(na.omit(log10(dat))), scale. = TRUE, center=TRUE)
  sdat <- summary(pca)$importance
  sdat <- as.data.frame(sdat)
  sdat <- sdat[2,]
  sdat = melt(sdat)
  
  #scree plot
  plot <- ggplot(data=sdat, aes(x=variable, y=value)) + geom_bar(stat="identity",color="black") +
    ylab('Proportion of Variance') +
    xlab('') +
    ggtitle('') +
    ## theme_bw() +  ##for white background
    theme(legend.title = element_blank()) +
    theme(text = element_text(size=8))
  ggsave(filename = paste(wd,'plot_Screeplot.png',sep = "/"),plot = plot)
  rm(plot)
  
  levels = data.frame(t1=c(1,2),t2=c(1,3),t3=c(2,3))
  levels
  
  for (i in 1:dim(levels)[2]) {
    ##i=1
    fn=levels[,i]
    print(fn)
    p <- data.frame(
      x = pca$x[, fn[1]],
      y = pca$x[, fn[2]],
      condition = condition,
      culture = culture,
      sample = sample)
    var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
    x1 = paste("PCA",fn[1],sep = "")
    y1 = paste("PCA",fn[2],sep = "")
    xx1 <- sprintf("(%5.1f%%)", var.perc[fn[1]])
    yy1 <- sprintf("(%5.1f%%)", var.perc[fn[2]])
    x1 = paste(x1,xx1)
    y1 = paste(y1,yy1)
    rm(xx1,yy1)
    point.size=1.2
    text.size=6
    label.size=2
    legend.text.size=7
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #cbPalette <- brewer.pal(n = 8, name = "RdBu")
    palette=cbPalette
    with.legend=TRUE
    g <- ggplot(p, aes_(x=~x, y=~y, label=~sample)) +
      theme_bw() +
      coord_cartesian(expand=TRUE) +
      geom_point(aes_(color=~condition,shape=culture), size=point.size) +
      ggrepel::geom_text_repel(size=label.size) +
      scale_color_manual(values=palette) +
      theme(text = element_text(size=text.size)) +
      labs(x=x1, y=y1)
    if(!with.legend) g <- g + theme(legend.position="none")
    ggsave(filename = paste(wd,"/plot_PCA_normalized_data_own_","PCA",fn[1],"_vs_","PCA",fn[2],".png",sep = ""),plot = g)
  }
  
  #set protein filter
  f1 <- kOverA(as.numeric(mt_sub[,2]),0,na.rm = FALSE)
  #f1 <- kOverA(4,0)
  ffun <- filterfun(f1)
  data <- p1n$tab
  data = log2(data)
  stopifnot(identical(colnames(data),colnames(p1$tab)))
  test = data ; data <- 0 #do this because I am not sure whether kOverA treats NA's as I expect
  test[is.na(test)] <- 0
  which <- genefilter(test, ffun)
  print(table(which))
  test <- test[which,]
  #set values back to NA because NA does not mean the same as 0 in Px
  test[test == 0] <- NA
  data <- test
  rm(test)
  #
  #DE analysis Limma way
  TS <- m$condition
  TS <- factor(TS)
  TS
  
  design <- model.matrix(~0+TS)
  colnames(design) <- levels(TS)
  #colnames(design)
  #use the self-written function to do that job
  x <- c(paste(colnames(design)[2],"-",colnames(design)[1],sep = ""))
  cont.matrix <- makeContrasts(contrasts=x,levels=design)
  colnames(cont.matrix) <- gsub("\\-","_vs_",colnames(cont.matrix))
  
  #
  fit1 <- lmFit(data, design) #Warning message:Partial NA coefficients for 131 probe(s) 
  fit1 <- contrasts.fit(fit1, cont.matrix)
  fit1 <- eBayes(fit1)
  nx = colnames(cont.matrix)
  nx
  #
  bg=0
  stats=0
  for (f in 1:length(nx)) {
    ##f=3
    print(paste("calculate contrast: ",nx[f],sep = ""))
    tmp = topTable(fit1,coef = f,number=Inf,adjust.method="BH")
    tmp = tmp[order(row.names(tmp)),]
    alpha = table(tmp$P.Value < 0.05)[2]
    beta = table(tmp$adj.P.Val < 0.05)[2]
    alpha = data.frame(id=nx[f],P.Value_0.05=alpha,adj.P.Val_0.05=beta)  
    alpha[is.na(alpha)] <- 0
    names(tmp) = paste(names(tmp),colnames(cont.matrix)[f],sep="_")
    tmp = tmp[,c(2,1,3:6)]
    tmp1 = tmp
    tmp1$protein = row.names(tmp1)
    tmp1$protein = sapply(strsplit(as.character(tmp1$protein),";"),"[[",1)
    tmp1 <- subset(tmp1, tmp1[,4] < 0.05)
    write.table(tmp1, file = paste(wd,"/DE_table_",nx[f],".txt",sep = ""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", 
                na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
    tmp$protein = row.names(tmp)
    bg = cbind(bg,tmp)
    stats = rbind(stats,alpha)
    rm(tmp,tmp1,alpha,beta)
  }
  stats <- stats[-1,]
  row.names(stats) <- NULL
  idx <- duplicated(t(bg))
  bg <- bg[, !idx]
  bg <- bg[,-1]
  #
  stopifnot(identical(row.names(bg),bg$protein))
  # add columns with number of good replicates
  ngood <- reshape2::dcast(p1n$stats, id ~ condition, fun.aggregate=sum, value.var="ngood")
  #ngood <- ngood[, c("id", condition)]
  names(ngood)[2:ncol(ngood)] <- paste0("ngood_", names(ngood)[2:ncol(ngood)])
  names(ngood)[1] = "protein"
  bg <- merge(bg, ngood, by="protein")
  #names(bg)[1] = "protein"
  ### (4) Protein annotations
  summary(p1n$proteins)
  #
  # Extract UniProt IDs_ edited version 
  luni3 <- lapply(as.character(bg$protein), function(prot) 
  {
    if(grepl(";", prot))                                         # ;
    {
      uniprot <- unlist(strsplit(prot, ";", fixed=TRUE))[1] # here only the UniProt IDs are listed, we just picked the first one every time!
      c(prot, uniprot)
    } else {
      c(prot, prot)
    }
  }
  )
  ids3 <- as.data.frame(do.call(rbind, luni3))                         
  names(ids3) <- c("protein", "uniprot") 
  head(ids3)
  #
  # Gene names and protein names from UniProt, ignore this because of use with biomart below
  #annot <- fetchFromUniProt(ids3$uniprot, verbose=TRUE,
  #                           batchsize = 100,columns = c("genes", "protein names"))
  #
  # #NCBI info
  # load("saureus_ensembl.Rdata")
  # attributes <- listAttributes(ensembl)
  # #ensembl = useDataset(mart = ensembl, "hsapiens_gene_ensembl")
  # ens <- getBM(attributes=c("uniprot_gn_id", "external_gene_name","entrezgene_id", "description", "chromosome_name"), 
  #              filters = "uniprot_gn_id", 
  #              values = ids3$uniprot, mart = ensembl, uniqueRows = TRUE,useCache = F)    
  # names(ens)[c(1,2)] = c("uniprot","gene_id")
  # names(ens)
  # annot <- merge(ids3,ens,by = "uniprot",all.x = TRUE)
  # annot <- annot[!duplicated(annot$uniprot),]
  # #
  # #final merge
  # bg <- merge(bg,annot,by="protein")
  # bg <- bg[!duplicated(bg$protein),]
  # rm(ids3,annot)
  #
  #bg <- bg[order(bg$adj.P.Val),]
  #
  
  #final merge
  bg <- merge(bg,ids3,by="protein")
  
  #loop through contrasts for plots
  row.names(bg) <- bg$protein
  rm(f)
  setwd(path)
  
  for (f in nx) {
    ##f=nx[1]
    print(f)
    sub <- bg[,grepl(paste(f,"ngood",sep = "|"),names(bg))]
    names(sub) <- gsub(paste("_",f,sep = ""),"",names(sub))
    png(paste(wd,'/tmp/plot_volcano-plot_',f,'.png',sep = ""), width = 5*300, height = 5*300, res = 300, pointsize = 8)
    print(plotVolcano(sub))
    dev.off()
    png(paste(wd,'/plot_enhanced_volcano-plot_1',f,'.png',sep = ""), width = 9*300, height = 5*300, res = 300, pointsize = 8)
    print(EnhancedVolcano(sub,
                  lab = row.names(sub),
                  x = 'logFC',
                  y = 'P.Value',
                  ylab = bquote(~-Log[10]~italic(p-value)),
                  xlab = bquote(~Log[2]~ "FC"),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  labSize = 6.0,
                  title = paste("Proteins: ",gsub("_"," ",f),sep=""),
                  ylim = c(-0.1, 10.0),
                  colAlpha = 2,
                  drawConnectors = F,
                  widthConnectors = 0.75,
                  subtitle = "Differential protein expression",
                  caption = paste('Total = ', nrow(sub), 'proteins. '),
                  legendLabels = c('NS', expression(Log[2]~FC),
                                   "P.Val < 0.05", "P.Val < 0.05 & logFC >1|<-1"),
                  selectLab = c("TEST") # this removes all the labels as we don't have a "TEST" gene name
  ))
  dev.off()
  
  png(paste(wd,'/plot_enhanced_volcano-plot_2',f,'.png',sep = ""), width = 9*300, height = 5*300, res = 300, pointsize = 8)
  print(EnhancedVolcano(sub,
                        lab = row.names(sub),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        ylab = bquote(~-Log[10]~italic(adj.P.Val)),
                        xlab = bquote(~Log[2]~ "FC"),
                        pCutoff = 0.05,
                        FCcutoff = 1,
                        labSize = 6.0,
                        title = paste("Proteins: ",gsub("_"," ",f),sep=""),
                        ylim = c(-0.1, 10.0),
                        colAlpha = 2,
                        drawConnectors = F,
                        widthConnectors = 0.75,
                        subtitle = "Differential protein expression",
                        caption = paste('Total = ', nrow(sub), 'proteins. '),
                        legendLabels = c('NS', expression(Log[2]~FC),
                                         "adj.P.Val < 0.05", "adj.P.Val < 0.05 & logFC >1|<-1"),
                        selectLab = c("TEST") # this removes all the labels as we don't have a "TEST" gene name
  ))
  dev.off()
  
  png(paste(wd,'/plot_p-value-distribution-plot_',f,'.png',sep = ""), width = 5*300, height = 5*300, res = 300, pointsize = 8)
  print(plotPdist(sub))
  dev.off()
  #
  #
  # Plots for top DE proteins
  setwd(wd)
  dir.create("plots_candidate_genes")
  genes <- sub[order(sub$adj.P.Val),]
  
  #as we have DE genes based on p.adj let's get them all plotted
  genes <- subset(genes,adj.P.Val < 0.05 )
  #### THIS ONLY GENERATES PLOTS FOR THE FIRST 5 GENES AT THE MOMENT! ####
  for (i in row.names(genes)[1:530])
  {
    #i=row.names(genes)[5]
    print(i)
    gene=subset(bg, row.names(bg) == i)
    uniprot <- as.character(gene$uniprot)
    #gene_id <- as.character(gene$uniprot)
    padj = which(names(gene) == paste("adj.P.Val_",f,sep = ""))
    padj <- round(gene[,padj],digits = 4)
    logFC = which(names(gene) == paste("logFC_",f,sep = ""))
    logFC <- round(gene[,logFC],digits = 4)
    plot_data <- subset(data, row.names(data) == i)
    plot_data <- as.data.frame(plot_data)
    plot_data <- melt(plot_data)      #Error in data.frame(indices, value = values) : arguments imply differing number of rows: 16, 0
    #### CUSTOM CODE TO DO THE SPLITTING SO THAT WE ACTUALLY OBTAIN THE QBIC IDS ####
    for (j in seq(1:length(plot_data$variable))) {
      splitter <- strsplit(as.character(plot_data$variable[j]),"_")
      len_split <- length(splitter[[1]])
      qbic_code <- splitter[[1]][len_split]
      plot_data$QBiC.Code[j] <- qbic_code
    }
    #plot_data$QBiC.Code <- gsub("-.*","",plot_data$variable)
    #plot_data$gene_id <- gene_id
    plot_data$uniprot <- uniprot
    m_set <- m[,grepl("QBiC.Code|Condition...",names(m))]
    plot_data <- merge(plot_data,m_set,by = "QBiC.Code")
    ##plot_data$uniprot <- ifelse(is.na(plot_data$uniprot),"",plot_data$uniprot)
    plot <- qplot(x=Condition..culture,y=value,data=plot_data, geom=c("boxplot", "jitter"),
                  fill=Condition..culture,main=paste("Uniprot:",uniprot,";logFC:",logFC,";padj:",padj,sep=""),
                  xlab="",ylab="log2 transformed expression") +
      facet_grid(~Condition..culture, scales="free",space = "free")
      theme(text = element_text(size=8))
    ggsave(filename=paste("plots_candidate_genes/",uniprot,".png",sep=""),
           width=7, height=7,plot=plot)
  }
  #write.table(genes, file = paste("plots_candidate_genes/","diff_expr_final_table.txt",sep = "/"), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
  system(paste("mv plots_candidate_genes ","plots_candidate_genes_",f,sep = ""))
  setwd(path)
  }
  
  setwd(path)
  # #
  
  
  #Upset Plot
  sub1 <- sub[,grepl("^ngood_",names(sub))]
  sub1 <- ifelse(sub1 == 0, 0,1)
  sub1 <- as.data.frame(sub1)
  upset(sub1,order.by = "freq",nsets = ncol(sub1),sets.bar.color = "#56B4E9")
  #### SOMEHOW HERE THE VARIABLE WD IS WRONG, THIS CURRENTLY FIXES IT FOR ME, PLEASE TRY OUT FOR YOU ####
  #png(paste(wd,'/plot_upset.png',sep = ""), width = 5*300, height = 5*300, res = 300, pointsize = 8)
  png(paste(path,'plot_upset.png',sep = ""), width = 5*300, height = 5*300, res = 300, pointsize = 8) 
  #print(upset(sub1,order.by = "freq",nsets = ncol(sub1),sets.bar.color = "#56B4E9"))
  upset(sub1,order.by = "freq",nsets = ncol(sub1),sets.bar.color = "#56B4E9")
  dev.off()
  
  #
  #
  #heatmaps, first adapt to better row.names
  #########################################################
  ### C) Customizing and plotting the heat map
  #########################################################
  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
  my_palette1 <- colorRampPalette(c("red", "#A6A3A1", "blue"))(n = 299)
  mypalette2 <- brewer.pal(11,"RdYlBu")
  my_palette2 <- colorRampPalette(mypalette2)(n=25)
  palettes <- list(my_palette,my_palette2)


  #heatmaps, all genes
  data1 <- data
  # datn <- sub[,c(1,2)]
  # row.names(datn) <- datn[,1]
  # datn[,1]=NULL
  # data1 <- merge(data1,datn,by=0)
  # row.names(data1) <- data1$uniprot
  # data1$uniprot <- NULL
  # data1[,1]=NULL
  # #
  rnames <- as.character(row.names(data1))                           # assign labels in column 1 to "rnames"
  mat_data <- data.matrix(data1[,1:ncol(data1)])  # transform column 2-5 into a matrix
  rownames(mat_data) <- rnames                  # assign row names
  colnames(mat_data) <- gsub("_Q.*","",colnames(mat_data))
  l <- dim(mat_data)[1]
  mat_data[is.na(mat_data)] <- 0
  #
  # creates a 5 x 5 inch image
  for (p in 1:length(palettes)) {
    fn = palettes[[p]]
    #### SAME AS BEFORE THE WD DOES NOT POINT TO THE CORRECT PATH WHERE THE FILE SHOULD BE SAVED ####
    # SO I REPLACED "WD" WITH "PATH"
    #png(paste(wd,"/plot_heatmaps_all_proteins_",p,".png",sep = ""),    # create PNG for the heat map
    png(paste(path,"/plot_heatmaps_all_proteins_",p,".png",sep = ""),    # create PNG for the heat map        
       width = 5*300,        # 5 x 300 pixels
        height = 5*300,
      res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size

    #par(cex.main=0.6)
    #par(oma=c(3,3,3,3))

    heatmap.2(mat_data,
              #cellnote = mat_data,  # same data set for cell labels
              main = paste(as.character(l)," , all proteins"), # heat map title", # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="density",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,5),     # widens margins around plot
              col=fn,       # use on color palette defined earlier
              #breaks=col_breaks,    # enable color transition at specified limits
              dendrogram="both",     # only draw a row dendrogram
              key.ylab="Density",
              key.xlab=paste("scaled",expression(log[2]("expression")),sep = " "),          # removes stupid color key x labelling
              keysize = 1,scale="c",key=TRUE, symkey=FALSE,
              #Colv="NA",            # turn off column clustering
              #sepwidth=c(0.05, 0.05),
              #sepcolor="white",
              #colsep=1:ncol(mat_data),
              #rowsep=1:nrow(M31_exprdata)),
              labRow = "",
              #adjCol = c(NA,0.5),
              srtCol=45,
              offsetCol = 1,
              #labCol=labCol,
              cexRow=0.5,
              cexCol=0.5)
    dev.off()               # close the PNG device
  }
  setwd(path)
  
  #write out final tables
  write.table(data, file = paste(wd,"/log2_transformed_expression.txt",sep = ""), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", 
              na = "NA", dec = ".", row.names = T,  col.names = T, qmethod = c("escape", "double"))
  write.table(bg, file = paste(wd,"final_table.txt",sep = "/"), append = FALSE, quote = FALSE, sep = "\t",eol = "\n", na = "NA", dec = ".", row.names = T,  col.names = NA, qmethod = c("escape", "double"))
  #
  #
  ##} 
  setwd(path)
  
  dir.create(paste(wd,"tables",sep = "/"))
  dir.create(paste(wd,"plots",sep = "/"))
  dir.create(paste(wd,"R",sep = "/"))
  
  setwd(wd)
  system("mv plots_candidate_genes* plots/")
  system("mv plot_* plots/")
  
  #end of script
  ####-------------save Sessioninfo
  fn <- paste("R_sessionInfo_",format(Sys.Date(), "%d_%m_%Y"),".txt",sep="")
  sink(fn)
  print(sessionInfo())
  sink()
  ####---------END----save Sessioninfo 
  system("mv R_* R/")   #mv: rename R_sessionInfo_13_10_2021.txt to R/: No such file or directory
  system("mv *.R R/")   #usage: mv [-f | -i | -n] [-v] source target
                        #mv [-f | -i | -n] [-v] source ... directory
  system("mv *.Rdata R/")  #mv: rename ensembl.Rdata to R/: No such file or directory
  system("mv *.txt tables/")
  #
  #setwd(path)
  # #create final folder   id
  fn1 <- paste("readme.txt",sep="")
  sink(fn1)
  fi <- paste("Processing of data: ",id_for_plots," normalized data + Filter based on ",filt,sep = "")
  print(fi)
  print("####")
  print("Contrasts analyzed: (diff=interaction term)")
  print(stats)
  print("####")
  sink()
  #Create an xtable for the object that you want to export to html: 
  stats <- xtable(stats,digits = 0)
  print.xtable(stats, type="html", file="stats.html")
  #system(paste("mv ",wd," ",wd,"_final_",ki,sep = ""))
  setwd(path)


setwd(path)
#system(paste("mv ",wd," ",wd,"_",tests,sep = ""))

##--end

##Improved Volcano Plot
##Monoculture vs. Biculture
final_table <- read.delim("Results_Proteus_LFQ_CBOLT/Tables/final_table.txt",sep = "\t",header = T)
final_table$logFC <- as.numeric(final_table$logFC_CBOLT_vs_Biculture)
final_table$adj.P.Val <- as.numeric(final_table$adj.P.Val_CBOLT_vs_Biculture)
p <- ggplot(data=final_table, aes(x=logFC, y=-log10(adj.P.Val))) + geom_point() + theme_bw()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# add a column of NAs
final_table$diffexpressed <- "NO"
# if log2Foldchange > 1.0 and pvalue < 0.05, set as "UP" 
final_table$diffexpressed[final_table$logFC > 1.0 & final_table$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
final_table$diffexpressed[final_table$logFC < -1.0 & final_table$adj.P.Val < 0.05] <- "DOWN"
# Re-plot but this time color the points with "diff expressed"
p <- ggplot(data=final_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + geom_point() + theme_bw()
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
## Change point color 
# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
#write the names of the proteins
ggplot(data=final_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=protein)) + 
  geom_point() + 
  theme_bw() +
  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
g3 <- ggplot(data=final_table, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=protein)) +
  geom_point() + 
  theme_bw() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-1.0, 1.0), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
ggsave(filename = paste(wd,'/plots/volcano_plot.png',sep = ""),plot = g3,width = 20,
       height =15)


