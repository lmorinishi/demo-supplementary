getCustomTrack <- function(dmrI=1, gen="hg19", samps, pal, cpgdata, rranges, groupmeth=targets[,14], txfilter=NULL,
                           startloc=NULL, endloc=NULL, txlab.just="left", pval=NULL) {

  # extract chromosome number and location from DMR results
  coords <- ranges(rranges[dmrI])
  chrom <- as.character(seqnames(rranges[dmrI]))
  if (is.null(startloc)) {
    startloc <- start(coords)
    endloc <- end(coords)
  }
  
  # add 25% extra space to plot
  minbase <- (startloc - (0.25*(endloc-startloc))) 
  maxbase <- (endloc + (0.25*(endloc-startloc)))
  
  # cislHMM=islHMM[islHMM$chr==chrom,]
  # islandData <- GRanges(seqnames=Rle(cislHMM[,1]),
  #                       ranges=IRanges(start=cislHMM[,2], end=cislHMM[,3]),
  #                       strand=Rle(strand(rep("*",nrow(cislHMM)))))
  
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  
  # extract data on CpGs in DMR
  cpgData <- subsetByOverlaps(cpgdata, rranges[dmrI])
  
  # Transcript track
  if (is.null(txfilter)) {
    txTr2 <- BiomartGeneRegionTrack(start=startloc, end=endloc, chromosome=chrom, genome=gen, strand="+-", showId=T, just.group=txlab.just, stacking="full", name="Transcripts", geneSymbols=T)
  } else {
    txTr <- BiomartGeneRegionTrack(start=startloc, end=endloc, chromosome=chrom, genome=gen, strand="+-", showId=T, just.group=txlab.just,
                                   stacking="full", name="Transcripts", geneSymbols=T, filter=list(with_protein_id=T))
    txTr2 <- txTr[transcript(txTr) %in% txfilter]
  }
  
  # Methylation data track
  if (is.null(pval)) {
    methTrack <- DataTrack(range=cpgData, groups=groupmeth, genome = gen,
                           chromosome=chrom, ylim=c(0.05,0.7), col=pal,
                           type=c("a","p"), name="DNAme",
                           background.panel="white", legend=F, cex.title=0.8,
                           cex.axis=0.8, cex.legend=0.8)
    methTrack <- methTrack[samps]
  } else {
    my_cpgs = cpgData[seqnames(cpgData)==chrom & start(cpgData) > minbase & start(cpgData) < maxbase]#& cpgData$ranges > minbase & cpgData$ranges < maxbase,]
    cpg_list = names(my_cpgs)
    elementMetadata(my_cpgs) = df[match(cpg_list, df$CpG),]
    methTrack <- DataTrack(range=my_cpgs, data=my_cpgs$adj.pval, name="DNAme", type=c("m", "p"),
                           col.mountain="black", fill.mountain=c("gray95", "gray95"), col.axis="black",
                           background.panel="white")
  }
  displayPars(methTrack) <- list(background.title="transparent", fontcolor.title="black", rotate.title=F, showAxis=T)
  
  tracks <-  list(iTrack, gTrack, methTrack, txTr2 )#, islandTrack)
  return(list("tracks"=tracks, "minbase"=minbase, "maxbase"=maxbase, "sizes"=rep(3, length(tracks))))
}

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

plotBothGene <- function(gene, allgenes, shape=T,
                         basecell=rep(c("KG1", "Thp1"), each=4), 
                         tet2=rep(c("WT", "WT", "mut", "mut")),
                         drug=c("None", "2HG"), pvals=NULL, normed=F, dfcols=3:26) {
  require(dplyr)

  df = data.frame(RNA_raw = unlist(elementMetadata(allgenes)[allgenes$external_gene_name==gene,dfcols]),
                  basecell = basecell,
                  tet2 = tet2,
                  drug = drug)
  df$treatment <- paste(df$basecell, df$tet2, sep="_")
  df$tet2 <- factor(df$tet2, levels=c("WT", "mut"))
  df$drug <- factor(df$drug, levels=c("None", "2HG", "AraC"))
  
  # not modular!!
  if (normed) {
    x = df %>% group_by(basecell) %>% filter(tet2=="WT") %>% mutate(mean_norm=mean(RNA_raw))
    df$RNA_raw[df$basecell=="KG1"] <- df$RNA_raw[df$basecell=="KG1"]/mean(df$RNA_raw[df$basecell=="KG1" & df$tet2=="WT"])
    df$RNA_raw[df$basecell=="Thp1"] <- df$RNA_raw[df$basecell=="Thp1"]/mean(df$RNA_raw[df$basecell=="Thp1" & df$tet2=="WT"])
  }
  
  if (shape) {
    p1 <- ggplot(df, aes(x=basecell, y=RNA_raw, color=tet2, shape=drug)) + 
      geom_point(position=position_jitterdodge(dodge.width=.9)) + scale_color_manual(values=c("gray55", "firebrick3")) + 
      theme_bw() + ggtitle(gene) + xlab("") + ylab("Normalized RNA count")
  } else {
    p1 <- ggplot(df, aes(x=basecell, y=RNA_raw, color=tet2)) + 
      geom_point(position=position_jitterdodge(dodge.width=.9)) + scale_color_manual(values=c("gray55", "firebrick3")) + 
      theme_bw() + ggtitle(gene) + xlab("") + ylab("Normalized RNA count")
  }
  
  if (!is.null(pvals)) {
    pval_K = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"padj_K"]
    log2FC_K = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"log2FoldChange_K"]
    pval_T = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"padj_T"]
    log2FC_T = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"log2FoldChange_T"]
    p1 <- p1 +
      annotate("text", x=1, y=max(df$RNA_raw, na.rm=T)*1.05, label=paste0("pval: ", formatC(pval_K, format = "e", digits = 2))) +
      annotate("text", x=1, y=max(df$RNA_raw, na.rm=T)*1.1, label=paste0("log2FC: ", round(log2FC_K,2))) +
      annotate("text", x=2, y=max(df$RNA_raw, na.rm=T)*1.05, label=paste0("pval: ", formatC(pval_T, format = "e", digits = 2))) +
      annotate("text", x=2, y=max(df$RNA_raw, na.rm=T)*1.1, label=paste0("log2FC: ", round(log2FC_T,2))) +
      ylim(-1e-5, max(df$RNA_raw)*1.1)
  }
  p1
}

plotGeneList <- function(genes, allgenes, shape=T, samples=data.frame(basecell=rep(c("KG1", "Thp1"), each=4), 
                                                                     tet2=rep(c("WT", "WT", "mut", "mut")),
                                                                     drug=c("None", "2HG")), 
                         gtype="boxplot", pvals=NULL, normed=F, dfcols=3:26, ylim=NULL) {
  require(dplyr)
  
  n_measurements = length(dfcols) * length(genes)
  df = data.frame(RNA_raw = as.vector(t(as.matrix(elementMetadata(allgenes)[match(genes, elementMetadata(allgenes)$external_gene_name),dfcols]))),
                  gene = rep(genes, each=length(dfcols)),
                  basecell = rep(samples$basecell, length.out=n_measurements),
                  tet2 = rep(samples$tet2, length.out=n_measurements),
                  drug = rep(samples$drug, length.out=n_measurements))
  df$gene <- factor(df$gene, levels=genes)
  df$treatment <- paste(df$basecell, df$tet2, sep="_")
  df$tet2 <- factor(df$tet2, levels=c("WT", "mut"))
  df$drug <- factor(df$drug, levels=c("None", "2HG", "AraC"))
  
  # not modular!!
  if (normed) {
    x = df %>% group_by(basecell, gene) %>% filter(tet2=="WT") %>% select_(.dots = c("basecell", "drug", "gene", "RNA_raw")) %>% summarize(mean=mean(RNA_raw))
    df <- merge(df, x, by=c("basecell", "gene"))
    df$display <- df$RNA_raw/df$mean
  } else {
    df$display <- df$RNA_raw
  }
  
  if (gtype=="boxplot") {
    p1 <- ggplot(df, aes(x=gene, y=display, color=tet2, fill=tet2)) +
      geom_boxplot(position=position_dodge(0.9), color='black', outlier.shape=NA) +
      geom_point(position=position_jitterdodge(dodge.width=.9), shape=21, color='black', size=2) +
      scale_color_manual(values=c("gray55", "firebrick3")) + scale_fill_manual(values=c("gray55", "firebrick3")) +
      theme_bw() + xlab("") + ylab("Normalized RNA count")
  } else if (gtype=="bar") {
    df_summary = df %>% group_by(gene, tet2) %>% summarize(mean=mean(display), sd=sd(display), sem=sd(display)/sqrt(n()))
    p1 <- ggplot(data=df, aes(x=gene, color=tet2, fill=tet2)) +
      geom_col(data=df_summary, position=position_dodge(0.5), width=0.5, aes(y=mean)) +
      geom_point(data=df, position=position_jitterdodge(dodge.width=.5), shape=21, color='black', size=2, aes(y=display)) +
      geom_errorbar(data=df_summary, width=0.2, position=position_dodge(0.5), color='black', aes(ymax=mean+sem, ymin=mean-sem)) +
      scale_color_manual(values=c("gray55", "firebrick3")) + scale_fill_manual(values=c("gray55", "firebrick3")) +
      theme_classic() + xlab("") + ylab("Normalized RNA count") + theme(legend.position="top")
    if (is.null(ylim)) p1 <- p1 + scale_y_continuous(expand=c(0,0)) else p1 <- p1 + scale_y_continuous(expand=c(0,0), limits=ylim)
  }
  
  if (!is.null(pvals)) {
    pval_K = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"padj_K"]
    log2FC_K = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"log2FoldChange_K"]
    pval_T = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"padj_T"]
    log2FC_T = elementMetadata(allgenes)[allgenes$external_gene_name==gene,"log2FoldChange_T"]
    p1 <- p1 +
      annotate("text", x=1, y=max(df$RNA_raw, na.rm=T)*1.05, label=paste0("pval: ", formatC(pval_K, format = "e", digits = 2))) +
      annotate("text", x=1, y=max(df$RNA_raw, na.rm=T)*1.1, label=paste0("log2FC: ", round(log2FC_K,2))) +
      annotate("text", x=2, y=max(df$RNA_raw, na.rm=T)*1.05, label=paste0("pval: ", formatC(pval_T, format = "e", digits = 2))) +
      annotate("text", x=2, y=max(df$RNA_raw, na.rm=T)*1.1, label=paste0("log2FC: ", round(log2FC_T,2))) +
      ylim(-1e-5, max(df$RNA_raw)*1.1)
  }
  p1
}























