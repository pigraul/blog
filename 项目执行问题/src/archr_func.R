#conda activate ArchR
library(AnnotationDbi, lib.loc = "/Personal/wangxiang/R_lib/")
library (readr)
library(ArchR)
library(pheatmap,lib.loc="/SGRNJ01/Public/Software/conda_env/ArchR/lib/R/library/")
library (Seurat)
library (tidyr)
library (dplyr)
library (chromVARmotifs)


color.use<-c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple","DarkSlateGray","Gold","DeepPink2","Red4","#4682B4","#FFDAB9","#708090","#836FFF","#CDC673","#CD9B1D","#FF6EB4","#CDB5CD","DarkGreen","#008B8B","#43CD80","#483D8B","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62","#68228B","#CDB7B5","#CD853F","#6B8E23","#696969","#7B68EE","#9F79EE","#B0C4DE","#7A378B","#66CDAA","#EEE8AA","#00FF00","#EEA2AD","#A0522D","#000080","#E9967A","#00CDCD","#8B4500","#DDA0DD","#EE9572","#EEE9E9","#8B1A1A","#8B8378","#EE9A49","#EECFA1","#8B4726","#8B8878","#EEB4B4","#C1CDCD","#8B7500","#0000FF","#EEEED1","#4F94CD","#6E8B3D","#B0E2FF","#76EE00","#A2B5CD","#548B54","#BBFFFF","#B4EEB4","#00C5CD","#008B8B","#7FFFD4","#8EE5EE","#43CD80","#68838B","#00FF00","#B9D3EE","#9ACD32","#00688B","#FFEC8B","#1C86EE","#CDCD00","#473C8B","#FFB90F","#EED5D2","#CD5555","#CDC9A5","#FFE7BA","#FFDAB9","#CD661D","#CDC5BF","#FF8C69","#8A2BE2","#CD8500","#B03060","#FF6347","#FF7F50","#CD0000","#F4A460","#FFB5C5","#DAA520","#CD6889","#32CD32","#FF00FF","#2E8B57","#CD96CD","#48D1CC","#9B30FF","#1E90FF","#CDB5CD","#191970","#E8E8E8","#FFDAB9")
source('/Public/Script/shouhou/SCRIPT/Seurat_Monocle_modify/color_protocol.R')
color.use <- c(color_protocol, color.use)
species_Builtin <- c("hg19", "hg38", "mm10", "mm9")

# ### ignore the chromosome prefixes
# #addArchRChrPrefix(chrPrefix = FALSE)
# 修正 fragments 染色体名 + 提取存在的染色体
fixFragmentChrNamesAndExtractChroms <- function(inputFiles,
                                                 validChroms = NULL,
                                                 chromSizes = NULL,
                                                 bgzip_path = "/SGRNJ/Public/Software/conda_env/celescope_atac1.3.0/bin/bgzip",
                                                 tabix_path = "/SGRNJ/Public/Software/conda_env/celescope_atac1.3.0/bin/tabix") {
  fixedFiles <- list()
  allChroms <- c()
  
  for (sample in names(inputFiles)) {
    fragPath <- inputFiles[[sample]]
    frags <- data.table::fread(fragPath, header = FALSE)
    
    # 修改染色体名
    #frags$V1 <- paste0("chr", frags$V1)
    frags <- frags %>%
        mutate(V1 = case_when(
            startsWith(V1, "chr") ~ V1,    # 已有 "chr" 则不变
            TRUE ~ paste0("chr", V1)       # 否则添加 "chr"
    ))
    frags$V1[frags$V1 == "chrMT"] <- "chrM"
    
    # 收集所有染色体
    allChroms <- union(allChroms, unique(frags$V1))
    
    # 按需保留 validChroms
    if (!is.null(validChroms)) {
      frags <- frags[frags$V1 %in% validChroms]
    }
    
    # 如果提供了 chromSizes，就去除超长片段
    if (!is.null(chromSizes)) {
      chrom_dt <- data.table::data.table(chr = names(chromSizes), chrom_len = as.numeric(chromSizes))
      setnames(frags, c("V1", "V2", "V3", "V4", "V5"), c("chr", "start", "end", "cell", "mapq"))
      frags <- merge(frags, chrom_dt, by = "chr", all.x = TRUE)
      
      before <- nrow(frags)
      frags <- frags[end <= chrom_len]
      after <- nrow(frags)
      message(sprintf("[%s] Removed %d out-of-bound fragments", sample, before - after))
      
      # 去掉 chrom_len 列
      frags[, chrom_len := NULL]
      
      # 恢复列名
      setnames(frags, c("chr", "start", "end", "cell", "mapq"), paste0("V", 1:5))
    }
    
    # 如果为空则跳过
    if (nrow(frags) == 0) {
      warning(sprintf("Sample %s has no fragments after filtering, skipping.", sample))
      next
    }
    
    # 输出文件名
    fixedPath_raw <- sub("\\.tsv\\.gz$", "_fixed.tsv", fragPath)
    data.table::fwrite(frags, fixedPath_raw, sep = "\t", col.names = FALSE)
    
    # 压缩并索引
    fixedPath <- paste0(fixedPath_raw, ".gz")
    system(paste(bgzip_path, "-f", fixedPath_raw))
    system(paste(tabix_path, "-f -s 1 -b 2 -e 3", fixedPath))
    
    fixedFiles[[sample]] <- fixedPath
  }
  
  return(list(files = fixedFiles, chroms = allChroms))
}


# 修剪 geneAnnotation 和 genomeAnnotation 对象
pruneAnnotations <- function(geneAnnotation, genomeAnnotation, chroms) {
    # 获取 genomeAnnotation 染色体名，支持 GRanges/Seqinfo 或命名向量
    genome_chroms <- if (inherits(genomeAnnotation$chromSizes, "GenomicRanges") || inherits(genomeAnnotation$chromSizes, "Seqinfo")) {
        seqlevels(genomeAnnotation$chromSizes)
    } else {
        names(genomeAnnotation$chromSizes)
    }

    valid_chroms <- intersect(
        chroms,
        intersect(seqlevels(geneAnnotation$genes), genome_chroms)
    )

    if (length(valid_chroms) == 0) {
        stop("没有任何有效的染色体用于注释，请检查 fragments 和注释文件的染色体是否一致。")
    }

    geneAnnotation$genes <- keepSeqlevels(geneAnnotation$genes, valid_chroms, pruning.mode = "coarse")

    # 只有在 chromSizes 是 GRanges 或 Seqinfo 时才修剪
    if (inherits(genomeAnnotation$chromSizes, "GenomicRanges") || inherits(genomeAnnotation$chromSizes, "Seqinfo")) {
        genomeAnnotation$chromSizes <- keepSeqlevels(genomeAnnotation$chromSizes, valid_chroms, pruning.mode = "coarse")
    } else {
        # 如果是命名向量，只保留 valid_chroms
        genomeAnnotation$chromSizes <- genomeAnnotation$chromSizes[valid_chroms]
    }

    return(list(geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation))
}



# Create Arrow file
createArrow <- function(species, inputFiles, minTSS, minFrags,
                genomeAnnotation = genomeAnnotation, geneAnnotation = geneAnnotation) {
    if (species %in% c("hg38","hg19","mm10","mm9")) {
        print(111)
        ArrowFiles <- createArrowFiles(
            inputFiles = inputFiles,
            sampleNames = names(inputFiles),
            minTSS = minTSS,
            minFrags = minFrags,
            addTileMat = TRUE,
            addGeneScoreMat = TRUE
        )
    } else {
        print(222)
        chromSizes <- seqlengths(genomeAnnotation$chromSizes)
        validChroms <- seqlevels(geneAnnotation$genes)
        # 1. 修正 fragments 并提取染色体
        fragFixResult <- fixFragmentChrNamesAndExtractChroms(inputFiles, validChroms, chromSizes)
        inputFiles_fixed <- c(fragFixResult$files)
        present_chroms <- fragFixResult$chroms

        # 修剪注释
        # annoFix <- pruneAnnotations(geneAnnotation, genomeAnnotation, present_chroms)
        # geneAnnotation <- annoFix$geneAnnotation
        # genomeAnnotation <- annoFix$genomeAnnotation
        genomeAnnotation$blacklist <- GRanges(seqinfo = seqinfo(genomeAnnotation$chromSizes))

        # 创建 ArrowFiles
        inputFiles_fixed_vec <- unlist(inputFiles_fixed)
        ArrowFiles <- createArrowFiles(
            inputFiles = inputFiles_fixed_vec,
            sampleNames = names(inputFiles_fixed_vec),
            minTSS = minTSS,
            minFrags = minFrags,
            addTileMat = TRUE,
            addGeneScoreMat = TRUE,
            geneAnnotation = geneAnnotation,
            genomeAnnotation = genomeAnnotation
        )
    }
    return(ArrowFiles)
}


get_genome_size <- function(projHeme1){
    genome_annotation <- getGenomeAnnotation(projHeme1)
    chrom_sizes <- seqlengths(genome_annotation$chromSizes)
    x <- sum(chrom_sizes)
    formatted_x <- format(x, scientific = TRUE, digits = 3)
    return(formatted_x)
}

step1_Percell_Quality_Control <- function(projHeme1, outdir){
    dip <- paste0(outdir,"/Plots/00.QCmetrics/")
    print(dip)
    dir.create (dip,showWarnings=T) 
    ## log10(Unique Fragments) vs TSS enrichment score
    df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
    p <- ggPoint(
        x = df[,1], 
        y = df[,2], 
        colorDensity = TRUE,
        continuousSet = "sambaNight",
        xlabel = "Log10 Unique Fragments",
        ylabel = "TSS Enrichment",
        xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
        ylim = c(0, quantile(df[,2], probs = 0.99))
    ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
    #plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)
    plotPDF(p, name = file.path("00.QCmetrics","TSS-vs-Frags.pdf"), ArchRProj = projHeme1, addDOC = FALSE)
    print(file.path(dip,"TSS-vs-Frags.png"))
    png(file.path(dip,"TSS-vs-Frags.png"))
    print(p)
    dev.off()

    #### Make a violin plot for each sample for the TSS enrichment scores
    p2 <- plotGroups(
        ArchRProj = projHeme1, 
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "TSSEnrichment",
        plotAs = "violin",
        alpha = 0.4,
        addBoxPlot = TRUE
    )
    plotPDF(p2, name = file.path("00.QCmetrics","sample-TSS-VlnPlot.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
    print(file.path(dip,"sample-TSS-VlnPlot.png"))
    png(paste0(dip,"/sample-TSS-VlnPlot.png"))
    print(p2)
    dev.off()

    #### Make a violin plot for each sample for the log10(unique nuclear fragments).
    p4 <- plotGroups(
        ArchRProj = projHeme1, 
        groupBy = "Sample", 
        colorBy = "cellColData", 
        name = "log10(nFrags)",
        plotAs = "violin",
        alpha = 0.4,
        addBoxPlot = TRUE
    )
    plotPDF(p4, name = file.path("00.QCmetrics","sample-Frags-VlnPlot.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
    print(file.path(dip,"sample-Frags-VlnPlot.png"))
    png(file.path(dip,"sample-Frags-VlnPlot.png"))
    print(p4)
    dev.off()

    #### Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
    p5 <- plotFragmentSizes(ArchRProj = projHeme1)
    plotPDF(p5, name = file.path("00.QCmetrics","Sample-FragSizes.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
    print(file.path(dip,"Sample-FragSizes.png"))
    png(file.path(dip,"Sample-FragSizes.png"))
    print(p5)
    dev.off()

    p6 <- plotTSSEnrichment(ArchRProj = projHeme1)
    plotPDF(p6, name = file.path("00.QCmetrics","Sample-TSSProfile.pdf"), ArchRProj = projHeme1, addDOC = FALSE, width = 5, height = 5)
    print(file.path(dip,"Sample-TSSProfile.png"))
    png(file.path(dip,"Sample-TSSProfile.png"))
    print(p6)
    dev.off()
    #### Saving and Loading an ArchRProject  
    saveArchRProject(ArchRProj = projHeme1, load = FALSE)
    file.copy(file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme1.rds"),overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    return(projHeme1)
}


step2_Dimensionality_Reduction <- function(projHeme1, outdir, prefix, resol=1.2, dims=20){
    outdir_dim <- "01.DimReduction"
    dir.create (file.path(outdir,"Plots",outdir_dim), showWarnings=T)

    ##添加判断:根据细胞数更改分辨率
    if(nrow(projHeme1@cellColData) < 1500){
        resol <- 0.8
    }
    print(paste0("resolution is: ",resol))

    #### LSI implementation
    projHeme2 <- addIterativeLSI(
        ArchRProj = projHeme1,
        useMatrix = "TileMatrix",
        name = "IterativeLSI", 
        iterations = 2, 
        clusterParams = list(     #See Seurat::FindClusters
            resolution = c(0.2), 
            sampleCells = 10000, 
            n.start = 10
        ), 
        varFeatures = 25000, 
        dimsToUse = 1:dims
    )

    #### Cluster with ArchR 
    projHeme2 <- addClusters(
        input = projHeme2,
        reducedDims = "IterativeLSI",
        method = "Seurat",
        name = "Clusters",
        resolution = resol
    )

    ## UMAP plot
    projHeme2 <- addUMAP(
        ArchRProj = projHeme2,
        reducedDims = "IterativeLSI",
        name = "UMAP",
        nNeighbors = 30,
        minDist = 0.5,
        metric = "cosine"
    )
    p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")   
    plotPDF(p1, name = file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_samples.png")))
    print (p1)
    dev.off()
    p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    plotPDF(p2, name =file.path(outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".UMAP_IterativeLSI_Cluster.png")))
    print (p2)
    dev.off()

    ## tSNE plot
    projHeme2 <- addTSNE(
        ArchRProj = projHeme2,
        reducedDims = "IterativeLSI",
        name = "TSNE",
        perplexity = 30
    )
    p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
    plotPDF(p3, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_samples.png")))
    print (p3)
    dev.off()
    p4 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
    plotPDF(p4, name = file.path(outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.pdf")), ArchRProj = projHeme2, addDOC = FALSE, width = 5, height = 5)
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".tSNE_IterativeLSI_Cluster.png")))
    print (p4)
    dev.off()

    ## cluster prop plot
    pal <- paletteDiscrete(values = projHeme2$Clusters)
    cluster <- table (projHeme2$Clusters,projHeme2$Sample)
    Cluster<-rownames(cluster)
    tt<-cluster
    colnames(tt)<-"num"
    tt<-cbind(Cluster,tt)
    write.table (tt,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPerCluster.xls")),col.names = T,row.names=F,sep='\t',quote = F)
    freq_table <- prop.table (x=table (projHeme2$Clusters,projHeme2$Sample),margin = 2)
    Cluster<-rownames(freq_table)
    tt<-freq_table
    colnames(tt)<-"prop"
    tt<-cbind(Cluster,tt)
    write.table (tt,file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPerCluster.xls")),col.names = T,row.names=F,sep='\t',quote =F)
    
    data <- cluster
    data <- as.data.frame(rowSums(data))
    label <- paste(rownames(data),data[,1],sep = ': ')
    freq <- prop.table(x = as.matrix(data), margin = 2)
    colnames(freq) <- "Fraction"
    label <- paste0(label,' (',round(freq[,1]*100,2),"%)")
    freq <- as.data.frame(freq)
    freq$label <- label
    freq$color <- color.use[1:nrow(freq)]
    freq <- freq[order(freq$Fraction, decreasing = T),]
    freq$label <- factor(freq$label,levels = freq$label)
    freq$ymax = cumsum(freq$Fraction)
    freq$ymin = c(0, head(freq$ymax, n = -1))

    p <- ggplot(data = freq, aes(fill = label, ymax = ymax, ymin = ymin, xmax = 4, xmin = 2)) +
        geom_rect(colour = "grey30", show.legend = T) +
        coord_polar(theta = "y") +
        labs(x = "", y = "") +
        xlim(c(0, 4)) +
        theme_bw() +
        theme(panel.grid=element_blank()) +
        theme(axis.text=element_blank()) +
        theme(axis.ticks=element_blank()) +
        theme(panel.border=element_blank()) +
        theme(plot.title = element_text(size = 15, face = "bold",hjust = .5))+
        scale_fill_manual(values = freq$color)+
        theme(legend.title = element_blank(), legend.position = 'right')
    pdf (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.pdf")))
    print(p)
    dev.off()
    png (file.path(outdir,"Plots",outdir_dim,paste0(prefix,".CellsPropPer.png")))
    print(p)
    dev.off()

    saveArchRProject(ArchRProj = projHeme2, load = FALSE)
    file.copy(file.path(outdir,"Save-ArchR-Project.rds"), file.path(outrds,"projHeme2-reDim.rds"), overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    return(projHeme2)
}


step3_Cell_Annotation <- function(projHeme2, outdir, prefix, scRNA="no"){
    outdir_cellanno <- "02.Cellannotation"    
    dianno <- file.path(outdir,"Plots",outdir_cellanno)
    dir.create (dianno, showWarnings=T)

    if (scRNA == "no"){
        projHeme3<-projHeme2
        #embeddings
        pal <- color.use[1:length(unique(projHeme2@cellColData$Clusters))]
        names (pal) <- unique(projHeme2@cellColData$Clusters)
        projHeme3@cellColData$Clusters2 <- projHeme3@cellColData$Clusters
        #### saveRDS 
        saveArchRProject(ArchRProj = projHeme3, load = FALSE)
        file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme3.anno.rds"),overwrite=TRUE)
        file.remove(file.path(outdir,"Save-ArchR-Project.rds"))
    }else{
        # read scRNA
        seRNA <- readRDS (scRNA)
        seRNA$barcode <- rownames (seRNA@meta.data)
        seRNA$BioClassification <- as.character (seRNA@active.ident)
        remapClust <- unique (seRNA$BioClassification)
        names (remapClust) <- unique (seRNA$BioClassification)
        pal <- color.use[1:length(levels(seRNA))]
        names (pal) <- levels(seRNA)

        ### Defining Cluster Identity with scRNA-seq (Unconstrained Integration)
        projHeme3 <- addGeneIntegrationMatrix(
            ArchRProj = projHeme2, 
            useMatrix = "GeneScoreMatrix",
            matrixName = "GeneIntegrationMatrix",
            reducedDims = "IterativeLSI",
            seRNA = seRNA,
            addToArrow = TRUE,
            groupRNA = "BioClassification",
            nameCell = "predictedCell_Un", 
            nameGroup = "predictedGroup_Un",
            nameScore = "predictedScore_Un",
            force = TRUE
        )
            
        ####  Labeling scATAC-seq clusters with scRNA-seq information
        cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup_Un)
        labelOld <- rownames(cM)
        labelNew <- colnames(cM)[apply(cM, 1, which.max)]
        remapClust <- remapClust[names(remapClust) %in% labelNew]
        labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
        projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
        
        #### saveRDS 
        saveArchRProject(ArchRProj = projHeme3, load = TRUE)
        file.copy(file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme3.anno.rds"),overwrite=TRUE)
        file.remove(file.path(outdir,"Save-ArchR-Project.rds"))
    }    

    #### UMAP plot 
    p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2",baseSize=8,size = 0.2,quantileCut = c(0.025,0.975),embedding = 'UMAP')
    plotPDF(p1, name = file.path(outdir_cellanno,paste0(prefix,"-UMAP-RNA-Integration.pdf")), ArchRProj = projHeme3, addDOC = FALSE)
    
    p2 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Sample",baseSize=8,size = 0.2,quantileCut = c(0.025,0.975),embedding = 'UMAP')
    plotPDF(p2, name = file.path(outdir_cellanno,paste0(prefix,"-sample-UMAP-RNA-Integration.pdf")), ArchRProj = projHeme3, addDOC = FALSE)
    
    ## uMAP replot
    umap.file <- projHeme3@embeddings$UMAP$df
    colnames (umap.file) <- c("UMAP_1","UMAP_2")
    write.table (umap.file,file.path(dianno,paste0(prefix,".UMAP_Dimension.file.xls")),col.names = NA,sep='\t',quote = F)
    data <- data.frame (barcode = rownames (projHeme3),celltype = projHeme3@cellColData$Clusters2)
    data.umap <- bind_cols (data,umap.file)
    clustcol <- pal[unique (data.umap$celltype)]

    p3 <- ggplot (data.umap,mapping = aes (x=UMAP_1,y=UMAP_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
            theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
            guides(color = guide_legend(override.aes = list(size = 2)))    
    pdf (file.path(dianno,paste0(prefix,".labumap.integration.pdf")))
    print(p3)
    dev.off()
    png (file.path(dianno,paste0(prefix,".labumap.integration.png")))
    print(p3)
    dev.off()

    ## TSNE plot 
    tsne.file <- projHeme3@embeddings$TSNE$df
    colnames (tsne.file) <- c("TSNE_1","TSNE_2")
    write.table (tsne.file,file.path(dianno,paste0(prefix,".TSNE_Dimension.file.xls")),col.names = NA,sep='\t',quote = F)
    data <- data.frame (barcode = rownames (projHeme3),celltype = projHeme3@cellColData$Clusters2)
    data.tsne <- bind_cols (data,tsne.file)
    clustcol <- pal[unique (data.tsne$celltype)]

    p3 <- ggplot (data.tsne,mapping = aes (x=TSNE_1,y=TSNE_2,color=celltype)) + geom_point(size = 0.2)  + scale_color_manual(values = clustcol) + 
            theme_classic() + theme(legend.title=element_blank(),legend.text = element_text (size = 10)) +
            guides(color = guide_legend(override.aes = list(size = 2)))
    pdf (file.path(dianno,paste0(prefix,".labtsne.integration.pdf")))
    print(p3)
    dev.off()
    png (file.path(dianno,paste0(prefix,".labtsne.integration.png")))
    print(p3)
    dev.off()

    ## cluster prop plot
    data <- data.frame (celltype = names (pal),pal) 
    data <- data[unique (projHeme3$Clusters2),]
    data <- data[order(data$celltype),]
    colnames(data)<-c("celltype","pal")
    cluster <- table (projHeme3$Clusters2,projHeme3$Sample)
    write.table (cluster,file.path(dianno,paste0(prefix,".CellsPerCellType.xls")),col.names = NA,sep='\t',quote = F)
    freq_table <- prop.table (x=table (projHeme3$Clusters2,projHeme3$Sample),margin = 2)
    write.table (freq_table,file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.xls")),col.names = NA,sep='\t',quote =F)

    data <- cluster
    data <- as.data.frame(rowSums(data))
    label <- paste(rownames(data),data[,1],sep = ': ')
    freq <- prop.table(x = as.matrix(data), margin = 2)
    colnames(freq) <- "Fraction"
    label <- paste0(label,' (',round(freq[,1]*100,2),"%)")
    freq <- as.data.frame(freq)
    freq$label <- label
    freq$color <- color.use[1:nrow(freq)]
    freq <- freq[order(freq$Fraction, decreasing = T),]
    freq$label <- factor(freq$label,levels = freq$label)
    freq$ymax = cumsum(freq$Fraction)
    freq$ymin = c(0, head(freq$ymax, n = -1))

    p <- ggplot(data = freq, aes(fill = label, ymax = ymax, ymin = ymin, xmax = 4, xmin = 2)) +
        geom_rect(colour = "grey30", show.legend = T) +
        coord_polar(theta = "y") +
        labs(x = "", y = "") +
        xlim(c(0, 4)) +
        theme_bw() +
        theme(panel.grid=element_blank()) +
        theme(axis.text=element_blank()) +
        theme(axis.ticks=element_blank()) +
        theme(panel.border=element_blank()) +
        theme(plot.title = element_text(size = 15, face = "bold",hjust = .5))+
        scale_fill_manual(values = freq$color)+
        theme(legend.title = element_blank(), legend.position = 'right')
    pdf (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.pdf")))
    print(p)
    dev.off()
    png (file.path(dianno,paste0(prefix,".CellsPropPerCellTypes.png")))
    print(p)
    dev.off()    

    cM <- confusionMatrix(paste0(projHeme3$Clusters2), paste0(projHeme3$Sample))
    freq_table <- as.data.frame (cM) 
    write.table (freq_table,file.path(outdir,"Plots",outdir_cellanno,paste0 (prefix,".CellsPerCellType.xls")),col.names=NA,sep="\t",quote=F)

    return(projHeme3)
}


step4_Call_Peaks <- function(projHeme3, outdir, prefix, genomeSize=NULL, minCells=50, testmethod=testmethod){
    projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2",force = TRUE)
    saveArchRProject(ArchRProj = projHeme4, load = FALSE)
    
    ##  Calling Peaks w/ Macs2
    message("[INFO] Filter cell number > ", minCells, ", can not find marker peak under this threshold.") 
    pathToMacs2 <- findMacs2()
    if (is.null(genomeSize)) {
        projHeme4 <- addReproduciblePeakSet(
            ArchRProj = projHeme4, 
            groupBy = "Clusters2", 
            pathToMacs2 = pathToMacs2,
            minCells = minCells
        )    
    }else {
        projHeme4 <- addReproduciblePeakSet(
            ArchRProj = projHeme4, 
            groupBy = "Clusters2", 
            pathToMacs2 = pathToMacs2,
            minCells = minCells,
            genomeSize = genomeSize
        )
    }

    getPeakSet(projHeme4)

    df1 <- DataFrame (iranges = getPeakSet(projHeme4)) 
    peak <- data.frame (celltype = rownames (df1),data.frame (getPeakSet(projHeme4)))
    dio_macscallpeak <- file.path(outdir,"Plots","03.1.MacsCallPeak")
    dir.create (dio_macscallpeak)
    for (ct in unique (peak$celltype)) {
        rs <- filter (peak,celltype == ct)
        rs <- dplyr::select (rs,-N)
        tt<- dplyr::select(rs,c("celltype","seqnames","start","end","width","distToGeneStart","nearestGene","peakType","distToTSS"))
        write.table (tt,file.path (dio_macscallpeak,paste0(ct,".macsCallpeak.xls")),sep='\t',quote =F,row.names = F,col.names = T)

    }
    write.table (tt,file.path (dio_macscallpeak,"Report.macsCallpeak.xls"),sep='\t',quote =F,row.names = F,col.names = T)

    ### Add Peak Matrix
    projHeme4 <- addPeakMatrix(projHeme4) 
    saveArchRProject(ArchRProj = projHeme4, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme4.peak.rds"),overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    ##  Identifying Marker Peaks with ArchR
    ### testMethod = "wilcoxon" ## "wilcoxon" "ttest", and "binomial".
    if (testmethod != "binomial") {
        markersPeaks <- getMarkerFeatures(
            ArchRProj = projHeme4, 
            useMatrix = "PeakMatrix", 
            groupBy = "Clusters2",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = testmethod,
            binarize = TRUE
        )             
    }else{
        markersPeaks <- getMarkerFeatures(
            ArchRProj = projHeme4, 
            useMatrix = "PeakMatrix", 
            groupBy = "Clusters2",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "binomial",
            binarize = TRUE
        )
    }

    ## markergene
    markersGS <- getMarkerFeatures(
        ArchRProj = projHeme4, 
        useMatrix = "GeneScoreMatrix", 
        groupBy = "Clusters2",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )
    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    mm <- as.data.frame (markerList) 
    markerGenes <- mm$name
    if (argv$species %in% species_Builtin) {
        all.genes <- getFeatures(projHeme4)
    } else {
        all.genes <- getFeatures(projHeme4) %>%
            as.data.frame() %>%
            tidyr::separate(col = ".", sep = ":", into = c("chr", "name"))
        all.genes <- all.genes$name
    }

    ### markerpeak 
    dio_markerpeak <- file.path(outdir,"Plots","03.2.Markerpeak")
    dir.create (dio_markerpeak)
    markerPeaksList_1 <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1")
    markerPeaksList_2 <- getMarkers(markersPeaks, cutOff = "FDR <= 1 & Log2FC >= 1") 
    fdr<-0.1
    for (ctp in names (markerPeaksList_1)) {
        print (ctp)
        markerpeak <- data.frame (markerPeaksList_1[ctp])
        markerpeak <- dplyr::select (markerpeak,-c("group","group_name"))
        markerpeak <- dplyr::select (markerpeak,c("idx","seqnames","start","end","Log2FC","FDR","MeanDiff"))
        if (nrow (markerpeak) > 0) {
            FDR<-0.1
            print(FDR)
            fdr<-c(FDR,fdr)
            write.table (markerpeak,file.path(dio_markerpeak,paste0(ctp,".markerpeak.xls")),sep='\t',quote=F,row.names = F,col.names = T)
            ### ### Marker Peak MA and Volcano Plots
            pma <- plotMarkers(seMarker = markersPeaks, name = ctp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
            pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.pdf")))
            print (pma)
            dev.off()
            png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.png")))
            print (pma)
            dev.off()

            pv <- plotMarkers(seMarker = markersPeaks, name = ctp, cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
            pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.pdf")))
            print (pv)
            dev.off()
            png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.png")))
            print (pv)
            dev.off()
            ### Marker Peak in Browser Tracks        
            p <- try({
                plotBrowserTrack(
                    ArchRProj = projHeme4, 
                    groupBy = "Clusters2", 
                    geneSymbol = sample(all.genes, 1),
                    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[ctp],
                    upstream = 10000,
                    downstream = 20000
                )
            }, silent = TRUE)
            if (inherits(p, "try-error")) {
                message("[WARNING] plotBrowserTrack Failed! Skip plotting for cluster ", ctp)
            } else {
                pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.pdf")))
                print (grid::grid.draw(p[[names(p)]]))
                dev.off()
                png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.png")))
                print (grid::grid.draw(p[[names(p)]]))
                dev.off() 
            }
        }else{
            ### not statistically significant markers
            markerpeak <- data.frame (markerPeaksList_2[ctp])
            markerpeak <- markerpeak[order(markerpeak$FDR),]
            markerpeak <- markerpeak[1,]
            FDR<-markerpeak$FDR[1]
            cutoff<-paste0("FDR <= ",eval(parse(text = "FDR"))," & Log2FC >= 1")
            print(cutoff)
            FDR_fdr<-markerpeak$FDR[1]
            print(FDR)
            fdr<-c(FDR_fdr,fdr)
            markerpeak <- dplyr::select (markerpeak,-c("group","group_name"))
            markerpeak <- dplyr::select (markerpeak,c("idx","seqnames","start","end","Log2FC","FDR","MeanDiff"))
            write.table (markerpeak,file.path(dio_markerpeak,paste0(ctp,".markerpeak.xls")),sep='\t',quote=F,row.names = F,col.names = T)
            ### ### Marker Peak MA and Volcano Plots
            pma <- plotMarkers(seMarker = markersPeaks, name = ctp, cutOff = cutoff, plotAs = "MA")
            pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.pdf")))
            print (pma)
            dev.off()
            png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.MA.png")))
            print (pma)
            dev.off()

            pv <- plotMarkers(seMarker = markersPeaks, name = ctp, cutOff = cutoff, plotAs = "Volcano")
            pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.pdf")))
            print (pv)
            dev.off()
            png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Volcano.png")))
            print (pv)
            dev.off()
            ### Marker Peak in Browser Tracks        
            p <- try({
                plotBrowserTrack(
                    ArchRProj = projHeme4, 
                    groupBy = "Clusters2", 
                    geneSymbol = sample(all.genes, 1),
                    features = getMarkers(markersPeaks, cutOff = cutoff, returnGR = TRUE)[ctp],
                    upstream = 10000,
                    downstream = 20000
                )
            }, silent = TRUE)
            if (inherits(p, "try-error")) {
                message("[WARNING] plotBrowserTrack Failed! Skip plotting for cluster ", ctp)
            } else {
                pdf (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.pdf")))
                print (grid::grid.draw(p[[names(p)]]))
                dev.off()
                png (file.path(dio_markerpeak,paste0(ctp,".MarkerPeak.Tracks.png")))
                print (grid::grid.draw(p[[names(p)]]))
                dev.off() 
            }
        }
    }

    ### Marker Peak Heatmaps 
    FDR<-max(fdr) 
    cutoff<-paste0("FDR <= ",eval(parse(text = "FDR"))," & Log2FC >= 1")
    print(cutoff)

    if (ncol(markersPeaks) <= 2){
            heatmapPeaks <- plotMarkerHeatmap(
            seMarker = markersPeaks, 
            plotLog2FC = TRUE,
            cutOff = cutoff,
            transpose = FALSE,
            returnMatrix = TRUE
        )
        heatmapPeaks <- as.matrix (heatmapPeaks)
        removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
        heatmapPeaks<-removeColsAllNa(heatmapPeaks)
        write.table (heatmapPeaks,file.path(dio_markerpeak,"MarkerPeakHeatmaps.xls"),sep='\t',quote=F,row.names =T,col.names = T)
        pdf (file.path(dio_markerpeak,"MarkerPeak.Heatmap.pdf"),width = 7, height = 20)
        print(pheatmap (heatmapPeaks,color =  paletteContinuous(set = "comet", n = 100),show_colnames=F,cluster_col = F))
        dev.off()
        png (file.path(dio_markerpeak,"MarkerPeak.Heatmap.png"),height = 950)
        print(pheatmap (heatmapPeaks,color =  paletteContinuous(set = "comet", n = 100),show_colnames=F,cluster_col = F))
        dev.off()
    }else{
        heatmapPeaks <- plotMarkerHeatmap(
            seMarker = markersPeaks, 
            cutOff = cutoff,
            transpose = TRUE,
            returnMatrix = TRUE
        )
        heatmapPeaks <- as.matrix (heatmapPeaks)
        write.table (heatmapPeaks,file.path(dio_markerpeak,"MarkerPeakHeatmaps.xls"),sep='\t',quote=F,row.names =T,col.names = T)
        pdf (file.path(dio_markerpeak,"MarkerPeak.Heatmap.pdf"),width = 7, height = 20)
        print(pheatmap (heatmapPeaks,color =  paletteContinuous(set = "comet", n = 100),show_colnames=F,cluster_col = F))
        dev.off()
        png (file.path(dio_markerpeak,"MarkerPeak.Heatmap.png"),height = 950)
        print(pheatmap (heatmapPeaks,color =  paletteContinuous(set = "comet", n = 100),show_colnames=F,cluster_col = F))
        dev.off()
    }

    return (list(archr_proj = projHeme4, fdr = FDR))
}


step5_Motif <- function(projHeme4, outdir, prefix, FDR, testmethod=testmethod){
    dir.create (file.path(outdir,"Plots","04.MarkerpeakMotifEnrichment"),showWarnings=T)
    outdir_motif <- file.path(outdir,"Plots","04.MarkerpeakMotifEnrichment")

    ## Chapter 14 Motif and Feature Enrichment with ArchR
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme4, motifSet = "JASPAR2020", name = "Motif")
    projHeme5 <- addImputeWeights(projHeme5)

    ## re get markerpeaks
    if (testmethod != "binomial") {
        markersPeaks <- getMarkerFeatures(
            ArchRProj = projHeme5, 
            useMatrix = "PeakMatrix", 
            groupBy = "Clusters2",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = testmethod,
            binarize = TRUE
        )             
    }else{
        markersPeaks <- getMarkerFeatures(
            ArchRProj = projHeme5, 
            useMatrix = "PeakMatrix", 
            groupBy = "Clusters2",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            testMethod = "binomial",
            binarize = TRUE
        )
    }
    cutoff<-paste0("FDR <= ",eval(parse(text = "FDR"))," & Log2FC >= 0.5")
    enrichMotifs <- peakAnnoEnrichment(
        seMarker = markersPeaks,
        ArchRProj = projHeme5,
        peakAnnotation = "Motif",
        cutOff = cutoff
    )
    
    mat <- assays(enrichMotifs)[["mlog10Padj"]]
    if(max(mat)>1){
        if(max(mat)>20){
            cutoff<-20
        }else{
            cutoff<-1
        }
        print(cutoff)

        heatmapEM.table <- plotEnrichHeatmap(enrichMotifs, n = 500, transpose = FALSE,returnMatrix = TRUE,cutOff = cutoff)
        heatmapEM.table <- data.frame (motiftpm = rownames (heatmapEM.table),heatmapEM.table)
        ht <- separate (heatmapEM.table,motiftpm,into=c("motif"),sep= "_")
        write.table (ht,file.path(outdir_motif,"MarkerPeakEnrichMotif.xls"),sep='\t',quote=F,row.names =F,col.names = T)

        heatmapEM.plot <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = FALSE,returnMatrix = TRUE,cutOff = cutoff)
        heatmapEM.plot <- data.frame (motiftpm = rownames (heatmapEM.plot),heatmapEM.plot)
        ht.plot <- separate (heatmapEM.plot,motiftpm,into=c("motif"),sep= "_")
        ###########ERROR HERE!!! repeat motifs
        ht.plot<-ht.plot[!duplicated(ht.plot$motif),]
        rownames(ht.plot) <- ht.plot$motif
        ht.plot <- ht.plot[,-1]
        write.table (ht.plot,file.path(outdir_motif,"MarkerPeakEnrichMotif.Plot.xls"),sep='\t',quote=F,row.names =T,col.names = NA)
        heatmapPeaks<-ht.plot[,which(colSums(ht.plot) > 0)]
        ht.plot <- as.matrix (heatmapPeaks)
        pdf (file.path(outdir_motif,"MarkerPeakEnrichMotif.pdf"),height = 13)
        print(pheatmap (ht.plot,color =  paletteContinuous(set = "comet", n = 100),cluster_row = F))
        dev.off()
        png (file.path(outdir_motif,"MarkerPeakEnrichMotif.png"),height = 800)
        print(pheatmap (ht.plot,color =  paletteContinuous(set = "comet", n = 100),cluster_row = F))
        dev.off()
    }

    saveArchRProject(ArchRProj = projHeme5, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme5.motif.rds"),overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    ### Chapter 15 ChromVAR Deviations Enrichment with ArchR
    dir.create (file.path(outdir,"Plots","05.ChromVAR.Variable.motif"),showWarnings=T)
    outdir_chromvar <- file.path(outdir,"Plots","05.ChromVAR.Variable.motif")

    if("Motif" %ni% names(projHeme5@peakAnnotation)){
        projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
    }

    projHeme5 <- addBgdPeaks(projHeme5)
    projHeme5 <- addDeviationsMatrix(
        ArchRProj = projHeme5, 
        peakAnnotation = "Motif",
        force = TRUE,
        matrixName = "MotifMatrix"
    )

    plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE,n=5)
    pdf (file.path(outdir_chromvar,"Variable-Motif-Deviation-Scores.pdf"))
    print(plotVarDev)
    dev.off()
    png (file.path(outdir_chromvar,"Variable-Motif-Deviation-Scores.png"))
    print(plotVarDev)
    dev.off()
    ## extract a subset of motifs for downstream analysis #(select top5 motif)
    varmotif <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = FALSE) %>% as.data.frame ()
    varmotif <- dplyr::select (varmotif,-c("seqnames","idx"))
    write.table (varmotif,file.path(outdir_chromvar,"ALL.Var-motif-Deviations.xls"),row.names=F,col.names=T,sep='\t',quote=F)
    varmotif <- separate(varmotif,name,into=c("motif_name"),sep="_",remove=F)
    motifs <- varmotif$motif_name
    #motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
    markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
    head(markerMotifs)

    ### MotifMatrix contains seqnames for both z-scores and deviations, shown above by by "z:" and "deviations:".
    ###  remove z:SREBF1_22
    markerMotifs <- grep("z:", markerMotifs, value = TRUE)[1:5]

    # markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
    # markerMotifs
    ## plot the distribution of chromVAR deviation scores for each cluster. 
    for (motif in markerMotifs) {
        marker <- unlist(strsplit(motif,':'))[2]
        p <- plotGroups(ArchRProj = projHeme5, 
            groupBy = "Clusters2", 
            colorBy = "MotifMatrix", 
            name = motif,
            imputeWeights = getImputeWeights(projHeme5)
        )
        pdf (file.path(outdir_chromvar,paste0(marker,".Plot-Groups-Deviations-w-Imputation.pdf")),width = 5,height = 5)
        print (p)
        dev.off()
        png (file.path(outdir_chromvar,paste0(marker,".Plot-Groups-Deviations-w-Imputation.png")))
        print (p)
        dev.off()
    }

    ### 1.overlay the z-scores on our UMAP embedding as we have done previously for gene scores.
    for (motif in sort(markerMotifs)) {
        marker <- unlist(strsplit(motif,':'))[2]
        p <- plotEmbedding(
            ArchRProj = projHeme5,
            colorBy = "MotifMatrix", 
            name = motif, 
            embedding = "UMAP",
            imputeWeights = getImputeWeights(projHeme5)
        )
        pdf (file.path(outdir_chromvar,paste0(marker,".motif.umap.pdf")),width = 5,height = 5)
        print (p)
        dev.off()
        png (file.path(outdir_chromvar,paste0(marker,".motif.umap.png")))
        print (p)
        dev.off()
    }

    ## 2. overlay the gene scores for each of these TFs on the UMAP embedding
    markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
    markerRNA <- markerRNA[markerRNA %in% motifs]
    markerRNA <- markerRNA[1:2]
    for (gene in sort(markerRNA)) {
        p <- plotEmbedding(
            ArchRProj = projHeme5, 
            colorBy = "GeneScoreMatrix", 
            name = gene, 
            embedding = "UMAP",
            imputeWeights = getImputeWeights(projHeme5)
        )
        pdf (file.path(outdir_chromvar,paste0(gene,".motif.link.genescore.umap.pdf")),width = 5,height = 5)
        print (p)
        dev.off()
        png (file.path(outdir_chromvar,paste0(gene,".motif.link.genescore.umap.png")))
        print (p)
        dev.off()  
    }

    # ## plot the linked gene expression for each of these TFs on the UMAP embedding 
    if (scRNA == "no"){
        markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
        markerRNA <- markerRNA[markerRNA %in% motifs]
        markerRNA <- markerRNA[1:2]
        for (gene in sort(markerRNA)) {
            p <- plotEmbedding(
                ArchRProj = projHeme5, 
                colorBy = "GeneScoreMatrix", 
                name = gene, 
                embedding = "UMAP",
                continuousSet = "blueYellow",
                imputeWeights = getImputeWeights(projHeme5)
            )
            pdf (file.path(outdir_chromvar,paste0(gene,".motif.link.geneScore.pdf")),width = 5,height = 5)
            print (p)
            dev.off()
            png (file.path(outdir_chromvar,paste0(gene,".motif.link.geneScore.png")))
            print (p)
            dev.off()  
        }
    }else{
        markerRNA <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "GeneIntegrationMatrix")
        markerRNA <- markerRNA[markerRNA %in% motifs]
        markerRNA <- markerRNA[1:2]
        for (gene in sort(markerRNA)) {
            p <- plotEmbedding(
                ArchRProj = projHeme5, 
                colorBy = "GeneIntegrationMatrix", 
                name = gene, 
                embedding = "UMAP",
                continuousSet = "blueYellow",
                imputeWeights = getImputeWeights(projHeme5)
            )
            pdf (file.path(outdir_chromvar,paste0(gene,".motif.link.geneIntegrationScore.pdf")),width = 5,height = 5)
            print (p)
            dev.off()
            png (file.path(outdir_chromvar,paste0(gene,".motif.link.geneIntegrationScore.png")))
            print (p)
            dev.off()  
        }
    }

    saveArchRProject(ArchRProj = projHeme5, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme5.motif.rds"),overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    return (projHeme5)
}


step6_Co_Accessibility <- function(projHeme5, outdir, prefix){
    dir.create (file.path (outdir,"Plots","06.Co-accessibility"))
    outdir_coaccess <- file.path (outdir,"Plots","06.Co-accessibility")

    projHeme6 <- addCoAccessibility(
        ArchRProj = projHeme5,
        reducedDims = "IterativeLSI"
    )
    ### retrieve this co-accessibility information from the ArchRProject  ( correlation column gives the numeric correlation of the accessibility between those two peaks)
    cA <- getCoAccessibility(
        ArchRProj = projHeme6,
        corCutOff = 0.5,
        resolution = 1000,
        returnLoops = FALSE
    )

    allpeak.tpm <- DataFrame (iranges = metadata(cA)[[1]])
    allpeak <- data.frame (celltype = rownames (allpeak.tpm),data.frame (allpeak.tpm$iranges)) 
    allpeak <- data.frame (allpeak,peak = paste0(allpeak$seqnames,":",allpeak$start,"-",allpeak$end))%>% dplyr::select (celltype,peak)
    coaccess <- data.frame (cA) %>% dplyr::select (queryHits,subjectHits,seqnames,correlation)
    copeak <- data.frame (query = allpeak[coaccess$queryHits,],subject = allpeak[coaccess$subjectHits,],coaccess)
    write.table (copeak,file.path(outdir_coaccess,"co-accessibility.xls"),col.names=T,row.names=F,sep='\t',quote=F)

    ### Plotting browser tracks of Co-accessibility
    markersGS <- getMarkerFeatures(
        ArchRProj = projHeme6, 
        useMatrix = "GeneScoreMatrix", 
        groupBy = "Clusters2",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
    )

    markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
    mm <- as.data.frame (markerList) 
    if ( nrow(mm) != 0) {
        percelltype_top1 <- group_by (mm,group) %>% sample_n (,size=1)
        markerGenes <- percelltype_top1$name[1:5]
        write.table (mm,file.path(outdir_coaccess,"celltype.markerGeneList.xls"),col.names=T,row.names = F,sep='\t',quote=F)
    }else{
        print("No Marker Genes found !!!")
        all.genes<-getFeatures(projHeme6)
        markerGenes<-sample(all.genes,5)
        print("Used genes are:")
        print(markerGenes)
    }

    p <- plotBrowserTrack(
        ArchRProj = projHeme6, 
        groupBy = "Clusters2",
        geneSymbol = markerGenes,
        upstream = 50000,
        downstream = 50000,
        loops = getCoAccessibility(projHeme6)
    )

    for(i in 1:length(markerGenes)){
        pdf(paste0(outdir_coaccess,"/",markerGenes[i],"-with-CoAccessibility.pdf"))
        grid::grid.draw(p[[markerGenes[i]]])
        dev.off()
        png(paste0(outdir_coaccess,"/",markerGenes[i],"-with-CoAccessibility.png"))
        grid::grid.draw(p[[markerGenes[i]]])
        dev.off()
        print(paste0(markerGenes[i]," PeakPlot DONE!"))
    }

    saveArchRProject(ArchRProj = projHeme6, load = FALSE)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme6.CoAccessibility.rds"),overwrite=TRUE)
    file.remove(file.path(outdir,"Save-ArchR-Project.rds"))

    ### Peak2GeneLinkage with ArchR  
    if (scRNA != "no"){
        dir.create (file.path(outdir,"Plots","07.Peak2GeneLinkage"),showWarnings=T)

        outdir_peak2gene <- file.path(outdir,"Plots","07.Peak2GeneLinkage")
        projHeme6 <- addPeak2GeneLinks(
            ArchRProj = projHeme6,
            reducedDims = "IterativeLSI"
        )
        
        ### Plotting a heatmap of peak-to-gene links 
        pdf (file.path (outdir_peak2gene,"Peak2GeneLinks-heatmap.pdf"))
        plotPeak2GeneHeatmap(ArchRProj = projHeme6, groupBy = "Clusters2",k = 5,returnMatrices = FALSE)
        dev.off()
        # get heatmap 
        p <- plotPeak2GeneHeatmap(ArchRProj = projHeme6, groupBy = "Clusters2",k = 5,returnMatrices = TRUE)
        hp <- data.frame (p$Peak2GeneLinks,kmeansId = p$ATAC$kmeansId)
        hp <- hp[order(hp$kmeansId),]
        write.table (hp,file.path(outdir_peak2gene,"peak2GeneLinks.xls"),col.names=T,row.names=F,sep='\t',quote=F)

        p2g.genes <- hp[order(hp$Correlation),] %>%filter (gene != "NA") %>% dplyr::select (gene) %>% unique()
        p2g.genes.use <- p2g.genes[1:5,]   
        p <- plotBrowserTrack(
            ArchRProj = projHeme6, 
            groupBy = "Clusters2", 
            geneSymbol = p2g.genes.use, 
            upstream = 10000,
            downstream = 20000,
            loops = getPeak2GeneLinks(projHeme6)
        )

        for(i in 1:length(p2g.genes.use)){
            pdf(paste0(outdir_peak2gene,"/",p2g.genes.use[i],"-Peak2GeneLinks.pdf"))
            print(grid::grid.draw(p[[p2g.genes.use[i]]]))
            dev.off()
            png(paste0(outdir_peak2gene,"/",p2g.genes.use[i],"-Peak2GeneLinks.png"))
            print(grid::grid.draw(p[[p2g.genes.use[i]]]))
            dev.off()
            print(paste0(p2g.genes.use[i]," PeakPlot DONE!"))
        }

    saveArchRProject(ArchRProj = projHeme6)
    file.copy (file.path(outdir,"Save-ArchR-Project.rds"),file.path(outrds,"projHeme6.coaccess.peak2gene.rds"),overwrite=TRUE)
    }

    return (projHeme6)
}


























