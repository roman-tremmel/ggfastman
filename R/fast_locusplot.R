#' Function to generate zoomed locus plots including linkage diequilibrium using r2 
#'
#' @param data A data.frame. Columns 'chr', 'pos','pvalue' and 'rsid' are required
#' @param build Genomic build. Currently only hg19 is supported.
#' @param snp Default: top snp with lowest pvalue is taken. Otherwise specify any rsid id included in data.
#' @param pop a 1000 Genomes Project population, (e.g. YRI or CEU), multiple allowed, default = "CEU".
#' @param show_MAF Show minor allele frequencies (MAF) using different sizes of points. Default=FALSE.  
#' @param show_regulom Show most significant RegulomeDB ranks. More information at https://regulomedb.org/regulome-help/. Default=FALSE.  
#' @param RDB_ranks Specify ranks to show. Default= 1 and 2 including subtypes a-f.  
#' @param color color scheme for R^2 values. Possible values are c("viridis","magma","inferno","plasma").
#' @param token !!!REQUIRED!!!: Register for token at https://ldlink.nci.nih.gov/?tab=apiaccess
#' @return A ggplot2 object/plot
#' @importFrom LDlinkR LDproxy
#' @importFrom Homo.sapiens Homo.sapiens
#' @importFrom ggbio tracks 
#' @importFrom ggbio geom_alignment 
#' @importFrom ggrepel geom_label_repel 
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges GRanges
#' @export
#' @examples
#' fast_locusplot(gwas_data, token = "mytoken")
fast_locusplot=function(data, build="hg19",snp="top",pop="CEU",show_MAF=FALSE,show_regulom=FALSE,RDB_ranks=c("1","2"),color="viridis",token=NULL){
  gwas_data <- data
  if (is.null(token)){
    stop('Register for token at https://ldlink.nci.nih.gov/?tab=apiaccess')
  }
  if (!all(c('chr','pos','pvalue','rsid') %in% colnames(data))){
    stop('data must have columns "chr", "pos" and "pvalue" as well as "rsid"')
  }
  if(is.numeric(data$chr)){
    print("Numeric 'chr' column detected. Tried to transform to character by prefixing 'chr'.")
    data$chr <- paste0("chr", data$chr)
  }
  
  if (!is.numeric(data$pos)) 
    stop("pos column should be numeric.")
  
  if (!is.numeric(data$pvalue)) 
    stop("pvalue column should be numeric.")
  
  # get first top SNP
  if(snp  == "top"){ snp <- head(data[which.min(gwas_data$pvalue),]$rsid, 1)}
  if(!any(snp == gwas_data$rsid)){
    stop("specified rsid is not detected in data!")
  }
  
  # LINKAGE -----------------------------------------------------------------
  df_proxies <- LDlinkR::LDproxy(snp,pop="CEU",r2d="r2", token = token)
  # claculate genomic range
  gwas_data_topsnp <- gwas_data[which(gwas_data$rsid == snp),]
  gwas_data_topsnp$pos_min <- gwas_data_topsnp$pos -50000
  gwas_data_topsnp$pos_max <- gwas_data_topsnp$pos +50000
  # and update SNP data
  dd <- gwas_data[gwas_data$chr == gwas_data_topsnp$chr & gwas_data$pos >= gwas_data_topsnp$pos_min & gwas_data$pos <= gwas_data_topsnp$pos_max,]
  dd <- base::merge(dd, df_proxies, all.x = T, by.x = "rsid", by.y = "RS_Number")
  dd$MAF[is.na(dd$MAF)] <- 0.1
 
  p1 <- ggplot2::ggplot(dd,ggplot2::aes(x=pos,y=-log10(pvalue), color=R2)) + 
        ggplot2::geom_segment(x = gwas_data_topsnp$pos,xend = gwas_data_topsnp$pos,y=-Inf, yend = -log10(gwas_data_topsnp$pvalue),linetype = 2, color = 1) +
        ggplot2::scale_color_viridis_c(direction = -1,breaks=seq(0,1,.2), na.value = "grey",option = color) +
        ggplot2::scale_shape_identity(guide = "none")  +
        ggplot2::scale_y_continuous(name=expression(-log[10](italic(p)))) +
        ggplot2::theme_classic(base_size = 14) 
        
  
  if(show_MAF){
    if(show_regulom){
      p1 <- p1 + 
        ggplot2::geom_hline(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(yintercept = -log10(pvalue)), linetype = 3) +
        ggplot2::geom_point(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(shape = ifelse(rsid == snp, 24,21)), size=6.5) +
        ggrepel::geom_text_repel(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(label = paste0(rsid, ";RegDB= ",RegulomeDB)), size=3) +
        ggplot2::scale_size_continuous(range = c(1,5))}
    
    p1 <- p1 +  ggplot2::geom_point(ggplot2::aes(size=MAF,shape = ifelse(rsid == snp, 17,16))) +
      ggrepel::geom_label_repel(data = dd[dd$rsid ==snp,], color = 1, ggplot2::aes(label = rsid), min.segment.length = ggplot2::unit(2,"mm"))+
      ggplot2::theme(legend.position = c(0.9, 0.7)) 
    
  }else{
   if(show_regulom){
    p1 <- p1 + 
         ggplot2::geom_hline(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(yintercept = -log10(pvalue)), linetype = 3)+
         ggplot2::geom_point(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(shape = ifelse(rsid == snp, 24,21)), size=5) +
         ggrepel::geom_text_repel(data = subset(dd, grepl(paste(RDB_ranks, collapse="|"), dd$RegulomeDB)), ggplot2::aes(label = paste0(rsid, ";RegDB= ",RegulomeDB)), size=3) 
   }
    p1 <- p1 +  ggplot2::geom_point(ggplot2::aes(shape = ifelse(rsid == snp, 17,16)), size=3) +
      ggrepel::geom_label_repel(data = dd[dd$rsid ==snp,], color = 1, ggplot2::aes(label = rsid), min.segment.length = ggplot2::unit(2,"mm"))+
      ggplot2::theme(legend.position = c(0.9, 0.7)) 
  }
  

  
  # GENOMIC CONTEXT ----------------------------------------------------------
  q=GenomicRanges::GRanges(seqnames=gwas_data_topsnp$chr,
                           ranges=IRanges::IRanges(start = gwas_data_topsnp$pos_min, end = gwas_data_topsnp$pos_max))
  p2 <- ggplot2::ggplot() + 
    ggbio::geom_alignment(Homo.sapiens::Homo.sapiens, which = IRanges::subsetByOverlaps(biovizBase::genesymbol, q))  +
    ggplot2::geom_vline(xintercept = gwas_data_topsnp$pos,linetype = 2) + 
    ggplot2::theme_classic(base_size = 14)
  # PLOT --------------------------------------------------------------------
   ggbio::tracks(p1, p2, heights = c(0.7,0.3), xlim = c(gwas_data_topsnp$pos_min, gwas_data_topsnp$pos_max))
}