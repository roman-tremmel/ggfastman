#' Function to generate zoomed locus plots including LD 
#'
#' @param data A data.frame. Columns 'chr', 'pos', and 'pvalue' are required
#' @param build Genomic build. Currently only hg19 is supported.
#' @param snp Default top snp with lowest pvalue or rsid
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
fast_locusplot=function(data, build="hg19",snp="top",token=NULL){
  if (is.null(token)){
    stop('Register for token at https://ldlink.nci.nih.gov/?tab=apiaccess')
  }
  # get first top SNP
  if(snp  == "top"){ snp <- head(gwas_data[which.min(gwas_data$pvalue),]$rsid, 1)}
  # LINKAGE -----------------------------------------------------------------
  df_proxies <- LDlinkR::LDproxy(snp,pop="CEU",r2d="r2", token = token)
  # claculate genomic range
  gwas_data_topsnp <- gwas_data[which(gwas_data$rsid == snp),]
  gwas_data_topsnp$pos_min <- gwas_data_topsnp$pos -50000
  gwas_data_topsnp$pos_max <- gwas_data_topsnp$pos +50000
  # and update SNP data
  dd <- gwas_data[gwas_data$chr == gwas_data_topsnp$chr & gwas_data$pos >= gwas_data_topsnp$pos_min & gwas_data$pos <= gwas_data_topsnp$pos_max,]
  dd <- base::merge(dd, df_proxies, by.x = "rsid", by.y = "RS_Number")
  
  p1 <- ggplot2::ggplot(dd,ggplot2::aes(x=pos,y=-log10(pvalue),color=R2, shape = ifelse(rsid == snp, 17,16))) + 
        ggplot2::geom_segment(x = gwas_data_topsnp$pos,xend = gwas_data_topsnp$pos,y=-Inf, yend = -log10(gwas_data_topsnp$pvalue),linetype = 2, color = 1) +
        ggplot2::geom_point(size=3) +
        ggrepel::geom_label_repel(data = dd[dd$rsid ==snp,], ggplot2::aes(label = rsid))+
        ggplot2::scale_color_viridis_c(direction = -1) +
        ggplot2::scale_shape_identity(guide = "none")  +
        ggplot2::scale_y_continuous(name=expression(-log[10](italic(p)))) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(legend.position = c(0.9, 0.8)) 
  
  # GENOMIC CONTEXT ----------------------------------------------------------
  q=GenomicRanges::GRanges(seqnames=gwas_data_topsnp$chr,
                           ranges=IRanges::IRanges(start = gwas_data_topsnp$pos_min, end = gwas_data_topsnp$pos_max))
  p2 <- ggplot2::ggplot() + 
    ggbio::geom_alignment(Homo.sapiens::Homo.sapiens, which = IRanges::subsetByOverlaps(biovizBase::genesymbol, q))  +
    ggplot2::geom_vline(xintercept = gwas_data_topsnp$pos,linetype = 2) + 
    ggplot2::theme_classic(base_size = 14)
  # PLOT --------------------------------------------------------------------
  
   ggbio::tracks(p1, p2,heights = c(0.7,0.3),xlim = c(gwas_data_topsnp$pos_min, gwas_data_topsnp$pos_max))
}