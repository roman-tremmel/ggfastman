#' A fast function for Manhattan plots
#' 
#' This is a very fast and easy-to-individualize plotting function for GWAS results e.g. pvalues based on \code{\link{ggplot2}} and \code{\link{scattermore}}. 
#'
#' @param gwas A data.frame. Columns "chr", "pos", and "pvalue" are required
#' @param build Genomic build. Currently hg18, hg19, and hg38 are supported. If you are not sure about, set the default "hg19".
#' @param color1 Color for odd-numbered chromosomes
#' @param color2 Color for even-numbered chromosomes
#' @param y_scale When an own y-scale is provided then set this to TRUE to avoid an error. 
#' @param log10p If TRUE (default), the -log10 transformed pvalues are plotted
#' @param speed The speed option. Fast and ultrafast use \code{\link{scattermore}} functionality
#' @param pointsize Only when using the 'fast' option you can increase pointsize. Default is 0. When using small pointsizes it could be that points are not shown in the RStudio Plots or Zoom window. But they will plotted when saving to pdf.
#' @param pixels Only when using the 'fast' option you can increase pixel width and height. Default is c(512, 512). 
#' @return A ggplot2 object/plot
#' @export
#' @examples
#' data("gwas_data")
#' head(gwas_data)
#' # slow
#' fast_manhattan(gwas_data, build='hg18')
#' # fast
#' fast_manhattan(gwas_data, build = "hg18", speed = "fast")
#' # ultrafast
#' fast_manhattan(gwas_data, build='hg18', speed = "ultrafast")
fast_manhattan=function(data,build="hg19",color1='black',color2='grey',y_scale = TRUE,log10p=TRUE,alpha = 1,
                   speed = "fast",pointsize=0, pixels=c(512, 512), ...){
  
  if (!all(c('chr','pos','pvalue') %in% colnames(data))){
    stop('data must have columns "chr", "pos" and "pvalue"')
  }
  
  if(is.numeric(data$chr)){
    print("Numeric 'chr' column detected. Tried to transform to character by prefixing 'chr'.")
    data$chr <- paste0("chr", data$chr)
  }
  
  if (!is.numeric(data$pos)) 
    stop("pos column should be numeric.")
  
  if (!is.numeric(data$pvalue)) 
    stop("pvalue column should be numeric.")
  
  if(log10p) {
    data$y <- -log10(data$pvalue)
  }
  else {
    data$y <- data$pvalue
  }
  
  build=match.arg(build, choices = c('hg18','hg19','hg38'))
  speed=match.arg(speed, choices = c("slow", "fast", "ultrafast"))
  chrom_lengths=get_chrom_lengths(build)[extract_which_chr(data)]
  data=add_cumulative_pos(data,chrom_lengths)
  data=add_color(data,color1 = color1, color2 = color2)
  data=add_shape(data,shape=16)
  x_breaks=get_x_breaks(chrom_lengths)
  
  color_map=unique(data$color)
  names(color_map)=unique(data$color)
  
  plot <- ggplot2::ggplot(data, ggplot2::aes(x=cumulative_pos,y=y,color=color,shape=shape))
  plot <- switch(speed,
                 slow = plot + ggplot2::geom_point(alpha =alpha),
                 fast  = plot + scattermore::geom_scattermore(pointsize=pointsize, pixels=pixels, alpha = alpha),
                 ultrafast = ggplot2::ggplot() + scattermore::geom_scattermost(data[, c("cumulative_pos", "y")])
  )
  
  if(y_scale) plot <- plot + ggplot2::scale_y_continuous(expand=c(0.01,0),name=expression(-log[10](italic(p))))
  
  return(
    plot + ggplot2::theme_classic()+
      ggplot2::scale_x_continuous(expand=c(0.01,0),breaks=x_breaks,
                                  labels=names(x_breaks),name='Chromosome')+
      ggplot2::scale_color_manual(values=color_map,guide='none')+
      ggplot2::scale_shape_identity() 
  )
}
