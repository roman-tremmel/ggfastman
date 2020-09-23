#' Add together two numbers
#'
#' @param gwas A GWAS study. Must have chrom, pos, and a y column reflecting -log10 pvalues
#' @param build Genomic build. Currently supports hg18, hg19, and hg38
#' @param color1 Color for odd-numbered chromosomes
#' @param color2 Color for even-numbered chromosomes
#' @return Ta ggplot object that makes a Manhattan plot.
#' @export
#' @examples
#' data(cad_gwas)
#' cad_gwas$y=-log10(cad_gwas$pval)
#' manhattan(cad_gwas,build='hg18')

manhattan=function(gwas,build=c('hg18','hg19','hg38'), color1='black',color2='grey', y_scale = TRUE,
                   speed = c("slow", "fast", "ultrafast"),
                   pointsize=0, pixels=c(512, 512), ...){
  data=gwas

  if (!all(c('chrom','pos','y')%in%colnames(data))){
    stop('data must have chrom, pos and y columns')
  }

  build=match.arg(build)
  data=add_cumulative_pos(data,build)
  data=add_color(data,color1 = color1, color2 = color2)
  data=add_shape(data,shape=16)
  data=add_fill(data)
  chrom_lengths=get_chrom_lengths(build)
  xmax=get_total_length(chrom_lengths)
  x_breaks=get_x_breaks(chrom_lengths)

  color_map=unique(data$color)
  names(color_map)=unique(data$color)




  plot <- ggplot2::ggplot(data, ggplot2::aes(x=cumulative_pos,y=y,color=color,shape=shape,fill=fill))
  plot <- switch(speed,
         slow = plot + ggplot2::geom_point(),
         fast  = plot + scattermore::geom_scattermore(pointsize=pointsize, pixels=pixels),
         ultrafast = ggplot2::ggplot() + scattermore::geom_scattermost(data[, c("cumulative_pos", "y")])
         )

  if(y_scale) plot <- plot + ggplot2::scale_y_continuous(expand=c(0.01,0),name=expression('-log10(P-value)'))

  return(
    plot + ggplot2::theme_classic()+
      ggplot2::scale_x_continuous(limits=c(0,xmax),expand=c(0.01,0),breaks=x_breaks,
                                  labels=names(x_breaks),name='Chromosome')+
      ggplot2::scale_color_manual(values=color_map,guide='none')+
      ggplot2::scale_fill_manual(values = color_map, guide = 'none')+
      ggplot2::scale_shape_identity())
}
