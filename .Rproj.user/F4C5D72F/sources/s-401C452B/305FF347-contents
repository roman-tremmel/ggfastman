#' A fast version of a quantile-quantile plot
#'
#' This is a very fast plotting function for quantile-quantile (QQ) plots using \code{\link{ggplot2}} and \code{\link{scattermore}} functionalities.
#'
#' @param pvalue A numeric vector of pvalues. 
#' @param speed The speed option. Fast and ultrafast use scattermore functionality
#' @param pointsize Only when using the 'fast' option you can increase pointsize. Default is 1.2
#' @param linecolor Color of the abline of the theoretical and sample distributions.
#' @param ci_color Color fof the confidence interval
#' @param color2 Fill color for even-numbered chromosomes
#' @param ci_alpha Alpha value (e.g. transparency) of the confidence's fill color
#' @param conf.alpha Alpha level of the confidence interval. Default value = 0.05
#' @param confint Print a confident interval. Default = TRUE
#' @param log10 Should the pvalues -log10 tranformed? Default = TRUE
#' @param inflation_method Method to calculate the genomic inflation factor. choices = c("median","regression"). The "median" method is the most common way,while the regression method to estimate lambda is adapted from the GenABEL package. 
#' @param title A title added to the plot 
#' @return A ggplot2 object/plot
#' @export
#' @examples
#' data("gwas_data")
#' head(gwas_data)
#' qq_fast(gwas_data$pvalue, inflation_method = "reg")
#' qq_fast(gwas_data$pvalue, inflation_method = "med")
fast_qq <- function(pvalue,speed ="fast",pointsize = 1.2,linecolor="deeppink",ci_color="steelblue",ci_alpha = 0.3,conf.alpha=.05,
                    confint=T,log10=T,inflation_method="median",title="QQ-plot"){
  speed=match.arg(speed,choices = c("slow", "fast", "ultrafast"))
  
  if (!is.numeric(data$pvalue)) stop("pvalue column should be numeric.")
  
  inflation_method <- match.arg(inflation_method, choices = c("median","regression"))
  # calculate expected and transform observed  
  qq_dat <- data.frame(pvalues = pvalue,
                       exp.pvalues = (rank(pvalue, ties.method="first")-.5)/(length(pvalue)+1))
  
  gif_out <- list()
  gif_out$es <- round(median(qchisq(pvalue, 1, lower.tail=FALSE))/qchisq(0.5,1), 2)
  gif_out$se <- ""
  
  if(inflation_method == "regression"){
    qq_dat_chi <- data.frame(pvalues = qchisq(pvalue, 1, lower.tail=FALSE))
    qq_dat_chi$exp.pvalues = (rank(qq_dat_chi$pvalue, ties.method="first")-.5)/(length(qq_dat_chi$pvalue)+1)
    qq_dat_chi$ppoints_chi_rank = qchisq(qq_dat_chi$exp.pvalues, df=1, lower.tail=FALSE)
    res <- round(summary(lm(sort(pvalues)~0+sort(ppoints_chi_rank), data = qq_dat_chi))$coeff,2)
    
    gif_out$es <- res[1, 1]
    gif_out$se <- paste(", se=",res[1, 2])
  }
  
  
  if(log10){
    qq_dat <- -log10(qq_dat)
  }
  
  # calculate confidence interval
  get_conf <- function(pvalue){
    conf.points=1000
    n <- length(pvalue)+1
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))}
    return(as.data.frame(mpts))
  }
  
  p <- ggplot2::ggplot(qq_dat, ggplot2::aes(exp.pvalues, pvalues))
  if(confint)  p <- p + ggplot2::geom_polygon(data=get_conf(pvalue), ggplot2::aes(V1, V2), fill = ci_color, alpha=ci_alpha) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = linecolor) + 
    ggplot2::xlab(expression(Expected ~ ~-log[10](italic(p)))) +  
    ggplot2::ylab(expression(Observed ~ ~-log[10](italic(p)))) +
    ggplot2::ggtitle(title, subtitle = substitute(paste(lambda~gc,"=", v1, v2), list(v1=gif_out$es, v2=gif_out$se))) +
    ggplot2::theme_bw()
  # plot
  switch (speed,
    slow = p + ggplot2::geom_point(),
    fast = p + scattermore::geom_scattermore(pointsize = pointsize),
    ultrafast = p + scattermore::geom_scattermost(qq_dat[,c("exp.pvalues", "pvalues")])
  )
  
}
