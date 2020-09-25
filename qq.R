# test this 
load("Y:/MIST+August/roman/gwas_data.rda")
qq_fast(gwas_data$pvalue, inflation_method = "reg")
qq_fast(gwas_data$pvalue, inflation_method = "med")
# or this
qq_fast(runif(nrow(gwas_data)))
qq_fast(runif(nrow(gwas_data)), inflation_method = "r")

# big data
big_gwas_data <-  do.call(rbind, replicate(15, gwas_data, simplify = FALSE)) 
qq_fast(big_gwas_data$pvalue)


# -----------------------------------------------------------------------
# QQ-plot function
# -----------------------------------------------------------------------
qq_fast <- function(pvalue,linecolor="deeppink",ci_color="steelblue",ci_alpha = 0.3,conf.alpha=.05,confint=T,log10=T,inflation_method="median"){

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

# function to calculate confidence interval
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
 if(confint)  p <- p + ggplot2::geom_polygon(data=get_conf(pvalue), ggplot2::aes(V1, V2), fill = ci_color, alpha=ci_alpha)
 
 p +
   ggplot2::geom_abline(slope = 1, intercept = 0, color = linecolor) + 
   scattermore::geom_scattermore(pointsize = 1.2) +      
   ggplot2::xlab(expression(Expected ~ ~-log[10](italic(p)))) +  
   ggplot2::ylab(expression(Observed ~ ~-log[10](italic(p)))) +
   ggplot2::ggtitle("QQ-plot", subtitle = paste0("lambda gc=",gif_out$es,gif_out$se)) +
   ggplot2::theme_bw()
 
  
}
