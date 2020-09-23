get_cumulative_length=function(chrom_lengths){
  n_chrom=length(chrom_lengths)
  cumulative_length=c(0,cumsum(x = unname(chrom_lengths))[1:(n_chrom-1)])
  names(cumulative_length)=names(chrom_lengths)
  return(cumulative_length)
}

get_total_length=function(chrom_lengths){
  sum(chrom_lengths)
}

get_x_breaks=function(chrom_lengths){
  cumulative_length=get_cumulative_length(chrom_lengths)
  x_breaks=cumulative_length+round(chrom_lengths/2)
  names(x_breaks)=gsub('chr','',names(x_breaks))
  names(x_breaks)[20]=''
  names(x_breaks)[22]=''
  return(x_breaks)
}

add_cumulative_pos=function(data,build=c('hg18','hg19','hg38')){
  build=match.arg(build)

  if (!all(c('chrom','pos','y')%in%colnames(data))){
    stop('data must have chrom, pos and y columns')
  }

  chrom_lengths=get_chrom_lengths(build)
  cumulative_length=get_cumulative_length(chrom_lengths)
for (i in paste0('chr',seq(1,22))){
    data[data$chrom==i,'cumulative_pos']=data[data$chrom==i,'pos']+cumulative_length[i]
  }

  return(data)
}

add_color=function(data,color1='black',color2='grey'){
  if ('color'%in%colnames(data)){
    user_color=data$color
  } else {
    user_color=rep(NA,nrow(data))
  }

  data$color=ifelse(data$chrom%in%paste0('chr',seq(1,22,2)),color1,color2)
  data$color=ifelse(is.na(user_color),data$color,user_color)
  return(data)
}

add_shape = function(data,shape=16){
  if ('shape'%in% colnames(data)){
    user_shape = data$shape
  } else {
    user_shape = rep(NA, nrow(data))
  }

  data$shape = shape
  data$shape = ifelse(is.na(user_shape),data$shape,user_shape)
  return(data)
}

add_fill = function(data){
  if ('fill'%in%colnames(data)){
    NULL
  } else {
    data$fill = data$color
  }
  return(data)
}
