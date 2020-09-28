extract_which_chr <- function(data){
  unique(data$chr)[order(match(unique(data$chr), c(paste0("chr",1:22), "chrX", "chrY")))]
}

get_chrom_lengths=function(build=c('hg18','hg19','hg38')){
  build=match.arg(build)
  
  chrom_lengths_hg18=c(chr1=247249719,chr2=242951149,chr3=199501827,
                       chr4=191273063,chr5=180857866,chr6=170899992,
                       chr7=158821424,chr8=146274826,chr9=140273252,
                       chr10=135374737,chr11=134452384,chr12=132349534,
                       chr13=114142980,chr14=106368585,chr15=100338915,
                       chr16=88827254,chr17=78774742,chr18=76117153,
                       chr19=63811651,chr20=62435964,chr21=46944323,
                       chr22=49691432,chrX=154913754,chrY=57772954)
  
  
  chrom_lengths_hg19=c(chr1=249250621,chr2=243199373,chr3=198022430,
                       chr4=191154276,chr5=180915260,chr6=171115067,
                       chr7=159138663,chr8=146364022,chr9=141213431,
                       chr10=135534747,chr11=135006516,chr12=133851895,
                       chr13=115169878,chr14=107349540,chr15=102531392,
                       chr16=90354753,chr17=81195210,chr18=78077248,
                       chr19=59128983,chr20=63025520,chr21=48129895,
                       chr22=51304566,chrX=155270560,chrY=59373566)
  
  # https://genome.ucsc.edu/goldenpath/help/hg38.chrom.sizes
  chrom_lengths_hg38=c(chr1=248956422,chr2=242193529,chr3=198295559,
                       chr4=190214555,chr5=181538259,chr6=170805979,
                       chr7=159345973,chr8=145138636,chr9=138394717,
                       chr10=133797422,chr11=135086622,chr12=133275309,
                       chr13=114364328,chr14=107043718,chr15=101991189,
                       chr16=90338345,chr17=83257441,chr18=80373285,
                       chr19=58617616,chr20=64444167,chr21=46709983,
                       chr22=50818468,chrX=156040895,chrY=57227415)
  
  switch(build,
         hg18 = chrom_lengths_hg18, 
         hg19 = chrom_lengths_hg19, 
         hg38 = chrom_lengths_hg38
  )
}

get_total_length=function(chrom_lengths){
  sum(chrom_lengths)
}

get_cumulative_length=function(chrom_lengths){
  cumulative_length <- 0
  
  if(length(chrom_lengths) > 1){
    cumulative_length=head(c(0,cumsum(x = unname(chrom_lengths))), -1)
    names(cumulative_length)=names(chrom_lengths)
  }
  return(cumulative_length)
}

get_x_breaks=function(chrom_lengths){
  cumulative_length=get_cumulative_length(chrom_lengths)
  x_breaks=cumulative_length+round(chrom_lengths/2)
  names(x_breaks)=gsub('chr','',names(x_breaks))
  if(length(chrom_lengths) == 21){
    names(x_breaks)[20]=''
  }
  if(length(chrom_lengths) > 21){
    names(x_breaks)[20]=''
    names(x_breaks)[22]=''
  }
  return(x_breaks)
}

add_cumulative_pos=function(data, chrom_lengths){
  cumulative_length=get_cumulative_length(chrom_lengths)
  
  tmp <- Map(function(x,y){x$cumulative_pos <-  x$pos + y; return(x)}, 
             split(data, data$chr)[names(chrom_lengths)], get_cumulative_length(chrom_lengths)) 
  return(do.call(rbind, tmp))
}  

add_color=function(data,color1='black',color2='grey'){
  if ('color'%in%colnames(data)){
    user_color=data$color
  } else {
    user_color=rep(NA,nrow(data))
  }
  data$color=ifelse(data$chr %in% extract_which_chr(data)[c(TRUE, FALSE)],color1,color2)
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

