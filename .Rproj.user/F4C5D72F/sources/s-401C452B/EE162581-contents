devtools::install_github("roman-tremmel/FASTGWASMAN")
library(FASTGWASMAN)
data("gwas_data")
head(gwas_data)

big_gwas_data <-  do.call(rbind, replicate(15, gwas_data, simplify = FALSE)) 
gwas_data$chr <- as.numeric(gsub("chr","",gwas_data$chr))
gwas_data[gwas_data$rsid=='rs602633','color']='green'
gwas_data[gwas_data$rsid=='rs1333045','color']='red'


FASTGWASMAN::fast_manhattan()



FASTGWASMAN::fast_manhattan(gwas_data[order(gwas_data$rsid),], build = "hg18", speed = "slow")
FASTGWASMAN::fast_manhattan(gwas_data, build = "hg18", speed = "fast")
FASTGWASMAN::fast_manhattan(gwas_data[1, ], build = "hg18", speed = "fast", pixels = c(512, 512))
FASTGWASMAN::fast_manhattan(gwas_data[,2:3 ], build = "hg18", speed = "fast", pixels = c(512, 512))
FASTGWASMAN::fast_manhattan(gwas_data[gwas_data$chr %in% c("chr1","chr9","chr21"),], build = "hg18", speed = "fast")

FASTGWASMAN::fast_manhattan(gwas_data, build = "hg18", speed = "fast", pixels = c(512, 512)*2, pointsize = 1.2)
FASTGWASMAN::manhattan(gwas_data, build='hg18', speed = "fast", color1 = "pink", color2 = "turquoise", pointsize = 3, pixels = c(1000, 500))
qqman::manhattan(gwas_data, bp = "pos", chr = "chrom", p = "pvalue")

FASTGWASMAN::fast_qq(gwas_data$pvalue)


FASTGWASMAN::manhattan(gwas_data, build = "hg18", speed = "fast", pixels = c(512, 512)*2)
ggsave("c:/Users/tremmel/Desktop/test.pdf")


# save(gwas_data, file = "data/gwas_data.rda")
devtools::document()
devtools::install(upgrade = "never")




# # code to show problem of lost points
# devtools::install_github("roman-tremmel/FASTGWASMAN")
# devtools::install_github('exaexa/scattermore', dependencies = F, force = T) # see the issue of the alpha values
# data("gwas_data")
# # this plot is fine
# FASTGWASMAN::manhattan(gwas_data, build = "hg18", speed = "fast", pixels = c(512, 512))
# # this plot loses its "skyscraper on x-axis value 9"
# FASTGWASMAN::manhattan(gwas_data, build = "hg18", speed = "fast", pixels = c(512, 512)*2)
# # here it is fine again
# FASTGWASMAN::manhattan(gwas_data, build = "hg18", speed = "fast",pointsize = 2, pixels = c(512, 512)*2)
# # here some points disappear e.g. look at y-values between 60-70
# FASTGWASMAN::manhattan(gwas_data, build = "hg18", speed = "fast",pointsize = 0.5, pixels = c(512, 512)*2)
# 
# 
