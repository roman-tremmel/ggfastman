library(devoid)
library(ggplot2)

bench_plot <- function(x){
  
  trash_dir <- tempfile()
  dir.create(trash_dir, recursive = TRUE)
  on.exit(unlink(trash_dir, recursive = TRUE))
  
  slow_plot <- function() {
    p <- FASTGWASMAN::manhattan(x, build = "hg18", speed = "slow")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".pdf"), plot = p)
    NULL
  }
  
  fast_plot <- function() {
    p <- FASTGWASMAN::manhattan(x, build = "hg18", speed = "fast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".pdf"), plot = p)
    NULL
  }
 
  ultrafast_plot <- function() {
    p <- FASTGWASMAN::manhattan(x, build = "hg18", speed = "ultrafast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".pdf"), plot = p)
    NULL
  }
  
  fastman_plot <- function(){
    pdf(tempfile(tmpdir = trash_dir, fileext =".pdf"))
    fastman::fastman(x, chr = "chr", ps = "pos", p = "pval")
    dev.off()
    NULL
  }
  
  qqman_plot <- function(){
    pdf(tempfile(tmpdir = trash_dir, fileext =".pdf"))
    qqman::manhattan(x, chr = "chr", bp = "pos", p = "pval", snp="rsid")
    dev.off()
    NULL
  }
  
  
  bench::mark(slow = slow_plot(), 
              fast = fast_plot(), 
              ultrafast = ultrafast_plot(), 
              fastman_plot = fastman_plot(),
              qqman_plot = qqman_plot(),
              min_iterations = 2) 
  
}


cad_gwas$y <- -log10(cad_gwas$pval)
cad_gwas$chr <- as.numeric(gsub("chr", "", cad_gwas$chrom))

res <- bench_plot(cad_gwas)
plot(res)


###
res <- bench::mark(slow = slow_plot(), fast = fast_plot(), min_iterations = 5)
benchplot(FASTGWASMAN::manhattan(big_cad_gwas, build = "hg18", speed = "slow"))
benchplot(FASTGWASMAN::manhattan(big_cad_gwas, build = "hg18", speed = "fast"))
