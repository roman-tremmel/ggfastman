library(bench)
library(ggplot2)
library(FASTGWASMAN)

bench_plot <- function(x){
  
  trash_dir <- tempfile()
  dir.create(trash_dir, recursive = TRUE)
  on.exit(unlink(trash_dir, recursive = TRUE))
  
  slow_plot <- function() {
    p <- FASTGWASMAN::fast_manhattan(x, build = "hg18", speed = "slow")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
  
  fast_plot <- function() {
    p <- FASTGWASMAN::fast_manhattan(x, build = "hg18", speed = "fast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
 
  ultrafast_plot <- function() {
    p <- FASTGWASMAN::fast_manhattan(x, build = "hg18", speed = "ultrafast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
  
  fastman_plot <- function(){
    png(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100,  res=300, units = "mm")
    fastman::fastman(x, chr = "chrom", ps = "pos", p = "pvalue")
    dev.off()
    NULL
  }
  
  qqman_plot <- function(){
    png(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100,  res=300, units = "mm")
    qqman::manhattan(x, chr = "chrom", bp = "pos", p = "pvalue", snp="rsid")
    dev.off()
    NULL
  }

  bench::mark(`FASTGWASMAN: slow` = slow_plot(), 
              `FASTGWASMAN: fast` = fast_plot(), 
              `FASTGWASMAN: ultrafast` = ultrafast_plot(), 
              `fastman: fastman` = fastman_plot(),
              `qqman: manhattan` = qqman_plot(),
              min_iterations = 10) 
  
}

plot_bench <- function(x){
  p <- plot(x)
  p + scale_x_discrete("Function", limits = rev(c("FASTGWASMAN: slow", 
                                      "FASTGWASMAN: fast",
                                      "FASTGWASMAN: ultrafast", 
                                      "fastman: fastman",
                                      "qqman: manhattan")))
}
  



