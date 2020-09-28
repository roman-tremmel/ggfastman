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
              min_iterations = 2) 
  
}

bench_qqplot <- function(x,y){
  
  trash_dir <- tempfile()
  dir.create(trash_dir, recursive = TRUE)
  on.exit(unlink(trash_dir, recursive = TRUE))
  
  slow_plot <- function() {
    p <- FASTGWASMAN::fast_qq(y, speed = "slow")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 100, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
  
  fast_plot <- function() {
    p <- FASTGWASMAN::fast_qq(y, speed = "fast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 100, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
  
  ultrafast_plot <- function() {
    p <- FASTGWASMAN::fast_qq(y, speed = "ultrafast")
    ggsave(tempfile(tmpdir = trash_dir, fileext =".png"), width = 100, height = 100, units = "mm", dpi = 300, plot = p)
    NULL
  }
  
  fastman_plot <- function(){
    png(tempfile(tmpdir = trash_dir, fileext =".png"), width = 100, height = 100,  res=300, units = "mm")
    fastman::fastqq(x, p = "pvalue")
    dev.off()
    NULL
  }
  
  qqman_plot <- function(){
    png(tempfile(tmpdir = trash_dir, fileext =".png"), width = 270, height = 100,  res=300, units = "mm")
    qqman::qq(y)
    dev.off()
    NULL
  }
  
  
  bench::mark(
              `FASTGWASMAN: slow` = slow_plot(),
              `FASTGWASMAN: fast` = fast_plot(),
              `FASTGWASMAN: ultrafast` = ultrafast_plot(),
              `fastman: fastman` = fastman_plot(),
              `qqman: manhattan` = qqman_plot(),
              min_iterations = 2) 
}
  


plot_bench <- function(x, funcs = c("FASTGWASMAN: slow", 
                                    "FASTGWASMAN: fast",
                                    "FASTGWASMAN: ultrafast", 
                                    "fastman: fastman",
                                    "qqman: manhattan")){
  require(ggsignif)
  p <- plot(x)
  p + scale_x_discrete("", limits = rev(funcs)) +
  ggsignif::geom_signif(comparisons = combn(funcs, 2, simplify = F), map_signif_level = T, step_increase = 0.2, color=1) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
}





