argv = commandArgs(trailingOnly=TRUE)

p.plot_grid <- function(input_file, output_path, ...){
  
  t <- read.csv(input_file,sep=' ',header=FALSE)
  colnames(t) = c('pos','hapMap','IBD', 'color')
  t$pos # "genomic position" "hapMap recomb rate" "IBDRecomb recomb rate" 
  library("ggpubr")
  ggscatter(t, x = "hapMap", y = "IBD", 
            add = "reg.line", #conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "deCode", ylab = "FastRecomb", color = "color", palette = c(blue = "blue", red = "red", black = "black"), alpha = 0.5)
  #dev.off()
  ggsave(output_path)
  
  
}

input_file = argv[1] 
p.plot_grid(input_file, paste(argv[1], '.png', sep = '')) 


