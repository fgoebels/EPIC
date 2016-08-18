#!/usr/local/bin/Rscript
library(gplots)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)


palette.breaks <- seq(0,1,0.1)
my_palette <- colorRampPalette(c("blue",  "red"))(10)

dataTable = read.table(args[1], row.names=1, header=T, check.names=F)
tmp = t(t(dataTable))

pdf(args[3])
heatmap.2(tmp, dendrogram="none", breaks = palette.breaks, col=my_palette, trace="none", density.info="none", Colv = "NA", Rowv = "NA", main = args[2])
dev.off()
