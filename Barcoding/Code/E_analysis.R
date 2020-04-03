# set directory
setwd("~/Desktop/FL_DC_paper/Barcoding/Barcode_data/inc_E_analysis")

# reading in raw count data
DL191 = read.table("DL191_pool.txt",header=TRUE)
DL198 = read.table("DL198_pool.txt",header=TRUE)

dt=rbind(DL191,DL198)

# colors
library(RColorBrewer)
colors = rev(brewer.pal(7,"RdYlBu"))

######################## plot heatmap ########################
library(pheatmap)

# specify what annotation to add per row
anno_row = dt["treatment"]

# change colors for annotation
my_colour = list(
  treatment = c(PBS = "#09A2D1", FL = "#F2AA4CFF")
)

# plot heatmap
pheatmap(data.matrix(dt[,c("cDC1","cDC2","pDC","mye","B","E")]),
         cluster_rows = T, cluster_cols = T, 
         annotation_row = anno_row,
         annotation_colors = my_colour,
         filename = "heatmap_E_2.pdf",
         width = 4,
         height = 5
)