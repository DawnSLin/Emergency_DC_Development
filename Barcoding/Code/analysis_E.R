# set directory
setwd("~/Desktop/FL_DC_paper/Barcoding/Barcode_data/inc_E_analysis")

# reading in raw count data
DL191 = read.table("DL191_pool.txt",header=TRUE)
DL198 = read.table("DL198_pool.txt",header=TRUE)

dt=rbind(DL191,DL198)
# order columns
dt = dt[c("cDC1","cDC2","pDC","mye","B","E","eos","mon","neu",
          "treatment","exp","mouse")]

# load packages
library(Rtsne)
library(ggplot2)
library(grid)
library(reshape2)
library(ggdendro)
library(ggbeeswarm)
library(dplyr)
library(gridExtra)
library(gplots)
library(pheatmap)
library(tidyverse)

# colors
library(RColorBrewer)
colors = rev(brewer.pal(7,"RdYlBu"))
my_colour = c("#09A2D1","#F2AA4CFF")

######################################## tSNE analysis ##############################################

colnames(dt)
set.seed(100)
tsne=Rtsne(dt[,c(1:6)],perplexity=30,check_duplicates=FALSE) 

# color code on cell type
cDC1=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"cDC1"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
  labs(title="cDC1",x="",y="") 
cDC2=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"cDC2"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="cDC2",x="",y="")
pDC=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"pDC"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="pDC",x="",y="") 
eos=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"eos"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="eos",x="",y="") 
mon=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"mon"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="mon",x="",y="") 
neu=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"neu"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="neu",x="",y="") 
mye=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"mye"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="mye",x="",y="") 
B=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"B"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
    labs(title="B",x="",y="") 

E=
  ggplot(data = NULL, aes(x=tsne$Y[,1],y=tsne$Y[,2],col=dt[,"E"]))+
  geom_point(size=0.5)+
  scale_color_gradientn(colours=colors)+
  theme_classic()+
  labs(title="E",x="",y="") 

pdf("tsne.pdf", height = 9,width = 15)
multiplot(cDC1,cDC2,pDC,mye,B,E, cols=3)
dev.off()

# add t-SNE coordinate to dataset
D = data.frame(x = tsne$Y[,1], y = tsne$Y[,2])
dt = cbind(D,dt)

  ggplot(dt,aes(x=dt$x,y=dt$y,col=dt$treatment))+
  geom_point(aes(x,y,colour = factor(treatment)),alpha=0.5,size=1,shape=21)+
  scale_colour_manual(values = my_colour)+
  theme_classic()+
  ggtitle("treatment")

####################################################### DBSCAN #######################################################

library(fpc)
library(dbscan)
library(factoextra)

set.seed(100)
# eps: the radius of neighborhood around a point x
# MinPts: minimum number of neighbors within “eps” radius
# function to determine optimal eps value for the correponding k, ie MinPts #
dbscan::kNNdistplot(tsne$Y, k =  20)
abline(h = 4, lty = 2) # add a line

# perform dbscan clustering
db=fpc::dbscan(tsne$Y, eps=2.5, MinPts=10)
print(db)

# plot clusters
fviz_cluster(db, data=tsne$Y,stand = FALSE,
             ellipse = FALSE, 
             #ellipse.type="norm", ellipse.level=0.9, ellipse.alpha=0.1,
             show.clust.cent = FALSE,
             aes=c(tsne$Y[,1],tsne$Y[,2]),
             xlab = FALSE, ylab = FALSE,
             geom = "point", shape=16, pointsize=1, 
             ggtheme = theme_classic())

# extract cluster ID from db, add to dataset
dt$cluster=db$cluster
write.table(dt, "tsne_db_clusters.txt", sep="\t", row.names=FALSE, col.names=TRUE)

# count #bc in each cluster
D=data.frame()
cluster=unique(dt$cluster)
for (c in cluster)
{ 
  if(c==0)next;
  #c=13
  dtc=subset(dt, dt$cluster==c)
  # count #bc from each mouse per cluster
  mouse=unique(dtc$mouse)
  a=data.frame()
  for (m in mouse)
  {
    #m=1127
    counts=sum(dtc$mouse==m) 
    b=data.frame(mouse=m,counts=counts)
    b$cluster=c
    a=rbind(a,b)
  }
  D=rbind(D,a)
}
write.table(D, "cluster_mouse_bc.txt", sep="\t", row.names=FALSE, col.names=TRUE)

# calculate proportional output to cell type (based on read count) per cluster
D=data.frame()
cluster=unique(dt$cluster)
for (c in cluster)
{ 
  if(c==0)next;
  #c=23
  dtc=subset(dt, dt$cluster==c)
  #dtc=subset(dtc, dtc$treatment=="FL")
  
  cDC1=sum(dtc["cDC1"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  cDC2=sum(dtc["cDC2"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  pDC=sum(dtc["pDC"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  mye=sum(dtc["mye"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  B=sum(dtc["B"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  E=sum(dtc["E"])/sum(dtc[,c("cDC1","cDC2","pDC","mye","B","E")])
  a=data.frame(cluster=c,cDC1=cDC1, cDC2=cDC2, pDC=pDC,mye=mye, B=B,E=E)
  D=rbind(D,a)
}
write.table(D, "proportion per cluster.txt", sep="\t", row.names=FALSE, col.names=TRUE)


# order columns
dt$cluster <- sub("^", "cluster ", dt$cluster)

# order rows
D = data.frame()
cluster = c("cluster 1",
            "cluster 14",
            "cluster 19",
            "cluster 23",
            "cluster 21",
            "cluster 16",
            "cluster 24",
            "cluster 20",
            "cluster 10",
            "cluster 9",
            "cluster 11",
            "cluster 7",
            "cluster 8",
            "cluster 17",
            "cluster 18",
            "cluster 4",
            "cluster 15",
            "cluster 13",
            "cluster 6",
            "cluster 12",
            "cluster 22",
            "cluster 5",
            "cluster 3",
            "cluster 2")
for (c in cluster)
{ 
  #c="cluster 1"
  dtc=subset(dt, dt$cluster==c)
  dtc=dtc[order(rowSums(dtc[,c(1:5)]),decreasing=T),]
  D=rbind(D,dtc)
}

dt_order = D
dt_order = dt_order[c("cDC1","cDC2","pDC","mye","B","E","eos","mon","neu",
          "treatment","exp","mouse","cluster","x","y")]

######################## plot heatmap ########################

# specify what annotation to add per row
anno_row = dt_order[c("treatment","cluster")]

# change colors for annotation
my_colour = list(
  treatment = c(PBS = "#09A2D1", FL = "#F2AA4CFF")
)

# plot heatmap
pheatmap(data.matrix(dt_order[,c(1:6)]),
         cluster_rows = F, cluster_cols = F, 
         annotation_row = anno_row,
         annotation_colors = my_colour,
         filename = "heatmap_E.pdf",
         width = 4,
         height = 8
)






