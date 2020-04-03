
setwd("/Users/lin.d/Desktop/FL_DC_paper/Barcoding/Barcode_data/IntraBM/")

DL197 = read.table("DL197_Counts.txt",header=TRUE)
DL197 = DL197[,!grepl("ctl", colnames(DL197))]
DL197 = DL197[,!grepl("ID", colnames(DL197))]

DL199 = read.table("DL199_Counts.txt",header=TRUE)
DL199 = DL199[,!grepl("ctl", colnames(DL199))]
DL199 = DL199[,!grepl("ID", colnames(DL199))]

dt = cbind(DL197,DL199)

# Quite a few samples in DL197 has low read counts and/or low correlation
# DL199 seems to have ok quality data
# Nevertheless, to combine both experiments
# barcode filtering is not applied before analysis

# Average between tech. replicates
Sample.names = unique(substr(names(dt),1,nchar(names(dt))-2))

d.avg=data.frame(row.names=row.names(dt))
new.names=vector()
for (s in Sample.names)
{
  e=as.data.frame(dt[,grep(s,names(dt))])
  if(ncol(e)<2)next;
  new.names=c(new.names,substr(names(e)[1],1,nchar(names(e)[1])-2))
  d.avg=cbind(d.avg,rowSums(e)/2)
}

names(d.avg)=new.names


# seperate data per mouse
pop=unlist(lapply(strsplit(names(d.avg),split ="_"),"[",4))
pop.names=vector()
pop.names=c(pop.names,unique(pop))
mouse=unique(unlist(lapply(strsplit(names(d.avg),split ="_"),"[",3)))

dt.pool = data.frame()
bc.count = data.frame()
for (m in mouse)
{ 
  #m=1101
  e=as.data.frame(d.avg[,grep(m,names(d.avg))])
  e=e[rowSums(e)>0,]
  t=unlist(lapply(strsplit(names(e),split ="_"),"[",2))  # extract treatment FL or PBS
  names(e)=pop.names
  e$mye=e$neu+e$mon+e$eos # add a extra column to pool reads from all mye samples
  e=sweep(e,2,colSums(e)/1e6,`/`) # normalize to 10^6
  e=asinh(e) # transformation
  e$mouse=m # add a extra column to record mouse ID
  e$treatment=t[1] # add a extra column to record treatment
  dt.pool=rbind(e,dt.pool)
  
  bc.total = nrow(e)
  bc.cDC1 = sum(e$cDC1>0)
  bc.cDC2 = sum(e$cDC2>0)
  bc.pDC = sum(e$pDC>0)
  bc.B = sum(e$B>0)
  bc.M = sum(e$mye>0)
  bc.E = sum(e$E>0)
  bc = data.frame(mouse=m, total=bc.total, cDC1=bc.cDC1, cDC2=bc.cDC2, pDC=bc.pDC, M=bc.M, B=bc.B, E=bc.E)
  bc.count = rbind(bc.count, bc)
}
write.table(bc.count, file=paste("IntraBM_barcode_count.txt", sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)



##### plot heatmap ##### 

library(pheatmap)

dt.pool = dt.pool[c("cDC1","cDC2","pDC","mye","B","E","eos","mon","neu",
                      "treatment","mouse","x","y")]
# specify what annotation to add per row
anno_row = dt.pool["treatment"]

# change colors for annotation
my_colour = list(
  treatment = c(PBS = "#09A2D1", FL = "#F2AA4CFF")
)

pheatmap(data.matrix(dt.pool[,c("cDC1","cDC2","pDC","mye","B","E")]),
         cluster_rows = T, cluster_cols = T, 
         annotation_row = anno_row,
         annotation_colors = my_colour,
         filename = "heatmap.pdf",
         width = 4,
         height = 5
)
