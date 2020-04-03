
setwd("/Users/lin.d/Desktop/FL_DC_paper/Barcoding/Barcode_data/Lowdose/")

DL188 = read.table("DL188_Counts.txt",header=TRUE)
DL188 = DL188[,!grepl("E", colnames(DL188))]
DL188 = DL188[,!grepl("B2", colnames(DL188))]
DL188 = DL188[,!grepl("ctl", colnames(DL188))]
DL188 = DL188[,!grepl("ID", colnames(DL188))]

DL208 = read.table("DL208_Counts.txt",header=TRUE)
DL208 = DL208[,!grepl("DL197", colnames(DL208))]
DL208 = DL208[,!grepl("Blank", colnames(DL208))]
DL208 = DL208[,!grepl("ID", colnames(DL208))]

dt = cbind(DL188,DL208)

# Due to low number barcode counts and low correlation between tech. replicates from these low dose exps, 
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
  bc = data.frame(mouse=m, total=bc.total, cDC1=bc.cDC1, cDC2=bc.cDC2, pDC=bc.pDC, B=bc.B, M=bc.M)
  bc.count = rbind(bc.count, bc)
}
write.table(bc.count, file=paste("low_dose_barcode_count.txt", sep="_"), sep="\t", row.names=FALSE, col.names=TRUE)

##### plot heatmap ##### 
library(pheatmap)
dt.pool = dt.pool[c("cDC1","cDC2","pDC","mye","B","eos","mon","neu",
                    "treatment","mouse","x","y")]
# specify what annotation to add per row
anno_row = dt.pool["treatment"]

# change colors for annotation
my_colour = list(
  treatment = c(PBS = "#09A2D1", FL = "#F2AA4CFF")
)

pheatmap(data.matrix(dt.pool[,c("cDC1","cDC2","pDC","mye","B")]),
         cluster_rows = T, cluster_cols = F, 
         annotation_row = anno_row,
         annotation_colors = my_colour,
         filename = "heatmap_2.pdf",
         width = 4,
         height = 4
)

