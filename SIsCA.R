library("vegan")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("ggfortify")
# vector with color scheme used in the figures, replace as required
# check out https://coolors.co/ for nice color palettes
#                 1         2       3         4         5         6         7         8         9         0
#               red       orange   green     lila       blue      cyan      brown   redgrey    grey    darkgreen
colorMT <- c("#FF0000","#FF9900","#99CC00","#6700CE","#0073FE","#00CC99","#996633","#B19B9B","#7F7F7F","#4C6700",
             "#FFCCCC","#FFEBCC","#F0FFC2","#E1C2FF","#CCE3FF","#C2FFF0","#F0E0D1","#EFEBEB","#F2F2F2","#EAFFAE",
             "#FF9999","#FFD699","#E0FF85","#C285FF","#99C7FF","#85FFE0","#E0C2A3","#E0D7D7","#D8D8D8","#D4FF5C",
             "#FF6666","#FFC266","#D1FF47","#A449FF","#65ABFF","#47FFD1","#D1A375","#D0C3C3","#BFBFBF","#BFFF0B",
             "#BF0000","#BF7300","#739900","#4D009B","#0056BE","#009973","#734D26","#8C6D6D","#2F2F2F","#394D00",
             "#7F0000","#7F4D00","#4D6600","#340067","#003A7F","#00664D","#4D3319","#5D4949","#000000","#263400")
# 1-10 normal color, 11-20 very bright, 21-30 bright, 31-40 slightly bright, 41-50 dark, 51-60 very dark

# Section 1 - data import (required)
# import protein-SIP data (see example data "SIsCA_example.csv" for reference)
SIsCA.raw <- read.csv(file='SIsCA_example.csv')
rownames(SIsCA.raw) <- sapply(1:nrow(SIsCA.raw),function(x) {paste(SIsCA.raw[x,1],"_",SIsCA.raw[x,4],SIsCA.raw[x,5],"_",x,sep="")})

# make vectors of peptide sequences, MAGs and/or taxonomy
SIsCA.MAG <- levels(SIsCA.raw$MAG)
SIsCA.tax <- levels(SIsCA.raw$taxonomy)
SIsCA.pept <- levels(SIsCA.raw$sequence)

# Section 2 - quality control (not required)
# visualize incorporation as heatmap (one heatmap per peptide sequence)
# (useful to find outliers/duplicates)
sapply(SIsCA.pept,function(x) {heatmap(as.matrix(SIsCA.raw[which(SIsCA.raw$sequence==x),11:31]), Colv=NA, margins=c(4,14))})

# compare incorporation in all peptides of a MAG (change to SIsCA.tax and SIsCA.raw$taxonomy for taxonomy-based comparison)
# (useful to find outliers)
# select MAG
x <- 3
# create one heatmap per peptide, numbers over plots indicate line numbers in SIsCA.raw
sapply(levels(droplevels(SIsCA.raw[which(SIsCA.raw$MAG==SIsCA.MAG[x]),]$sequence)), function(p) {
  heatmap(t(sapply(c("C","B","A"), function(r) {sapply(1:3, function(t) {as.numeric(SIsCA.raw[which(SIsCA.raw$sequence==p&SIsCA.raw$replicate==r&SIsCA.raw$time.point==t),11:31])})}))
  , Colv=NA, Rowv=NA, main=p)
  text(0.2,0.9,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="A"&SIsCA.raw$time.point==1), collapse=" "))
  text(0.45,0.9,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="A"&SIsCA.raw$time.point==2), collapse=" "))
  text(0.7,0.9,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="A"&SIsCA.raw$time.point==3), collapse=" "))
  text(0.2,0.53,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="B"&SIsCA.raw$time.point==1), collapse=" "))
  text(0.45,0.53,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="B"&SIsCA.raw$time.point==2), collapse=" "))
  text(0.7,0.53,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="B"&SIsCA.raw$time.point==3), collapse=" "))
  text(0.2,0.16,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="C"&SIsCA.raw$time.point==1), collapse=" "))
  text(0.45,0.16,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="C"&SIsCA.raw$time.point==2), collapse=" "))
  text(0.7,0.16,paste(which(SIsCA.raw$sequence==p&SIsCA.raw$replicate=="C"&SIsCA.raw$time.point==3), collapse=" "))
  })
SIsCA.MAG[x]

# Section 3 - creating time series (required)
# merging of replicates and time points
SIsCA.concat <- sapply(SIsCA.pept,function(p) {sapply(1:3,function(t) {apply(SIsCA.raw[which(SIsCA.raw$sequence==p&SIsCA.raw$time.point==t),11:31],2,mean)})})
# for color-coding of MAGs
SIsCA.MAG.col <- merge(distinct(SIsCA.raw[,c(1,6)]),cbind(MAG=SIsCA.MAG,color=colorMT[c(1:4)]), by="MAG")
# create heatmap where each row is one peptide, columns represent R² values, grey bars on top show time points
heatmap(t(SIsCA.concat), 
        Colv=NA, 
        labRow=as.vector(SIsCA.MAG.col[order(SIsCA.MAG.col$sequence),]$MAG),
        cexRow=1.0, 
        labCol="",
        ColSideColors = c(rep("#AAAAAA",21),rep("#888888",21),rep("#555555",21)), 
        RowSideColors = as.vector(SIsCA.MAG.col[order(SIsCA.MAG.col$sequence),]$color))


# Section 4 - PCA of peptides (not required)
# replace missing values with 0
SIsCA.concat.0 <- SIsCA.concat
SIsCA.concat.0[is.na(SIsCA.concat.0)] <- 0
# PCA
SIsCA.concat.pca <- rda(t(SIsCA.concat.0))
# calculate the % variances
SIsCA.concat.var <- sprintf("%.1f %%",SIsCA.concat.pca$CA$eig/SIsCA.concat.pca$tot.chi*100)  
# plot
ggplot(scores(SIsCA.concat.pca)$sites, aes(x=PC1, y=PC2, label=as.vector(SIsCA.MAG.col[order(SIsCA.MAG.col$sequence),]$MAG))) +
  geom_point(color=as.vector(SIsCA.MAG.col[order(SIsCA.MAG.col$sequence),]$color),size=3) +
  geom_text_repel(size=3)+
  xlab(bquote("PC1 ("*.(SIsCA.concat.var[1])*" variance)"))+
  ylab(bquote("PC2 ("*.(SIsCA.concat.var[2])*" variance)"))+
  theme(panel.background=element_rect(fill="#FFFFFF"),axis.line=element_line(size=0.5, linetype="solid", colour="#000000"), axis.text=element_text(size=12, colour="#000000"), axis.title=element_text(size=12))

# Section 5 - merging peptides of each MAG (required)
# average all peptides of a MAG
SIsCA.MAG.concat <- sapply(SIsCA.MAG,function(m) {sapply(1:3,function(t) {apply(SIsCA.raw[which(SIsCA.raw$MAG==m&SIsCA.raw$time.point==t),11:31],2,mean)})})
# heatmap where each row represents one MAG, columns represent R² values, grey bars on top show time points
heatmap(t(SIsCA.MAG.concat), 
        Colv=NA, 
        labRow=SIsCA.MAG,
        cexRow=1,
        labCol="",
        ColSideColors = c(rep("#AAAAAA",21),rep("#888888",21),rep("#555555",21)), 
        RowSideColors = colorMT[c(1:4)])

# Section 6 - PCA of MAGs (required)
# replace missinv values with 0
SIsCA.MAG.concat.0 <- SIsCA.MAG.concat
SIsCA.MAG.concat.0[is.na(SIsCA.MAG.concat.0)] <- 0
# PCA
SIsCA.MAG.concat.pca <- rda(t(SIsCA.MAG.concat.0))
# calculate the % variances 
SIsCA.MAG.concat.var <- sprintf("%.1f %%",SIsCA.MAG.concat.pca$CA$eig/SIsCA.MAG.concat.pca$tot.chi*100) 
# plot
ggplot(scores(SIsCA.MAG.concat.pca)$sites, aes(x=PC1, y=PC2, label=SIsCA.MAG)) +
  geom_point(color=colorMT[c(1:4)],size=3) +
  geom_text_repel(size=3)+
  xlab(bquote("PC1 ("*.(SIsCA.MAG.concat.var[1])*" variance)"))+
  ylab(bquote("PC2 ("*.(SIsCA.MAG.concat.var[2])*" variance)"))+
  theme(panel.background=element_rect(fill="#FFFFFF"),axis.line=element_line(size=0.5, linetype="solid", colour="#000000"), axis.text=element_text(size=12, colour="#000000"), axis.title=element_text(size=12))

# Section 7 - MAG heatmap (not required)
sapply(SIsCA.MAG, function(MAG) {heatmap(t(cbind(T3=SIsCA.MAG.concat[43:63,MAG],T2=SIsCA.MAG.concat[22:42,MAG],T1=SIsCA.MAG.concat[1:21,MAG])), Colv=NA, Rowv=NA, labCol=NA, main=MAG)})
