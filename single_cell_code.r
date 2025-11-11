rm(list = ls())
options(stringsAsFactors = F)
setwd("/home/xiaoyu_li/rproj/singleDrEndo/")
library(Seurat)
library(ggplot2)
endos = readRDS("/home/public/bbb/selected_endothelial_SeuratObject.rds")


DimPlot(endos,label = T)

DimPlot(endos, group.by = "samples")
DimPlot(endos, label = T, reduction = "umap", group.by = "cell_type")

endos <- SetIdent(endos, value = "cell_type")
DotPlot(endos,scale.by = "size",
        features = c(
          "mki67", "top2a","pclaf", "ccnb1",
          "syne1b","gpr182","cxcl12a","prox1a","cdh6",
          "mrc1a","lyve1b","dab2", "stab1", "stab2",  
          "cldn5a","slc2a1a","slc7a5","abcc2",
          "hey1", "hey2", "cxcr4a", "cxcr4b","efnb2a",
          "plvapa", "plvapb","rgcc", "apln","dll4"
        )
        )+
  scale_color_gradientn(colors = rev(viridis::viridis(10)))+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust = 1, vjust = 1))

amkss <- FindAllMarkers(endos, only.pos = T)
colnames(amkss)
topn <- amkss%>%group_by(cluster)%>%top_n(50, -p_val_adj)%>%top_n(25, avg_log2FC)
write.xlsx(topn, file = "./Supplementary Table S2-1.xlsx")


avgc = as.data.frame(AverageExpression(endos, assays = "RNA", group.by = "cell_type"))
colnames(avgc) = gsub("RNA.", "", colnames(avgc))

avgcs = avgc[rownames(avgc)%in%topn$gene,]
avgcs = avgcs[topn$gene,]
pheatmap::pheatmap(avgcs, scale = "row", cluster_col = F,
                   cluster_row = F, border_color = F,
                   color = colorRampPalette(colors = c("#87c6fa","white","#fabb87"))(100))

avgct = as.data.frame(AverageExpression(endos, assays = "RNA", group.by = c("cell_type","samples")))
colnames(avgct) = gsub("RNA.", "", colnames(avgct))
avgcts = avgct[rownames(avgct)%in%topn$gene,]
avgcts = avgcts[topn$gene,]
pheatmap::pheatmap(avgcts, scale = "row", cluster_col = F,
                   cluster_row = F, border_color = F,
                   color = colorRampPalette(colors = c("#87c6fa","white","#fabb87"))(100))



###################################################################################
###################################################################################
###################################################################################

abc = read.table("./reanalysisOutput/abc family.txt", header = F)
abc = as.vector(abc$V1)
abc = rownames(scs)[rownames(scs) %in% abc]
length(abc)
slc = read.table("./reanalysisOutput/slc family.txt", header = F)
slc = as.vector(slc$V1)
slc = rownames(scs)[rownames(scs) %in% slc]
tj = read.table("./reanalysisOutput/tjz.txt", header = F)
tj = as.vector(tj$V1)
tj =  rownames(scs)[rownames(scs) %in% tj]
scs = endos
geneSet = list(
 "ABC Transporter" = abc,
  "SLC Transporter" = slc,
  "Tight Junction" = tj
)
scs = AddModuleScore(scs, geneSet, name = "BBB_TERMS_")
colnames(scs@meta.data)[grep("^BBB_TERMS_", colnames(scs@meta.data))] = names(geneSet)
data = scs@meta.data
colnames(data)
dat = data[,c(1,6,7:ncol(data))]
dat = reshape2::melt(dat)
head(dat)

dat$cell_type = factor(dat$cell_type, levels = c("CapEC", "AEC","AngEC", "VEC", "LymEC", "MEC"))
i = 1
ggplot(data = dat[dat$variable == unique(dat$variable)[i],])+labs(x = NULL, y = "Module Score")+
  geom_boxplot(aes(x = cell_type, y = value), outliers = F,linewidth = 0.35, fill = "skyblue")+
  scale_y_continuous(expand = c(0,0.01,0,0.1))+
  ggtitle(label = unique(dat$variable)[i])+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))



###################################################################################
###################################################################################
###################################################################################


mkcap <- FindMarkers(endos, ident.1 = "CapEC",logfc.threshold = 0.1,min.pct = 0.01)
mkcap$gene = rownames(mkcap)
mkcaps <- mkcap[mkcap$avg_log2FC > 0.2 & mkcap$p_val_adj < 0.05,]
nrow(mkcaps)
abc = read.table("./reanalysisOutput/abc family.txt", header = F)
abc = as.vector(abc$V1)
abc = rownames(mkcaps)[rownames(mkcaps) %in% abc]
length(abc)
slc = read.table("./reanalysisOutput/slc family.txt", header = F)
slc = as.vector(slc$V1)
slc = rownames(mkcaps)[rownames(mkcaps) %in% slc]
length(slc)

tj = read.table("./reanalysisOutput/tjz.txt", header = F)
tj = as.vector(tj$V1)
tj =  rownames(mkcaps)[rownames(mkcaps) %in% tj]

geneSet = list(
  "ABC Transporter" = abc,
  "SLC Transporter" = slc,
  "Tight Junction" = tj
)
scs = endos
scs = AddModuleScore(scs, geneSet, name = "BBB_TERMS_")
colnames(scs@meta.data)[grep("^BBB_TERMS_", colnames(scs@meta.data))] = names(geneSet)
data = scs@meta.data
colnames(data)
dat = data[,c(1,6,7:ncol(data))]
dat = reshape2::melt(dat)
head(dat)

dats = dat[dat$cell_type %in% c('CapEC'),]
i = 3 # i = c(1,2,3)
ggplot(data = dats[dats$variable==unique(dats$variable)[i],])+
  labs(x = NULL, y = "Module Score")+
  geom_boxplot(aes(x = samples, y = value), outliers = F,linewidth = 0.35, fill = "#e7ba90")+
  scale_y_continuous(expand = c(0,0.02,0,0.5))+
  ggtitle(label = unique(dats$variable)[i])+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5))


###########################################################################################
###########################################################################################
###########################################################################################

load("/home/public/bbb/mouse_e.Rdata")
sts = st
abc = read.table("./reanalysisOutput/Abc family.txt", header = F)
abc = as.vector(abc$V1)
abc = rownames(st)[rownames(st) %in% abc]
length(abc)
slc = read.table("./reanalysisOutput/Slc family.txt", header = F)
slc = as.vector(slc$V1)
slc = rownames(st)[rownames(st) %in% slc]
tj = read.table("./reanalysisOutput/Tjp.txt", header = F)
tj = as.vector(tj$V1)
tj =  rownames(st)[rownames(st) %in% tj]
geneSet = list(
  "ABC Transporter" = abc,
  "SLC Transporter" = slc,
  "Tight Junction" = tj
)
st = sts
st = AddModuleScore(st, features = geneSet, name = "BBB_TERMS_")
colnames(st@meta.data)[grep("^BBB_TERMS_", colnames(st@meta.data))] = names(geneSet)

st$orig.ident.2 = gsub(".[12]$", "", st$orig.ident)

data = st@meta.data
colnames(data)
dat = data[,c("orig.ident.2", "harmony_clusters", names(geneSet))]
dat = reshape2::melt(dat)
head(dat)

dats = dat[dat$harmony_clusters == 18,]
i = 1
ggplot(data = dats[dats$variable == unique(dats$variable)[i],])+labs(fill = "samples", x = NULL, y = "Score")+
  geom_boxplot(aes(x = orig.ident.2, y = value), outliers = F,linewidth = 0.35, fill = "#e7ba90")+
  ggtitle(label = unique(dats$variable)[i])+
  scale_y_continuous(expand = c(0,0.01,0,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


###########################################################################################
###########################################################################################
###########################################################################################

sth = readRDS("/home/public/bbb/GSE256493_Fetal_CNS_sorted_endothelial_cells_seurat_object.rds")
age_2 = data.frame(
  age_1 = unique(sth$age),
  age_2 = c( "GW14.4-16.4", "GW15", "GW18","GW9")
)
sth$age_2 = ""
sth$age_2 = age_2[match(sth$age,age_2$age_1),]$age_2
sth$age_2 = factor(sth$age_2, levels = c("GW9","GW14.4-16.4", "GW15", "GW18"))
sth = SetIdent(sth , value = "ECclusters")
mkcap = FindMarkers(sth, ident.1 = "Capillary")
mkcaps = mkcap[mkcap$avg_log2FC > 0 & mkcap$pct.1 > 0.1 & mkcap$p_val_adj < 0.05,]

abc = read.table("./reanalysisOutput/ABC family.txt", header = F)
abc = as.vector(abc$V1)
abc = rownames(sth)[rownames(sth) %in% abc]
length(abc)
slc = read.table("./reanalysisOutput/SLC family.txt", header = F)
slc = as.vector(slc$V1)
slc = rownames(sth)[rownames(sth) %in% slc]
tj = read.table("./reanalysisOutput/TJ.txt", header = F)
tj = as.vector(tj$V1)
tj =  rownames(sth)[rownames(sth) %in% tj]
geneSet = list(
  "ABC Transporter" = abc,
  "SLC Transporter" = slc,
  "Tight Junction" = tj
)

sths = sth
sths = AddModuleScore(sths, features = geneSet, name = "BBB_TERMS_")
colnames(sths@meta.data)[grep("^BBB_TERMS_", colnames(sths@meta.data))] = names(geneSet)

data = sths@meta.data
colnames(data)
dat = data[,c("age_2", "ECclusters", names(geneSet))]
dat = reshape2::melt(dat)
head(dat)
dats = dat[dat$ECclusters == "Capillary",]
head(dats)
i = 3
ggplot(data = dats[dats$variable == unique(dats$variable)[i],])+labs(fill = "samples", x = NULL, y = "Score")+
  geom_boxplot(aes(x = age_2, y = value), outliers = F,linewidth = 0.35, fill = "#e7ba90")+
  ggtitle(label = unique(dats$variable)[i])+
  scale_y_continuous(expand = c(0,0.01,0,0.1))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

