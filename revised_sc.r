################################################################################
sc1 = readRDS('/home/public/bbb/obj_integrated_anno.rds')
sc1$ct = as.character(sc1$class)
sc1$ct[which(sc1$ct %in% c('VEC_1','VEC_2'))] = 'VEC'
sc1$ct[which(sc1$ct %in% c('MEC_1','MEC_2'))] = 'MEC'
sc1$ct[which(sc1$ct %in% c('CarEC'))] = 'HEC'
sc1@meta.data = sc1@meta.data[c('orig.ident','nCount_RNA','nFeature_RNA','ct')]
st1 = sc1
###
sc2 = '/home/public/bbb/Zbe_6dpf/filtered_feature_bc_matrix/'
sc3 = '/home/public/bbb/Zbe_11dpf/filtered_feature_bc_matrix/'
#
sce2 = CreateSeuratObject(counts = Read10X(sc2),min.cells = 3,min.features = 100,assay = "RNA")
sce2$orig.ident = '6dpf_2'
sce3 = CreateSeuratObject(counts = Read10X(sc3),min.cells = 3,min.features = 100,assay = "RNA")
sce3$orig.ident = '11dpf_4'
#########
st = sce2
###
st$percent.mt <- PercentageFeatureSet(object = st, pattern = "^mt")
Idents(st) = 'ec'
st <- subset(st, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
#
st <- NormalizeData(st)
st <- FindVariableFeatures(st)
st <- ScaleData(st,features = rownames(st))
st <- RunPCA(st)
#
st <- FindNeighbors(st, reduction = "pca", dims = 1:30)
st <- FindClusters(st, resolution = 1, cluster.name = "pca_clusters")
st <- RunUMAP(st, reduction = "pca", dims = 1:30, reduction.name = "umap")
#
#DimPlot(st,reduction = "umap.harmony",group.by = c("harmony_clusters"))
#
DotPlot(st,features = c('kdrl','pecam1','cldn5a','slc7a5','cxcr4b','nrp1a','prox1a','cdh6','pcna','mki67','stab1','stab2','plvapb','rgcc','hand2','vcam1b'))
st$ct = as.character(st$seurat_clusters)
st$ct[which(st$ct %in% c(3,9,15,27,29))] = 'CapEC'
st$ct[which(st$ct %in% c(1))] = 'AEC'
st$ct[which(st$ct %in% c(5))] = 'LymEC'
st$ct[which(st$ct %in% c(7,11,14,20))] = 'MEC'
st$ct[which(st$ct %in% c(5,28))] = 'VEC'
st$ct[which(st$ct %in% c(2,13))] = 'AngEC'
st$ct[which(st$ct %in% c(12))] = 'HEC'
st = st[,which(st$ct %in% c('CapEC','AEC','LymEC','MEC','VEC','AngEC','HEC'))]
#
st2 = st
########
#########
st = sce3
###
st$percent.mt <- PercentageFeatureSet(object = st, pattern = "^mt")
Idents(st) = 'ec'
st <- subset(st, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
#
st <- NormalizeData(st)
st <- FindVariableFeatures(st)
st <- ScaleData(st,features = rownames(st))
st <- RunPCA(st)
#
st <- FindNeighbors(st, reduction = "pca", dims = 1:30)
st <- FindClusters(st, resolution = 1, cluster.name = "pca_clusters")
st <- RunUMAP(st, reduction = "pca", dims = 1:30, reduction.name = "umap")
#
#DimPlot(st,reduction = "umap.harmony",group.by = c("harmony_clusters"))
#
DotPlot(st,features = c('kdrl','pecam1','cldn5a','slc7a5','cxcr4b','nrp1a','prox1a','cdh6','pcna','mki67','stab1','stab2','plvapb','rgcc','hand2','vcam1b'))
st$ct = as.character(st$seurat_clusters)
st$ct[which(st$ct %in% c(3,7,13,20))] = 'CapEC'
st$ct[which(st$ct %in% c(5,6,19))] = 'AEC'
st$ct[which(st$ct %in% c(23))] = 'LymEC'
st$ct[which(st$ct %in% c(15,16))] = 'MEC'
st$ct[which(st$ct %in% c(9,11))] = 'VEC'
st$ct[which(st$ct %in% c(0,2))] = 'AngEC'
st$ct[which(st$ct %in% c(8))] = 'HEC'
st = st[,which(st$ct %in% c('CapEC','AEC','LymEC','MEC','VEC','AngEC','HEC'))]
#
st3 = st
########
rm(list = c('sc1','sce2','sce3','st'))
st1@meta.data = st1@meta.data[c('orig.ident','nCount_RNA','nFeature_RNA','ct')]
st2@meta.data = st2@meta.data[c('orig.ident','nCount_RNA','nFeature_RNA','ct')]
st3@meta.data = st3@meta.data[c('orig.ident','nCount_RNA','nFeature_RNA','ct')]
st1 = st1[,which(st1$ct %in% c('CapEC','AEC','LymEC','MEC','VEC','AngEC','HEC'))]
#
st = merge(st1,y = list(st2,st3),add.cell.ids = c('orig','new_6','new_11'))
st$dpf = substr(st$orig.ident,1,(nchar(st$orig.ident) - 2))

################################################################################
################################################################################
################################################################################
################################################################################

st = readRDS('/home/public/bbb/selected_endothelial_SeuratObject.rds')
DefaultAssay(st) = 'RNA'
st$ct = st$cell_type
st = st[,which(st$ct %in% c('CapEC','AEC','LymEC','MEC','VEC','AngEC'))]
Idents(st) = st$ct
#st = JoinLayers(st)
sm <- FindAllMarkers(object = st, only.pos = T, min.pct = 0.50)
#
save(sm,file = 'sm_gene_pct0.50.Rdata')


topn = sc.markers %>% 
  filter( avg_log2FC > 0 ) %>%
  group_by(cluster) %>% 
  top_n(25, avg_log2FC)

save(topn, file = 'top_25.Rdata')

################################################################################
################################################################################
################################################################################
################################################################################
#
#
#
#
#fig 3a
#zebrafish for correlation
#stt = readRDS('/home/public/bbb/selected_endothelial_SeuratObject.rds')
#DefaultAssay(stt) = 'RNA'
#stz = stt
#stz$ct1 = 'others'
#stz$ct1[which(stz$cell_type %in% c('CapEC'))] = 'CapEC'
#stz$ct1[which(stz$cell_type %in% c('AEC'))] = 'AEC'
#stz$ct1[which(stz$cell_type %in% c('VEC'))] = 'VEC'
#ae0 <- AverageExpression(object = stz,assays = "RNA",slot = "data",group.by = "ct1")
#ae0 <- as.data.frame(ae0)
#human
sth = readRDS('/home/ccwu/mouse_bbb/human_bbb/GSE256493_Fetal_CNS_sorted_endothelial_cells_seurat_object.rds')
stn = sth
stn$ct1 = 'others'
#stn$ct1[which(stn$ECclusters %in% c('Capillary'))] = 'CapEC'
stn$ct1[which(stn$ECclusters %in% c('Artery'))] = 'AEC'
#stn$ct1[which(stn$ECclusters %in% c('Vein'))] = 'VEC'
#ae1 <- AverageExpression(object = stn,assays = "RNA",slot = "data",group.by = "ct1")
#ae1 <- as.data.frame(ae1)
fm1 = FindMarkers(stn, ident.1 = 'AEC', group.by = 'ct1',logfc.threshold = 0,min.pct = 0,min.diff.pct = 0)
#Capillary Artery Vein
#mouse
load('/home/ccwu/mouse_bbb/mouse_bbb/mouse_e.Rdata')#st
stm = st
stm$ct1 = 'others'
#stm$ct1[which(stm$harmony_clusters == 18)] = 'CapEC'
stm$ct1[which(stm$harmony_clusters == 19)] = 'AEC'
#stm$ct1[which(stm$harmony_clusters == 16)] = 'VEC'
#ae2 <- AverageExpression(object = stm,assays = "RNA",slot = "data",group.by = "ct1")
#ae2 <- as.data.frame(ae2)
fm2 = FindMarkers(stm, ident.1 = 'AEC', group.by = 'ct1',logfc.threshold = 0,min.pct = 0,min.diff.pct = 0)
#
#
#
#
load('/home/ccwu/bcl6b/sc/analysis/homo_gene.Rdata')#h2m" "h2z" "m2h" "z2h" "z2m
#gene = read.table('fig_2_cell_marker.txt',head = T)
#load('heatmatp.specific.rdata') #avgcs
load('sm_gene_pct0.50.Rdata')
# avgcs = sm[which(sm$pct.2 <= 0.20),]
avgcs = sm#[which(sm$pct.1 >= 0.60),]
ec = unique(avgcs$cluster)
ge = avgcs
ge$ct = ge$cluster
for (i in 1:6){
  # temp = ge[which(ge$ct == ec[i]),]
  gene = ge[which(ge$ct == ec[i]),]
  gene = gene[1:35,]
  assign(paste0('gene_',i),gene)
}
# AEC   AngEC CapEC LymEC MEC   VEC  
#
i=1
###
gene = get(paste0('gene_',i))
colnames(gene)[7] = 'zebrafish'
gene_h = merge(gene,z2h,by = 'zebrafish',sort = F)
temp = as.data.frame(table(gene_h$zebrafish))
temp2 = temp[which(temp$Freq == 1),]
colnames(temp2)[1] = 'zebrafish'
temp3 = merge(gene_h,temp2,by = 'zebrafish',sort = F)
fe1 = unique(temp3$human)
fe0 = unique(temp3$zebrafish)
###
gene = get(paste0('gene_',i))
colnames(gene)[7] = 'zebrafish'
gene_m = merge(gene,z2m,by = 'zebrafish',sort = F)
temp = as.data.frame(table(gene_m$zebrafish))
temp2 = temp[which(temp$Freq == 1),]
colnames(temp2)[1] = 'zebrafish'
temp3 = merge(gene_m,temp2,by = 'zebrafish',sort = F)
fe2 = unique(temp3$mouse)
###
fe1 = fe1[-c(2,7,10)]
fe2 = fe2[-2]
#fe1 = fe1[-2]
#fe2 = fe2
re1 = fm1[fe1,]
re2 = fm2[fe2,]
r1 = length(which(sign(re1[,2]) == 1))/nrow(re0)
r2 = length(which(sign(re2[,2]) == 1))/nrow(re0)
#
p = DotPlot(stn,features = rev(fe1), group.by = c('ct1')) + coord_flip()+ theme(legend.position = 'none')
assign(paste0('p',i),p)
q = DotPlot(stm,features = rev(fe2), group.by = c('ct1')) + coord_flip()+ theme(legend.position = 'none')
assign(paste0('q',i),q)
#
library(patchwork)
p3+q3+p1+q1+p6+q6+plot_layout(nrow = 2)

x1 = fm1[fe1,]
x2 = fm2[fe2,]

write.csv(x1, file = 'aec_human.csv')
write.csv(x2, file = 'aec_mouse.csv')



