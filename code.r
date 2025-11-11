#sc = readRDS("./SeuratObject/obj_integrated_anno.rds") #original data
sc = readRDS('/home/public/bbb/selected_endothelial_SeuratObject.rds')
DefaultAssay(sc) = "RNA"

sc$Class = sc$class
sc$Class = gsub("VEC_1", "AngEC", sc$Class)
sc$Class = gsub("_[12]", "", sc$Class)
sc$Class = factor(sc$Class, unique(sc$Class))
sc = SetIdent(sc, value = "Class")
sc = sc[,!sc$Class %in% c("CarEC", "UnEC")]
DefaultAssay(sc)="integrated"
sc = RunTSNE(sc, dims = 1:18)
DimPlot(sc, group.by = "Class", reduction = "tsne")

sc = SetIdent(sc, value = "Class")
DefaultAssay(sc) = "RNA"
avgc = as.data.frame(AverageExpression(sc, assays = "RNA", group.by = "Class"))
colnames(avgc) <- gsub("RNA.", "", colnames(avgc))
amk = FindAllMarkers(sc, only.pos = T)
amks = amk[amk$p_val_adj < 0.01 & amk$pct.1 > .2,]
top50 = amks%>%group_by(cluster)%>%top_n(25, -p_val_adj)
top50 = top50%>%group_by(cluster)%>%top_n(25, avg_log2FC)
avgcs = avgc[top50$gene,]
pheatmap::pheatmap(t(avgcs), cluster_col = F,scale = "column",
                   cluster_row = F, border_color = F,
                   color = colorRampPalette(colors = c("#87c6fa","grey91","#fabb87"))(100))

amkd = amk[amk$pct.1 > 0.2 & amk$p_val_adj < 0.05, ]
g = amkd%>%group_by(gene)%>%reframe(gene_num = n())
g = g[g$gene_num >= 3,]
# amkd = amkd[amkd$gene %in% g$gene,]
# amkd = amkd[!amkd$gene %in% top50$gene,]
avgcd = avgc[g$gene,]
avgcd$rm = rowMeans(avgcd)
avgcd = avgcd[order(avgcd$rm, decreasing = F),]
avgcd = avgcd %>% top_n(150,rm)
avgcd[avgcd > 20] = 20
# avgcd[avgcd < 1] = 0
avgcd = avgcd[,1:6]
pheatmap::pheatmap(log2(t(avgcd) + 1), cluster_col = F,# scale = "row",
                   cluster_row =F, border_color = "grey",
                   color = colorRampPalette(colors = c("#87c6fa","grey91","#fabb87"))(100))
save(amk, amks, amkd, avgc, avgcs, avgcd, file = "./abcd.rda")


