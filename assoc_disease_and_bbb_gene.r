#
setwd("/home/ccwu/mouse_bbb/gwas_nealelab/bbb_gene")
a  = read.table('hg19.tss.bed')
ge = c('SLC16A1','SLC3A2','CLDN5','SLC7A5','SLC2A1','ABCC2','ZIC2','LEF1','SLCO1C1')
#
di = c('AD','DM','EP','MS','PD','I60','I61','I62','I63','I64','I65','I67')
fin = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(di)){
	#
	temp = read.table(paste0('../tmp',i,'.tsv'))
	#
	for (j in 1:length(ge)){
		#
		gen = ge[j]
		aa = unique(a[which(a$V5 == gen),c(1,2)])
		pos = aa[1,2]
		chr = aa[1,1]
		chr = substr(chr,4,nchar(chr))
		pos1 = pos - 1000000
		pos2 = pos + 1000000
		#
		tem = na.omit(temp[which(temp$V1 == chr),])
		num1 = which(tem$V2 > pos1)
		num2 = which(tem$V2 < pos2)
		num3 = intersect(num1,num2)
		te = tem[num3,]
		tt = min(te$V15)
		ds = abs(te$V2[which(te$V15 == tt)] - pos)
		mat = as.data.frame(array(,dim = c(1,4)))
		mat[1,1] = di[i]
		mat[1,2] = gen
		mat[1,3] = tt
		mat[1,4] = ds
		fin = rbind(fin,mat)
	}
}
save(fin, file = 'fin_v1.Rdata')
fin$V5 = fin$V1
fin$V5[which(fin$V5 == 'AD')] = '01'
fin$V5[which(fin$V5 == 'DM')] = '02'
fin$V5[which(fin$V5 == 'EP')] = '03'
fin$V5[which(fin$V5 == 'MS')] = '04'
fin$V5[which(fin$V5 == 'PD')] = '05'
fin$V5[which(fin$V5 == 'I60')] = '06'
fin$V5[which(fin$V5 == 'I61')] = '07'
fin$V5[which(fin$V5 == 'I62')] = '08'
fin$V5[which(fin$V5 == 'I63')] = '09'
fin$V5[which(fin$V5 == 'I64')] = '10'
fin$V5[which(fin$V5 == 'I65')] = '11'
fin$V5[which(fin$V5 == 'I67')] = '12'
fin$V6 = fin$V2
fin$V6[which(fin$V6 == 'SLC16A1')] = 9
fin$V6[which(fin$V6 == 'SLC3A2')] = 8
fin$V6[which(fin$V6 == 'CLDN5')] = 7
fin$V6[which(fin$V6 == 'SLC7A5')] = 6
fin$V6[which(fin$V6 == 'SLC2A1')] = 5
fin$V6[which(fin$V6 == 'ABCC2')] = 4
fin$V6[which(fin$V6 == 'ZIC2')] = 3
fin$V6[which(fin$V6 == 'LEF1')] = 2
fin$V6[which(fin$V6 == 'SLCO1C1')] = 1
fin$V7 = 'a'
fin$V7[which(fin$V4 < 1000000)] = 'a'
fin$V7[which(fin$V4 < 100000)] = 'b'
fin$V7[which(fin$V4 < 10000)] = 'c'
fin$V7[which(fin$V7 == 'a')] = 1
fin$V7[which(fin$V7 == 'b')] = 2
fin$V7[which(fin$V7 == 'c')] = 3
# fin$V7[which(fin$V7 < 1000)] = 1
fin$V3[which(fin$V3 < 1e-40)] = 1e-40
library(RColorBrewer)
# cl = brewer.pal(11,'Spectral')
# cl2 = rev(cl)
ggplot(fin,aes(x = V5,y = V6)) + geom_point(aes(color = -log(V3,10),size = V7)) + scale_color_gradientn(values = seq(0,1,0.01),colors = colorRampPalette(colors = c("#F57D15","#9F2A63","#280B54"))(201)) #+ scale_color_gradientn(values = seq(0,1,0.01),colors = colorRampPalette(colors = cl2[6:1])(201))




########################################################################################################################################################################################
hs1 = readRDS('/home/ccwu/monkey_public/GSE256493_Adult_control_brain_temporal_lobe_unsorted_endothelial_and_perivascular_cells_seurat_object.rds')
hs2 = readRDS('/home/ccwu/monkey_public/GSE256493_Adult_control_brain_temporal_lobe_sorted_endothelial_cells_seurat_object.rds')
hs1$clusters = hs1$Cellclusters
hs1@meta.data <- hs1@meta.data[,c('clusters','Patient','sex','age')]
hs2$clusters = hs2$ECclusters
hs2@meta.data <- hs2@meta.data[,c('clusters','Patient','sex','age')]
hs = merge(hs1,y = hs2)
sm <- FindMarkers(object = hs, ident.1 = 'Capillary', group.by = 'clusters', assay = 'RNA',slot = 'counts',logfc.threshold = 0,min.pct = 0)
#
ge1 = rownames(sm)[which(sm$pct.1 > 0.5 & sm$pct.2 < 0.5)]
ge2 = rownames(sm)[which(sm$pct.1 < 0.05 & sm$pct.2 > 0.1)]
save(ge1,ge2,file = '/home/ccwu/mouse_bbb/gwas_nealelab/bbb_gene/EC_gene.Rdata')
#
setwd("/home/ccwu/mouse_bbb/gwas_nealelab/bbb_gene")
a  = read.table('hg19.tss.bed')
a  = a[which(a$V1 %in% c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22')),]
load('EC_gene.Rdata')
ge1 = intersect(ge1,a$V5)
ge2 = intersect(ge2,a$V5)
# ge = c('SLC16A1','SLC3A2','CLDN5','SLC7A5','SLC2A1','ABCC2','ZIC2','LEF1','SLCO1C1')
#
di = c('AD','DM','EP','MS','PD','I60','I61','I62','I63','I64','I65','I67')
fin1 = as.data.frame(array(,dim = c(0,0)))
fin2 = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(di)){
	#
	temp = read.table(paste0('../tmp',i,'.tsv'))
	#
	for (j in 1:length(ge1)){
		#
		gen = ge1[j]
		aa = unique(a[which(a$V5 == gen),c(1,2)])
		pos = aa[1,2]
		chr = aa[1,1]
		chr = substr(chr,4,nchar(chr))
		pos1 = pos - 1000000
		pos2 = pos + 1000000
		#
		tem = na.omit(temp[which(temp$V1 == chr),])
		num1 = which(tem$V2 > pos1)
		num2 = which(tem$V2 < pos2)
		num3 = intersect(num1,num2)
		te = tem[num3,]
		tt = min(te$V15)
		ds = min(abs(te$V2[which(te$V15 == tt)] - pos))
		mat = as.data.frame(array(,dim = c(1,4)))
		mat[1,1] = di[i]
		mat[1,2] = gen
		mat[1,3] = tt
		mat[1,4] = ds
		fin1 = rbind(fin1,mat)
	}
	#
	for (j in 1:length(ge2)){
		#
		gen = ge2[j]
		aa = unique(a[which(a$V5 == gen),c(1,2)])
		pos = aa[1,2]
		chr = aa[1,1]
		chr = substr(chr,4,nchar(chr))
		pos1 = pos - 1000000
		pos2 = pos + 1000000
		#
		tem = na.omit(temp[which(temp$V1 == chr),])
		num1 = which(tem$V2 > pos1)
		num2 = which(tem$V2 < pos2)
		num3 = intersect(num1,num2)
		te = tem[num3,]
		tt = min(te$V15)
		ds = min(abs(te$V2[which(te$V15 == tt)] - pos))
		mat = as.data.frame(array(,dim = c(1,4)))
		mat[1,1] = di[i]
		mat[1,2] = gen
		mat[1,3] = tt
		mat[1,4] = ds
		fin2 = rbind(fin2,mat)
	}
	print(i)
	print(di[i])
}
save(fin1, fin2, file = 'fin_v2.Rdata')
#
load('fin_v2.Rdata')
di = c('AD','DM','EP','MS','PD','I60','I61','I62','I63','I64','I65','I67')
##
for (i in 1:length(di)){
	#
	temp = di[i]
	tem1 = fin1[which(fin1$V1 == temp & fin1$V3 <= 5e-8),]
	tem2 = fin2[which(fin2$V1 == temp & fin2$V3 <= 5e-8),]
	tem1$type = 'EC-associated'
	tem2$type = 'non-associated'
	tem = rbind(tem1,tem2)
	#
	tem$pv = -log(tem$V3,10)
	# tem$score[which(tem$score > 15)] = 15
	p1 = ggplot(tem,aes(x = type, y = V4)) + geom_boxplot()
	p2 = ggplot(tem,aes(x = type, y = pv)) + geom_boxplot() + ylim(0,20)
	assign(paste0('x',i),p1)
	assign(paste0('y',i),p2)
}
library(patchwork)
x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+plot_layout(nrow = 3) 
# y1+y2+y3+y4+y5+y6+y7+y8+y9+plot_layout(nrow = 3) 




















