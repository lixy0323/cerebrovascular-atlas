###################################################################################################
#3dpf
#细胞表达信息 细胞注释信息
g3 = read.csv('3-6-11dpf-final/vas_cell_matrix_annotation_F314.csv')
g3 = g3[,-c(40:41)]
#血管注释信息
v3 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F314.csv')
v3 = v3[,-1]
###
cg3 = g3
###
cg3$cell_type[which(cg3$cell_type == 'VEC_1')] = 'VEC'
cg3$cell_type[which(cg3$cell_type == 'VEC_2')] = 'VEC'
###
###################################################################################################
###################################################################################################
#6dpf
#细胞表达信息 细胞注释信息
g6 = read.csv('3-6-11dpf-final/vas_cell_matrix_annotation_F610.csv')
g6 = g6[,-c(40:41)]
#血管注释信息
v6 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F610.csv')
v6 = v6[,-1]
###
cg6 = g6
###
cg6$cell_type[which(cg6$cell_type == 'VEC_1')] = 'VEC'
cg6$cell_type[which(cg6$cell_type == 'VEC_2')] = 'VEC'
###
###################################################################################################
###################################################################################################
#11dpf
#细胞表达信息 细胞注释信息
g11 = read.csv('3-6-11dpf-final/vas_cell_matrix_annotation_F1109.csv')
g11 = g11[,-c(40:41)]
#血管注释信息
v11 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F1109.csv')
v11 = v11[,-1]
###
cg11 = g11
###
cg11$cell_type[which(cg11$cell_type == 'VEC_1')] = 'VEC'
cg11$cell_type[which(cg11$cell_type == 'VEC_2')] = 'VEC'
###
###################################################################################################
####
# cap_gene = c('cldn5a','slc16a1slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','plvapb','adgrl1a','slc2a1a','bcl6b')
# ge = read.csv('D:/工作/KEGG/kegg_zebrafish.csv')
#
# for (i in 1:length(cap_gene)){
	#
	# gene = cap_gene[i]
	# num = grep(gene,ge$geneID)
	# pw = ge$Description[num]
	# print(c(i,gene))
	# print(pw)
# }
###################################################################
###################################################################
###################################################################
colnames(cg3)[which(colnames(cg3) == 'zgc.158423')]='zgc:158423'
colnames(cg6)[which(colnames(cg6) == 'zgc158423')]='zgc:158423'
colnames(cg11)[which(colnames(cg11) == 'zgc.158423')]='zgc:158423'
v3$region_name[which(v3$region_name == 'willis')]='Willis'
v3$region_name[which(v3$region_name == 'dorsal')]='Dorsal'

vn1 = unique(v3$region_name)
vn2 = unique(v6$region_name)
vn3 = unique(v11$region_name)
vn = intersect(c(vn1,vn2),vn3)
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
vn = intersect(vname$vascular_name,vn)
#
# xge1 = c('slc16a1a','zgc:158423','slc2a1a','slc7a5')
# xge2 = c('slco1d1','abcb4','abcg2a','abcc2')
# xge = c('cldn5a','slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','bcl6b','plvapb','adgrl1a','slc2a1a')
# xge = c('cldn5a')
xge = colnames(cg3)[c(2:13,15:20,22:39)]
for (x in 1:length(xge)){
	mat1 = as.data.frame(array(,dim = c(length(vn),0)))
	mat2 = as.data.frame(array(,dim = c(length(vn),0)))
	mat3 = as.data.frame(array(,dim = c(length(vn),0)))
	mat1[,1] = vn
	mat2[,1] = vn
	mat3[,1] = vn
	mat1[,2] = '3dpf'
	mat2[,2] = '6dpf'
	mat3[,2] = '9dpf'
	###
	mk = xge[x]
	for (i in 1:length(vn)){
		###
		tmp = vn[i]
		###
		ri1 = v3$region_id[which(v3$region_name == tmp)]
		ri2 = v6$region_id[which(v6$region_name == tmp)]
		ri3 = v11$region_id[which(v11$region_name == tmp)]
		###
		tmp1 = cg3[which(cg3$blood %in% ri1 & cg3$cell_type == 'CapEC'),]
		tmp2 = cg6[which(cg6$blood %in% ri2 & cg6$cell_type == 'CapEC'),]
		tmp3 = cg11[which(cg11$blood %in% ri3 & cg11$cell_type == 'CapEC'),]
		###
		x1 = which(colnames(tmp1) == mk)
		tem1 = tmp1[,x1]
		x2 = which(colnames(tmp2) == mk)
		tem2 = tmp2[,x2]
		x3 = which(colnames(tmp3) == mk)
		tem3 = tmp3[,x3]
		##
		mat1[i,3] = mean(tem1)
		mat1[i,4] = length(tem1)
		mat2[i,3] = mean(tem2)
		mat2[i,4] = length(tem2)
		mat3[i,3] = mean(tem3)
		mat3[i,4] = length(tem3)
	}
	###
	vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
	# library(forcats)
	library(RColorBrewer)
	cl = brewer.pal(11,'Spectral')[6:1]
	mat = rbind(mat1,mat2,mat3)
	mat$V1 = factor(mat$V1, levels = rev(vname$vascular_name))
	mat[which(is.na(mat) == T,arr.ind = T)] = 0
	mat$region = rep(vname$location,3)
	mat$V3[which(mat$V3 > 0.4)] = 0.4
	p = ggplot(mat,aes(x = V2,y = V1)) + geom_point(aes(size = V4,color = V3)) + 
	scale_color_gradientn(limits = c(0,0.4), colours = cl) + #limits = c(0,0.5))
	coord_fixed() +
	theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())
	assign(paste0('p',x),p)
}
library(patchwork)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20+p21+p22+p23+p24+p25+p26+p27+p28+p29+p30+p31+p32+p33+p34+p35+p36+plot_layout(nrow = 2, guides = "collect")
#
###############
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
x=27
mat1 = as.data.frame(array(,dim = c(length(vn),0)))
mat2 = as.data.frame(array(,dim = c(length(vn),0)))
mat3 = as.data.frame(array(,dim = c(length(vn),0)))
mat1[,1] = vn
mat2[,1] = vn
mat3[,1] = vn
mat1[,2] = '3dpf'
mat2[,2] = '6dpf'
mat3[,2] = '9dpf'
###
mk = xge[x]
for (i in 1:length(vn)){
	###
	tmp = vn[i]
	###
	ri1 = v3$region_id[which(v3$region_name == tmp)]
	ri2 = v6$region_id[which(v6$region_name == tmp)]
	ri3 = v11$region_id[which(v11$region_name == tmp)]
	###
	tmp1 = cg3[which(cg3$blood %in% ri1),]
	tmp2 = cg6[which(cg6$blood %in% ri2),]
	tmp3 = cg11[which(cg11$blood %in% ri3),]
	###
	x1 = which(colnames(tmp1) == mk)
	tem1 = tmp1[,x1]
	x2 = which(colnames(tmp2) == mk)
	tem2 = tmp2[,x2]
	x3 = which(colnames(tmp3) == mk)
	tem3 = tmp3[,x3]
	##
	mat1[i,3] = mean(tem1)
	mat1[i,4] = length(tem1)
	mat2[i,3] = mean(tem2)
	mat2[i,4] = length(tem2)
	mat3[i,3] = mean(tem3)
	mat3[i,4] = length(tem3)
}
###
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')




