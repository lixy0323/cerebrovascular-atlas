###################################################################################################
#3dpf
#细胞表达信息 细胞注释信息
a3 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f314-1230.csv')
g3 = read.csv('3-6-11dpf-final/vas_cell_matrix_annotation_F314.csv')
g3 = g3[,-c(40:41)] #wrong position
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
a6 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f610-1119.csv')
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
a11 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f1109-0312.csv')
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
#cg3
gg3 = cg3[,2:39]
fin3 = apply(gg3,2,sum)
fin3 = as.data.frame(fin3)
fin3$gene = rownames(fin3)
fin1 = fin3[order(fin3[,1],decreasing = F),]
xx = fin1[1:37,1]
names(xx) = fin1$gene[1:37]
barplot(xx,las=2,horiz = T,xlim = c(0,5000))
#
###################################################################################################
#
#a3
# a = a3
# a = a[-which(a$cell_type == 'UnEC'),]
# a$cell_type[which(a$cell_type == 'VEC_1')] = 'VEC'
# a$cell_type[which(a$cell_type == 'VEC_2')] = 'VEC'
# num = unique(a$zstack)
# for (i in 1:num){
	
	# temp = a[which(a$zstack == i),]
# }
###################################################################################################
###################################cg3 cg6 cg11
ct = c('AEC','VEC','AngEC','CapEC')
#
cg = cg3
for (i in 1:length(ct)){
	#
	temp = ct[i]
	temp1 = cg[which(cg$cell_type == temp),c(2:13,15:20,22:39)]
	temp2 = apply(temp1,2,mean)
	
	
	
	
	
	



}































