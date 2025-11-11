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
xge = c('cldn5a','slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','bcl6b','plvapb','adgrl1a','slc2a1a')
# xge = c('cldn5a')
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
		mat3[i,3] = mean(tem2)
		mat3[i,4] = length(tem2)
	}
	###
	vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
	# library(forcats)
	library(RColorBrewer)
	cl = brewer.pal(11,'Spectral')[6:1]
	mat = rbind(mat1,mat2,mat3)
	mat$V1 = factor(mat$V1, levels = rev(vname$vascular_name))
	mat = na.omit(mat)
	mat$region = rep(vname$location,3)
	p = ggplot(mat,aes(x = V2,y = V1)) + geom_point(aes(size = V4,color = V3)) + 
	scale_color_gradientn(colours =cl) + #limits = c(0,0.5))
	coord_fixed() 
	assign(paste0('p',x),p)
}
####################################################################################
####################################################################################
####################################################################################
####################################################################################
############################draw with diff brain region#############################
colnames(cg3)[which(colnames(cg3) == 'zgc.158423')]='zgc:158423'
colnames(cg6)[which(colnames(cg6) == 'zgc158423')]='zgc:158423'
colnames(cg11)[which(colnames(cg11) == 'zgc.158423')]='zgc:158423'
v3$region_name[which(v3$region_name == 'willis')]='Willis'
v3$region_name[which(v3$region_name == 'dorsal')]='Dorsal'
v11$region_name[which(v11$region_name == 'Wiliis')]='Willis'
vn1 = unique(v3$region_name)
vn2 = unique(v6$region_name)
vn3 = unique(v11$region_name)
vn = intersect(c(vn1,vn2),vn3)
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
vn = intersect(vname$vascular_name,vn)
vm = vname$location
#
# xge = c('slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','slc2a1a','abcb4','abcg2a','abcc2')
xge = c('cldn5a')
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
	vname = read.csv('vascular_name.csv')
	#
	mat1$location = vm
	mat2$location = vm
	mat3$location = vm
	uvm = unique(vm)
	matt1 = as.data.frame(array(,dim = c(length(uvm),4)))
	matt2 = as.data.frame(array(,dim = c(length(uvm),4)))
	matt3 = as.data.frame(array(,dim = c(length(uvm),4)))
	for (i in 1:length(uvm)){
		temp = uvm[i]
		tem = mat1[which(mat1$location == temp),]
		matt1[i,1] = temp
		matt1[i,2] = '3dpf'
		matt1[i,3] = mean(tem[,3])
		matt1[i,4] = sum(tem[,4])
		tem = mat2[which(mat2$location == temp),]
		matt2[i,1] = temp
		matt2[i,2] = '6dpf'
		matt2[i,3] = mean(tem[,3])
		matt2[i,4] = sum(tem[,4])
		tem = mat3[which(mat3$location == temp),]
		matt3[i,1] = temp
		matt3[i,2] = '9dpf'
		matt3[i,3] = mean(tem[,3])
		matt3[i,4] = sum(tem[,4])
	}
	# library(forcats)
	library(RColorBrewer)
	cl = brewer.pal(11,'Spectral')[6:1]
	mat = rbind(matt1,matt2,matt3)
	mat$V1 = factor(mat$V1, levels = rev(uvm))
	mat = na.omit(mat)
	mat$V3 = scale(mat$V3)
	p = ggplot(mat,aes(x = V2,y = V1)) + geom_point(aes(size = V4,color = V3)) + 
	#scale_color_gradientn(colours =cl,limits = c(0,xx)) +
	scale_color_gradientn(colours =cl) +	
	coord_fixed() 
	assign(paste0('p',x),p)
}
p1+p2+p3+p4+p5+p6+p7+p8+p9+plot_layout(ncol = 5)
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#slc2a1a
colnames(cg3)[which(colnames(cg3) == 'zgc.158423')]='zgc:158423'
colnames(cg6)[which(colnames(cg6) == 'zgc158423')]='zgc:158423'
v3$region_name_pri[which(v3$region_name_pri == 'willis')]='Willis'
v3$region_name_pri[which(v3$region_name_pri == 'dorsal')]='Dorsal'
vn1 = unique(v3$region_name_pri)
vn2 = unique(v6$region_name_pri)
# vn = intersect(vn1,vn2)
vname = read.csv('vascular_name.csv')
vn = vname$vascular_name
vm = vname$location
#
xge = c('slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','slc2a1a')
x = 10
###
d3 = as.data.frame(array(,dim = c(0,0)))
d6= as.data.frame(array(,dim = c(0,0)))
###
mk = xge[x]
te = c('AMsCtA','MMsCtA','PMsCtA')
for (i in 1:length(te)){
	###
	tmp = te[i]
	###
	ri1 = v3$region_id[which(v3$region_name_pri == tmp)]
	ri2 = v6$region_id[which(v6$region_name_pri == tmp)]
	###
	# tmp1 = cg3[which(cg3$region_id %in% ri1),]
	# tmp2 = cg6[which(cg6$region_id %in% ri2),]
	tmp1 = cg3[which(cg3$region_id %in% ri1 & cg3$cell_type == 'CapEC'),]
	tmp2 = cg6[which(cg6$region_id %in% ri2 & cg6$cell_type == 'CapEC'),]
	###
	d3 = rbind(d3,tmp1)
	d6 = rbind(d6,tmp2)
}
###
ds3 = d3[which(d3$slc2a1a > 0),]
ds3[(nrow(ds3) + 1),] = ds3[nrow(ds3),]
ds3[nrow(ds3),c(4,5,6)] = c(0,0,0)
ds3[(nrow(ds3) + 1),] = ds3[nrow(ds3),]
ds3[nrow(ds3),c(4,5,6)] = c(2200,3500,22)
ds3$zstack = ds3$zstack - 2
ds6 = d6[which(d6$slc2a1a > 0),]
ds6[(nrow(ds6) + 1),] = ds6[nrow(ds6),]
ds6[nrow(ds6),c(4,5,6)] = c(0,0,0)
ds6[(nrow(ds6) + 1),] = ds6[nrow(ds6),]
ds6[nrow(ds6),c(4,5,6)] = c(2200,3500,18)
plot_ly(data = ds3,x=~x,y=~y,z=~zstack,type = 'scatter3d',mode = 'markers',marker = list(size = 1),color = ds3$slc2a1a)
plot_ly(data = ds6,x=~x,y=~y,z=~zstack,type = 'scatter3d',mode = 'markers',marker = list(size = 1),color = ds6$slc2a1a)
#不可信 未对准
































###################################################################
###################################################################
###################################################################
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
colnames(cg3)[which(colnames(cg3) == 'zgc.158423')]='zgc:158423'
colnames(cg6)[which(colnames(cg6) == 'zgc158423')]='zgc:158423'
colnames(cg11)[which(colnames(cg11) == 'zgc.158423')]='zgc:158423'
v3$region_name[which(v3$region_name == 'willis')]='Willis'
v3$region_name[which(v3$region_name == 'dorsal')]='Dorsal'
#cg
xge = c('cldn5a','slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','bcl6b','plvapb','adgrl1a','slc2a1a')
#
mat = as.data.frame(array(,dim = c(45,4)))
mat[,1] = rep(xge,3)
mat[,2] = c(rep('d3',15),rep('d6',15),rep('d9',15))
t4 = c()
t5 = c()
for (z in c(3,6,11)){
	temp = get(paste0('cg',z))
	temp2=temp[which(temp$cell_type == 'CapEC'),]
	temp3 = temp2[,xge]
	temp4 = apply(temp3,2,sum)/nrow(temp3)
	temp5 = apply(temp3,2,function(x) length(which(x > 0))/length(x))
	#
	t4 = c(t4,temp4)
	t5 = c(t5,temp5)
}
mat[,3] = t4
mat[,4] = t5
ggplot(mat,aes(x = V2, y = V1)) + geom_point(aes(color = V3,size = V4))
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
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
# cap_gene = c('cldn5a','slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','plvapb','adgrl1a','slc2a1a','bcl6b')
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
xge1 = c('slc16a1a','zgc:158423','slc2a1a','slc7a5')
xge2 = c('slco1d1','abcb4','abcg2a','abcc2')
#xge = c('cldn5a','slc16a1a','zgc:158423','slco1c1','slc7a5','slco1d1','abcb4','abcg2a','abcc2','abcb4','zic2b','bcl6b','plvapb','adgrl1a','slc2a1a')
# xge = c('cldn5a')
# for (x in 1:length(xge)){
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
mk = xge2
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
	x1 = which(colnames(tmp1) %in% mk)
	tem1 = tmp1[,x1]
	x2 = which(colnames(tmp2) %in% mk)
	tem2 = tmp2[,x2]
	x3 = which(colnames(tmp3) %in% mk)
	tem3 = tmp3[,x3]
	##
	mat1[i,3] = sum(tem1)
	mat1[i,4] = nrow(tem1)
	# mat1[i,5] = mat1[i,3]/mat1[i,4]
	mat2[i,3] = sum(tem2)
	mat2[i,4] = nrow(tem2)
	# mat2[i,5] = mat2[i,3]/mat2[i,4]
	mat3[i,3] = sum(tem3)
	mat3[i,4] = nrow(tem3)
	# mat3[i,5] = mat3[i,3]/mat3[i,4]
}
mat1[,5] = mat1[,3]/sum(mat1[,3])
mat2[,5] = mat2[,3]/sum(mat2[,3])
mat3[,5] = mat3[,3]/sum(mat3[,3])
###
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
# library(forcats)
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')[6:1]
mat = rbind(mat1,mat2,mat3)
mat$V1 = factor(mat$V1, levels = rev(vname$vascular_name))
mat = na.omit(mat)
mat$region = rep(vname$location,3)
p = ggplot(mat,aes(x = V2,y = V1)) + geom_point(size = 5,aes(color = V5)) + 
scale_color_gradientn(colours =cl) + #limits = c(0,0.5))
coord_fixed() 
assign(paste0('p',2),p)
# }

