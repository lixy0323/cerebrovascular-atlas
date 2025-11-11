vname = read.csv('vascular_name.csv')
####
c3 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f314-1230.csv')
##
c3$cell_type[which(c3$cell_type == 'VEC_1')] = 'VEC'
c3$cell_type[which(c3$cell_type == 'VEC_2')] = 'VEC'
zs = unique(c3$zstack)
mat = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(zs)){
	tmp = c3[which(c3$zstack == zs[i]),]
	vs = unique(tmp$region_id)
	matt = as.data.frame(array(,dim = c(length(vs),7)))
	matt[,1] = zs[i]
	matt[,2] = vs
	###
	for (j in 1:length(vs)){
		#
		tmp2 = tmp[which(tmp$region_id == vs[j]),]
		matt[j,3] = length(which(tmp2$cell_type == 'AEC'))
		matt[j,4] = length(which(tmp2$cell_type == 'VEC'))
		matt[j,5] = length(which(tmp2$cell_type == 'AngEC'))
		matt[j,6] = length(which(tmp2$cell_type == 'CapEC'))
		matt[j,7] = length(which(tmp2$cell_type == 'UnEC'))
	}
	mat = rbind(mat,matt)
}
colnames(mat)[1:7] = c('zstack','region_id','AEC','VEC','AngEC','CapEC','UnEC')
mat$count = apply(mat[,3:6],1,function(x) length(which(x>0)))
#
# length(which(mat$AEC == 0 & mat$AngEC == 0 & mat$VEC == 0 & mat$CapEC == 0))
#
#vascular name  vname
v3 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F314.csv')
x3 = unique(v3[,c(2,3)])
matx = merge(mat,x3,by = c('region_id'),sort = F)
colnames(matx)[9] = 'vascular_name'
matx$vascular_name = toupper(matx$vascular_name)
vname$vascular_name = toupper(vname$vascular_name)
matv = merge(matx,vname,by = c('vascular_name'),sort = F)
##
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),5)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	mat1[i,1] = length(which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,2] = length(which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,3] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0))
	mat1[i,4] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0))
	mat1[i,5] = nrow(tmp) - mat1[i,1] - mat1[i,2] - mat1[i,3] - mat1[i,4]
	# print(c(num1,num2,num3,num4,num5,num6))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:28,1:28,1:28,1:28,1:28)
p1=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
############################################################################################
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),3)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	num1 = which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num2 = which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num3 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0)
	num4 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0)
	num5 = c(num1,num2,num3,num4)
	tmp2 = tmp[-num5,]
	mat1[i,1] = length(which(tmp2$VEC > tmp2$CapEC))
	mat1[i,2] = length(which(tmp2$VEC == tmp2$CapEC))
	mat1[i,3] = length(which(tmp2$VEC < tmp2$CapEC))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:28,1:28,1:28)
p2=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
vname = read.csv('vascular_name.csv')
####
c3 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f610-1119.csv')
##
c3$cell_type[which(c3$cell_type == 'VEC_1')] = 'VEC'
c3$cell_type[which(c3$cell_type == 'VEC_2')] = 'VEC'
zs = unique(c3$zstack)
mat = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(zs)){
	tmp = c3[which(c3$zstack == zs[i]),]
	vs = unique(tmp$region_id)
	matt = as.data.frame(array(,dim = c(length(vs),7)))
	matt[,1] = zs[i]
	matt[,2] = vs
	###
	for (j in 1:length(vs)){
		#
		tmp2 = tmp[which(tmp$region_id == vs[j]),]
		matt[j,3] = length(which(tmp2$cell_type == 'AEC'))
		matt[j,4] = length(which(tmp2$cell_type == 'VEC'))
		matt[j,5] = length(which(tmp2$cell_type == 'AngEC'))
		matt[j,6] = length(which(tmp2$cell_type == 'CapEC'))
		matt[j,7] = length(which(tmp2$cell_type == 'UnEC'))
	}
	mat = rbind(mat,matt)
}
colnames(mat)[1:7] = c('zstack','region_id','AEC','VEC','AngEC','CapEC','UnEC')
mat$count = apply(mat[,3:6],1,function(x) length(which(x>0)))
#
# length(which(mat$AEC == 0 & mat$AngEC == 0 & mat$VEC == 0 & mat$CapEC == 0))
#
#vascular name  vname
v3 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F610.csv')
x3 = unique(v3[,c(2,3)])
matx = merge(mat,x3,by = c('region_id'),sort = F)
colnames(matx)[9] = 'vascular_name'
matx$vascular_name = toupper(matx$vascular_name)
vname$vascular_name = toupper(vname$vascular_name)
matv = merge(matx,vname,by = c('vascular_name'),sort = F)
##
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),5)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	mat1[i,1] = length(which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,2] = length(which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,3] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0))
	mat1[i,4] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0))
	mat1[i,5] = nrow(tmp) - mat1[i,1] - mat1[i,2] - mat1[i,3] - mat1[i,4]
	# print(c(num1,num2,num3,num4,num5,num6))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:25,1:25,1:25,1:25,1:25)
p3=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
############################################################################################
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),3)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	num1 = which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num2 = which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num3 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0)
	num4 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0)
	num5 = c(num1,num2,num3,num4)
	tmp2 = tmp[-num5,]
	mat1[i,1] = length(which(tmp2$VEC > tmp2$CapEC))
	mat1[i,2] = length(which(tmp2$VEC == tmp2$CapEC))
	mat1[i,3] = length(which(tmp2$VEC < tmp2$CapEC))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:25,1:25,1:25)
p4=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
####
vname = read.csv('vascular_name.csv')
####
c3 = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f1109-0312.csv')
##
c3$cell_type[which(c3$cell_type == 'VEC_1')] = 'VEC'
c3$cell_type[which(c3$cell_type == 'VEC_2')] = 'VEC'
zs = unique(c3$zstack)
mat = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(zs)){
	tmp = c3[which(c3$zstack == zs[i]),]
	vs = unique(tmp$region_id)
	matt = as.data.frame(array(,dim = c(length(vs),7)))
	matt[,1] = zs[i]
	matt[,2] = vs
	###
	for (j in 1:length(vs)){
		#
		tmp2 = tmp[which(tmp$region_id == vs[j]),]
		matt[j,3] = length(which(tmp2$cell_type == 'AEC'))
		matt[j,4] = length(which(tmp2$cell_type == 'VEC'))
		matt[j,5] = length(which(tmp2$cell_type == 'AngEC'))
		matt[j,6] = length(which(tmp2$cell_type == 'CapEC'))
		matt[j,7] = length(which(tmp2$cell_type == 'UnEC'))
	}
	mat = rbind(mat,matt)
}
colnames(mat)[1:7] = c('zstack','region_id','AEC','VEC','AngEC','CapEC','UnEC')
mat$count = apply(mat[,3:6],1,function(x) length(which(x>0)))
#
# length(which(mat$AEC == 0 & mat$AngEC == 0 & mat$VEC == 0 & mat$CapEC == 0))
#
#vascular name  vname
v3 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F1109.csv')
x3 = unique(v3[,c(2,3)])
matx = merge(mat,x3,by = c('region_id'),sort = F)
colnames(matx)[9] = 'vascular_name'
matx$vascular_name = toupper(matx$vascular_name)
vname$vascular_name = toupper(vname$vascular_name)
matv = merge(matx,vname,by = c('vascular_name'),sort = F)
##
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),5)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	mat1[i,1] = length(which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,2] = length(which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0))
	mat1[i,3] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0))
	mat1[i,4] = length(which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0))
	mat1[i,5] = nrow(tmp) - mat1[i,1] - mat1[i,2] - mat1[i,3] - mat1[i,4]
	# print(c(num1,num2,num3,num4,num5,num6))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:25,1:25,1:25,1:25,1:25)
p5=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
############################################################################################
zs = unique(c3$zstack)
mat1 = as.data.frame(array(,dim = c(length(zs),3)))
for (i in 1:length(zs)){
	###
	tmp = matv[which(matv$zstack == zs[i]),]
	num1 = which(tmp$AEC != 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num2 = which(tmp$AEC == 0 & tmp$VEC != 0 & tmp$CapEC == 0 & tmp$UnEC == 0)
	num3 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC != 0 & tmp$UnEC == 0)
	num4 = which(tmp$AEC == 0 & tmp$VEC == 0 & tmp$CapEC == 0 & tmp$UnEC != 0)
	num5 = c(num1,num2,num3,num4)
	tmp2 = tmp[-num5,]
	mat1[i,1] = length(which(tmp2$VEC > tmp2$CapEC))
	mat1[i,2] = length(which(tmp2$VEC == tmp2$CapEC))
	mat1[i,3] = length(which(tmp2$VEC < tmp2$CapEC))
}
#
library(reshape2)
mat1 = mat1/apply(mat1,1,sum)
mat2 = melt(mat1)
mat2$zstack = c(1:25,1:25,1:25)
p6=ggplot(mat2,aes(x = zstack,y = value)) + geom_bar(stat = 'identity',aes(fill = variable))
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################
########################################################################################################################################################################################







































#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)

num = unique(c3$zstack)
	i=1
	ii = num[i]
	b = c3[which(c3$zstack == ii),]
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),4] = 3500
	b[nrow(b),5] = 2200
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),4] = 0
	b[nrow(b),5] = 0
	ggplot(b,aes(x = (2200 - y),y = (3500 - x))) + geom_point(aes(color = cell_type),size = 1) +# scale_color_manual(values = cc) + 
		coord_fixed() + 
		#theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		#theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  #theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'black',linewidth = 1),
			  #panel.background = element_rect(fill = "transparent",colour = NA),
			  #plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + annotate("text", x = 200, y = 3300, size = 15, label = i) 























