###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#3d
a = read.csv('gene_coord_with_zstack_f314-1230.csv')
#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)

num = unique(a$zstack)

for (i in 1:length(num)){
	ii = num[i]
	b = a[which(a$zstack == ii),]
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 2200
	b[nrow(b),4] = 3500
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 0
	b[nrow(b),4] = 0
	ggplot(b,aes(x = y_location,y = (3500-x_location))) + geom_point(aes(color = feature_name),size = 0.2)+ scale_color_manual(values = cc) + 
		coord_fixed() + 
		#theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		#theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'white',linewidth = 1),
			  #panel.background = element_rect(fill = "transparent",colour = NA),
			  #plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + annotate("text", x = 200, y = 3300, size = 15, label = i) 

	ggsave(paste0(i,'.png'), width = 4.4, height = 7, units = "in")
}
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#6d
a = read.csv('gene_coord_with_zstack_f610-1119.csv')
#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)

num = unique(a$zstack)

for (i in 1:length(num)){
	ii = num[i]
	b = a[which(a$zstack == ii),]
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 3500
	b[nrow(b),4] = 2200
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 0
	b[nrow(b),4] = 0
	ggplot(b,aes(x = (2200 - x_location),y = y_location)) + geom_point(aes(color = feature_name),size = 0.2)+ scale_color_manual(values = cc) + 
		coord_fixed() + 
		#theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		#theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'white',size = 1),
			  #panel.background = element_rect(fill = "transparent",colour = NA),
			  #plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + annotate("text", x = 200, y = 3300, size = 15, label = i) 

	ggsave(paste0(i,'.png'), width = 4.4, height = 7, units = "in")
}
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
#11d
a = read.csv('gene_coord_with_zstack_F1109.csv')
#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)

num = unique(a$zstack)

for (i in 1:length(num)){
	ii = num[i]
	b = a[which(a$zstack == ii),]
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),2] = 2200
	b[nrow(b),3] = 3500
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),2] = 0
	b[nrow(b),3] = 0
	ggplot(b,aes(x = y_location,y = x_location)) + geom_point(aes(color = feature_name),size = 0.1)+ scale_color_manual(values = cc) + 
		coord_fixed() + 
		theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'white',size = 1),
			  panel.background = element_rect(fill = "transparent",colour = NA),
			  plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) #+ annotate("text", x = 200, y = 3300, size = 15, label = i) 

	ggsave(paste0(i,'.png'), width = 4.4, height = 7, units = "in")
}
##############################################################################################################################################################################
###cell position
# a = read.csv('3-6-11dpf-final/all_cell_coord_with_zstack_f1109-0312.csv')
a = read.csv('3-6-11dpf-final/vascular_cell_coord_with_zstack_f1109-0312.csv')
a = a[,-1]
#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)
#
num = unique(a$zstack)
#血管注释
v11 = read.csv('3-6-11dpf-final/region_annotation_with_feature_F1109.csv')
v11 = v11[,c(-1,-4)]
v11 = unique(v11)
vname = read.csv('3-6-11dpf-final/vascular_name_and_location.csv')
colnames(vname)[2] = 'region_name'
###
for (i in 1:length(num)){
	ii = num[i]
	b = a[which(a$zstack == ii),]
	#
	b = merge(b,v11,by = 'region_id',all.x = T)
	b = merge(b,vname,by = 'region_name',all.x = T)
	#
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),5] = 3500
	b[nrow(b),6] = 2200
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),5] = 0
	b[nrow(b),6] = 0
	#
	
	#
	ggplot(b,aes(x = y,y = x)) + geom_point(aes(color = location),size = 1)+ #scale_color_manual(values = cc) + 
		coord_fixed() + 
		theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  #theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'white',size = 1),
			  panel.background = element_rect(fill = "transparent",colour = NA),
			  plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) #+ annotate("text", x = 200, y = 3300, size = 15, label = i) 

	ggsave(paste0(i,'.png'), width = 4.4, height = 7, units = "in")
}
###




















###################################################################
a = read.csv('gene_coord_with_zstack_f314-1230.csv')
#
library(RColorBrewer)
cl = brewer.pal(11,'Spectral')
cc = colorRampPalette(cl)(38)
#
num = unique(a$zstack)
#
for (i in 1:length(num)){
	ii = num[i]
	b = a[which(a$zstack == ii),]
	#
	xm = mean(b$y_location)
	ym = mean(b$x_location)
	#
	
	
	
	
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 2200
	b[nrow(b),4] = 3500
	b[(nrow(b) + 1),] = b[nrow(b),]
	b[nrow(b),3] = 0
	b[nrow(b),4] = 0
	
	
	
	
	
	ggplot(b,aes(x = y_location,y = (3500-x_location))) + geom_point(aes(color = feature_name),size = 0.2)+ scale_color_manual(values = cc) + 
		coord_fixed() + 
		#theme_classic() + 
		theme(panel.background = element_rect(fill='white'), strip.background = element_rect(colour=NA, fill=NA), panel.border = element_rect(fill = NA, color = NA)) + 
		#theme_bw() +
	  theme(axis.ticks =element_blank()) +
	  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
	  theme(axis.title.x =element_blank()) +
	  theme(axis.title.y =element_blank()) + 
	  theme(legend.position = 'none') + 
	  theme(panel.grid.major=element_line(colour=NA),
	          panel.border = element_rect(color = 'white',linewidth = 1),
			  #panel.background = element_rect(fill = "transparent",colour = NA),
			  #plot.background = element_rect(fill = "transparent",colour = NA),
			  panel.grid.minor = element_blank()) + 
		scale_y_continuous(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + annotate("text", x = 200, y = 3300, size = 15, label = i) 

	ggsave(paste0(i,'.png'), width = 4.4, height = 7, units = "in")
}
######














