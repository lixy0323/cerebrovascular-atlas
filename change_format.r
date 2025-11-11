######
vname = read.csv('vascular_name_and_location_v.2.csv')
###
a = read.table('pheatmap_data.txt',head = T)
#
va = unique(a$cell_type)
#
for (i in 1:length(va)){
	#
	temp = va[i]
	b = a[which(a$cell_type == temp),]
	#
	dp = c('3dpf','6dpf','11dpf')
	#
	mat = as.data.frame(array(,dim = c(12,3)))
	for (j in 1:3){
		#
		temp2 = dp[j]
		d = b[which(b$group == temp2),]
		mat[,j] = d$norm12
	}
	rownames(mat) = a$region_name[1:12]
	colnames(mat) = dp
	mat$vascular_name = rownames(mat)
	matb = merge(vname,mat,by = c('vascular_name'),sort = F)
	matd = matb[,c(3:5)]
	rownames(matd) = matb$vascular_name
	pheatmap(t(matd),cluster_row = F,cluster_col =F,color = colorRampPalette(colors = c("white","red"))(100))

}

#########################
vname = read.csv('vascular_name_and_location_v.2.csv')
a = read.csv('marker_gene_expression.csv')
#
gene = unique(a$variable)
#
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#
for (i in 1:length(gene)){
	#
	temp = gene[i]
	b = a[which(a$variable == temp),]
	#
	dp = c('3dpf','6dpf','11dpf')
	#
	mat = as.data.frame(array(,dim = c(12,3)))
	for (j in 1:3){
		#
		temp2 = dp[j]
		d = b[which(b$group == temp2),]
		mat[,j] = d$mean_expr
	}
	rownames(mat) = b$region_name[1:12]
	colnames(mat) = dp
	mat$vascular_name = rownames(mat)
	matb = merge(vname,mat,by = c('vascular_name'),sort = F)
	matd = matb[,c(3:5)]
	rownames(matd) = matb$vascular_name
	xx = pheatmap(t(matd),cluster_row = F,cluster_col =F,color = colorRampPalette(colors = c("white","red"))(100))
	save_pheatmap_pdf(xx,paste0('gene/',temp,'.pdf'),10,3)	
}
##################################################################################################################















