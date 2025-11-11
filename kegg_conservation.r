a = read.csv('capEC_kegg_1.csv')
a = a[which(a$pvalue < 0.05),]
b = read.csv('capEC_kegg_2.csv')
b = b[which(b$pvalue < 0.05),]
d = read.csv('capEC_kegg_3.csv')
d = d[which(d$pvalue < 0.05),]
#
#
pw = intersect(intersect(a$Description,b$Description),d$Description)
#
mat = as.data.frame(array(,dim = c(0,0)))
for (i in 1:length(pw)){
	#
	temp = pw[i]
	temp1 = a[which(a$Description == temp),c(2,5,8,12)]
	temp2 = b[which(b$Description == temp),c(2,5,8,12)]
	temp3 = d[which(d$Description == temp),c(2,5,8,12)]
	matt = as.data.frame(array(,dim = c(3,4)))
	matt[1,] = temp1
	matt[2,] = temp2 
	matt[3,] = temp3
	matt[,5] = c('H','M','Z')
	mat = rbind(mat,matt)
}











