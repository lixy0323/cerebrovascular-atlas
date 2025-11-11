#
load('gene_list.rda')
ge = gene_list
#
mat = read.table('bbb_expression_gene.txt',head = T)
#
mat2 = mat[,2:4]
mat2[which(mat2 > 10, arr.ind = T)] = 10
rownames(mat2) = mat[,1]
#
mat3 = mat2[which(rownames(mat2) %in% ge$ns2ns),]


