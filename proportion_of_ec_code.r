rm(list = ls())

vas_local = read.csv("./localjupyter/vascular_name_and_location.csv")
colnames(vas_local) = c("location", "region_name", "index")



vas_cell_3dpf = read.csv("./localjupyter/vas_cell_matrix_annotation_F314.csv", row.names = 1)
dim(vas_cell_3dpf)
vas_cell_3dpf = vas_cell_3dpf[!vas_cell_3dpf$cell_type %in% c("other", "UnEC"),]
dim(vas_cell_3dpf)
colnames(vas_cell_3dpf) = gsub("blood", "region_id", colnames(vas_cell_3dpf))
colnames(vas_cell_3dpf) = gsub("zgc.158423", "zgc:158423", colnames(vas_cell_3dpf))

vas_anno_3dpf = read.csv("./localjupyter/region_annotation_with_feature_F314.csv", row.names = 1)
vas_anno_3dpf = na.omit(vas_anno_3dpf)
vas_anno_3dpf[vas_anno_3dpf$region_name == "willis",]$region_name = "Willis"
vas_anno_3dpf[vas_anno_3dpf$region_name == "dorsal",]$region_name = "Dorsal"
vas_anno_3dpf = vas_anno_3dpf[!duplicated(vas_anno_3dpf$region_id), c("region_id", "region_name")]
vas_anno_3dpf = vas_anno_3dpf[vas_anno_3dpf$region_name %in% vas_local$region_name,]
vas_cell_3dpf = vas_cell_3dpf[vas_cell_3dpf$region_id %in% vas_anno_3dpf$region_id,]
dim(vas_cell_3dpf)


vas_cell_3dpf$region_name = vas_anno_3dpf[match(vas_cell_3dpf$region_id, vas_anno_3dpf$region_id),]$region_name
vas_cell_3dpf$location = vas_local[match(vas_cell_3dpf$region_name, vas_local$region_name),]$location
vas_cell_3dpf[vas_cell_3dpf$cell_type %in% c("VEC_1", "VEC_2"),]$cell_type = "VEC"
dim(vas_cell_3dpf)

vas_cell_6dpf = read.csv("./localjupyter/vas_cell_matrix_annotation_F610.csv", row.names = 1)
dim(vas_cell_6dpf)
vas_cell_6dpf = vas_cell_6dpf[!vas_cell_6dpf$cell_type %in% c("other", "UnEC"),]
dim(vas_cell_6dpf)
colnames(vas_cell_6dpf) = gsub("blood", "region_id", colnames(vas_cell_6dpf))
colnames(vas_cell_6dpf) = gsub("zgc158423", "zgc:158423", colnames(vas_cell_6dpf))

vas_anno_6dpf = read.csv("./localjupyter/region_annotation_with_feature_F610.csv", row.names = 1)
vas_anno_6dpf = na.omit(vas_anno_6dpf)
vas_anno_6dpf[vas_anno_6dpf$region_name == "DorsalS",]$region_name = "Dorsal"
vas_anno_6dpf[vas_anno_6dpf$region_name == "PrCmS",]$region_name = "PrCm"
unique(vas_anno_6dpf$region_name)
vas_anno_6dpf = vas_anno_6dpf[!duplicated(vas_anno_6dpf$region_id), c("region_id", "region_name")]
vas_anno_6dpf = vas_anno_6dpf[vas_anno_6dpf$region_name %in% vas_local$region_name,]
vas_cell_6dpf = vas_cell_6dpf[vas_cell_6dpf$region_id %in% vas_anno_6dpf$region_id,]
dim(vas_cell_6dpf)
vas_cell_6dpf$region_name = vas_anno_6dpf[match(vas_cell_6dpf$region_id, vas_anno_6dpf$region_id),]$region_name
vas_cell_6dpf$location = vas_local[match(vas_cell_6dpf$region_name, vas_local$region_name),]$location
vas_cell_6dpf[vas_cell_6dpf$cell_type %in% c("VEC_1", "VEC_2"),]$cell_type = "VEC"



vas_cell_11dpf = read.csv("./localjupyter/vas_cell_matrix_annotation_F1109.csv", row.names = 1)
dim(vas_cell_11dpf)
vas_cell_11dpf = vas_cell_11dpf[!vas_cell_11dpf$cell_type %in% c("other", "UnEC"),]
dim(vas_cell_11dpf)
colnames(vas_cell_11dpf) = gsub("blood", "region_id", colnames(vas_cell_11dpf))
colnames(vas_cell_11dpf) = gsub("zgc.158423", "zgc:158423", colnames(vas_cell_11dpf))

vas_anno_11dpf = read.csv("./localjupyter/region_annotation_with_feature_F1109.csv", row.names = 1)
vas_anno_11dpf = na.omit(vas_anno_11dpf)
vas_anno_11dpf[vas_anno_11dpf$region_name == "Wiliis",] = "Willis"
unique(vas_anno_11dpf$region_name)

vas_anno_11dpf = vas_anno_11dpf[!duplicated(vas_anno_11dpf$region_id), c("region_id", "region_name")]
vas_anno_11dpf = vas_anno_11dpf[vas_anno_11dpf$region_name %in% vas_local$region_name,]
vas_cell_11dpf = vas_cell_11dpf[vas_cell_11dpf$region_id %in% vas_anno_11dpf$region_id,]

vas_cell_11dpf$region_name = vas_anno_11dpf[match(vas_cell_11dpf$region_id, vas_anno_11dpf$region_id),]$region_name
vas_cell_11dpf$location = vas_local[match(vas_cell_11dpf$region_name, vas_local$region_name),]$location
vas_cell_11dpf[vas_cell_11dpf$cell_type %in% c("VEC_1", "VEC_2"),]$cell_type = "VEC"
dim(vas_cell_11dpf)





vas_cell_3dpf$group = "3dpf"
vas_cell_6dpf$group = "6dpf"
vas_cell_11dpf$group = "11dpf"

reg_3 = read.csv("./localjupyter/region_area_F314.csv", row.names = 1)
reg_6 = read.csv("./localjupyter/region_area_F610.csv", row.names = 1)
reg_11 = read.csv("./localjupyter/region_area_F1109.csv", row.names = 1)

vas_cell_3dpf$region_area = reg_3[match(vas_cell_3dpf$region_id, reg_3$region_id),]$region_area
vas_cell_6dpf$region_area = reg_6[match(vas_cell_6dpf$region_id, reg_6$region_id),]$region_area
vas_cell_11dpf$region_area = reg_11[match(vas_cell_11dpf$region_id, reg_11$region_id),]$region_area

vas_cell = rbind(vas_cell_3dpf, vas_cell_6dpf, vas_cell_11dpf)
vas_cell = vas_cell[vas_cell$location %in% vas_local$location,]
cell_in_vas = vas_cell %>% group_by(location,cell_type, group) %>% summarise(count = n(), 
                                                                             count_ =  length(cell_type), 
                                                                             .groups = "keep")

cell_in_vas = cell_in_vas %>% group_by(group, location) %>% reframe(
  cell_type = cell_type,
  count = count, 
  count_ = count_,
  percent = count_ / sum(count_)
)

cell_in_vas$group = factor(cell_in_vas$group, levels = c("3dpf", "6dpf", "11dpf"))
unique(cell_in_vas$location)
cell_in_vas$location = factor(cell_in_vas$location, 
                              levels = c("Prosencephalon","Mesencephalon","Metenbcephalon","Lateral","Ventral"))
write.csv(cell_in_vas, "./折线图.csv")
ggplot(data = cell_in_vas)+labs(y = "Frequence", x = NULL, color = NULL)+
  geom_line(aes(x = group, y = percent, group = cell_type),width = 0.75,color = "black", linewidth = 0.25)+
  geom_point(aes(x = group, y = percent, color = cell_type), size=  2.5)+
  facet_grid(rows = vars(location))+
  scale_y_continuous(expand = c(0.05,0,0.05,0), limits = c(0,1))+
  scale_x_discrete(expand = c(0,0.1,0,0.1))+
  scale_color_manual(values = c( "#e7ba90", "#90e7ba","#ba90e7","#90bae7"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        panel.border = element_rect(fill = "transparent", color = "black"),
        panel.grid = element_blank(),
        legend.position = "top")
  
# vas_cells = vas_cell[,c("slc7a5", "slc2a1a", "stab1","stab2", 'region_name','cell_type_t', 'group')]
# 
# vas_cells = vas_cells[(vas_cells$slc7a5 > 0 | vas_cells$slc2a1a > 0) |  
#                         (vas_cells$stab2 > 0 | vas_cells$stab2 > 0),]
# vas_cells$status = ""
# vas_cells[(vas_cells$slc7a5 > 0 & vas_cells$slc2a1a > 0) & (vas_cells$stab2 > 0 | vas_cells$stab2 > 0), ]$status = "VEN-CapEC"
# vas_cells[(vas_cells$slc7a5 > 0 & vas_cells$slc2a1a > 0) & (vas_cells$stab2 == 0 | vas_cells$stab2 == 0), ]$status = "CapEC"
# vas_cells[(vas_cells$slc7a5 == 0 | vas_cells$slc2a1a == 0) & (vas_cells$stab2 > 0 & vas_cells$stab2 > 0), ]$status = "VEC"
cell_in_vas_3 = vas_cell %>% group_by(region_name,cell_type, group) %>% summarise(count = n(), 
                                                                               count_ = length(cell_type),
                                                                                .groups = "keep")
# 
# cell_in_vas_3 = cell_in_vas_3[cell_in_vas_3$status != "",]

cell_in_vas_3$group = factor(cell_in_vas_3$group, levels = c("3dpf", "6dpf", "11dpf"))
cell_in_vas_3 = cell_in_vas_3[cell_in_vas_3$region_name %in% vas_local$region_name,]
cell_in_vas_3$region_name = factor(cell_in_vas_3$region_name, levels = unique(vas_local$region_name))
# cell_in_vas_3 = cell_in_vas_3[cell_in_vas_3$cell_type %in% c("VEC-CapEC", "CapEC", "VEC", "VEC-AngEC"),]
cell_in_vas_3$cell_type = factor(cell_in_vas_3$cell_type, levels = c("AEC","AngEC","VEC","CapEC"))
cell_in_vas_3 = cell_in_vas_3[!cell_in_vas_3$region_name %in% c("CaDI", "PrA", "PrV",'Dorsal'),]
write.csv(cell_in_vas_3, "./source_data_for_Fig4c.csv")
ggplot(data = cell_in_vas_3)+labs(y = "Frequence", x = NULL,fill = "Cell Type")+
  geom_bar(aes(x = group, y = count, fill = cell_type),width = 0.75,
           stat = "identity", color = "black", linewidth = 0.25, position = "fill")+
  facet_grid(~region_name)+
  scale_y_continuous(expand = c(0,0,0,0))+
  scale_fill_manual(values = c( "#ecc7a5","#f2cbd4","#cbe7f2","#e1cff5"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 8),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))