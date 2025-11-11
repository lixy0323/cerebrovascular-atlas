rm(list = ls())

options(stringAsFactors = F)
setwd("~/rproj/st")

library(openxlsx)




#############################################################################################
#############################################################################################
#############################################################################################
datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 1)
data = datas[datas$group == "length",]

colnames(data)
data$p_slc = ""
for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,8:13])
  # x = t.test(sample(x = data[i,2:6], size = 3, replace = F),sample(x = data[i,7:11], size = 3, replace = F))
  data$p_slc[i] = as.numeric(x$p.value)
}

data$p_zgc = ""

for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,14:20])
  data$p_zgc[i] = as.numeric(x$p.value)
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,21:26])
  data$p_cldc[i] = as.numeric(x$p.value)
}

data_length_3dpf = data
dat = data_length_3dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Length of vessels (μm)")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = .2)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), 
                width = .2, position=position_dodge(0.75), size = 0.5)+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,4000))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )
  
#############################################################################################

datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 1)
data = datas[datas$group == "number",]

data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}
data_number_3dpf = data

dat = data_number_3dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Number of vascular segments")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = .2)+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.5)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0, 120))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )


#############################################################################################


datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 1)
data = datas[datas$group == "straightness",]

data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}
data_straightness_3dpf = data

dat = data_straightness_3dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Vascular straightness")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.5)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,1.25))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )


#############################################################################################


datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 2)
data = datas
data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}

data_3dpf = data

dat = data_3dpf[,1:26]
colnames(dat)[1] = "para"
colnames(dat)
dat = reshape2::melt(dat)
head(dat)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)

dats = dat[dat$para == unique(dat$para)[1],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Total length of vessels (μm)")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,20000))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_blank()
    # axis.text.y = element_text(color = "black")
  )


dats = dat[dat$para == unique(dat$para)[2],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Number of vascular segments")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,500))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    # axis.text.y = element_text(color = "black")
    axis.text.y = element_blank()
  )


dats = dat[dat$para == unique(dat$para)[3],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Vascular straightness")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,1))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    # axis.text.y = element_text(color = "black")
    axis.text.y = element_blank()
  )





#############################################################################################
#############################################################################################
#############################################################################################
datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 3)
data = datas[datas$group == "length",]

colnames(data)
data$p_slc = ""
for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,8:13])
  # x = t.test(sample(x = data[i,2:6], size = 3, replace = F),sample(x = data[i,7:11], size = 3, replace = F))
  data$p_slc[i] = as.numeric(x$p.value)
}

data$p_zgc = ""

for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,14:20])
  data$p_zgc[i] = as.numeric(x$p.value)
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  x = t.test(data[i,2:7], data[i,21:26])
  data$p_cldc[i] = as.numeric(x$p.value)
}

data_length_6dpf = data
dat = data_length_6dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Length of vessels (μm)")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = .2)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), 
                width = .2, position=position_dodge(0.75), size = 0.5)+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,8000))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )

#############################################################################################

datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 3)
data = datas[datas$group == "number",]

data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}
data_number_6dpf = data

dat = data_number_6dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Number of vascular segements")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = .2)+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.5)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0, 200))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )


#############################################################################################


datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 3)
data = datas[datas$group == "straightness",]

data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}
data_straightness_6dpf = data

dat = data_straightness_6dpf[,1:26]
dat = reshape2::melt(dat)
head(dat)
dat$vas_name = factor(dat$vas_name, levels = unique(dat$vas_name))
table(dat$vas_name)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)
errbar = dat %>% group_by(vas_name, group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dat)+labs(x = NULL, y = "Vascular straightness")+
  geom_bar(aes(x = vas_name, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.65, color = "black", linewidth = .2)+
  geom_errorbar(data = errbar, aes(x = vas_name,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.5)+
  geom_point(aes(x = vas_name, y = value, fill = group),size = .75,
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,1.25))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )


#############################################################################################


datas = read.xlsx("./crispant-reconstruction.xlsx", sheet = 4)
data = datas
data$p_slc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,8:13]))
  if (any(sd == 0)){
    data$p_slc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,8:13])
    data$p_slc[i] = as.numeric(x$p.value)
  }
}


data$p_zgc = ""
for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,14:20]))
  if (any(sd == 0)){
    data$p_zgc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,14:20])
    data$p_zgc[i] = as.numeric(x$p.value)
  }
}

data$p_cldc = ""

for(i in 1:nrow(data)){
  sd = c(sd(data[i,2:7]), sd(data[i,21:26]))
  if (any(sd == 0)){
    data$p_cldc[i] = 1
  }else{
    x = t.test(data[i,2:7], data[i,21:26])
    data$p_cldc[i] = as.numeric(x$p.value)
  }
}

data_6dpf = data

dat = data_6dpf[,1:26]
colnames(dat)[1] = "para"
colnames(dat)
dat = reshape2::melt(dat)
head(dat)
dat$group = gsub("_\\d", "", dat$variable)
dat$group = factor(dat$group, levels = c("control", "slc16a1a", "zgc158423", "cldn1"))
head(dat)

dats = dat[dat$para == unique(dat$para)[1],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Total length of vessels (μm)")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,40000))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    # axis.text.y = element_text(color = "black")
    axis.text.y = element_blank()
  )


dats = dat[dat$para == unique(dat$para)[2],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Number of vascular segments")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,800))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    # axis.text.y = element_text(color = "black")
    axis.text.y = element_blank()
  )


dats = dat[dat$para == unique(dat$para)[3],]
head(dats)
errbar = dats %>% group_by(group) %>% reframe(mean = mean(value),sd = sd(value))
errbar$lower = errbar$mean - errbar$sd
errbar$upper = errbar$mean + errbar$sd

ggplot(data = dats)+labs(x = NULL, y = "Vascular straightness")+
  geom_bar(aes(x = group, y = value, fill = group), stat = "summary", fun = "mean", 
           position = position_dodge(0.75),width = 0.75, color = "black", linewidth = 0.2)+
  geom_errorbar(data = errbar, aes(x = group,ymin = lower, ymax = upper, group = group), width = .2, 
                position=position_dodge(0.75), size = 0.75)+
  geom_point(aes(x = group, y = value, fill = group),size = 1,color = "black",
             position=position_jitterdodge(jitter.width = 0.1))+
  scale_y_continuous(expand = c(0,0,0,0), limits = c(0,1))+
  scale_fill_manual(values = c("#c7c9c8", "#89bbe5", "#f6d4a9", "#df9799"))+
  theme_classic()+
  theme(
    axis.text.x = element_text(angle =45, hjust = 1, vjust = 1, color = "black"),
    # axis.text.y = element_text(color = "black")
    axis.text.y = element_blank()
  )
