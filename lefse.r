#!/usr/bin/env Rscript
# Feng jiarong, shenzhen, 2021.03

args <- commandArgs(T)

cat("R script for lefse input formating")
cat("Example usage: Rscipts OTU_table_for_biom.txt sample.metaData.csv lefse.out")


if (length(args) != 3) {
  stop("Rscipts *.r OTU_table_for_biom.txt sample.metaData.csv lefse.out \n")
}

options(repos=structure(c(CRAN=“https://mirrors.tuna.tsinghua.edu.cn/CRAN/”)))
##--功能函数#---------
 # 功能函数
 ### 提取OTU表格
 vegan_otu <- function(physeq){
 OTU <- otu_table(physeq)
 if(taxa_are_rows(OTU)){
 OTU <- t(OTU)
 }
 return(as(OTU,"matrix"))
 }
 ### 添加OTU注释信息
 vegan_tax <- function(physeq){
 tax <- tax_table(physeq)

 return(as(tax,"matrix"))
 }
#

otu.table <- args[1]

tax.table <- 
map       <- read.delim(args[3], row.names =1)
out.path  <- args[4]

result = tryCatch({
    otu = read.delim(input.table,row.names = 1) %>% as.matrix()
    taxonomy <- reshape2::colsplit(otu$taxonomy,";",c("Kingdom","Phylum","Class","Order","Family","Genus","Species") )
    otu <- otu[!colnames(otu) %in% c("taxonomy") ]
}, error = function(e) {
    ps = readRDS(args[1])
}, finally = {
    cat("Please input a csv format table or phyloseq object")
}

if(! require("tidyverse") ) install.packages("tidyverse")
suppressMessages( library("tidyverse") )

if(! require("BiocManager") ) install.packages("BiocManager")
suppressMessages( library("BiocManager") )

if(! require("devtools") ) install.packages("devtools")
suppressMessages( library("devtools") )


if(! require("tidyverse") ) BiocManager::install("microbiomeViz")
suppressMessages( library("microbiomeViz"))

if(! require(tidyverse) ) install.packages("tidyverse")
suppressMessages( library("ggtree"))

if(! require(tidyverse) ) install.packages("tidyverse")
suppressMessages( library("tidyverse"))

if(! require(tidyverse) ) install.packages("tidyverse")
suppressMessages( library("phyloseq"))

if(! require(tidyverse) ) install.packages("tidyverse")
suppressMessages( library("dplyr"))
 library("tidyverse")

# 整理LEFSE分析的数据
# 这里如果大家熟悉phyloseq对象，直接用这个phyloseq对象就是了
# 导入数据

# 如果不熟悉phyloseq可以使用普通数据导入
# 导入otu表格

ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE),
 sample_data(map) ,
 tax_table(tax)
 # phy_tree(tree)
)

# 没有不要用全部OTU，按照序列数过滤，当然可以不过滤
# otu数量很多，所以选择一部分展示
ps = filter_taxa(ps, function(x) sum(x ) > 400 , TRUE)

# 提取注释文件
# 这一步是我我们准备Lefse物种注释表格整理的开始。
#-提取otu表和tax表格#--------
otu_table = as.data.frame(t(vegan_otu(ps)))
tax_table = (vegan_tax(ps))

# 整理物种注释文件
# 首先去除未注释出来的物种，这里为：Unassigned。
# 其次将物种不同级别之间的分隔符号修改为|。
# 再有将每一级的名称都加上前面级别的全部注释名称，形成结构。
# 最后按照不同级别数据合并otu表格，然后将注释文件和otu文件合并。
# 初步整理，去除Unassigned
tax_table[tax_table == "Unassigned"] = NA

tax_table = as.data.frame(
 tax_table
)
# 物种数据再整理
# 其次将物种不同级别之间的分隔符号修改为|。
# design = as.data.frame(sample_data(ps))
#合并otu表格和tax表格#---------
otu_tax = merge(otu_table,tax_table,by = "row.names",all = F)
# head(otu_tax)

 #--------添加门类标记，未注释的结果去掉#----
 otu_tax$Kingdom = paste("k",otu_tax$Kingdom,sep = "__")
 otu_tax$Kingdom[otu_tax$Kingdom == "k__NA"] =""

 otu_tax$Phylum = paste("|p",otu_tax$Phylum,sep = "__")
 otu_tax$Phylum[otu_tax$Phylum == "|p__NA"] = ""

 otu_tax$Class = paste("|c",otu_tax$Class,sep = "__")
 otu_tax$Class[otu_tax$Class == "|c__NA"] = ""

 otu_tax$Order = paste("|o",otu_tax$Order,sep = "__")
 otu_tax$Order[otu_tax$Order == "|o__NA"] = ""

 otu_tax$Family = paste("|f",otu_tax$Family,sep = "__")
 otu_tax$Family [otu_tax$Family == "|f__NA"] = ""

 otu_tax$Genus = paste("|g",otu_tax$Genus,sep = "__")
 otu_tax$Genus[ otu_tax$Genus == "|g__NA"] = ""

 otu_tax$Species = paste("|s",otu_tax$Species,sep = "__")
 otu_tax$Species[otu_tax$Species == "|s__NA"] = ""

 otu_tax$Row.names = paste("|t",otu_tax$Row.names,sep = "__")
 #合并得到结合全部门类的OTU名称#----
 OTU_name = paste(otu_tax$Kingdom,otu_tax$Phylum,otu_tax$Class,otu_tax$Order,otu_tax$Family,
 otu_tax$Genus,sep = "")

# 再有将每一级的名称都加上前面级别的全部注释名称，形成结构。这一步，我们去除了一些特殊的符号。主要是怕影响美观
 #----不同等级注释修改#----
 otu_tax$Phylum=paste(otu_tax$Kingdom,otu_tax$Phylum,sep = "")
 otu_tax$Class=paste(otu_tax$Phylum,otu_tax$Class,sep = "")
 otu_tax$Order=paste(otu_tax$Class,otu_tax$Order,sep = "")
 otu_tax$Family=paste(otu_tax$Order,otu_tax$Family,sep = "")
 otu_tax$Genus=paste(otu_tax$Family,otu_tax$Genus,sep = "")
 otu_tax$Species=paste(otu_tax$Genus,otu_tax$Species,sep = "")


 #替换两个括号等特殊符号#--------
 library("tidyverse")
 OTU_name = str_replace(OTU_name,"[(]","")
 OTU_name = str_replace(OTU_name,"[)]","")
 # as.character(OTU_name )
 # OTU_name = gsub("(","",OTU_name[311])
 # row.names(otu_tax) = OTU_name
 # #将otu表格和tax文件分离#-----
 # otu_table = otu_tax[2:(ncol(otu_table)+1)]
 # tax_table = otu_tax[(ncol(otu_table)+2):(ncol(otu_table)+2+6)]
 #
 # 再整理
 # 最后按照不同级别数据合并otu表格，然后将注释文件和otu文件合并。

 #-------将不同分类登记的也添加上去
 HA = otu_table
 # 按Kingdom合并
 grp <- otu_tax[rownames(otu_tax), "Kingdom", drop=F]
 merge=cbind(HA, grp)
 HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
 colnames(HA_Kingdom)[1]="Class"

 # 按Phylum合并
 grp <- otu_tax[rownames(otu_tax), "Phylum", drop=F]
 merge=cbind(HA, grp)
 HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
 colnames(HA_Phylum)[1]="Class"

 # 按Class合并
 grp <- otu_tax[rownames(otu_tax), "Class", drop=F]
 merge=cbind(HA, grp)
 HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
 colnames(HA_Class)[1]="Class"

 # 按Order合并
 grp <- otu_tax[rownames(otu_tax), "Order", drop=F]
 merge=cbind(HA, grp)
 HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
 colnames(HA_Order)[1]="Class"

 # 按Family合并
 grp <- otu_tax[rownames(otu_tax), "Family", drop=F]
 merge=cbind(HA, grp)
 HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
 colnames(HA_Family)[1]="Class"

 # 按Genus合并
 grp <- otu_tax[rownames(otu_tax), "Genus", drop=F]
 merge=cbind(HA, grp)
 HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
 colnames(HA_Genus)[1]="Class"
 # colnames(otu_tax)
 # # 按Species合并
 # grp <- otu_tax[rownames(otu_tax), "Species", drop=F]
 # merge=cbind(HA, grp)
 # HA_Species = merge %>% group_by(Species) %>% summarise_all(sum)
 # colnames(HA_Species)[1]="Class"
 ## OTU水平
 # merge=cbind(HA, OTU_name)
 # head(otu_table)
 # HA_OTU = otu_table

 # HA_OTU = data.frame(Class = row.names(otu_table),otu_table)
 # colnames(HA_OTU )
 # head(HA_OTU)
 # 合并6个分类级
 all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)
 dim(all)


# head(OTU_name)
# 因为我们去除了未注释出来的，所以会有重复
# 举个例子：
# 1|2|3|4|5|6|7

# 这里只需要简单去除就可以了。我们还构建了phyloseq对象，也就是说基于这个群落我们做其他分析也是是十分方便。
# 去除重复
all = distinct(all, Class, .keep_all = TRUE)
dim(all)

#-------------lefse分析构建
all1 = all
row.names(all1) = all1$Class
all1$Class = NULL
all1 = as.matrix(all1)
# head(all1)

#-构建phylose对象
ps_G_graphlan = phyloseq(otu_table(as.matrix(all1),taxa_are_rows = TRUE),
 sample_data(ps))
ps_G_graphlan

#----提取OTU表格
otu_table = as.data.frame((vegan_otu(ps_G_graphlan)))
otu_table[otu_table==0] <- 1

# row.names(otu_table)

# head(design)

design = as.data.frame(sample_data(ps_G_graphlan))
# 如果我们要使用功能python那个脚本运行，保存就可以使用以下代码
# 这里我就不做运行了，之前的帖子相信大家都看过了，我是做了对比的，发现算法都是一致。
# filename = paste(lefsepath,"/LEFSE_to_run_G_level.txt",sep = "")
# write.table(otu_table, filename,append = F, quote = F,col.names= F,sep = "\t")
# 文件预处理
# format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res -l 4
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf

# # 注意这里 –c用来指定分组信息-u 1指定样品信息
# #文件分析,这里-l设置LDA阈值，默认为2，我们使用4 会更加严格
# ~/src/nsegata-lefse/run_lefse.py pri_lefse.in pri_lefse_2.res -l 4
# #柱状图绘制
# plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
# #树状图绘制
# plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
# #做每个差异的柱状图
# mkdir biomarkers_raw_images
## plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/

#-------------R语言进行lefse分析
# ps_sub <- subset_samples(ps_G_graphlan,Group %in% c("July-Orchard","July-cropland"));ps_sub
ps_sub = ps_G_graphlan

otu_table = as.data.frame((vegan_otu(ps_sub)))
otu_table[otu_table==0] <- 1

# head(otu_table)

map = as.data.frame(sample_data(ps_sub))
otu = (otu_table)
claslbl= map$Group

# 随机数种子
set.seed(12345);
#KW rank sum test
dat3t = otu_table
# head(otu)
# dat3t = as.matrix(dat3t)
# str(dat3t)
rawpvalues <- apply(dat3t, 2, function(x) kruskal.test(x, claslbl)$p.value);
#--得到计算后得到的p值
ord.inx <- order(rawpvalues);
rawpvalues <- rawpvalues[ord.inx];
clapvalues <- p.adjust(rawpvalues, method ="fdr");
# p.adjust
dat3t <- dat3t[,ord.inx];
dim(dat3t)

wil_datadf <- as.data.frame(dat3t);
# head(wil_datadf)
#if no subclass within classes then no wilcoxon rank sum test
#Linear Discriminant analysis (LDA)
library( MASS)
ldares <- lda(claslbl~ .,data = wil_datadf);
# ldares
ldamean <- as.data.frame(t(ldares$means));
# ldamean
class_no <<- length(unique(claslbl));
ldamean$max <- apply(ldamean[,1:class_no],1,max);
ldamean$min <- apply(ldamean[,1:class_no],1,min);

#---计算LDA
ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);

# head(ldamean)

a = rep("A",length(ldamean$max))
i = 1
for (i in 1:length(ldamean$max)) {
 name =colnames(ldamean[,1:class_no])
 a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
}
ldamean$class = a
# head(ldamean)

ldamean$Pvalues <- signif(rawpvalues,digits=5);
ldamean$FDR <- signif(clapvalues,digits=5);
resTable <- ldamean;
# head(ldamean)
# it seems lda add ` around names containing dash "-", need to strip this off
rawNms <- rownames(resTable);
rownames(resTable) <- gsub("`", '', rawNms);
# head(resTable)
# if(pvalOpt == "raw"){
# de.Num <- sum(rawpvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
# }else{
p.lvl =0.05
lda.lvl = 2
de.Num <- sum(clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
# }

if(de.Num == 0){
 current.msg <<- "No significant features were identified with given criteria.";
}else{
 current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
}
current.msg
# sort by p value
ord.inx <- order(resTable$Pvalues, resTable$LDAscore);
resTable <- resTable[ord.inx, ,drop=FALSE];
#p-values column to appear first; then FDR and then others
resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];

write.csv(resTable,"cs.csv")
# R语言绘制物种分类树
# 这部分我首先进行了颜色的映射数据整理，然后使用microbiomeviz做了圈图的输出。
taxtree = resTable[clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl,]

delimiter = "\\|"
tax_split <- strsplit(row.names(taxtree), delimiter)
row.names(taxtree)<- vapply(tax_split, tail, n = 1, "")

# head(taxtree)

#-提取所需要的颜色
colour = c('darkgreen','red',"blue")

selececol = colour[1:length(levels(as.factor(taxtree$class)))]

names(selececol) = levels(as.factor(taxtree$class))
A = rep("a",length(row.names(taxtree)))

i = 1
for (i in 1:length(row.names(taxtree))) {
 A[i] = selececol [taxtree$class[i]]
}

# A

lefse_lists = data.frame(node=row.names(taxtree),
 color=A,
 stringsAsFactors = FALSE
)

# str(all)

## 计算均值用于呈现结点大小
dat <- data.frame(V1=all[,1], V2=rowMeans(all[,-1]), stringsAsFactors = FALSE)

# head(dat)
# dim(dat)
# write.csv(dat,"./tree_tax.csv")
# dat$V2 = NULL
# 用物种和丰度生成树骨架
tr <- parseMetaphlanTSV(dat, node.size.offset=2, node.size.scale=0.8)

# tr

# 构造树#------

p =tree.backbone(tr, size=1,layout= 'circular')
p
# ?tree.backbone

# 注释树
p <- clade.anno(p, lefse_lists, alpha=0.3)
p

ggsave("./cs.pdf",p,width = 10,height = 10)