####转成scanpy格式####
#法1 SeuratDisk  常规转换比较慢
#R转为Python
library(qs)
qsave(sce_allcell,file= 'sce_allcell_harmony.qs')
sce_allcell<-qread('sce_allcell_harmony.qs')
library(SeuratDisk)
SaveH5Seurat(sce_allcell, filename = "data.h5Seurat",overwrite = TRUE) 
Convert("data.h5Seurat", dest = "h5ad",overwrite = TRUE)
#Python转为R--搞不好--转成法4  zellkonverter
Convert("/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5ad", dest = "h5seurat", overwrite = TRUE) 
sce_Epithelial_T2 <- LoadH5Seurat("/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5seurat",meta.data = FALSE, misc = FALSE)

#法2  zellkonverter  Python转R
BiocManager::install("zellkonverter")
library(zellkonverter)
sce1=readH5AD("Epithelial_T_afterumap.h5ad", verbose = TRUE)
sce_Epithelial_T2 <- as.Seurat(sce1, counts = "X", data = NULL)
qsave(sce_Epithelial_T2,file= 'sce_Epithelial_T2.qs') #目前07-24最新

#qs文件
##allcell
sce_allcell<-qread('/home/wling_32/NGmetadata/sce_allcell_harmony.qs')
##原版3w基因经过infercnv筛选的纯Tumor上皮  
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
##经过python筛选2000高辨基因的纯Tumor上皮 
sce_Epithelial_T <- qread(file = "/home/wling_32/jupyter_home/NGmetadata/sce_Epithelial_T2.qs")

#h5ad文件
##allcell
"/home/wling_32NGmetadata/Rawdata/allcell.h5ad"
#原版3w基因上皮（未经infercnv） 
"/home/wling_32/jupyter_home/NGmetadata/Epithelial.h5ad"
##原版3w基因经过infercnv筛选的纯Tumor上皮  
"/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T.h5ad"
##经过python筛选2000高辨基因的纯Tumor上皮 
"/home/wling_32/jupyter_home/NGmetadata/cNMF/Epithelial_T_afterumap.h5ad"
##最新版加了fr聚类的纯Tumor上皮
"/home/wling_32/jupyter_home/NGmetadata/Pseudotime/Epithelial_T_afterumap.h5ad"


####infercnvpy后续分析####
sce_Epithelial <- qread(file = "~/NGmetadata/InferCNV/sce_Epithelial.qs")
qsave(sce_Epithelial,file= 'sce_Epithelial.qs')

##把python的obs导出到excel中，然后添加到seurat里
epithelial_obs <- read_excel("~/jupyter_home/NGmetadata/epithelial_obs.xlsx")
identical(epithelial_obs$index,colnames(sce_Epithelial))
sce_Epithelial <- AddMetaData(sce_Epithelial,epithelial_obs$Sample_type,col.name = "Sample_type")
sce_Epithelial <- AddMetaData(sce_Epithelial,epithelial_obs$cnv_leiden,col.name = "cnv_leiden")
sce_Epithelial <- AddMetaData(sce_Epithelial,epithelial_obs$cnv_score,col.name = "cnv_score")
range(sce_Epithelial$cnv_score)


##按照cnvscore区间分析---数据不合适
start_value <- 0.00284852
end_value <- 0.06112092
segments <- seq(from = start_value, to = end_value, length.out = 26)
sce_Epithelial$cnv_score_segments <- cut(sce_Epithelial$cnv_score, 
                                         breaks = segments, 
                                         labels = seq(from=1, to=25), 
                                         right = FALSE)

table(sce_Epithelial$cnv_score_segments)

ggplot(data=sce_Epithelial@meta.data, mapping=aes(x=cnv_score_segments,fill=Sample_type))+
  geom_bar(stat="count",width=0.5,position='stack')+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  geom_text(stat='count',aes(label=..count..), color="white", size=3.5,position=position_stack(0.5))+
  theme_minimal()

##按照cnvscore和leiden聚类看图-用了
CNV_score<- sce_Epithelial@meta.data[,c("Sample_type","cnv_leiden","cnv_score")]
#CNV_score$cnv_leiden <- as.factor(CNV_score$cnv_leiden,levels=c(0:48))

ggplot(CNV_score,aes(cnv_leiden,cnv_score))+
  geom_violin(aes(fill=cnv_leiden))+
  scale_fill_manual(values = rainbow(50) )+
  theme_bw()
#结果没呈现小提琴图  因为cnvscore都是几种值的重复  最后能看出那些点 0 1 25最低

sce_Epithelial$cnv_status <- ifelse(sce_Epithelial$Sample_type=="Tumor cell" & sce_Epithelial$cnv_leiden %in% c(0,1,25),"Normal",
                                    ifelse(sce_Epithelial$Sample_type=="Normal cell","Normal","Tumor"))

table(sce_Epithelial$cnv_status)

library(qs)
library(Seurat)
qsave(sce_Epithelial,file= '/home/wling_32/NGmetadata/sce_Epithelial.qs')
sce_Epithelial <- qread(file = "~/NGmetadata/sce_Epithelial.qs")

sce_Epithelial_T <- subset(sce_Epithelial,cnv_status=="Tumor")
qsave(sce_Epithelial_T,file= '/home/wling_32/NGmetadata/sce_Epithelial_T.qs')


####cNMF后续分析-runGSEA####
library(qs)
library(Seurat)
qsave(sce_Epithelial,file= '/home/wling_32/NGmetadata/sce_Epithelial.qs')
sce_Epithelial <- qread(file = "~/NGmetadata/sce_Epithelial.qs")
sce_Epithelial_T <- qread(file = "~/NGmetadata/sce_Epithelial_T.qs")

sce_Epithelial_T <- subset(sce_Epithelial,cnv_status=="Tumor")
qsave(sce_Epithelial_T,file= '/home/wling_32/NGmetadata/sce_Epithelial_T.qs')

topgenes <- read_excel("~/jupyter_home/NGmetadata/cNMF/topgenes.xlsx")
topgenes <- topgenes[,c(2:8)]
colnames(topgenes) <- c("cNMF_Cluster1","cNMF_Cluster2","cNMF_Cluster3","cNMF_Cluster4","cNMF_Cluster5",
                        "cNMF_Cluster6","cNMF_Cluster7")

cNMF_clusterscore <- read.delim("~/jupyter_home/NGmetadata/cNMF/Epithelial_T_cNMF/Epithelial_T_cNMF.usages.k_7.dt_0_15.consensus.txt")
rownames(cNMF_clusterscore) <- cNMF_clusterscore$index
cNMF_clusterscore <- cNMF_clusterscore[,c(-1)]
colnames(cNMF_clusterscore) <- c("cNMF_Cluster1","cNMF_Cluster2","cNMF_Cluster3","cNMF_Cluster4","cNMF_Cluster5",
                                 "cNMF_Cluster6","cNMF_Cluster7")
identical(rownames(cNMF_clusterscore),colnames(sce_Epithelial_T))
sce_Epithelial_T <- AddMetaData(sce_Epithelial_T,metadata = cNMF_clusterscore)

##runGSEA确定每个cNMFcluster的特征
install.packages("GeneNMF")
install.packages("GeneNMF", version = "0.2.0")
library(GeneNMF)
library(msigdbr)
library(fgsea)
top_p <- lapply(topgenes, function(program) {
  runGSEA(program, universe=rownames(sce_Epithelial_T), category = "H")
})


#去给每个细胞归属cNMF_cluster
metadata <- sce_Epithelial_T@meta.data
metadata$cNMF_Cluster <- apply(metadata[, c("cNMF_Cluster1", "cNMF_Cluster2", "cNMF_Cluster3", "cNMF_Cluster4", "cNMF_Cluster5", "cNMF_Cluster6", "cNMF_Cluster7")], 1, function(row) {
  names(which.max(row))
})
table(metadata$cNMF_Cluster)
identical(rownames(metadata),colnames(sce_Epithelial_T))
sce_Epithelial_T$cNMF_Cluster <- metadata$cNMF_Cluster

DimPlot(sce_Epithelial_T,reduction = "umap", group.by = 'cNMF_Cluster')  + 
  NoAxes() + 
  theme(plot.title = element_blank())


head(top_p_H$cNMF_Cluster7,10)

C2 <- top_p_C2$cNMF_Cluster2




####cNMF后续-GSVA####
setwd("cNMF")
library(clusterProfiler)
library(Seurat)
library(GSVA)
library(msigdbr)
library(qs)

sce_Epithelial_T <- qread(file = "sce_Epithelial_T2.qs")
Idents(sce_Epithelial_T) <- "cNMF_Cluster"  
expr <- AverageExpression(sce_Epithelial_T, assays = "originalexp", slot = "data")[[1]] 
expr <- expr[rowSums(expr)>0,]  #选取非零基因 
expr <- as.matrix(expr) 
head(expr)

# Prepare gene set
#reactome/kegg/go/hallmarkers
geneset.c2.cp <- read.gmt("/home/wling_32/jupyter_home/NGmetadata/cNMF/GSVA/c2.cp.v2023.2.Hs.symbols.gmt")  
geneset.c5.all <- read.gmt("/home/wling_32/jupyter_home/NGmetadata/cNMF/GSVA/c5.all.v2023.2.Hs.symbols.gmt")  
geneset.h.all <- read.gmt("/home/wling_32/jupyter_home/NGmetadata/cNMF/GSVA/h.all.v2023.2.Hs.symbols.gmt")  

geneset <- rbind(geneset.c2.cp,geneset.c5.all,geneset.h.all)
genesets = split(geneset$gene, geneset$term)

##分组gsva
tumorcell_gsva.res <- gsva(expr, genesets, method="gsva")  
tumorcell_gsva.df <- data.frame(Genesets=rownames(tumorcell_gsva.res), tumorcell_gsva.res, check.names = F)
save(tumorcell_gsva.res,tumorcell_gsva.df, file="tumorcell_gsva.Rdata") 


##绘图方法1
pathway <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
             "HALLMARK_P53_PATHWAY",
             "KEGG_TGF_BETA_SIGNALING_PATHWAY",  #C1
             "GOBP_REGULATION_OF_COLLAGEN_FIBRIL_ORGANIZATION",#ecm
             "REACTOME_DNA_METHYLATION","GOBP_MRNA_METHYLATION",#METHYLATION
             "REACTOME_METABOLISM_OF_RNA","GOBP_DETOXIFICATION_OF_NITROGEN_COMPOUND",
             "GOBP_GLUTAMATE_BIOSYNTHETIC_PROCESS","GOBP_GLYCOLYTIC_PROCESS_THROUGH_FRUCTOSE_6_PHOSPHATE",#metabolism
             "GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL","GOBP_NATURAL_KILLER_CELL_CYTOKINE_PRODUCTION",
             "GOBP_CD8_POSITIVE_ALPHA_BETA_T_CELL_DIFFERENTIATION",#immue
             "GOBP_REGULATION_OF_RESPONSE_TO_DRUG","REACTOME_PD_1_SIGNALING" #drug
)
gsva.res <- as.matrix(tumorcell_gsva.res[pathway,])

annotation_col <- data.frame(Cell= c("Stromal",
                                     rep("Myeloids",2),
                                     rep("T cell",3),
                                     rep("B cell",3)
))
annotation_row <- data.frame(Pathway= c(rep("ECM-related",4),
                                        rep("Methylation-related",2),
                                        rep("Metabolism-related",4),
                                        rep("Immune-related",3),
                                        rep("Thearpy-related",2)))
row.names(annotation_col) <- colnames(gsva.res)
row.names(annotation_row) <- pathway


cellcolor <- c("#e85d14", "#993c00", "#f5c563","#006f80")
names(cellcolor) <- c("Stromal","Myeloids","T cell","B cell") #类型颜色
pathwaycolor <- c("#3C5488","#00A087","#4DBBD5",'#E64B35',"#EB9412") 
names(pathwaycolor) <- c("ECM-related","Methylation-related","Metabolism-related","Immune-related","Thearpy-related") #类型颜色
ann_colors <- list(Cell=cellcolor, Pathway=pathwaycolor ) #颜色设置

pheatmap::pheatmap(gsva.res, show_colnames = T,  
                   cluster_rows = F,cluster_cols = F,
                   #scale = "row",
                   angle_col = "90",  
                   fontsize_row = 8,
                   fontsize_col = 10,
                   annotation_colors = ann_colors,  #图例颜色注释信息
                   annotation_col = annotation_col,
                   annotation_row = annotation_row,
                   row_split = annotation_row$Pathway,
                   cellwidth = 12,cellheight = 12,
                   column_split = annotation_col$Cell,border_color = "black",
                   color = colorRampPalette(c("#106894", "white", "#BA0D0D"))(50))


##绘图方法2
library(dplyr)
library(tibble)
library(paletteer)
library(ComplexHeatmap)
library(circlize)

new_row_names <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_P53_PATHWAY",
  "REACTOME_SIGNALING_BY_EGFR",#C1
  "GOBP_MIRNA_PROCESSING",
  "GOBP_POSITIVE_REGULATION_OF_MRNA_PROCESSING", #C2
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "REACTOME_COLLAGEN_FORMATION", #C3
  "GOBP_CELL_CYCLE",
  "HALLMARK_G2M_CHECKPOINT",
  "GOBP_MEIOTIC_CELL_CYCLE_PHASE_TRANSITION", #C4
  "GOBP_REGULATION_OF_TRANSLATIONAL_INITIATION_IN_RESPONSE_TO_STRESS",
  "GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE",
  "GOBP_RESPONSE_TO_LAMINAR_FLUID_SHEAR_STRESS", #C5
  "KEGG_FATTY_ACID_METABOLISM",
  "REACTOME_METABOLISM_OF_PORPHYRINS",
  "KEGG_NITROGEN_METABOLISM",   #C6
  "GOBP_GLUTAMATE_SECRETION_NEUROTRANSMISSION",
  "GOBP_ENTEROENDOCRINE_CELL_DIFFERENTIATION" #C7
)

new_names <- c( "Oncogene-drive", "Oncogene-drive", "Oncogene-drive", 
                "Mixed", "Mixed", 
                "p-EMT", "p-EMT",
                "Cycling", "Cycling", "Cycling", 
                "Stress", "Stress", "Stress",  
                "Metabolism","Metabolism","Metabolism",
                "Secretory","Secretory") 

new_names<-factor(x = new_names, levels =unique(c( "Oncogene-drive", "Oncogene-drive", "Oncogene-drive", 
                                                   "Mixed", "Mixed", 
                                                   "p-EMT", "p-EMT",
                                                   "Cycling", "Cycling", "Cycling", 
                                                   "Stress", "Stress", "Stress",  
                                                   "Metabolism","Metabolism","Metabolism",
                                                   "Secretory","Secretory")) )
column_names <- paste0("cNMF_Cluster", 1:7) 
set.seed(123)  # 设置随机种子以确保结果可重复
new_data <- tumorcell_gsva.res[new_row_names,]

# 创建数据框并赋予行名和列名
new_df <- as.data.frame(new_data)
rownames(new_df) <- new_row_names
colnames(new_df) <- column_names

# # 定义第一个左侧注释，带有背景框和文本
left_annot <- rowAnnotation(
  text = anno_text(new_row_names, gp = gpar(col = c(
    "black",
    "black",
    "black",
    "black",
    "black"
  ),
  border = "black"),
  #                       location = 0.5,just = "center",
  width = unit(12, "cm")
  ),
  block2 = anno_block(
    gp = gpar(fill = c("#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB")),
    labels_rot=0,
    labels = unique(new_names),
    width = max_text_width(unique(new_names))*1.05,
    labels_gp = gpar(col = "black", fontsize = 16)
  )
)

#cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y,
            width = width, 
            height = height, 
            gp = gpar(col = "black", fill = "white"))
  grid.circle(x =x, y = y, r =  unit(3 * data[i,j], "mm"), 
              gp = gpar(fill = fill, col = NA))
}

options(repr.plot.width = 15, repr.plot.height = 15)
ht<-Heatmap(
  new_data, 
  name = "matrix",
  cluster_rows = FALSE, 
  cluster_columns = FALSE,
  show_heatmap_legend=FALSE,
  top_annotation = HeatmapAnnotation(barplot = anno_empty(height = unit(1.5, "cm"))),
  row_title = NULL, 
  col=c("#004b8f", "#FFFFFF", "#B80909"),
  show_row_names=FALSE,
  column_names_rot=45,
  width = unit(7, "cm"),
  height = unit(18, "cm"),
  row_split=factor(new_names),
  heatmap_legend_param = list(title = "scaled score",
                              title_gp = gpar(fontsize = 16), # 设置标题字体大小和粗体
                              labels_gp = gpar(fontsize = 16)),
  column_names_gp = gpar(fontsize = 18),
  right_annotation = left_annot, # 组合两个左侧注释
  #cell_fun = cell_fun
)

ht
breaks <- c(-1, 0, 1)
colors <- c("#004b8f", "#FFFFFF", "#B80909")
color_mapping <- colorRamp2(breaks, colors)

custom_legend2 <- Legend(title = "Scale",
                         at = breaks,
                         labels = c("-1", "0", "1"),
                         col_fun = color_mapping,
                         legend_width = unit(4, "cm"),
                         legend_height = unit(2.5, "cm"),
                         title_gp = gpar(fontsize = 16), # 设置标题字体大小和粗体
                         labels_gp = gpar(fontsize = 16  ),
                         direction = "horizontal")
draw(ht, annotation_legend_list = list(custom_legend2))


# add cNMF_Cluster ratio
prop.table(table(sce_Epithelial_T$cNMF_Cluster))
ratio=c(0.354317074 ,  0.226496523  , 0.128140235  , 0.095372566 ,  0.109944854 ,  0.082318779 ,  0.003409969)
decorate_annotation("barplot", {
  #     pushViewport(viewport(xscale = c(0:1), yscale = c(0, max(p)*1.1)))
  grid.rect(x = 1:7, y = 0, width = 0.5, height = ratio, just = "bottom",
            gp = gpar(fill = c("#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB")), default.units = "native")
  grid.yaxis()
  grid.text("cNMF_Cluster ratio", x = unit(-1.5, "cm"),rot = 0, just = "bottom")
  popViewport()
})


####AUCcell基因评分####
library(AUCell)
library(clusterProfiler)
load("~/Glycolysis/单细胞分cluster/sceT.Rdata")
sce_allcell<-qread('/home/wling_32/NGmetadata/sce_allcell_harmony.qs')
cells_rankings <- qread('/home/wling_32/NGmetadata/基因集评分/cells_rankings_sce_allcell.qs')

cells_rankings <- AUCell_buildRankings(sce_allcell@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(sceT@assays$RNA@data)
save(cells_rankings,file = "AUCcell_cells_rankings_sceT.rda")
qsave(cells_rankings,file= 'cells_rankings_sce_allcell.qs')

gene <- intersect(rownames(sceT@assays$RNA),gene)
markers <- list()
markers$neddylation <- gene
cells_AUC <- AUCell_calcAUC(markers, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
##可视化
geneSet <- "neddylation"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sce_allcell$AUC  <- aucs
library(ggraph)
ggplot(data.frame(sce_allcell@meta.data, sce_allcell@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=AUC)
) + geom_point( size=1.5
) + scale_color_viridis(option="A")  + theme_light(base_size = 15)+labs(title = "Neddylation")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  theme(plot.title = element_text(hjust = 0.5))

###正常和肿瘤上皮Neddylation评分绘图
library(ggplot2)
library(readxl)
epithelial_obs <- read_excel("epithelial_obs.xlsx")
df <- epithelial_obs

library("ggpubr")
my_comparisons=list(c("Normal cell","Tumor cell"))
ggviolin(df, x="Sample_type", y="Neddylation_aucell", fill = "Sample_type",
         palette = c('nejm'),
         add = "boxplot",
         add.params = list(fill="white"),
         title ="Neddylation AUCell score",
         xlab = "Sample type",
         ylab = "Score",
)+
  scale_y_continuous(limits = c(0.05, 0.35))+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")+
  theme_bw()+theme(plot.margin = margin(t = 2,  # 顶部边缘距离
                                        r = 3,  # 右边边缘距离
                                        b = 2,  # 底部边缘距离
                                        l = 3,  # 左边边缘距离
                                        unit = "cm"))#设置单位为cm

##肿瘤cNMF亚群nedd评分
library(ggpubr)
identical(Epithelial_T_obs$index,Epithelial_T_nedd$index)
Epithelial_T_obs$Neddylation_aucell <- Epithelial_T_nedd$Neddylation_aucell
Epithelial_T_obs$cNMF_Cluster <- factor(Epithelial_T_obs$cNMF_Cluster,levels=c("cNMF_Cluster4",
                                                                               "cNMF_Cluster1",
                                                                               "cNMF_Cluster2",
                                                                               "cNMF_Cluster3",
                                                                               "cNMF_Cluster5",
                                                                               "cNMF_Cluster6",
                                                                               "cNMF_Cluster7"))
p2 <- ggboxplot(Epithelial_T_obs, x="cNMF_Cluster", y="Neddylation_aucell", width = 0.6, 
                color = "black",#轮廓颜色
                fill="cNMF_Cluster",#填充
                #add = "boxplot",
                #add.params = list(fill="white"),"#D55640","#65BA8E","#F3AE63","#7BBDDB","#73558B","#E2BECB","#DBBD99"
                palette = c("#65BA8E","#F3AE63","#D55640","#DBBD99","#E2BECB","#73558B","#7BBDDB"),#"#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB"
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=0.5, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") #图例放右边 
###指定组比较
my_comparisons <- list(c("cNMF_Cluster4", "cNMF_Cluster1"), c("cNMF_Cluster4", "cNMF_Cluster2"),c("cNMF_Cluster4", "cNMF_Cluster7"),
                       c("cNMF_Cluster4", "cNMF_Cluster6"),c("cNMF_Cluster4", "cNMF_Cluster5"),c("cNMF_Cluster4", "cNMF_Cluster3"))
p2+stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")


####拟时序分析-cytotrace2.0####
setwd("Pseudotime")
devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r") 
library(CytoTRACE2)
library(tidyverse)
library(Seurat)
library(qs)
##输入原版3w基因的sce_epi_T
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
cytotrace2_result_sce <- cytotrace2(sce_Epithelial_T, 
                                    is_seurat = TRUE, 
                                    slot_type = "counts", 
                                    species = 'human',
                                    seed = 1234)
cytotrace2_result_sce
qsave(cytotrace2_result_sce,file="cytotrace2_result_sce.qs")

##可视化
# making an annotation dataframe that matches input requirements for plotData function
library(readxl)
Epithelial_T_obs <- read_excel("~/jupyter_home/NGmetadata/cNMF/Epithelial_T_obs.xlsx")
identical(colnames(sce_Epithelial_T),Epithelial_T_obs$index)
sce_Epithelial_T$cNMF_Cluster <- Epithelial_T_obs$cNMF_Cluster
annotation <- data.frame(phenotype = sce_Epithelial_T@meta.data$cNMF_Cluster) %>% 
  set_rownames(., colnames(sce_Epithelial_T))

# plotting
plots <- plotData(cytotrace2_result = cytotrace2_result_sce, 
                  annotation = annotation, 
                  is_seurat = TRUE)
#save(plots,file="cytotrace2_plots.rda")太大了不存了
# 绘制CytoTRACE2_Potency的umap图
p1 <- plots$CytoTRACE2_UMAP
# 绘制CytoTRACE2_Potency的umap图
p2 <- plots$CytoTRACE2_Potency_UMAP
# 绘制CytoTRACE2_Relative的umap图 ，v1 
p3 <- plots$CytoTRACE2_Relative_UMAP 
# 绘制各细胞类型CytoTRACE2_Score的箱线图
p4 <- plots$CytoTRACE2_Boxplot_byPheno

library(patchwork)
(p1+p2+p3+p4) + plot_layout(ncol = 2)
p1
cytotracedata <- cytotrace2_result_sce@meta.data
write.csv(cytotracedata,file="cytotrace_metadata.csv")

##量化差异绘箱线图
identical(Epithelial_T_obs$index,rownames(cytotracedata))
cytotracedata$cNMF_Cluster <- Epithelial_T_obs$cNMF_Cluster
cytotracedata$cNMF_Cluster <- factor(cytotracedata$cNMF_Cluster,levels = c("cNMF_Cluster2",
                                                                           "cNMF_Cluster4",
                                                                           "cNMF_Cluster1",
                                                                           "cNMF_Cluster7",
                                                                           "cNMF_Cluster6",
                                                                           "cNMF_Cluster5",
                                                                           "cNMF_Cluster3"))
library(ggpubr)
p1 <- ggboxplot(cytotracedata, x="cNMF_Cluster", y="CytoTRACE2_Relative", width = 0.6, 
                color = "black",#轮廓颜色
                fill="cNMF_Cluster",#填充
                #add = "boxplot",
                #add.params = list(fill="white"),"#D55640","#65BA8E","#F3AE63","#7BBDDB","#73558B","#E2BECB","#DBBD99"
                palette = c("#D55640","#65BA8E","#F3AE63","#7BBDDB","#73558B","#E2BECB","#DBBD99"),#"#F3AE63","#D55640","#DBBD99","#65BA8E","#E2BECB","#73558B","#7BBDDB"
                xlab = F, #不显示x轴的标签
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                size=0.5, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") #图例放右边 
###指定组比较
my_comparisons <- list(c("cNMF_Cluster2", "cNMF_Cluster1"), c("cNMF_Cluster2", "cNMF_Cluster4"),c("cNMF_Cluster2", "cNMF_Cluster7"),
                       c("cNMF_Cluster2", "cNMF_Cluster6"),c("cNMF_Cluster2", "cNMF_Cluster5"),c("cNMF_Cluster2", "cNMF_Cluster3"))
p1+stat_compare_means(comparisons = my_comparisons,
                      method = "wilcox.test")



####拟时序分析-PAGA后续分析）####
##绘图--定义早期晚期cNMF_Cluster
#读入PAGA网络细胞坐标节点数据，修改一下列名
cell_emb <- read.csv("clust_embedding.csv", header = T)
colnames(cell_emb)[1] <- "cNMF_Cluster"
cell_emb$cNMF_Cluster <- as.factor(cell_emb$cNMF_Cluster)
cell_emb$cNMF_Cluster <- paste0("cNMF_Cluster",cell_emb$cNMF_Cluster)

#读入PAGA网络互作连接信息，修改一下列名
cell_edge <- read.csv("edge.csv", header = T)
colnames(cell_edge)[1] <- "From"

#加载R包作图
library(igraph)
library(ggraph)
library(tidyverse)
library(tidygraph)
#转化数据，我们使用的是ggraph作图，所以数据需要转化为他识别的格式
graph <- as_tbl_graph(cell_edge) %>% 
  mutate(interaction = centrality_degree(mode='out'))

#导出布局，进行替换
#导出布局之后，其实就是我们作图的数据。这里布局自动生成了每个细胞节点的坐标，我们更具细胞类型替换它在PAGA里面的坐标信息
gglayout <- create_layout(graph, layout = "kk")
gglayout
gglayout$x <- cell_emb$X
gglayout$y <- cell_emb$Y

#添加节点大小数据，细胞类型细胞数的多少用来表示节点大小
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
cell_number <- as.data.frame(table(sce_Epithelial_T$cNMF_Cluster))
gglayout$Freq <- cell_number$Freq


#作图，这里我们使用ggunchull包添加底色。使用annotate添加注释文字
# devtools::install_github("sajuukLyu/ggunchull", type = "source")
library(ggunchull)
ggraph(gglayout) + 
  #stat_unchull(data=gglayout[gglayout$name %in% c("cNMF_Cluster7","cNMF_Cluster6","cNMF_Cluster5","cNMF_Cluster3"),],
  #             aes(x,y),
  #             fill="#B35656",
  #             show.legend = F,
  #             nbin = 50, 
  #             nsm = 30,
  #             qval = 0.8,
  #             sfac = 14)+
  #stat_unchull(data=gglayout[gglayout$name %in% c("cNMF_Cluster2","cNMF_Cluster4","cNMF_Cluster1"),],
  #             aes(x,y),
  #             fill="#D4EAC4",
  #             show.legend = F,
  #             nbin = 50, 
  #             nsm = 30,
  #             qval = 0.8,
  #             sfac = 12)+
  geom_edge_link(aes(edge_width=Conn),color='black') + 
  geom_node_point(aes(size = Freq,color=name),shape=16) +
  geom_node_text(aes(label=name),size=3) +
  scale_edge_width(range=c(0,4))+
  scale_size(range = c(0,15))+
  theme_classic()+
  NoLegend()+
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  annotate("text", label="Late_cNMF_Cluster",x=2,y=0.5, color='#F97B72', size=6, fontface="bold")+
  annotate("text", label="Early_cNMF_Cluster",x=-0.5,y=-2, color='#87C55F', size=6, fontface="bold")



####高低nedd分组分析-KEGG/GO####
####KEGG/GO分析
setwd("Neddgroup_Analysis")
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(Seurat)
library(qs)
library(clusterProfiler)
library(dplyr)
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
library(readxl)
Epithelial_T_obs <- read_excel("~/jupyter_home/NGmetadata/cNMF/Epithelial_T_obs.xlsx")
Epithelial_T_nedd <- read_excel("~/jupyter_home/NGmetadata/Neddylation_AUcell/Epithelial_T_nedd.xlsx")
cytotrace_metadata <- read.csv("~/jupyter_home/NGmetadata/Pseudotime/cytotrace_metadata.csv", row.names=1)

identical(Epithelial_T_obs$index,Epithelial_T_nedd$index)
identical(Epithelial_T_obs$index,rownames(cytotrace_metadata))

Epithelial_T_obs$Neddylation_aucell <- Epithelial_T_nedd$Neddylation_aucell
Epithelial_T_obs <- cbind(Epithelial_T_obs,cytotrace_metadata[,c(24:28)])

write.csv(Epithelial_T_obs,file="/home/wling_32/jupyter_home/NGmetadata/Epithelial_T_obs.csv")

###区分nedd高低组
Epithelial_T_obs$Neddylation_Group <- ifelse(Epithelial_T_obs$Neddylation_aucell>=median(Epithelial_T_obs$Neddylation_aucell),"High","Low")
table(Epithelial_T_obs$Neddylation_Group)
identical(colnames(sce_Epithelial_T),Epithelial_T_obs$index)
sce_Epithelial_T$Neddylation_Group <- Epithelial_T_obs$Neddylation_Group
Idents(sce_Epithelial_T) <- sce_Epithelial_T$Neddylation_Group

##跑差异基因分析
diff <- FindMarkers(sce_Epithelial_T, ident.1 = "High", ident.2 = "Low",
                    group.by = "Neddylation_Group", logfc.threshold = 0.25,min.pct = 0.25)
save(diff,file="neddgroup_diffgene.rda")
#筛选logFC绝对值大于等于1的
diff1 <- diff[which(diff$avg_log2FC>=1 | diff$avg_log2FC<=-1),]
diff <- diff1
diff$gene <- rownames(diff)
diff$group <- ""
diff$group <- ifelse(diff$avg_log2FC>=1,"up",'down')

#富集分析,我们这里就以KEGG为例子
group <- data.frame(gene=diff$gene,group=diff$group)#分组情况
#gene转化为ID
Gene_ID <- bitr(diff$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
diff_KEGG <- compareCluster(ENTREZID~group,
                            data=data,
                            fun = "enrichKEGG",#函数选择什么定义什么分析 #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"    
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01,
                            qvalueCutoff = 0.01,
                            organism= "hsa")#物种

diff_GO <- compareCluster(ENTREZID~group,
                          data=data,
                          fun = "enrichGO",#函数选择什么定义什么分析 #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"    
                          #pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          OrgDb='org.Hs.eg.db')

#将gene ID转化为gene symbol
diff_KEGG = setReadable(diff_KEGG,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
diff_GO = setReadable(diff_GO,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
#获取富集分析表格文件
diff_enrich <- diff_KEGG@compareClusterResult
diff_enrich <- diff_GO@compareClusterResult

#diff_KEGG2 <- diff_KEGG %>% group_by(group) %>% top_n(n = 5, wt = -qvalue) 
diff_enrich <- diff_enrich[which(diff_enrich$Description %in% c("Lysosome",
                                                          "Proteasome",
                                                          "Cell cycle",
                                                          "Ubiquitin mediated proteolysis",
                                                          "p53 signaling pathway",
                                                          "RNA degradation",
                                                          "Colorectal cancer")),]
diff_enrich <- diff_enrich[which(diff_enrich$Description %in% c("chaperone binding",
                                                                "ubiquitin-like protein ligase binding",
                                                                "ubiquitin protein ligase binding",
                                                                "ubiquitin binding",
                                                                "ubiquitin-like protein conjugating enzyme activity",
                                                                "ubiquitin-like protein binding")),]
#排序
diff_enrich$group <- factor(diff_enrich$group, levels = c("up","down"))
# 使用排序索引重新排列数据框
diff_enrich <- diff_enrich[order(diff_enrich$group), ]
#terms因子顺序
diff_enrich$Description <- factor(diff_enrich$Description, levels = diff_enrich$Description)

#展示的基因，我们选择每个terms展示5个基因，实际情况可以展示自己关注的基因
diff_enrich$geneID  <- sapply(strsplit(diff_enrich$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))


ggplot(diff_enrich, aes(x = -log10(qvalue), y = rev(Description), fill = group))+
  geom_bar(stat = "identity", width = 0.5)+
  geom_text(aes(x=0.1,y=rev(Description),label = Description),size=4, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(colour = 'black', size = 12),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none")+
  scale_x_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("#CF9999","#9DC3D9"))+  #"#CB5640","#65B0C6" 
  geom_text(data = diff_enrich,
            aes(x = 0.1, y = rev(Description), label = geneID, color = group),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3)+
  scale_color_manual(values = c("#CB5640","#65B0C6"))+
  scale_y_discrete(expand = c(0.1,0))+
  labs(title = "Enrichment of genes",
       y=c("Down                                     Up"))


####GSEA分析
KS_GSEA <- function(gene,
                    LogFC,
                    analysis=c('GO',"KEGG"),
                    package=c('clusterProfiler','fgsea'),
                    OrgDb=c("org.Hs.eg.db", "org.Mm.eg.db")){
  if(OrgDb=="org.Hs.eg.db"){
    organism = "hsa"
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    organism = "mmu"
    species = "Mus musculus"
  }
  
  if(package=="clusterProfiler"){
    entrezID <- bitr(gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = OrgDb)#genesymbol转化为ID
    
    genelist <- LogFC
    names(genelist) <- gene
    genelist <- genelist[names(genelist) %in% entrezID[,1]]
    names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
    genelist <- sort(genelist,decreasing = T)
    
    #install.packages('R.utils')
    R.utils::setOption( "clusterProfiler.download.method",'auto')
    
    if(analysis=="KEGG"){
      KEGG_gesa <- gseKEGG(geneList = genelist,
                           organism = organism,#不同物种查询：https://www.genome.jp/kegg/catalog/org_list.html
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = FALSE,
                           eps = 0)
      
      KEGG_gesa_result = setReadable(KEGG_gesa, OrgDb = OrgDb, keyType = "ENTREZID")
      
      KEGG_gesa_table <- KEGG_gesa_result@result
      write.csv(KEGG_gesa_table, file = './KEGG_gesa_table.csv')
      
      
      return(KEGG_gesa_result)
      
    }
    else{
      GO_gesa <- gseGO(geneList = genelist,
                       ont = "BP",
                       OrgDb=OrgDb,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
      
      GO_gesa_result = setReadable(GO_gesa,OrgDb = OrgDb, keyType = "ENTREZID")
      
      GO_gesa_table <- GO_gesa_result@result 
      write.csv(GO_gesa_table, file = './GO_gesa_table.csv')
      
      return(GO_gesa_result)
      
    }
    
    
  }
  
  if(package=="fgsea"){
    if(analysis=="KEGG"){
      
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    ## GSEA分析
    GSEA_df <- fgsea(pathways = GSEA_geneset, 
                     stats = ranks,
                     minSize=10,
                     maxSize=500,
                     eps=0.0)
    library(data.table)
    fwrite(GSEA_df, file="./GSEA_df.txt", sep="\t", sep2=c("", " ", ""))
    
    return(GSEA_df)
    
  }
  
}
KS_GSEA_plot <- function(inputType=c('clusterProfiler','fgsea'),
                         analysis=c('GO',"KEGG"),
                         data,
                         term,
                         gene,
                         LogFC,
                         OrgDb
){
  
  if(OrgDb=="org.Hs.eg.db"){
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    species = "Mus musculus"
  }
  
  
  
  if(inputType=='clusterProfiler'){
    #clusterprofile
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    
    # define function
    gsInfo <- function(object, terms) {
      geneList <- object@geneList
      
      gsea_result_table <- object@result
      site=which(gsea_result_table$Description==terms, arr.ind = TRUE)
      genesetid=gsea_result_table$ID[site]
      
      geneSetID <- object@result[genesetid, "ID"]
      
      geneSet <- object@geneSets[[geneSetID]]
      exponent <- object@params[["exponent"]]
      df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
      df$ymin <- 0
      df$ymax <- 0
      pos <- df$position == 1
      h <- diff(range(df$runningScore))/20
      df$ymin[pos] <- -h
      df$ymax[pos] <- h
      df$geneList <- geneList
      
      df$Description <- object@result[geneSetID, "Description"]
      return(df)
    }
    
    gsea_plotData = gsInfo(data, terms = term)
    
    gsea_plotData <- gsea_plotData %>%
      mutate("gene_name" = data@gene2Symbol) %>%
      filter(position == 1) 
    
    colnames(gsea_plotData)[2] <- 'y'
    
    test = gsea_plotData
    
    #NES和adjust Pvalue
    gesa_table <- data@result
    terms = gesa_table[gesa_table$Description==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$qvalue, scientific = T,digits = 2)
    
    
    
  }
  
  if(inputType=='fgsea'){
    if(analysis=="KEGG"){
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
      
      
    }
    
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    
    
    ToPlot  <- sapply(data$pathway, function(Pathway){
      pathway <- GSEA_geneset[[Pathway]]
      stats <- ranks
      rnk <- rank(-stats)
      ord <- order(rnk)
      statsAdj <- stats[ord]
      statsAdj <- sign(statsAdj)*(abs(statsAdj)^1)
      statsAdj <- statsAdj/max(abs(statsAdj))
      pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
      pathway <- sort(pathway)
      gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
      bottoms <- gseaRes$bottoms
      tops <- gseaRes$tops
      n <- length(statsAdj)
      xs <- as.vector(rbind(pathway - 1, pathway))
      ys <- as.vector(rbind(bottoms, tops))
      toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0), pathway = Pathway)
      return(list(toPlot))
    })
    
    #点图数据
    gsea_plotData <- ToPlot[[term]]
    
    test <-  gsea_plotData[1:floor(nrow(gsea_plotData)/length(1)),]
    #NES和adjust Pvalue
    terms = data[data$pathway==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$padj, scientific = T,digits = 2)
    
  }
  
  
  
  test$xend = test$x+1
  test$xend2 = ''
  test$xend3 = ''
  
  for (i in 1:nrow(test)) {
    test$xend2[i] = test$xend[i+1] - test$xend[i]
  }
  
  for (i in 1:ncol(test)){
    test[,i][is.na(test[,i])] <- 5
  } 
  
  test$xend2 <- as.numeric(test$xend2)
  
  for (i in 1:nrow(test)) {
    test$xend3[i] = (test$x[i] + test$xend2[i])
  }
  
  test$xend3 <- as.numeric(test$xend3)
  
  
  if(analysis=="KEGG"){
    title_terms = paste0("KEGG:",term)
  }
  
  if(analysis=="GO"){
    title_terms = paste0("GO:",term)
  }
  
  
  if(labels_NES>0){
    fill_color = "#90191B"
  }
  
  if(labels_NES<0){
    fill_color = "#3F90C9"
  }
  
  
  p=ggplot(test) + 
    geom_rect(aes(xmin = x-1,xmax = xend3, ymin = -0.04 , ymax = 0.04, fill=x), lwd=4)+
    scale_fill_gradientn(colours = colorRampPalette(c("#90191B","white","#3F90C9"))(100))+
    geom_rect(aes(xmin = x,xmax = xend, ymin = -0.04, ymax = 0.04, fill=x), color="black", size=0.5)+
    geom_point(data=gsea_plotData, aes(x = x, y = y),fill=fill_color, shape=21, size=4) + 
    geom_hline(yintercept = 0, linetype=3, lwd = 1) +
    scale_x_continuous(expand = c(0.01,0.01))+
    ylab("Enrichment Score") +
    xlab('')+
    labs(title=title_terms)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 15,colour = 'black'),
          legend.position="none",
          plot.title = element_text(size=15,hjust =0.5, face = 'bold'),
          axis.ticks.x = element_blank())
  
  if(labels_NES>0){
    lable_y = max(gsea_plotData$y)-0.2
    
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/1.8)]
    
    ylim = expand_limits(y=c(-0.2, NA))
  }
  
  if(labels_NES<0){
    
    lable_y = min(gsea_plotData$y)+0.2
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/2)]
    ylim = expand_limits(y=c(NA, 0.2))
    
  }
  
  p+ylim+annotate(geom = 'text', label=paste("NES =",labels_NES, "\n", "FDR =", labels_FDR), 
                  x=lable_x, y=lable_y, size=4)
  
  
}

library(Seurat)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(ggplot2)

load("~/jupyter_home/NGmetadata/Neddgroup_Analysis/neddgroup_diffgene.rda")
df <- diff
df$gene = rownames(diff)
df <- df[-which(is.infinite(df$avg_log2FC)==TRUE),]

#高低nedd分组分析-GSEA分析-----------------------------------------------------------------------
##最后都用了fgsea，因为得出结果比较多
GSEA_CP <- KS_GSEA(gene = df$gene,
                   LogFC = df$avg_log2FC,
                   analysis = "KEGG",
                   package = 'clusterProfiler',
                   OrgDb = 'org.Hs.eg.db')
GSEA_F <- KS_GSEA(gene = df$gene,
                  LogFC = df$avg_log2FC,
                  analysis = "KEGG",
                  package = 'fgsea',
                  OrgDb = 'org.Hs.eg.db')
save(GSEA_CP,GSEA_F,file="neddgroup_GSEA.rda")  

p1=KS_GSEA_plot(inputType = "clusterProfiler",
                analysis = "KEGG",
                data = GSEA_CP,
                term = 'Proteasome',
                gene = df$gene,
                LogFC  = df$avg_log2FC,
                OrgDb = 'org.Hs.eg.db')
p1

p2=KS_GSEA_plot(inputType = "fgsea",
                analysis = "KEGG",
                data = GSEA_F,
                term = 'Cell cycle',
                gene = df$gene,
                LogFC = df$avg_log2FC,
                OrgDb = 'org.Hs.eg.db')
p2







####hdWGCNA####
setwd("hdWGCNA")
devtools::install_github('smorabit/hdWGCNA', ref='dev')
install.packages("harmony")
library(hdWGCNA)
library(WGCNA)
library(Seurat) 
library(tidyverse) 
library(igraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
library(qs)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

sce_Epithelial_T2 <- qread(file = "/home/wling_32/jupyter_home/NGmetadata/cNMF/sce_Epithelial_T2.qs")
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
#DefaultAssay(sce_Epithelial_T2) <- 'RNA'
##把Epithelial_T(3w gene版) 进行suerat的umap，然后添加scanpy的X_umap(必须先seruat的umap后才能添加)
sce_Epithelial_T <- sce_Epithelial_T|>
  Seurat::NormalizeData() |>
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) |>
  ScaleData() |>
  RunPCA()
library("harmony")
sce_Epithelial_T <- RunHarmony(sce_Epithelial_T, group.by.vars = "orig.ident")  

ElbowPlot(sce_Epithelial_T, ndims = 50)
pc.num=1:30
sce_Epithelial_T <- RunUMAP(sce_Epithelial_T, reduction="harmony", dims=pc.num) 
sce_Epithelial_T@reductions$X_umap<-sce_Epithelial_T@reductions$umap
Epithelial_T_umap <- read.csv("~/jupyter_home/NGmetadata/Epithelial_T_umap.csv")
sce_Epithelial_T@reductions$X_umap@cell.embeddings<-as.matrix(Epithelial_T_umap)


##开整hdWGCNA
sce_Epithelial_T3 <- SetupForWGCNA(sce_Epithelial_T,
                         wgcna_name = "Epithelial_T",
                         #Select genes that will be used for co-expression network analysis
                         gene_select = "fraction",
                         fraction = 0.05
)
length(sce_Epithelial_T3@misc$Epithelial_T$wgcna_genes)

Epithelial_T_obs <- read.csv("~/jupyter_home/NGmetadata/Epithelial_T_obs.csv", row.names=1)
sce_Epithelial_T3 <- AddMetaData(sce_Epithelial_T3,metadata = Epithelial_T_obs[,c(25:32)])

sce_Epithelial_T3 <- MetacellsByGroups(
  seurat_obj =  sce_Epithelial_T3,
  group.by = c("cNMF_Cluster","Neddylation_Group"),
  k = 25,#需要聚合的细胞数，一般20-75
  max_shared=20,
  reduction = 'pca',
  assay = "RNA",
  slot = "counts",
  ident.group = 'cNMF_Cluster')
qsave(sce_Epithelial_T3,file= 'sce_Epithelial_T3.qs')

sce_Epithelial_T3 <- NormalizeMetacells(sce_Epithelial_T3)

sce_Epithelial_T3 <- SetDatExpr(
  sce_Epithelial_T3,
  group_name = "cNMF_Cluster4", # the name of the group of interest in the group.by column
  group.by='cNMF_Cluster', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

#Select soft-power threshold选择软阈值
#其他参数参照WGCNA包的pickSoftThreshold函数。
sce_Epithelial_T3 <- TestSoftPowers(sce_Epithelial_T3, networkType = 'signed') # you can also use "unsigned" or "signed hybrid"


# plot the results:
plot <- PlotSoftPowers(sce_Epithelial_T3)

# assemble with patchwork
#可视化软阈值
wrap_plots(plot, ncol=2)

#选择阈值
power_table <- GetPowerTable(sce_Epithelial_T3)
head(power_table)
a=power_table$SFT.R.sq
i = 1
for (b in a){
  print(b)
  if (b<0.8) {
    i=i+1
    print(i)
  }else if(b>0.8){
    break
  }
}
select_soft_power = power_table$Power[i]

#Construct co-expression network
#构建网络,其他参数参照WGCNA包的blockwiseConsensusModules
sce_Epithelial_T3 <- ConstructNetwork(
  sce_Epithelial_T3,
  soft_power = 12,
  tom_name = "cNMF_Cluster4_Test",
  setDatExpr=F)

#Module identification
#可视化tree
PlotDendrogram(sce_Epithelial_T3, main='scRNA hdWGCNA Dendrogram')

#inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(sce_Epithelial_T3)

# Compute harmonized module eigengenes
sce_Epithelial_T3 <- ScaleData(sce_Epithelial_T3)
sce_Epithelial_T3 <- ModuleEigengenes(sce_Epithelial_T3, group.by.vars = 'orig.ident')
qsave(sce_Epithelial_T3,file= 'sce_Epithelial_T3_2.qs')
## Harmonize module eigengenes:
hMEs <- GetMEs(sce_Epithelial_T3)
#hMEs <- GetMEs(sce_Epithelial_T4)
# module eigengenes:
MEs <- GetMEs(sce_Epithelial_T3, harmonized=FALSE)


#Compute module connectivity
sce_Epithelial_T3 <- ModuleConnectivity(sce_Epithelial_T3)

# rename the modules
sce_Epithelial_T3 <- ResetModuleNames(sce_Epithelial_T3,new_name = "M")


# get the module table
modules <- GetModules(sce_Epithelial_T3)

#修改module颜色

# get a table of just the module and it's unique color
mod_color_df <- GetModules(sce_Epithelial_T3) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

#不要grey模块
n_mods <- nrow(mod_color_df) - 1

sce_Epithelial_T4 <- sce_Epithelial_T3
# reset the module colors(我们这里有个module)
newcolor <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e",
              "#4aef7b", "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
              "#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92",
              "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
              "#edd05e", "#6f25e8")
sce_Epithelial_T4 <- ResetModuleColors(sce_Epithelial_T4, newcolor)


#重新getmodule
modules <- GetModules(sce_Epithelial_T4)
write.csv(modules, file = 'modules.csv')


#重新做一下基因聚类图
PlotDendrogram(sce_Epithelial_T4, main='scRNA hdWGCNA Dendrogram')


# get hub genes
hub_df <- GetHubGenes(sce_Epithelial_T4, n_hubs = 200)
write.csv(hub_df, file = 'hub_df.csv')

library(UCell)
sce_Epithelial_T4 <- ModuleExprScore(
  sce_Epithelial_T4,
  n_genes = 25,
  method='UCell'
)

#保存文件
qsave(sce_Epithelial_T4,file= 'sce_Epithelial_T3_3.qs')##更新了module的颜色和Ucell评分


###plot1--可放补图
sce_Epithelial_T4 <- qread("sce_Epithelial_T3_3.qs")
# plot genes ranked by kME for each module
#按照kME排序，可视化每个module的基因
PlotKMEs(sce_Epithelial_T4, ncol=6)


###plot2--太多了没用
#Plot module eigengenes--using hMEs
plot_hMEs <- ModuleFeaturePlot(
  sce_Epithelial_T4,
  reduction = "pca",
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  raster = T
)
plot_hMEs

pdf("5_ModuleFeaturePlot_hMEs.pdf",width=10,height=8)

wrap_plots(plot_hMEs, ncol=3)

dev.off()


###plot3--可以用

## add MEvalues to metadata
sce_Epithelial_T4@meta.data <- cbind(
  sce_Epithelial_T4@meta.data,
  GetMEs(sce_Epithelial_T4, harmonized=TRUE)
)

MEs <- GetMEs(sce_Epithelial_T4, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# plot with Seurat's DotPlot function
mods=c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10",
       "M11","M12","M13","M14","M15","M16","M17","M18","M19","M20",
       "M21","M22","M23","M24","M25","M26","M27","M28","M29","M30",
       "M31","M32","M33","M34","M35","M36")
DotPlot(sce_Epithelial_T4, features=mods, group.by = 'Neddylation_Group',
        dot.scale=4)+
  coord_flip()+
  theme_bw()+
  theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_color_gradientn(colours = c("#264D59","#43978D", "white", "#F9AD6A", "#D46C4E"))

#选择6 7 9 14 21 29 31这几个M群基因
##筛选最终基因
hubgene <- hub_df[which(hub_df$module%in%c("M6","M7","M9","M14","M21","M29","M31")),]
neddgene <- read.gmt("/home/wling_32/jupyter_home/NGmetadata/Neddylation_AUcell/REACTOME_NEDDYLATION.v2023.1.Hs.gmt")  
gene <- intersect(hubgene$gene_name,neddgene$gene)
gene
save(hubgene,neddgene,gene,file="Selected_Genes.rda")

###plot4--可以用

p1 <- VlnPlot(sce_Epithelial_T4,
              features = c('M9'),
              group.by = 'Neddylation_Group',
              pt.size = 0)
p1= p1+geom_boxplot(width=.25, fill='white')

## Change axis labels and remove legend:
p1 <- p1 + xlab('') + ylab('hME') + NoLegend()
## Plot output
p1

###plot5--可以用

ModuleNetworkPlot(sce_Epithelial_T4,outdir = "./ModuleNetworks_final")



####WGCNA-普通转录组####
setwd("WGCNA")
library(qs)
library(WGCNA)
library(ggplot2)
library(ggrepel)
enableWGCNAThreads(8)#R语言运行比较慢，可设置线程（量力而行），示例而已

#==============================================================================#
#                 1、读入表达数据和表型数据，例如这里Bulk RNA FPKM             #
#==============================================================================#
#读入数据并检查
sce_Epithelial_T3 <- qread("/home/wling_32/jupyter_home/NGmetadata/hdWGCNA/sce_Epithelial_T3.qs")
wgcnagene <- sce_Epithelial_T3@misc$Epithelial_T$wgcna_genes
rm(sce_Epithelial_T3)
load("~/jupyter_home/NGmetadata/转录组作图/特定基因云雨图/all_clin_1829.rda")
load("~/jupyter_home/NGmetadata/转录组作图/特定基因云雨图/ComBat_data(all).rda")
comgene <- intersect(wgcnagene,rownames(combat_edata))
combat_edata <- combat_edata[comgene,]
human_data <- t(combat_edata)#分析数据需要转置，行为样本，列为基因/蛋白
human_data <- data.frame(human_data)
  
#如果有整理好的表型数据或者临床数据直接读入即可，和表达数据的样本分组对应即可
#这里我们的数据比较简单，所以直接构建一个即可
datTraits <- all_clin_1829
# write.csv(datTraits, file = 'datTraits.csv')


#检查下数据，检查数据中的缺失、离群
gsg = goodSamplesGenes(human_data, verbose = 3)
gsg$allOK
gsg$goodSamples

#本示例数据数据ok，如果数据有不好的，可通过以下代码过滤
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(human_data)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(human_data)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  human_data = human_data[gsg$goodSamples, gsg$goodGenes]
}


#做一下聚类图看看样本情况，是否有离群值
#binarizeCategoricalColumns函数构建一个trait数据，这个函数挺有用，这里是为了聚类指示样本
datTraits1 <- datTraits[,"Dataset"] 
datTraits1 <- as.data.frame(datTraits1)
rownames(datTraits1) <- datTraits$id

pheno <- binarizeCategoricalColumns(datTraits1, dropFirstLevelVsAll = F,minCount=0)
colnames(pheno)
colnames(pheno) <- c('GSE14333','GSE17538','GSE38832','GSE39582','TCGA')


#样本聚类
sampleTree = hclust(dist(human_data), method = "average")
traitColors = numbers2colors(pheno,signed=T)

# 做一个聚类图，下面显示样本
pdf("1.traits_samples.pdf", width = 8, height = 4)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(pheno),
                    main = "Sample dendrogram and trait heatmap",
                    addTextGuide=T,cex.colorLabels=1,
                    cex.dendroLabels = 1)
dev.off()


#==============================================================================#
#                               2、确定软阈值                                  #
#==============================================================================#
powers = c(c(1:10), seq(from = 11, to=20, by=1))
sft = pickSoftThreshold(human_data,#转置的表达矩阵
                        powerVector = powers,
                        verbose = 5)
#作图
pdf("2.softpower1.pdf", width = 8, height = 5)
par(mfrow = c(1,2))

plot(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
abline(h=0.9, col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

#--
#当然了，这个作图你需要修饰的话，也可以提取数据用ggplot完成
data_sft = data.frame(x=sft$fitIndices[,1],
                      y=-sign(sft$fitIndices[,3])*sft$fitIndices[,2])

p1 <- ggplot(data_sft, aes(x,y))+
  geom_point(color='#2f5688', size=5)+
  labs(x='Soft Threshold (power)',y='Scale Free Topology Model Fit,signed R^2',
       title = 'Scale independence')+
  theme_bw()+
  theme(axis.text = element_text(colour = 'black', size = 10),
        legend.position = "none")+
  geom_text(aes(label=x, color='red'), size=3)+
  geom_hline(yintercept = 0.9, linewidth=0.5)

data_mean <- data.frame(x=sft$fitIndices[,1], y=sft$fitIndices[,5])
p2 <- ggplot(data_mean, aes(x,y))+
  geom_point(color='#2f5688', size=5)+
  labs(x='Soft Threshold (power)',y='Mean Connectivity',
       title = 'Mean connectivity')+
  theme_bw()+
  theme(axis.text = element_text(colour = 'black', size = 10),
        legend.position = "none")+
  geom_text(aes(label=x, color='red'), size=3)

p1+p2

#sft$powerEstimate返回的是最佳的阈值，但是具体采用哪个，需要实际运行探索

#软阈值的数值确定一般是选择Scale independence靠近0.9的数字。或者
#Mean connectivity中最靠近0的数字。可以看到9可能就是选择的合适阈值，
#其实用ggplot用点表示出来之后，看阈值也更加清晰。

#==============================================================================#
#                               3、网络构建（1）                               #
#==============================================================================#
#power选择12,构建邻接矩阵
softPower = 12
adjacency = adjacency(human_data, power = softPower,type="signed")
#“signed”表示将  不连接强负相关的基因表达谱，原文中使用signed
#邻接矩阵转为拓扑重叠
TOM = TOMsimilarity(adjacency, TOMType = "signed")
save(TOM,file="TOM.rda")
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


#使用动态混合树切割算法来切割层次聚类树
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;#按照自己的需求设定大小，最小模块大小
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize, 
                            method = 'tree')

dynamicColors = labels2colors(dynamicMods)#为每个module赋予颜色
# table(dynamicColors)#查看基因分成了多少module、每个module的基因

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#合并模块（对相关性比较高的模块）
#计算特征基因eigengenes
MEList = moduleEigengenes(human_data, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
#对模块eigengenes进行聚类
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.3#cut为0.3，对应相关>=0.7及以上的module合并,自行设置
abline(h=MEDissThres, col = "red")


# Call an automatic merging function
merge = mergeCloseModules(human_data, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
table(mergedColors)
#作图--看一下合并前后的对比
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#接下来可视化合并的模块，并将每个样本也添加在上面，看看效果
#单独看merge module
pdf("6.module_merge.pdf", width = 6, height = 5)
plotDendroAndColors(geneTree, mergedColors,c("colors"),
                    dendroLabels = FALSE, hang = 0.01,
                    addGuide = TRUE, guideHang = 0.01)
dev.off()


#==============================================================================#
#                          5、Module与样本性状的相关性                         #
#==============================================================================#
#这里就是很常规的了，大多数文献常见的内容，看一下module与性状的相关性
#使用merge后的每个module的eigengene与表型数据进行相关分析
ME_merge <- merge$newMEs
save(ME_merge,file = "ME_merge.rda")

##用老版本方法
robustY=T
corType="pearsoon"
all_clin_1829$Neddylation_group <- ifelse(all_clin_1829$Neddylation_group=="Low",0,1)
traitData <- all_clin_1829[,c(4,5,20,26)]
if (corType=="pearsoon") {
  moduleTraitCor = cor(ME_merge, traitData, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(ME_merge))
} else {
  moduleTraitCorP = bicorAndPvalue(ME_merge, traitData, robustY=robustY)
  moduleTraitCor = moduleTraitCorP$bicor
  moduleTraitPvalue   = moduleTraitCorP$p
}

#将模块相关性评分大于0.5的提取出来
textCor <- signif(moduleTraitCor,2)#保留小数点后两位
textCor[textCor < 0.5] <- "" 
textMatrix = paste(textCor, "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(ME_merge), 
               cex.lab = 1.2, 
               ySymbols = colnames(ME_merge), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

###做气泡图
library(tidyr)
library(reshape2) 
save(moduleTraitCor,moduleTraitPvalue,file = "moduleTraitCorP.rda")
moduleTraitCor <- as.data.frame(moduleTraitCor)
moduleTraitCor <- rownames_to_column(moduleTraitCor,"Module")
moduleTraitCor$rowname <- rep("Neddylation",8)
moduleTraitCor$name <- moduleTraitCor$Module
moduleTraitCor$value <- moduleTraitCor$Neddylation
cor1 <- moduleTraitCor[,c(6,7,8)]
moduleTraitCor$rowname <- rep("Cell_cycle",8)
moduleTraitCor$name <- moduleTraitCor$Module
moduleTraitCor$value <- moduleTraitCor$`cell cycle`
cor2 <- moduleTraitCor[,c(6,7,8)]
cor <- rbind(cor1,cor2)


moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
moduleTraitPvalue <- rownames_to_column(moduleTraitPvalue,"Module")
moduleTraitPvalue$rowname <- rep("Neddylation",8)
moduleTraitPvalue$name <- moduleTraitPvalue$Module
moduleTraitPvalue$value <- moduleTraitPvalue$Neddylation
p1 <- moduleTraitPvalue[,c(6,7,8)]
moduleTraitPvalue$rowname <- rep("Cell_cycle",8)
moduleTraitPvalue$name <- moduleTraitPvalue$Module
moduleTraitPvalue$value <- moduleTraitPvalue$`cell cycle`
p2 <- moduleTraitPvalue[,c(6,7,8)]
p <- rbind(p1,p2)
p <- p%>%mutate(p_value=case_when(
  value > 1e-3 ~ "*",
  value >1e-10 & value <= 1e-3 ~ "**",
  value > 1e-100 & value <= 1e-10 ~ "***",
  value <= 1e-100 ~ "****"))

library(paletteer)
p$p_value <- factor(p$p_value,levels=c("*","**","***","****"))

p$cor <- cor$value
ggplot()+
  theme(legend.key = element_rect(colour="black"),
        axis.text.x = element_text(angle = 45,hjust=0.5,vjust=0.5,color = "red"))+
  coord_equal()+
  geom_point(data=p,
             aes(x=rowname,y=name,
                 size=p_value,
                 color=cor))+
  scale_color_gradientn(colours = c("#1B8713", "#E7C258", "#E53B59"))

#用这个
ggplot()+
  geom_point(data=p,
             aes(x=rowname,y=name,
                 size=p_value,
                 color=cor))+
  coord_flip()+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_color_gradientn(colours = c("#264D59","#43978D", "white", "#F9AD6A", "#D46C4E"))


#==============================================================================#
#                          6、提取module基因                  
#==============================================================================#
# A = as.data.frame(moduleColors)
# A = A[geneTree$order,]
# 提取module基因
moduleColors = mergedColors
colorOrder = c("grey", unique(c(standardColors(50), unique(mergedColors))))
moduleLabels = match(moduleColors, colorOrder)-1

# Modules
module_dataframe <- data.frame(gene_id=colnames(human_data), 
                               module_name=paste0('module_', moduleLabels), 
                               module_color=moduleColors)

write.csv(module_dataframe, file = "module_dataframe_gene.csv")


#==============================================================================#
#                          8. 感兴趣性状的模块的具体基因分析                  
#==============================================================================#
####8.感兴趣性状的模块的具体基因分析
### 计算模块与基因的相关性矩阵
dataExpr <- human_data
nSamples <- nrow(ME_merge)
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, ME_merge, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, ME_merge, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "cyan"
pheno = "cell cycle"
modNames = substring(colnames(ME_merge), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.25,v=0.7,col="red",lwd=1.5)
##提取该模块的所有基因
module = "cyan"
column = match(module, modNames)
moduleGenes = moduleColors==module
table(moduleGenes)
cyan_module<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][moduleGenes]) 
names(cyan_module)="genename"
save(cyan_module,file = "cyan_module.rda")

##根据条件筛选hub基因
library(tidyverse)
a1 <- rownames_to_column(geneModuleMembership,"gene")
a2 <- rownames_to_column(geneTraitCor,"gene")
geneModuleMembership <- as.data.frame(geneModuleMembership)
hub<- abs(geneModuleMembership$MEcyan)>0.5 & abs(geneTraitCor$`cell cycle`)>0.25
table(hub)
hub
hubgene_cyan <- dimnames(data.frame(dataExpr))[[2]][hub]
hubgene_cyan#选择了0.5和0.25两个阈值
save(hubgene_cyan,geneModuleMembership,geneTraitCor,dataExpr,file = "hubgene_cyangene.rda")


####如何筛选最终基因
load("~/jupyter_home/NGmetadata/hdWGCNA/Selected_Genes.rda")
hubgene_cyan
#library(dplyr)
#hubgene2 <- hubgene %>%group_by(module) %>%slice(1:150)  #选每个模块的前10个
gene3 <- intersect(hubgene_cyan,gene)
gene3

gene <- intersect(hubgene$gene_name,neddgene$gene)
gene <- intersect(hubgene$gene_name,hubgene_cyan)
gene <- intersect(neddgene$gene,hubgene_cyan)   
gene <- intersect(gene,hubgene$gene_name)
##最终选择交集15个基因
#"COMMD4" "GPS1"   "PSMA5"  "PSMB1"  "PSMB2"  "PSMB3"  "PSMB5"  "PSMC3"  "PSMC4"  "PSMC5"  "PSMD13"
#"PSMD2"  "PSMD3"  "PSMD8"  "UBE2M" 

write.csv(hubgene,file="hdWGCNA_hubgene.csv")
write.csv(hubgene_cyan,file="WGCNA_hubgene_cyan.csv")
write.csv(neddgene,file="neddgene.csv")




####空转（GSE226997）####
setwd("/home/wling_32/Spatial/GSE226997_RAW")
install.packages("tidydr")
library(Seurat)
library(ggplot2)
library(hdf5r)
library(tidydr)
P1 <- Load10X_Spatial(data.dir = './GSM7089855_Ajou_Visium_P1/',
                                             filename = "filtered_feature_bc_matrix.h5",
                                             assay = "Spatial")
P2 <- Load10X_Spatial(data.dir = './GSM7089856_Ajou_Visium_P2/',
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial")
P3 <- Load10X_Spatial(data.dir = './GSM7089857_Ajou_Visium_P3/',
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial")
P4 <- Load10X_Spatial(data.dir = './GSM7089858_Ajou_Visium_P4/',
                      filename = "filtered_feature_bc_matrix.h5",
                      assay = "Spatial")
P1$orig.ident <- "P1"
P1@project.name <- "P1"
P2$orig.ident <- "P2"
P2@project.name <- "P2"
P3$orig.ident <- "P3"
P3@project.name <- "P3"
P4$orig.ident <- "P4"
P4@project.name <- "P4"

#merge到一起, 过滤下spots，当然也可以单个单个过滤
Spatial_list <- list(P1,P2,P3,P4)
names(Spatial_list) <- c("P1","P2","P3","P4")
Spatial_merge <-  Reduce(function(x,y) merge(x,y) , Spatial_list) 
Spatial_merge <- Spatial_merge[,Spatial_merge$nCount_Spatial >=5]
Spatial_merge <- Spatial_merge[,Spatial_merge$nFeature_Spatial >=10]

#这里使用SCTransform标准化，使用seurat中的FindIntegrationAnchors进行整合
DefaultAssay(Spatial_merge) <- 'Spatial'
object_splitlist <- SplitObject(Spatial_merge, split.by = "orig.ident")

for (i in names(object_splitlist)) {
  object_splitlist[[i]] <- SCTransform(object_splitlist[[i]], verbose = T, assay = 'Spatial')
}

#这就和做scRNA的没啥区别了，nfeatures可以自行选择调整
Integration.features <- SelectIntegrationFeatures(object.list = object_splitlist, nfeatures = 2000)
object_splitlist <- PrepSCTIntegration(object.list = object_splitlist, anchor.features = Integration.features, verbose = T)


#integration，因为是cca整合，速度可能会稍微慢一点，耐心等待
integration.anchors <- FindIntegrationAnchors(object.list = object_splitlist, normalization.method = "SCT",
                                              anchor.features = Integration.features, verbose = T)
Spatial_integrated <- IntegrateData(anchorset =integration.anchors, normalization.method = "SCT")

#接下来就是降维聚类
# DefaultAssay(Spatial_integrated)
Spatial_integrated <- RunPCA(object = Spatial_integrated, verbose = T)
Spatial_integrated <- FindNeighbors(Spatial_integrated, dims = 1:30)
Spatial_integrated <- FindClusters(Spatial_integrated, resolution = 0.8)#resolution可设置多个，自行选择
Spatial_integrated <- FindClusters(Spatial_integrated, resolution = 0.4)
Spatial_integrated <- RunUMAP(Spatial_integrated, dims = 1:30)
#Spatial_integrated <- RunUMAP(Spatial_integrated, dims = 1:30, umap.method = 'umap-learn', metric = 'correlation')
save(Spatial_integrated, file = "Spatial_integrated.RData")

#因为整合，导致image槽出现多个，只保留对应样本的4个即可（注意选择的规律2.2，3.3这样的）
Spatial_integrated@images[-which(names(Spatial_integrated@images) %in% c('slice1', 'slice1.P2.2', 'slice1.P3.3', 'slice1.P4.4'))] <- NULL

names(Spatial_integrated@images) <- c("P1","P2","P3","P4")


#基本可视化
cols= c("#EDB931","#eb6841","#cc2a36","#00a0b0","#7A989A", "#849271", "#CF9546", "#C67052", "#C1AE8D",
        "#3F6F76", "#C65840", "#62496F", "#69B7CE","#91323A", "#3A4960", "#6D7345", "#D7C969",
        "#C1395E", "#AEC17B", "#E07B42", "#89A7C2", "#F0CA50","#a53e1f", "#457277", "#8f657d", "#8dcee2",
        "#E69253", "#EDB931", "#E4502E", "#4378A0", "#272A2A","#3F6148", "#A4804C", "#4B5F80", "#DBD3A4")

#转录组UMAP降维结果
#Figure1-2
DimPlot(Spatial_integrated, label = F, cols = cols,pt.size = 0.1)+
  theme_dr()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())


#展示空间位置
#Figure3
#SpatialDimPlot(Spatial_integrated, pt.size.factor = 0,ncol = 4) #只展示组织图片
#展示空间图片与cluster spots
SpatialDimPlot(Spatial_integrated, stroke=0.1,ncol=4)& 
  scale_fill_manual(values = cols) &
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

#因为我们是合并数据，如果要展示单个image得图形，只需要选择images参数是对应得slice即可
SpatialDimPlot(Spatial_integrated, images = 'P1', 
                pt.size.factor=1.2, stroke=.1)& 
  scale_fill_manual(values = cols) &
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

#展示基因expression
#figure4
DefaultAssay(Spatial_integrated) <- 'SCT'
SpatialFeaturePlot(Spatial_integrated, 'UBE2M', ncol = 4)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

##提取markers
devtools::install_github("mahmoudibrahim/genesorteR")
library(genesorteR)
remotes::install_github(repo = 'genecell/COSGR')
library(COSG)
library(qs)
#恶性上皮细胞（用了COSG筛选marker）
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
Epithelial_T_obs <- read.csv("~/jupyter_home/NGmetadata/Epithelial_T_obs.csv", row.names=1)
sce_Epithelial_T <- AddMetaData(sce_Epithelial_T,metadata = Epithelial_T_obs[,c(25:32)])
Idents(sce_Epithelial_T) <- sce_Epithelial_T$cNMF_Cluster
#gs <- sortGenes(sce_Epithelial_T@assays$RNA@data, Idents(sce_Epithelial_T))
#gs_marker <- getMarkers(gs, quant = 0.99)
marker_cosg500 <- cosg(
  sce_Epithelial_T,
  groups='all', #考虑全部分组
  assay='RNA',
  slot='data',
  mu=1,         #惩罚项参数，值越大
  remove_lowly_expressed=TRUE,   #是否过滤低表达基因
  expressed_pct=0.1,             #设置低表达的阈值
  n_genes_user=500      #每个cluster定义Top-N个marker gene
)
cNMFcluster <- marker_cosg500$names
cNMFcluster4 <- cNMFcluster$cNMF_Cluster4
##全部细胞
sce_allcell<-qread('/home/wling_32/NGmetadata/sce_allcell_harmony.qs')
Idents(sce_allcell) <- sce_allcell$cell.type
marker_cosg500all <- cosg(
  sce_allcell,
  groups='all', #考虑全部分组
  assay='RNA',
  slot='data',
  mu=1,         #惩罚项参数，值越大
  remove_lowly_expressed=TRUE,   #是否过滤低表达基因
  expressed_pct=0.1,             #设置低表达的阈值
  n_genes_user=500      #每个cluster定义Top-N个marker gene
)
allmarker <- marker_cosg500all$names
epi_marker <- allmarker$Epithelial



###AddmoduleScore评分法(用的这个)
#注意，基因集是list
neddylation_gene <- read.csv("~/jupyter_home/NGmetadata/Neddylation_AUcell/neddylation_gene.csv", sep="")
gene <- intersect(rownames(Spatial_integrated@assays$RNA),neddylation_gene$x)
neddgene <- list(gene)
cNMFcluster4 <- list(cNMFcluster4)
epi_marker <- list(epi_marker)
#计算评分
Spatial_integrated <- AddModuleScore(object = Spatial_integrated, 
                                     features = neddgene, 
                                     name = 'Neddylation')
Spatial_integrated <- AddModuleScore(object = Spatial_integrated, 
                                     features = cNMFcluster4, 
                                     name = 'cNMFcluster4')
Spatial_integrated <- AddModuleScore(object = Spatial_integrated, 
                                     features = epi_marker, 
                                     name = 'Epithelial')

###AUcell评分法（没用，亮度不高）
library(AUCell)
library(clusterProfiler)
cells_rankings <- AUCell_buildRankings(Spatial_integrated@assays$SCT@data)
qsave(cells_rankings,file= 'cells_rankings_Spatial_integrated.qs')

markers <- list()
markers$neddylation <- neddylation_gene$x
markers$cNMFcluster4 <- cNMFcluster$cNMF_Cluster4
markers$Epithelial <- allmarker$Epithelial
cells_AUC <- AUCell_calcAUC(markers, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)

geneSet <- "neddylation"
geneSet <- "cNMFcluster4"
aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
Spatial_integrated$Neddylation_auc  <- aucs
Spatial_integrated$cNMFcluster4_auc  <- aucs


#plot==figure9
SpatialFeaturePlot(Spatial_integrated, 
                   features = c("Epithelial1"), 
                   alpha = c(0.1,1),
                   ncol = 4,
                   min.cutoff = 0,
                   max.cutoff = 1)&
  theme_bw()&
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())


#小提琴图可视化各个区域不同样本评分变化-figure10
VlnPlot(Spatial_integrated, 'Neddylation1', 
        group.by = 'orig.ident', split.by = 'integrated_snn_res.0.4') + facet_grid(~split)+
  scale_fill_manual(values = cols)

qsave(Spatial_integrated,file= 'Spatial_integrated2.qs')


####细胞周期####
setwd("细胞周期相关")
##seurat包分析--用了这个
library(qs)
library(Seurat)
sce_Epithelial_T <- qread(file = "/home/wling_32/NGmetadata/sce_Epithelial_T.qs")
sce_Epithelial_T = CellCycleScoring(sce_Epithelial_T,
                        s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes)
table(sce_Epithelial_T$Phase)

Epithelial_T_obs <- read.csv("~/jupyter_home/NGmetadata/Epithelial_T_obs.csv", row.names=1)
sce_Epithelial_T <- AddMetaData(sce_Epithelial_T,metadata = Epithelial_T_obs[,c(25:32)])
table(sce_Epithelial_T$Phase,sce_Epithelial_T$cNMF_Cluster)
table(sce_Epithelial_T$Phase,sce_Epithelial_T$Neddylation_Group)

Epithelial_T_obs$S.Score <- sce_Epithelial_T$S.Score 
Epithelial_T_obs$G2M.Score <- sce_Epithelial_T$G2M.Score 
Epithelial_T_obs$Phase <- sce_Epithelial_T$Phase 
write.csv(Epithelial_T_obs,file="/home/wling_32/jupyter_home/NGmetadata/Epithelial_T_obs.csv")

##绘图
#线性关系
library(ggstatsplot)
library(ggplot2)
ggscatterstats(Epithelial_T_obs, 
               y ='S.Score', 
               x ='Neddylation_aucell',
               title = "Correlation betwen Neddylation and G2M phase",
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#866AA3",
               yfill = "#D95B5B", 
               marginal.type = "densigram",
               bf.message = F,results.subtitle = F,
               ggtheme=ggplot2::theme_bw())

ggscatterstats(Epithelial_T_obs, 
               y ='cNMF_Cluster4', 
               x ='Neddylation_aucell',
               title = "Correlation betwen Neddylation and cNMF_Cluster4",
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#866AA3",
               yfill = "#F49600", 
               marginal.type = "densigram",
               smooth.line.args = list(size = 1.5, color = "blue", method = "lm", formula = y ~ x,
                                       na.rm = TRUE),
               bf.message = F,
               ggtheme=ggplot2::theme_bw())

##比例图
a <- ggbarstats(
  data             = Epithelial_T_obs,
  x                = Neddylation_Group,
  y                = Phase,
  title            = "Proportion of Cell Cycle Phase in Neddylation group",
  xlab             = "Neddylation group",
  ylab             = "Relative percent (%)",
  proportion.test  = TRUE,
  bf.message  = F,results.subtitle = F,
) + scale_fill_manual(values = c("#7BBEEB","#D69494"))
a


##scran包分析--没seurat好
library(scran)
hs.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
names(hs.pairs)

mt_count = sce_Epithelial_T@assays$RNA@counts
## ID转换
library(org.Hs.eg.db)
library(tidyverse)
gene_ids = AnnotationDbi::select(org.Hs.eg.db, 
                                 keys=rownames(mt_count), 
                                 columns=c("ENSEMBL"), 
                                 keytype="SYMBOL")
gene_ids = gene_ids %>% 
  na.omit() %>% 
  dplyr::distinct(SYMBOL, .keep_all = T) %>% 
  dplyr::distinct(ENSEMBL, .keep_all = T) 
mt_count = mt_count[gene_ids$SYMBOL,]
rownames(mt_count) = gene_ids$ENSEMBL
dim(mt_count)
##周期分析
assignments = scran::cyclone(mt_count,
                             pairs = hs.pairs)
names(assignments)
dim(assignments$normalized.scores)
assignments$normalized.scores[1:3,]
table(assignments$phases)
scores = as.data.frame(assignments$scores)
phases = as.data.frame(assignments$phases)

Epithelial_T_obs <- read.csv("~/jupyter_home/NGmetadata/Epithelial_T_obs.csv", row.names=1)
sce_Epithelial_T <- AddMetaData(sce_Epithelial_T,metadata = Epithelial_T_obs[,c(25:32)])
sce_Epithelial_T <- AddMetaData(sce_Epithelial_T,metadata = scores)
sce_Epithelial_T$phases <- assignments$phases

table(sce_Epithelial_T$phases,sce_Epithelial_T$cNMF_Cluster)
table(sce_Epithelial_T$phases,sce_Epithelial_T$Neddylation_Group)


####转录组作图####
##这里是转录组neddylation的特定基因在nedd分组的表达-云雨图##
install.packages("ggrain")
library("ggrain")
load("~/jupyter_home/NGmetadata/转录组作图/特定基因云雨图/all_clin_1829.rda")
load("~/jupyter_home/NGmetadata/转录组作图/特定基因云雨图/ComBat_data(all).rda")
gene <- as.data.frame(t(combat_edata[c("NEDD8","CDH1","CDH2",
                                       "NAE1","CUL1","CUL2","CUL3","CUL4","CUL5",
                                       "UBE2D1","UBE2M","CCNB1","CCNB2","CDK1","TP53"),]))
identical(row.names(gene),row.names(all_clin_1829))
all_clin <- cbind(all_clin_1829,gene)
#my_comparisons=list(c("High","Low"))
ggplot(all_clin,aes(Neddylation_group,NEDD8,fill=Neddylation_group))+
  geom_rain(point.args=list(alpha = .2,size=0.2),
            boxplot.args.pos=list(
              width=0.05,position=position_nudge(x=0.13)),
            violin.args.pos=list(
              side="r",alpha=1,
              width=0.7,position=position_nudge(x=0.2)))+
  scale_fill_manual(values=c("#D69494","#7BBEEB"))+
  scale_y_continuous(limits = c(6,12))+
  guides(fill='none',color='none')+
  coord_flip()+ 
  ggsignif::geom_signif(textsize = 8,
    comparisons=list(c("High","Low")),
    map_signif_level=TRUE)+
  #stat_compare_means(comparisons = my_comparisons,method = "t.test")+
  theme_bw()+theme(plot.margin = margin(t = 2.5,  # 顶部边缘距离
                              r = 4,  # 右边边缘距离
                              b = 3,  # 底部边缘距离
                              l = 4,  # 左边边缘距离
                              unit = "cm"),
                   axis.text.x = element_text(size = 15),  
                   axis.text.y = element_text(size = 15),
                   axis.title.x = element_text(size = 14), 
                   axis.title.y = element_text(size = 14),)


####这里是UBE2M和GPS1与CCNB1、CCNB2、CDK1和NEDD8的相关性热图
library(corrplot)
library(RColorBrewer)
library(psych)
load("~/jupyter_home/NGmetadata/转录组作图/特定基因云雨图/ComBat_data(all).rda")
cordata <- as.data.frame(t(combat_edata[c("UBE2M","GPS1","CCNB1","CCNB2","CDK1","NEDD8"),]))
corr <- cor(cordata)
corr_results <- corr.test(cordata, method = "pearson")
corr_results
corrplot(corr, method = 'square', order = 'AOE', type = 'upper', 
         col=colorRampPalette(c("#264D59","#43978D", "white", "#F9AD6A", "#D46C4E"))(56),
         diag = F,
         is.corr=T,
         addCoef.col = "black")


