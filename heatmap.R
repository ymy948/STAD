# library(pheatmap)
# 
# data1<- read.csv("E:/sirebrowser/STAD/miR/分析/rpm筛选.csv",head = T,row.names=1)
# 
# data1<- read.table("E:/sirebrowser/3-COXPH/42-CESC_diff_miRNA.txt",head = T,row.names=1,sep="\t")
# range(data1)
# dim(data1)
# data1<-as.matrix(data1)
# 
# bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
# 
# a<-pheatmap(data1 , scale = "none", clustering_method = "ward",
#             clustering_distance_rows = "euclidean",
#             cluster_row=TRUE,cluster_cols=TRUE,
#             labels_col = "", labels_row = "",
#             color = c(colorRampPalette(colors = c("MediumBlue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
#             legend_breaks=seq(-8,8,2),
#             breaks=bk
#            #color = colorRampPalette(colors = c("blue","white","red"))(100)
#            #color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#            )
# 
# dev.off()         
# 
# order_row = a$tree_row$order  #记录热图的行排序
# order_col = a$tree_col$order    #记录热图的列排序
# datat = data.frame(data1[order_row,order_col])   # 按照热图的顺序，重新排原始数据
# datat = data.frame(rownames(datat),datat,check.names =F)  # 将行名加到表格数据中
# colnames(datat)[1] = "Sample" 
# write.csv(datat,file="E:/sirebrowser/STAD/miR/分析/rpm筛选排序.csv",row.names=FALSE,quote = FALSE)  #输出结果，按照热图中的顺序


library(gplots)
data1<- read.csv("E:/sirebrowser/STAD/miR/分析/rpm筛选排序.csv",head = T,row.names=1,stringsAsFactors = F)
data1<-as.matrix(data1)
colors <- colorRampPalette(c("MediumBlue",'white', "red"))(100)

heatmap.2(data1,trace="none",col=colors,cexRow = 0.4, cexCol = 0.1, scale="row")
heatmap.2(data1,  col=redblue(100),  key=TRUE, keysize=1, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


