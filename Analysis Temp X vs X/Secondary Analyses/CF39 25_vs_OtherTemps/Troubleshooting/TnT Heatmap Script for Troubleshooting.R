
library(readr)
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(dendextend)

TnT_List_Analysis_Shrunk2 <- read_csv("../TnT CF39 thermoregulation Shrunk2.csv")
rawNamesOfColumns <- c("locus", "FoldChange_25vs29", "FoldChange_25vs33", "FoldChange_25vs37", "FoldChange_25vs41", "Gene", "Product", "PAO1_ID", "Protein Description")
colnames(TnT_List_Analysis_Shrunk2) <- rawNamesOfColumns

FoldChangeData <- data.frame(TnT_List_Analysis_Shrunk2$locus, TnT_List_Analysis_Shrunk2$FoldChange_25vs29, TnT_List_Analysis_Shrunk2$FoldChange_25vs33, TnT_List_Analysis_Shrunk2$FoldChange_25vs37, TnT_List_Analysis_Shrunk2$FoldChange_25vs41)
NamesOfColumns <- c("locus", "FoldChange_25vs29", "FoldChange_25vs33", "FoldChange_25vs37", "FoldChange_25vs41")
colnames(FoldChangeData) <- NamesOfColumns
rownames(FoldChangeData) <- FoldChangeData$locus
FoldChangeData <- FoldChangeData[,-1]



FoldChangeDataNoNA <- FoldChangeData
FoldChangeDataNoNA[is.na(FoldChangeDataNoNA)] <- 0
FoldChangeDataNoNAMat <- data.matrix(FoldChangeDataNoNA)

mat <- FoldChangeDataNoNAMat
colorSetRedBlue <- colorRamp2(c(-10, 0, 10), c(rgb(0, 114/255, 178/255), "white", rgb(213/255, 94/255, 0))) #First color is Blue and second is Vermillion from the Nature Colorblind palette
ColOrder = order(as.numeric(gsub("FoldChange_", "", colnames(mat))))


row_dend = as.dendrogram(hclust(dist(mat, method = "maximum")), method = "complete")
row_dend = color_branches(row_dend, k = 4) # `color_branches()` returns a dendrogram object

split = data.frame(cutree(hclust(dist(mat, method = "maximum")), h=4))


InCellTextSize = 7.75
RowLabelTextCell = 9
EmptyLableColumn = array("", nrow = 1, ncol = length(TnT_List_Analysis_Shrunk2$Gene))


NameDataFrame = data.frame(TnT_List_Analysis_Shrunk2$locus, TnT_List_Analysis_Shrunk2$Gene)
NameDataArray = array(data = c(TnT_List_Analysis_Shrunk2$locus, TnT_List_Analysis_Shrunk2$Gene), dim = c(10, 2))
EmptyColumn = cbind((1:length(TnT_List_Analysis_Shrunk2$locus)), " ")
EmptyColumn = EmptyColumn[, 2]


HeatmapTest <- Heatmap(
  matrix = mat, 
  cluster_rows = color_branches(as.dendrogram(hclust(dist(mat, method = "canberra")), method = "complete"), k = 5),
  row_split = 5,
  row_gap = unit(1.75, 'mm'),
  border = FALSE, 
  
  show_row_names = FALSE,
  row_labels = TnT_List_Analysis_Shrunk2$PAO1_ID,
  
  
  row_names_max_width = unit(50, "cm"), #Does not appear to be doing anything
  right_annotation = rowAnnotation(
    Col1 = anno_text(TnT_List_Analysis_Shrunk2$locus, location = unit(2, "mm"), gp = gpar(fontsize = RowLabelTextCell)), 
    Col2 = anno_text(TnT_List_Analysis_Shrunk2$Gene, location = unit(7, "mm"), gp = gpar(fontsize = RowLabelTextCell)),
    Col3 = anno_text(TnT_List_Analysis_Shrunk2$PAO1_ID, location = unit(9, "mm"), gp = gpar(fontsize = RowLabelTextCell))#,
  ),
  
  
  row_title = "Transcribed Gene Loci", 
  row_names_gp = gpar(fontsize = RowLabelTextCell),
  
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(mat[i, j] < 0 || mat[i, j] > 0)
      grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = InCellTextSize))},
  
  
  column_dend_reorder = FALSE, 
  column_order = ColOrder, 
  column_names_rot = 45,
  
  
  column_title = "Conditions tested (CF39S vs CF39)", 
  column_title_side = "bottom",
  
  
  row_dend_width = unit(3.1, "cm"),
  
  
  col = colorSetRedBlue,
  height = unit(10, 'cm')
  , raster_quality = 10
)
#HeatmapTest <- draw(HeatmapTest, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


HeatmapTest