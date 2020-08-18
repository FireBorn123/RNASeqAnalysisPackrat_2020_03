---
title: "ComplexHeatMap Test Script"
author: "Trevor Randall"
date: "12/05/2020"
output: html_document
---

Analysis from: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html




#Clear Global Environment

```r
remove(list = ls())
```



```r
#Set Library Directory
#PC Path
.libPaths(c("L:/RStudios/RPackRatLibLocations", "L:/RStudios/RPackRat_2019_04_DESEQLibs"))

#Set working directory
#PC Path
setwd("E:/Dropbox/Dropbox/Harrison Lab - Trevor Randall/RNASeq Analysis/RNASeqAnlyPkrat_2020_03/Analysis Temp X vs X/Secondary Analyses/CF39 25_vs_OtherTemps/TnT Heatmaps")
#Laptop Path
#Set Library Directory
#.libPaths(c(""))
#setwd("")
#sink(file = "./RSessionRawRun.txt") 
```


#Libraries required

```r
library(readr)
library(ComplexHeatmap)
```

```
## Loading required package: grid
```

```
## ========================================
## ComplexHeatmap version 2.5.3
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite:
## Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##   genomic data. Bioinformatics 2016.
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(ComplexHeatmap))
## ========================================
```

```r
library(circlize)
```

```
## Warning: package 'circlize' was built under R version 3.5.3
```

```
## ========================================
## circlize version 0.4.8
## CRAN page: https://cran.r-project.org/package=circlize
## Github page: https://github.com/jokergoo/circlize
## Documentation: http://jokergoo.github.io/circlize_book/book/
## 
## If you use it in published research, please cite:
## Gu, Z. circlize implements and enhances circular visualization 
##   in R. Bioinformatics 2014.
## ========================================
```

```r
library(cluster)
library(dendextend)
```

```
## 
## ---------------------
## Welcome to dendextend version 1.13.4
## Type citation('dendextend') for how to cite the package.
## 
## Type browseVignettes(package = 'dendextend') for the package vignette.
## The github page is: https://github.com/talgalili/dendextend/
## 
## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
## Or contact: <tal.galili@gmail.com>
## 
## 	To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
## ---------------------
```

```
## 
## Attaching package: 'dendextend'
```

```
## The following object is masked from 'package:stats':
## 
##     cutree
```

```r
#library(seriation)
```

PTA = Pure Temperature Analysis
#Dataset

```r
TnT_List_Analysis_Shrunk2 <- read_csv("../TnT CF39 thermoregulation Shrunk2-2.csv")
```

```
## Parsed with column specification:
## cols(
##   locus = col_character(),
##   `CF39 FC 25vs29` = col_double(),
##   `CF39 FC 25vs33` = col_double(),
##   `CF39 FC 25vs37` = col_double(),
##   `CF39 FC 25vs41` = col_double(),
##   gene = col_character(),
##   product = col_character(),
##   `PAO1 ID (Pseudomonas DB Protein search)` = col_character(),
##   `Protein Description` = col_character()
## )
```

```r
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
```


```r
mat <- FoldChangeDataNoNAMat
colorSetRedBlue <- colorRamp2(c(-10, 0, 10), c(rgb(0, 114/255, 178/255), "white", rgb(213/255, 94/255, 0))) #First color is Blue and second is Vermillion from the Nature Colorblind palette
ColOrder = order(as.numeric(gsub("FoldChange_", "", colnames(mat))))
```

```
## Warning in order(as.numeric(gsub("FoldChange_", "", colnames(mat)))): NAs
## introduced by coercion
```

```r
RowOrder = sort(rownames(mat))
```



#Possible clustering options that i could see

        clustering_distance_rows = "pearson", #Pre-defined distance method (1 - pearson) #Seems to look the most split up along with the single distancing
        clustering_distance_rows = "spearman",
        clustering_distance_rows = "kendall",
        clustering_distance_rows = function(mat) dist(mat), #A function that calculates distance between martix points
        clustering_distance_rows = function(x, y) 1 - cor(x, y), #A function that calculates pairwise distance between martix points
        clustering_method_rows = "single",
        clustering_method_rows = "complete",
        cluster_rows = diana(mat), #Cluster used here
        cluster_rows = function(mat) as.dendrogram(diana(mat)),
        
##TODO Look into Seration package to improve pattern visibility

#Coloration of row dendrigram

```r
row_dend = as.dendrogram(hclust(dist(mat, method = "maximum")), method = "complete")
row_dend = color_branches(row_dend, k = 4) # `color_branches()` returns a dendrogram object

split = data.frame(cutree(hclust(dist(mat, method = "maximum")), h=4))
```

  cluster_rows = function(mat) color_branches(as.dendrogram(diana(mat)), k = 10), #Same as above dendrogram
  cluster_rows = function(mat) color_branches(as.dendrogram(diana(mat)), k = 10), #Looks good
  cluster_rows = color_branches(as.dendrogram(hclust(dist(mat, method = "euclidean")), method = "complete"), k = 10) #Looks good
  cluster_rows = color_branches(as.dendrogram(hclust(dist(mat, method = "canberra")), method = "complete"), k = 10) #Looks very logical
  cluster_rows = color_branches(as.dendrogram(hclust(dist(mat, method = "maximum")), method = "complete"), k = 10) #Looks good










