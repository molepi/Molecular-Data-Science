# Download expression dataset #

To download the data use the following code from within `R`: 

```r
file <- url("https://raw.githubusercontent.com/molepi/FOS2017/master/QuadraticModelsAge/GSE55164_small_RNA_Seq_normalized.txt")
data <- read.table(file, sep="\t", header=TRUE)

file <- url("https://raw.githubusercontent.com/molepi/FOS2017/master/QuadraticModelsAge/FOS_brain_training.txt")
data <- read.table(file, sep="\t", header=TRUE)
```
