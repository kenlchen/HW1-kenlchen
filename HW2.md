HW2
========================================================

Assignment: **This is your second homework due Tuesday Feb 18.**

Reproduce the results presented in Figure 2 of the following paper: 
Qian, F., Bolen, C. R., Jing, C., Wang, X., Zheng, W., Zhao, H., et al. (2013). Impaired toll-like receptor 3-mediated immune responses from macrophages of patients chronically infected with hepatitis C virus. Clinical and Vaccine Immunology : CVI, 20(2), 146–155. doi:10.1128/CVI.00530-12

You will have to:

1. Get the data from GEO
2. Normalize the data (if necessary)
3. Use limma to test for differential expression
4. Display the results using a heatmap [Hint: Use the pheatmap package]

Prior to running code, please ensure that the following packages are installed:
GEOmetadb
GEOquery
Biobase
lumi
data.table
pheatmap
RSQLite

If not, run the following chunk of code by deleting "eval=FALSE"
Packages already downloaded can be removed from this operation by deleting their names from the input inside the parentheses


```r
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GEOmetadb", "GEOquery", "lumi"))
install.packages(c("data.table", "pheatmap", "RSQLite"))
```


Additionally, please download and save (to the planned working directory) the GEO expression set list associated with the GEO accession number GSE40812. Load this data set into a variable named "dataSet".

The following code will download and save the required data. The data must be saved as an object name "dataSet" with the file name "GSE40812_getGEO" (or renamed such after) for the rest of the code to work


```r
library(GEOquery)  #load needed library
```

```
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Setting options('download.file.method.GEOquery'='auto')
```

```r
dataSet = getGEO("GSE40812")  #download GEO eset (as a list of esets)
```

```
## Found 1 file(s)
## GSE40812_series_matrix.txt.gz
## File stored at: 
## /var/folders/bf/j3ypc9090f19v32b9qrjsmj00000gn/T//Rtmp1VFNRr/GPL10558.soft
```

```r
# save(dataSet, file = 'GSE40812_getGEO')
```


Subset the expression set to data from macrophages, clean data labels, and retain only those necessary: patient ID, treatment, 

```r
library(Biobase)
# load('GSE40812_getGEO')
eset = dataSet[[1]]
pd = pData(eset)
pdMacrophage = pd[pd$source_name_ch1 == "Monocyte-derived Macrophage", ]
pdMacrophage$cell = gsub(".* ", "", pdMacrophage$source_name_ch1)  #create clean cell source column
pdMacrophage$HCV = gsub(".*: ", "", pdMacrophage$characteristics_ch1)
pdMacrophage$HCV = ifelse(pdMacrophage$HCV == "Neg", "-", "+")  #create clean HCV status column
pdMacrophage$treatment = ifelse(pdMacrophage$characteristics_ch1.2 == "treatment: Mock", 
    "mock", "polyic")  #create clean perturbation/treatment column

# create clean patient ID column by splitting the 'title' column
title = pdMacrophage$title
titleDF = data.frame(strsplit(as.character(title), "_"))
pdMacrophage$patient = t(titleDF)[, 2]

# create cleaned pData data frame retaining only relevant columns
pdMacrophageCleaned = pdMacrophage[, c("patient", "cell", "HCV", "treatment")]

# subset eset and assign cleaned pData to new eset
esetMacrophage = eset[, rownames(pdMacrophageCleaned)]
pData(esetMacrophage) = pdMacrophageCleaned
```


This data has already been normalized. However, normalize using lumi library for good measure.


```r
library(lumi)
```

```
## Warning: replacing previous import by 'graphics::image' when loading
## 'methylumi'
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
```

```
## Warning: replacing previous import by 'nleqslv::nleqslv' when loading
## 'lumi'
```

```r
esetMacrophage = lumiN(esetMacrophage)
```

```
## Perform quantile normalization ...
```


Use Limma to test for differential expression 
First find the set of genes (probes) that respond to the poly IC treatment (ie. are differentially expressed between mock and poly IC)


```r
library(limma)
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
mmTreatment = model.matrix(~treatment, esetMacrophage)  #design matrix
fitTreatment = lmFit(esetMacrophage, mmTreatment)
ebayTreatment = eBayes(fitTreatment)
ttTreatment = topTable(ebayTreatment, coef = "treatmentpolyic", number = Inf, 
    sort.by = "none")  #calculate statistics on differential expression
ttTreatmentDiff = ttTreatment[ttTreatment$adj.P.Val < 0.05 & abs(ttTreatment$logFC) > 
    log2(1.5), ]  #subset to only include probes below meeting the threshold of FDR of 0.05 and fold change > 1.5
dim(ttTreatmentDiff)  #should have 1146 rows to match 1146 probes mentioned in supplemental data
```

```
## [1] 1146   36
```

```r
esetDiff = esetMacrophage[rownames(ttTreatmentDiff), ]  #subset expression set to only these probes
```


Now, find the probes for which the response to poly IC is different between HCV + and - patients. This chunk of code does not yield correct results. I'm not sure why not.

```r
mmHCVTreat = model.matrix(~treatment + HCV, esetDiff)
fitHCVTreat = lmFit(esetDiff, mmHCVTreat)
ebayHCVTreat = eBayes(fitHCVTreat)
ttHCVTreat = topTable(ebayHCVTreat, coef = c("HCV+", "treatmentpolyic"), number = Inf, 
    sort.by = "none")
ttHCVTreatDiff = ttHCVTreat[ttHCVTreat$P.Value < 0.1, ]
esetFinal = esetDiff[rownames(ttHCVTreatDiff)]
```


Create an expression set for which the rows are probes and the columns are patients. Each entry is the fold change between the macrophage response to the poly IC treatment and the mock treatment. 


```r
library(reshape)
```

```
## Loading required package: plyr
## 
## Attaching package: 'reshape'
## 
## The following objects are masked from 'package:plyr':
## 
##     rename, round_any
```

```r
pdPremelt = pdMacrophageCleaned
rownames(pdPremelt) = rownames(pdMacrophageCleaned)
pdPremelt$treatment = ifelse(pdPremelt$treatment == "mock", -1, 1)  #replace values in treatment column with 1 for polyic and -1 for mock
pdMelted = melt(pdPremelt)  #'melt' the data frame. 'melt' is somewhat like reversing the creation of a pivot table in excel
```

```
## Using patient, cell, HCV as id variables
```

```r
pdMelted$GSM = rownames(pdPremelt)  #add a column for the sample numbers
sampleByPatient = cast(pdMelted, GSM ~ patient ~ variable)  #'cast' is somewhat like creating a pivot table. In this case the rows are sample numbers (GSM), the columns are patients, and the values in the cell denoted by a sample number and patient ID is the value originally associated with the treatment column corresponding to the original row with that sample/patient combo.
sampleByPatient[which(is.na(sampleByPatient))] = 0  #because samples are specific to patients (ie. there were no rows originally corresponding to most of the possible sample id-patient number pairings), most of the entries are NA. Change them to zero
sbpMatrix = matrix(sampleByPatient, nrow = 40, ncol = 20)

# create new expression set with expression matrix of probes v patients
esetDiff2 = esetDiff
foldChange = exprs(esetDiff) %*% sbpMatrix
exprs(esetDiff2) = exprs(esetDiff) %*% sbpMatrix
rownames(exprs(esetDiff2)) = rownames(exprs(esetDiff))

# Create a table of patient ID and HCV status, assign to pData
patientStatus = unique(pdMacrophageCleaned[, c("patient", "HCV")])
patientStatus = patientStatus[with(patientStatus, order(patient)), ]
rownames(patientStatus) = NULL
pData(esetDiff2) = patientStatus
rownames(pData(esetDiff2)) = patientStatus$patient
pData(esetDiff2)$patient = NULL

# change the column names of the expression set expression matrix to reflect
# the patient IDs
colnames(exprs(esetDiff2)) = patientStatus$patient

# create model matrix, linear model, and eBayes
mmHCV = model.matrix(~HCV, esetDiff2)
# mmHCV2 = model.matrix(~HCV, patientStatus) # this should be the same
fitHCV = lmFit(esetDiff2, mmHCV)
# fitHCV2 = lmFit(foldChange, mmHCV2)
ebayHCV = eBayes(fitHCV)
ttHCV = topTable(ebayHCV, coef = "HCV+", number = Inf, sort.by = "none")
ttHCVDiff = ttHCV[ttHCV$P.Value < 0.1, ]
length(ttHCVDiff$ID)  #should equal 43 according to paper
```

```
## [1] 43
```


Draw Heatmap


```r
library(pheatmap)
esetHM = esetDiff[rownames(ttHCVDiff), ]
a = pdMacrophageCleaned[order(pdMacrophageCleaned$patient), ]
b = a[order(a$HCV), ]
c = b[order(b$treatment), ]
esetHM = esetHM[, rownames(c)]
matrixHM = exprs(esetHM)
mu = rowMeans(matrixHM)
sigma = apply(matrixHM, 1, sd)
matrixZ = (matrixHM - mu)/sigma
colnames(matrixZ) = c$patient

hm = pheatmap(matrixZ, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 
