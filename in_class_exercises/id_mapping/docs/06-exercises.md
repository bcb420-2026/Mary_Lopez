# Exercises

These are intended to be done **after** completing the worked examples.




```r
library(readr)
library(biomaRt) 
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:biomaRt':
## 
##     select
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```


## Exercise 1 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119732

Using **GSE119732**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Count how many contain a `.`.
3. Create a new column with versions stripped.
4. Map the identifiers to HGNC symbols.


```r
# helper function
strip_ensembl_version <- function(x) sub("\\..*$", "", x)

source("./fetch_geo_supp.R")

params_119732 <- list(
  gse       = "GSE119732",
  file_grep = "count_table_RNA_seq"
)

fetch_geo_supp(gse = params_119732$gse)
```

```
## Setting options('download.file.method.GEOquery'='auto')
```

```
## Setting options('GEOquery.inmemory.gpl'=FALSE)
```

```r
path_119732 <- file.path("data", params_119732$gse)
files_119732 <- list.files(
  path       = path_119732,
  pattern    = params_119732$file_grep,
  full.names = TRUE,
  recursive  = TRUE
)

files_119732  # should show GSE119732_count_table_RNA_seq.txt.gz
```

```
## [1] "data/GSE119732/GSE119732_count_table_RNA_seq.txt.gz"
```

```r
# read the count table
x1 <- read_tsv(files_119732[1], show_col_types = FALSE)

# the ID column is the first column
id_col1 <- names(x1)[1]
ids1    <- as.character(x1[[id_col1]])

## 1. Extract the first 20 IDs.
first20_1 <- ids1[1:20]
first20_1
```

```
##  [1] "ENSG00000223972.5" "ENSG00000227232.5" "ENSG00000278267.1"
##  [4] "ENSG00000243485.4" "ENSG00000237613.2" "ENSG00000268020.3"
##  [7] "ENSG00000240361.1" "ENSG00000186092.4" "ENSG00000238009.6"
## [10] "ENSG00000239945.1" "ENSG00000233750.3" "ENSG00000268903.1"
## [13] "ENSG00000269981.1" "ENSG00000239906.1" "ENSG00000241860.6"
## [16] "ENSG00000222623.1" "ENSG00000241599.1" "ENSG00000279928.1"
## [19] "ENSG00000279457.3" "ENSG00000273874.1"
```

```r
## 2. Count how many contain a "."
sum(grepl("\\.", first20_1))
```

```
## [1] 20
```

```r
## 3. Create a new column with versions stripped.
x1$ID_stripped <- strip_ensembl_version(ids1)

head(x1[, c(id_col1, "ID_stripped")], 10)
```

```
## # A tibble: 10 × 2
##    gene_id           ID_stripped    
##    <chr>             <chr>          
##  1 ENSG00000223972.5 ENSG00000223972
##  2 ENSG00000227232.5 ENSG00000227232
##  3 ENSG00000278267.1 ENSG00000278267
##  4 ENSG00000243485.4 ENSG00000243485
##  5 ENSG00000237613.2 ENSG00000237613
##  6 ENSG00000268020.3 ENSG00000268020
##  7 ENSG00000240361.1 ENSG00000240361
##  8 ENSG00000186092.4 ENSG00000186092
##  9 ENSG00000238009.6 ENSG00000238009
## 10 ENSG00000239945.1 ENSG00000239945
```

```r
## 4. Map the identifiers to HGNC symbols.

ensembl_ids1 <- unique(x1$ID_stripped)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map1 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids1,
  mart       = ensembl
)

expr_mapped_1 <- x1 %>%
  mutate(ensembl_gene_id = ID_stripped) %>%
  left_join(
    map1,
    by = "ensembl_gene_id",
    relationship = "many-to-many"  # acknowledge non 1:1 mapping
  )

head(expr_mapped_1[, 1:8])
```

```
## # A tibble: 6 × 8
##   gene_id              A1    A2    A3    A4    B1    B2    B3
##   <chr>             <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
## 1 ENSG00000223972.5     0     0     0     0     0     0     0
## 2 ENSG00000227232.5    79   119    84    50    80    72    50
## 3 ENSG00000278267.1    17    10    22    19    19    22    16
## 4 ENSG00000243485.4     0     0     0     0     0     0     0
## 5 ENSG00000237613.2     0     0     0     0     0     0     0
## 6 ENSG00000268020.3     0     0     0     0     0     0     0
```


## Exercise 2 — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122380

Using **GSE122380**, confirm whether the ID column contains Ensembl IDs with version suffixes.

1. Extract the first 20 IDs.
2. Create a new column with versions stripped.
3. Map the identifiers to HGNC symbols.
4. What is different about this file? 



```r
params_122380 <- list(
  gse       = "GSE122380",
  file_grep = "raw_counts"
)

fetch_geo_supp(gse = params_122380$gse)

path_122380 <- file.path("data", params_122380$gse)

files_122380 <- list.files(
  path       = path_122380,
  pattern    = params_122380$file_grep,
  full.names = TRUE,
  recursive  = TRUE
)

files_122380  # should be "data/GSE122380/GSE122380_raw_counts.txt.gz"
```

```
## [1] "data/GSE122380/GSE122380_raw_counts.txt.gz"
```

```r
# read the raw counts

x2 <- read_table(
  files_122380[1],
  comment = "#",    
  col_types = cols()  
)

# check the structure
dim(x2)
```

```
## [1] 16319   298
```

```r
head(x2[, 1:5])
```

```
## # A tibble: 6 × 5
##   Gene_id         `18489_0` `18489_10` `18489_11` `18489_12`
##   <chr>               <dbl>      <dbl>      <dbl>      <dbl>
## 1 ENSG00000000419       825        549        576        599
## 2 ENSG00000000457       175        238        279        254
## 3 ENSG00000000460       337        103        120        124
## 4 ENSG00000000938         7          1          0          0
## 5 ENSG00000000971         0          3          5         10
## 6 ENSG00000001036      2251        860        819        732
```

```r
names(x2)
```

```
##   [1] "Gene_id"  "18489_0"  "18489_10" "18489_11" "18489_12" "18489_13"
##   [7] "18489_14" "18489_15" "18489_1"  "18489_2"  "18489_3"  "18489_4" 
##  [13] "18489_5"  "18489_6"  "18489_7"  "18489_8"  "18489_9"  "18499_0" 
##  [19] "18499_10" "18499_11" "18499_12" "18499_13" "18499_14" "18499_15"
##  [25] "18499_1"  "18499_2"  "18499_3"  "18499_4"  "18499_5"  "18499_6" 
##  [31] "18499_7"  "18499_8"  "18499_9"  "18505_0"  "18505_10" "18505_11"
##  [37] "18505_12" "18505_13" "18505_14" "18505_15" "18505_1"  "18505_2" 
##  [43] "18505_3"  "18505_4"  "18505_5"  "18505_6"  "18505_7"  "18505_8" 
##  [49] "18505_9"  "18508_0"  "18508_10" "18508_11" "18508_12" "18508_13"
##  [55] "18508_14" "18508_15" "18508_1"  "18508_2"  "18508_3"  "18508_4" 
##  [61] "18508_5"  "18508_6"  "18508_7"  "18508_8"  "18508_9"  "18511_0" 
##  [67] "18511_10" "18511_11" "18511_12" "18511_13" "18511_14" "18511_15"
##  [73] "18511_1"  "18511_2"  "18511_3"  "18511_4"  "18511_5"  "18511_6" 
##  [79] "18511_7"  "18511_8"  "18511_9"  "18517_0"  "18517_10" "18517_11"
##  [85] "18517_12" "18517_13" "18517_14" "18517_15" "18517_1"  "18517_2" 
##  [91] "18517_3"  "18517_4"  "18517_5"  "18517_6"  "18517_7"  "18517_8" 
##  [97] "18517_9"  "18520_0"  "18520_10" "18520_11" "18520_12" "18520_13"
## [103] "18520_14" "18520_15" "18520_1"  "18520_3"  "18520_4"  "18520_5" 
## [109] "18520_6"  "18520_7"  "18520_8"  "18520_9"  "18855_0"  "18855_10"
## [115] "18855_11" "18855_12" "18855_13" "18855_14" "18855_15" "18855_1" 
## [121] "18855_3"  "18855_4"  "18855_5"  "18855_6"  "18855_7"  "18855_8" 
## [127] "18855_9"  "18858_0"  "18858_10" "18858_11" "18858_12" "18858_13"
## [133] "18858_14" "18858_15" "18858_1"  "18858_2"  "18858_3"  "18858_4" 
## [139] "18858_5"  "18858_6"  "18858_7"  "18858_8"  "18858_9"  "18870_0" 
## [145] "18870_10" "18870_11" "18870_12" "18870_13" "18870_14" "18870_15"
## [151] "18870_1"  "18870_2"  "18870_3"  "18870_5"  "18870_6"  "18870_7" 
## [157] "18870_8"  "18870_9"  "18907_0"  "18907_10" "18907_11" "18907_12"
## [163] "18907_13" "18907_14" "18907_15" "18907_1"  "18907_2"  "18907_3" 
## [169] "18907_4"  "18907_5"  "18907_6"  "18907_7"  "18907_8"  "18907_9" 
## [175] "18912_0"  "18912_10" "18912_11" "18912_12" "18912_13" "18912_14"
## [181] "18912_15" "18912_1"  "18912_2"  "18912_3"  "18912_4"  "18912_5" 
## [187] "18912_6"  "18912_7"  "18912_8"  "18912_9"  "19093_0"  "19093_10"
## [193] "19093_11" "19093_12" "19093_13" "19093_14" "19093_15" "19093_1" 
## [199] "19093_2"  "19093_3"  "19093_4"  "19093_5"  "19093_6"  "19093_7" 
## [205] "19093_8"  "19093_9"  "19108_0"  "19108_10" "19108_11" "19108_12"
## [211] "19108_13" "19108_14" "19108_15" "19108_1"  "19108_2"  "19108_3" 
## [217] "19108_5"  "19108_6"  "19108_7"  "19108_8"  "19108_9"  "19127_0" 
## [223] "19127_10" "19127_11" "19127_12" "19127_13" "19127_14" "19127_15"
## [229] "19127_1"  "19127_2"  "19127_3"  "19127_4"  "19127_5"  "19127_6" 
## [235] "19127_7"  "19127_8"  "19127_9"  "19159_0"  "19159_10" "19159_11"
## [241] "19159_12" "19159_13" "19159_14" "19159_15" "19159_1"  "19159_2" 
## [247] "19159_3"  "19159_4"  "19159_5"  "19159_6"  "19159_7"  "19159_8" 
## [253] "19159_9"  "19190_0"  "19190_10" "19190_11" "19190_12" "19190_13"
## [259] "19190_14" "19190_15" "19190_1"  "19190_2"  "19190_3"  "19190_4" 
## [265] "19190_5"  "19190_6"  "19190_7"  "19190_8"  "19190_9"  "19193_0" 
## [271] "19193_10" "19193_11" "19193_12" "19193_13" "19193_14" "19193_15"
## [277] "19193_1"  "19193_2"  "19193_3"  "19193_5"  "19193_6"  "19193_7" 
## [283] "19193_8"  "19193_9"  "19209_0"  "19209_10" "19209_11" "19209_12"
## [289] "19209_14" "19209_15" "19209_1"  "19209_3"  "19209_4"  "19209_5" 
## [295] "19209_6"  "19209_7"  "19209_8"  "19209_9"
```

```r
## 1) Extract the first 20 IDs

id_col2 <- names(x2)[1]      # gene ID column
ids2    <- as.character(x2[[id_col2]])

first20_2 <- ids2[1:20]
first20_2
```

```
##  [1] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" "ENSG00000000938"
##  [5] "ENSG00000000971" "ENSG00000001036" "ENSG00000001084" "ENSG00000001167"
##  [9] "ENSG00000001460" "ENSG00000001461" "ENSG00000001561" "ENSG00000001617"
## [13] "ENSG00000001626" "ENSG00000001629" "ENSG00000001630" "ENSG00000001631"
## [17] "ENSG00000002016" "ENSG00000002330" "ENSG00000002549" "ENSG00000002587"
```

```r
## 2) Create a new column with versions stripped.

x2$ID_stripped <- strip_ensembl_version(ids2)
head(x2[, c(id_col2, "ID_stripped")], 10)
```

```
## # A tibble: 10 × 2
##    Gene_id         ID_stripped    
##    <chr>           <chr>          
##  1 ENSG00000000419 ENSG00000000419
##  2 ENSG00000000457 ENSG00000000457
##  3 ENSG00000000460 ENSG00000000460
##  4 ENSG00000000938 ENSG00000000938
##  5 ENSG00000000971 ENSG00000000971
##  6 ENSG00000001036 ENSG00000001036
##  7 ENSG00000001084 ENSG00000001084
##  8 ENSG00000001167 ENSG00000001167
##  9 ENSG00000001460 ENSG00000001460
## 10 ENSG00000001461 ENSG00000001461
```

```r
## 3) Map the identifiers to HGNC symbols.

ensembl_ids2 <- unique(x2$ID_stripped)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map2 <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids2,
  mart       = ensembl
)

str(map2$ensembl_gene_id)  
```

```
##  chr [1:15956] "ENSG00000000419" "ENSG00000000457" "ENSG00000000460" ...
```

```r
expr_mapped_2 <- x2 %>%
  mutate(ensembl_gene_id = ID_stripped) %>%
  left_join(map2, by = "ensembl_gene_id")  

head(expr_mapped_2[, 1:8])
```

```
## # A tibble: 6 × 8
##   Gene_id       `18489_0` `18489_10` `18489_11` `18489_12` `18489_13` `18489_14`
##   <chr>             <dbl>      <dbl>      <dbl>      <dbl>      <dbl>      <dbl>
## 1 ENSG00000000…       825        549        576        599        607        582
## 2 ENSG00000000…       175        238        279        254        327        274
## 3 ENSG00000000…       337        103        120        124        154        151
## 4 ENSG00000000…         7          1          0          0          0          0
## 5 ENSG00000000…         0          3          5         10          4         43
## 6 ENSG00000001…      2251        860        819        732       1071        860
## # ℹ 1 more variable: `18489_15` <dbl>
```

```r
head(expr_mapped_2[, c("Gene_id", "ID_stripped", "hgnc_symbol")], 20)
```

```
## # A tibble: 20 × 3
##    Gene_id         ID_stripped     hgnc_symbol
##    <chr>           <chr>           <chr>      
##  1 ENSG00000000419 ENSG00000000419 DPM1       
##  2 ENSG00000000457 ENSG00000000457 SCYL3      
##  3 ENSG00000000460 ENSG00000000460 FIRRM      
##  4 ENSG00000000938 ENSG00000000938 FGR        
##  5 ENSG00000000971 ENSG00000000971 CFH        
##  6 ENSG00000001036 ENSG00000001036 FUCA2      
##  7 ENSG00000001084 ENSG00000001084 GCLC       
##  8 ENSG00000001167 ENSG00000001167 NFYA       
##  9 ENSG00000001460 ENSG00000001460 STPG1      
## 10 ENSG00000001461 ENSG00000001461 NIPAL3     
## 11 ENSG00000001561 ENSG00000001561 ENPP4      
## 12 ENSG00000001617 ENSG00000001617 SEMA3F     
## 13 ENSG00000001626 ENSG00000001626 CFTR       
## 14 ENSG00000001629 ENSG00000001629 ANKIB1     
## 15 ENSG00000001630 ENSG00000001630 CYP51A1    
## 16 ENSG00000001631 ENSG00000001631 KRIT1      
## 17 ENSG00000002016 ENSG00000002016 RAD52      
## 18 ENSG00000002330 ENSG00000002330 BAD        
## 19 ENSG00000002549 ENSG00000002549 LAP3       
## 20 ENSG00000002587 ENSG00000002587 HS3ST1
```
4. what is different about this file is that the version numbers are already stripped so my stripped helper function does nothing, none of the ids contain a dot.
For GSE122380, the raw counts file is whitespace delimited, so i had to change the way i read the table using read table instead of read tsv which is the first way i tried.



## Exercise 3 - 
 
Can you use the worked example to process the above two GEO records?  How?

i used the worked example to practice processing the geo records. i was inspired by the workflow and i used the same steps. i downloaded the supplementary data. i used read tsv for the first exercise, however i had to use read table for the second exercise as it was whitespace delimited. then i mapped the ensembl ids the same as the worked example. i took the first column of Ensembl IDs then i used biomart to map to HGNC symbols. for 6.1 i removed the version suffix like in the worked example however it was not necessary for 6.2.



