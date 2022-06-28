# Gene Biotypes

6-27-2022

Mira Sohn

## Overview

This workflow is written to demonstrate how to explore [gene biotypes](https://useast.ensembl.org/info/genome/genebuild/biotypes.html) in DESeq2-mediated RNA-seq workflow. This analysis used [GSE206057](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206057) available on [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) as an input dataset. This workflow will minimize to talk about theoretical and technical background. Also, many important steps performed in the `DESeq2`-mediated RNA-seq workflow will be skipped in order to focus on analyzing gene biotypes. Please visit manual pages listed below if you are analyzing real-world datasets.


## Analysis

### Loading libraries

```r
library(tidyverse)
library(data.table)
library(ggplot2)
library(DESeq2)
library(AnnotationHub)
library(ensembldb)
```

In addition to tools for data importing, cleaning, and plotting, this workflow uses following tools for genomic analysis.

* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): RNA-seq pipeline
* [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html), [ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html): used to retrieve gene biotypes from Ensembl DB

### Assigning input/output file paths

```r

# Assign input file path
input.path <- "GSE206057_Results.genes.deseq.counts.xlsx"

# Assign output file directory
output.path <- "plots"

# Create the output directory if absent
if (!dir.exists(output.path)) { dir.create(output.path) }

```

This step is optional but makes the analysis way more error-free than hard-coding. Users can use a separate R or yaml file instead.


### Loading input dataset

```r
# Import read count matrix (returns a data frame)
read.count <- readxl::read_excel(input.path, sheet=1)
```

The input count matrix was provided as an excel file. Given that, this workflow imports the input matrix using the package [`readxl`](https://readxl.tidyverse.org/). 



### Data cleaning

In the current workflow, data cleaning is aimed at converting the input data frame to `DESeq2`-loadable data shape.

#### Exploring input data frame

The functions `head()` and `str()` were used to explore the input data frame.

```r
# > head(read.count)
# # A tibble: 6 × 10
#   ...1            DMSO1  DMSO2   DMSO3   LCX1   LCX2   LCX3   HCX1   HCX2   HCX3
#   <chr>           <dbl>  <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
# 1 0610005C13Rik    0       0      0    2.02e0 0         0   0         0   0
# 2 0610009B22Rik  931.   1010.   951.   8.25e2 8.84e2  868.  9.75e2  893.  8.98e2
# 3 0610009E02Rik    3.43   14.2    9.81 1.08e1 7.12e0   17.3 7.87e0   14.8 9.71e0
# 4 0610009L18Rik   34.5    31.7   40.2  3.97e1 1.46e1   38.1 3.35e1   24.0 2.51e1
# 5 0610010B08Rik  232.    463.   283.   3.04e2 2.66e2  264.  5.20e2  407.  5.36e2
# 6 0610010F05Rik 2785.   2477.  2831.   2.86e3 2.77e3 2692.  2.79e3 2694.  2.66e3

# > str(read.count)
# tibble [24,972 × 10] (S3: tbl_df/tbl/data.frame)
#  $ ...1 : chr [1:24972] "0610005C13Rik" "0610009B22Rik" "0610009E02Rik" "0610009L18Rik" ...
#  $ DMSO1: num [1:24972] 0 930.63 3.43 34.52 231.98 ...
#  $ DMSO2: num [1:24972] 0 1010.5 14.2 31.7 462.7 ...
#  $ DMSO3: num [1:24972] 0 950.98 9.81 40.23 282.86 ...
#  $ LCX1 : num [1:24972] 2.02 825.06 10.83 39.74 303.52 ...
#  $ LCX2 : num [1:24972] 0 884.36 7.12 14.55 266.08 ...
#  $ LCX3 : num [1:24972] 0 868.2 17.3 38.1 263.9 ...
#  $ HCX1 : num [1:24972] 0 974.88 7.87 33.48 520.18 ...
#  $ HCX2 : num [1:24972] 0 892.9 14.8 24 407.1 ...
#  $ HCX3 : num [1:24972] 0 898.49 9.71 25.15 536.01 ...
```

#### Reformatting dataset

The `DESeq2` requires specific input data format accepted by the arguments below:

- `countData`: consisting of columns per sample and rows per gene
- `colData`: consisting of columns per condition (e.g. samplename, genotype, treatment, biological replicates, etc) and rows per sample

Check out the [DESeq2 tutorial](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) for more info.


```r
# Rename the first column
colnames(read.count)[1] <- "geneid"


# Clean the count matrix and conver the object from data frame to matrix
read.count <- read.count %>%
    dplyr::filter(!(geneid %in% c("44621", "44622"))) %>%
    column_to_rownames("geneid") %>%
    as.matrix()

# convert the counts to integers
read.count <- apply(read.count, 2, round)

# dim(read.count)
# [1] 24968     9

# Filter out undetected genes
read.count <- read.count[rowSums(read.count) > 0, ]

# > dim(read.count)
# [1] 17181     9



# > head(read.count)
#               DMSO1 DMSO2 DMSO3 LCX1 LCX2 LCX3 HCX1 HCX2 HCX3
# 0610005C13Rik     0     0     0    2    0    0    0    0    0
# 0610009B22Rik   931  1010   951  825  884  868  975  893  898
# 0610009E02Rik     3    14    10   11    7   17    8   15   10
# 0610009L18Rik    35    32    40   40   15   38   33   24   25
# 0610010B08Rik   232   463   283  304  266  264  520  407  536
# 0610010F05Rik  2785  2477  2831 2859 2770 2692 2790 2694 2658

# Build colData
colData <- data.frame(samplename=colnames(read.count)) %>%
    mutate(treatment=c(rep("DMSO", 3), rep("LCX", 3), rep("HCX", 3)),
           replicate=rep(1:3, 3))




# Assign rownames
rownames(colData) <- colData$samplename

# > colData
#       samplename treatment replicate
# DMSO1      DMSO1      DMSO         1
# DMSO2      DMSO2      DMSO         2
# DMSO3      DMSO3      DMSO         3
# LCX1        LCX1       LCX         1
# LCX2        LCX2       LCX         2
# LCX3        LCX3       LCX         3
# HCX1        HCX1       HCX         1
# HCX2        HCX2       HCX         2
# HCX3        HCX3       HCX         3
```


### Building DESeq2 object


In addition to the arguments `countData` and `colData`, the `design` is also required for building `DESeq2` object named `dds`. Note that it has to match one of the column names in the `colData`. This demonstration will compare the `treatment` conditions.

```r
# Build dds obj
# NOTE: Row names of the colData and column names of the countData should match.
dds <- DESeqDataSetFromMatrix(countData=read.count,
                              colData=colData,
                              design=~treatment)
```


### Extracting normalized read counts

For convenience, QC is skipped in this workflow. Check out the [DESeq2 doc](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) for details about QC. For real-world analysis, QC using the function `vst()` is highly recommended in DESeq2-mediated RNA-seq.


In DESeq2 workflow, raw read counts are normalized by being divided by size factors which are computed per sample. The size factors reflect sequencing depth. Therefore, computing size factors using the function `estimateSizeFactors()` is performed prior to extracting the normalized counts. Run `sizeFactors(dds)` if you'd like to explore your size factors.

```r

# Estimate size factors
dds <- estimateSizeFactors(dds)


# Extract normalized read counts as a data frame
# NOTE: `counts(dds, normalized=TRUE)` gives counts scaled by size or normalization factors.
norm.counts <- as.data.frame(counts(dds, normalized=TRUE))

# Clean the normalized count data frame
# 1. Remove zero-count rows
# 2. Move rownames to column
norm.counts <- norm.counts[rowSums(norm.counts) > 0,] %>%
    rownames_to_column("Genename")

# > head(norm.counts)
#        Genename       DMSO1      DMSO2       DMSO3        LCX1        LCX2
# 1 0610005C13Rik    0.000000    0.00000    0.000000    2.016367    0.000000
# 2 0610009B22Rik  911.637375  989.35418  929.840851  831.751290  890.065586
# 3 0610009E02Rik    2.937607   13.71382    9.777506   11.090017    7.048031
# 4 0610009L18Rik   34.272082   31.34588   39.110025   40.327335   15.102923
# 5 0610010B08Rik  227.174942  453.53563  276.703429  306.487748  267.825165
# 6 0610010F05Rik 2727.078507 2426.36664 2768.012039 2882.396290 2789.006417
#         LCX3        HCX1       HCX2       HCX3
# 1    0.00000    0.000000    0.00000    0.00000
# 2  874.87565  982.828433  900.34353  904.01896
# 3   17.13466    8.064233   15.12335   10.06703
# 4   38.30101   33.264962   24.19736   25.16757
# 5  266.09121  524.175164  410.34694  539.59261
# 6 2713.32401 2812.401363 2716.15395 2675.81559


```


### Exploring gene biotypes

The [Ensembl](https://useast.ensembl.org/info/genome/genebuild/biotypes.html) database defines **biotype** as a gene or transcript classification.

#### Retrieving annotation data from ensembl DB

This workflow retrieves gene biotypes per gene from Ensembl DB using the package `AnnotationHub`.

References:
- [Bioconductor AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
- [AnnotationHub doc1](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub-HOWTO.html)
- [AnnotationHub doc2](https://bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub.html)
- [Bioconductor ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html)
- [ensembldb doc](https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html)

```r

# Create an AnnotationHub obj
ah <- AnnotationHub()

# Query Ensembl DB for mouse
ensdb <- query(ah, pattern=c("Ensdb", "Mus musculus"))

# > ensdb
# AnnotationHub with 19 records
# # snapshotDate(): 2021-10-20
# # $dataprovider: Ensembl
# # $species: Mus musculus
# # $rdataclass: EnsDb
# # additional mcols(): taxonomyid, genome, description,
# #   coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
# #   rdatapath, sourceurl, sourcetype
# # retrieve records with, e.g., 'object[["AH53222"]]'
#
#             title
#   AH53222 | Ensembl 87 EnsDb for Mus Musculus
#   AH53726 | Ensembl 88 EnsDb for Mus Musculus
#   AH56691 | Ensembl 89 EnsDb for Mus Musculus
#   AH57770 | Ensembl 90 EnsDb for Mus Musculus
#   AH60788 | Ensembl 91 EnsDb for Mus Musculus
#   ...       ...
#   AH83247 | Ensembl 101 EnsDb for Mus musculus
#   AH89211 | Ensembl 102 EnsDb for Mus musculus
#   AH89457 | Ensembl 103 EnsDb for Mus musculus
#   AH95775 | Ensembl 104 EnsDb for Mus musculus
#   AH98078 | Ensembl 105 EnsDb for Mus musculus

# > mcols(ensdb) %>% rownames()
#  [1] "AH53222" "AH53726" "AH56691" "AH57770" "AH60788" "AH60992" "AH64461"
#  [8] "AH64944" "AH67971" "AH69210" "AH73905" "AH75036" "AH78811" "AH79718"
# [15] "AH83247" "AH89211" "AH89457" "AH95775" "AH98078"



# Retrieve the newest ensdb
newest.ensdb <- length(mcols(ensdb) %>% rownames())
ensdb <- ah[[ (mcols(ensdb) %>% rownames())[newest.ensdb] ]]

# NOTE: If you'd like to keep the version constantly, assign a specific name (e.g. AH98078) instead of using the newest db.
# ensdb <- ah[["AH98078"]]
```

#### Adding gene biotypes to the normalized count data

Normalized read counts per biotype per treatment is calculated.

```r

# Retrieve gene biotypes matching the input genes
biotype.df <- select(ensdb,
                     keys=norm.counts$Genename,
                     keytype="GENENAME",
                     columns=c("GENENAME", "GENEID", "GENEBIOTYPE"))


# > head(biotype.df)
#        GENENAME             GENEID    GENEBIOTYPE
# 1 0610005C13Rik ENSMUSG00000109644         lncRNA
# 2 0610009B22Rik ENSMUSG00000007777 protein_coding
# 3 0610009E02Rik ENSMUSG00000086714         lncRNA
# 4 0610009L18Rik ENSMUSG00000043644         lncRNA
# 5 0610010K14Rik ENSMUSG00000020831 protein_coding
# 6 0610012G03Rik ENSMUSG00000107002 protein_coding


# > unique(biotype.df$GENEBIOTYPE)
#  [1] "lncRNA"                             "protein_coding"
#  [3] "TEC"                                "transcribed_processed_pseudogene"
#  [5] "transcribed_unprocessed_pseudogene" "processed_pseudogene"
#  [7] "pseudogene"                         "transcribed_unitary_pseudogene"
#  [9] "unprocessed_pseudogene"             "polymorphic_pseudogene"
# [11] "miRNA"                              "ribozyme"
# [13] "snRNA"                              "scaRNA"
# [15] "snoRNA"                             "misc_RNA"


# Add gene biotypes and sample metadata by joining data tables
cleaned.df <- norm.counts %>%
    inner_join(biotype.df, by=c("Genename"="GENENAME")) %>%  # join retrieved biotype data frame to the count data frame
    gather('Sample', 'Count', colnames(norm.counts)[-1]) %>%  # reshape the data frame by gathering all the counts from individual samples to a single column
    inner_join(colData, by=c("Sample"="samplename")) %>%   # add sample metadata by joining the data frame `colData`
    group_by(treatment, GENEBIOTYPE) %>%   # reshape the data frame to contain counts per treatment and biotype
    summarize(Count=mean(Count))  # Count is mean of biological replicates

# head(cleaned.df)
# A tibble: 6 × 3
# Groups:   treatment [1]
# treatment GENEBIOTYPE               Count
# <chr>     <chr>                     <dbl>
# DMSO      lncRNA                  145.
# DMSO      miRNA                     0.457
# DMSO      misc_RNA                  2.29
# DMSO      polymorphic_pseudogene  511.
# DMSO      processed_pseudogene   1744.
# DMSO      protein_coding         2286.


# Rename columns
colnames(cleaned.df) <- c('Treatment', 'Gene_biotype', 'Count')

# head(cleaned.df)
# A tibble: 6 × 3
# Groups:   Treatment [1]
# Treatment Gene_biotype              Count
# <chr>     <chr>                     <dbl>
# DMSO      lncRNA                  145.
# DMSO      miRNA                     0.457
# DMSO      misc_RNA                  2.29
# DMSO      polymorphic_pseudogene  511.
# DMSO      processed_pseudogene   1744.
# DMSO      protein_coding         2286.



```




#### Plotting count of gene biotypes

A stacked bar plot is created using `ggplot2` package to present normalized read counts per biotype across the treatment conditions.

```r
# Plot
p <- ggplot(cleaned.df, aes(x=Treatment, y=Count, fill=Gene_biotype)) +
    geom_bar(position="stack", stat="identity", width=0.5, aes(fill=Gene_biotype)) +
    theme_bw() +
    ylab("Normalized Read Counts")

# Print
print(p)


```

![biotype_count.png](https://github.com/Mira0507/genebiotype/blob/master/plots/biotype_count.png)


```r
# Save
# NOTE: Assign `device=pdf` and `filename="biotype_count.pdf" if you'd like to save as pdf
ggsave(filename=file.path(output.path, "biotype_count.png"),
       plot=p,
       device="png")

```


#### Plotting proportion of gene biotypes

Proportion: 100 x (normalized counts per biotype per treatment condition) / (normalized counts per treatment condition)

```r

# Clean the data frame and calculate percentage
proportion.df <- cleaned.df %>%
    # Reshape the data frame to have summed counts per treatment
    group_by(Treatment) %>%
    summarize(Total_count=sum(Count)) %>%
    # Add normalized counts per gene biotype by joining another `cleaned.df`
    inner_join(cleaned.df, by='Treatment') %>%
    # Add a new column containing percentage counts
    mutate(Proportion=round(100*Count/Total_count, 2))



# head(proportion.df)
# A tibble: 6 × 5
# Treatment Total_count Gene_biotype              Count Proportion
# <chr>           <dbl> <chr>                     <dbl>      <dbl>
# DMSO            7120. lncRNA                  145.          2.04
# DMSO            7120. miRNA                     0.457       0.01
# DMSO            7120. misc_RNA                  2.29        0.03
# DMSO            7120. polymorphic_pseudogene  511.          7.18
# DMSO            7120. processed_pseudogene   1744.         24.5
# DMSO            7120. protein_coding         2286.         32.1

# Plot
p <- ggplot(proportion.df, aes(x=Treatment, y=Proportion, fill=Gene_biotype)) +
    geom_bar(position="stack", stat="identity", width=0.5, aes(fill=Gene_biotype)) +
    theme_bw() +
    ylab("Proportion (%)")


# Print
print(p)


```

![biotype_proportion.png](https://github.com/Mira0507/genebiotype/blob/readme/plots/biotype_proportion.png)




```r
# Save
ggsave(filename=file.path(output.path, "biotype_proportion.png"),
       plot=p,
       device="png")

```


#### Plotting complexity of gene biotypes

Complexity indicates total number of unique gene IDs per biotype.


```r

# Clean the data frame
cleaned.df <- norm.counts %>%
    # join retrieved biotype data frame to the count data frame
    inner_join(biotype.df, by=c("Genename"="GENENAME")) %>%
    # reshape the data frame by gathering all the counts from individual samples to a single column
    gather('Sample', 'Count', colnames(norm.counts)[-1]) %>%
    # add sample metadata by joining the data frame `colData`
    inner_join(colData, by=c("Sample"="samplename")) %>%
    # reshape the data frame to contain counts per treatment and biotype
    group_by(treatment, GENEBIOTYPE) %>%
    # count is mean of biological replicates and total number of unique genes per biotype per condition
    summarize(Count=mean(Count), nGene=length(unique(GENEID))) %>%
    # add computed proportion data by joining data frames
    inner_join(proportion.df,
               by=c("treatment"="Treatment", "GENEBIOTYPE"="Gene_biotype", "Count"="Count"))


# head(cleaned.df)
# A tibble: 6 × 6
# Groups:   treatment [1]
# treatment GENEBIOTYPE               Count nGene Total_count Proportion
# <chr>     <chr>                     <dbl> <int>       <dbl>      <dbl>
# DMSO      lncRNA                  145.     1030       7120.       2.04
# DMSO      miRNA                     0.457     5       7120.       0.01
# DMSO      misc_RNA                  2.29      1       7120.       0.03
# DMSO      polymorphic_pseudogene  511.       16       7120.       7.18
# DMSO      processed_pseudogene   1744.       80       7120.      24.5
# DMSO      protein_coding         2286.    15348       7120.      32.1

# Rename the columns
colnames(cleaned.df) <- c("Treatment",
                          "Gene_biotype",
                          "Count",
                          "nGene",
                          "Total_count",
                          "Proportion")

# head(cleaned.df)
# A tibble: 6 × 6
# Groups:   treatment [1]
# treatment GENEBIOTYPE               Count nGene Total_count Proportion
# <chr>     <chr>                     <dbl> <int>       <dbl>      <dbl>
# DMSO      lncRNA                  145.     1030       7120.       2.04
# DMSO      miRNA                     0.457     5       7120.       0.01
# DMSO      misc_RNA                  2.29      1       7120.       0.03
# DMSO      polymorphic_pseudogene  511.       16       7120.       7.18
# DMSO      processed_pseudogene   1744.       80       7120.      24.5
# DMSO      protein_coding         2286.    15348       7120.      32.1

```

A scatter plot is created using `ggplot2` to display total number of unique gene IDs per biotype (differentiated by color) per treatment condition (differentiated by shape) on x-axis along with proportion of each biotype on y-axis.

```r
# Plot
p <- ggplot(cleaned.df,
            aes(x=log10(nGene), y=Proportion, shape=Treatment, color=Gene_biotype)) +
    geom_point(size=5, alpha=0.7) +
    theme_bw() +
    theme(legend.title=element_blank()) +
    xlab("log10 (Number of Unique Gene IDs Expressed") +
    ylab("Proportion (%)")

# Print
```


![biotype_complexity.png](https://github.com/Mira0507/genebiotype/blob/master/plots/biotype_complexity.png)
