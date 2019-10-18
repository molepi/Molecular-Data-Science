This final practical introduces Bioconductor and some basic analysis and
annotations.

-----

## Part 1: Bioconductor

[Bioconductor](https://www.bioconductor.org/) is a collection of more
than 1,500 packages for analysing and annotating high-throughput genomic
data.

It is useful for a wide range of applications including bulk and
single-cell RNA sequencing, copy number analysis, flow cytometry, and
methylation microarray data.

-----

Unlike CRAN packages, Bioconductor packages are installed differently
using `BiocManager`. You can check if the package you want to use is
available:

``` r
BiocManager::available("GenomicRanges")
```

    ## [1] "GenomicRanges"

This package imports an object class called `GenomicRanges`, which
contains genomic data, such as SNPs, gene expression, or DNA
methylation.

-----

All Bioconductor packages have obligatory documentation, called
‘vignettes’. You can view these from within R using:

``` r
browseVignettes("GenomicRanges")
```

This takes you to a list of files that can help you when using this
package in R. In addition to this, help is readily available on the
[Bioconductor](https://support.bioconductor.org) forums.

-----

## Part 2: Genomic Ranges

Another type of package in R contains highly curated data, as resources
for teaching and learning. `DMRcatedata` is one such package. We can
load this package into R and import the data, which is an object of
class `RangedSummarizedExperiment`

``` r
library(DMRcatedata)
data(dmrcatedata)
```

This loads in a series of objects into R, including several
`GenomicRanges` class objects. tx.hg19 is a GRanges object containing
simulated WGBS data.

``` r
tx.hg19
```

    ## GRanges object with 215170 ranges and 4 metadata columns:
    ##                                   seqnames              ranges strand |
    ##                                      <Rle>           <IRanges>  <Rle> |
    ##   ENST00000000233                     chr7 127228399-127231759      + |
    ##   ENST00000000412                    chr12     9092961-9102551      - |
    ##   ENST00000000442                    chr11   64073050-64084210      + |
    ##   ENST00000001008                    chr12     2904119-2913124      + |
    ##   ENST00000001146                     chr2   72356367-72375167      - |
    ##               ...                      ...                 ...    ... .
    ##   ENST00000610276                    chr21   33108045-33108720      + |
    ##   ENST00000610277 chrHSCHR19LRC_LRC_I_CTG1   54677109-54693666      - |
    ##   ENST00000610278                    chr22   21335650-21336044      - |
    ##   ENST00000610279                    chr10   69609283-69610504      + |
    ##   ENST00000610280                    chr11   58059298-58060237      - |
    ##                         gene_name      gene_type         gene_id
    ##                       <character>    <character>     <character>
    ##   ENST00000000233            ARF5 protein_coding ENSG00000004059
    ##   ENST00000000412            M6PR protein_coding ENSG00000003056
    ##   ENST00000000442           ESRRA protein_coding ENSG00000173153
    ##   ENST00000001008           FKBP4 protein_coding ENSG00000004478
    ##   ENST00000001146         CYP26B1 protein_coding ENSG00000003137
    ##               ...             ...            ...             ...
    ##   ENST00000610276      AP000255.6        lincRNA ENSG00000273091
    ##   ENST00000610277          MBOAT7 protein_coding ENSG00000273130
    ##   ENST00000610278 XXbac-B135H6.18        lincRNA ENSG00000272829
    ##   ENST00000610279    RP11-57G10.8        lincRNA ENSG00000272892
    ##   ENST00000610280         OR10Q2P     pseudogene ENSG00000272900
    ##                               tx_name
    ##                           <character>
    ##   ENST00000000233            ARF5-001
    ##   ENST00000000412            M6PR-001
    ##   ENST00000000442           ESRRA-002
    ##   ENST00000001008           FKBP4-001
    ##   ENST00000001146         CYP26B1-001
    ##               ...                 ...
    ##   ENST00000610276      AP000255.6-001
    ##   ENST00000610277          MBOAT7-001
    ##   ENST00000610278 XXbac-B135H6.18-001
    ##   ENST00000610279    RP11-57G10.8-001
    ##   ENST00000610280         OR10Q2P-001
    ##   -------
    ##   seqinfo: 265 sequences from an unspecified genome; no seqlengths

-----

For convenience, let’s work only with the ‘standard’ autosomal
chromosomes.

``` r
tx.hg19 <- keepStandardChromosomes(tx.hg19, pruning.mode="coarse")
tx.hg19
```

    ## GRanges object with 196317 ranges and 4 metadata columns:
    ##                   seqnames              ranges strand |       gene_name
    ##                      <Rle>           <IRanges>  <Rle> |     <character>
    ##   ENST00000000233     chr7 127228399-127231759      + |            ARF5
    ##   ENST00000000412    chr12     9092961-9102551      - |            M6PR
    ##   ENST00000000442    chr11   64073050-64084210      + |           ESRRA
    ##   ENST00000001008    chr12     2904119-2913124      + |           FKBP4
    ##   ENST00000001146     chr2   72356367-72375167      - |         CYP26B1
    ##               ...      ...                 ...    ... .             ...
    ##   ENST00000610274     chr2 171197675-171207618      - |      AC012594.1
    ##   ENST00000610276    chr21   33108045-33108720      + |      AP000255.6
    ##   ENST00000610278    chr22   21335650-21336044      - | XXbac-B135H6.18
    ##   ENST00000610279    chr10   69609283-69610504      + |    RP11-57G10.8
    ##   ENST00000610280    chr11   58059298-58060237      - |         OR10Q2P
    ##                        gene_type         gene_id             tx_name
    ##                      <character>     <character>         <character>
    ##   ENST00000000233 protein_coding ENSG00000004059            ARF5-001
    ##   ENST00000000412 protein_coding ENSG00000003056            M6PR-001
    ##   ENST00000000442 protein_coding ENSG00000173153           ESRRA-002
    ##   ENST00000001008 protein_coding ENSG00000004478           FKBP4-001
    ##   ENST00000001146 protein_coding ENSG00000003137         CYP26B1-001
    ##               ...            ...             ...                 ...
    ##   ENST00000610274      antisense ENSG00000231898      AC012594.1-029
    ##   ENST00000610276        lincRNA ENSG00000273091      AP000255.6-001
    ##   ENST00000610278        lincRNA ENSG00000272829 XXbac-B135H6.18-001
    ##   ENST00000610279        lincRNA ENSG00000272892    RP11-57G10.8-001
    ##   ENST00000610280     pseudogene ENSG00000272900         OR10Q2P-001
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome; no seqlengths

-----

There are two parts to a `GenomicRanges` object.

The first part is the `seqnames`, which must consist at least of the
start and end coordinates along with the strand. Here, this part also
shows the chromosome.

The second part are additional elements, like `gene_name`, `gene_type`,
and `gene_id`.

-----

Components can be accessed when needed:

``` r
head(start(tx.hg19))
```

    ## [1] 127228399   9092961  64073050   2904119  72356367  37458809

``` r
head(width(tx.hg19))
```

    ## [1]  3361  9591 11161  9006 18801 17495

-----

You can `subset()` these objects, for example if you want to focus on a
particular set of chromosomes.

``` r
subset(tx.hg19, seqnames %in% c("chr1", "chr2"))
```

    ## GRanges object with 31510 ranges and 4 metadata columns:
    ##                   seqnames              ranges strand |     gene_name
    ##                      <Rle>           <IRanges>  <Rle> |   <character>
    ##   ENST00000001146     chr2   72356367-72375167      - |       CYP26B1
    ##   ENST00000002125     chr2   37458809-37476303      + |       NDUFAF7
    ##   ENST00000003583     chr1   24683495-24740215      - |         STPG1
    ##   ENST00000003912     chr1   24742293-24799466      + |        NIPAL3
    ##   ENST00000005756     chr2 158958382-158992478      + |          UPP2
    ##               ...      ...                 ...    ... .           ...
    ##   ENST00000610254     chr2 162102956-162111141      - |    AC009299.2
    ##   ENST00000610262     chr2   86422713-86423172      + | RP11-301O19.1
    ##   ENST00000610265     chr2 145276024-145277950      + |      ZEB2-AS1
    ##   ENST00000610272     chr1 179850742-179851730      - | RP11-533E19.7
    ##   ENST00000610274     chr2 171197675-171207618      - |    AC012594.1
    ##                        gene_type         gene_id           tx_name
    ##                      <character>     <character>       <character>
    ##   ENST00000001146 protein_coding ENSG00000003137       CYP26B1-001
    ##   ENST00000002125 protein_coding ENSG00000003509       NDUFAF7-001
    ##   ENST00000003583 protein_coding ENSG00000001460         STPG1-001
    ##   ENST00000003912 protein_coding ENSG00000001461        NIPAL3-001
    ##   ENST00000005756 protein_coding ENSG00000007001          UPP2-001
    ##               ...            ...             ...               ...
    ##   ENST00000610254      antisense ENSG00000235724    AC009299.2-005
    ##   ENST00000610262      antisense ENSG00000273080 RP11-301O19.1-001
    ##   ENST00000610265      antisense ENSG00000238057      ZEB2-AS1-010
    ##   ENST00000610272        lincRNA ENSG00000272906 RP11-533E19.7-001
    ##   ENST00000610274      antisense ENSG00000231898    AC012594.1-029
    ##   -------
    ##   seqinfo: 24 sequences from an unspecified genome; no seqlengths

-----

There are also packages on Bioconductor containing Annotation data, such
as the TxDb.\* family of packages. These contain information on the
genomic coordinates of exons, genes, transcripts, and so on.

``` r
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
tx <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
tx
```

    ## GRanges object with 82960 ranges and 2 metadata columns:
    ##                 seqnames        ranges strand |     tx_id     tx_name
    ##                    <Rle>     <IRanges>  <Rle> | <integer> <character>
    ##       [1]           chr1   11874-14409      + |         1  uc001aaa.3
    ##       [2]           chr1   11874-14409      + |         2  uc010nxq.1
    ##       [3]           chr1   11874-14409      + |         3  uc010nxr.1
    ##       [4]           chr1   69091-70008      + |         4  uc001aal.1
    ##       [5]           chr1 321084-321115      + |         5  uc001aaq.2
    ##       ...            ...           ...    ... .       ...         ...
    ##   [82956] chrUn_gl000237        1-2686      - |     82956  uc011mgu.1
    ##   [82957] chrUn_gl000241   20433-36875      - |     82957  uc011mgv.2
    ##   [82958] chrUn_gl000243   11501-11530      + |     82958  uc011mgw.1
    ##   [82959] chrUn_gl000243   13608-13637      + |     82959  uc022brq.1
    ##   [82960] chrUn_gl000247     5787-5816      - |     82960  uc022brr.1
    ##   -------
    ##   seqinfo: 93 sequences (1 circular) from hg19 genome

-----

As above, we want to work only with autosomal chromosomes.

``` r
tx <- keepStandardChromosomes(tx, pruning.mode="coarse")
tx
```

    ## GRanges object with 78827 ranges and 2 metadata columns:
    ##           seqnames        ranges strand |     tx_id     tx_name
    ##              <Rle>     <IRanges>  <Rle> | <integer> <character>
    ##       [1]     chr1   11874-14409      + |         1  uc001aaa.3
    ##       [2]     chr1   11874-14409      + |         2  uc010nxq.1
    ##       [3]     chr1   11874-14409      + |         3  uc010nxr.1
    ##       [4]     chr1   69091-70008      + |         4  uc001aal.1
    ##       [5]     chr1 321084-321115      + |         5  uc001aaq.2
    ##       ...      ...           ...    ... .       ...         ...
    ##   [78823]     chrM    7587-15888      - |     78823  uc022bqs.1
    ##   [78824]     chrM    8367-14149      - |     78824  uc022bqt.1
    ##   [78825]     chrM   10761-14149      - |     78825  uc031tgb.1
    ##   [78826]     chrM   14675-14698      - |     78826  uc022bqv.1
    ##   [78827]     chrM   15960-16024      - |     78827  uc022bqx.1
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg19 genome

-----

## Part 3: Finding Overlaps

It is very useful to count overlaps in two distinct `GenomicRanges`
objects.

``` r
olaps <- countOverlaps(tx, tx.hg19)
xtabs( ~olaps)
```

    ## olaps
    ##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
    ## 3099 5925 4878 4867 5055 5065 5018 4860 4383 4231 3943 3377 3124 2603 2390 
    ##   15   16   17   18   19   20   21   22   23   24   25   26   27   28   29 
    ## 2003 1714 1535 1298 1152 1093 1057  772  676  551  532  396  471  319  270 
    ##   30   31   32   33   34   35   36   37   38   39   40   41   42   43   44 
    ##  281  264  205  173  170   84  100   54   81   78   78   73   57   68   35 
    ##   45   46   47   48   49   50   51   52   53   54   55   56   57   58   59 
    ##   15   16   49   60   34   10    3    7   18    5   17    6   16    3    7 
    ##   60   61   62   63   64   65   66   67   68   69   71   74   77   78   79 
    ##   12    4    3    2    3    2    9    1    8    7    1    1    1   14    4 
    ##   80   82   89   92   96   98  103  104  107  111  116  119  127  132  207 
    ##    1   11    3    1    1    1    5    1    1    1    1    1    1    1    1

-----

You can add output from `countOverlaps()` into your `GenomicRanges`,
coupling derived data with the original annotation.

``` r
tx$Overlaps <- countOverlaps(tx, tx.hg19)
tx
```

    ## GRanges object with 78827 ranges and 3 metadata columns:
    ##           seqnames        ranges strand |     tx_id     tx_name  Overlaps
    ##              <Rle>     <IRanges>  <Rle> | <integer> <character> <integer>
    ##       [1]     chr1   11874-14409      + |         1  uc001aaa.3         4
    ##       [2]     chr1   11874-14409      + |         2  uc010nxq.1         4
    ##       [3]     chr1   11874-14409      + |         3  uc010nxr.1         4
    ##       [4]     chr1   69091-70008      + |         4  uc001aal.1         1
    ##       [5]     chr1 321084-321115      + |         5  uc001aaq.2         4
    ##       ...      ...           ...    ... .       ...         ...       ...
    ##   [78823]     chrM    7587-15888      - |     78823  uc022bqs.1         0
    ##   [78824]     chrM    8367-14149      - |     78824  uc022bqt.1         0
    ##   [78825]     chrM   10761-14149      - |     78825  uc031tgb.1         0
    ##   [78826]     chrM   14675-14698      - |     78826  uc022bqv.1         0
    ##   [78827]     chrM   15960-16024      - |     78827  uc022bqx.1         0
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg19 genome

-----

After this, you can perform co-ordinated actions, like subsetting for
transcripts satisfying particular conditions.

``` r
subset(tx, Overlaps > 10)
```

    ## GRanges object with 27503 ranges and 3 metadata columns:
    ##           seqnames            ranges strand |     tx_id     tx_name
    ##              <Rle>         <IRanges>  <Rle> | <integer> <character>
    ##       [1]     chr1     762971-794826      + |        15  uc001abp.2
    ##       [2]     chr1     762971-794826      + |        16  uc009vjn.2
    ##       [3]     chr1     762971-794826      + |        17  uc009vjo.2
    ##       [4]     chr1     762971-794826      + |        18  uc021oem.2
    ##       [5]     chr1     763178-794826      + |        19  uc031pjk.1
    ##       ...      ...               ...    ... .       ...         ...
    ##   [27499]     chrY 15434982-15592550      - |     78731  uc022cny.1
    ##   [27500]     chrY 15435435-15592550      - |     78732  uc022cnz.1
    ##   [27501]     chrY 15447443-15592550      - |     78733  uc022coa.1
    ##   [27502]     chrY 24026223-24329089      - |     78764  uc010nxc.1
    ##   [27503]     chrY 24049765-24329089      - |     78769  uc011nbg.2
    ##            Overlaps
    ##           <integer>
    ##       [1]        11
    ##       [2]        11
    ##       [3]        11
    ##       [4]        11
    ##       [5]        11
    ##       ...       ...
    ##   [27499]        12
    ##   [27500]        12
    ##   [27501]        12
    ##   [27502]        19
    ##   [27503]        15
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg19 genome

## Part 4: Summarized Experiments

Another class of object often used by Bioconductor packages is the
\``SummarizedExperiment`. This is used to store matrices of experimental
results, commonly produced by sequencing and microarray experiments.

Each object stores multiple samples along with meta-data, which
describes features and phenotypes. The anatomy of a
`SummarizedExperiment` object is shown here:

![Figure
1](https://github.com/molepi/Molecular-Data-Science/blob/master/RIntro_practical/summarizedexperiment.svg)

-----

Another Bioconductor package that contains teaching data is `airway`.
Lets load it and then read in the data.

``` r
library(airway)
data(airway)
airway
```

    ## class: RangedSummarizedExperiment 
    ## dim: 64102 8 
    ## metadata(1): ''
    ## assays(1): counts
    ## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
    ## rowData names(0):
    ## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    ## colData names(9): SampleName cell ... Sample BioSample

-----

Phenotype and feature information is stored in the `colData`, which is a
data frame class object. For example, `cell` represents which cell line
was used.

``` r
colData(airway)
```

    ## DataFrame with 8 rows and 9 columns
    ##            SampleName     cell      dex    albut        Run avgLength
    ##              <factor> <factor> <factor> <factor>   <factor> <integer>
    ## SRR1039508 GSM1275862   N61311    untrt    untrt SRR1039508       126
    ## SRR1039509 GSM1275863   N61311      trt    untrt SRR1039509       126
    ## SRR1039512 GSM1275866  N052611    untrt    untrt SRR1039512       126
    ## SRR1039513 GSM1275867  N052611      trt    untrt SRR1039513        87
    ## SRR1039516 GSM1275870  N080611    untrt    untrt SRR1039516       120
    ## SRR1039517 GSM1275871  N080611      trt    untrt SRR1039517       126
    ## SRR1039520 GSM1275874  N061011    untrt    untrt SRR1039520       101
    ## SRR1039521 GSM1275875  N061011      trt    untrt SRR1039521        98
    ##            Experiment    Sample    BioSample
    ##              <factor>  <factor>     <factor>
    ## SRR1039508  SRX384345 SRS508568 SAMN02422669
    ## SRR1039509  SRX384346 SRS508567 SAMN02422675
    ## SRR1039512  SRX384349 SRS508571 SAMN02422678
    ## SRR1039513  SRX384350 SRS508572 SAMN02422670
    ## SRR1039516  SRX384353 SRS508575 SAMN02422682
    ## SRR1039517  SRX384354 SRS508576 SAMN02422673
    ## SRR1039520  SRX384357 SRS508579 SAMN02422683
    ## SRR1039521  SRX384358 SRS508580 SAMN02422677

-----

Inside these classes of objects, `GenomicRanges` store genomic
information. These can be accessed using `rowRanges()`

``` r
rowRanges(airway)
```

    ## GRangesList object of length 64102:
    ## $ENSG00000000003 
    ## GRanges object with 17 ranges and 2 metadata columns:
    ##        seqnames            ranges strand |   exon_id       exon_name
    ##           <Rle>         <IRanges>  <Rle> | <integer>     <character>
    ##    [1]        X 99883667-99884983      - |    667145 ENSE00001459322
    ##    [2]        X 99885756-99885863      - |    667146 ENSE00000868868
    ##    [3]        X 99887482-99887565      - |    667147 ENSE00000401072
    ##    [4]        X 99887538-99887565      - |    667148 ENSE00001849132
    ##    [5]        X 99888402-99888536      - |    667149 ENSE00003554016
    ##    ...      ...               ...    ... .       ...             ...
    ##   [13]        X 99890555-99890743      - |    667156 ENSE00003512331
    ##   [14]        X 99891188-99891686      - |    667158 ENSE00001886883
    ##   [15]        X 99891605-99891803      - |    667159 ENSE00001855382
    ##   [16]        X 99891790-99892101      - |    667160 ENSE00001863395
    ##   [17]        X 99894942-99894988      - |    667161 ENSE00001828996
    ## 
    ## ...
    ## <64101 more elements>
    ## -------
    ## seqinfo: 722 sequences (1 circular) from an unspecified genome

-----

Assay information can be accessed using `assay()`. FOr example, here we
see a table of read counts for different exons in each sample.

``` r
head(assay(airway))
```

    ##                 SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ## ENSG00000000003        679        448        873        408       1138
    ## ENSG00000000005          0          0          0          0          0
    ## ENSG00000000419        467        515        621        365        587
    ## ENSG00000000457        260        211        263        164        245
    ## ENSG00000000460         60         55         40         35         78
    ## ENSG00000000938          0          0          2          0          1
    ##                 SRR1039517 SRR1039520 SRR1039521
    ## ENSG00000000003       1047        770        572
    ## ENSG00000000005          0          0          0
    ## ENSG00000000419        799        417        508
    ## ENSG00000000457        331        233        229
    ## ENSG00000000460         63         76         60
    ## ENSG00000000938          0          0          0

-----

Subsetting a `SummarizedExperiment` object is quite simple. Suppose we
want to look only at samples treated with dexamethasone.

``` r
subset(airway, , dex=="trt")
```

    ## class: RangedSummarizedExperiment 
    ## dim: 64102 4 
    ## metadata(1): ''
    ## assays(1): counts
    ## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
    ## rowData names(0):
    ## colnames(4): SRR1039509 SRR1039513 SRR1039517 SRR1039521
    ## colData names(9): SampleName cell ... Sample BioSample

-----

We can also extract interested summaries from these objects. For
instance, we can calculate the total number of reads overlapping genes
in each sample and store them in our `SummarizedExperiment` class
object.

``` r
airway$libSize <- colSums(assay(airway))
airway$libSize
```

    ## SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 
    ##   20637971   18809481   25348649   15163415   24448408   30818215 
    ## SRR1039520 SRR1039521 
    ##   19126151   21164133

-----

## Part 5: Down-stream analysis

It is important to look at how skills learned here translate to working
with objects in Bioconductor packages.

`DESeq2` is a package for analysing bulk RNAseq differential expression
data. Load it into R.

``` r
library(DESeq2)
```

This package utilizes a class that combines `SummarizedExperiment` type
count data with `formula` objects that describe the experimental design.

-----

We may want to include cell line as a covariate, and investigate the
effect of dexamethasone treatment.

``` r
dSeqData <- DESeqDataSet(airway, design = ~ cell + dex)
dSeqData
```

    ## class: DESeqDataSet 
    ## dim: 64102 8 
    ## metadata(2): '' version
    ## assays(1): counts
    ## rownames(64102): ENSG00000000003 ENSG00000000005 ... LRG_98 LRG_99
    ## rowData names(0):
    ## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    ## colData names(10): SampleName cell ... BioSample libSize

This class of object can be manipulated in a very similar way to
`SummarizedExperiment` objects.

-----

The `DESeq` workflow is summarized by a single function call, which
performs advanced statistical analysis on a `DESeqDataSet` class object.

``` r
dSeqRes <- DESeq(dSeqData)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

-----

We can now extract a table of differential expression, which we can
visualise or manipulate further.

``` r
head(results(dSeqRes))
```

    ## log2 fold change (MLE): dex untrt vs trt 
    ## Wald test p-value: dex untrt vs trt 
    ## DataFrame with 6 rows and 6 columns
    ##                          baseMean      log2FoldChange             lfcSE
    ##                         <numeric>           <numeric>         <numeric>
    ## ENSG00000000003  708.602169691234   0.381253982523046 0.100654411867397
    ## ENSG00000000005                 0                  NA                NA
    ## ENSG00000000419  520.297900552084  -0.206812601085046  0.11221864650836
    ## ENSG00000000457  237.163036796015 -0.0379204312050866  0.14344465493375
    ## ENSG00000000460  57.9326331250967  0.0881681758708956 0.287141822933389
    ## ENSG00000000938 0.318098378392895    1.37822703433214  3.49987280259267
    ##                              stat               pvalue                padj
    ##                         <numeric>            <numeric>           <numeric>
    ## ENSG00000000003  3.78775232451127 0.000152016273081236 0.00128362968201804
    ## ENSG00000000005                NA                   NA                  NA
    ## ENSG00000000419 -1.84294328545158   0.0653372915124159   0.196546085664247
    ## ENSG00000000457 -0.26435583272587     0.79150574160785   0.911459479590442
    ## ENSG00000000460 0.307054454729672    0.758801923977344   0.895032784968828
    ## ENSG00000000938 0.393793463954221    0.693733530342183                  NA

-----

## Feedback

Hopefully, this 2 day R course was informative and helpful. Please send
any questions or comments to <l.j.sinke@lumc.nl>.

Feedback is very much appreciated. Enjoy the rest of this Molecular Data
Science FOS course\! :)

-----
