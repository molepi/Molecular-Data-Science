-----

## Part 1: Bioconductor

Load `IRanges`

``` r
library(IRanges)
```

-----

## Part 2: IRanges

Create `IRanges` object

``` r
ir <- IRanges(5,10)
ir
```

    ## IRanges object with 1 range and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5        10         6

Create `IRanges` object with multiple ranges

``` r
ir <- IRanges(c(4, 5, 8, 15, 19, 28, 40), width=c(15, 6, 7, 12, 9, 3, 6))
ir
```

    ## IRanges object with 7 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         4        18        15
    ##   [2]         5        10         6
    ##   [3]         8        14         7
    ##   [4]        15        26        12
    ##   [5]        19        27         9
    ##   [6]        28        30         3
    ##   [7]        40        45         6

Give the ranges names

``` r
names(ir) = letters[1:7]
ir[c("b", "d", "e")]
```

    ## IRanges object with 3 ranges and 0 metadata columns:
    ##         start       end     width
    ##     <integer> <integer> <integer>
    ##   b         5        10         6
    ##   d        15        26        12
    ##   e        19        27         9

Annotate `IRanges` object

``` r
mcols(ir) <- mtcars[1:7, 1:3]
ir
```

    ## IRanges object with 7 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   b         5        10         6 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258
    ##   e        19        27         9 |      18.7         8       360
    ##   f        28        30         3 |      18.1         6       225
    ##   g        40        45         6 |      14.3         8       360

#### Answers to Question 1

What is the mean `disp` for ranges whose `start` is less than 20?
**209.2**

``` r
mean(mcols(ir[start(ir) < 20])$disp)
```

    ## [1] 209.2

What is the maximum `width` for ranges whose `cyl` is 6? **15**

``` r
max(width(ir[mcols(ir)$cyl==6]))
```

    ## [1] 15

-----

## Part 2

Creating two objects

``` r
ir1 <- ir[c(2,5,7)]
ir2 <- ir[-c(2,5,7)]
ir1
```

    ## IRanges object with 3 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   b         5        10         6 |        21         6       160
    ##   e        19        27         9 |      18.7         8       360
    ##   g        40        45         6 |      14.3         8       360

``` r
ir2
```

    ## IRanges object with 4 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258
    ##   f        28        30         3 |      18.1         6       225

Find overlaps

``` r
olaps <- findOverlaps(ir1, ir2)
olaps
```

    ## Hits object with 3 hits and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         1           1
    ##   [2]         1           2
    ##   [3]         2           3
    ##   -------
    ##   queryLength: 3 / subjectLength: 4

The interpretation of the findOverlaps() output is as follows:

  - The 1st range from ir1 overlaps with the 1st range in ir2 - \[5,10\]
    overlaps with \[4,18\]
  - The 1st range from ir1 overlaps with the 2nd range in ir2 - \[5,10\]
    overlaps with \[8,14\]
  - The 2nd range from ir1 overlaps with the 3rd range in ir2 -
    \[19,27\] overlaps with \[15,26\]

Show ranges using `subjectHits` or `queryHits`

``` r
ir2[subjectHits(olaps)]
```

    ## IRanges object with 3 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258

``` r
ir1[queryHits(olaps)]
```

    ## IRanges object with 3 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   b         5        10         6 |        21         6       160
    ##   b         5        10         6 |        21         6       160
    ##   e        19        27         9 |      18.7         8       360

Count overlaps

``` r
nolaps <- countOverlaps(ir2, ir1)
nolaps
```

    ## a c d f 
    ## 1 1 1 0

``` r
nolaps <- countOverlaps(ir1, ir2)
nolaps
```

    ## b e g 
    ## 2 1 0

Add overlaps to annotation

``` r
mcols(ir2)$Overlaps <- countOverlaps(ir2, ir1)
ir2
```

    ## IRanges object with 4 ranges and 4 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258
    ##   f        28        30         3 |      18.1         6       225
    ##      Overlaps
    ##     <integer>
    ##   a         1
    ##   c         1
    ##   d         1
    ##   f         0

``` r
mcols(ir1)$Overlaps <- countOverlaps(ir1, ir2)
ir1
```

    ## IRanges object with 3 ranges and 4 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   b         5        10         6 |        21         6       160
    ##   e        19        27         9 |      18.7         8       360
    ##   g        40        45         6 |      14.3         8       360
    ##      Overlaps
    ##     <integer>
    ##   b         2
    ##   e         1
    ##   g         0

Subsetting

``` r
subset(ir2, Overlaps > 0)
```

    ## IRanges object with 3 ranges and 4 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258
    ##      Overlaps
    ##     <integer>
    ##   a         1
    ##   c         1
    ##   d         1

#### Answers to Question 2

Load in the example `IRanges`

``` r
load(url("https://raw.githubusercontent.com/molepi/Molecular-Data-Science/master/RIntro_practical/practical3_data_iranges.RData"))
```

How many overlaps are there between `ir1` and `ir2`? **26**

``` r
nolaps <- countOverlaps(ir2, ir1)
sum(nolaps)
```

    ## [1] 26

What is the name of the range in `ir2` with the most overlaps with
`ir1`? **r has 6 overlaps**

``` r
nolaps[nolaps == max(nolaps)]
```

    ## r 
    ## 6

Subset `ir2` to show only ranges that have more than 2 overlaps with
`ir1`

``` r
mcols(ir2)$Overlaps <- nolaps
ir2[mcols(ir2)$Overlaps > 2]
```

    ## IRanges object with 4 ranges and 1 metadata column:
    ##         start       end     width |  Overlaps
    ##     <integer> <integer> <integer> | <integer>
    ##   q        11        15         5 |         3
    ##   r        19        25         7 |         6
    ##   t        27        35         9 |         3
    ##   w        36        48        13 |         4

-----

## Part 4: Genomic Ranges

``` r
library(Homo.sapiens)
hg <- genes(Homo.sapiens)
hg
```

    ## GRanges object with 23056 ranges and 1 metadata column:
    ##         seqnames              ranges strand |       GENEID
    ##            <Rle>           <IRanges>  <Rle> | <FactorList>
    ##       1    chr19   58858172-58874214      - |            1
    ##      10     chr8   18248755-18258723      + |           10
    ##     100    chr20   43248163-43280376      - |          100
    ##    1000    chr18   25530930-25757445      - |         1000
    ##   10000     chr1 243651535-244006886      - |        10000
    ##     ...      ...                 ...    ... .          ...
    ##    9991     chr9 114979995-115095944      - |         9991
    ##    9992    chr21   35736323-35743440      + |         9992
    ##    9993    chr22   19023795-19109967      - |         9993
    ##    9994     chr6   90539619-90584155      + |         9994
    ##    9997    chr22   50961997-50964905      - |         9997
    ##   -------
    ##   seqinfo: 93 sequences (1 circular) from hg19 genome

Remove non-standard chromosomes

``` r
hg <- keepStandardChromosomes(hg, pruning.mode="coarse")
hg[order(seqnames(hg))]
```

    ## GRanges object with 23033 ranges and 1 metadata column:
    ##             seqnames              ranges strand |       GENEID
    ##                <Rle>           <IRanges>  <Rle> | <FactorList>
    ##       10000     chr1 243651535-244006886      - |        10000
    ##   100126331     chr1 117637265-117637350      + |    100126331
    ##   100126346     chr1 154166141-154166219      - |    100126346
    ##   100126348     chr1   94312388-94312467      + |    100126348
    ##   100126349     chr1 166123980-166124035      - |    100126349
    ##         ...      ...                 ...    ... .          ...
    ##       90655     chrY     3447126-3448082      + |        90655
    ##       90665     chrY     6778727-6959724      + |        90665
    ##        9081     chrY   24636544-24660784      + |         9081
    ##        9086     chrY   22737611-22755040      + |         9086
    ##        9087     chrY   15815447-15817902      + |         9087
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg19 genome

-----

#### Answers to Question 3:

How many base pairs does chromosomes 15 span? **102.5 million**

``` r
seqinfo(hg)["chr15"]
```

    ## Seqinfo object with 1 sequence from hg19 genome:
    ##   seqnames seqlengths isCircular genome
    ##   chr15     102531392         NA   hg19

Can you make a table of the number of genes in all of the standard human
chromosomes?

``` r
table(seqnames(hg))
```

    ## 
    ##  chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9 chr10 chr11 chr12 
    ##  2326  1464  1271   873  1022  1000  1108   818   945   903  1439  1173 
    ## chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22  chrX  chrY 
    ##   449   779   791   938  1337   341  1607   647   296   535   918    53 
    ##  chrM 
    ##     0

-----

## Part 5: Overlap with GWAS Hits

Load `gwascat`

``` r
library(gwascat)
```

Load data

``` r
data(ebicat37)
ebicat37
```

    ## gwasloc instance with 22688 records and 36 attributes per record.
    ## Extracted:  2016-01-18 
    ## Genome:  GRCh37 
    ## Excerpt:
    ## GRanges object with 5 ranges and 3 metadata columns:
    ##       seqnames    ranges strand |                  DISEASE/TRAIT
    ##          <Rle> <IRanges>  <Rle> |                    <character>
    ##   [1]    chr11  41820450      * | Post-traumatic stress disorder
    ##   [2]    chr15  35060463      * | Post-traumatic stress disorder
    ##   [3]     chr8  97512977      * | Post-traumatic stress disorder
    ##   [4]     chr9 100983826      * | Post-traumatic stress disorder
    ##   [5]    chr15  54715642      * | Post-traumatic stress disorder
    ##              SNPS   P-VALUE
    ##       <character> <numeric>
    ##   [1]  rs10768747     5e-06
    ##   [2]  rs12232346     2e-06
    ##   [3]   rs2437772     6e-06
    ##   [4]   rs7866350     1e-06
    ##   [5]  rs73419609     6e-06
    ##   -------
    ##   seqinfo: 23 sequences from GRCh37 genome

Change genome

``` r
genome(ebicat37) <- "hg19"
```

Genes that overlap GWAS hits

``` r
goh <- findOverlaps(hg, ebicat37)
goh
```

    ## Hits object with 12934 hits and 0 metadata columns:
    ##           queryHits subjectHits
    ##           <integer>   <integer>
    ##       [1]         2         140
    ##       [2]         4       10865
    ##       [3]         4       16003
    ##       [4]         5        6746
    ##       [5]         5        7962
    ##       ...       ...         ...
    ##   [12930]     23013       11148
    ##   [12931]     23013       11149
    ##   [12932]     23013       22003
    ##   [12933]     23025       19309
    ##   [12934]     23027       21478
    ##   -------
    ##   queryLength: 23033 / subjectLength: 22688

Genes that overlap hits

``` r
length(unique(queryHits(goh)))
```

    ## [1] 4720

Estimate proportion of GWAS SNPs in genes / exons

``` r
mean(reduce(ebicat37) %over% hg)
```

    ## [1] 0.511465

``` r
mean(reduce(ebicat37) %over% exons(Homo.sapiens))
```

    ## [1] 0.07059754
