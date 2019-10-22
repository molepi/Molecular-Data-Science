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

Before moving on, make sure to set your library path.

``` r
.libPaths("C:/fos_2019/library")
```

-----

Unlike CRAN packages, Bioconductor packages are installed differently
using `BiocManager`. You can check if the package you want to use is
available:

``` r
BiocManager::available("IRanges")
```

    ## [1] "IRanges"

All Bioconductor packages have obligatory documentation, called
‘vignettes’. You can view these from within R using:

``` r
browseVignettes("IRanges")
```

This takes you to a list of files that can help you when using this
package in R. In addition to this, help is readily available on the
[Bioconductor](https://support.bioconductor.org) forums.

Now load in the `IRanges` package.

``` r
library(IRanges)
```

-----

## Part 2: IRanges

This package provides efficient low-level and highly reusable classes
for storing and manipulating **annotated ranges of integers**. It also
contains useful functions, such as to identify overlaps.

You can imagine how this kind of functionality would be useful when
handling sequences. `IRanges` forms the basis for one of the most
important packages handling genomic data, `GenomicRanges`, which we will
explore later.

`IRanges` objects can be created by specifying an integer range.

``` r
ir <- IRanges(5,10)
ir
```

    ## IRanges object with 1 range and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5        10         6

We have created an integer range from 5 to 10, and the width is listed
as 6 since both the start and end integer contribute to the width of an
`IRanges` object.

We can access this information using the `start()`, `end()`, and
`width()` functions.

``` r
start(ir)
```

    ## [1] 5

``` r
end(ir)
```

    ## [1] 10

``` r
width(ir)
```

    ## [1] 6

-----

`IRanges` can be shifted or resized. Here, `shift()` performs an intra
range transformation and reduces the start and end integer by 2, keeping
the width the same.

``` r
shift(ir, -2)
```

    ## IRanges object with 1 range and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         3         8         6

The `resize()` function adjusts the width to a new specified value. Its
default behaviour is to keep the `start` the same, and adjust the width
using the `end` integer.

``` r
resize(ir, 3)
```

    ## IRanges object with 1 range and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5         7         3

-----

You can input multiple ranges at once using vectors to create a
multi-range object.

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

These can also be resized or shifted, and this will be performed on each
range individually.

``` r
shift(ir,1)
```

    ## IRanges object with 7 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         5        19        15
    ##   [2]         6        11         6
    ##   [3]         9        15         7
    ##   [4]        16        27        12
    ##   [5]        20        28         9
    ##   [6]        29        31         3
    ##   [7]        41        46         6

-----

The `reduce()` function will merge any overlapping regions and create an
`IRanges` object of only distinct ranges. Here, there was a lot of
overlap, with the only gap between 30 and 40.

``` r
reduce(ir)
```

    ## IRanges object with 2 ranges and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]         4        30        27
    ##   [2]        40        45         6

Any gaps can also be displayed easily. None of our `IRanges` object
overlapped with integers 31 to 39.

``` r
gaps(ir)
```

    ## IRanges object with 1 range and 0 metadata columns:
    ##           start       end     width
    ##       <integer> <integer> <integer>
    ##   [1]        31        39         9

-----

Each row in an `IRanges` object can be given a name, and they can be
subsetted similarly to vectors or data frames.

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

Associating each range with a collection of attributes is very useful
for genomic annotation. This metadata is assigned using the `mcols()`
function.

R has some standard built-in data sets. One, called `mtcars`, is used
here to demonstrate annotation.

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

This annotation will remain the same through any `shift()` or `resize()`
commands. Specific `mcol` variables can be accessed using the `$`
operator, similar to in data frames.

There are over 300 functions that can be used to manipulate `IRanges`
objects. In the next section, their role in genome-scale analysis will
become clearer.

-----

#### Question 1: Using the above created `IRanges` object:

  - What is the mean `disp` for ranges whose `start` is less than 20?
  - What is the maximum `width` for ranges whose `cyl` is 6?

-----

## Part 3: Finding Overlaps

It is very useful to count overlaps in two distinct `IRanges` objects.
We can subset the `ir` object above to create two new `IRanges` objects.

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

-----

We can then use `findOverlaps()` to identify overlapping ranges in the
`query` (first argument) and `subject` (second argument) objects.

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

The interpretation of the `findOverlaps()` output is as follows:

  - The 1st range from `ir1` overlaps with the 1st range in `ir2` -
    \[5,10\] overlaps with \[4,18\]
  - The 1st range from `ir1` overlaps with the 2nd range in `ir2` -
    \[5,10\] overlaps with \[8,14\]
  - The 2nd range from `ir1` overlaps with the 3rd range in `ir2` -
    \[19,27\] overlaps with \[15,26\]

-----

You can list the ranges with overlaps in a particular `IRanges` object
after using `findOverlaps()` by utilising subsetting and `queryHits()`
or `subjectHits()`.

So, we can show the ranges within `ir2` that overlap with ranges in
`ir1`:

``` r
ir1[queryHits(olaps)]
```

    ## IRanges object with 3 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   b         5        10         6 |        21         6       160
    ##   b         5        10         6 |        21         6       160
    ##   e        19        27         9 |      18.7         8       360

``` r
ir2[subjectHits(olaps)]
```

    ## IRanges object with 3 ranges and 3 metadata columns:
    ##         start       end     width |       mpg       cyl      disp
    ##     <integer> <integer> <integer> | <numeric> <numeric> <numeric>
    ##   a         4        18        15 |        21         6       160
    ##   c         8        14         7 |      22.8         4       108
    ##   d        15        26        12 |      21.4         6       258

-----

You can also `countOverlaps()`.

``` r
nolapsir1 <- countOverlaps(ir1, ir2)
nolapsir1
```

    ## b e g 
    ## 2 1 0

``` r
nolapsir2 <- countOverlaps(ir2, ir1)
nolapsir2
```

    ## a c d f 
    ## 1 1 1 0

So, here we see that in `ir1`, `b` overlaps with two ranges in `ir2`,
`e` with one, and `g` with none.

Also, in `ir2`, `a`, `c`, and `d` ranges in `ir2` overlap with ranges in
`ir1` once each, and `f` does not overlap with any ranges in `ir1`.

-----

You can add output from `countOverlaps()` into your `IRanges`, coupling
derived data with the original annotation.

``` r
mcols(ir1)$Overlaps <- nolapsir1
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

``` r
mcols(ir2)$Overlaps <- nolapsir2
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

-----

After this, you can perform co-ordinated actions, like subsetting for
transcripts satisfying particular conditions.

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

-----

#### Question 2: Load the example `IRanges` into R using the following command:

`load(url("https://raw.githubusercontent.com/molepi/Molecular-Data-Science/master/RIntro_practical/practical3_data_iranges.RData"))`

  - How many overlaps are there between `ir1` and `ir2`?
  - What is the name of the range in `ir2` with the most overlaps with
    `ir1`?
  - Subset `ir2` to show only ranges that have more than 2 overlaps with
    `ir1`

-----

## Part 4: Genomic Ranges

As an extension to `IRanges`, `GenomicRanges` objects contain obligatory
metadata called `seqnames`. This gives the chromosome occupied by the
sequence whose `start` and `end` position are modelled by the associated
`IRanges`.

The package `Homo.sapiens` is an annotation package from Bioconductor.
It contains information on genes in the human genome that can be
accessed using the `genes()` function.

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

This `GenomicRanges` object functions similarly to an `IRanges` object.
You can explore `start()`, `end()`, and `width()` in the same way, and
`subset()` it as before.

There are two parts to a `GenomicRanges` object. The first part is the
`seqnames`, which must consist at least of the start and end coordinates
along with the strand. The second part are additional elements or
meta-data, like `GENEID`.

-----

You can order a `GenomicRanges` object to make exploring it more
intuitive. Here, we sort the ranges by chromosome using `seqnames()`.

``` r
hg[order(seqnames(hg))]
```

    ## GRanges object with 23056 ranges and 1 metadata column:
    ##                   seqnames              ranges strand |       GENEID
    ##                      <Rle>           <IRanges>  <Rle> | <FactorList>
    ##       10000           chr1 243651535-244006886      - |        10000
    ##   100126331           chr1 117637265-117637350      + |    100126331
    ##   100126346           chr1 154166141-154166219      - |    100126346
    ##   100126348           chr1   94312388-94312467      + |    100126348
    ##   100126349           chr1 166123980-166124035      - |    100126349
    ##         ...            ...                 ...    ... .          ...
    ##      283788 chrUn_gl000219         53809-99642      - |       283788
    ##   100507412 chrUn_gl000220        97129-126696      + |    100507412
    ##   100288687 chrUn_gl000228       112605-124171      + |    100288687
    ##   100653046 chrUn_gl000228         86060-90745      + |    100653046
    ##      728410 chrUn_gl000228        79448-113887      + |       728410
    ##   -------
    ##   seqinfo: 93 sequences (1 circular) from hg19 genome

-----

An important component of `GenomicRanges`is `seqinfo()`. This shows how
many base pairs are in each chromosome as `seqlengths` and which genome
they are from.

``` r
seqinfo(hg)
```

    ## Seqinfo object with 93 sequences (1 circular) from hg19 genome:
    ##   seqnames       seqlengths isCircular genome
    ##   chr1            249250621       <NA>   hg19
    ##   chr2            243199373       <NA>   hg19
    ##   chr3            198022430       <NA>   hg19
    ##   chr4            191154276       <NA>   hg19
    ##   chr5            180915260       <NA>   hg19
    ##   ...                   ...        ...    ...
    ##   chrUn_gl000245      36651       <NA>   hg19
    ##   chrUn_gl000246      38154       <NA>   hg19
    ##   chrUn_gl000247      36422       <NA>   hg19
    ##   chrUn_gl000248      39786       <NA>   hg19
    ##   chrUn_gl000249      38502       <NA>   hg19

-----

You can inspect `seqinfo()` for a specific chromosome. The human X
chromosome spans more than 155 million base pairs.

``` r
seqinfo(hg)["chrX"]
```

    ## Seqinfo object with 1 sequence from hg19 genome:
    ##   seqnames seqlengths isCircular genome
    ##   chrX      155270560         NA   hg19

-----

You can `subset()` a `GenomicRanges` object to return only ranges in a
specific chromosome, `seqnames()` alongside the `%in%` operator. Here,
we see there are 296 genes on chromosome 21.

``` r
subset(hg, seqnames %in% "chr21")
```

    ## GRanges object with 296 ranges and 1 metadata column:
    ##             seqnames            ranges strand |       GENEID
    ##                <Rle>         <IRanges>  <Rle> | <FactorList>
    ##   100129027    chr21 47247755-47256333      - |    100129027
    ##   100131902    chr21 31661463-31661832      - |    100131902
    ##   100132288    chr21   9907189-9968593      - |    100132288
    ##   100133286    chr21 37441940-37498938      - |    100133286
    ##   100151643    chr21 31992946-31993169      + |    100151643
    ##         ...      ...               ...    ... .          ...
    ##        9619    chr21 43619799-43717354      + |         9619
    ##        9875    chr21 33683330-33765312      - |         9875
    ##        9946    chr21 34961648-35014160      - |         9946
    ##        9980    chr21 37536839-37666572      + |         9980
    ##        9992    chr21 35736323-35743440      + |         9992
    ##   -------
    ##   seqinfo: 93 sequences (1 circular) from hg19 genome

-----

A useful function for `GenomicRanges` is the `keepStandardChromosomes()`
function. This removes unmapped or non-standard chromosomes for the
species of interest from the object.

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

If you compare this with the above output, you can see that ranges from
chromosomes like `chrUn_gl000228` have been removed.

-----

#### Question 3:

  - How many base pairs does chromosome 15 span?
  - Can you make a `table()` of the number of genes in all of the
    standard human chromosomes?

-----

## Part 5: Overlaps with GWAS Hits

Genome-wide association studies (GWAS) are systematically represented in
the EMBL-EBI GWAS catalog. The `gwascat` Bioconductor package allows us
to retrieve a version of this.

``` r
library(gwascat)
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

Per GWAS finding, there is one range, and each corresponds to a SNP that
has beem replicated as significantly associated with the given
phenotype.

The genome here is listed as `GRCh37`, which is very similar to `hg19`.
To avoid an error, we need to rename it.

``` r
genome(ebicat37) <- "hg19"
```

-----

You can view the top traits in the catalog using `topTraits()`.

``` r
topTraits(ebicat37)
```

    ## 
    ##  Obesity-related traits                  Height       IgG glycosylation 
    ##                     957                     822                     699 
    ##         Type 2 diabetes    Rheumatoid arthritis         Crohn's disease 
    ##                     340                     294                     249 
    ##           Schizophrenia Blood metabolite levels         HDL cholesterol 
    ##                     248                     245                     220 
    ##           Breast cancer 
    ##                     199

If we are interested in a specific trait, like LDL cholesterol, we can
`subsetByTraits()`.

``` r
subsetByTraits(ebicat37, tr="LDL cholesterol")
```

    ## gwasloc instance with 184 records and 36 attributes per record.
    ## Extracted:  2016-01-18 
    ## Genome:  hg19 
    ## Excerpt:
    ## GRanges object with 5 ranges and 3 metadata columns:
    ##       seqnames    ranges strand |   DISEASE/TRAIT        SNPS   P-VALUE
    ##          <Rle> <IRanges>  <Rle> |     <character> <character> <numeric>
    ##   [1]     chr2  44072576      * | LDL cholesterol   rs4299376     4e-72
    ##   [2]     chr1 150958836      * | LDL cholesterol    rs267733     5e-09
    ##   [3]     chr2  21263900      * | LDL cholesterol   rs1367117    1e-182
    ##   [4]    chr19  45422946      * | LDL cholesterol   rs4420638    2e-178
    ##   [5]    chr17  64210580      * | LDL cholesterol   rs1801689     1e-11
    ##   -------
    ##   seqinfo: 23 sequences from hg19 genome

-----

We can use `findOverlaps()` as above to determine which genes have
annotated intervals that overlap GWAS hits.

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

To make more sense of this, we can count the number of genes that
overlap with one or more GWAS hits using `length()` coupled with
`unique()`.

``` r
length(unique(queryHits(goh)))
```

    ## [1] 4720

-----

Lastly, implementing the `reduce()` function, we can estimate the
proportion of GWAS SNPs that lie within gene regions or, more
specifically, exons.

``` r
mean(reduce(ebicat37) %over% hg)
```

    ## [1] 0.511465

``` r
mean(reduce(ebicat37) %over% exons(Homo.sapiens))
```

    ## [1] 0.07059754

We can see that the vast majority of GWAS hits lie in non-coding regions
and almost half are intergenic. This is the focus of a lot of research,
since the regulatory mechanisms of non-coding regions is considerably
complex and elusive.

-----

## Feedback

Hopefully, this 2 day R course was informative and helpful. Please send
any questions or comments to <l.j.sinke@lumc.nl>.

It is obligatory to hand in to the Blackboard one file, but if you have
several please turn in one and e-mail me the others. Thank you.

Feedback is very much appreciated. Enjoy the rest of this Molecular Data
Science FOS course\! :)

-----
