These introductory practicals are designed to teach you the basics of R.
At the end, we expect you to be able to follow along with the other
practicals of this course. Save your R-script or paste your answers in a
document turn it in on the blackboard.

-----

## Part 1: Finding your way around RStudio

Open RStudio

After opening a new R-script with `ctrl + shift + N`, you should have a
layout similar to the image below.

![Figure
1](https://github.com/molepi/Molecular-Data-Science/blob/master/RIntro_practical/figure.png?raw=true)

You can type equations in the *Console* (see Figure) and R will
calculate the result. Try this by typing the following equation in the
Console:

``` r
1 + 1     
```

    ## [1] 2

-----

You can also type your code in the *R-script* (see Figure) and run a
line or a selection using `ctrl + enter`. This is preferable to typing
in the *Console* directly, since you can save the scripts for later use
or share them with collaborators. Additionally, scripts can be commented
using `#` to improve their readability.

Copy the following equations in your *R-script*, place your cursor on
the first line, and run the equations line-by-line. You can also select
multiple lines and run them together.

    5 - 2         # subtraction
    2 * 2         # multiplication
    6 / 2         # division
    3 ^ 2         # exponentiation
    sqrt(9)       # square root
    abs(-5)       # absolute value

-----

Using the `up arrow` and `down arrow`, you can scroll through historical
commands. These can also be viewed in the **History** tab of the
*Environment* section (see Figure).

In the **Environment** tab of the *Environment* section, you can view
any stored variables available for use.

You can store values in variables using the assignment operator, `<-`.
Note that variable names in R are case-sensitive.

``` r
x <- 1 + 1
x
```

    ## [1] 2

You can list all currently saved variables, or remove them from memory.
Remember, `rm()` is irreversible. If you want to use the variable again,
you will need to redefine it.

``` r
ls()        # list all stored variables
```

    ## [1] "x"

``` r
rm(x)       # removing x
ls()
```

    ## character(0)

``` r
x <- 1 + 1  # redefining x
```

-----

In the *bottom right* section of RStudio, are some very useful tabs.
**Files** will display the content of the current working directory. You
can set this using `setwd()` and view the path using `getwd()`.

Any **Plots** that you create can also be viewed here, and you can
browse through your plot history with the `Forward` and `Back` buttons.

The `tidyverse` package is a collection of popular R packages that can
be installed together. It contains packages like `ggplot2`, meantioned
in the lecture, for creating graphics, and `dplyr`, which we will use
today for data manipulation.

A list of currently installed packages is shown in the **Packages** tab.
The functions contained within these can be loaded into R with the
`library()` function.

Let’s load `tidyverse`, so that we have access to it later:

``` r
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------------------------------------------------------------ tidyverse 1.2.1 --

    ## v ggplot2 3.2.1     v purrr   0.3.2
    ## v tibble  2.1.3     v dplyr   0.8.3
    ## v tidyr   1.0.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.4.0

    ## -- Conflicts --------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

We will use this package later in the practical. After loading a package
into R, the functions become available for use.

Lastly, the **Help** tab can be used to view documentation on any loaded
function. `Tab` can be used as autocomplete in R, but can also be used
to browse options for completion. Find help on the function, `mean`.

    ?mean

#### Question 1:

  - What does the `na.rm` option of the `mean` function do?
  - What classes of input are accepted by the `sqrt` function?

-----

## Part 2: Vectors

You can store multiple values in a vector using the concatenate
function, `c()`.

``` r
y <- c(5, 6, 7)
y
```

    ## [1] 5 6 7

You can perform calculations on elements of a vector, using the same
operators and functions you saw above.

``` r
x + y
```

    ## [1] 7 8 9

``` r
y * 3
```

    ## [1] 15 18 21

``` r
y ^ 2
```

    ## [1] 25 36 49

-----

You can also sum the elements of a vector using `sum()`, view the
maximum value using `max()`, or return the number of elements using
`length()`. Additionally, vectors can be reversed using `rev()`.

``` r
sum(y)
```

    ## [1] 18

``` r
max(y)
```

    ## [1] 7

``` r
length(y)
```

    ## [1] 3

``` r
rev(y)
```

    ## [1] 7 6 5

-----

You can select specific values from a vector using square brackets. A
sequence of values can be defined with `:`. Values can also be added to
currently defined vectors using `c()`.

``` r
y[1]        # Return the first value of y
```

    ## [1] 5

``` r
y[2:3]      # Return the second and third value of y
```

    ## [1] 6 7

``` r
c(y, 8)     # Add an 8 as the last element of y
```

    ## [1] 5 6 7 8

``` r
y[-3]       # Return all except the third value of y
```

    ## [1] 5 6

-----

Some functions can help you explore the nature of stored variables.

``` r
str(y)
```

    ##  num [1:3] 5 6 7

``` r
class(y)
```

    ## [1] "numeric"

-----

In this case, `y` is a numeric class object. A vector can also contain
character data, defined using quotation marks.

``` r
z <- c("a", "b", "c")
z
```

    ## [1] "a" "b" "c"

``` r
c(z, "d")
```

    ## [1] "a" "b" "c" "d"

``` r
class(z)
```

    ## [1] "character"

-----

Another class of data in R is the logical class, where `T` represents
TRUE and `F` represents FALSE.

``` r
w <- c(T,T,F)
w
```

    ## [1]  TRUE  TRUE FALSE

``` r
class(w)
```

    ## [1] "logical"

-----

In addition to the above methods, there are other ways to subset vectors
in R. These can either use numerical or logical vectors defined by `c()`
or logical operators, such as `>`, `<=`, or `!=`.

``` r
y[c(1, 3)]    # Return the first and third value of y
```

    ## [1] 5 7

``` r
z[w]          # Uses a logical vector to subset z
```

    ## [1] "a" "b"

``` r
y[y == 6]     # Return all values of y equal to 6
```

    ## [1] 6

``` r
z[z != "a"]   # Return all values of z not equal to "a"
```

    ## [1] "b" "c"

``` r
y[y > 5]      # Return all values of y greater than 5
```

    ## [1] 6 7

-----

#### Question 2: Obtain the following vectors by adding, removing or subsetting vectors y or z. There are multiple ways to do this, but you only need to give one of them in your answers.

  - `7`
  - `"a", "c"`
  - `5, 7, 9`
  - `"c", "b", "a"`

-----

## Part 3: Data Frames

R also has a data frame class which is comparable to a spreadsheet.

You can click on a data frame in the **Environment** to view it.
However, when the data frame is large it is better to use functions such
as `str`, `summary`, `head` and `dim`.

``` r
df <- data.frame(y, z)  # Create a data frame from vectors y and z
df
```

    ##   y z
    ## 1 5 a
    ## 2 6 b
    ## 3 7 c

``` r
class(df)
```

    ## [1] "data.frame"

``` r
str(df)
```

    ## 'data.frame':    3 obs. of  2 variables:
    ##  $ y: num  5 6 7
    ##  $ z: Factor w/ 3 levels "a","b","c": 1 2 3

``` r
summary(df)
```

    ##        y       z    
    ##  Min.   :5.0   a:1  
    ##  1st Qu.:5.5   b:1  
    ##  Median :6.0   c:1  
    ##  Mean   :6.0        
    ##  3rd Qu.:6.5        
    ##  Max.   :7.0

-----

You can view the dimensions of a data frame.

``` r
nrow(df)
```

    ## [1] 3

``` r
ncol(df)
```

    ## [1] 2

``` r
dim(df)
```

    ## [1] 3 2

-----

You can add row names and column names to a data frame.

``` r
rownames(df) <- c("row1", "row2", "row3")
colnames(df) <- c("col1", "col2")
df
```

    ##      col1 col2
    ## row1    5    a
    ## row2    6    b
    ## row3    7    c

-----

You can subset data frames in a similar manner to vectors, or using the
`$` operator.

``` r
df[1, ]           # Return the first row
```

    ##      col1 col2
    ## row1    5    a

``` r
df[, 1]           # Return the first column
```

    ## [1] 5 6 7

``` r
df[1, 1]          # Return the first value in the first column
```

    ## [1] 5

``` r
df$col1           # Return column 'col1'
```

    ## [1] 5 6 7

``` r
df[df$col1 > 5, ] # Return rows where column 'col1' is greater than 5
```

    ##      col1 col2
    ## row2    6    b
    ## row3    7    c

-----

Extra columns can be added to a data frame using the `$` operator.

``` r
df$col3 <- c(T, T, F) 
df
```

    ##      col1 col2  col3
    ## row1    5    a  TRUE
    ## row2    6    b  TRUE
    ## row3    7    c FALSE

Adding extra rows requires the `rbind()` function, and for values of the
new row to be of the correct class.

``` r
df2 <- data.frame(col1 = 8, col2 = "d", col3 = F) 
rownames(df2) <- "row4" 
df3 <- rbind(df, df2) 
df3
```

    ##      col1 col2  col3
    ## row1    5    a  TRUE
    ## row2    6    b  TRUE
    ## row3    7    c FALSE
    ## row4    8    d FALSE

-----

#### Question 3:

  - Add a column with the name `col4` to data frame `df` with the values
    `1, 1, 1`.
  - Make a new data frame `df4` which is a subset of data frame `df` but
    only contains columns `col1` and `col2`, and rows `row2` and `row3`.
  - Calculate the sum of `col1` of data frame `df3`.
  - Change the row names of data frame `df4` to `rowX` and `rowY`.

-----

## Part 4: Real Data

We are now going to load a data set from the Framingham heart study
[BioLINCC](https://biolincc.nhlbi.nih.gov/home/). This data is stored as
a .csv file, so we read it into R using `read.csv()`.

``` r
fhs <- read.csv(url("https://raw.githubusercontent.com/molepi/Molecular-Data-Science/master/RIntro_practical/data.csv"))
```

The first few lines of the data frame include:

``` r
summary(fhs)
```

    ##      SEX          TOTCHOL         AGE            SYSBP      
    ##  Female:2490   Min.   :107   Min.   :32.00   Min.   : 83.5  
    ##  Male  :1944   1st Qu.:206   1st Qu.:42.00   1st Qu.:117.5  
    ##                Median :234   Median :49.00   Median :129.0  
    ##                Mean   :237   Mean   :49.93   Mean   :132.9  
    ##                3rd Qu.:264   3rd Qu.:57.00   3rd Qu.:144.0  
    ##                Max.   :696   Max.   :70.00   Max.   :295.0  
    ##                NA's   :52                                   
    ##      DIABP        CURSMOKE        BMI        DIABETES    BPMEDS    
    ##  Min.   : 48.00   No :2253   Min.   :15.54   No :4313   No  :4229  
    ##  1st Qu.: 75.00   Yes:2181   1st Qu.:23.09   Yes: 121   Yes : 144  
    ##  Median : 82.00              Median :25.45              NA's:  61  
    ##  Mean   : 83.08              Mean   :25.85                         
    ##  3rd Qu.: 90.00              3rd Qu.:28.09                         
    ##  Max.   :142.50              Max.   :56.80                         
    ##                              NA's   :19                            
    ##     GLUCOSE            EDUC         MI      
    ##  Min.   : 40.00   Min.   :1.000   No :3703  
    ##  1st Qu.: 72.00   1st Qu.:1.000   Yes: 731  
    ##  Median : 78.00   Median :2.000             
    ##  Mean   : 82.19   Mean   :1.976             
    ##  3rd Qu.: 87.00   3rd Qu.:3.000             
    ##  Max.   :394.00   Max.   :4.000             
    ##  NA's   :397      NA's   :113

``` r
head(fhs)
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1   Male     195  39 106.0    70       No 26.97       No     No      77
    ## 2 Female     250  46 121.0    81       No 28.73       No     No      76
    ## 3   Male     245  48 127.5    80      Yes 25.34       No     No      70
    ## 4 Female     225  61 150.0    95      Yes 28.58       No     No     103
    ## 5 Female     285  46 130.0    84      Yes 23.10       No     No      85
    ## 6 Female     228  43 180.0   110       No 30.30       No     No      99
    ##   EDUC  MI
    ## 1    4 Yes
    ## 2    2  No
    ## 3    1  No
    ## 4    3  No
    ## 5    3  No
    ## 6    2 Yes

From looking at the above, we can see that this data frame contains 4434
observations of 12 variables. The first 11 variables represent values
measured at baseline (e.g. sex, age), and the last variable, MI,
represents whether of not the individual had suffered a myocardial
infarction at 24 years after baseline.

-----

There are 642 missing values in the data. While there are various
methods to handle these, we will look only at complete cases in this
instance by removing observations with missing values.

``` r
sum(is.na(fhs))
```

    ## [1] 642

``` r
fhs <- fhs[complete.cases(fhs), ]
sum(is.na(fhs))
```

    ## [1] 0

-----

A very useful concept in R is the `formula`. This is specified as `lhs ~
rhs`, where the left-hand side is typically a dependent variable, and
the right-hand side represents the independent variables, combined using
`+`.

``` r
f <- y ~ x + z
class(f)
```

    ## [1] "formula"

-----

We make use of formulae and create tables of variables in our data using
`xtabs()` and `prop.table()`.

``` r
xtabs(~ MI, fhs)              # Table of MI in FHS
```

    ## MI
    ##   No  Yes 
    ## 3208  643

``` r
prop.table(xtabs(~ MI, fhs))  # Proportion table of MI in FHS
```

    ## MI
    ##        No       Yes 
    ## 0.8330304 0.1669696

``` r
xtabs(~ MI + SEX, fhs)        # Table of MI by SEX
```

    ##      SEX
    ## MI    Female Male
    ##   No    1905 1303
    ##   Yes    203  440

``` r
xtabs(~ MI, fhs[fhs$AGE>50, ])# Table of MI in individuals over 50
```

    ## MI
    ##   No  Yes 
    ## 1392  372

-----

#### Question 4:

  - Make a table of current smoking status
  - Make a proportion table of BP medication
  - Make a table showing the number of smokers for each education level
  - Cross tabulate sex and diabetes status

-----

## Part 5: Data Manipulation

As mentioned above, `dplyr` is a `tidyverse` package that is frequently
used for data summary and transformation. `dplyr` code is often more
intuitive to read and use that base R packages for data manipulation.

This package contains a number of functions, which we will introduce one
by one and then combine.

-----

The `select()` function allows you to subset specific columns out of a
dataset, similar to subsetting introduced above.

``` r
fhsBP <- select(fhs, SYSBP, DIABP, BPMEDS)
fhsBP
```

    ##      SYSBP DIABP BPMEDS
    ## 1    106.0  70.0     No
    ## 2    121.0  81.0     No
    ## 3    127.5  80.0     No
    ## 4    150.0  95.0     No
    ## 5    130.0  84.0     No
    ## 6    180.0 110.0     No
    ## 7    138.0  71.0     No
    ## 8    100.0  71.0     No
    ## 9    141.5  89.0     No
    ## 10   162.0 107.0     No
    ## 11   133.0  76.0     No
    ## 12   131.0  88.0     No
    ## 13   142.0  94.0     No
    ## 14   124.0  88.0    Yes
    ## 16   140.0  90.0     No
    ## 17   138.0  90.0     No
    ## 18   112.0  78.0     No
    ## 19   122.0  84.5     No
    ## 20   139.0  88.0     No
    ## 21   108.0  70.5     No
    ## 23   148.0  78.0     No
    ## 24   132.0  82.0     No
    ## 25   137.5  90.0     No
    ## 26   102.0  68.0     No
    ## 28   132.0  91.0     No
    ## 29   182.0 121.0     No
    ## 30   130.0  88.0     No
    ## 31   102.0  68.0     No
    ## 32   115.0  85.5     No
    ## 33   125.0  90.0     No
    ## 34   150.0  85.0     No
    ## 36   125.0  74.5     No
    ## 37   147.0  74.0     No
    ## 38   124.5  92.5     No
    ## 40   160.0  98.0     No
    ## 41   153.0 101.0     No
    ## 42   111.0  73.0     No
    ## 43   116.5  80.0     No
    ## 44   122.0  78.0     No
    ## 46   132.0  83.5     No
    ## 47   206.0  92.0    Yes
    ## 48    96.0  63.0     No
    ## 49   179.5 114.0     No
    ## 50   119.0  77.5     No
    ## 51   116.0  69.0     No
    ## 53   156.5  92.5     No
    ## 54   112.0  66.0     No
    ## 55   130.0  78.0     No
    ## 56   145.0  82.5     No
    ## 58   116.0  71.0     No
    ## 59   114.0  76.0     No
    ## 60   143.5  81.0     No
    ## 61   115.0  69.0     No
    ## 62   158.0 102.0     No
    ## 63   121.0  79.0     No
    ## 64   157.0  89.0     No
    ## 65   123.5  75.0     No
    ## 66   126.5  80.0     No
    ## 67   136.0  84.0     No
    ## 68   154.0  87.0     No
    ## 69   190.0  99.0     No
    ## 70   107.0  73.0     No
    ## 71   112.5  60.0     No
    ## 72   110.0  67.5     No
    ## 74   116.0  72.5     No
    ## 76   150.0 106.0     No
    ## 77   110.0  76.0     No
    ## 78   138.5  85.0     No
    ## 79   155.0  85.0     No
    ## 81   151.0 101.0     No
    ## 82   152.0  90.0     No
    ## 83   179.0  94.0     No
    ## 84   155.0 110.0     No
    ## 85   138.0  86.5     No
    ## 86   124.0  78.0     No
    ## 87   182.0 101.0     No
    ## 88   113.0  72.5     No
    ## 89   138.0  82.0     No
    ## 90   200.0 104.0     No
    ## 91   124.0  86.0     No
    ## 92   117.5  80.0     No
    ## 93   121.0  61.5     No
    ## 94   138.0  89.0     No
    ## 95   132.5  80.0     No
    ## 96   102.0  71.5     No
    ## 97   154.0  80.0     No
    ## 98   136.0  96.0     No
    ## 99   126.0  79.0     No
    ## 100  123.0  76.5     No
    ## 101  134.0  80.0     No
    ## 103  180.0  90.0     No
    ## 104  121.0  75.0     No
    ## 105  132.5  87.0     No
    ## 106  141.0  84.0     No
    ## 107  110.0  64.0     No
    ## 108  141.0  82.0     No
    ## 109  145.0  77.0     No
    ## 110  100.0  63.0     No
    ## 111  135.0  82.0     No
    ## 112  115.0  79.0     No
    ## 113  138.0  90.0     No
    ## 114  187.0  88.0     No
    ## 115  123.0  69.0     No
    ## 116  139.0  80.0     No
    ## 118  138.0  88.5     No
    ## 119  158.0 105.0     No
    ## 121  127.0  82.0     No
    ## 122  160.5  96.0     No
    ## 123  100.0  68.0     No
    ## 124  112.0  85.5     No
    ## 125  157.0  88.0     No
    ## 126  105.0  75.0     No
    ## 127  112.0  70.0     No
    ## 128  109.0  71.0     No
    ## 129  135.0  97.0     No
    ## 130  158.0  86.0     No
    ## 131  112.5  80.0     No
    ## 132  128.0  77.0     No
    ## 133  118.0  77.0     No
    ## 134  115.0  77.5     No
    ## 135  112.0  73.0     No
    ## 136  106.0  67.5     No
    ## 138  134.0  89.0     No
    ## 139  128.0  64.0     No
    ## 140  121.0  84.0     No
    ## 141  131.0  81.0     No
    ## 142  131.0  81.0     No
    ## 143  124.0  70.0     No
    ## 144  154.0 100.0     No
    ## 145  127.5  81.5     No
    ## 146  138.0  82.0     No
    ## 147  100.0  72.5     No
    ## 148  117.5  72.5     No
    ## 149  151.0  95.0     No
    ## 150  136.0  94.0     No
    ## 151  149.0 100.0     No
    ## 152  141.0  86.0     No
    ## 153  153.0  98.0     No
    ## 154  123.0  81.0     No
    ## 155  180.5 106.5     No
    ## 156  130.0  77.0     No
    ## 157  142.0  93.0     No
    ## 158  136.5  76.0     No
    ## 159  212.0 104.0     No
    ## 161  134.0  82.0     No
    ## 162  153.0  80.5     No
    ## 163  150.0  77.5     No
    ## 164  191.0 124.5    Yes
    ## 165  108.0  81.0     No
    ## 166  121.5  73.0     No
    ## 167  139.0  86.0    Yes
    ## 168  180.0  97.5     No
    ## 169  118.0  86.5     No
    ## 170  102.0  61.0     No
    ## 171  117.5  80.0     No
    ## 172  173.0  89.0     No
    ## 173  109.0  78.0     No
    ## 174  110.0  78.0     No
    ## 175  134.0  83.0     No
    ## 176  109.0  77.0     No
    ## 177  144.0  78.0     No
    ## 178  129.5  83.0     No
    ## 179  122.0  80.5     No
    ## 180  110.0  74.0     No
    ## 181  111.0  64.0     No
    ## 182  124.0  76.0     No
    ## 183  158.0  86.5     No
    ## 184  151.0 101.0     No
    ## 185  182.0 106.0     No
    ## 186  122.0  83.0     No
    ## 187  117.5  71.0     No
    ## 188  110.0  67.0     No
    ## 189  112.0  60.0     No
    ## 190  128.0  82.0     No
    ## 191  118.0  76.0     No
    ## 192  128.0  77.0     No
    ## 194  158.0  98.0     No
    ## 195  117.0  74.5     No
    ## 196  138.0  97.0     No
    ## 197  125.0  86.5     No
    ## 198  144.5  83.5     No
    ## 199  133.0  92.0     No
    ## 200  109.0  70.0     No
    ## 201  127.0  81.0     No
    ## 203  128.0  87.0     No
    ## 204  106.0  71.0     No
    ## 205  117.0  76.0     No
    ## 206  157.0  82.5     No
    ## 207  170.0  86.0     No
    ## 208  137.0  84.0     No
    ## 209  111.0  68.0     No
    ## 210  110.0  67.5     No
    ## 211   94.0  66.5     No
    ## 213  135.0  91.0     No
    ## 215  144.0  78.0     No
    ## 216  155.0  86.0     No
    ## 217  117.5  72.5     No
    ## 218  128.0  68.0     No
    ## 219  125.0  65.0     No
    ## 220  124.0  84.0     No
    ## 222  166.0  98.0     No
    ## 224  126.5  93.0     No
    ## 227  114.0  85.5     No
    ## 228  123.0  72.0     No
    ## 229   96.0  70.0     No
    ## 230  117.5  77.5     No
    ## 231  124.0  66.0     No
    ## 232  130.0  70.0     No
    ## 233  106.0  82.0     No
    ## 234  160.0  85.0     No
    ## 235  134.0  60.0     No
    ## 236  142.0  89.0     No
    ## 237  135.0  92.0     No
    ## 238  174.0  92.5    Yes
    ## 239  177.5  98.0     No
    ## 240  131.0  93.0     No
    ## 241  136.5  99.5     No
    ## 242  128.0  91.0     No
    ## 243  129.0  70.0     No
    ## 244  112.5  79.0     No
    ## 245  159.0  90.0     No
    ## 246  128.0  74.0     No
    ## 247  108.0  74.0     No
    ## 248  155.0  85.0     No
    ## 249  124.0  85.0     No
    ## 250  122.0  84.0     No
    ## 251  129.0  88.0     No
    ## 252  139.0  88.0     No
    ## 253  148.0  90.0     No
    ## 254  134.0  81.0     No
    ## 255  178.0 127.0     No
    ## 256  124.0  76.0     No
    ## 257  130.5  82.0     No
    ## 261  107.5  73.0     No
    ## 262  200.0 122.5     No
    ## 264  144.0  84.0     No
    ## 265  109.0  70.0     No
    ## 266  124.0  74.0     No
    ## 267  125.0  80.0     No
    ## 268  189.0 110.0    Yes
    ## 269  168.0  98.0     No
    ## 270  139.0  90.0     No
    ## 271  108.0  72.5     No
    ## 272  197.5 105.0     No
    ## 273  146.0  89.0     No
    ## 274  137.0  91.0     No
    ## 275  150.0  81.0     No
    ## 277  134.0  97.0     No
    ## 278  107.5  72.5     No
    ## 279  156.5  86.0     No
    ## 280  140.0  93.0     No
    ## 281  134.0  98.0     No
    ## 282  125.0  85.0     No
    ## 283  139.0  98.0     No
    ## 284  174.0  84.5     No
    ## 285  129.0  83.0     No
    ## 287  136.0  87.0     No
    ## 288  127.5  86.0     No
    ## 289  129.0  80.0     No
    ## 290  147.0  92.0     No
    ## 291  127.5  91.0     No
    ## 293  132.0  78.0     No
    ## 294  123.0  79.0     No
    ## 296  100.0  74.5     No
    ## 297  142.0  85.0     No
    ## 298   98.0  64.0     No
    ## 299  105.0  57.0     No
    ## 300  105.0  69.0     No
    ## 301  111.0  68.0     No
    ## 302  126.0  72.0     No
    ## 303  131.5  76.0     No
    ## 304  116.0  72.0     No
    ## 305  125.0  76.0     No
    ## 308  117.5  80.0     No
    ## 310  101.0  63.0     No
    ## 311  130.0  80.0     No
    ## 312  130.0  80.0    Yes
    ## 313  121.0  86.0     No
    ## 316  158.5  81.0     No
    ## 317   97.0  65.0     No
    ## 320  101.0  61.0     No
    ## 321   97.0  64.0     No
    ## 323  114.0  77.0     No
    ## 324  151.5  88.0     No
    ## 325  110.0  65.0     No
    ## 326  140.0  94.0     No
    ## 329  168.0  98.0     No
    ## 330  157.0 101.0     No
    ## 331  125.0  78.0     No
    ## 333  118.0  81.0     No
    ## 334  148.0  85.0     No
    ## 335  131.0  82.0     No
    ## 336  112.0  74.0     No
    ## 337  133.0  78.0     No
    ## 338  127.0  70.0     No
    ## 339   97.5  57.5     No
    ## 340  125.0  79.0     No
    ## 341  189.0 111.0     No
    ## 343  120.0  72.0     No
    ## 344  101.0  67.0     No
    ## 345  110.0  78.5     No
    ## 346  204.0  94.0    Yes
    ## 347  145.0  99.0     No
    ## 348  126.0  87.0     No
    ## 349  126.5  78.5     No
    ## 350  143.5  89.0     No
    ## 352  160.0  82.0     No
    ## 353  158.0 101.0     No
    ## 354  105.0  75.0     No
    ## 355  111.0  60.0     No
    ## 358  136.0  66.0     No
    ## 360  144.0  98.0    Yes
    ## 361  128.0  86.0     No
    ## 362  130.0  86.0     No
    ## 363  110.0  69.0     No
    ## 364  157.5 104.5     No
    ## 365  110.0  67.0     No
    ## 366  106.0  71.0     No
    ## 367  121.0  83.0     No
    ## 368  137.0  89.5     No
    ## 370  140.5  90.0     No
    ## 371  127.5  83.0     No
    ## 372  118.0  79.0     No
    ## 373  171.0 112.0     No
    ## 374  146.0  95.0     No
    ## 375  100.0  70.0     No
    ## 376  118.0  78.0     No
    ## 377   98.0  57.0     No
    ## 378  106.0  66.0     No
    ## 379  130.0  89.0     No
    ## 380  133.0  80.0     No
    ## 381  114.0  67.0     No
    ## 382  128.0  85.0     No
    ## 383  215.0 105.0     No
    ## 384  122.0  85.0     No
    ## 385  159.0 102.0     No
    ## 386  145.0  89.0     No
    ## 387   95.0  55.0     No
    ## 388  156.0  93.0    Yes
    ## 389  155.0  82.5     No
    ## 390  116.0  87.0     No
    ## 391  122.0  74.0     No
    ## 392  109.0  73.0     No
    ## 393  170.0  98.0     No
    ## 394  123.0  89.0     No
    ## 395  125.0  80.0     No
    ## 396  154.0 100.0    Yes
    ## 398  135.0  85.0     No
    ## 399  122.5  66.5     No
    ## 400  120.0  80.0     No
    ## 401  178.0  88.5     No
    ## 402  124.0  87.0     No
    ## 403  127.5  86.5     No
    ## 404  111.0  67.0     No
    ## 405  146.5  83.0     No
    ## 406  117.0  76.0     No
    ## 407  124.0  74.0     No
    ## 408  134.0  86.5     No
    ## 409  149.0  96.0     No
    ## 411  136.5  79.0     No
    ## 412  122.5  77.5     No
    ## 413  132.5  66.0     No
    ## 414  159.0  82.0     No
    ## 416  113.0  68.0     No
    ## 417  113.5  75.0     No
    ## 418  116.0  67.0     No
    ## 419  119.0  75.0     No
    ## 420  138.0  84.0     No
    ## 421  141.0  91.0     No
    ## 423  176.0 119.0     No
    ## 424  137.0  86.0     No
    ## 425  177.5 120.0     No
    ## 426  131.0  86.0     No
    ## 427  156.0 100.0     No
    ## 430  197.0 118.0     No
    ## 431  221.0 118.0    Yes
    ## 432  130.0  80.0     No
    ## 433  111.0  72.0     No
    ## 434  148.0 106.0     No
    ## 436  119.0  82.0     No
    ## 437   90.0  60.0     No
    ## 439  112.5  70.0     No
    ## 440  160.0 120.0     No
    ## 441  110.0  80.0     No
    ## 442  147.0  90.0     No
    ## 444  135.0  74.0     No
    ## 448  109.5  72.0     No
    ## 449   96.0  64.0     No
    ## 450  117.5  80.0     No
    ## 452  120.0  80.0     No
    ## 453  165.0 105.0     No
    ## 455  141.0  89.0     No
    ## 456  124.0  76.0     No
    ## 457  113.0  80.0     No
    ## 458  116.0  81.0     No
    ## 459  140.0  87.0     No
    ## 460  107.0  68.0     No
    ## 461   95.5  59.0     No
    ## 462  146.0  91.0     No
    ## 463  130.0  70.0     No
    ## 464  209.0 133.0     No
    ## 465  112.5  65.0     No
    ## 467  115.0  75.0     No
    ## 468  112.0  76.0     No
    ## 469  132.0  82.0     No
    ## 470  107.5  80.0     No
    ## 471  127.5  83.5     No
    ## 472  136.5  77.5     No
    ## 475  129.0  81.0     No
    ## 476  150.0  85.0    Yes
    ## 477  109.0  71.0     No
    ## 478  124.0  83.0     No
    ## 479  110.0  70.0     No
    ## 480  135.0  89.0     No
    ## 481  146.0 104.0     No
    ## 482  152.0  99.0     No
    ## 483  154.0  95.5     No
    ## 484  125.0  74.0     No
    ## 486  110.0  64.0     No
    ## 487  162.5  92.5     No
    ## 488  134.0  86.0     No
    ## 489  177.5  95.0    Yes
    ## 491  107.0  66.0     No
    ## 492  150.0  96.5     No
    ## 493  116.0  74.5     No
    ## 494  127.0  74.0     No
    ## 495  125.0  80.0     No
    ## 496  138.0  90.0     No
    ## 497  145.0 100.0     No
    ## 498  103.0  70.0     No
    ## 499  145.0  85.0     No
    ## 500  110.0  75.0     No
    ## 501  295.0 135.0     No
    ## 502  158.0  90.0     No
    ## 503  125.0  80.0     No
    ## 504  126.0  81.0     No
    ## 505  145.0  90.0     No
    ## 506  120.0  80.0     No
    ## 508  150.0 101.0     No
    ## 509   98.0  64.5     No
    ## 511  130.0  72.0     No
    ## 512  134.0  78.0     No
    ## 513  135.0  93.0     No
    ## 514  134.0  89.5     No
    ## 515  114.0  81.0     No
    ## 516  137.0  76.0     No
    ## 517  151.0  95.0     No
    ## 518  121.0  76.0     No
    ## 519  111.0  68.5     No
    ## 521  103.0  71.0     No
    ## 522  155.0  93.0     No
    ## 523  132.0  82.0     No
    ## 525  144.0  88.0     No
    ## 526  108.0  71.0     No
    ## 527  160.5  98.5     No
    ## 528  139.0  81.0     No
    ## 529  134.5  89.0     No
    ## 530  115.5  62.0     No
    ## 531  132.0  84.0     No
    ## 533  155.0  84.0     No
    ## 534  137.5  89.5     No
    ## 535  133.0  84.0     No
    ## 536  135.0  95.0     No
    ## 537  130.0  78.0     No
    ## 538  155.0 100.0     No
    ## 539  122.0  72.5     No
    ## 541  117.5  65.0     No
    ## 542  129.0  80.0     No
    ## 543  148.0  96.0     No
    ## 544  132.0  84.5     No
    ## 545  173.0 117.0     No
    ## 546  119.0  80.0     No
    ## 547  130.0  82.0     No
    ## 548  122.5  76.0     No
    ## 549  105.0  74.0     No
    ## 550  120.0  73.0     No
    ## 551  125.0  83.0     No
    ## 552  189.0 121.0     No
    ## 553  113.0  75.0     No
    ## 554  123.0  83.0     No
    ## 555  145.0 102.0     No
    ## 556  105.0  59.5     No
    ## 557  130.0  85.0     No
    ## 558  171.0 101.0     No
    ## 559  114.0  79.0     No
    ## 560  125.0  83.0     No
    ## 562  131.0  88.0    Yes
    ## 563  115.0  72.5     No
    ## 564  145.0  74.0     No
    ## 565  110.0  71.0     No
    ## 566  132.0  92.0     No
    ## 567  117.0  65.0     No
    ## 568  112.0  83.0     No
    ## 569  121.5  78.0     No
    ## 570  144.0  96.5     No
    ## 571  114.0  80.0     No
    ## 572  128.0  94.0     No
    ## 574  142.5  90.0     No
    ## 575  174.5 103.0     No
    ## 576  140.0  82.0     No
    ## 577  110.0  75.0     No
    ## 578  129.0  76.0     No
    ## 579  130.0  85.0     No
    ## 580  136.5  83.5     No
    ## 581  120.0  79.0     No
    ## 582  122.5  75.0     No
    ## 583  113.5  75.5     No
    ## 584  128.0  93.0     No
    ## 586  164.0 102.5     No
    ## 587  116.0  72.0     No
    ## 588  125.0  80.0     No
    ## 589  165.0  77.0     No
    ## 591  108.0  66.0     No
    ## 592  135.0  89.0     No
    ## 593  121.0  82.0     No
    ## 594  133.0  91.0     No
    ## 595  120.0  74.0     No
    ## 596  150.0  96.5     No
    ## 597  148.0 100.0     No
    ## 598  122.0  76.0     No
    ## 601  121.5  82.5     No
    ## 602  163.0  78.0     No
    ## 603  139.0  92.0     No
    ## 604  117.0  80.0     No
    ## 605  146.0  94.0     No
    ## 606  125.0  74.0     No
    ## 607  118.5  85.5     No
    ## 608  142.0  87.0     No
    ## 609   98.0  65.0     No
    ## 610  116.0  73.0     No
    ## 611  185.0 114.0     No
    ## 612  126.0  82.0     No
    ## 613  220.0 118.0    Yes
    ## 615  115.0  72.5     No
    ## 616  123.0  88.0     No
    ## 617  125.0  73.5     No
    ## 618  119.0  82.0     No
    ## 619  121.5  74.0     No
    ## 620  144.0  78.0     No
    ## 621  120.0  70.0     No
    ## 622  123.0  90.0     No
    ## 623  113.0  75.0     No
    ## 624  103.0  73.0     No
    ## 625  132.0  91.0     No
    ## 626  140.0  73.0     No
    ## 627  164.0  88.0     No
    ## 628  153.0 108.0     No
    ## 629  145.0  79.0     No
    ## 630  120.5  78.0     No
    ## 631  145.0  91.0     No
    ## 632  116.0  86.0     No
    ## 634  155.0  92.5     No
    ## 635   98.5  69.5     No
    ## 636  107.0  61.0     No
    ## 637  108.0  70.0     No
    ## 638  127.5  72.5     No
    ## 639  120.5  75.0     No
    ## 640  117.5  77.5     No
    ## 641  161.0  96.0     No
    ## 642  148.0  98.0     No
    ## 643  139.5  89.0     No
    ## 644  127.5  80.0     No
    ## 645  127.5  87.5     No
    ## 647  105.0  68.0     No
    ## 648  180.0 101.0     No
    ## 649  123.0  92.5     No
    ## 650  115.0  71.0     No
    ## 651   96.0  67.0     No
    ## 652  130.5  85.5     No
    ## 653  117.0  75.0     No
    ## 654  168.5 108.0     No
    ## 655  103.5  70.0     No
    ## 656  129.5  93.0     No
    ## 657  148.0 108.0     No
    ## 658  134.5  78.5     No
    ## 659  117.5  67.5     No
    ## 660  176.0  99.0     No
    ## 661  111.0  73.0     No
    ## 662  109.0  71.0     No
    ## 663  125.0  80.0     No
    ## 664  141.0  93.0     No
    ## 665  131.0  85.5     No
    ## 666  119.0  77.0     No
    ## 667  123.5  81.0     No
    ## 668  163.5 103.0     No
    ## 669  117.0  90.0     No
    ## 671  127.5  81.0     No
    ## 672  128.5  83.0     No
    ## 673  128.5  74.5     No
    ## 674  134.0  84.0     No
    ## 675  134.0  81.0     No
    ## 676  111.0  70.0     No
    ## 677  128.5  82.0     No
    ## 678  148.0  91.0     No
    ## 679  107.0  68.5     No
    ## 680  129.0  73.0     No
    ## 681  130.0  86.5     No
    ## 682  121.0  75.0     No
    ## 683  112.0  76.0     No
    ## 684  127.0  78.0     No
    ## 685  126.5  92.5    Yes
    ## 687  138.0  91.0     No
    ## 688  167.0  89.5     No
    ## 689  205.5 104.5     No
    ## 690  150.0  92.0    Yes
    ## 691  136.0  92.0     No
    ## 692  118.5  73.0     No
    ## 693  112.0  68.5     No
    ## 694  145.0  81.0     No
    ## 695  138.0  79.0     No
    ## 696  115.0  70.0     No
    ## 697  110.0  78.0     No
    ## 699  119.5  85.5     No
    ## 701  119.0  80.0     No
    ## 702  167.5  95.0     No
    ## 703  130.0  78.0     No
    ## 705  136.5  85.5     No
    ## 706  148.0  92.5     No
    ## 707  166.0 102.0     No
    ## 708  144.0  80.0     No
    ## 709  137.0  79.0     No
    ## 710  152.5  82.0     No
    ## 711  127.0  81.0     No
    ## 712  119.0  83.0     No
    ## 713  125.0  74.0     No
    ## 714  172.5  75.0     No
    ## 716  126.0  84.0     No
    ## 717  176.0  98.0     No
    ## 718  145.0  95.0     No
    ## 719  146.0  93.5     No
    ## 721  132.5  91.0     No
    ## 722  105.0  70.0     No
    ## 723  117.0  68.0     No
    ## 724  186.0 102.0     No
    ## 725  183.0  93.0    Yes
    ## 726  118.0  73.0     No
    ## 727  110.0  67.5     No
    ## 728  136.0  86.0     No
    ## 729  109.0  72.5     No
    ## 730  135.0  85.5     No
    ## 732  140.5  84.0     No
    ## 733  114.0  79.0     No
    ## 734  148.0  84.0     No
    ## 735  120.0  80.0     No
    ## 736  134.0  80.0     No
    ## 737  134.0  90.5     No
    ## 738  107.5  67.5     No
    ## 739  126.0  78.0     No
    ## 741  127.5  80.0     No
    ## 742  141.0  99.0     No
    ## 743  124.0  84.0     No
    ## 744  141.0  79.0     No
    ## 746  120.0  80.0     No
    ## 747  153.5 105.0     No
    ## 748  134.0  79.0     No
    ## 749  152.0  99.0    Yes
    ## 750  163.5  93.0     No
    ## 751  135.0  87.0     No
    ## 752  130.0  77.5     No
    ## 753  103.0  68.0     No
    ## 754  117.0  74.0     No
    ## 755  107.0  73.0     No
    ## 756  122.0  78.0     No
    ## 757  155.0  92.5     No
    ## 758  117.0  73.0     No
    ## 759  120.0  75.0     No
    ## 760  131.0  74.0     No
    ## 761  144.0  84.0     No
    ## 762  147.5  92.5     No
    ## 763  212.0 103.0     No
    ## 765  183.0 114.5    Yes
    ## 766  120.0  74.0     No
    ## 767  178.0  91.0     No
    ## 768  128.0  76.0     No
    ## 769  126.0  80.0     No
    ## 770  132.0  78.0     No
    ## 771  127.0  76.5     No
    ## 772  143.5  92.0     No
    ## 773  175.0  94.0     No
    ## 774  160.0 108.0     No
    ## 775  119.0  82.0     No
    ## 776  142.5  90.0     No
    ## 777  100.0  73.0     No
    ## 778  192.0 102.0    Yes
    ## 779  127.0  75.0     No
    ## 780  110.0  62.5     No
    ## 781  129.0  85.0     No
    ## 782  117.5  75.0     No
    ## 784  134.5  83.0     No
    ## 786  131.0  87.0     No
    ## 787  160.0  92.0     No
    ## 789   96.5  71.0     No
    ## 790  107.0  66.5     No
    ## 791  157.0  94.0     No
    ## 792  125.0  89.0     No
    ## 794  180.0 108.0     No
    ## 795  185.0 100.0     No
    ## 796  147.0  94.0     No
    ## 798  123.0  69.0     No
    ## 799  148.0  85.5     No
    ## 800  141.0 102.0     No
    ## 801  133.0  83.0     No
    ## 802  118.0  84.0     No
    ## 803  160.0 100.0     No
    ## 804  114.0  66.5     No
    ## 808  118.0  82.0     No
    ## 809  112.0  62.5     No
    ## 811  150.0 101.0     No
    ## 812  137.0  81.0     No
    ## 813  118.5  82.5     No
    ## 814  126.0  82.0    Yes
    ## 815  121.5  74.0     No
    ## 816  133.0  76.0     No
    ## 817  140.5  94.0     No
    ## 818  159.5  91.0     No
    ## 819  118.5  69.0     No
    ## 820  164.0 107.0     No
    ## 821  118.5  82.0     No
    ## 822  129.0  90.0     No
    ## 823  130.0  92.0     No
    ## 824  177.0 110.0    Yes
    ## 825  156.0 104.0     No
    ## 827  116.0  78.5     No
    ## 828  152.0  82.0    Yes
    ## 829  117.0  80.0     No
    ## 830  127.0  79.0     No
    ## 831  150.0  90.0     No
    ## 832  162.0  87.5     No
    ## 833  124.0  86.5     No
    ## 834  132.0  91.0     No
    ## 835  130.0  82.0     No
    ## 836  120.0  70.0     No
    ## 837  142.0  93.0     No
    ## 838  125.0  72.5     No
    ## 839  164.0  94.5     No
    ## 840  112.5  82.5     No
    ## 841  108.0  66.0     No
    ## 844  142.0  82.0     No
    ## 845  135.0  79.0     No
    ## 846  111.0  70.0     No
    ## 847  128.0  82.0     No
    ## 848  124.0  79.0     No
    ## 850  129.0  86.0     No
    ## 851  168.0  82.0     No
    ## 853  102.5  65.0     No
    ## 854  110.0  78.5     No
    ## 855  134.0  87.0     No
    ## 856  127.5  85.0     No
    ## 857  122.0  70.0     No
    ## 858  126.0  85.0     No
    ## 860  105.0  84.0     No
    ## 861  147.0  90.0     No
    ## 863  200.0 140.0     No
    ## 865  127.0  75.0     No
    ## 866  128.0  91.0     No
    ## 867  137.0  96.5     No
    ## 868  138.0  83.0     No
    ## 869  133.0  80.0     No
    ## 870  120.0  73.0     No
    ## 871  122.0  82.0     No
    ## 872  121.5  78.0     No
    ## 873  116.0  79.0     No
    ## 874  132.0  84.0     No
    ## 875  131.5  78.5     No
    ## 876  118.0  66.0     No
    ## 877  118.5  80.5     No
    ## 878  114.0  70.5     No
    ## 879  102.0  69.0     No
    ## 880  131.0  82.0     No
    ## 881  119.0  85.0     No
    ## 882  107.0  76.5     No
    ## 883  122.5  67.5     No
    ## 884  131.0  80.0     No
    ## 885  114.0  83.0     No
    ## 886  111.0  62.0     No
    ## 887  118.0  86.0     No
    ## 888  146.0  88.0     No
    ## 889  141.0  81.0     No
    ## 890  107.0  73.0     No
    ## 891  128.0  83.0     No
    ## 892  163.0  85.0     No
    ## 893  136.0  99.0     No
    ## 894  244.0 124.0    Yes
    ## 895  134.0  92.0     No
    ## 896  113.0  77.0     No
    ## 897  113.5  74.0     No
    ## 898  104.0  76.0     No
    ## 899  135.0  82.0     No
    ## 900  114.0  81.0     No
    ## 901  131.0  79.0     No
    ## 904  127.5  75.0     No
    ## 905  122.0  73.0     No
    ## 906  128.0  86.5     No
    ## 907  127.5  83.5     No
    ## 908  153.0  75.0     No
    ## 909  168.0 102.0     No
    ## 910  127.5  81.0     No
    ## 911  134.0  91.0     No
    ## 912  113.5  73.0     No
    ## 913  156.0  98.0     No
    ## 915  112.0  69.0     No
    ## 916  164.0  85.5     No
    ## 917  176.0  97.0     No
    ## 918  121.0  80.0     No
    ## 919  141.0  83.0     No
    ## 920  148.0  93.0     No
    ## 921  123.0  73.5     No
    ## 922  152.5  90.0    Yes
    ## 923  108.0  72.0     No
    ## 924  111.0  61.0     No
    ## 925  126.0  85.0     No
    ## 926  141.0  92.0     No
    ## 927  130.0  82.5     No
    ## 929  108.0  77.0     No
    ## 930  168.0 102.0     No
    ## 931  146.0  89.0     No
    ## 932  128.0  77.0     No
    ## 933  134.5  87.0     No
    ## 934  124.5  82.0     No
    ## 935  213.0  94.5     No
    ## 936  115.0  81.0     No
    ## 938  127.0  76.0     No
    ## 939  129.0  92.0     No
    ## 940  126.0  82.0     No
    ## 941  167.5  91.5     No
    ## 942  108.0  66.0     No
    ## 943  125.0  75.0     No
    ## 944  153.0  94.0     No
    ## 945  117.0  71.0     No
    ## 946  152.0 103.0     No
    ## 947  118.0  74.0     No
    ## 948  154.0  88.5     No
    ## 950  114.0  84.0     No
    ## 952  120.0  90.0     No
    ## 953  111.0  67.5     No
    ## 954  116.0  82.0     No
    ## 955  115.0  80.0     No
    ## 957  159.0 115.0     No
    ## 958  134.0  84.0     No
    ## 959  206.0  98.0     No
    ## 960  123.0  71.0     No
    ## 961  157.0  99.0     No
    ## 962  142.0  74.0     No
    ## 963  106.0  70.5     No
    ## 964  139.0  88.0     No
    ## 965  242.0 141.0    Yes
    ## 966  111.0  70.0     No
    ## 967  151.0  96.0     No
    ## 968  199.0  83.0     No
    ## 969  184.0 109.0     No
    ## 970  129.5  87.0     No
    ## 971  127.5  80.0     No
    ## 973  171.5  83.0     No
    ## 974  129.0  90.0     No
    ## 975  133.0  86.0     No
    ## 977  167.5 102.5     No
    ## 978  160.0  85.0     No
    ## 981  127.5  75.0     No
    ## 982  115.5  68.0     No
    ## 983  128.0  85.0     No
    ## 984  118.0  77.0     No
    ## 985  198.0 107.0     No
    ## 986  117.5  67.5     No
    ## 987  127.0  86.0     No
    ## 989  114.5  62.5     No
    ## 990  137.5  77.5     No
    ## 991  151.0  98.0     No
    ## 992  172.5  88.0     No
    ## 993  151.0  89.5    Yes
    ## 994  158.5 107.0     No
    ## 995  118.0  84.0     No
    ## 996  160.0  97.0     No
    ## 997  137.0  82.5     No
    ## 998  131.5  92.0     No
    ## 999  123.5  84.0     No
    ## 1000 140.0  82.0     No
    ## 1001  98.0  66.0     No
    ## 1002 109.5  67.0     No
    ## 1003 115.5  65.5     No
    ## 1004 106.0  79.0     No
    ## 1006 122.0  81.0     No
    ## 1007 137.0  88.0     No
    ## 1008 160.0 109.0     No
    ## 1009 134.0  98.0     No
    ## 1010 125.5  86.0     No
    ## 1012 120.5  67.5     No
    ## 1013 131.0  74.0     No
    ## 1014 159.0 109.0     No
    ## 1015 170.0 100.0    Yes
    ## 1016 124.0  90.0     No
    ## 1017 154.0 100.0     No
    ## 1018 132.0  75.0     No
    ## 1020 111.0  72.0     No
    ## 1021 110.0  77.0     No
    ## 1022 149.0  81.0     No
    ## 1023 135.0  82.5     No
    ## 1024 127.0  89.0     No
    ## 1025 121.0  67.0     No
    ## 1026 115.0  82.0     No
    ## 1027 139.0  88.0     No
    ## 1028 122.0  80.0     No
    ## 1029 113.0  80.0     No
    ## 1030 124.0  86.0     No
    ## 1031 139.0  92.0     No
    ## 1032 118.0  77.5     No
    ## 1033 122.0  76.0     No
    ## 1035 107.5  71.0     No
    ## 1036 125.0  85.0     No
    ## 1039 102.0  64.0     No
    ## 1040 124.0  85.0     No
    ## 1042 206.0 104.0     No
    ## 1043 117.0  67.0     No
    ## 1044 145.0  94.0     No
    ## 1045 107.0  73.5     No
    ## 1046 142.0  94.0     No
    ## 1047 167.0  96.5     No
    ## 1048 142.0  82.5     No
    ## 1049 111.5  70.0     No
    ## 1050 148.0  81.0     No
    ## 1051 133.0  88.0     No
    ## 1052 100.0  71.0     No
    ## 1053 132.0  86.0     No
    ## 1054 139.0  84.0     No
    ## 1055 119.5  74.0     No
    ## 1056 127.0  83.0     No
    ## 1057 147.0  71.5     No
    ## 1059 136.5  95.0     No
    ## 1060 149.0  90.0     No
    ## 1061 155.0  74.0     No
    ## 1062 133.0  84.0     No
    ## 1063 131.5  84.0     No
    ## 1064 117.5  77.5     No
    ## 1065 131.0  88.0     No
    ## 1068 105.5  67.0     No
    ## 1069 138.0  86.0     No
    ## 1070 126.0  84.0     No
    ## 1071 145.0  90.0     No
    ## 1072 122.0  83.0     No
    ## 1073 141.5  91.0     No
    ## 1074 111.0  71.0     No
    ## 1075 142.0  90.0     No
    ## 1076 131.0  88.0     No
    ## 1078 117.0  80.0     No
    ## 1079 143.0  98.0     No
    ## 1080 165.0 115.0    Yes
    ## 1081 161.5  95.0     No
    ## 1082 122.0  66.5     No
    ## 1084 164.5  93.5     No
    ## 1085 120.0  80.0     No
    ## 1086 141.5  98.0     No
    ## 1087 124.0  82.0     No
    ## 1088 112.0  77.0     No
    ## 1089 137.0  97.0     No
    ## 1090 163.0 105.0     No
    ## 1091 123.0  76.5     No
    ## 1092 116.0  62.0     No
    ## 1093 171.5 105.5     No
    ## 1094 150.0  97.0     No
    ## 1095 126.0  68.0     No
    ## 1096 132.5  90.0     No
    ## 1097 146.0  88.0     No
    ## 1098 177.0 103.5     No
    ## 1099 105.0  74.0     No
    ## 1100 108.5  63.5     No
    ## 1101 110.0  60.0     No
    ## 1102 102.0  72.0     No
    ## 1105 129.5  92.5     No
    ## 1106 152.5  91.0     No
    ## 1107 132.0  85.0     No
    ## 1108 131.0  80.0     No
    ## 1109 146.5  82.0     No
    ## 1110 114.0  77.0     No
    ## 1112 124.0  70.0     No
    ## 1113 160.0  93.0    Yes
    ## 1115 153.0  93.0     No
    ## 1116 109.5  69.0     No
    ## 1118 132.0  67.0     No
    ## 1119 131.0  79.0     No
    ## 1120 128.0  80.5     No
    ## 1121 125.0  83.0    Yes
    ## 1123 201.0  93.0    Yes
    ## 1124 120.5  84.0     No
    ## 1125 145.0  90.0     No
    ## 1126 144.0  86.0     No
    ## 1127 127.0  70.0     No
    ## 1128 122.5  82.5     No
    ## 1129 131.0  81.0     No
    ## 1130 118.0  79.0     No
    ## 1131 112.5  76.5     No
    ## 1132 135.0  91.0     No
    ## 1134 118.0  78.0     No
    ## 1135 123.0  75.0     No
    ## 1136 112.0  66.0     No
    ## 1137 137.0  87.0     No
    ## 1138 114.0  81.0     No
    ## 1139 137.0  76.0     No
    ## 1140 125.0  71.0     No
    ## 1141 113.5  66.5     No
    ## 1142 113.0  76.5     No
    ## 1143 114.0  78.0     No
    ## 1144 110.0  72.5     No
    ## 1145 125.0  87.0     No
    ## 1147 131.0  87.0     No
    ## 1148 149.0  73.0     No
    ## 1151 137.0  99.0     No
    ## 1152 213.0 103.0    Yes
    ## 1153 105.0  65.0     No
    ## 1154 130.0  77.5     No
    ## 1155 145.0  90.0     No
    ## 1156 148.5  90.0    Yes
    ## 1157 159.5  94.0     No
    ## 1158 143.5  89.0     No
    ## 1159 140.0  80.0     No
    ## 1160 119.0  75.0     No
    ## 1162 105.5  64.5     No
    ## 1163 127.5  79.5     No
    ## 1166 114.0  71.0     No
    ## 1167 101.0  72.0     No
    ## 1169 133.0  72.5     No
    ## 1170 138.0  92.0     No
    ## 1171 173.0 102.0     No
    ## 1172 112.5  64.5     No
    ## 1173 109.0  74.0     No
    ## 1174 189.0 103.0     No
    ## 1176 141.5  88.5     No
    ## 1177 159.0  64.0     No
    ## 1178 146.0  77.0     No
    ## 1179 140.0  84.0     No
    ## 1180 122.0  87.5     No
    ## 1181 130.0  80.0     No
    ## 1182 125.0  80.5     No
    ## 1183 132.0  79.0     No
    ## 1184 159.0  89.0     No
    ## 1186 144.0  99.0     No
    ## 1187 133.0  88.0     No
    ## 1188 147.0  96.0     No
    ## 1190 141.0  73.5     No
    ## 1191 156.0  90.0     No
    ## 1192 174.0 112.0     No
    ## 1193 150.0  84.0     No
    ## 1194 140.0  88.0     No
    ## 1195 129.0  82.0     No
    ## 1196 130.0  84.0     No
    ## 1197 100.0  60.0     No
    ## 1198 113.0  68.0     No
    ## 1199 128.0  84.5     No
    ## 1200 145.0  77.0     No
    ## 1201 172.0  84.0     No
    ## 1202 167.0 107.0     No
    ## 1204 116.0  69.0     No
    ## 1205 122.0  73.0     No
    ## 1206 107.5  65.0     No
    ## 1207 144.0  91.0     No
    ## 1208 188.0 107.0    Yes
    ## 1209 119.0  73.0     No
    ## 1211 100.0  65.0     No
    ## 1212 130.0  80.0     No
    ## 1213 116.0  79.0     No
    ## 1214 147.0  86.0     No
    ## 1215 140.0  90.0     No
    ## 1216 114.0  81.0     No
    ## 1217 153.5 103.5     No
    ## 1218 127.5  67.5     No
    ## 1219 126.0  86.0     No
    ## 1220 138.0  76.0     No
    ## 1221 175.0 107.5     No
    ## 1223 123.0  69.0     No
    ## 1224 165.0  91.0     No
    ## 1226 131.0  83.0     No
    ## 1227 135.0  94.5     No
    ## 1229 131.0  84.0     No
    ## 1231 145.0  87.5     No
    ## 1233 130.0  94.0     No
    ## 1234 184.0 106.0     No
    ## 1235 135.0  80.0     No
    ## 1236 118.0  76.5     No
    ## 1237 160.5 109.0     No
    ## 1238 114.0  64.0     No
    ## 1239 177.0  96.0     No
    ## 1240 243.0 142.5     No
    ## 1242 178.0 103.0     No
    ## 1243 139.0  80.0     No
    ## 1244 127.0  77.0     No
    ## 1245 184.0  90.0     No
    ## 1246 145.5  87.5     No
    ## 1247 133.0  87.0     No
    ## 1248 187.5  85.5    Yes
    ## 1249 119.0  76.0     No
    ## 1250 125.0  86.0     No
    ## 1251 165.0  96.0     No
    ## 1252 138.0  92.0     No
    ## 1253 126.0  84.0     No
    ## 1255 178.0 106.0     No
    ## 1256 138.0 105.0     No
    ## 1257 123.0  85.5     No
    ## 1259 110.0  69.0     No
    ## 1260 122.0  82.0     No
    ## 1261 119.0  76.0     No
    ## 1262 118.0  72.0     No
    ## 1263 110.0  78.0     No
    ## 1265 155.0  85.0     No
    ## 1266 111.0  81.0     No
    ## 1268 150.0  88.0     No
    ## 1269 116.0  75.0     No
    ## 1270 117.0  76.0     No
    ## 1271 173.0  84.0     No
    ## 1272  99.0  67.0     No
    ## 1273 115.0  77.0     No
    ## 1274 181.0  90.0     No
    ## 1275 115.0  75.0     No
    ## 1276 174.0 110.0     No
    ## 1277  96.0  61.0     No
    ## 1278 122.0  72.0     No
    ## 1279 120.5  85.5     No
    ## 1280 127.0  82.0     No
    ## 1281 102.5  67.5     No
    ## 1283 130.0  80.0     No
    ## 1284 124.0  81.0     No
    ## 1285 180.0 109.5     No
    ## 1286 112.0  74.5     No
    ## 1287 133.5  76.0     No
    ## 1288 100.5  66.0     No
    ## 1290 120.0  86.0     No
    ## 1291 124.0  70.0     No
    ## 1292 156.0  92.0     No
    ## 1293 156.5 105.0     No
    ## 1295 140.0  88.0     No
    ## 1296 128.0  69.0     No
    ## 1297 119.0  84.0     No
    ## 1298 130.0  68.0     No
    ## 1299 118.0  78.0     No
    ## 1300 156.0  90.0     No
    ## 1301 148.0 100.0     No
    ## 1302 113.0  77.0     No
    ## 1303 159.0  88.0     No
    ## 1304 127.0  93.0     No
    ## 1306 143.0 101.0     No
    ## 1307 125.0  92.0     No
    ## 1308 144.0  86.0     No
    ## 1310 122.0  70.0     No
    ## 1311 129.0  89.0     No
    ## 1312 108.5  71.5     No
    ## 1313 162.0  99.0     No
    ## 1314 135.5  80.0     No
    ## 1316 144.0  85.0     No
    ## 1317 123.0  71.0     No
    ## 1318 128.0  76.0    Yes
    ## 1319 126.5  67.5     No
    ## 1320 131.5  89.0     No
    ## 1321 120.0  74.0     No
    ## 1322 143.0 104.0     No
    ## 1323 128.0  74.0     No
    ## 1324 160.0  96.0     No
    ## 1325 122.0  74.0     No
    ## 1326  97.0  63.0     No
    ## 1327 146.0  95.0     No
    ## 1328 155.0 102.0     No
    ## 1329 128.0  90.0     No
    ## 1330 114.0  85.0     No
    ## 1332 175.0  95.0     No
    ## 1333 125.0  94.0     No
    ## 1334 121.0  85.5     No
    ## 1335 116.0  77.0     No
    ## 1336 117.5  80.0     No
    ## 1337 108.0  64.0     No
    ## 1338 132.5  87.5     No
    ## 1339 108.0  70.0     No
    ## 1341 108.5  75.0     No
    ## 1343 110.0  70.0     No
    ## 1345 140.0  85.0     No
    ## 1346 110.0  73.0     No
    ## 1347 103.0  61.0     No
    ## 1348 122.5  80.0     No
    ## 1349 115.0  81.0     No
    ## 1350 163.0 105.0     No
    ## 1351  95.0  58.0     No
    ## 1352 134.0  90.0     No
    ## 1354 173.0  59.0     No
    ## 1355 148.0  99.0     No
    ## 1356 131.0  87.0     No
    ## 1358 172.5  98.0     No
    ## 1359 174.0 100.0     No
    ## 1360 125.0  89.0     No
    ## 1361 126.0  82.0     No
    ## 1362 120.0  72.5     No
    ## 1363 137.0  75.0     No
    ## 1364 164.0  81.0     No
    ## 1365 121.0  74.0     No
    ## 1366 132.0  87.0     No
    ## 1368 115.0  70.0     No
    ## 1369 128.0  87.0     No
    ## 1370 108.0  67.0     No
    ## 1372 114.0  70.0     No
    ## 1373 136.0  70.0     No
    ## 1375 127.5  85.0     No
    ## 1376 106.0  78.0     No
    ## 1377 127.0  81.0     No
    ## 1378 126.0  70.0     No
    ## 1379 158.0  97.0     No
    ## 1380 154.0  98.0     No
    ## 1381 130.0  86.0     No
    ## 1382 124.0  78.0     No
    ## 1383 112.5  73.5     No
    ## 1384 135.0  90.0     No
    ## 1386 140.0 108.0     No
    ## 1387 133.5  87.5     No
    ## 1388 122.0  76.0     No
    ## 1389 120.0  84.0     No
    ## 1390 113.0  80.0     No
    ## 1391 114.0  72.0     No
    ## 1392 128.0  94.0     No
    ## 1394 151.0  94.0     No
    ## 1395 124.0  83.0     No
    ## 1396 115.0  82.0     No
    ## 1397 120.0  77.0     No
    ## 1398 102.0  66.5     No
    ## 1399 121.0  73.0     No
    ## 1400 127.5  80.0     No
    ## 1401 111.0  79.0     No
    ## 1402 112.5  70.0     No
    ## 1403 103.5  66.5     No
    ## 1404 142.5  97.5     No
    ## 1406 114.0  80.0     No
    ## 1407 110.0  73.0     No
    ## 1408 114.0  73.0     No
    ## 1409 123.0  80.0     No
    ## 1410 114.0  84.0     No
    ## 1411 124.0  83.0     No
    ## 1412 141.5  95.0     No
    ## 1414 119.0  73.0     No
    ## 1415 115.0  83.0     No
    ## 1416 163.0  86.0     No
    ## 1417 120.0  68.0     No
    ## 1419 145.0  89.0     No
    ## 1421 137.0  81.0     No
    ## 1422 127.0  76.0     No
    ## 1423 144.0  88.0     No
    ## 1424 138.0  87.0     No
    ## 1425 149.5  78.0     No
    ## 1426 146.0  89.0     No
    ## 1428 101.0  70.0     No
    ## 1429 133.0  83.0     No
    ## 1430 140.0  83.0     No
    ## 1431 128.0  77.0     No
    ## 1432 106.0  67.5     No
    ## 1433 113.5  74.0     No
    ## 1435 138.0  88.0     No
    ## 1437 135.0  76.0     No
    ## 1438 138.0  86.0     No
    ## 1440 182.5  88.0     No
    ## 1441 109.0  73.0     No
    ## 1442 132.0  86.0     No
    ## 1443 172.5 116.5     No
    ## 1444 121.0  79.0     No
    ## 1445 126.0  73.0     No
    ## 1446 112.0  82.0     No
    ## 1447 121.0  81.0     No
    ## 1448 134.0  87.0    Yes
    ## 1449 154.0  90.0     No
    ## 1450 125.0  85.0     No
    ## 1451 109.0  72.0     No
    ## 1452 117.0  79.0     No
    ## 1453 105.0  70.0     No
    ## 1456 115.0  76.0     No
    ## 1457 149.0 100.0     No
    ## 1458 123.0  78.0     No
    ## 1459 149.0  90.0     No
    ## 1461 122.5  82.5     No
    ## 1462 162.0  98.0     No
    ## 1463 129.0  86.0     No
    ## 1464 144.0  90.0     No
    ## 1465 137.0  85.0     No
    ## 1466 122.0  94.0     No
    ## 1467 113.0  79.0     No
    ## 1469 175.0  67.5     No
    ## 1470 143.0  87.0     No
    ## 1471 133.0  88.0     No
    ## 1473 100.0  64.0     No
    ## 1474 115.0  77.0     No
    ## 1475 142.0  85.5     No
    ## 1476 142.0  88.0     No
    ## 1477 132.0  70.0     No
    ## 1478 139.0  86.0     No
    ## 1479 131.0  85.0     No
    ## 1480 146.0  86.0     No
    ## 1481 123.0  87.0     No
    ## 1482 118.0  78.0     No
    ## 1483  95.0  57.0     No
    ## 1484 158.5 100.5     No
    ## 1485 142.0  92.0     No
    ## 1486  95.0  65.0     No
    ## 1487 115.0  74.0     No
    ## 1488 115.0  76.0     No
    ## 1489 131.0  89.0     No
    ## 1490 136.0  90.0     No
    ## 1492 123.0  72.0     No
    ## 1493 130.0  72.0     No
    ## 1494 127.5  77.5     No
    ## 1495 130.0  82.0     No
    ## 1496 121.0  82.0     No
    ## 1497 135.0  80.0     No
    ## 1498 145.0  82.0     No
    ## 1499 125.0  80.0     No
    ## 1500 136.5  83.0     No
    ## 1501 115.0  77.5     No
    ## 1502 117.0  81.0     No
    ## 1503 116.0  66.0     No
    ## 1504 165.0  80.0     No
    ## 1505 141.0  98.0     No
    ## 1506 120.0  74.0     No
    ## 1507 110.0  70.0     No
    ## 1508 162.5 105.0     No
    ## 1509 114.0  64.0     No
    ## 1511 127.0  80.0     No
    ## 1513 115.0  69.0     No
    ## 1514 110.0  72.0     No
    ## 1515 112.0  65.0     No
    ## 1516 118.0  82.0     No
    ## 1517 151.0  97.5     No
    ## 1518 109.0  81.0     No
    ## 1519 119.0  76.0     No
    ## 1520 135.0  70.0     No
    ## 1521 135.0  89.0     No
    ## 1522 100.0  60.0     No
    ## 1523 117.5  70.5     No
    ## 1524 139.0  93.0     No
    ## 1525 132.0  75.5     No
    ## 1526 125.5  94.0     No
    ## 1527 121.0  81.0     No
    ## 1528 160.0  95.0     No
    ## 1529 111.0  80.0     No
    ## 1530 127.0  79.0     No
    ## 1531 134.0  95.0     No
    ## 1532 148.0  93.0     No
    ## 1534 107.5  70.0     No
    ## 1535 125.0  82.0     No
    ## 1536 124.0  90.0     No
    ## 1537 156.0  96.0     No
    ## 1538 128.0  71.0     No
    ## 1539 132.0  94.0     No
    ## 1540 137.5  92.0     No
    ## 1542 137.0  92.5     No
    ## 1543 115.0  70.0     No
    ## 1544 131.5  85.0     No
    ## 1545 117.5  72.0     No
    ## 1546 129.0  78.0     No
    ## 1547 112.0  66.5     No
    ## 1548 132.0  80.0     No
    ## 1549 141.0  92.0     No
    ## 1550 130.5  98.0     No
    ## 1551 122.5  82.5     No
    ## 1552 109.0  70.0     No
    ## 1553 141.0  75.0     No
    ## 1554 126.0  66.0     No
    ## 1555 135.0  88.0     No
    ## 1557 123.0  69.0     No
    ## 1558 105.0  70.0     No
    ## 1559 121.0  88.0     No
    ## 1560 123.0  68.0     No
    ## 1561 132.5  85.5     No
    ## 1562 122.0  70.0     No
    ## 1563 127.5  87.0     No
    ## 1564 100.0  68.0     No
    ## 1565 126.0  85.0     No
    ## 1566 129.0  61.0     No
    ## 1567 111.0  63.0     No
    ## 1568 142.5  79.0     No
    ## 1569 122.5  77.0     No
    ## 1571 143.0  94.5     No
    ## 1572 112.5  80.0     No
    ## 1573 107.0  73.0     No
    ## 1574 142.5  93.5     No
    ## 1575 125.0  85.0     No
    ## 1576 146.0  78.0     No
    ## 1577 143.5  90.0     No
    ## 1578 129.0  82.0     No
    ## 1579 132.5  85.0     No
    ## 1580 121.0  66.0     No
    ## 1581 112.5  60.0     No
    ## 1582 129.0  81.0     No
    ## 1583 142.0 108.0     No
    ## 1584 199.0 106.0     No
    ## 1586 134.0  92.0     No
    ## 1587 112.0  80.0     No
    ## 1589 179.0 107.0     No
    ## 1590 125.0  85.0     No
    ## 1591 179.0 100.0     No
    ## 1592 101.0  67.0     No
    ## 1593 104.0  57.0     No
    ## 1594 112.0  71.0     No
    ## 1595 104.0  66.0     No
    ## 1596 138.0  92.0     No
    ## 1597 130.0  85.0     No
    ## 1598 120.0  72.0     No
    ## 1599 124.0  72.0     No
    ## 1600 133.0  69.0     No
    ## 1602 108.0  81.0     No
    ## 1603 129.0  80.0     No
    ## 1605 118.0  82.0     No
    ## 1606 186.5  97.0     No
    ## 1607 131.0  99.0     No
    ## 1608 134.0  75.0     No
    ## 1609 116.0  83.0     No
    ## 1610 138.5  87.5     No
    ## 1611 120.0  72.0     No
    ## 1612 138.0  91.0     No
    ## 1613 127.5  85.5     No
    ## 1614 145.0 100.5     No
    ## 1615 150.0  90.0    Yes
    ## 1616 121.0  71.0     No
    ## 1617 120.0  80.0     No
    ## 1618 117.0  78.0     No
    ## 1619 155.0  90.0     No
    ## 1620 164.0 102.0     No
    ## 1621 122.5  76.5     No
    ## 1622  96.0  62.0     No
    ## 1623 170.0 103.0    Yes
    ## 1624 173.0  89.0     No
    ## 1625 134.0  80.0     No
    ## 1626 186.0 101.0     No
    ## 1627 154.0 106.0     No
    ## 1628 102.0  65.0     No
    ## 1629 160.0  94.0     No
    ## 1630 103.0  73.0     No
    ## 1631 115.0  84.0     No
    ## 1633 204.0 118.0     No
    ## 1634 114.0  78.0     No
    ## 1636 127.5  90.0     No
    ## 1637 116.0  82.0     No
    ## 1638 131.5  85.0     No
    ## 1642 158.0  94.0     No
    ## 1643 130.0  80.0     No
    ## 1644 123.5  84.5     No
    ## 1645 120.0  80.0     No
    ## 1646 111.0  95.5     No
    ## 1648 112.5  75.5     No
    ## 1649 118.5  84.5     No
    ## 1650 112.5  70.0     No
    ## 1651 140.0  93.0     No
    ## 1652 139.0  90.0     No
    ## 1653 124.0  79.0     No
    ## 1654 110.0  65.0     No
    ## 1655 217.0 112.0     No
    ## 1656 196.0 116.0     No
    ## 1657 130.0  75.5     No
    ## 1658 146.0  89.0     No
    ## 1659 123.0  83.0     No
    ## 1660 175.0 110.0     No
    ## 1661 114.0  77.0     No
    ## 1663 176.0  89.0     No
    ## 1664 140.0  94.0     No
    ## 1665 138.0  82.0     No
    ## 1666 111.0  78.0     No
    ## 1667 156.5  85.0     No
    ## 1668 113.0  66.5     No
    ## 1669 142.0  76.0     No
    ## 1670 167.0  94.0    Yes
    ## 1671 113.0  70.0     No
    ## 1672 193.0 104.0     No
    ## 1673 187.0  95.5     No
    ## 1675 107.5  75.0     No
    ## 1676 147.0  97.0     No
    ## 1677 112.5  72.5     No
    ## 1678 196.0 120.0     No
    ## 1679 189.0  87.0     No
    ## 1680 168.0  98.0     No
    ## 1681 130.0  85.0     No
    ## 1682 126.5  85.0     No
    ## 1683 136.0  76.0     No
    ## 1684 196.0 119.0     No
    ## 1685 112.0  79.0     No
    ## 1688 116.0  82.0     No
    ## 1689 170.0 118.0     No
    ## 1690 190.0  97.0     No
    ## 1691 120.0  70.5     No
    ## 1693 125.0  88.0     No
    ## 1695 158.0 102.0     No
    ## 1696 104.0  73.0     No
    ## 1697 107.5  70.0     No
    ## 1698 114.0  68.0     No
    ## 1699 125.5  82.5     No
    ## 1700 117.5  80.0     No
    ## 1701 131.0  85.0     No
    ## 1702 161.0 100.5     No
    ## 1703 114.0  79.0     No
    ## 1704 110.5  66.0     No
    ## 1705 132.5  82.5     No
    ## 1706 122.0  81.0     No
    ## 1707 130.0  85.5     No
    ## 1708 121.5  72.5     No
    ## 1709 129.0  85.0     No
    ## 1710 151.5  96.5     No
    ## 1711 122.0  82.0     No
    ## 1712 185.0  86.0     No
    ## 1714 110.0  74.0     No
    ## 1715 113.5  72.5     No
    ## 1716 119.0  76.0     No
    ## 1718 154.0  98.0     No
    ## 1720 123.0  82.0     No
    ## 1721 155.5  99.5     No
    ## 1722 122.0  84.0     No
    ## 1723 111.5  74.0     No
    ## 1724 119.0  84.0     No
    ## 1725 101.0  62.0     No
    ## 1727 131.0  80.0     No
    ## 1728 130.0  86.0     No
    ## 1729 127.0  86.0     No
    ## 1731 143.5 109.0     No
    ## 1732 159.0  91.5     No
    ## 1733 121.5  69.0     No
    ## 1734 139.0  75.0     No
    ## 1735 155.0  81.0     No
    ## 1736 115.0  61.0     No
    ## 1737 124.0  78.0     No
    ## 1738 109.0  71.0     No
    ## 1739 145.0 100.0     No
    ## 1740 104.0  69.0     No
    ## 1741 128.0  92.5     No
    ## 1742 114.5  80.0     No
    ## 1744 130.0  94.0     No
    ## 1745 122.0  82.0     No
    ## 1746 115.0  80.0     No
    ## 1748 126.5  90.0     No
    ## 1749 144.0  90.0     No
    ## 1750 120.0  76.0     No
    ## 1751 117.5  73.5     No
    ## 1752 143.5  90.0     No
    ## 1754 103.0  65.0     No
    ## 1755 115.0  71.0     No
    ## 1756 112.0  63.0     No
    ## 1757 146.0  91.0     No
    ## 1758  92.0  69.0     No
    ## 1759 133.0  86.0     No
    ## 1761 193.0 132.0     No
    ## 1762 121.0  86.0     No
    ## 1763 119.0  72.0     No
    ## 1764 141.5  93.0     No
    ## 1765 120.0  87.0     No
    ## 1766 169.0  85.0     No
    ## 1768 134.0  90.0     No
    ## 1770 110.0  76.0     No
    ## 1771 131.0  87.0     No
    ## 1772 117.5  71.0     No
    ## 1773 140.0  93.0     No
    ## 1774 147.5  95.0     No
    ## 1775 110.0  69.0     No
    ## 1776 120.0  70.0     No
    ## 1777 120.0  78.0     No
    ## 1778 133.5  81.5     No
    ## 1779 146.5  77.5     No
    ## 1780 110.0  71.5     No
    ## 1782 120.0  81.0     No
    ## 1783 107.0  74.0     No
    ## 1784 130.0  85.0     No
    ## 1785 146.5  80.0     No
    ## 1786 102.0  64.0     No
    ## 1787 126.5  84.0     No
    ## 1788 146.5  97.0     No
    ## 1789 118.0  76.0     No
    ## 1790 138.5  99.0    Yes
    ## 1792 166.0  93.0    Yes
    ## 1793 145.5  82.0     No
    ## 1794 116.5  82.0     No
    ## 1797 142.0  54.0     No
    ## 1799 122.0  70.0     No
    ## 1800 110.0  71.0     No
    ## 1801 119.0  85.5     No
    ## 1802 153.0  75.0     No
    ## 1803 166.5 107.0     No
    ## 1804 135.0  80.0     No
    ## 1806 124.0  75.0     No
    ## 1807 151.0 102.0     No
    ## 1808 124.0  69.0     No
    ## 1809 119.0  78.0     No
    ## 1810 106.0  72.0     No
    ## 1811 148.0  89.0     No
    ## 1812 124.0  87.0     No
    ## 1813 124.0  78.0     No
    ## 1814 141.0  83.5     No
    ## 1815 122.0  70.0     No
    ## 1816 105.5  57.5     No
    ## 1817 127.0  74.0     No
    ## 1818 120.0  77.5     No
    ## 1819 159.0  95.0     No
    ## 1821 123.0  86.0     No
    ## 1822 130.0  80.0     No
    ## 1825 166.5 106.5     No
    ## 1826 123.0  75.0     No
    ## 1827 104.0  65.0     No
    ## 1828 202.0 132.0     No
    ## 1829 131.5  82.0     No
    ## 1830 110.0  70.0     No
    ## 1831 100.0  64.0     No
    ## 1832 116.0  79.5     No
    ## 1833 132.0  82.0     No
    ## 1834 111.0  79.0     No
    ## 1836 173.0 102.0     No
    ## 1837 177.0 124.0     No
    ## 1838 161.0  97.0     No
    ## 1839 159.0 105.0    Yes
    ## 1840 117.0  70.0     No
    ## 1841 154.0  82.0     No
    ## 1842 148.0  74.0     No
    ## 1843 123.0  90.0     No
    ## 1844 136.5  86.5     No
    ## 1845 135.0  83.0    Yes
    ## 1846 109.0  66.5     No
    ## 1848 104.0  64.0     No
    ## 1849 170.0 105.0     No
    ## 1850 150.5  98.0     No
    ## 1851 138.0  78.0     No
    ## 1853 129.0  87.0     No
    ## 1854 124.0  77.0     No
    ## 1856 103.5  60.0     No
    ## 1857 134.0  84.0     No
    ## 1858 158.0  93.0     No
    ## 1859 171.0  89.0     No
    ## 1860 134.0  84.0     No
    ## 1861 110.0  73.0     No
    ## 1864 147.0  88.0     No
    ## 1865 119.5  73.0     No
    ## 1866 150.0  99.0     No
    ## 1867 140.0  88.0     No
    ## 1868 132.0  94.0     No
    ## 1869 123.5  83.0     No
    ## 1870 120.0  70.0     No
    ## 1871 129.0  80.0     No
    ## 1872 119.0  80.0     No
    ## 1873 138.0  90.0     No
    ## 1874 155.5 100.5     No
    ## 1875 105.0  65.0     No
    ## 1876 132.0  82.0     No
    ## 1877 119.0  83.0     No
    ## 1878 129.0  83.5     No
    ## 1881 123.0  79.0     No
    ## 1882  97.0  64.0     No
    ## 1883 158.0  98.0    Yes
    ## 1884 155.0  89.0     No
    ## 1885 112.5  75.5     No
    ## 1886 105.5  67.5     No
    ## 1887 135.0  83.0     No
    ## 1888 140.0  91.0     No
    ## 1889 107.0  67.5     No
    ## 1891 135.0  80.0     No
    ## 1892 141.0  87.0     No
    ## 1893 132.5  85.0     No
    ## 1894 131.0  92.0     No
    ## 1895 135.0  88.0     No
    ## 1896 125.0  79.0     No
    ## 1897 136.5  89.0     No
    ## 1899 120.0  70.0     No
    ## 1900 158.0  70.0     No
    ## 1901 128.0  71.0     No
    ## 1902 130.0  83.0     No
    ## 1903 100.0  70.0     No
    ## 1904 144.0  79.0     No
    ## 1906 120.0  82.5     No
    ## 1907 129.0  72.0     No
    ## 1909 122.0  69.0     No
    ## 1910 141.0  92.0     No
    ## 1911 102.5  60.0     No
    ## 1912 124.5  99.0     No
    ## 1913 163.0  82.0     No
    ## 1914 115.0  79.5     No
    ## 1916 128.0  72.0     No
    ## 1917 154.0  94.0     No
    ## 1918 129.0  85.5     No
    ## 1919 102.0  60.0     No
    ## 1920 180.0 100.0     No
    ## 1921 115.0  70.0     No
    ## 1922 150.0  86.0     No
    ## 1923 154.0 110.0     No
    ## 1924 142.5 104.5     No
    ## 1925 121.0  78.0     No
    ## 1926 140.0  88.0     No
    ## 1927 107.0  73.0     No
    ## 1928 156.0  91.0     No
    ## 1929 127.0  83.0     No
    ## 1930 107.5  66.5     No
    ## 1932 121.5  74.5     No
    ## 1933 120.0  67.5     No
    ## 1934 126.5  84.0     No
    ## 1935 104.0  73.5     No
    ## 1936 145.0  87.0     No
    ## 1937 134.0  82.0     No
    ## 1938 105.0  72.5     No
    ## 1939 123.0  78.0     No
    ## 1940 195.0  90.0    Yes
    ## 1941 151.0  88.0     No
    ## 1942 120.0  67.0    Yes
    ## 1944 119.5  79.0     No
    ## 1945 167.0 114.0     No
    ## 1946 157.5  78.0     No
    ## 1949 137.0  91.0     No
    ## 1950 121.0  78.0     No
    ## 1951 122.0  78.0     No
    ## 1952 114.0  87.0     No
    ## 1953 155.0  93.5     No
    ## 1955 124.0  77.0     No
    ## 1956 119.0  82.0     No
    ## 1958 127.0  76.0     No
    ## 1959 128.5  82.0     No
    ## 1960 126.0  80.0     No
    ## 1961 130.0  81.0     No
    ## 1962 137.0  85.0     No
    ## 1963 119.0  75.5     No
    ## 1964 124.0  76.0     No
    ## 1965 127.5  91.5     No
    ## 1966 200.0 120.0     No
    ## 1968 147.5  87.5     No
    ## 1969 135.0  76.5     No
    ## 1970 126.0  76.0     No
    ## 1972 144.0  91.0     No
    ## 1974 113.0  70.0     No
    ## 1975 125.0  90.0     No
    ## 1976 146.0  82.0     No
    ## 1977 146.0  94.5     No
    ## 1978 132.0  81.0     No
    ## 1979 120.0  85.0     No
    ## 1980 152.0 104.0     No
    ## 1981 127.0  90.0     No
    ## 1982 136.5  92.0     No
    ## 1983 112.5  67.0     No
    ## 1984 135.5  90.0     No
    ## 1985 128.0  78.0     No
    ## 1986 113.0  77.0     No
    ## 1987 127.0  82.0     No
    ## 1988 137.5  87.0     No
    ## 1989 112.0  74.0     No
    ## 1990 141.0  78.0     No
    ## 1991 180.0 114.0    Yes
    ## 1992 107.0  76.0     No
    ## 1993 113.5  70.0     No
    ## 1994 145.0  85.0     No
    ## 1995 110.0  74.0     No
    ## 1996 133.0  78.0     No
    ## 1997 113.0  65.0     No
    ## 1998 119.0  73.5     No
    ## 1999 107.5  65.0     No
    ## 2000 126.0  73.0     No
    ## 2001 118.0  84.0     No
    ## 2003 108.0  73.0     No
    ## 2005 111.5  72.5     No
    ## 2006 127.5  70.0     No
    ## 2007 147.0  96.0     No
    ## 2008 170.0 110.0     No
    ## 2009 130.0  75.0     No
    ## 2010 132.0  82.0     No
    ## 2011 143.0  95.0     No
    ## 2012 116.5  77.5     No
    ## 2014 105.0  67.5     No
    ## 2015 146.0  85.0     No
    ## 2017 133.5  80.0     No
    ## 2018 129.0  91.5     No
    ## 2019 112.0  82.0     No
    ## 2020 158.0  90.0     No
    ## 2021 150.0  80.0     No
    ## 2023 114.0  68.0     No
    ## 2024 135.0  82.0     No
    ## 2025 111.0  79.5     No
    ## 2026 112.5  70.0     No
    ## 2027 173.0 112.0     No
    ## 2028 142.0  79.0     No
    ## 2029 129.0  88.0     No
    ## 2030 142.0  84.0     No
    ## 2032 118.0  84.0     No
    ## 2033 134.0  72.0     No
    ## 2034 155.0  92.5     No
    ## 2035 108.0  65.5     No
    ## 2036 137.5 101.5     No
    ## 2037 122.0  82.0     No
    ## 2038 119.0  81.0     No
    ## 2039 182.0  92.0    Yes
    ## 2040 150.0  92.0     No
    ## 2041 145.0  90.0     No
    ## 2042 145.0  94.0     No
    ## 2043 157.0  71.0     No
    ## 2045 134.0  84.0     No
    ## 2046 120.5  73.5     No
    ## 2047 128.0  83.5     No
    ## 2048 121.0  85.5     No
    ## 2049 175.0  85.0     No
    ## 2050 111.0  70.0     No
    ## 2052 111.0  85.0     No
    ## 2053 142.0  83.0     No
    ## 2054 119.5  69.5     No
    ## 2055 105.0  72.5     No
    ## 2056 125.0  83.0     No
    ## 2057 145.5  99.0     No
    ## 2058  97.0  67.0     No
    ## 2059 119.0  73.0     No
    ## 2060 139.0  80.0     No
    ## 2061 126.0  89.5     No
    ## 2062 112.0  75.0     No
    ## 2063 171.0  87.0     No
    ## 2064 139.0  88.0     No
    ## 2065 132.0  82.0     No
    ## 2067 116.0  86.0     No
    ## 2068 112.5  75.0     No
    ## 2069 151.0  95.0     No
    ## 2072 126.0  79.0     No
    ## 2073 108.0  82.0     No
    ## 2075 126.0  72.0     No
    ## 2077 116.0  74.0     No
    ## 2078 105.0  70.0     No
    ## 2084 115.5  65.0     No
    ## 2085 232.0 136.0     No
    ## 2086 128.5  73.5     No
    ## 2087 119.0  80.0     No
    ## 2088 122.5  80.0     No
    ## 2089 134.0  88.0     No
    ## 2090 173.0  96.0    Yes
    ## 2091 132.0  85.0    Yes
    ## 2092 111.5  77.0     No
    ## 2093  85.5  51.0     No
    ## 2095 132.0  81.0     No
    ## 2096 110.0  80.0     No
    ## 2097 134.0  80.0     No
    ## 2100 112.0  62.0     No
    ## 2101 102.0  66.5     No
    ## 2103 162.5 100.0     No
    ## 2104 191.0  81.0    Yes
    ## 2106 133.0  77.0     No
    ## 2108 168.0  94.0     No
    ## 2109 148.0  99.0     No
    ## 2110 118.0  80.0     No
    ## 2112 108.0  70.5     No
    ## 2113 128.0  82.0     No
    ## 2115 110.0  80.0     No
    ## 2116 126.0  91.0     No
    ## 2117 141.5  93.0     No
    ## 2119 142.0  91.0     No
    ## 2120 132.0  82.0     No
    ## 2122 150.5  97.0     No
    ## 2123 184.5  83.0     No
    ## 2124 175.0  80.0     No
    ## 2126  95.0  59.0     No
    ## 2127 122.0  74.5     No
    ## 2128 135.0  80.0     No
    ## 2129 145.0  92.0     No
    ## 2130 144.0  82.0     No
    ## 2131 100.0  80.0     No
    ## 2132 118.5  88.0     No
    ## 2133 113.0  79.0     No
    ## 2134 130.0  87.0     No
    ## 2135 150.0  93.0     No
    ## 2136 139.0  82.5     No
    ## 2137 103.0  72.5     No
    ## 2139 126.0  91.0     No
    ## 2140 100.0  70.0     No
    ## 2141 157.0  97.0     No
    ## 2142 156.0  95.0     No
    ## 2143 134.5  92.5     No
    ## 2144 130.0  87.0     No
    ## 2145 112.5  77.5     No
    ## 2147 116.0  70.0     No
    ## 2148 105.0  59.0     No
    ## 2151 122.0  81.5     No
    ## 2152 107.5  74.0     No
    ## 2153 120.0  80.0     No
    ## 2154 105.0  72.5     No
    ## 2155 130.0  80.0     No
    ## 2156 146.5  82.0     No
    ## 2157 126.0  78.0     No
    ## 2158 114.0  75.0     No
    ## 2159 112.5  75.0     No
    ## 2160 117.0  77.0     No
    ## 2162 138.5  90.0     No
    ## 2163 120.0  73.0     No
    ## 2164 122.0  74.0     No
    ## 2166 149.0  86.0     No
    ## 2167 120.0  83.5     No
    ## 2170 117.0  84.0     No
    ## 2171 117.5  72.5     No
    ## 2172 124.0  77.0     No
    ## 2173 127.0  79.0     No
    ## 2174 116.0  72.0     No
    ## 2176 112.0  73.5     No
    ## 2177 121.0  70.0     No
    ## 2178 145.0  86.0     No
    ## 2179  95.5  64.0     No
    ## 2181 117.0  78.0     No
    ## 2183 130.0  72.0     No
    ## 2184 152.0  98.0     No
    ## 2185 134.0  87.0     No
    ## 2187 155.0  84.0     No
    ## 2188 141.0  90.0     No
    ## 2189 171.0 120.0     No
    ## 2190 102.0  69.0     No
    ## 2191 113.0  69.0     No
    ## 2192 186.0 109.0    Yes
    ## 2194 123.0  74.0     No
    ## 2195 188.0 128.0     No
    ## 2196 161.0 104.0    Yes
    ## 2197 107.0  76.0     No
    ## 2198 124.0  81.0     No
    ## 2200 162.5  93.5     No
    ## 2201 104.0  72.0     No
    ## 2202 157.0  98.0     No
    ## 2203 114.0  80.0     No
    ## 2204 119.0  74.0     No
    ## 2205 123.0  94.5     No
    ## 2208 130.0  86.0     No
    ## 2209 120.0  83.0     No
    ## 2214 144.5  88.5     No
    ## 2215 126.5  76.0     No
    ## 2216 119.0  85.0     No
    ## 2217 107.5  75.0     No
    ## 2218 146.0  76.0     No
    ## 2219 113.0  73.0     No
    ## 2220 127.0  79.0     No
    ## 2222 146.0  98.5     No
    ## 2224 137.0  93.0     No
    ## 2225 104.0  61.0     No
    ## 2226 135.0  87.0     No
    ## 2227 167.5 102.5     No
    ## 2228 136.0  88.0     No
    ## 2229 154.0  96.0     No
    ## 2230 141.0  98.0     No
    ## 2231 146.0  77.0     No
    ## 2232 133.0  83.0     No
    ## 2233 134.0  85.0     No
    ## 2234 125.0  82.0     No
    ## 2235 205.0  92.5     No
    ## 2236 169.5 104.5     No
    ## 2238 125.5  80.0     No
    ## 2239 185.0 107.5     No
    ## 2240 120.5  76.0     No
    ## 2242  96.5  67.0     No
    ## 2243 152.0  97.0     No
    ## 2244 130.0  88.0     No
    ## 2245 132.0  85.5     No
    ## 2246 127.0  86.0     No
    ## 2247 129.0  94.0     No
    ## 2248 149.5  85.0     No
    ## 2249 110.0  69.0     No
    ## 2250 124.0  74.0     No
    ## 2251 148.0  92.0     No
    ## 2252 106.0  59.5     No
    ## 2253 126.0  87.0     No
    ## 2254 121.0  79.5     No
    ## 2255 110.5  69.0     No
    ## 2256 146.5  80.0     No
    ## 2257 116.0  69.0     No
    ## 2258 130.0  95.0     No
    ## 2259 110.0  78.0     No
    ## 2260 116.0  77.0     No
    ## 2261 129.0  76.5     No
    ## 2262 188.0  92.0     No
    ## 2263 115.0  72.0     No
    ## 2266 220.0  96.0     No
    ## 2267 116.0  79.0     No
    ## 2268 116.0  83.0     No
    ## 2269 124.0  82.0     No
    ## 2270 146.0  90.0     No
    ## 2271 127.0  90.0     No
    ## 2272 111.0  72.0     No
    ## 2273 134.0  83.0     No
    ## 2274 182.0 111.0     No
    ## 2275 118.0  79.5     No
    ## 2276 124.0  86.0     No
    ## 2279 150.0  88.0     No
    ## 2280 107.5  67.5     No
    ## 2282 124.0  79.0     No
    ## 2283 210.0 120.0     No
    ## 2284 106.0  64.0     No
    ## 2285 109.0  75.0     No
    ## 2287 120.0  74.0     No
    ## 2289 132.0  81.0     No
    ## 2290 130.0  80.0     No
    ## 2291 151.0  92.0     No
    ## 2292 184.0 102.0    Yes
    ## 2293 105.0  73.0     No
    ## 2294 130.0  91.0     No
    ## 2295 134.5  91.0     No
    ## 2296 115.0  80.0     No
    ## 2297 127.0  86.0     No
    ## 2298  97.0  62.0     No
    ## 2299 128.0  83.0     No
    ## 2300 124.0  81.0     No
    ## 2301 114.0  70.0     No
    ## 2302 193.0  63.0     No
    ## 2303 121.5  79.0     No
    ## 2304 120.0  80.0     No
    ## 2305 154.0  98.0     No
    ## 2306 128.0  79.0     No
    ## 2307 115.0  69.0     No
    ## 2308 120.0  75.0     No
    ## 2309 158.0  96.0     No
    ## 2310 120.0  80.0     No
    ## 2311 181.5 102.5     No
    ## 2312 109.0  73.0     No
    ## 2313 115.0  71.0     No
    ## 2315 130.0  75.0     No
    ## 2316 185.0 105.0     No
    ## 2317 150.0 101.0     No
    ## 2318 154.0  84.0     No
    ## 2319 139.0  80.0     No
    ## 2320 188.5 106.5     No
    ## 2321 127.5  72.5     No
    ## 2323 106.0  80.0     No
    ## 2324 132.5  85.0     No
    ## 2326 120.0  83.0     No
    ## 2327 122.5  68.5     No
    ## 2328 165.0  95.0     No
    ## 2329 116.0  75.0     No
    ## 2330 113.0  68.0     No
    ## 2331 122.5  80.0     No
    ## 2332 142.5  74.5     No
    ## 2333 170.0  92.0     No
    ## 2334 132.5  80.0     No
    ## 2335 118.0  76.0     No
    ## 2336 119.0  83.5     No
    ## 2337 110.0  79.0     No
    ## 2338 112.0  70.0     No
    ## 2339 107.5  68.0     No
    ## 2340 125.0  82.0     No
    ## 2341 111.5  67.0     No
    ## 2342 139.0  80.0     No
    ## 2343 182.5  97.5     No
    ## 2344 165.0  84.0     No
    ## 2345 145.0  92.5     No
    ## 2346 108.0  75.0     No
    ## 2347 138.0  97.0     No
    ## 2348 116.0  68.0     No
    ## 2349 138.0  72.0     No
    ## 2350 108.0  74.0     No
    ## 2351 108.0  62.0     No
    ## 2352 164.0 111.0     No
    ## 2353 148.0 108.0     No
    ## 2354 131.5  91.0     No
    ## 2355  96.0  67.0     No
    ## 2356 120.0  60.0     No
    ## 2358 114.0  73.5     No
    ## 2359 148.0 103.0     No
    ## 2360 140.0  86.0     No
    ## 2362 120.0  78.0     No
    ## 2363 130.5  90.0     No
    ## 2364 128.5  87.5     No
    ## 2365 116.0  67.0     No
    ## 2366 120.0  80.0     No
    ## 2367 100.0  70.0     No
    ## 2368 109.5  72.0     No
    ## 2371 132.0  86.0     No
    ## 2372 115.0  75.0     No
    ## 2373 139.0  80.0     No
    ## 2375 148.0  75.0     No
    ## 2378 133.0  96.0     No
    ## 2379 124.5  66.5     No
    ## 2381 134.0  88.0     No
    ## 2382 174.0  90.0    Yes
    ## 2383 110.0  70.0     No
    ## 2384 120.0  70.0     No
    ## 2385 192.0 105.0    Yes
    ## 2386 109.0  77.0     No
    ## 2387 138.0  79.0     No
    ## 2388 146.0  92.0     No
    ## 2389 176.5 115.0     No
    ## 2390 110.0  71.0     No
    ## 2391 101.0  69.0     No
    ## 2392 148.0  91.0     No
    ## 2393 110.0  78.0     No
    ## 2394 137.0  84.0     No
    ## 2395 154.5  83.0     No
    ## 2396 127.5  62.5     No
    ## 2397 110.5  69.0     No
    ## 2398 153.0 100.0     No
    ## 2399 177.0  97.0     No
    ## 2400 127.5  83.5     No
    ## 2401 119.0  76.0     No
    ## 2402 140.0  81.0     No
    ## 2404 132.0  92.0     No
    ## 2405 133.0  77.5     No
    ## 2406 160.0  96.0     No
    ## 2407 113.5  72.0     No
    ## 2408 149.0  95.0     No
    ## 2409 183.0 108.0    Yes
    ## 2410 129.0  86.0     No
    ## 2411 116.0  85.5     No
    ## 2412 123.0  92.0     No
    ## 2413 114.0  68.0     No
    ## 2414 199.0 114.0     No
    ## 2415 159.0  92.0     No
    ## 2416 171.0  97.0     No
    ## 2417 197.5 125.0     No
    ## 2418 190.0 100.0     No
    ## 2419 116.0  88.5     No
    ## 2420 116.0  81.0     No
    ## 2421 112.5  80.0     No
    ## 2422 120.0  83.5     No
    ## 2423 116.0  77.5     No
    ## 2426 174.0 100.0     No
    ## 2427 133.0  84.0     No
    ## 2428 107.5  72.0     No
    ## 2429 110.0  66.0     No
    ## 2430 133.0  87.0     No
    ## 2431 116.0  83.0     No
    ## 2432 195.0 108.0     No
    ## 2433 140.0  85.0     No
    ## 2434 145.0  85.0     No
    ## 2435 135.0  75.0     No
    ## 2436 101.0  71.0     No
    ## 2437 143.0  81.0     No
    ## 2438 145.0  88.0     No
    ## 2439 122.0  78.0     No
    ## 2440 111.0  71.0     No
    ## 2441 146.0  78.0     No
    ## 2442 122.0  87.5     No
    ## 2443 120.0  80.0     No
    ## 2444 150.0  93.0     No
    ## 2445 180.0 100.0     No
    ## 2446 109.0  70.0     No
    ## 2447 115.5  79.0     No
    ## 2448 121.0  78.5     No
    ## 2449 132.0  85.0     No
    ## 2453 102.0  67.0     No
    ## 2454 110.0  72.5     No
    ## 2455 115.0  80.0     No
    ## 2456 136.5  85.0     No
    ## 2458 118.0  72.0     No
    ## 2459 144.0  79.0     No
    ## 2460 175.0  78.0     No
    ## 2461 120.0  81.0     No
    ## 2462 146.0 106.0     No
    ## 2463 150.0 101.0     No
    ## 2464 132.0  81.0     No
    ## 2465 125.0  75.0     No
    ## 2467 139.0  92.5     No
    ## 2468 165.0 106.0     No
    ## 2469 126.0  86.0     No
    ## 2471 105.0  72.5     No
    ## 2472 128.0  81.5     No
    ## 2473 129.0  74.0     No
    ## 2474 132.0  96.0     No
    ## 2475 128.0  83.0     No
    ## 2476 113.0  81.0     No
    ## 2477 131.0  96.0     No
    ## 2478 134.5  84.0     No
    ## 2479 111.0  70.0     No
    ## 2480 118.0  84.0     No
    ## 2482 150.0  94.0     No
    ## 2483 142.5  87.0     No
    ## 2484 148.0  92.5     No
    ## 2485 146.0  80.0     No
    ## 2486 176.0  87.0     No
    ## 2487 158.5  94.5     No
    ## 2489 115.0  75.0     No
    ## 2490 144.0  88.0     No
    ## 2491 133.0  87.0     No
    ## 2493 137.0  80.0     No
    ## 2494 135.0  86.0     No
    ## 2495 145.0  72.5     No
    ## 2496 166.0  96.0     No
    ## 2497 136.0  96.0     No
    ## 2498 121.0  88.0     No
    ## 2499 137.0  88.0     No
    ## 2500 146.0  91.0     No
    ## 2501 132.5  92.0     No
    ## 2502 142.0  85.0     No
    ## 2503 114.0  76.0     No
    ## 2504 108.0  62.0     No
    ## 2505 129.0  89.0     No
    ## 2506 112.0  80.0     No
    ## 2507 113.0  61.0     No
    ## 2508 137.0  95.0     No
    ## 2509 137.5  88.5     No
    ## 2510 102.0  61.0     No
    ## 2511 208.0 104.0     No
    ## 2512 137.0  83.5     No
    ## 2513 122.0  95.0     No
    ## 2514 118.5  76.0     No
    ## 2515 110.0  68.0     No
    ## 2516 136.0  88.0     No
    ## 2518 124.0  78.0     No
    ## 2519 105.0  69.0     No
    ## 2520 120.0  87.0     No
    ## 2521 130.0  85.0     No
    ## 2522 136.5  97.0     No
    ## 2523 125.0  82.0     No
    ## 2524 123.0  69.0     No
    ## 2525 130.0  87.0     No
    ## 2526 123.5  78.0     No
    ## 2527 150.0  89.0     No
    ## 2529 118.0  70.0     No
    ## 2530 102.5  64.5     No
    ## 2531 113.0  79.0     No
    ## 2532 129.0  81.0     No
    ## 2533  94.0  62.0     No
    ## 2534 163.0 102.0     No
    ## 2535 210.0 130.0     No
    ## 2537 132.0  91.0     No
    ## 2538 112.5  75.0     No
    ## 2539 184.5 110.5     No
    ## 2540 141.0  82.5     No
    ## 2541  92.5  70.0     No
    ## 2542 126.0  85.0     No
    ## 2543 172.0 105.0     No
    ## 2544 120.0  82.0     No
    ## 2545 128.0  87.0     No
    ## 2546 110.0  73.0     No
    ## 2547 102.0  68.0     No
    ## 2548 115.0  79.0     No
    ## 2551 132.0  76.0     No
    ## 2552 128.0  74.0     No
    ## 2553 125.0  88.0     No
    ## 2554 127.5  77.5     No
    ## 2555 145.0  89.0     No
    ## 2556 165.0  99.0     No
    ## 2557 144.0  98.0     No
    ## 2558 127.0  82.5     No
    ## 2559 126.0  82.5     No
    ## 2560 117.5  77.5     No
    ## 2561 160.0  97.5     No
    ## 2562 136.5  87.0     No
    ## 2563 138.0  91.0     No
    ## 2564 127.0  68.5     No
    ## 2565 103.0  76.5     No
    ## 2566 172.0  95.0    Yes
    ## 2567 129.5  82.5     No
    ## 2568 176.0 113.0     No
    ## 2569 125.0  94.0     No
    ## 2570 116.0  82.0     No
    ## 2571 162.0  85.0     No
    ## 2572 102.0  69.0     No
    ## 2573 136.0  86.0     No
    ## 2574 124.0  72.5     No
    ## 2575 127.0  78.0     No
    ## 2576 115.0  76.0     No
    ## 2577 109.0  75.0     No
    ## 2578 128.0  91.0     No
    ## 2579 202.5  85.0     No
    ## 2580 120.5  78.0     No
    ## 2581 140.0 100.0     No
    ## 2583 136.5  85.0     No
    ## 2584 145.0  88.0    Yes
    ## 2585 112.0  67.0     No
    ## 2586 102.0  72.0     No
    ## 2587 156.0  69.0     No
    ## 2588 140.0  97.5    Yes
    ## 2590 155.0  85.0    Yes
    ## 2591 135.0  90.0     No
    ## 2592 120.0  72.0     No
    ## 2593 136.0  94.0     No
    ## 2594 139.0  74.0     No
    ## 2595 154.0  97.0    Yes
    ## 2597 142.5  83.5     No
    ## 2598 166.0 101.0     No
    ## 2599 128.0  85.0     No
    ## 2600 122.5  78.5     No
    ## 2601 146.0  92.0     No
    ## 2602 136.5  87.5     No
    ## 2603 112.0  66.0     No
    ## 2604  98.0  74.0     No
    ## 2605 130.0  84.0     No
    ## 2606 122.0  80.0     No
    ## 2607 112.5  68.0     No
    ## 2608 131.0  84.0     No
    ## 2609 150.0  95.5     No
    ## 2610 141.0  84.5     No
    ## 2611 166.0  85.0     No
    ## 2612 120.0  81.0     No
    ## 2613 141.0  77.5     No
    ## 2614 120.0  81.0     No
    ## 2615 175.0  82.0    Yes
    ## 2616 146.0  98.5     No
    ## 2617 164.0 104.0     No
    ## 2618 138.0  89.0     No
    ## 2619 104.0  72.0     No
    ## 2620 130.0  80.0     No
    ## 2621 137.5  72.5     No
    ## 2622 121.0  74.0     No
    ## 2623 115.0  70.0     No
    ## 2624 140.0  84.0     No
    ## 2625 145.0  81.0     No
    ## 2626 141.0 105.0     No
    ## 2627 124.0  85.0     No
    ## 2628 124.0  78.0     No
    ## 2629 158.0  94.0     No
    ## 2630 117.0  83.0     No
    ## 2632 159.0 102.0    Yes
    ## 2633  96.0  59.0     No
    ## 2634 169.0 117.0     No
    ## 2635 102.5  65.0     No
    ## 2636 149.0  82.0    Yes
    ## 2637 119.0  80.0     No
    ## 2638 100.0  61.5     No
    ## 2639 133.5  80.0     No
    ## 2640 113.0  70.0     No
    ## 2641 114.0  82.0     No
    ## 2643 125.0  79.0     No
    ## 2644 128.0  82.0     No
    ## 2645 110.0  79.0     No
    ## 2646 110.0  71.0     No
    ## 2649 120.0  83.5     No
    ## 2650 110.0  70.0     No
    ## 2651 107.5  80.0     No
    ## 2652 122.0  81.0     No
    ## 2655 112.5  85.0     No
    ## 2656 149.5  93.0     No
    ## 2657 117.5  80.0     No
    ## 2658 159.0  95.0     No
    ## 2659 121.0  83.0     No
    ## 2660 113.0  82.5     No
    ## 2661 140.0  86.0     No
    ## 2662 126.0  73.0     No
    ## 2663 122.0  86.0     No
    ## 2664 110.0  71.0     No
    ## 2665 155.0 100.0     No
    ## 2667 144.0  90.0     No
    ## 2668 122.0  81.0     No
    ## 2669 123.0  80.0     No
    ## 2670 119.0  78.5     No
    ## 2671 182.5  97.0     No
    ## 2672 132.5  85.5     No
    ## 2673 166.0  90.0     No
    ## 2674 126.0  84.0     No
    ## 2675 137.5  94.0     No
    ## 2676 131.5  89.0     No
    ## 2677 123.0  76.5     No
    ## 2679 141.0  81.5     No
    ## 2680 147.5  87.5     No
    ## 2681 154.5 104.0     No
    ## 2682 125.0  85.0     No
    ## 2683 105.0  71.0     No
    ## 2684 122.0  87.0    Yes
    ## 2685 135.0  85.0     No
    ## 2686 142.5  95.0     No
    ## 2687 119.0  72.0     No
    ## 2688 146.0  98.0     No
    ## 2689 160.0  90.0    Yes
    ## 2690 109.0  75.0     No
    ## 2692 112.0  78.0     No
    ## 2693 180.0 108.0     No
    ## 2695 115.0  80.0     No
    ## 2696 139.0  96.0     No
    ## 2698 147.0 102.0     No
    ## 2699 175.0  85.0     No
    ## 2700 112.5  85.0     No
    ## 2701 123.0  77.5     No
    ## 2702 116.5  83.0     No
    ## 2703 112.0  77.0     No
    ## 2704 131.0  71.0     No
    ## 2705 122.0  78.0     No
    ## 2706 115.0  84.0     No
    ## 2708 129.0  86.0     No
    ## 2710 165.0  99.0     No
    ## 2711 102.5  72.5     No
    ## 2712 142.0  92.0     No
    ## 2714 114.0  85.0     No
    ## 2715 119.0  75.0     No
    ## 2716 108.0  74.0     No
    ## 2718 105.0  70.0     No
    ## 2719 177.0  75.0     No
    ## 2721 107.5  66.0     No
    ## 2722 155.0  99.0     No
    ## 2723 144.0  82.0     No
    ## 2724 135.0  95.0     No
    ## 2725 182.0 102.0     No
    ## 2726 129.0  84.0     No
    ## 2727 130.0  82.5     No
    ## 2728 125.0  72.5     No
    ## 2729 137.5  82.5     No
    ## 2730 162.0  91.0     No
    ## 2731 170.0 104.0     No
    ## 2732 133.0  77.0     No
    ## 2735 110.0  72.0     No
    ## 2736 128.0  89.0     No
    ## 2738 129.0  73.0     No
    ## 2739 166.0  71.0     No
    ## 2741 103.0  78.0     No
    ## 2742 141.0  86.0    Yes
    ## 2743 139.5  87.0     No
    ## 2744 116.0  72.0     No
    ## 2745 191.0 106.0     No
    ## 2746 134.0  88.0     No
    ## 2747 113.0  84.0     No
    ## 2748 110.0  70.0     No
    ## 2749 133.0  89.0     No
    ## 2750 152.5  85.0     No
    ## 2752 132.0  80.0     No
    ## 2754 165.0  90.0     No
    ## 2755 116.0  69.0     No
    ## 2756 123.0  70.0     No
    ## 2757 139.0  67.0     No
    ## 2758 131.0  79.0     No
    ## 2759 121.5  86.5     No
    ## 2760 182.5 103.0     No
    ## 2761 134.0  87.0     No
    ## 2762 118.0  83.0     No
    ## 2763 105.0  70.0     No
    ## 2764 114.0  72.5     No
    ## 2766 127.0  79.5     No
    ## 2768 121.0  61.0     No
    ## 2769 147.0  86.5     No
    ## 2770 125.0  85.0     No
    ## 2771 118.0  79.0     No
    ## 2773 160.0  82.0     No
    ## 2774 120.0  73.5     No
    ## 2775 114.0  78.0     No
    ## 2776 152.0  88.0     No
    ## 2777 147.5  93.0     No
    ## 2778 143.5  90.0     No
    ## 2779 131.0  74.0     No
    ## 2780 122.5  77.5     No
    ## 2781 111.0  81.0     No
    ## 2783 147.0  98.0    Yes
    ## 2784 119.0  80.0     No
    ## 2785 190.0 130.0     No
    ## 2787 145.0  88.0     No
    ## 2788 119.0  82.0     No
    ## 2789 120.0  75.0     No
    ## 2790 151.5  79.0     No
    ## 2791 127.5  90.0     No
    ## 2792 125.5  80.5     No
    ## 2793  83.5  58.0     No
    ## 2794  96.0  67.0     No
    ## 2795 152.0  73.0     No
    ## 2796 161.0  90.0    Yes
    ## 2798 159.0  87.0    Yes
    ## 2799 144.0  74.0     No
    ## 2800 142.0  84.0     No
    ## 2801 137.5  81.0     No
    ## 2802 115.0  75.0     No
    ## 2803 152.5  88.0     No
    ## 2804 131.5  81.0     No
    ## 2805 115.0  69.0     No
    ## 2806 128.0  87.0     No
    ## 2807  98.0  73.0     No
    ## 2808 118.0  80.0     No
    ## 2809 127.0  68.0     No
    ## 2810 146.5  71.0     No
    ## 2811 150.0  94.0     No
    ## 2812 210.0 135.0     No
    ## 2813 117.5  77.5     No
    ## 2814 131.5  83.0     No
    ## 2815 171.0  84.0     No
    ## 2816 158.0 109.0     No
    ## 2817 121.0  75.0     No
    ## 2819 116.0  69.0     No
    ## 2820 110.0  80.0     No
    ## 2821 115.0  63.5     No
    ## 2822 122.0  76.0     No
    ## 2823 115.0  81.0     No
    ## 2826 150.0  89.0     No
    ## 2827 148.0  95.0     No
    ## 2829 124.0  75.0     No
    ## 2830  99.0  60.0     No
    ## 2831 136.5  84.0     No
    ## 2832 122.0  83.0     No
    ## 2833 128.0  90.0     No
    ## 2834 150.0  95.0     No
    ## 2835 144.0  85.0    Yes
    ## 2836 106.5  75.0     No
    ## 2837 132.0  90.0     No
    ## 2838 121.0  82.0     No
    ## 2841 113.0  74.0     No
    ## 2842 176.0  95.0     No
    ## 2843 155.0 105.0     No
    ## 2844 130.0  80.0     No
    ## 2845 121.0  78.0     No
    ## 2846 152.0  95.0     No
    ## 2847 143.0  76.0     No
    ## 2848 148.0  70.0     No
    ## 2850 101.0  68.0     No
    ## 2851 142.0  80.0     No
    ## 2852 123.0  77.0     No
    ## 2853 130.0  77.5     No
    ## 2854 177.0 101.0     No
    ## 2855 103.0  73.0     No
    ## 2856 116.0  68.0     No
    ## 2859 140.5  92.5     No
    ## 2860 120.0  84.0     No
    ## 2861 165.0 112.0     No
    ## 2862 140.0  92.0     No
    ## 2863 141.0 115.0    Yes
    ## 2864 170.5 100.5     No
    ## 2865 127.0  84.0     No
    ## 2866 142.0  82.5     No
    ## 2867 119.0  72.0     No
    ## 2869 108.0  76.0     No
    ## 2870 117.5  83.5     No
    ## 2871 132.0  85.0     No
    ## 2872 131.0  80.0     No
    ## 2873 107.0  74.0     No
    ## 2874 111.0  76.0     No
    ## 2875 156.0  86.0     No
    ## 2876 121.0  82.0     No
    ## 2877 108.0  73.0     No
    ## 2879 126.5  67.0     No
    ## 2880 132.0  94.0     No
    ## 2881 120.0  77.5     No
    ## 2883 107.0  68.0     No
    ## 2884 114.0  74.0     No
    ## 2885 175.0 114.0    Yes
    ## 2887 150.0  85.0     No
    ## 2888 138.0  78.0     No
    ## 2889 113.5  65.5     No
    ## 2890 132.5  87.0     No
    ## 2891 109.0  73.0     No
    ## 2892 126.0  80.0     No
    ## 2894 134.0  87.0     No
    ## 2895 107.0  65.0     No
    ## 2896 165.0  94.5     No
    ## 2897 111.0  60.0     No
    ## 2898 123.0  82.0     No
    ## 2899  93.0  71.0     No
    ## 2900 123.5  79.0     No
    ## 2902 115.0  72.0     No
    ## 2904 101.0  71.0     No
    ## 2905 171.0 118.0     No
    ## 2906 132.0  87.0     No
    ## 2907 118.5  72.0     No
    ## 2908 133.0  88.0     No
    ## 2909 181.0  97.5    Yes
    ## 2910 143.0  93.0     No
    ## 2911 102.5  71.0     No
    ## 2912 142.5 100.0     No
    ## 2913 126.0  85.5     No
    ## 2915 113.5  61.0     No
    ## 2916 116.5  72.5     No
    ## 2917 128.0  76.0     No
    ## 2918 138.0  98.0     No
    ## 2919 148.0  94.0     No
    ## 2920 156.0 100.0     No
    ## 2921 142.5  83.5     No
    ## 2923 148.0  98.0     No
    ## 2924 138.0  78.0     No
    ## 2925 149.0 103.5     No
    ## 2926 166.0  98.0    Yes
    ## 2927 141.0  82.5     No
    ## 2928 111.0  80.0     No
    ## 2929 155.0  71.0     No
    ## 2931 197.0 109.0     No
    ## 2932 164.0  94.0     No
    ## 2933 156.0  91.0     No
    ## 2934 116.0  74.0     No
    ## 2936 146.0  89.0     No
    ## 2937 118.0  85.0     No
    ## 2938 110.0  70.0     No
    ## 2939 162.5  93.5    Yes
    ## 2940 122.5  67.5     No
    ## 2941 144.0  79.0     No
    ## 2942 105.0  60.0     No
    ## 2943 121.5  76.5     No
    ## 2944 114.0  80.0     No
    ## 2945 119.0  69.0     No
    ## 2946 102.5  69.0     No
    ## 2947 116.0  86.0     No
    ## 2948 123.0  86.0     No
    ## 2949 114.0  80.0     No
    ## 2950 181.0  74.0     No
    ## 2951 132.0  80.0     No
    ## 2952 129.0  90.0     No
    ## 2953 145.0  88.0     No
    ## 2954 116.0  81.0     No
    ## 2955 107.5  70.0     No
    ## 2957 117.5  77.5     No
    ## 2958 139.0  82.0     No
    ## 2959 135.0  87.5     No
    ## 2960 131.0  90.0     No
    ## 2961 108.0  73.0     No
    ## 2962 138.0  76.0     No
    ## 2963 138.0  85.0     No
    ## 2964 102.0  59.0     No
    ## 2965 125.0  84.0     No
    ## 2966 126.0  93.0     No
    ## 2968 198.0 108.0    Yes
    ## 2969 123.0  82.0     No
    ## 2970 135.0  97.5    Yes
    ## 2971 134.0  79.0     No
    ## 2972 163.0  87.0     No
    ## 2973 119.0  82.0     No
    ## 2974 115.0  75.0     No
    ## 2975 141.0  85.0     No
    ## 2976 114.5  76.0     No
    ## 2977 124.0  83.0     No
    ## 2980 163.0  83.0     No
    ## 2981 113.5  80.0     No
    ## 2982 131.0  79.0     No
    ## 2983 135.0  82.0     No
    ## 2984 150.0  88.0     No
    ## 2985 126.0  76.0     No
    ## 2986 107.0  82.0     No
    ## 2987 111.0  72.0     No
    ## 2988 135.0  77.5     No
    ## 2989 190.0 110.0     No
    ## 2990 107.0  82.0     No
    ## 2992 136.0  85.5     No
    ## 2993 106.0  75.0     No
    ## 2994 112.0  71.0     No
    ## 2995 125.0  80.0     No
    ## 2996 122.5  82.0     No
    ## 2997 175.5 113.0    Yes
    ## 2998  96.5  72.5     No
    ## 2999 132.0  95.0     No
    ## 3000 114.0  64.0     No
    ## 3001 130.0  79.0     No
    ## 3002 119.0  86.0     No
    ## 3003 172.0  82.0     No
    ## 3004 118.0  71.0     No
    ## 3006 149.0  88.0     No
    ## 3007 122.0  85.0     No
    ## 3008 126.5  82.0     No
    ## 3009 142.5  86.5     No
    ## 3010 160.0  99.0     No
    ## 3011 103.0  64.5     No
    ## 3012 135.0  87.0     No
    ## 3013 128.0  89.0     No
    ## 3014 140.0  88.0     No
    ## 3016 181.0 107.0    Yes
    ## 3017 163.5 102.0     No
    ## 3019 122.0  78.0     No
    ## 3020 124.0  73.0     No
    ## 3021 166.0  89.0     No
    ## 3022 102.0  66.5     No
    ## 3023 103.0  66.0     No
    ## 3024 130.5  82.0     No
    ## 3025 113.0  62.0     No
    ## 3026 102.5  66.5     No
    ## 3027  90.0  70.0     No
    ## 3029 109.0  70.0     No
    ## 3030 124.0  72.5     No
    ## 3031 132.5  73.0     No
    ## 3032 171.0  95.0     No
    ## 3033 112.5  73.5     No
    ## 3036 120.0  62.0     No
    ## 3039 133.5  89.5     No
    ## 3040 149.0  89.0     No
    ## 3041 127.5  76.0     No
    ## 3042 204.0  96.0    Yes
    ## 3043 129.0  79.0     No
    ## 3045 125.0  86.0     No
    ## 3046 120.0  76.0     No
    ## 3047 118.0  70.5     No
    ## 3048 132.5  87.0     No
    ## 3049 120.0  72.0     No
    ## 3050 142.0  89.0     No
    ## 3051 122.0  82.5     No
    ## 3052 137.0  90.0     No
    ## 3053 151.5 103.0     No
    ## 3054 134.0  89.0     No
    ## 3055 132.5  69.0     No
    ## 3056 112.0  78.0     No
    ## 3057 107.5  73.5     No
    ## 3059 127.0  80.0     No
    ## 3060 140.5  89.0     No
    ## 3063 207.5 118.0     No
    ## 3064 115.0  75.0     No
    ## 3065 115.0  81.0     No
    ## 3066 127.0  80.0     No
    ## 3067 129.5  85.5     No
    ## 3068 102.0  59.5     No
    ## 3069 181.0 101.0     No
    ## 3070 191.0  97.0     No
    ## 3071 138.0  88.0     No
    ## 3072 145.0  91.0     No
    ## 3073 136.0  87.0     No
    ## 3074 110.0  79.0     No
    ## 3076 113.0  83.0     No
    ## 3079 150.5  94.0     No
    ## 3081 132.0  88.0     No
    ## 3082 116.0  72.0     No
    ## 3084 104.0  73.0     No
    ## 3085 121.0  70.0     No
    ## 3086 195.0 110.0     No
    ## 3087 130.5  86.0     No
    ## 3088 137.5  80.0     No
    ## 3089 110.0  77.0     No
    ## 3091 131.0  81.5     No
    ## 3092 147.5 100.0     No
    ## 3093 165.0  88.0     No
    ## 3094 142.0  88.0     No
    ## 3095 177.0 111.0     No
    ## 3096 108.0  75.0     No
    ## 3097 124.0  72.0     No
    ## 3098 145.0  67.0     No
    ## 3099 108.0  74.0     No
    ## 3100 155.0  92.5     No
    ## 3102 105.0  72.5     No
    ## 3103 147.5  97.5    Yes
    ## 3104 127.5  80.0     No
    ## 3106 136.5  87.0     No
    ## 3107 134.0  85.5     No
    ## 3108 174.5 101.5     No
    ## 3109 126.0  88.0     No
    ## 3110 126.0  77.0     No
    ## 3111 133.5  86.0     No
    ## 3112 125.0  86.0     No
    ## 3113 148.5  88.0     No
    ## 3114 105.0  85.0     No
    ## 3115 158.0  88.0     No
    ## 3116 130.0  84.0     No
    ## 3117 198.0 106.0    Yes
    ## 3118 152.5 105.0     No
    ## 3119 153.0  85.0    Yes
    ## 3122 130.0  90.0     No
    ## 3123 107.5  75.0     No
    ## 3124 150.0  84.0     No
    ## 3125 160.0 105.0     No
    ## 3126  98.0  53.0     No
    ## 3127 122.0  83.0     No
    ## 3129 126.0  75.0     No
    ## 3130 121.0  85.5     No
    ## 3131 103.0  67.5     No
    ## 3132 135.0  88.0     No
    ## 3133 127.0  80.0     No
    ## 3134 125.0  86.0     No
    ## 3136 125.0  80.0     No
    ## 3137 150.0 105.0     No
    ## 3138  94.0  62.0     No
    ## 3139 163.5  87.0     No
    ## 3140 126.0  80.0     No
    ## 3142 132.0  86.0     No
    ## 3143 118.5  81.0     No
    ## 3144 105.0  65.0     No
    ## 3145 131.0  87.0     No
    ## 3147 118.0  80.0     No
    ## 3148 131.0  94.0     No
    ## 3149 106.0  64.0     No
    ## 3150 126.0  86.0     No
    ## 3151 122.0  88.5     No
    ## 3152 125.0  79.0     No
    ## 3153 147.0  89.0     No
    ## 3154 148.0 103.0     No
    ## 3155 130.0  86.0     No
    ## 3156 120.5  80.5     No
    ## 3157 137.0  79.0     No
    ## 3158 109.0  70.0     No
    ## 3159 122.0  82.0     No
    ## 3160 197.0  91.0     No
    ## 3161 147.5  97.0     No
    ## 3162 115.0  71.0     No
    ## 3163 131.0  84.0     No
    ## 3164 125.0  80.0     No
    ## 3165 186.5  99.0     No
    ## 3166 141.0  76.0     No
    ## 3167 119.0  80.0     No
    ## 3170 108.5  73.0     No
    ## 3171 112.5  77.5     No
    ## 3172 142.0  93.0     No
    ## 3173 131.0  81.0     No
    ## 3174 114.0  79.0     No
    ## 3175 121.0  77.0     No
    ## 3176 149.0  86.0     No
    ## 3177 124.5  80.0     No
    ## 3178 126.0  80.0     No
    ## 3179 123.0  78.0     No
    ## 3180 107.0  77.0     No
    ## 3183 147.0  98.0     No
    ## 3185 130.0  72.5     No
    ## 3186 153.0  89.0     No
    ## 3188 128.5  75.0     No
    ## 3189 127.5  83.0     No
    ## 3190 193.0 109.0    Yes
    ## 3192 106.0  65.0     No
    ## 3193 164.0 113.0     No
    ## 3194 154.0  99.0     No
    ## 3195 123.0  83.0     No
    ## 3196 131.0  92.0     No
    ## 3197 164.0 120.0    Yes
    ## 3198 148.5 100.0     No
    ## 3200 127.5  81.5     No
    ## 3201 136.0  75.0     No
    ## 3202 130.0  86.5     No
    ## 3203 116.0  64.0     No
    ## 3204 130.0  82.0     No
    ## 3205 114.0  75.0     No
    ## 3206 138.5  77.5     No
    ## 3207 123.0  72.0     No
    ## 3208 127.0  83.0     No
    ## 3209  98.0  67.0     No
    ## 3210 133.0  92.0     No
    ## 3211 196.0 101.0     No
    ## 3212 105.0  70.0     No
    ## 3213 139.0  96.0     No
    ## 3214 111.0  73.0     No
    ## 3215 117.5  90.0     No
    ## 3216 141.0  92.0     No
    ## 3217 158.0  80.0     No
    ## 3218 111.0  73.0     No
    ## 3219 127.5  92.0     No
    ## 3220 147.5  92.0     No
    ## 3221 131.5  82.5     No
    ## 3222 152.0  93.0     No
    ## 3223 127.0  69.0     No
    ## 3224 122.5  77.0     No
    ## 3225 159.0  90.0    Yes
    ## 3228 152.0  74.0     No
    ## 3229 133.0  89.0     No
    ## 3230 141.5 108.5     No
    ## 3231 172.0 111.0     No
    ## 3232 114.0  80.0     No
    ## 3233 139.0  83.0     No
    ## 3234 115.0  78.0     No
    ## 3236 127.0  72.0     No
    ## 3237 146.0  87.0     No
    ## 3238 150.5  95.0     No
    ## 3239 162.5  99.5     No
    ## 3240 115.0  71.0     No
    ## 3241 120.0  79.5     No
    ## 3243 132.0  88.0     No
    ## 3244 130.0  73.5     No
    ## 3245 121.0  82.0     No
    ## 3246 124.0  79.0     No
    ## 3249 163.5  97.0     No
    ## 3251 117.5  85.0     No
    ## 3252 132.0  86.0     No
    ## 3253 129.5  83.0     No
    ## 3254 163.0  97.0     No
    ## 3256 116.0  86.0     No
    ## 3257 114.5  70.0     No
    ## 3258 179.5  97.0     No
    ## 3259 125.0  87.0     No
    ## 3261 142.5  80.0     No
    ## 3262 126.0  80.0     No
    ## 3263 120.0  78.0     No
    ## 3264 126.0  76.0     No
    ## 3265 134.0  84.0     No
    ## 3266 118.5  74.5     No
    ## 3267 125.0  76.0     No
    ## 3268 160.0  90.0     No
    ## 3269 155.0 105.0     No
    ## 3270 133.0  82.0     No
    ## 3271 137.0  81.0     No
    ## 3272 144.0  91.0     No
    ## 3273 177.0 101.0     No
    ## 3274 136.0  84.0    Yes
    ## 3275 114.0  83.0     No
    ## 3276 138.0  99.0     No
    ## 3277 140.0  94.0     No
    ## 3278 116.0  79.0     No
    ## 3279 133.0  94.0     No
    ## 3280 140.0  78.0     No
    ## 3281 170.0  89.0    Yes
    ## 3282 124.0  72.5     No
    ## 3284 136.0  84.0     No
    ## 3287 151.0 106.5     No
    ## 3288 132.0  86.0     No
    ## 3289 107.0  73.0     No
    ## 3291 125.0  83.0     No
    ## 3292 131.0  87.0     No
    ## 3293 122.0  82.0     No
    ## 3294 123.5  77.0     No
    ## 3297 122.5  76.5     No
    ## 3298 139.0  74.0     No
    ## 3300 154.0  92.0     No
    ## 3301 103.0  71.0     No
    ## 3302 212.0 116.0     No
    ## 3303 123.0  73.0     No
    ## 3305 153.0 100.0     No
    ## 3306 140.0  87.0     No
    ## 3307 108.0  72.0     No
    ## 3308 167.0 109.0     No
    ## 3309 132.0  80.0     No
    ## 3310 115.0  81.0     No
    ## 3311 115.0  81.0     No
    ## 3312 122.0  76.5     No
    ## 3313 140.5  89.0     No
    ## 3314 107.0  73.5     No
    ## 3315 117.5  76.0     No
    ## 3316 136.5  81.0     No
    ## 3317 113.0  81.0     No
    ## 3318 162.0  90.0     No
    ## 3319 121.0  72.0     No
    ## 3320 131.0  93.0     No
    ## 3321 119.5  87.0     No
    ## 3322 128.0  84.0     No
    ## 3323 164.0  98.0     No
    ## 3324 152.0  89.0     No
    ## 3325 113.0  81.0     No
    ## 3326 117.5  72.0     No
    ## 3327 151.5  95.0     No
    ## 3328 120.0  81.0     No
    ## 3329 126.0  79.0     No
    ## 3330 148.0  90.0     No
    ## 3331 120.0  66.5     No
    ## 3334  95.5  70.0     No
    ## 3335 116.0  81.5     No
    ## 3336 105.0  69.0     No
    ## 3337 140.0  94.0     No
    ## 3338 124.0  80.0     No
    ## 3339 130.0  87.0     No
    ## 3340 132.0  80.0     No
    ## 3341 104.0  74.0     No
    ## 3342 142.0  91.0     No
    ## 3343 111.0  72.0     No
    ## 3344 120.0  80.0     No
    ## 3345 128.0  80.0     No
    ## 3346 144.0  92.0     No
    ## 3347 105.0  86.0     No
    ## 3348 140.0  88.0     No
    ## 3349 122.0  81.0     No
    ## 3350 115.0  83.0     No
    ## 3352 137.5  88.5     No
    ## 3353 127.5  82.5     No
    ## 3354 124.0  80.0     No
    ## 3355 170.0 107.0     No
    ## 3356 138.5  85.5     No
    ## 3357 154.0  91.0     No
    ## 3358 170.0 107.5     No
    ## 3359 199.5 107.0     No
    ## 3360 141.0  84.0     No
    ## 3361 161.0 103.0     No
    ## 3362 126.0  85.0     No
    ## 3363 111.0  71.0     No
    ## 3364 170.0  81.0     No
    ## 3366 152.0 102.0    Yes
    ## 3367 109.0  69.0     No
    ## 3368 128.0  82.0     No
    ## 3369 123.0  70.0     No
    ## 3370 110.0  68.0    Yes
    ## 3371 124.0  92.0     No
    ## 3372 152.0  92.0     No
    ## 3374 120.0  78.0     No
    ## 3375 108.0  73.0     No
    ## 3376 132.0  81.0     No
    ## 3377 136.5  84.0     No
    ## 3378 165.0  85.0     No
    ## 3379 128.5  73.0     No
    ## 3380 193.0  95.0     No
    ## 3382 119.0  80.0     No
    ## 3383 117.0  77.0     No
    ## 3384 115.0  80.0     No
    ## 3385 134.0  93.0     No
    ## 3386 131.0  85.0     No
    ## 3387 134.0  92.0     No
    ## 3388 132.5  97.5     No
    ## 3389 133.0  93.0     No
    ## 3390 129.0  84.0     No
    ## 3391 146.0  92.0     No
    ## 3392 135.0  80.0     No
    ## 3393 115.5  85.0     No
    ## 3394 122.5  73.0     No
    ## 3395 123.0  75.0     No
    ## 3396 131.0  80.0     No
    ## 3397 119.0  81.0     No
    ## 3398 129.0  81.0     No
    ## 3399 131.0  81.0     No
    ## 3400 114.0  78.0     No
    ## 3401 110.0  70.0     No
    ## 3402 168.0 100.0     No
    ## 3403 145.0  85.0     No
    ## 3404 143.0  84.0     No
    ## 3406 108.0  75.0     No
    ## 3407 121.0  82.0     No
    ## 3408 122.0  80.0     No
    ## 3409 117.0  73.0     No
    ## 3410 124.5  84.5     No
    ## 3411 127.0  81.0     No
    ## 3412 106.0  72.0     No
    ## 3413 166.0  88.0     No
    ## 3414 105.5  74.0     No
    ## 3415 196.0 103.0     No
    ## 3416 135.0  95.0     No
    ## 3417 113.5  77.0     No
    ## 3419 102.0  70.5     No
    ## 3420 118.5  73.0     No
    ## 3421 121.0  85.0     No
    ## 3422 117.0  78.0     No
    ## 3423 115.0  70.0     No
    ## 3424 116.5  87.0     No
    ## 3425  99.0  62.0     No
    ## 3426 115.0  79.0     No
    ## 3427 122.0  75.0     No
    ## 3428 151.0  74.0     No
    ## 3429 142.0  94.0     No
    ## 3430 112.0  66.0     No
    ## 3431 117.0  78.0     No
    ## 3433 123.0  85.0     No
    ## 3434 112.0  74.5     No
    ## 3435 112.5  77.5     No
    ## 3437 164.5 102.0     No
    ## 3438 130.0  80.0     No
    ## 3439 127.0  88.5     No
    ## 3440 142.5  85.0     No
    ## 3444 120.0  80.0     No
    ## 3445 137.0  82.0     No
    ## 3446 124.5  72.0     No
    ## 3447 112.5  60.0     No
    ## 3449 128.0  86.5     No
    ## 3450 119.0  65.0     No
    ## 3451 120.0  77.5     No
    ## 3452 136.5  87.0     No
    ## 3453 101.0  59.0     No
    ## 3454 129.0  85.0     No
    ## 3455 139.0  79.0     No
    ## 3456 127.0  76.5     No
    ## 3457 119.0  75.0     No
    ## 3458 114.0  80.0     No
    ## 3459 110.0  60.0     No
    ## 3461 127.5  76.0     No
    ## 3463 119.0  82.5     No
    ## 3464 134.5  87.0     No
    ## 3466 124.0  75.5     No
    ## 3467 159.0 100.0     No
    ## 3469 108.5  73.5     No
    ## 3470 119.0  62.5     No
    ## 3471 110.0  67.0     No
    ## 3472 181.0 112.5    Yes
    ## 3473 143.5  85.0     No
    ## 3474 177.5 120.0     No
    ## 3475 119.0  80.0     No
    ## 3476 125.0  75.0     No
    ## 3477 157.0  96.0     No
    ## 3478 102.0  74.5     No
    ## 3479 111.0  78.0     No
    ## 3480 128.0  77.5     No
    ## 3481 131.0  52.0     No
    ## 3482 122.0  68.0     No
    ## 3484 124.0  80.0     No
    ## 3486 125.0  75.0     No
    ## 3487 141.0  82.0     No
    ## 3488 141.0  87.0     No
    ## 3489 111.0  67.0     No
    ## 3490 107.5  70.0     No
    ## 3491 159.0 100.0     No
    ## 3492 143.0  87.0     No
    ## 3494 160.0  87.0    Yes
    ## 3495 154.5  93.0     No
    ## 3496 133.5  92.0     No
    ## 3497 112.0  83.0     No
    ## 3498 147.5  92.5     No
    ## 3499 123.0  76.0     No
    ## 3500 135.0  88.0     No
    ## 3501 113.0  77.0     No
    ## 3502 168.0 103.0     No
    ## 3503 117.5  82.5     No
    ## 3504 138.0  80.0     No
    ## 3505 195.0 110.0     No
    ## 3506 125.0  83.0     No
    ## 3507 159.0 102.0    Yes
    ## 3508 126.0  84.0     No
    ## 3509 140.0  94.0     No
    ## 3510 134.5  80.0     No
    ## 3511 115.0  78.0     No
    ## 3512 110.0  80.0     No
    ## 3513 144.5  91.5     No
    ## 3514 149.5  86.0     No
    ## 3515 142.0  90.0     No
    ## 3516 114.0  82.0     No
    ## 3517 150.0  74.0     No
    ## 3519 128.0  86.5     No
    ## 3520 165.0  84.0     No
    ## 3521 102.5  66.0     No
    ## 3522 112.0  75.5     No
    ## 3523 127.5  75.0     No
    ## 3524 143.5  93.0     No
    ## 3525 143.0  90.0     No
    ## 3526 121.5  81.5     No
    ## 3528 120.0  72.0     No
    ## 3530 141.0  92.0     No
    ## 3531 133.0  93.0     No
    ## 3532 102.0  73.0     No
    ## 3533 116.5  81.0     No
    ## 3535 158.0 105.0     No
    ## 3536 168.0  92.0     No
    ## 3538 139.0  84.0     No
    ## 3539 125.0  72.0     No
    ## 3541 112.0  70.0     No
    ## 3542 138.0  80.0     No
    ## 3543 110.0  65.0     No
    ## 3544 113.0  78.0     No
    ## 3545 182.0  86.0     No
    ## 3546 121.0  84.0     No
    ## 3547 130.0  87.0     No
    ## 3548 140.0  94.0     No
    ## 3549 122.5  87.0     No
    ## 3550 111.0  80.0     No
    ## 3551 160.5 106.5     No
    ## 3552 103.0  67.0     No
    ## 3553 118.5  71.0     No
    ## 3554 165.0 108.0     No
    ## 3555 140.0  89.5     No
    ## 3556 132.5  55.0    Yes
    ## 3557 123.0  76.0     No
    ## 3558 162.0  97.5     No
    ## 3559 132.0  88.0     No
    ## 3560 112.0  63.5     No
    ## 3561 150.5  87.0     No
    ## 3562 157.0  91.0     No
    ## 3564 153.0 106.0     No
    ## 3565 112.0  83.0     No
    ## 3566 113.0  82.0     No
    ## 3567 192.5 110.0     No
    ## 3568 108.0  70.5     No
    ## 3569 142.0  92.0     No
    ## 3570 153.0  87.0     No
    ## 3571 130.0  80.0     No
    ## 3572 108.0  65.0     No
    ## 3573 123.0  75.0     No
    ## 3574 167.0 100.0     No
    ## 3575 228.0 130.0    Yes
    ## 3576 113.0  70.0     No
    ## 3577 125.0  77.5     No
    ## 3578 129.0  83.0     No
    ## 3579 167.0  96.0     No
    ## 3580 151.0 102.0     No
    ## 3581 127.0  79.5     No
    ## 3582 128.0  90.0     No
    ## 3583 122.0  80.0     No
    ## 3584 119.0  66.0     No
    ## 3585 145.5  92.5     No
    ## 3586 168.5  97.0     No
    ## 3587 165.0  86.0     No
    ## 3588 133.0  85.5     No
    ## 3589 163.0  91.0     No
    ## 3590 114.0  72.0     No
    ## 3591 115.0  67.0     No
    ## 3592 132.0  89.5     No
    ## 3593 135.0  86.0     No
    ## 3594 101.5  67.0     No
    ## 3595 124.0  80.0     No
    ## 3596 122.0  70.0     No
    ## 3597 112.5  60.0     No
    ## 3598 118.0  80.0     No
    ## 3599 112.0  61.0     No
    ## 3600 148.0  89.0     No
    ## 3601 147.0  94.0     No
    ## 3602 146.0  92.0     No
    ## 3603 147.5  90.0     No
    ## 3604 136.0  90.0     No
    ## 3605 129.0  86.0     No
    ## 3606 174.0  97.0    Yes
    ## 3607 128.0  87.0     No
    ## 3608 108.0  77.0     No
    ## 3609 126.0  81.0     No
    ## 3610 159.5  93.5     No
    ## 3611 138.0  96.0     No
    ## 3614 125.0  72.0     No
    ## 3615 132.0  77.0     No
    ## 3616 120.0  72.5     No
    ## 3617 156.0 105.0     No
    ## 3618 142.0  66.0     No
    ## 3619 123.0  79.0     No
    ## 3620 135.0  85.0     No
    ## 3621 108.0  81.0     No
    ## 3623 128.0  83.5     No
    ## 3624 135.0  85.0     No
    ## 3625 118.0  68.0     No
    ## 3626 138.0  72.0     No
    ## 3628 163.0  89.0     No
    ## 3629 128.5  87.5     No
    ## 3630 112.0  78.0     No
    ## 3631 144.5  95.0     No
    ## 3632 158.0 108.0     No
    ## 3633 110.0  74.0     No
    ## 3634 115.0  80.0     No
    ## 3635 127.0  79.0     No
    ## 3636 127.5  75.5     No
    ## 3637 151.5 110.0     No
    ## 3639 159.0  90.0     No
    ## 3640 158.0  74.0     No
    ## 3641 123.0  92.0     No
    ## 3642 100.5  69.0     No
    ## 3643 166.0 107.0     No
    ## 3644 161.0 105.0    Yes
    ## 3646 111.0  79.0     No
    ## 3647  83.5  55.0     No
    ## 3649 248.0 130.0    Yes
    ## 3650 126.0  73.0     No
    ## 3651 121.0  79.0     No
    ## 3652 111.0  72.5     No
    ## 3653 152.0  70.0     No
    ## 3654 155.0  90.0     No
    ## 3655 108.0  70.0     No
    ## 3656 172.5 112.5     No
    ## 3657 122.0  81.0     No
    ## 3658 163.0  94.0     No
    ## 3659 123.0  81.0     No
    ## 3660 108.0  73.0     No
    ## 3661 112.5  62.5     No
    ## 3662 107.5  71.0     No
    ## 3663 143.0  82.5     No
    ## 3664 129.0  87.0     No
    ## 3665  98.0  60.0     No
    ## 3666 116.0  76.0     No
    ## 3667 114.0  72.0     No
    ## 3668 124.0  89.0     No
    ## 3669 196.0 109.0     No
    ## 3670 160.0  98.0     No
    ## 3671 122.5  85.0     No
    ## 3672 156.0  95.0     No
    ## 3673 126.5  75.5     No
    ## 3675 108.0  75.0     No
    ## 3676 122.0  81.0     No
    ## 3677 176.5  92.0    Yes
    ## 3678 115.0  65.0     No
    ## 3679 121.0  69.0     No
    ## 3680 131.0  79.0     No
    ## 3681 127.0  77.0     No
    ## 3682 161.5  96.0     No
    ## 3683 100.5  66.0     No
    ## 3684 131.5  83.0     No
    ## 3685 131.5  77.0     No
    ## 3687 135.5  77.0     No
    ## 3688 125.0  76.0     No
    ## 3689 131.0  83.0     No
    ## 3690 143.5  77.5     No
    ## 3691 136.0  90.0     No
    ## 3692 123.5  78.0     No
    ## 3693 164.0 119.0     No
    ## 3694 145.0  87.0     No
    ## 3696 133.0  82.0     No
    ## 3697 139.0  81.0     No
    ## 3698 124.0  84.0     No
    ## 3699 139.0  81.5     No
    ## 3700 113.0  87.0     No
    ## 3701 116.0  77.0     No
    ## 3702 182.0 110.0    Yes
    ## 3703 113.0  74.0     No
    ## 3704 126.0  71.0     No
    ## 3705 112.0  67.0     No
    ## 3706 130.0  74.0     No
    ## 3707 104.0  76.0     No
    ## 3708 106.0  77.0     No
    ## 3709 177.5 110.0     No
    ## 3710 121.0  78.0     No
    ## 3711 134.0  78.0     No
    ## 3712 110.0  77.5     No
    ## 3713 128.5  80.0     No
    ## 3714 115.0  72.0     No
    ## 3715 154.0  96.0     No
    ## 3716 202.0 124.0    Yes
    ## 3717 151.0  85.0     No
    ## 3718 154.0  87.0     No
    ## 3719 162.0 109.0     No
    ## 3721 117.0  78.5     No
    ## 3723 123.0  82.0     No
    ## 3724 110.0  83.0     No
    ## 3725 117.0  72.0     No
    ## 3726 122.0  85.0     No
    ## 3727 154.0  80.0     No
    ## 3728 106.0  65.0     No
    ## 3729 114.5  77.0     No
    ## 3730  99.0  62.0     No
    ## 3731  99.5  66.0     No
    ## 3732 185.0  95.0     No
    ## 3734 108.0  70.0     No
    ## 3736 121.0  82.0    Yes
    ## 3737 116.0  67.0     No
    ## 3739 158.0  89.0     No
    ## 3740  85.0  70.0     No
    ## 3742 169.0 111.0     No
    ## 3743 116.0  66.0     No
    ## 3744 130.0  84.5     No
    ## 3745 117.0  86.0     No
    ## 3746 105.5  67.0     No
    ## 3747 135.0  88.0     No
    ## 3748 146.0  89.0     No
    ## 3749 165.0 100.0     No
    ## 3750 143.0  96.0     No
    ## 3752 176.0 109.0     No
    ## 3753 118.0  74.0     No
    ## 3756 136.0  94.0     No
    ## 3757 117.0  81.5     No
    ## 3758 130.0  94.0     No
    ## 3759 135.0  84.0     No
    ## 3760 122.5  75.0     No
    ## 3762 135.0  97.0     No
    ## 3763 107.0  69.0     No
    ## 3764 167.5  92.5     No
    ## 3765 122.5  80.0     No
    ## 3766 133.0  72.0     No
    ## 3767 138.0  87.0     No
    ## 3768 151.0  91.0     No
    ## 3770 101.0  75.0     No
    ## 3771 109.0  78.5     No
    ## 3772 132.0  78.0     No
    ## 3774 123.0  66.0     No
    ## 3775 119.0  67.0     No
    ## 3776 137.0  79.0     No
    ## 3777 130.0  89.0     No
    ## 3778 128.0  78.0     No
    ## 3780 161.0  85.0     No
    ## 3781 120.0  75.0     No
    ## 3782 230.0 110.0     No
    ## 3783 118.0  86.5     No
    ## 3785 143.0  75.0     No
    ## 3786 140.0  71.0     No
    ## 3787 155.0  79.0     No
    ## 3788 153.0 102.5     No
    ## 3789 146.5  79.0     No
    ## 3791 137.5  82.5     No
    ## 3792 138.0  72.0     No
    ## 3793 147.5  95.0     No
    ## 3794 129.0  80.0     No
    ## 3795 101.0  66.0     No
    ## 3797 162.0  85.0     No
    ## 3798 197.0  72.0     No
    ## 3799 133.0  89.0     No
    ## 3800 164.0  89.0     No
    ## 3801 189.0 121.0     No
    ## 3802 132.0  74.0     No
    ## 3803 125.5  82.0     No
    ## 3804 134.0  93.0     No
    ## 3805 118.0  92.0     No
    ## 3806 142.5  85.0     No
    ## 3807 124.0  82.0     No
    ## 3808 126.0  52.0     No
    ## 3809 136.0  95.5     No
    ## 3810 112.0  78.0     No
    ## 3812 106.0  48.0     No
    ## 3813 118.0  71.0     No
    ## 3815 147.0 100.0     No
    ## 3817 128.0  78.0     No
    ## 3818 142.0  84.0     No
    ## 3819 140.0  92.5     No
    ## 3820 113.0  68.0     No
    ## 3821 151.0  89.0     No
    ## 3823 144.0  88.0     No
    ## 3824 153.0 102.5     No
    ## 3825 146.0  86.0     No
    ## 3826 141.0  93.0     No
    ## 3829 112.5  71.0     No
    ## 3830 113.0  79.0     No
    ## 3831 100.0  76.0     No
    ## 3832 140.0  90.0     No
    ## 3833 114.0  82.0     No
    ## 3834 130.0  77.5     No
    ## 3835 155.0  92.5     No
    ## 3836 128.0  82.0     No
    ## 3838  97.0  63.0     No
    ## 3839 112.0  73.0     No
    ## 3841 167.0  92.0    Yes
    ## 3842 113.0  73.0     No
    ## 3844 166.0  90.0    Yes
    ## 3845 214.0  94.0    Yes
    ## 3846 143.0  88.0     No
    ## 3847 134.5  90.5     No
    ## 3848 132.0  91.0     No
    ## 3849 124.0  80.0     No
    ## 3850 159.5  82.5     No
    ## 3852 125.0  75.0     No
    ## 3853 123.0  73.0     No
    ## 3854 117.5  75.0     No
    ## 3856 119.0  75.0     No
    ## 3857 131.0  81.0     No
    ## 3858 123.0  74.0     No
    ## 3859 117.5  72.5     No
    ## 3860 140.0  84.0     No
    ## 3861 116.0  70.0     No
    ## 3862 114.0  75.0     No
    ## 3863 119.0  65.0     No
    ## 3864 125.0  79.0     No
    ## 3865 147.5  88.0     No
    ## 3866 124.5  86.5     No
    ## 3867 111.0  56.0     No
    ## 3868 141.0  84.0     No
    ## 3869 105.0  77.0     No
    ## 3870 127.0  83.0     No
    ## 3871 126.5  79.0     No
    ## 3872 136.5  78.0     No
    ## 3873 145.0  95.0     No
    ## 3874 134.0  86.5     No
    ## 3875 135.0  86.0     No
    ## 3876 145.0  92.0     No
    ## 3877 138.0  94.0     No
    ## 3878 149.0  96.0     No
    ## 3879 170.0 118.0     No
    ## 3880 141.0  90.0     No
    ## 3881 141.0  92.0     No
    ## 3884 142.5  82.0     No
    ## 3885 120.0  79.0     No
    ## 3886 142.0  88.0    Yes
    ## 3887 172.5  85.0     No
    ## 3888 128.0  81.0     No
    ## 3889 109.5  69.0     No
    ## 3890 104.0  78.0     No
    ## 3891 144.5  88.0     No
    ## 3892 118.0  73.0     No
    ## 3893 150.0  84.0     No
    ## 3894 158.0 102.0     No
    ## 3895 148.0  91.5     No
    ## 3896 120.0  84.0     No
    ## 3897 155.0  90.0     No
    ## 3898 108.0  70.0     No
    ## 3899 196.0 102.0     No
    ## 3900 115.0  85.0     No
    ## 3901 116.0  83.0     No
    ## 3902 167.5 110.0     No
    ## 3903 143.0  93.0     No
    ## 3904  96.5  62.5     No
    ## 3905 133.5  85.5     No
    ## 3906 132.0  81.0     No
    ## 3907 124.0  84.0     No
    ## 3909 146.5  81.0     No
    ## 3910 153.0 100.0     No
    ## 3911 122.0  78.0     No
    ## 3912 140.0  95.0     No
    ## 3913 143.0  92.0     No
    ## 3915 145.5  94.0     No
    ## 3916 163.0 112.5     No
    ## 3918 122.0  86.0     No
    ## 3919 149.0  85.0     No
    ## 3920 134.0  74.0     No
    ## 3921 110.0  67.5     No
    ## 3922 163.0 101.0     No
    ## 3924 113.0  75.5     No
    ## 3925 115.0  77.0     No
    ## 3926 100.5  66.0     No
    ## 3928 158.0  97.0     No
    ## 3929 118.0  71.0     No
    ## 3930 110.0  76.0     No
    ## 3931 122.5  84.0     No
    ## 3932 130.0  71.5     No
    ## 3933 160.0  92.0     No
    ## 3934 118.0  75.5     No
    ## 3937  99.0  59.0     No
    ## 3938 125.0  88.0     No
    ## 3939 106.0  73.0     No
    ## 3940 115.0  83.0     No
    ## 3941 126.0  82.0     No
    ## 3942 120.0  89.0     No
    ## 3944 102.0  64.5     No
    ## 3946 131.0  87.0     No
    ## 3948 113.0  83.0     No
    ## 3949 141.0  81.0     No
    ## 3950 139.0  85.0     No
    ## 3951 101.0  59.0     No
    ## 3953 100.5  62.0     No
    ## 3954 133.0  90.0     No
    ## 3955 131.0  81.0     No
    ## 3956 129.0  81.0     No
    ## 3957 115.0  90.0     No
    ## 3958 137.0  86.0     No
    ## 3960 119.0  78.0     No
    ## 3961 154.0  92.0     No
    ## 3962 136.0  87.0     No
    ## 3964 102.0  68.0     No
    ## 3965 164.0  97.0    Yes
    ## 3966 170.0 107.0    Yes
    ## 3967 152.5  90.0     No
    ## 3968 116.0  71.0     No
    ## 3969 113.0  68.0     No
    ## 3970 127.5  82.5     No
    ## 3971 141.5  91.0     No
    ## 3972 110.0  84.0     No
    ## 3973 125.0  80.0     No
    ## 3974 115.0  75.0     No
    ## 3975 137.0  87.0     No
    ## 3976 173.0  89.0     No
    ## 3977 162.5 104.0     No
    ## 3978 105.0  67.5     No
    ## 3979 100.0  70.0     No
    ## 3980 146.0  86.0     No
    ## 3981 148.5  94.0     No
    ## 3982 135.0  86.5     No
    ## 3983 128.0  82.0     No
    ## 3984 119.0  77.0     No
    ## 3985 129.0  93.0     No
    ## 3986 109.0  74.0     No
    ## 3987 115.5  82.5     No
    ## 3989 141.0  93.0     No
    ## 3990 129.0  94.0     No
    ## 3991 170.0 113.0     No
    ## 3992 121.0  81.0     No
    ## 3993 124.0  78.0     No
    ## 3994 100.0  72.0     No
    ## 3995 158.0  90.0     No
    ## 3996 109.0  61.0     No
    ## 3997  99.0  59.0     No
    ## 3998 120.5  80.0     No
    ## 3999 133.0  87.0     No
    ## 4000 130.0  71.0     No
    ## 4001 132.5  87.0     No
    ## 4002 156.5  93.0     No
    ## 4003 131.0  81.0     No
    ## 4004 108.0  70.0     No
    ## 4005 127.0  66.5     No
    ## 4006 111.0  72.0     No
    ## 4007 124.0  69.0     No
    ## 4008 167.0 101.0     No
    ## 4009 111.0  72.0     No
    ## 4010 168.5  98.5    Yes
    ## 4011 155.5  98.0     No
    ## 4012 107.0  74.0     No
    ## 4013 111.0  79.0     No
    ## 4014 123.5  80.5     No
    ## 4015 120.0  80.0     No
    ## 4016 215.0 110.0    Yes
    ## 4017 110.0  70.0     No
    ## 4018 192.5 113.0    Yes
    ## 4019 135.0  82.5     No
    ## 4020 111.0  74.0     No
    ## 4021 106.0  60.0     No
    ## 4022 149.0  91.0     No
    ## 4023 110.0  70.0     No
    ## 4024 137.0  75.0     No
    ## 4025 172.5 100.5     No
    ## 4026 130.0  62.0     No
    ## 4027 130.0  80.0     No
    ## 4028 172.0 108.0     No
    ## 4029 124.0  76.0     No
    ## 4030 102.0  56.0     No
    ## 4031 104.5  65.0     No
    ## 4032 133.0  86.5     No
    ## 4033 120.5  76.0     No
    ## 4034 151.0  90.0     No
    ## 4035 169.0  82.0    Yes
    ## 4036 111.0  80.0     No
    ## 4038 157.0  98.0     No
    ## 4039 140.0  86.5    Yes
    ## 4040 147.5  87.5     No
    ## 4041 123.0  88.0     No
    ## 4043 172.0  87.0     No
    ## 4044 148.0  89.5     No
    ## 4045 144.0  96.0     No
    ## 4048 110.0  61.0     No
    ## 4049 147.0  65.0     No
    ## 4050 160.0  85.0     No
    ## 4051 140.0  99.0     No
    ## 4053 120.0  77.0     No
    ## 4054 111.5  67.0     No
    ## 4055 120.0  80.0     No
    ## 4056 117.0  73.0     No
    ## 4058 152.5 108.0     No
    ## 4059 112.0  71.0     No
    ## 4060 106.0  73.0     No
    ## 4061 127.0  85.0     No
    ## 4062 128.0  81.0     No
    ## 4063 165.0  95.0     No
    ## 4064 105.0  70.0     No
    ## 4065 164.0 100.0     No
    ## 4066 119.0  73.0     No
    ## 4067 143.5  92.0     No
    ## 4069 121.0  81.0     No
    ## 4070 110.0  68.0     No
    ## 4071 139.0  75.0     No
    ## 4072 156.0  90.0     No
    ## 4073 110.0  65.0     No
    ## 4074 128.0  88.0     No
    ## 4075 128.0  79.0     No
    ## 4076 150.0  98.0     No
    ## 4077 116.5  74.0     No
    ## 4079 155.0  80.5     No
    ## 4080 119.0  75.5     No
    ## 4081 128.0  65.0     No
    ## 4083 173.0  88.0     No
    ## 4084 130.0  94.0     No
    ## 4085 137.0  79.0     No
    ## 4086 112.0  65.0     No
    ## 4087 187.0  97.0     No
    ## 4088 116.0  72.0     No
    ## 4089 139.0  85.0     No
    ## 4090 141.0  89.0     No
    ## 4091 127.5  81.0     No
    ## 4093 173.0  98.5     No
    ## 4094 129.0  97.0     No
    ## 4095 138.5  85.0     No
    ## 4096 130.0  86.0     No
    ## 4097  92.5  66.5     No
    ## 4098 143.0  96.0     No
    ## 4099 114.0  80.0     No
    ## 4100 110.0  76.0     No
    ## 4101 131.5  77.5     No
    ## 4102 151.5  71.0     No
    ## 4103 157.0 112.0     No
    ## 4105 131.5  84.0     No
    ## 4106 194.0 111.0    Yes
    ## 4107 170.0 100.0     No
    ## 4108 130.5  82.0     No
    ## 4109 153.0  80.0     No
    ## 4110 132.0  95.0     No
    ## 4111 150.0  99.0     No
    ## 4112 133.0  86.0     No
    ## 4113 130.5  63.0     No
    ## 4114 118.0  79.0     No
    ## 4115 144.0  80.0     No
    ## 4116 121.5  74.5     No
    ## 4117 150.0  91.0     No
    ## 4118 116.0  84.0     No
    ## 4119 126.5  81.0     No
    ## 4120  93.5  58.0     No
    ## 4121 126.5  76.5     No
    ## 4123 152.0  96.5     No
    ## 4124 132.0  88.0     No
    ## 4125 152.5  92.5     No
    ## 4127 111.0  79.0     No
    ## 4128 124.5  75.5     No
    ## 4129 111.5  72.0     No
    ## 4130 126.0  86.0     No
    ## 4131 137.0  83.0    Yes
    ## 4132 174.0 130.0     No
    ## 4133 173.0  67.0     No
    ## 4134 110.0  68.0     No
    ## 4135 150.0  82.0     No
    ## 4136 118.0  81.0     No
    ## 4138 177.5 100.0     No
    ## 4139 117.0  76.0     No
    ## 4142 161.0  90.0     No
    ## 4143 125.0  87.0     No
    ## 4144 125.0  87.0     No
    ## 4145 144.0  87.0     No
    ## 4146 102.5  63.5     No
    ## 4147 129.5 100.0     No
    ## 4148 120.5  80.0     No
    ## 4149  96.0  66.0     No
    ## 4150 112.5  76.5     No
    ## 4151 176.0  78.0     No
    ## 4152 120.0  87.0     No
    ## 4153 146.0  82.5     No
    ## 4154 147.0  90.0     No
    ## 4155 126.0  79.0     No
    ## 4156 169.0  91.0     No
    ## 4157 157.0  94.0     No
    ## 4158 110.0  80.0     No
    ## 4160 125.0  80.0     No
    ## 4161 207.0 122.5     No
    ## 4162 136.5  85.0     No
    ## 4163 145.0  99.0     No
    ## 4164 134.0  88.0     No
    ## 4165 113.0  68.0     No
    ## 4166 100.0  66.5     No
    ## 4167 129.0  77.0     No
    ## 4171 131.0  76.0     No
    ## 4172 129.0  84.0     No
    ## 4173 101.0  68.5     No
    ## 4174 122.0  78.0     No
    ## 4175 153.0  80.0     No
    ## 4176 118.0  87.0     No
    ## 4177 119.0  86.0     No
    ## 4178 132.5  82.5     No
    ## 4179 111.0  60.5     No
    ## 4180 133.0  84.0     No
    ## 4181 110.0  72.0     No
    ## 4182 122.5  75.0     No
    ## 4183 165.0  95.0     No
    ## 4184 122.5  77.5     No
    ## 4186 126.0  67.0     No
    ## 4187 142.0  87.0     No
    ## 4188 115.0  74.0     No
    ## 4189 135.0 100.0     No
    ## 4191 141.0  76.0     No
    ## 4192 118.0  73.0     No
    ## 4193 123.0  82.0     No
    ## 4195 118.0  74.0     No
    ## 4196 140.0  88.0     No
    ## 4197 108.0  72.0     No
    ## 4198 129.0  85.0     No
    ## 4200 125.0  71.0     No
    ## 4201 140.0  92.0     No
    ## 4202 107.0  70.0     No
    ## 4203 136.0  73.0     No
    ## 4204 157.0  94.0     No
    ## 4205 147.0 101.0     No
    ## 4206 133.0  90.0     No
    ## 4207 162.5  92.5     No
    ## 4208 130.0  85.0     No
    ## 4209 143.0  87.0     No
    ## 4211 173.0 106.0    Yes
    ## 4212 108.5  68.0     No
    ## 4213 130.0  86.5     No
    ## 4214 112.0  67.0     No
    ## 4215 112.5  80.0     No
    ## 4216 130.0  80.0     No
    ## 4217 159.0 110.0     No
    ## 4218 185.5 115.5     No
    ## 4219 129.0  83.0     No
    ## 4220 103.0  62.0     No
    ## 4221 126.0  75.0     No
    ## 4222 213.0 133.0     No
    ## 4223 160.0 105.0     No
    ## 4224 175.0 101.0     No
    ## 4225 111.5  74.0     No
    ## 4226 129.5  93.5     No
    ## 4227 124.0  77.0     No
    ## 4228 110.0  73.0     No
    ## 4229 118.0  79.5     No
    ## 4230 148.0  95.0     No
    ## 4232  97.5  62.5     No
    ## 4233 124.0  89.0     No
    ## 4234 104.0  64.0     No
    ## 4235 107.0  69.0     No
    ## 4236 127.0  81.0     No
    ## 4237 115.0  80.0     No
    ## 4238 130.0  80.0     No
    ## 4239 157.0  80.0    Yes
    ## 4240 192.5 115.0     No
    ## 4241 135.0  80.0     No
    ## 4242 176.0 109.0     No
    ## 4243 121.0  82.0     No
    ## 4245 112.0  72.0     No
    ## 4246 142.5  95.0     No
    ## 4247 107.0  81.0     No
    ## 4248 116.0  71.5     No
    ## 4249 117.0  78.5     No
    ## 4250 130.0  86.0     No
    ## 4251 133.5  97.5     No
    ## 4252 135.0  83.0     No
    ## 4253 158.0 103.0     No
    ## 4255 205.5 123.5     No
    ## 4257 140.0  93.0     No
    ## 4258 109.0  79.0     No
    ## 4259 118.0  79.0     No
    ## 4260 136.0  80.5     No
    ## 4261 192.5 125.0     No
    ## 4262 131.0  84.0     No
    ## 4263 200.0 125.0    Yes
    ## 4264 115.0  72.5     No
    ## 4265 120.0  81.0     No
    ## 4266 104.0  69.0     No
    ## 4267 124.0  70.0     No
    ## 4268 116.0  76.0     No
    ## 4269 128.0  81.0     No
    ## 4271 115.0  71.0     No
    ## 4272 187.0  95.0     No
    ## 4273 121.0  81.0     No
    ## 4275 126.0  74.0     No
    ## 4276 143.0  88.0     No
    ## 4277 149.5  84.0     No
    ## 4278 151.0  87.5     No
    ## 4279 116.0  74.0     No
    ## 4280  93.0  62.5     No
    ## 4281 125.0  80.0     No
    ## 4282 121.0  79.0     No
    ## 4283 140.5  83.0     No
    ## 4284 136.0  80.0     No
    ## 4285 152.0  82.0     No
    ## 4287 122.0  67.5     No
    ## 4288 123.0  77.0     No
    ## 4289 133.0  97.0     No
    ## 4290 134.0  75.0     No
    ## 4291 130.0  84.0     No
    ## 4292 123.0  80.0     No
    ## 4293 142.0  86.0     No
    ## 4294 110.0  74.5     No
    ## 4295 176.0  98.0    Yes
    ## 4296 104.0  79.0     No
    ## 4297 158.0 102.5    Yes
    ## 4298 146.5  92.0     No
    ## 4299 139.0  74.0     No
    ## 4300 141.0  80.0     No
    ## 4301 108.0  68.0     No
    ## 4302 142.0  68.0     No
    ## 4303 115.0  73.0     No
    ## 4304 190.0  88.0     No
    ## 4305 149.0  83.0     No
    ## 4306 146.5  97.5     No
    ## 4307 137.0  82.0     No
    ## 4308 122.0  72.0     No
    ## 4309 140.0  77.0    Yes
    ## 4310 127.0  81.0     No
    ## 4313 206.0 116.0     No
    ## 4314 108.5  70.5     No
    ## 4315 164.0 102.0     No
    ## 4316 124.0  84.0     No
    ## 4317 155.0  80.0     No
    ## 4318 114.0  80.0     No
    ## 4319 147.5  87.5     No
    ## 4320 125.0  82.0     No
    ## 4321 120.0  77.5     No
    ## 4322 170.0 101.0     No
    ## 4323  97.5  60.0     No
    ## 4324 134.0  96.0     No
    ## 4325 115.0  78.0     No
    ## 4326 128.0  86.0     No
    ## 4327 118.0  80.0     No
    ## 4330 126.0  73.0     No
    ## 4331 135.0  82.0     No
    ## 4332 118.0  81.0     No
    ## 4333 162.0 110.0     No
    ## 4334 105.0  82.0     No
    ## 4335 131.5  84.0     No
    ## 4336 144.0  96.5     No
    ## 4337 103.0  70.0     No
    ## 4338 128.0  74.0     No
    ## 4339 108.0  78.0     No
    ## 4340 136.0  93.0     No
    ## 4341 127.5  75.0     No
    ## 4342 103.0  72.5     No
    ## 4343 147.0  85.0     No
    ## 4345 157.5  83.0     No
    ## 4346 126.0  84.0     No
    ## 4347  99.0  62.0     No
    ## 4348 125.5  85.5     No
    ## 4349 143.0  81.0    Yes
    ## 4350 135.0  97.5     No
    ## 4352 109.0  72.0     No
    ## 4354 118.5  77.5     No
    ## 4355 136.0  86.0     No
    ## 4357 130.0  81.5     No
    ## 4358 133.0  78.0     No
    ## 4359 106.0  58.0     No
    ## 4360 180.0 108.0    Yes
    ## 4361 101.0  70.0     No
    ## 4363 125.0  87.0     No
    ## 4364 154.5  88.0     No
    ## 4365 210.0 127.5     No
    ## 4366 129.5  80.0     No
    ## 4367 135.0  94.0     No
    ## 4368 107.0  75.0     No
    ## 4369 111.0  70.0     No
    ## 4370 121.5  57.0     No
    ## 4371 132.5  85.0     No
    ## 4372 137.5  89.0     No
    ## 4373 129.0  95.0     No
    ## 4374 121.0  83.0     No
    ## 4375 157.0  86.0     No
    ## 4376 182.0  99.0     No
    ## 4377 118.5  70.0     No
    ## 4378 137.5  87.5     No
    ## 4380 106.0  63.0     No
    ## 4381 143.0  79.0     No
    ## 4382 110.0  70.0     No
    ## 4383 117.0  74.0     No
    ## 4384 150.0  89.0     No
    ## 4385 101.0  69.0     No
    ## 4386 129.0  80.0     No
    ## 4387 195.0 105.0     No
    ## 4388 179.0  96.0    Yes
    ## 4389 124.0  78.0     No
    ## 4390 116.0  80.0     No
    ## 4391 170.0 105.0     No
    ## 4392 157.5 105.0     No
    ## 4393 133.0  83.0     No
    ## 4394 115.0  60.0     No
    ## 4395 130.0  80.0     No
    ## 4396 160.5 100.0     No
    ## 4397 146.0  84.0     No
    ## 4398 142.0  72.0     No
    ## 4399 136.0  77.0     No
    ## 4400 103.0  67.5     No
    ## 4401 124.0  76.5     No
    ## 4403 135.0  80.0     No
    ## 4404 126.5  88.0     No
    ## 4405 105.0  70.0     No
    ## 4406 144.0  88.0     No
    ## 4407 141.0  95.0     No
    ## 4408 123.0  78.5     No
    ## 4409 155.0  82.0     No
    ## 4410 125.0  80.0     No
    ## 4411 167.0  94.0     No
    ## 4412 137.5  84.5     No
    ## 4413 125.0  84.5     No
    ## 4414 128.0  82.0     No
    ## 4415 119.0  74.0     No
    ## 4416 188.0 110.0     No
    ## 4417 149.0  98.0    Yes
    ## 4418 120.0  80.0     No
    ## 4419 137.5  85.0     No
    ## 4420 122.0  84.0     No
    ## 4421 125.5  84.0     No
    ## 4422 129.5  88.0     No
    ## 4423 190.0 130.0     No
    ## 4426 141.0  81.0     No
    ## 4427 168.0  97.0     No
    ## 4428 179.0  92.0     No
    ## 4429 126.5  80.0     No
    ## 4432 133.5  83.0     No
    ## 4433 141.0  98.0     No
    ## 4434 133.0  86.0     No

You can also perform negative indexing, using the `-` operator.

``` r
fhsCOV <- select(fhs, -MI)
fhsCOV
```

    ##         SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1      Male     195  39 106.0  70.0       No 26.97       No     No      77
    ## 2    Female     250  46 121.0  81.0       No 28.73       No     No      76
    ## 3      Male     245  48 127.5  80.0      Yes 25.34       No     No      70
    ## 4    Female     225  61 150.0  95.0      Yes 28.58       No     No     103
    ## 5    Female     285  46 130.0  84.0      Yes 23.10       No     No      85
    ## 6    Female     228  43 180.0 110.0       No 30.30       No     No      99
    ## 7    Female     205  63 138.0  71.0       No 33.11       No     No      85
    ## 8    Female     313  45 100.0  71.0      Yes 21.68       No     No      78
    ## 9      Male     260  52 141.5  89.0       No 26.36       No     No      79
    ## 10     Male     225  43 162.0 107.0      Yes 23.61       No     No      88
    ## 11   Female     254  50 133.0  76.0       No 22.91       No     No      76
    ## 12   Female     247  43 131.0  88.0       No 27.64       No     No      61
    ## 13     Male     294  46 142.0  94.0      Yes 26.31       No     No      64
    ## 14   Female     332  41 124.0  88.0       No 31.31       No    Yes      84
    ## 16   Female     221  38 140.0  90.0      Yes 21.35       No     No      70
    ## 17     Male     232  48 138.0  90.0      Yes 22.37       No     No      72
    ## 18   Female     291  46 112.0  78.0      Yes 23.38       No     No      89
    ## 19   Female     195  38 122.0  84.5      Yes 23.24       No     No      78
    ## 20     Male     195  41 139.0  88.0       No 26.88       No     No      65
    ## 21   Female     190  42 108.0  70.5      Yes 21.59       No     No      85
    ## 23   Female     234  52 148.0  78.0       No 34.17       No     No     113
    ## 24   Female     215  52 132.0  82.0      Yes 25.11       No     No      75
    ## 25     Male     270  44 137.5  90.0      Yes 21.96       No     No      83
    ## 26     Male     294  47 102.0  68.0      Yes 24.18       No     No      66
    ## 28     Male     225  35 132.0  91.0      Yes 26.09       No     No      83
    ## 29   Female     272  61 182.0 121.0       No 32.80       No     No      65
    ## 30   Female     247  60 130.0  88.0       No 30.36       No     No      74
    ## 31     Male     295  36 102.0  68.0      Yes 28.15       No     No      63
    ## 32     Male     226  43 115.0  85.5      Yes 27.57       No     No      75
    ## 33     Male     227  40 125.0  90.0       No 23.62       No     No      79
    ## 34   Female     209  59 150.0  85.0       No 20.77       No     No      88
    ## 36     Male     235  62 125.0  74.5       No 25.90       No     No      73
    ## 37     Male     214  54 147.0  74.0      Yes 24.71       No     No      87
    ## 38     Male     225  37 124.5  92.5       No 38.53       No     No      83
    ## 40     Male     178  52 160.0  98.0       No 40.11      Yes     No     225
    ## 41   Female     233  42 153.0 101.0      Yes 28.93       No     No      90
    ## 42     Male     180  36 111.0  73.0       No 27.78       No     No      80
    ## 43   Female     243  43 116.5  80.0      Yes 26.87       No     No      78
    ## 44   Female     237  41 122.0  78.0      Yes 23.28       No     No      74
    ## 46     Male     195  54 132.0  83.5       No 26.21       No     No     100
    ## 47   Female     311  53 206.0  92.0       No 21.51      Yes    Yes     215
    ## 48   Female     208  49  96.0  63.0       No 20.68       No     No      98
    ## 49   Female     252  65 179.5 114.0       No 30.47       No     No      87
    ## 50     Male     261  46 119.0  77.5      Yes 23.59       No     No      74
    ## 51   Female     179  63 116.0  69.0      Yes 22.15       No     No      75
    ## 53   Female     267  63 156.5  92.5      Yes 27.10       No     No      79
    ## 54     Male     216  51 112.0  66.0       No 23.47       No     No      95
    ## 55   Female     237  47 130.0  78.0      Yes 19.66       No     No      75
    ## 56   Female     240  62 145.0  82.5       No 28.27       No     No      75
    ## 58   Female     250  46 116.0  71.0      Yes 20.35       No     No      94
    ## 59   Female     266  54 114.0  76.0      Yes 17.61      Yes     No      55
    ## 60     Male     255  49 143.5  81.0      Yes 25.65       No     No      80
    ## 61     Male     185  44 115.0  69.0       No 22.29       No     No      82
    ## 62   Female     205  40 158.0 102.0      Yes 25.45       No     No      87
    ## 63     Male     270  56 121.0  79.0      Yes 23.58       No     No      93
    ## 64   Female     254  67 157.0  89.0       No 24.25       No     No      74
    ## 65     Male     220  53 123.5  75.0      Yes 19.64       No     No      73
    ## 66   Female     235  57 126.5  80.0      Yes 24.88       No     No      72
    ## 67     Male     220  57 136.0  84.0       No 26.84       No     No      64
    ## 68   Female     252  63 154.0  87.0       No 28.60       No     No      45
    ## 69   Female     212  62 190.0  99.0       No 29.64      Yes     No     202
    ## 70     Male     223  38 107.0  73.0      Yes 23.01       No     No      78
    ## 71   Female     300  47 112.5  60.0      Yes 20.13       No     No      83
    ## 72   Female     302  52 110.0  67.5       No 23.51       No     No      87
    ## 74     Male     175  42 116.0  72.5      Yes 28.61       No     No      95
    ## 76   Female     189  41 150.0 106.0      Yes 33.80       No     No      75
    ## 77   Female     221  44 110.0  76.0      Yes 22.16       No     No      83
    ## 78   Female     258  59 138.5  85.0       No 34.55       No     No     103
    ## 79   Female     202  44 155.0  85.0       No 24.04       No     No      68
    ## 81   Female     183  45 151.0 101.0       No 45.80       No     No      63
    ## 82     Male     274  41 152.0  90.0      Yes 30.58       No     No      65
    ## 83     Male     170  60 179.0  94.0      Yes 26.52       No     No      83
    ## 84     Male     285  39 155.0 110.0       No 32.51       No     No      70
    ## 85   Female     210  53 138.0  86.5       No 22.49       No     No      87
    ## 86   Female     170  52 124.0  78.0      Yes 26.03       No     No      82
    ## 87   Female     210  61 182.0 101.0       No 29.35       No     No      83
    ## 88   Female     197  36 113.0  72.5      Yes 22.73       No     No      65
    ## 89   Female     261  62 138.0  82.0       No 23.89       No     No      77
    ## 90   Female     326  61 200.0 104.0      Yes 38.46       No     No      78
    ## 91     Male     252  41 124.0  86.0      Yes 28.56       No     No      70
    ## 92     Male     286  62 117.5  80.0      Yes 31.56      Yes     No     195
    ## 93     Male     274  41 121.0  61.5       No 25.42       No     No      76
    ## 94     Male     188  53 138.0  89.0      Yes 18.23       No     No      75
    ## 95     Male     256  39 132.5  80.0      Yes 24.80       No     No      97
    ## 96   Female     244  51 102.0  71.5       No 27.38       No     No      77
    ## 97   Female     311  66 154.0  80.0       No 28.55       No     No     104
    ## 98     Male     288  52 136.0  96.0      Yes 22.70       No     No      77
    ## 99     Male     243  60 126.0  79.0      Yes 28.57       No     No      65
    ## 100  Female     193  65 123.0  76.5       No 29.33       No     No      96
    ## 101  Female     239  63 134.0  80.0      Yes 26.64      Yes     No     126
    ## 103  Female     296  56 180.0  90.0       No 23.72       No     No     120
    ## 104  Female     269  56 121.0  75.0      Yes 22.36       No     No      66
    ## 105  Female     220  47 132.5  87.0      Yes 27.98       No     No      75
    ## 106  Female     275  60 141.0  84.0       No 29.66       No     No     105
    ## 107  Female     268  45 110.0  64.0      Yes 20.68       No     No      71
    ## 108    Male     219  62 141.0  82.0       No 31.01       No     No      90
    ## 109  Female     265  48 145.0  77.0       No 24.23       No     No      64
    ## 110  Female     173  42 100.0  63.0      Yes 23.25       No     No      99
    ## 111  Female     273  63 135.0  82.0       No 26.76       No     No      56
    ## 112  Female     250  42 115.0  79.0       No 26.93       No     No      79
    ## 113    Male     290  40 138.0  90.0      Yes 27.54       No     No      73
    ## 114  Female     278  66 187.0  88.0       No 40.52       No     No      84
    ## 115    Male     197  49 123.0  69.0      Yes 29.62       No     No      60
    ## 116  Female     264  67 139.0  80.0       No 25.75       No     No      87
    ## 118    Male     233  48 138.0  88.5       No 23.62       No     No      68
    ## 119  Female     282  64 158.0 105.0       No 24.37       No     No      71
    ## 121    Male     257  50 127.0  82.0       No 32.23       No     No     117
    ## 122    Male     278  60 160.5  96.0       No 26.40       No     No      75
    ## 123  Female     185  37 100.0  68.0      Yes 18.38       No     No      72
    ## 124    Male     210  36 112.0  85.5       No 21.93       No     No      77
    ## 125    Male     175  50 157.0  88.0      Yes 25.09       No     No      85
    ## 126  Female     241  39 105.0  75.0       No 26.12       No     No      87
    ## 127  Female     237  46 112.0  70.0      Yes 20.20       No     No      62
    ## 128    Male     288  66 109.0  71.0       No 29.29       No     No      80
    ## 129  Female     200  46 135.0  97.0       No 36.81       No     No      97
    ## 130  Female     183  55 158.0  86.0       No 24.45       No     No     102
    ## 131    Male     223  41 112.5  80.0       No 29.44       No     No      58
    ## 132    Male     244  44 128.0  77.0       No 25.95       No     No      90
    ## 133    Male     235  38 118.0  77.0       No 25.87       No     No      82
    ## 134    Male     210  43 115.0  77.5      Yes 25.10       No     No      68
    ## 135  Female     213  41 112.0  73.0       No 24.81       No     No      74
    ## 136    Male     244  53 106.0  67.5      Yes 21.84       No     No      65
    ## 138    Male     303  59 134.0  89.0      Yes 24.07       No     No      75
    ## 139  Female     246  56 128.0  64.0      Yes 25.54       No     No      92
    ## 140  Female     150  50 121.0  84.0      Yes 28.69       No     No      88
    ## 141    Male     270  46 131.0  81.0      Yes 26.40       No     No      83
    ## 142  Female     266  43 131.0  81.0       No 24.03       No     No      94
    ## 143    Male     246  61 124.0  70.0       No 25.63       No     No      78
    ## 144  Female     187  41 154.0 100.0      Yes 20.50       No     No      78
    ## 145    Male     256  49 127.5  81.5      Yes 28.21       No     No      85
    ## 146  Female     286  55 138.0  82.0       No 24.27       No     No      90
    ## 147    Male     240  43 100.0  72.5      Yes 24.30       No     No      65
    ## 148    Male     154  40 117.5  72.5       No 26.82       No     No      87
    ## 149  Female     266  57 151.0  95.0       No 38.39       No     No     109
    ## 150  Female     279  56 136.0  94.0      Yes 32.99       No     No     102
    ## 151    Male     293  48 149.0 100.0       No 31.61       No     No      76
    ## 152  Female     259  59 141.0  86.0      Yes 25.97       No     No      86
    ## 153  Female     219  47 153.0  98.0       No 22.02       No     No      92
    ## 154    Male     210  64 123.0  81.0       No 26.49       No     No      75
    ## 155  Female     230  54 180.5 106.5       No 28.92       No     No      64
    ## 156  Female     320  64 130.0  77.0       No 26.24       No     No      74
    ## 157  Female     220  61 142.0  93.0       No 23.37       No     No      79
    ## 158  Female     312  61 136.5  76.0       No 31.13       No     No      85
    ## 159  Female     214  66 212.0 104.0       No 25.32       No     No      84
    ## 161    Male     209  39 134.0  82.0      Yes 28.34       No     No      75
    ## 162  Female     195  58 153.0  80.5       No 23.36       No     No      73
    ## 163  Female     265  49 150.0  77.5       No 21.83       No     No     107
    ## 164  Female     254  49 191.0 124.5       No 28.35       No    Yes      54
    ## 165  Female     165  46 108.0  81.0       No 24.19       No     No      72
    ## 166  Female     159  36 121.5  73.0      Yes 20.41       No     No      75
    ## 167    Male     244  41 139.0  86.0      Yes 30.77       No    Yes      67
    ## 168    Male     168  54 180.0  97.5      Yes 22.99       No     No     132
    ## 169  Female     174  47 118.0  86.5       No 26.15       No     No      86
    ## 170    Male     242  57 102.0  61.0      Yes 23.50       No     No      82
    ## 171    Male     301  50 117.5  80.0       No 28.04       No     No      72
    ## 172  Female     266  62 173.0  89.0       No 42.00       No     No      75
    ## 173    Male     167  39 109.0  78.0       No 21.19       No     No      71
    ## 174  Female     205  37 110.0  78.0      Yes 24.43       No     No      75
    ## 175    Male     265  41 134.0  83.0      Yes 23.35       No     No     117
    ## 176    Male     210  48 109.0  77.0      Yes 24.12       No     No      79
    ## 177  Female     197  46 144.0  78.0      Yes 22.51       No     No      60
    ## 178  Female     235  52 129.5  83.0       No 28.86       No     No      79
    ## 179    Male     220  38 122.0  80.5       No 21.66       No     No      77
    ## 180  Female     265  50 110.0  74.0      Yes 25.26       No     No      88
    ## 181  Female     200  39 111.0  64.0      Yes 19.24       No     No      60
    ## 182  Female     308  55 124.0  76.0      Yes 27.23       No     No      68
    ## 183    Male     245  62 158.0  86.5      Yes 26.51       No     No      74
    ## 184  Female     259  67 151.0 101.0       No 21.67       No     No      86
    ## 185  Female     325  60 182.0 106.0       No 27.61       No     No      77
    ## 186  Female     229  61 122.0  83.0      Yes 25.45       No     No      61
    ## 187  Female     236  39 117.5  71.0      Yes 27.27       No     No      74
    ## 188  Female     214  42 110.0  67.0      Yes 22.54       No     No      75
    ## 189  Female     300  37 112.0  60.0       No 23.67       No     No      75
    ## 190    Male     225  42 128.0  82.0      Yes 26.79       No     No      85
    ## 191    Male     215  36 118.0  76.0      Yes 18.99       No     No      97
    ## 192  Female     225  53 128.0  77.0      Yes 23.95       No     No      78
    ## 194  Female     216  52 158.0  98.0       No 24.53       No     No      86
    ## 195  Female     224  45 117.0  74.5      Yes 16.75       No     No      87
    ## 196    Male     215  56 138.0  97.0       No 30.76       No     No      69
    ## 197    Male     216  42 125.0  86.5      Yes 23.25       No     No      57
    ## 198  Female     245  55 144.5  83.5       No 28.96       No     No      65
    ## 199    Male     253  38 133.0  92.0      Yes 28.82       No     No      63
    ## 200    Male     241  56 109.0  70.0      Yes 20.12       No     No      87
    ## 201    Male     235  61 127.0  81.0       No 28.63       No     No      90
    ## 203  Female     464  42 128.0  87.0       No 22.90       No     No      72
    ## 204  Female     226  49 106.0  71.0      Yes 22.89       No     No      57
    ## 205    Male     308  48 117.0  76.0      Yes 30.85       No     No      54
    ## 206  Female     248  55 157.0  82.5      Yes 22.91       No     No      83
    ## 207  Female     215  58 170.0  86.0      Yes 29.06       No     No      98
    ## 208    Male     240  60 137.0  84.0       No 29.51       No     No      88
    ## 209  Female     171  38 111.0  68.0       No 18.76       No     No      83
    ## 210    Male     189  53 110.0  67.5      Yes 23.59       No     No      63
    ## 211    Male     240  52  94.0  66.5      Yes 22.93       No     No      88
    ## 213  Female     186  37 135.0  91.0      Yes 21.48       No     No      84
    ## 215    Male     227  42 144.0  78.0       No 23.75       No     No      97
    ## 216  Female     273  64 155.0  86.0       No 27.53       No     No      91
    ## 217  Female     212  56 117.5  72.5      Yes 27.30       No     No      75
    ## 218  Female     249  67 128.0  68.0       No 25.81       No     No      87
    ## 219    Male     185  43 125.0  65.0      Yes 20.65       No     No      76
    ## 220    Male     176  45 124.0  84.0      Yes 20.27       No     No      75
    ## 222  Female     285  66 166.0  98.0       No 26.04       No     No     132
    ## 224  Female     254  52 126.5  93.0       No 26.79       No     No      65
    ## 227    Male     237  49 114.0  85.5      Yes 28.57       No     No      92
    ## 228  Female     175  57 123.0  72.0       No 22.37       No     No      74
    ## 229  Female     196  48  96.0  70.0       No 22.72       No     No      68
    ## 230  Female     150  36 117.5  77.5      Yes 23.71       No     No      74
    ## 231  Female     180  38 124.0  66.0      Yes 29.29       No     No      68
    ## 232  Female     242  60 130.0  70.0       No 29.17       No     No      84
    ## 233  Female     173  40 106.0  82.0       No 23.05       No     No      83
    ## 234  Female     257  43 160.0  85.0       No 21.95       No     No      84
    ## 235  Female     194  62 134.0  60.0       No 30.23       No     No      85
    ## 236    Male     193  53 142.0  89.0       No 29.56       No     No      78
    ## 237    Male     243  54 135.0  92.0      Yes 31.30       No     No      65
    ## 238  Female     320  58 174.0  92.5       No 31.77       No    Yes      79
    ## 239    Male     239  56 177.5  98.0       No 29.44       No     No     105
    ## 240  Female     273  50 131.0  93.0       No 27.61       No     No      94
    ## 241    Male     220  46 136.5  99.5       No 27.78       No     No      70
    ## 242  Female     303  53 128.0  91.0       No 27.35       No     No      77
    ## 243  Female     310  59 129.0  70.0       No 23.29       No     No      70
    ## 244    Male     232  48 112.5  79.0       No 28.62       No     No     100
    ## 245  Female     164  48 159.0  90.0       No 26.73       No     No      63
    ## 246    Male     245  54 128.0  74.0      Yes 24.85       No     No      75
    ## 247  Female     135  36 108.0  74.0      Yes 22.53       No     No      75
    ## 248    Male     265  59 155.0  85.0      Yes 27.06       No     No      75
    ## 249  Female     291  55 124.0  85.0       No 26.89       No     No      80
    ## 250  Female     273  55 122.0  84.0      Yes 27.15       No     No      97
    ## 251    Male     238  40 129.0  88.0      Yes 26.32       No     No      60
    ## 252    Male     207  65 139.0  88.0       No 24.04       No     No      73
    ## 253    Male     219  65 148.0  90.0      Yes 29.35       No     No      97
    ## 254  Female     246  45 134.0  81.0       No 21.99       No     No      76
    ## 255  Female     310  62 178.0 127.0      Yes 31.77       No     No      79
    ## 256  Female     197  40 124.0  76.0      Yes 18.06       No     No      69
    ## 257    Male     260  63 130.5  82.0      Yes 20.12       No     No      72
    ## 261    Male     342  55 107.5  73.0      Yes 21.97       No     No      74
    ## 262  Female     180  60 200.0 122.5      Yes 44.27      Yes     No     150
    ## 264  Female     287  58 144.0  84.0       No 21.81       No     No      68
    ## 265  Female     182  39 109.0  70.0       No 20.59       No     No      66
    ## 266    Male     224  61 124.0  74.0       No 21.90       No     No      75
    ## 267    Male     238  51 125.0  80.0      Yes 19.36       No     No      66
    ## 268  Female     252  60 189.0 110.0       No 28.77       No    Yes      70
    ## 269    Male     212  36 168.0  98.0       No 29.77       No     No      75
    ## 270    Male     300  39 139.0  90.0      Yes 30.96       No     No     107
    ## 271  Female     165  42 108.0  72.5       No 24.85       No     No      87
    ## 272  Female     352  60 197.5 105.0       No 36.29       No     No      95
    ## 273    Male     240  52 146.0  89.0       No 24.59      Yes     No      73
    ## 274    Male     284  47 137.0  91.0      Yes 27.33       No     No      61
    ## 275  Female     219  46 150.0  81.0       No 25.43       No     No      93
    ## 277    Male     176  57 134.0  97.0       No 38.14       No     No      94
    ## 278    Male     203  45 107.5  72.5       No 22.32       No     No      83
    ## 279    Male     265  45 156.5  86.0      Yes 24.15       No     No      76
    ## 280  Female     219  55 140.0  93.0       No 27.78       No     No      70
    ## 281    Male     245  41 134.0  98.0       No 24.26       No     No      78
    ## 282  Female     312  50 125.0  85.0       No 27.07       No     No      94
    ## 283  Female     273  54 139.0  98.0       No 29.06       No     No      73
    ## 284  Female     254  57 174.0  84.5      Yes 24.22       No     No      76
    ## 285    Male     273  62 129.0  83.0       No 25.49       No     No      70
    ## 287  Female     254  46 136.0  87.0       No 31.24       No     No      80
    ## 288  Female     262  53 127.5  86.0      Yes 24.11       No     No      73
    ## 289    Male     215  52 129.0  80.0       No 29.29       No     No      87
    ## 290  Female     302  67 147.0  92.0       No 25.23       No     No      87
    ## 291  Female     185  42 127.5  91.0       No 23.78       No     No      59
    ## 293  Female     284  51 132.0  78.0       No 21.94       No     No      94
    ## 294  Female     325  60 123.0  79.0       No 25.82       No     No      85
    ## 296  Female     260  39 100.0  74.5       No 20.51       No     No      66
    ## 297    Male     164  68 142.0  85.0       No 30.28      Yes     No     120
    ## 298  Female     270  54  98.0  64.0       No 22.02       No     No      75
    ## 299  Female     186  38 105.0  57.0       No 21.10       No     No      68
    ## 300    Male     258  39 105.0  69.0       No 24.10       No     No      79
    ## 301  Female     228  36 111.0  68.0      Yes 23.86       No     No      68
    ## 302    Male     155  36 126.0  72.0      Yes 25.14       No     No      70
    ## 303  Female     270  44 131.5  76.0       No 22.19       No     No     113
    ## 304    Male     323  48 116.0  72.0       No 26.22       No     No      99
    ## 305  Female     242  54 125.0  76.0      Yes 22.16       No     No      87
    ## 308    Male     310  41 117.5  80.0      Yes 26.74       No     No      78
    ## 310  Female     206  34 101.0  63.0      Yes 21.50       No     No      66
    ## 311    Male     173  41 130.0  80.0       No 28.39       No     No      61
    ## 312  Female     283  62 130.0  80.0       No 35.20       No    Yes      83
    ## 313    Male     319  38 121.0  86.0       No 29.77       No     No      77
    ## 316  Female     282  55 158.5  81.0       No 30.24       No     No      70
    ## 317  Female     188  46  97.0  65.0      Yes 21.17       No     No      60
    ## 320    Male     221  49 101.0  61.0      Yes 23.94       No     No      79
    ## 321  Female     269  39  97.0  64.0      Yes 23.09       No     No      67
    ## 323  Female     250  47 114.0  77.0       No 24.16       No     No      93
    ## 324    Male     194  62 151.5  88.0      Yes 21.61      Yes     No     105
    ## 325    Male     258  39 110.0  65.0      Yes 19.97       No     No      85
    ## 326    Male     272  45 140.0  94.0      Yes 27.87       No     No      79
    ## 329  Female     260  63 168.0  98.0      Yes 21.05       No     No      73
    ## 330    Male     270  41 157.0 101.0       No 33.11       No     No      75
    ## 331    Male     340  44 125.0  78.0      Yes 26.41       No     No      90
    ## 333    Male     197  44 118.0  81.0      Yes 17.44       No     No      75
    ## 334    Male     240  49 148.0  85.0      Yes 29.69       No     No      68
    ## 335  Female     240  53 131.0  82.0      Yes 24.22       No     No      80
    ## 336  Female     220  51 112.0  74.0      Yes 31.23       No     No      66
    ## 337    Male     200  43 133.0  78.0       No 26.72       No     No      71
    ## 338  Female     328  60 127.0  70.0       No 22.36       No     No      63
    ## 339    Male     222  39  97.5  57.5      Yes 23.22       No     No      64
    ## 340    Male     246  56 125.0  79.0      Yes 29.64       No     No      85
    ## 341  Female     246  59 189.0 111.0       No 19.88       No     No      85
    ## 343  Female     186  43 120.0  72.0      Yes 24.33       No     No      86
    ## 344  Female     197  48 101.0  67.0      Yes 21.35       No     No     100
    ## 345  Female     311  61 110.0  78.5       No 29.04       No     No      77
    ## 346  Female     368  55 204.0  94.0       No 25.20       No    Yes      81
    ## 347    Male     220  39 145.0  99.0      Yes 26.50       No     No      90
    ## 348  Female     218  42 126.0  87.0      Yes 22.50       No     No      73
    ## 349    Male     241  51 126.5  78.5      Yes 29.91       No     No      70
    ## 350  Female     219  54 143.5  89.0       No 28.47       No     No      96
    ## 352  Female     312  64 160.0  82.0       No 27.59       No     No      94
    ## 353  Female     200  58 158.0 101.0       No 23.06       No     No      77
    ## 354  Female     237  42 105.0  75.0      Yes 23.85       No     No      87
    ## 355  Female     207  41 111.0  60.0      Yes 18.48       No     No      76
    ## 358  Female     238  63 136.0  66.0      Yes 20.20       No     No      92
    ## 360  Female     250  53 144.0  98.0       No 28.78       No    Yes      78
    ## 361    Male     238  41 128.0  86.0      Yes 26.09       No     No      72
    ## 362  Female     276  53 130.0  86.0      Yes 27.09       No     No      56
    ## 363  Female     339  59 110.0  69.0       No 26.89       No     No      73
    ## 364  Female     231  45 157.5 104.5       No 22.86       No     No      92
    ## 365    Male     187  49 110.0  67.0      Yes 19.26       No     No      85
    ## 366    Male     198  41 106.0  71.0      Yes 21.51       No     No      84
    ## 367    Male     196  39 121.0  83.0      Yes 25.38       No     No      68
    ## 368    Male     225  50 137.0  89.5      Yes 25.77       No     No      91
    ## 370    Male     202  53 140.5  90.0       No 25.82       No     No      74
    ## 371  Female     231  54 127.5  83.0      Yes 21.31      Yes     No     115
    ## 372  Female     210  40 118.0  79.0      Yes 21.21       No     No      84
    ## 373    Male     230  51 171.0 112.0      Yes 25.08       No     No      88
    ## 374    Male     283  40 146.0  95.0      Yes 23.27       No     No      64
    ## 375  Female     203  38 100.0  70.0      Yes 22.73       No     No      80
    ## 376  Female     200  48 118.0  78.0       No 24.53       No     No      98
    ## 377  Female     244  46  98.0  57.0      Yes 24.01       No     No      95
    ## 378    Male     212  59 106.0  66.0      Yes 26.46       No     No     117
    ## 379    Male     218  37 130.0  89.0       No 22.70       No     No      88
    ## 380    Male     185  40 133.0  80.0      Yes 21.86       No     No      58
    ## 381    Male     282  59 114.0  67.0       No 28.04       No     No      79
    ## 382    Male     198  46 128.0  85.0      Yes 29.00       No     No      77
    ## 383  Female     248  67 215.0 105.0       No 22.91       No     No      97
    ## 384  Female     290  40 122.0  85.0       No 24.06       No     No      70
    ## 385    Male     239  55 159.0 102.0       No 32.35       No     No      71
    ## 386  Female     245  52 145.0  89.0       No 25.02       No     No      80
    ## 387  Female     200  42  95.0  55.0      Yes 23.68       No     No      83
    ## 388  Female     201  54 156.0  93.0       No 27.91       No    Yes      70
    ## 389  Female     290  58 155.0  82.5      Yes 29.50       No     No     113
    ## 390  Female     285  45 116.0  87.0      Yes 23.85       No     No      55
    ## 391    Male     203  67 122.0  74.0       No 15.54       No     No      79
    ## 392    Male     250  68 109.0  73.0      Yes 24.68       No     No      66
    ## 393  Female     290  61 170.0  98.0       No 26.98       No     No      84
    ## 394  Female     258  55 123.0  89.0       No 31.35       No     No      84
    ## 395  Female     246  42 125.0  80.0       No 29.02       No     No      98
    ## 396  Female     296  59 154.0 100.0       No 24.04       No    Yes      73
    ## 398    Male     220  39 135.0  85.0      Yes 27.17       No     No      81
    ## 399    Male     195  40 122.5  66.5      Yes 19.98       No     No      72
    ## 400    Male     235  39 120.0  80.0       No 27.23       No     No      87
    ## 401  Female     223  63 178.0  88.5       No 27.18       No     No      91
    ## 402  Female     308  55 124.0  87.0       No 31.82       No     No      84
    ## 403  Female     167  39 127.5  86.5       No 28.28       No     No      82
    ## 404  Female     185  39 111.0  67.0      Yes 23.87       No     No      87
    ## 405    Male     240  51 146.5  83.0      Yes 25.62       No     No     140
    ## 406  Female     245  54 117.0  76.0       No 26.64       No     No      76
    ## 407  Female     310  48 124.0  74.0       No 25.94       No     No      73
    ## 408    Male     300  51 134.0  86.5      Yes 24.76       No     No      81
    ## 409  Female     256  65 149.0  96.0       No 29.75       No     No      80
    ## 411    Male     202  39 136.5  79.0       No 24.35       No     No      60
    ## 412  Female     155  47 122.5  77.5      Yes 21.34       No     No      78
    ## 413    Male     207  59 132.5  66.0      Yes 26.84       No     No      76
    ## 414    Male     276  66 159.0  82.0       No 31.42       No     No      65
    ## 416    Male     215  41 113.0  68.0       No 25.13       No     No      87
    ## 417  Female     277  41 113.5  75.0       No 29.73       No     No      90
    ## 418  Female     179  41 116.0  67.0       No 18.58       No     No      68
    ## 419    Male     226  45 119.0  75.0      Yes 25.34       No     No      66
    ## 420    Male     245  37 138.0  84.0      Yes 27.45       No     No      76
    ## 421  Female     240  43 141.0  91.0       No 29.25       No     No      65
    ## 423  Female     288  51 176.0 119.0       No 28.70      Yes     No      82
    ## 424    Male     218  41 137.0  86.0       No 25.52       No     No      88
    ## 425    Male     229  44 177.5 120.0      Yes 39.88       No     No      78
    ## 426    Male     174  36 131.0  86.0       No 26.23       No     No      64
    ## 427  Female     285  60 156.0 100.0       No 23.02       No     No      85
    ## 430  Female     323  55 197.0 118.0      Yes 27.51       No     No     112
    ## 431  Female     185  63 221.0 118.0       No 24.96       No    Yes      76
    ## 432  Female     206  61 130.0  80.0      Yes 21.93       No     No      82
    ## 433    Male     246  45 111.0  72.0      Yes 21.79       No     No     118
    ## 434  Female     304  45 148.0 106.0       No 22.98       No     No      72
    ## 436  Female     235  52 119.0  82.0      Yes 24.25       No     No      79
    ## 437  Female     232  46  90.0  60.0       No 19.20       No     No      82
    ## 439  Female     177  42 112.5  70.0      Yes 20.62       No     No      83
    ## 440  Female     188  58 160.0 120.0       No 35.58       No     No      85
    ## 441  Female     176  38 110.0  80.0      Yes 24.03       No     No     113
    ## 442    Male     276  63 147.0  90.0      Yes 31.80       No     No      75
    ## 444  Female     268  64 135.0  74.0       No 30.18       No     No      83
    ## 448  Female     195  40 109.5  72.0       No 22.36       No     No      85
    ## 449  Female     197  57  96.0  64.0       No 18.59       No     No      77
    ## 450  Female     268  48 117.5  80.0      Yes 36.11       No     No      67
    ## 452    Male     212  40 120.0  80.0      Yes 24.01       No     No      57
    ## 453    Male     220  43 165.0 105.0      Yes 31.07       No     No     114
    ## 455    Male     240  45 141.0  89.0      Yes 25.01       No     No      76
    ## 456    Male     226  36 124.0  76.0      Yes 25.84       No     No      70
    ## 457    Male     237  63 113.0  80.0       No 27.61       No     No      71
    ## 458  Female     199  33 116.0  81.0      Yes 21.61       No     No      93
    ## 459  Female     206  56 140.0  87.0       No 27.72       No     No      85
    ## 460  Female     265  43 107.0  68.0      Yes 21.08       No     No      95
    ## 461  Female     260  60  95.5  59.0       No 25.94      Yes     No     160
    ## 462    Male     243  44 146.0  91.0       No 26.72       No     No     104
    ## 463    Male     258  61 130.0  70.0      Yes 24.35       No     No      78
    ## 464    Male     182  61 209.0 133.0       No 30.77       No     No      75
    ## 465    Male     209  46 112.5  65.0      Yes 27.48       No     No      78
    ## 467  Female     266  48 115.0  75.0      Yes 31.16       No     No      90
    ## 468  Female     210  55 112.0  76.0       No 20.53       No     No      65
    ## 469  Female     292  50 132.0  82.0      Yes 22.54      Yes     No     110
    ## 470  Female     190  46 107.5  80.0       No 25.13       No     No      75
    ## 471  Female     234  47 127.5  83.5       No 27.65       No     No      89
    ## 472  Female     173  44 136.5  77.5       No 26.62       No     No      72
    ## 475    Male     227  55 129.0  81.0       No 26.57       No     No      74
    ## 476  Female     186  47 150.0  85.0      Yes 22.53       No    Yes      93
    ## 477  Female     195  48 109.0  71.0      Yes 21.10       No     No      65
    ## 478  Female     284  40 124.0  83.0      Yes 27.90       No     No      71
    ## 479  Female     177  35 110.0  70.0      Yes 25.71       No     No      84
    ## 480  Female     305  49 135.0  89.0       No 25.04       No     No     117
    ## 481    Male     195  38 146.0 104.0      Yes 29.60       No     No      80
    ## 482    Male     227  50 152.0  99.0       No 29.46       No     No     115
    ## 483    Male     245  38 154.0  95.5      Yes 30.02       No     No      87
    ## 484  Female     250  57 125.0  74.0      Yes 21.08       No     No      72
    ## 486  Female     152  44 110.0  64.0      Yes 25.71       No     No      83
    ## 487  Female     222  47 162.5  92.5       No 23.45       No     No      80
    ## 488    Male     243  56 134.0  86.0       No 25.94       No     No      66
    ## 489    Male     270  61 177.5  95.0       No 28.15      Yes    Yes     123
    ## 491    Male     252  56 107.0  66.0       No 25.11       No     No      72
    ## 492    Male     235  46 150.0  96.5      Yes 27.22       No     No      78
    ## 493    Male     161  53 116.0  74.5      Yes 19.82       No     No      90
    ## 494    Male     168  49 127.0  74.0       No 27.38       No     No      85
    ## 495  Female     247  46 125.0  80.0      Yes 21.51       No     No      80
    ## 496  Female     250  37 138.0  90.0      Yes 19.56       No     No      74
    ## 497    Male     219  50 145.0 100.0       No 26.26       No     No     108
    ## 498    Male     266  61 103.0  70.0      Yes 23.88       No     No      93
    ## 499  Female     320  51 145.0  85.0      Yes 24.03       No     No      98
    ## 500    Male     222  37 110.0  75.0      Yes 26.55       No     No      78
    ## 501  Female     253  64 295.0 135.0       No 38.82       No     No      70
    ## 502    Male     319  44 158.0  90.0      Yes 29.15       No     No     100
    ## 503    Male     176  37 125.0  80.0       No 26.75       No     No      65
    ## 504  Female     239  62 126.0  81.0       No 29.35       No     No      70
    ## 505    Male     230  54 145.0  90.0      Yes 25.72       No     No      85
    ## 506  Female     190  39 120.0  80.0       No 27.16       No     No      85
    ## 508  Female     181  44 150.0 101.0      Yes 23.74       No     No      86
    ## 509  Female     220  51  98.0  64.5       No 21.14       No     No      79
    ## 511    Male     207  63 130.0  72.0       No 26.76       No     No      81
    ## 512  Female     269  46 134.0  78.0      Yes 26.80       No     No     104
    ## 513    Male     258  40 135.0  93.0      Yes 31.69       No     No      57
    ## 514    Male     340  56 134.0  89.5      Yes 21.91       No     No      72
    ## 515    Male     277  44 114.0  81.0      Yes 27.51       No     No      76
    ## 516    Male     245  51 137.0  76.0      Yes 22.26       No     No      73
    ## 517    Male     233  45 151.0  95.0      Yes 29.17       No     No      95
    ## 518    Male     315  40 121.0  76.0      Yes 26.32       No     No      85
    ## 519    Male     258  56 111.0  68.5       No 24.04       No     No      70
    ## 521  Female     210  40 103.0  71.0      Yes 24.40       No     No      68
    ## 522  Female     248  52 155.0  93.0       No 23.09       No     No      70
    ## 523  Female     251  64 132.0  82.0       No 28.87       No     No      82
    ## 525  Female     258  50 144.0  88.0       No 24.19       No     No      83
    ## 526  Female     201  41 108.0  71.0      Yes 20.47       No     No      75
    ## 527    Male     303  57 160.5  98.5       No 25.84       No     No     100
    ## 528  Female     248  51 139.0  81.0      Yes 31.16       No     No      95
    ## 529  Female     296  50 134.5  89.0       No 24.91       No     No      76
    ## 530  Female     218  46 115.5  62.0      Yes 23.48       No     No      77
    ## 531  Female     235  53 132.0  84.0       No 30.10       No     No      63
    ## 533  Female     240  65 155.0  84.0       No 29.93       No     No      91
    ## 534  Female     234  58 137.5  89.5       No 24.25       No     No      72
    ## 535    Male     278  56 133.0  84.0      Yes 22.67       No     No      96
    ## 536    Male     269  46 135.0  95.0      Yes 26.55       No     No      92
    ## 537    Male     203  64 130.0  78.0       No 28.66       No     No      70
    ## 538    Male     243  43 155.0 100.0      Yes 24.89       No     No      63
    ## 539    Male     231  35 122.0  72.5      Yes 22.78       No     No      93
    ## 541  Female     271  59 117.5  65.0       No 19.77       No     No      89
    ## 542  Female     261  50 129.0  80.0       No 23.06       No     No      90
    ## 543    Male     203  62 148.0  96.0      Yes 30.84       No     No      95
    ## 544  Female     210  53 132.0  84.5      Yes 27.08       No     No      84
    ## 545    Male     206  52 173.0 117.0      Yes 29.63       No     No      77
    ## 546  Female     301  57 119.0  80.0       No 24.79       No     No      73
    ## 547  Female     241  46 130.0  82.0      Yes 34.84       No     No      93
    ## 548    Male     215  49 122.5  76.0      Yes 27.17       No     No      61
    ## 549    Male     261  47 105.0  74.0      Yes 20.57       No     No      75
    ## 550  Female     266  65 120.0  73.0       No 24.33       No     No      69
    ## 551  Female     237  41 125.0  83.0       No 24.56       No     No      61
    ## 552  Female     217  61 189.0 121.0       No 37.41       No     No     100
    ## 553  Female     176  41 113.0  75.0      Yes 22.29       No     No      55
    ## 554  Female     370  53 123.0  83.0       No 24.64       No     No      74
    ## 555    Male     182  36 145.0 102.0       No 30.92       No     No      70
    ## 556    Male     245  43 105.0  59.5       No 30.55       No     No      77
    ## 557    Male     292  37 130.0  85.0      Yes 28.07       No     No      81
    ## 558  Female     246  62 171.0 101.0       No 23.88       No     No      78
    ## 559    Male     193  64 114.0  79.0       No 16.59       No     No      64
    ## 560  Female     247  40 125.0  83.0      Yes 22.55       No     No      80
    ## 562  Female     300  49 131.0  88.0       No 26.85       No    Yes      73
    ## 563  Female     213  46 115.0  72.5       No 19.98       No     No     107
    ## 564    Male     439  47 145.0  74.0      Yes 22.42       No     No      90
    ## 565  Female     171  38 110.0  71.0       No 21.80       No     No      78
    ## 566    Male     227  39 132.0  92.0      Yes 26.47       No     No      74
    ## 567    Male     207  43 117.0  65.0      Yes 24.42       No     No     100
    ## 568  Female     145  43 112.0  83.0       No 22.36       No     No      70
    ## 569    Male     210  50 121.5  78.0      Yes 26.29       No     No      77
    ## 570    Male     245  48 144.0  96.5       No 32.58       No     No      77
    ## 571    Male     258  45 114.0  80.0       No 26.60       No     No      68
    ## 572  Female     233  55 128.0  94.0      Yes 36.62       No     No      95
    ## 574    Male     310  63 142.5  90.0      Yes 22.63       No     No      70
    ## 575    Male     240  53 174.5 103.0       No 29.82       No     No      74
    ## 576  Female     213  50 140.0  82.0      Yes 22.18       No     No      72
    ## 577  Female     242  40 110.0  75.0       No 16.69       No     No      68
    ## 578  Female     214  50 129.0  76.0       No 26.39       No     No      75
    ## 579    Male     241  56 130.0  85.0       No 25.79       No     No      78
    ## 580  Female     250  44 136.5  83.5      Yes 21.33       No     No      95
    ## 581  Female     240  37 120.0  79.0      Yes 23.09       No     No      80
    ## 582  Female     188  42 122.5  75.0       No 24.56       No     No      68
    ## 583  Female     275  54 113.5  75.5       No 19.63       No     No      78
    ## 584  Female     252  39 128.0  93.0       No 30.36       No     No      63
    ## 586    Male     178  62 164.0 102.5       No 28.35       No     No      83
    ## 587    Male     253  42 116.0  72.0       No 21.96       No     No      88
    ## 588  Female     278  55 125.0  80.0       No 26.27       No     No      78
    ## 589    Male     232  60 165.0  77.0       No 29.23       No     No      57
    ## 591  Female     189  34 108.0  66.0       No 20.81       No     No      88
    ## 592  Female     242  44 135.0  89.0      Yes 23.29       No     No      77
    ## 593  Female     247  49 121.0  82.0      Yes 29.07       No     No      69
    ## 594    Male     230  50 133.0  91.0       No 25.74       No     No      70
    ## 595  Female     223  47 120.0  74.0       No 20.27       No     No      76
    ## 596  Female     263  63 150.0  96.5       No 24.85       No     No      75
    ## 597  Female     277  39 148.0 100.0      Yes 24.12       No     No      72
    ## 598    Male     215  43 122.0  76.0      Yes 26.84       No     No      74
    ## 601    Male     210  41 121.5  82.5      Yes 23.27       No     No      89
    ## 602    Male     300  59 163.0  78.0       No 28.83       No     No      95
    ## 603  Female     287  51 139.0  92.0       No 37.48       No     No      74
    ## 604    Male     234  62 117.0  80.0      Yes 26.97       No     No      67
    ## 605  Female     176  51 146.0  94.0      Yes 27.42       No     No      85
    ## 606  Female     240  46 125.0  74.0      Yes 22.89       No     No      76
    ## 607    Male     300  37 118.5  85.5       No 25.83       No     No      82
    ## 608  Female     258  60 142.0  87.0       No 32.53      Yes     No     145
    ## 609    Male     215  52  98.0  65.0      Yes 20.87       No     No      76
    ## 610    Male     197  40 116.0  73.0      Yes 24.01       No     No      83
    ## 611  Female     210  52 185.0 114.0      Yes 27.01       No     No      83
    ## 612  Female     220  45 126.0  82.0       No 23.87       No     No      90
    ## 613  Female     294  60 220.0 118.0      Yes 24.22       No    Yes      59
    ## 615    Male     232  39 115.0  72.5      Yes 30.22       No     No     105
    ## 616  Female     220  46 123.0  88.0      Yes 32.49       No     No      79
    ## 617  Female     205  40 125.0  73.5      Yes 20.68       No     No      99
    ## 618    Male     198  44 119.0  82.0       No 28.04       No     No      66
    ## 619    Male     208  38 121.5  74.0       No 27.05       No     No      90
    ## 620    Male     276  60 144.0  78.0       No 26.98       No     No      88
    ## 621  Female     216  42 120.0  70.0      Yes 21.93       No     No      88
    ## 622    Male     254  44 123.0  90.0       No 24.44       No     No      72
    ## 623  Female     246  47 113.0  75.0      Yes 21.66       No     No      68
    ## 624  Female     330  55 103.0  73.0       No 24.50       No     No      67
    ## 625    Male     170  48 132.0  91.0      Yes 27.61       No     No      57
    ## 626  Female     245  61 140.0  73.0      Yes 30.74       No     No      91
    ## 627    Male     286  62 164.0  88.0      Yes 19.53       No     No     126
    ## 628  Female     215  53 153.0 108.0      Yes 27.31       No     No      70
    ## 629  Female     250  64 145.0  79.0      Yes 25.16       No     No      86
    ## 630    Male     270  46 120.5  78.0      Yes 22.54       No     No      83
    ## 631  Female     287  54 145.0  91.0       No 23.81       No     No      83
    ## 632    Male     218  42 116.0  86.0      Yes 17.81       No     No      69
    ## 634    Male     222  55 155.0  92.5      Yes 28.35       No     No      68
    ## 635  Female     170  35  98.5  69.5       No 19.64       No     No      77
    ## 636  Female     197  40 107.0  61.0       No 23.65       No     No      80
    ## 637    Male     157  44 108.0  70.0      Yes 20.56       No     No      82
    ## 638  Female     244  52 127.5  72.5      Yes 24.29       No     No     118
    ## 639  Female     252  63 120.5  75.0       No 24.30       No     No      77
    ## 640    Male     203  39 117.5  77.5       No 27.29       No     No      60
    ## 641  Female     398  51 161.0  96.0      Yes 23.63       No     No      83
    ## 642    Male     207  59 148.0  98.0       No 25.63       No     No      93
    ## 643  Female     264  51 139.5  89.0      Yes 29.38       No     No      76
    ## 644    Male     214  67 127.5  80.0       No 22.11       No     No      84
    ## 645  Female     272  47 127.5  87.5      Yes 22.35       No     No      72
    ## 647    Male     200  50 105.0  68.0      Yes 23.30       No     No      68
    ## 648    Male     269  63 180.0 101.0      Yes 24.42       No     No      84
    ## 649  Female     250  39 123.0  92.5       No 29.23       No     No      71
    ## 650    Male     219  38 115.0  71.0      Yes 24.58       No     No      93
    ## 651  Female     187  42  96.0  67.0      Yes 24.23       No     No      84
    ## 652    Male     243  38 130.5  85.5      Yes 19.53       No     No      79
    ## 653  Female     190  52 117.0  75.0      Yes 21.48       No     No      67
    ## 654    Male     267  60 168.5 108.0      Yes 26.67       No     No      83
    ## 655    Male     352  61 103.5  70.0       No 27.86       No     No     105
    ## 656    Male     218  41 129.5  93.0      Yes 27.80       No     No      83
    ## 657  Female     195  41 148.0 108.0      Yes 18.21       No     No      69
    ## 658  Female     266  56 134.5  78.5      Yes 30.78       No     No      84
    ## 659  Female     229  38 117.5  67.5      Yes 23.47       No     No      80
    ## 660  Female     264  63 176.0  99.0       No 23.62       No     No      65
    ## 661  Female     273  42 111.0  73.0      Yes 19.27       No     No      89
    ## 662  Female     162  40 109.0  71.0       No 20.99       No     No      85
    ## 663    Male     250  51 125.0  80.0       No 26.98       No     No     108
    ## 664  Female     314  60 141.0  93.0      Yes 25.23       No     No      94
    ## 665    Male     219  40 131.0  85.5      Yes 31.96       No     No      74
    ## 666    Male     250  34 119.0  77.0      Yes 29.04       No     No      80
    ## 667  Female     230  40 123.5  81.0       No 27.91       No     No      65
    ## 668  Female     195  54 163.5 103.0       No 26.89       No     No      80
    ## 669    Male     194  36 117.0  90.0       No 27.08       No     No      87
    ## 671    Male     224  38 127.5  81.0       No 20.39       No     No      90
    ## 672    Male     166  37 128.5  83.0       No 26.81       No     No     108
    ## 673  Female     160  39 128.5  74.5      Yes 20.56       No     No      83
    ## 674  Female     238  60 134.0  84.0       No 27.49       No     No      66
    ## 675  Female     281  58 134.0  81.0      Yes 22.54       No     No      74
    ## 676  Female     292  56 111.0  70.0      Yes 23.17       No     No      74
    ## 677  Female     203  50 128.5  82.0      Yes 18.99       No     No      84
    ## 678    Male     250  57 148.0  91.0      Yes 27.60       No     No     103
    ## 679  Female     240  41 107.0  68.5      Yes 23.47       No     No      83
    ## 680    Male     193  40 129.0  73.0      Yes 19.11       No     No      80
    ## 681    Male     243  51 130.0  86.5      Yes 29.86       No     No      74
    ## 682    Male     292  43 121.0  75.0      Yes 21.73       No     No      82
    ## 683  Female     200  36 112.0  76.0       No 26.03       No     No      62
    ## 684  Female     295  64 127.0  78.0       No 22.89       No     No      73
    ## 685    Male     255  50 126.5  92.5      Yes 25.35       No    Yes      74
    ## 687  Female     285  49 138.0  91.0       No 25.94       No     No      84
    ## 688  Female     206  59 167.0  89.5      Yes 25.83       No     No      75
    ## 689  Female     300  62 205.5 104.5       No 32.19       No     No     117
    ## 690  Female     289  56 150.0  92.0      Yes 25.68       No    Yes      84
    ## 691  Female     287  57 136.0  92.0       No 26.24       No     No      71
    ## 692    Male     253  48 118.5  73.0      Yes 26.25       No     No      70
    ## 693  Female     199  38 112.0  68.5      Yes 23.88       No     No      67
    ## 694  Female     262  49 145.0  81.0       No 23.10       No     No      83
    ## 695  Female     355  65 138.0  79.0       No 28.38       No     No     108
    ## 696    Male     223  42 115.0  70.0      Yes 24.58       No     No      81
    ## 697    Male     307  45 110.0  78.0      Yes 28.57       No     No      69
    ## 699    Male     220  38 119.5  85.5      Yes 31.39       No     No      85
    ## 701  Female     240  48 119.0  80.0       No 31.67       No     No      79
    ## 702    Male     234  58 167.5  95.0      Yes 27.73       No     No      66
    ## 703    Male     239  60 130.0  78.0      Yes 28.36       No     No      99
    ## 705  Female     225  55 136.5  85.5       No 20.56       No     No      90
    ## 706    Male     267  54 148.0  92.5       No 26.58       No     No      98
    ## 707    Male     227  45 166.0 102.0      Yes 29.32       No     No      85
    ## 708  Female     300  64 144.0  80.0       No 25.81       No     No     102
    ## 709  Female     248  49 137.0  79.0      Yes 21.60       No     No      74
    ## 710  Female     215  45 152.5  82.0       No 25.92       No     No      75
    ## 711    Male     215  48 127.0  81.0       No 24.87       No     No      76
    ## 712  Female     156  45 119.0  83.0      Yes 22.02       No     No      78
    ## 713  Female     289  57 125.0  74.0       No 18.64       No     No      69
    ## 714    Male     219  64 172.5  75.0      Yes 29.29       No     No      91
    ## 716    Male     329  50 126.0  84.0      Yes 21.45       No     No      65
    ## 717  Female     317  59 176.0  98.0       No 28.70       No     No      63
    ## 718  Female     270  58 145.0  95.0       No 23.17       No     No     120
    ## 719    Male     267  40 146.0  93.5      Yes 27.47       No     No      89
    ## 721  Female     272  54 132.5  91.0      Yes 23.09       No     No      78
    ## 722  Female     195  39 105.0  70.0       No 26.97       No     No      64
    ## 723  Female     203  46 117.0  68.0      Yes 21.50       No     No      85
    ## 724  Female     212  60 186.0 102.0      Yes 23.06       No     No      60
    ## 725    Male     190  62 183.0  93.0       No 28.96       No    Yes      77
    ## 726  Female     222  60 118.0  73.0       No 24.48       No     No      90
    ## 727    Male     206  40 110.0  67.5       No 25.88       No     No      80
    ## 728  Female     262  54 136.0  86.0      Yes 23.28       No     No      69
    ## 729  Female     294  47 109.0  72.5       No 28.59       No     No      77
    ## 730    Male     251  58 135.0  85.5      Yes 21.24       No     No     103
    ## 732  Female     160  51 140.5  84.0       No 26.28       No     No     122
    ## 733    Male     143  47 114.0  79.0      Yes 26.59       No     No      72
    ## 734  Female     245  35 148.0  84.0      Yes 23.74       No     No      73
    ## 735    Male     200  44 120.0  80.0      Yes 31.44       No     No      74
    ## 736    Male     215  50 134.0  80.0      Yes 29.98       No     No      75
    ## 737    Male     185  57 134.0  90.5      Yes 27.77       No     No     103
    ## 738    Male     212  39 107.5  67.5      Yes 27.43       No     No      77
    ## 739    Male     255  41 126.0  78.0      Yes 25.48       No     No      76
    ## 741  Female     211  45 127.5  80.0       No 27.05       No     No      68
    ## 742    Male     212  49 141.0  99.0      Yes 25.94       No     No      70
    ## 743  Female     200  38 124.0  84.0      Yes 20.67       No     No      75
    ## 744  Female     278  59 141.0  79.0       No 26.45       No     No      94
    ## 746  Female     240  39 120.0  80.0      Yes 24.79       No     No      75
    ## 747    Male     298  59 153.5 105.0      Yes 25.05       No     No      84
    ## 748    Male     271  64 134.0  79.0       No 24.95       No     No      90
    ## 749  Female     283  51 152.0  99.0       No 31.63       No    Yes      73
    ## 750    Male     241  57 163.5  93.0       No 31.68       No     No      67
    ## 751  Female     251  40 135.0  87.0       No 31.60       No     No      80
    ## 752    Male     292  62 130.0  77.5      Yes 24.75       No     No      72
    ## 753  Female     245  53 103.0  68.0       No 21.80       No     No      63
    ## 754  Female     168  37 117.0  74.0      Yes 21.51       No     No      77
    ## 755  Female     197  48 107.0  73.0       No 19.78       No     No      76
    ## 756  Female     152  40 122.0  78.0       No 18.92       No     No      80
    ## 757  Female     225  51 155.0  92.5       No 23.84       No     No      63
    ## 758  Female     296  48 117.0  73.0      Yes 24.59       No     No      78
    ## 759    Male     211  64 120.0  75.0      Yes 18.70       No     No      61
    ## 760    Male     334  44 131.0  74.0      Yes 28.82       No     No      77
    ## 761  Female     241  64 144.0  84.0       No 30.93       No     No      66
    ## 762    Male     196  41 147.5  92.5      Yes 22.72       No     No      71
    ## 763  Female     202  62 212.0 103.0       No 24.94      Yes     No     105
    ## 765    Male     215  49 183.0 114.5       No 29.25       No    Yes      90
    ## 766    Male     188  39 120.0  74.0      Yes 26.48       No     No      80
    ## 767  Female     290  61 178.0  91.0      Yes 28.87       No     No      80
    ## 768  Female     229  58 128.0  76.0       No 32.49       No     No      75
    ## 769    Male     196  52 126.0  80.0       No 22.32       No     No      74
    ## 770  Female     249  48 132.0  78.0      Yes 23.10       No     No     137
    ## 771  Female     210  45 127.0  76.5      Yes 21.67       No     No      72
    ## 772    Male     235  46 143.5  92.0      Yes 28.87       No     No     115
    ## 773  Female     232  66 175.0  94.0       No 29.84       No     No      67
    ## 774    Male     205  44 160.0 108.0      Yes 22.92       No     No      76
    ## 775    Male     245  60 119.0  82.0       No 28.56       No     No      89
    ## 776    Male     249  42 142.5  90.0       No 26.14       No     No      82
    ## 777    Male     160  42 100.0  73.0       No 22.56       No     No      93
    ## 778  Female     251  67 192.0 102.0       No 44.09       No    Yes      62
    ## 779  Female     207  63 127.0  75.0       No 22.26       No     No      77
    ## 780    Male     214  62 110.0  62.5       No 23.80       No     No      95
    ## 781    Male     296  50 129.0  85.0       No 25.86       No     No     106
    ## 782  Female     250  63 117.5  75.0       No 25.88       No     No      91
    ## 784  Female     193  41 134.5  83.0      Yes 22.28       No     No     127
    ## 786  Female     278  53 131.0  87.0      Yes 33.38       No     No      74
    ## 787    Male     254  53 160.0  92.0       No 25.61       No     No      80
    ## 789  Female     192  42  96.5  71.0      Yes 26.03       No     No      68
    ## 790    Male     184  47 107.0  66.5      Yes 16.87       No     No      70
    ## 791  Female     267  58 157.0  94.0       No 33.32      Yes     No     205
    ## 792    Male     162  45 125.0  89.0       No 27.98       No     No      83
    ## 794    Male     220  55 180.0 108.0      Yes 23.59       No     No      90
    ## 795    Male     204  67 185.0 100.0      Yes 22.69       No     No     114
    ## 796    Male     280  56 147.0  94.0      Yes 28.30       No     No      85
    ## 798    Male     180  38 123.0  69.0      Yes 22.70       No     No      70
    ## 799  Female     210  63 148.0  85.5       No 24.01       No     No      88
    ## 800    Male     246  40 141.0 102.0      Yes 26.17       No     No      67
    ## 801  Female     262  45 133.0  83.0       No 22.19       No     No      92
    ## 802  Female     193  54 118.0  84.0       No 24.90       No     No      82
    ## 803  Female     218  67 160.0 100.0       No 20.50       No     No      71
    ## 804  Female     152  39 114.0  66.5       No 24.56       No     No      71
    ## 808  Female     253  46 118.0  82.0      Yes 19.70       No     No      70
    ## 809    Male     175  35 112.0  62.5      Yes 21.03       No     No      69
    ## 811    Male     177  47 150.0 101.0      Yes 28.96       No     No      60
    ## 812  Female     199  43 137.0  81.0      Yes 21.85       No     No      72
    ## 813    Male     211  52 118.5  82.5      Yes 24.83       No     No      80
    ## 814    Male     258  45 126.0  82.0       No 27.18       No    Yes      70
    ## 815    Male     256  42 121.5  74.0      Yes 23.59       No     No     115
    ## 816  Female     207  56 133.0  76.0       No 23.42       No     No      85
    ## 817    Male     225  42 140.5  94.0      Yes 31.14       No     No      73
    ## 818    Male     260  63 159.5  91.0       No 27.01       No     No      66
    ## 819  Female     272  52 118.5  69.0       No 18.98       No     No      75
    ## 820    Male     208  46 164.0 107.0      Yes 20.63       No     No      70
    ## 821    Male     222  49 118.5  82.0      Yes 28.47       No     No      87
    ## 822  Female     219  44 129.0  90.0       No 33.47       No     No      73
    ## 823    Male     228  42 130.0  92.0       No 24.86       No     No      76
    ## 824  Female     230  64 177.0 110.0      Yes 28.91       No    Yes     113
    ## 825  Female     167  59 156.0 104.0      Yes 15.96       No     No      45
    ## 827    Male     175  54 116.0  78.5       No 25.82       No     No      67
    ## 828  Female     237  67 152.0  82.0       No 29.40       No    Yes      91
    ## 829    Male     280  57 117.0  80.0      Yes 26.56       No     No      60
    ## 830  Female     242  53 127.0  79.0      Yes 19.64       No     No      74
    ## 831  Female     231  35 150.0  90.0      Yes 23.09       No     No      72
    ## 832  Female     238  46 162.0  87.5       No 26.95       No     No      93
    ## 833  Female     314  40 124.0  86.5       No 26.79       No     No      66
    ## 834    Male     219  47 132.0  91.0       No 27.93       No     No      80
    ## 835    Male     219  49 130.0  82.0      Yes 24.35       No     No      65
    ## 836  Female     326  58 120.0  70.0       No 24.69       No     No      68
    ## 837    Male     255  54 142.0  93.0      Yes 22.17       No     No     118
    ## 838  Female     191  47 125.0  72.5      Yes 23.81       No     No      85
    ## 839  Female     239  60 164.0  94.5      Yes 25.01       No     No      89
    ## 840    Male     292  53 112.5  82.5       No 25.04       No     No      82
    ## 841    Male     233  66 108.0  66.0       No 25.16       No     No      72
    ## 844    Male     237  36 142.0  82.0       No 27.50       No     No      87
    ## 845    Male     231  37 135.0  79.0      Yes 28.46       No     No      67
    ## 846    Male     260  49 111.0  70.0      Yes 24.24       No     No      87
    ## 847    Male     225  44 128.0  82.0      Yes 26.79       No     No      85
    ## 848    Male     210  41 124.0  79.0       No 25.26       No     No      91
    ## 850    Male     205  37 129.0  86.0       No 27.27       No     No      76
    ## 851  Female     305  55 168.0  82.0       No 26.45       No     No      78
    ## 853  Female     178  36 102.5  65.0      Yes 20.87       No     No      94
    ## 854  Female     217  37 110.0  78.5      Yes 32.26       No     No      84
    ## 855  Female     328  47 134.0  87.0      Yes 22.34       No     No      99
    ## 856    Male     232  50 127.5  85.0       No 25.09       No     No      79
    ## 857    Male     261  68 122.0  70.0      Yes 26.37       No     No      60
    ## 858    Male     217  44 126.0  85.0      Yes 28.49       No     No      68
    ## 860  Female     197  42 105.0  84.0       No 25.75       No     No      67
    ## 861  Female     173  50 147.0  90.0       No 24.06       No     No      85
    ## 863  Female     248  53 200.0 140.0       No 43.30      Yes     No     130
    ## 865    Male     205  57 127.0  75.0      Yes 20.55       No     No      65
    ## 866  Female     234  47 128.0  91.0      Yes 25.59       No     No      93
    ## 867    Male     240  49 137.0  96.5      Yes 23.38       No     No     118
    ## 868    Male     266  58 138.0  83.0       No 30.64       No     No     100
    ## 869  Female     244  55 133.0  80.0       No 25.01       No     No      70
    ## 870  Female     241  59 120.0  73.0       No 23.76       No     No      88
    ## 871    Male     207  55 122.0  82.0       No 24.15       No     No      74
    ## 872  Female     214  43 121.5  78.0       No 26.11       No     No      82
    ## 873  Female     225  46 116.0  79.0       No 29.21       No     No      70
    ## 874    Male     214  43 132.0  84.0       No 26.77       No     No     102
    ## 875  Female     304  46 131.5  78.5       No 21.02       No     No     112
    ## 876    Male     238  49 118.0  66.0      Yes 26.91       No     No      82
    ## 877  Female     241  42 118.5  80.5       No 32.36       No     No      75
    ## 878  Female     255  61 114.0  70.5       No 24.79       No     No     107
    ## 879  Female     212  36 102.0  69.0       No 33.36       No     No      71
    ## 880    Male     207  39 131.0  82.0      Yes 24.25       No     No      83
    ## 881  Female     222  48 119.0  85.0       No 30.46       No     No      80
    ## 882    Male     220  39 107.0  76.5       No 22.45       No     No      78
    ## 883  Female     240  59 122.5  67.5      Yes 25.40       No     No      81
    ## 884  Female     231  48 131.0  80.0       No 21.14       No     No      85
    ## 885    Male     268  44 114.0  83.0       No 31.16       No     No      76
    ## 886    Male     256  42 111.0  62.0      Yes 21.67       No     No      74
    ## 887  Female     195  44 118.0  86.0       No 23.09       No     No      75
    ## 888  Female     247  52 146.0  88.0       No 27.26       No     No      63
    ## 889    Male     166  60 141.0  81.0      Yes 19.42       No     No     101
    ## 890  Female     248  35 107.0  73.0       No 20.64       No     No      80
    ## 891    Male     165  40 128.0  83.0      Yes 24.71       No     No      60
    ## 892  Female     283  64 163.0  85.0      Yes 21.17       No     No      68
    ## 893  Female     261  53 136.0  99.0      Yes 21.02       No     No      94
    ## 894  Female     264  59 244.0 124.0       No 19.61       No    Yes     120
    ## 895  Female     273  55 134.0  92.0       No 32.17       No     No      67
    ## 896  Female     234  58 113.0  77.0      Yes 20.68       No     No      67
    ## 897    Male     197  56 113.5  74.0      Yes 21.03       No     No      81
    ## 898  Female     202  40 104.0  76.0      Yes 19.93       No     No      62
    ## 899    Male     163  42 135.0  82.0       No 25.75       No     No      77
    ## 900    Male     287  61 114.0  81.0       No 27.53       No     No      68
    ## 901    Male     247  48 131.0  79.0       No 22.12       No     No      83
    ## 904    Male     251  65 127.5  75.0       No 26.46       No     No      74
    ## 905    Male     205  52 122.0  73.0       No 22.73       No     No      85
    ## 906    Male     225  39 128.0  86.5      Yes 25.13       No     No     100
    ## 907  Female     250  55 127.5  83.5       No 30.61       No     No      73
    ## 908  Female     255  64 153.0  75.0       No 23.39       No     No      74
    ## 909    Male     239  51 168.0 102.0      Yes 30.38       No     No      68
    ## 910    Male     264  49 127.5  81.0       No 25.16       No     No      70
    ## 911    Male     232  53 134.0  91.0      Yes 25.13       No     No      75
    ## 912  Female     232  47 113.5  73.0       No 28.78       No     No      77
    ## 913    Male     225  52 156.0  98.0      Yes 30.93       No     No     100
    ## 915  Female     221  50 112.0  69.0      Yes 24.07       No     No      79
    ## 916    Male     242  58 164.0  85.5      Yes 18.84       No     No     106
    ## 917  Female     194  64 176.0  97.0       No 33.19       No     No      89
    ## 918    Male     238  37 121.0  80.0       No 28.95       No     No      67
    ## 919    Male     206  57 141.0  83.0      Yes 27.98       No     No      69
    ## 920    Male     261  63 148.0  93.0       No 30.03       No     No      75
    ## 921    Male     273  64 123.0  73.5       No 22.30       No     No      84
    ## 922  Female     290  66 152.5  90.0       No 23.63       No    Yes      76
    ## 923    Male     290  35 108.0  72.0       No 22.91       No     No      84
    ## 924    Male     180  38 111.0  61.0      Yes 21.51       No     No      75
    ## 925  Female     233  39 126.0  85.0      Yes 22.89       No     No      87
    ## 926  Female     199  42 141.0  92.0      Yes 43.69       No     No      60
    ## 927    Male     255  45 130.0  82.5      Yes 28.56       No     No      94
    ## 929  Female     200  57 108.0  77.0      Yes 18.55       No     No      87
    ## 930  Female     244  62 168.0 102.0      Yes 26.39       No     No     105
    ## 931    Male     255  58 146.0  89.0      Yes 27.47       No     No      73
    ## 932  Female     257  51 128.0  77.0       No 24.94       No     No      88
    ## 933    Male     176  59 134.5  87.0      Yes 31.76       No     No      93
    ## 934  Female     217  44 124.5  82.0      Yes 22.36       No     No      68
    ## 935  Female     282  60 213.0  94.5       No 28.58      Yes     No      78
    ## 936  Female     318  54 115.0  81.0       No 25.84       No     No      76
    ## 938    Male     220  53 127.0  76.0       No 24.27       No     No      74
    ## 939  Female     201  43 129.0  92.0      Yes 24.54       No     No      63
    ## 940    Male     206  64 126.0  82.0      Yes 24.35       No     No      97
    ## 941  Female     259  61 167.5  91.5       No 29.53       No     No      85
    ## 942    Male     224  39 108.0  66.0       No 28.57       No     No      97
    ## 943  Female     154  35 125.0  75.0      Yes 23.10       No     No      75
    ## 944  Female     273  66 153.0  94.0       No 25.27       No     No      76
    ## 945  Female     303  53 117.0  71.0       No 22.01       No     No      85
    ## 946    Male     229  40 152.0 103.0      Yes 32.73       No     No      72
    ## 947  Female     253  46 118.0  74.0      Yes 26.42       No     No      64
    ## 948  Female     214  55 154.0  88.5       No 28.20       No     No      90
    ## 950    Male     219  46 114.0  84.0      Yes 21.89       No     No      60
    ## 952    Male     152  49 120.0  90.0       No 23.03       No     No      93
    ## 953  Female     194  42 111.0  67.5      Yes 21.34       No     No      47
    ## 954  Female     353  60 116.0  82.0      Yes 22.66       No     No      71
    ## 955  Female     298  56 115.0  80.0       No 31.11       No     No      77
    ## 957    Male     253  39 159.0 115.0      Yes 32.66       No     No      74
    ## 958    Male     199  51 134.0  84.0       No 23.04       No     No     120
    ## 959  Female     318  62 206.0  98.0       No 27.23       No     No      87
    ## 960  Female     196  45 123.0  71.0      Yes 20.56       No     No      76
    ## 961    Male     360  61 157.0  99.0      Yes 28.74       No     No      73
    ## 962    Male     230  55 142.0  74.0      Yes 23.65       No     No      82
    ## 963    Male     243  58 106.0  70.5      Yes 23.72       No     No      80
    ## 964  Female     247  53 139.0  88.0      Yes 23.71       No     No      53
    ## 965  Female     206  62 242.0 141.0      Yes 39.86       No    Yes      94
    ## 966  Female     229  37 111.0  70.0      Yes 20.24       No     No      70
    ## 967    Male     280  36 151.0  96.0      Yes 25.35       No     No      94
    ## 968  Female     335  60 199.0  83.0       No 24.61       No     No      90
    ## 969  Female     219  53 184.0 109.0       No 22.73       No     No      73
    ## 970    Male     168  62 129.5  87.0      Yes 20.56       No     No      80
    ## 971  Female     250  59 127.5  80.0       No 29.16       No     No     108
    ## 973  Female     177  62 171.5  83.0       No 30.98       No     No      81
    ## 974    Male     260  43 129.0  90.0       No 25.29       No     No      62
    ## 975  Female     232  47 133.0  86.0      Yes 20.15       No     No      74
    ## 977  Female     284  53 167.5 102.5       No 31.50       No     No      87
    ## 978  Female     247  47 160.0  85.0      Yes 27.05       No     No      77
    ## 981    Male     150  48 127.5  75.0      Yes 26.60       No     No     112
    ## 982  Female     202  43 115.5  68.0       No 23.33       No     No      73
    ## 983  Female     215  47 128.0  85.0      Yes 20.89       No     No      90
    ## 984  Female     229  43 118.0  77.0       No 25.32       No     No     103
    ## 985    Male     285  56 198.0 107.0      Yes 24.87       No     No      97
    ## 986  Female     220  34 117.5  67.5      Yes 20.79       No     No      86
    ## 987  Female     274  47 127.0  86.0       No 21.93       No     No      83
    ## 989    Male     234  66 114.5  62.5      Yes 28.62      Yes     No     216
    ## 990  Female     170  39 137.5  77.5      Yes 27.35       No     No      70
    ## 991  Female     268  58 151.0  98.0       No 20.34       No     No      60
    ## 992    Male     243  66 172.5  88.0      Yes 25.98       No     No      90
    ## 993    Male     239  61 151.0  89.5       No 25.02       No    Yes      88
    ## 994    Male     257  57 158.5 107.0      Yes 27.10       No     No      67
    ## 995  Female     252  63 118.0  84.0       No 25.31       No     No      82
    ## 996  Female     285  53 160.0  97.0       No 31.31       No     No      65
    ## 997    Male     209  38 137.0  82.5       No 26.69       No     No      84
    ## 998    Male     262  57 131.5  92.0      Yes 28.30       No     No      78
    ## 999  Female     261  49 123.5  84.0       No 20.94       No     No      75
    ## 1000 Female     202  47 140.0  82.0       No 22.88       No     No      95
    ## 1001 Female     160  36  98.0  66.0      Yes 25.07       No     No      73
    ## 1002 Female     218  42 109.5  67.0      Yes 23.48       No     No      71
    ## 1003 Female     284  39 115.5  65.5      Yes 20.39       No     No      78
    ## 1004   Male     265  52 106.0  79.0      Yes 26.48       No     No     163
    ## 1006 Female     210  61 122.0  81.0       No 22.48       No     No      97
    ## 1007 Female     266  54 137.0  88.0       No 29.76       No     No      80
    ## 1008 Female     256  61 160.0 109.0       No 42.53       No     No      79
    ## 1009   Male     199  52 134.0  98.0      Yes 27.78       No     No      89
    ## 1010 Female     171  38 125.5  86.0       No 26.57       No     No      71
    ## 1012 Female     205  37 120.5  67.5      Yes 22.89       No     No     113
    ## 1013   Male     225  52 131.0  74.0      Yes 24.54       No     No      77
    ## 1014   Male     272  50 159.0 109.0      Yes 31.27       No     No      68
    ## 1015   Male     269  60 170.0 100.0      Yes 29.59       No    Yes      83
    ## 1016   Male     160  39 124.0  90.0       No 26.82       No     No      67
    ## 1017   Male     158  57 154.0 100.0      Yes 24.07       No     No      70
    ## 1018 Female     239  64 132.0  75.0       No 20.67       No     No      80
    ## 1020 Female     258  45 111.0  72.0      Yes 26.24       No     No      65
    ## 1021 Female     201  67 110.0  77.0       No 24.53       No     No      77
    ## 1022 Female     340  40 149.0  81.0       No 28.46       No     No      99
    ## 1023 Female     192  48 135.0  82.5      Yes 32.67       No     No      69
    ## 1024 Female     300  53 127.0  89.0      Yes 25.46       No     No      70
    ## 1025 Female     177  50 121.0  67.0      Yes 22.02       No     No      77
    ## 1026   Male     194  47 115.0  82.0       No 28.23       No     No      76
    ## 1027   Male     266  42 139.0  88.0       No 23.61       No     No      78
    ## 1028   Male     287  53 122.0  80.0       No 26.26       No     No      85
    ## 1029   Male     247  39 113.0  80.0       No 28.59       No     No      82
    ## 1030   Male     199  39 124.0  86.0      Yes 23.39       No     No      72
    ## 1031 Female     205  64 139.0  92.0       No 32.32       No     No      90
    ## 1032 Female     310  61 118.0  77.5      Yes 24.03       No     No      70
    ## 1033 Female     270  46 122.0  76.0       No 21.35       No     No      88
    ## 1035   Male     199  49 107.5  71.0      Yes 26.62       No     No      68
    ## 1036 Female     255  55 125.0  85.0       No 22.89       No     No      81
    ## 1039 Female     168  41 102.0  64.0       No 23.64       No     No      75
    ## 1040   Male     216  46 124.0  85.0       No 29.91       No     No     103
    ## 1042 Female     263  64 206.0 104.0       No 26.15       No     No      91
    ## 1043 Female     175  43 117.0  67.0       No 22.36       No     No      70
    ## 1044   Male     240  50 145.0  94.0       No 28.86       No     No      68
    ## 1045 Female     220  38 107.0  73.5       No 23.09       No     No      80
    ## 1046 Female     310  56 142.0  94.0       No 31.10       No     No      65
    ## 1047   Male     186  53 167.0  96.5      Yes 25.09       No     No     113
    ## 1048 Female     220  51 142.0  82.5      Yes 21.02       No     No      78
    ## 1049   Male     232  42 111.5  70.0      Yes 28.30       No     No      80
    ## 1050   Male     240  58 148.0  81.0      Yes 25.67       No     No      78
    ## 1051 Female     187  54 133.0  88.0      Yes 31.82       No     No      77
    ## 1052   Male     213  53 100.0  71.0       No 23.85       No     No      75
    ## 1053   Male     210  58 132.0  86.0       No 28.92       No     No      74
    ## 1054   Male     267  60 139.0  84.0      Yes 28.76       No     No     112
    ## 1055   Male     332  62 119.5  74.0       No 28.50       No     No      92
    ## 1056   Male     236  59 127.0  83.0       No 26.53       No     No      86
    ## 1057 Female     232  53 147.0  71.5       No 25.45       No     No      74
    ## 1059   Male     242  35 136.5  95.0       No 24.43       No     No      88
    ## 1060 Female     224  50 149.0  90.0       No 29.94       No     No      85
    ## 1061   Male     167  36 155.0  74.0       No 19.42       No     No      81
    ## 1062 Female     277  57 133.0  84.0       No 36.21       No     No      74
    ## 1063   Male     227  34 131.5  84.0       No 25.41       No     No      87
    ## 1064   Male     217  46 117.5  77.5       No 32.40       No     No      83
    ## 1065   Male     200  42 131.0  88.0      Yes 25.09       No     No      88
    ## 1068 Female     197  46 105.5  67.0       No 23.14       No     No      69
    ## 1069   Male     279  45 138.0  86.0      Yes 30.63      Yes     No     144
    ## 1070 Female     236  60 126.0  84.0       No 20.34       No     No      76
    ## 1071   Male     218  45 145.0  90.0      Yes 26.65       No     No      70
    ## 1072 Female     187  44 122.0  83.0       No 30.40       No     No      75
    ## 1073   Male     222  39 141.5  91.0      Yes 27.06       No     No      73
    ## 1074 Female     183  42 111.0  71.0       No 19.78       No     No      81
    ## 1075 Female     264  62 142.0  90.0      Yes 31.78       No     No      97
    ## 1076   Male     270  48 131.0  88.0       No 27.13       No     No      55
    ## 1078   Male     207  45 117.0  80.0      Yes 20.17       No     No      59
    ## 1079   Male     192  58 143.0  98.0      Yes 29.01       No     No      68
    ## 1080 Female     285  56 165.0 115.0      Yes 24.25       No    Yes     116
    ## 1081 Female     250  55 161.5  95.0       No 27.76       No     No      83
    ## 1082   Male     211  57 122.0  66.5      Yes 20.19       No     No      57
    ## 1084   Male     201  61 164.5  93.5       No 27.73       No     No      95
    ## 1085 Female     280  49 120.0  80.0      Yes 22.33       No     No      75
    ## 1086   Male     164  42 141.5  98.0       No 32.52       No     No      76
    ## 1087   Male     220  50 124.0  82.0       No 24.54       No     No      83
    ## 1088   Male     270  42 112.0  77.0      Yes 24.77       No     No      85
    ## 1089   Male     346  46 137.0  97.0      Yes 29.11       No     No      82
    ## 1090   Male     266  53 163.0 105.0      Yes 28.04       No     No      79
    ## 1091 Female     236  59 123.0  76.5      Yes 30.67       No     No     107
    ## 1092   Male     169  44 116.0  62.0       No 19.44       No     No      87
    ## 1093 Female     237  54 171.5 105.5       No 34.25       No     No     104
    ## 1094   Male     250  58 150.0  97.0      Yes 32.00       No     No      65
    ## 1095   Male     245  67 126.0  68.0      Yes 29.04       No     No      94
    ## 1096   Male     271  44 132.5  90.0       No 24.06       No     No      95
    ## 1097   Male     176  41 146.0  88.0      Yes 24.04       No     No      83
    ## 1098 Female     240  52 177.0 103.5       No 24.39       No     No      75
    ## 1099 Female     187  40 105.0  74.0       No 23.26       No     No      69
    ## 1100 Female     252  55 108.5  63.5      Yes 25.23       No     No     121
    ## 1101 Female     140  41 110.0  60.0      Yes 23.38       No     No      82
    ## 1102   Male     185  40 102.0  72.0       No 24.08       No     No      83
    ## 1105   Male     324  41 129.5  92.5       No 34.99       No     No     103
    ## 1106   Male     309  52 152.5  91.0      Yes 30.12       No     No     103
    ## 1107 Female     273  48 132.0  85.0       No 26.03       No     No      78
    ## 1108 Female     245  52 131.0  80.0       No 32.04       No     No      81
    ## 1109   Male     229  53 146.5  82.0      Yes 27.80      Yes     No     172
    ## 1110   Male     200  52 114.0  77.0       No 28.28       No     No      84
    ## 1112 Female     276  65 124.0  70.0       No 25.61       No     No      75
    ## 1113 Female     233  64 160.0  93.0       No 27.60       No    Yes     111
    ## 1115   Male     181  48 153.0  93.0       No 29.34       No     No      88
    ## 1116 Female     170  45 109.5  69.0      Yes 17.38       No     No      66
    ## 1118 Female     190  58 132.0  67.0       No 23.08       No     No      70
    ## 1119 Female     225  54 131.0  79.0      Yes 25.91       No     No      62
    ## 1120   Male     183  38 128.0  80.5       No 33.32       No     No      72
    ## 1121 Female     273  56 125.0  83.0       No 24.48       No    Yes      66
    ## 1123 Female     263  67 201.0  93.0       No 30.04       No    Yes      78
    ## 1124   Male     290  55 120.5  84.0      Yes 25.05       No     No      90
    ## 1125 Female     273  59 145.0  90.0       No 23.94       No     No      82
    ## 1126 Female     265  49 144.0  86.0      Yes 25.57       No     No      68
    ## 1127 Female     279  54 127.0  70.0       No 23.48       No     No      79
    ## 1128   Male     232  45 122.5  82.5       No 27.55       No     No      59
    ## 1129 Female     346  55 131.0  81.0      Yes 22.69       No     No      77
    ## 1130 Female     253  46 118.0  79.0       No 26.61       No     No      83
    ## 1131 Female     216  40 112.5  76.5      Yes 27.22       No     No      77
    ## 1132   Male     220  44 135.0  91.0      Yes 27.23       No     No      88
    ## 1134 Female     239  40 118.0  78.0      Yes 23.48       No     No      75
    ## 1135 Female     261  37 123.0  75.0       No 26.72       No     No      85
    ## 1136   Male     210  47 112.0  66.0      Yes 24.58       No     No      84
    ## 1137   Male     253  47 137.0  87.0      Yes 24.50       No     No      81
    ## 1138   Male     263  43 114.0  81.0      Yes 25.68       No     No      74
    ## 1139 Female     265  53 137.0  76.0       No 25.46       No     No      84
    ## 1140   Male     219  51 125.0  71.0      Yes 21.19       No     No      75
    ## 1141 Female     246  52 113.5  66.5      Yes 19.47       No     No      60
    ## 1142   Male     200  52 113.0  76.5       No 27.22       No     No      73
    ## 1143   Male     220  50 114.0  78.0       No 26.26       No     No      83
    ## 1144 Female     279  37 110.0  72.5      Yes 24.89       No     No      70
    ## 1145 Female     274  50 125.0  87.0       No 21.67       No     No      73
    ## 1147   Male     229  43 131.0  87.0      Yes 23.31       No     No      74
    ## 1148 Female     352  60 149.0  73.0      Yes 25.96       No     No      79
    ## 1151   Male     230  53 137.0  99.0      Yes 25.77       No     No      79
    ## 1152 Female     283  55 213.0 103.0       No 28.62       No    Yes      76
    ## 1153 Female     209  42 105.0  65.0      Yes 23.80       No     No      64
    ## 1154   Male     227  61 130.0  77.5       No 26.18       No     No      91
    ## 1155 Female     159  40 145.0  90.0       No 20.33       No     No      83
    ## 1156 Female     201  59 148.5  90.0       No 25.85       No    Yes      83
    ## 1157 Female     600  52 159.5  94.0       No 28.27      Yes     No     140
    ## 1158 Female     241  63 143.5  89.0       No 26.45       No     No      80
    ## 1159 Female     238  66 140.0  80.0       No 26.69       No     No      83
    ## 1160   Male     220  53 119.0  75.0      Yes 26.28       No     No      67
    ## 1162 Female     175  39 105.5  64.5      Yes 25.83       No     No      72
    ## 1163   Male     272  58 127.5  79.5      Yes 26.37       No     No      66
    ## 1166 Female     267  56 114.0  71.0       No 24.81       No     No      90
    ## 1167   Male     241  54 101.0  72.0      Yes 28.35       No     No      70
    ## 1169   Male     251  64 133.0  72.5       No 24.28      Yes     No      86
    ## 1170 Female     210  44 138.0  92.0       No 23.13       No     No      92
    ## 1171 Female     259  46 173.0 102.0       No 27.22       No     No      75
    ## 1172   Male     202  40 112.5  64.5      Yes 22.85       No     No     103
    ## 1173 Female     253  42 109.0  74.0      Yes 24.38       No     No      60
    ## 1174 Female     236  63 189.0 103.0      Yes 27.91       No     No      74
    ## 1176 Female     292  61 141.5  88.5       No 23.95       No     No      70
    ## 1177 Female     215  52 159.0  64.0       No 24.56       No     No     124
    ## 1178 Female     206  51 146.0  77.0       No 23.58       No     No      90
    ## 1179 Female     293  64 140.0  84.0       No 34.56       No     No      76
    ## 1180 Female     262  64 122.0  87.5       No 24.77       No     No      85
    ## 1181   Male     282  40 130.0  80.0      Yes 23.90       No     No      68
    ## 1182 Female     245  44 125.0  80.5      Yes 24.58       No     No      80
    ## 1183 Female     192  45 132.0  79.0      Yes 24.53       No     No     112
    ## 1184   Male     160  58 159.0  89.0       No 27.71       No     No      85
    ## 1186   Male     217  44 144.0  99.0      Yes 25.16       No     No      60
    ## 1187   Male     249  50 133.0  88.0      Yes 28.50       No     No      75
    ## 1188   Male     256  56 147.0  96.0       No 30.42       No     No      75
    ## 1190 Female     251  58 141.0  73.5       No 20.97       No     No      96
    ## 1191 Female     315  63 156.0  90.0       No 25.92       No     No      74
    ## 1192 Female     224  57 174.0 112.0       No 22.73       No     No      82
    ## 1193 Female     196  64 150.0  84.0       No 25.98       No     No      93
    ## 1194   Male     225  57 140.0  88.0      Yes 24.71       No     No      78
    ## 1195   Male     270  55 129.0  82.0      Yes 27.63       No     No      60
    ## 1196 Female     265  55 130.0  84.0       No 29.66       No     No      69
    ## 1197   Male     193  41 100.0  60.0      Yes 29.69       No     No      69
    ## 1198 Female     392  46 113.0  68.0      Yes 23.35       No     No      63
    ## 1199 Female     239  38 128.0  84.5       No 33.49       No     No      80
    ## 1200 Female     261  56 145.0  77.0       No 26.67       No     No      95
    ## 1201 Female     242  67 172.0  84.0      Yes 19.81       No     No     111
    ## 1202 Female     222  45 167.0 107.0       No 25.56       No     No      88
    ## 1204 Female     242  50 116.0  69.0      Yes 21.65       No     No      73
    ## 1205 Female     257  41 122.0  73.0       No 24.17       No     No     104
    ## 1206 Female     291  46 107.5  65.0      Yes 24.10       No     No      78
    ## 1207 Female     229  60 144.0  91.0       No 24.96       No     No      82
    ## 1208 Female     283  63 188.0 107.0       No 33.40       No    Yes      78
    ## 1209 Female     223  41 119.0  73.0      Yes 24.22       No     No      77
    ## 1211   Male     214  61 100.0  65.0       No 30.18       No     No      66
    ## 1212   Male     239  68 130.0  80.0       No 23.25       No     No      95
    ## 1213 Female     250  45 116.0  79.0       No 28.59       No     No      87
    ## 1214 Female     334  52 147.0  86.0      Yes 29.01      Yes     No      63
    ## 1215   Male     322  56 140.0  90.0       No 29.47       No     No      64
    ## 1216 Female     291  53 114.0  81.0       No 26.21       No     No      79
    ## 1217 Female     214  51 153.5 103.5       No 23.45       No     No      71
    ## 1218 Female     308  65 127.5  67.5       No 28.82      Yes     No     366
    ## 1219   Male     200  47 126.0  86.0      Yes 26.84       No     No      96
    ## 1220   Male     256  67 138.0  76.0       No 22.81       No     No     100
    ## 1221   Male     221  49 175.0 107.5       No 25.97       No     No      78
    ## 1223 Female     185  36 123.0  69.0       No 18.98       No     No      75
    ## 1224 Female     304  62 165.0  91.0       No 26.97       No     No      81
    ## 1226 Female     268  53 131.0  83.0       No 21.03       No     No      78
    ## 1227   Male     248  36 135.0  94.5      Yes 36.52       No     No      85
    ## 1229   Male     221  45 131.0  84.0      Yes 28.58       No     No      72
    ## 1231 Female     242  47 145.0  87.5       No 22.01       No     No      73
    ## 1233   Male     170  38 130.0  94.0      Yes 23.90       No     No      75
    ## 1234 Female     233  57 184.0 106.0       No 38.88       No     No      40
    ## 1235 Female     315  64 135.0  80.0      Yes 25.23       No     No      89
    ## 1236   Male     205  46 118.0  76.5       No 23.48       No     No      77
    ## 1237   Male     267  49 160.5 109.0       No 28.33       No     No      75
    ## 1238   Male     261  50 114.0  64.0      Yes 22.32       No     No      71
    ## 1239 Female     306  63 177.0  96.0       No 32.51       No     No     126
    ## 1240 Female     169  48 243.0 142.5       No 28.49       No     No      77
    ## 1242   Male     260  60 178.0 103.0      Yes 24.62       No     No      79
    ## 1243 Female     242  41 139.0  80.0      Yes 19.68       No     No      60
    ## 1244 Female     280  64 127.0  77.0       No 30.39       No     No      78
    ## 1245 Female     231  62 184.0  90.0       No 26.03       No     No      75
    ## 1246   Male     270  62 145.5  87.5      Yes 23.88       No     No      67
    ## 1247 Female     216  66 133.0  87.0      Yes 30.06       No     No      91
    ## 1248 Female     320  59 187.5  85.5       No 25.33       No    Yes      68
    ## 1249 Female     214  36 119.0  76.0      Yes 21.67       No     No      75
    ## 1250   Male     304  40 125.0  86.0      Yes 30.07       No     No      84
    ## 1251   Male     248  53 165.0  96.0      Yes 24.45       No     No      71
    ## 1252 Female     246  49 138.0  92.0       No 28.23       No     No      91
    ## 1253   Male     240  48 126.0  84.0      Yes 24.35       No     No      73
    ## 1255   Male     253  64 178.0 106.0      Yes 24.68       No     No      76
    ## 1256   Male     256  46 138.0 105.0       No 26.97       No     No     100
    ## 1257 Female     245  58 123.0  85.5       No 24.49       No     No      76
    ## 1259 Female     205  41 110.0  69.0       No 25.99       No     No      67
    ## 1260 Female     199  40 122.0  82.0      Yes 22.16       No     No      77
    ## 1261 Female     304  42 119.0  76.0      Yes 32.52       No     No      80
    ## 1262 Female     281  53 118.0  72.0       No 19.83       No     No      64
    ## 1263   Male     206  52 110.0  78.0       No 33.03       No     No      62
    ## 1265 Female     309  49 155.0  85.0       No 23.06       No     No      63
    ## 1266 Female     234  52 111.0  81.0      Yes 22.35       No     No      77
    ## 1268   Male     305  46 150.0  88.0       No 26.82       No     No      75
    ## 1269 Female     208  40 116.0  75.0       No 18.52       No     No      82
    ## 1270   Male     299  66 117.0  76.0       No 27.44       No     No      83
    ## 1271   Male     171  65 173.0  84.0      Yes 28.70       No     No      70
    ## 1272   Male     210  35  99.0  67.0      Yes 22.39       No     No      75
    ## 1273 Female     271  42 115.0  77.0       No 28.68       No     No      82
    ## 1274   Male     264  58 181.0  90.0       No 24.49       No     No      71
    ## 1275 Female     209  43 115.0  75.0       No 27.99       No     No      90
    ## 1276   Male     248  44 174.0 110.0       No 31.74       No     No     100
    ## 1277 Female     155  38  96.0  61.0      Yes 24.19       No     No      68
    ## 1278 Female     149  39 122.0  72.0      Yes 21.30       No     No      75
    ## 1279   Male     255  41 120.5  85.5       No 30.85       No     No      79
    ## 1280   Male     137  40 127.0  82.0       No 27.04       No     No      71
    ## 1281 Female     160  38 102.5  67.5      Yes 21.16       No     No      68
    ## 1283 Female     226  41 130.0  80.0      Yes 25.25       No     No      73
    ## 1284 Female     197  42 124.0  81.0      Yes 21.50       No     No      63
    ## 1285   Male     279  58 180.0 109.5      Yes 26.04       No     No      82
    ## 1286   Male     195  42 112.0  74.5       No 23.37       No     No     110
    ## 1287 Female     231  47 133.5  76.0      Yes 25.77       No     No      73
    ## 1288   Male     229  59 100.5  66.0       No 25.18       No     No      81
    ## 1290   Male     242  55 120.0  86.0      Yes 26.77       No     No      73
    ## 1291   Male     230  68 124.0  70.0       No 22.85       No     No      60
    ## 1292   Male     183  61 156.0  92.0      Yes 24.69       No     No      79
    ## 1293 Female     208  36 156.5 105.0       No 33.82      Yes     No     186
    ## 1295 Female     261  45 140.0  88.0      Yes 21.44       No     No      78
    ## 1296   Male     246  51 128.0  69.0      Yes 27.57       No     No      72
    ## 1297   Male     172  45 119.0  84.0       No 28.25       No     No      98
    ## 1298 Female     216  35 130.0  68.0       No 25.94       No     No      90
    ## 1299 Female     222  49 118.0  78.0       No 21.18       No     No      73
    ## 1300 Female     298  50 156.0  90.0      Yes 24.24       No     No     100
    ## 1301   Male     208  40 148.0 100.0       No 32.84       No     No     102
    ## 1302 Female     213  43 113.0  77.0      Yes 29.34       No     No      73
    ## 1303 Female     266  49 159.0  88.0      Yes 20.66       No     No      84
    ## 1304   Male     220  47 127.0  93.0      Yes 30.70       No     No      57
    ## 1306 Female     210  58 143.0 101.0       No 31.34       No     No      80
    ## 1307   Male     252  38 125.0  92.0       No 24.72       No     No      69
    ## 1308   Male     290  56 144.0  86.0      Yes 22.26       No     No      90
    ## 1310 Female     285  47 122.0  70.0      Yes 23.48       No     No      82
    ## 1311   Male     175  60 129.0  89.0       No 22.16       No     No      75
    ## 1312 Female     225  45 108.5  71.5       No 25.74       No     No      80
    ## 1313   Male     213  58 162.0  99.0      Yes 28.30       No     No      70
    ## 1314   Male     252  63 135.5  80.0       No 28.78       No     No      79
    ## 1316   Male     318  52 144.0  85.0       No 27.66       No     No      64
    ## 1317 Female     230  47 123.0  71.0       No 26.98       No     No      73
    ## 1318 Female     248  46 128.0  76.0      Yes 28.87       No    Yes      77
    ## 1319 Female     270  49 126.5  67.5       No 26.56       No     No      77
    ## 1320   Male     245  36 131.5  89.0       No 26.33       No     No      62
    ## 1321   Male     196  58 120.0  74.0      Yes 20.12       No     No      73
    ## 1322   Male     220  58 143.0 104.0       No 29.85       No     No      87
    ## 1323 Female     190  46 128.0  74.0       No 23.01       No     No      78
    ## 1324   Male     218  61 160.0  96.0       No 28.89      Yes     No     223
    ## 1325 Female     184  59 122.0  74.0      Yes 24.66       No     No      67
    ## 1326 Female     205  40  97.0  63.0       No 26.56       No     No      80
    ## 1327 Female     272  54 146.0  95.0       No 23.66       No     No      80
    ## 1328 Female     186  56 155.0 102.0       No 24.38       No     No      75
    ## 1329   Male     218  45 128.0  90.0       No 32.15       No     No      77
    ## 1330 Female     234  51 114.0  85.0       No 28.68       No     No      84
    ## 1332 Female     270  59 175.0  95.0      Yes 29.69       No     No      76
    ## 1333   Male     154  40 125.0  94.0      Yes 29.29       No     No      84
    ## 1334   Male     253  53 121.0  85.5      Yes 28.52       No     No      68
    ## 1335   Male     214  64 116.0  77.0      Yes 22.48       No     No      71
    ## 1336   Male     226  48 117.5  80.0      Yes 26.18       No     No      66
    ## 1337 Female     187  41 108.0  64.0       No 23.63       No     No      70
    ## 1338   Male     241  43 132.5  87.5       No 32.02       No     No      83
    ## 1339   Male     195  56 108.0  70.0      Yes 22.92       No     No     117
    ## 1341 Female     233  40 108.5  75.0       No 28.30       No     No      73
    ## 1343 Female     218  45 110.0  70.0       No 20.24       No     No      78
    ## 1345 Female     238  60 140.0  85.0       No 28.41       No     No      88
    ## 1346 Female     244  40 110.0  73.0      Yes 21.84       No     No      67
    ## 1347 Female     222  55 103.0  61.0      Yes 23.18       No     No      75
    ## 1348   Male     225  42 122.5  80.0      Yes 25.54       No     No      90
    ## 1349   Male     197  47 115.0  81.0       No 22.41       No     No      66
    ## 1350 Female     240  50 163.0 105.0       No 31.37       No     No      75
    ## 1351 Female     222  45  95.0  58.0      Yes 21.68       No     No      77
    ## 1352   Male     310  39 134.0  90.0      Yes 35.11       No     No      88
    ## 1354   Male     175  40 173.0  59.0       No 27.99       No     No      75
    ## 1355   Male     312  52 148.0  99.0      Yes 26.73       No     No      65
    ## 1356   Male     253  50 131.0  87.0       No 26.54       No     No      76
    ## 1358 Female     226  63 172.5  98.0       No 26.47       No     No      81
    ## 1359   Male     180  55 174.0 100.0       No 26.83      Yes     No      78
    ## 1360   Male     259  46 125.0  89.0      Yes 24.80       No     No      57
    ## 1361   Male     216  62 126.0  82.0       No 21.18       No     No      69
    ## 1362   Male     260  41 120.0  72.5      Yes 26.36       No     No      88
    ## 1363   Male     233  38 137.0  75.0      Yes 20.55       No     No     107
    ## 1364   Male     280  59 164.0  81.0      Yes 29.76       No     No      68
    ## 1365 Female     170  34 121.0  74.0       No 20.82       No     No      83
    ## 1366   Male     212  50 132.0  87.0      Yes 25.90       No     No      83
    ## 1368   Male     206  59 115.0  70.0      Yes 24.79       No     No      76
    ## 1369   Male     230  44 128.0  87.0      Yes 26.02       No     No      73
    ## 1370 Female     158  33 108.0  67.0       No 19.84       No     No      69
    ## 1372   Male     183  60 114.0  70.0       No 23.56       No     No      59
    ## 1373 Female     211  55 136.0  70.0       No 34.40       No     No      83
    ## 1375   Male     238  50 127.5  85.0      Yes 27.82       No     No      69
    ## 1376 Female     186  35 106.0  78.0      Yes 24.73       No     No      70
    ## 1377 Female     239  57 127.0  81.0       No 21.85       No     No      87
    ## 1378 Female     277  60 126.0  70.0       No 25.13       No     No     100
    ## 1379 Female     220  55 158.0  97.0       No 25.11       No     No      86
    ## 1380 Female     235  53 154.0  98.0      Yes 26.91       No     No      65
    ## 1381   Male     207  43 130.0  86.0      Yes 24.96       No     No     103
    ## 1382   Male     256  56 124.0  78.0      Yes 26.67       No     No      80
    ## 1383 Female     230  36 112.5  73.5       No 22.36       No     No      65
    ## 1384   Male     199  38 135.0  90.0       No 32.19       No     No      85
    ## 1386   Male     271  49 140.0 108.0       No 27.66       No     No      77
    ## 1387 Female     250  49 133.5  87.5       No 29.45       No     No      84
    ## 1388 Female     152  42 122.0  76.0      Yes 21.26       No     No      78
    ## 1389   Male     300  35 120.0  84.0      Yes 24.69       No     No      99
    ## 1390 Female     287  38 113.0  80.0       No 24.56       No     No      75
    ## 1391   Male     237  62 114.0  72.0       No 25.65      Yes     No      89
    ## 1392   Male     214  44 128.0  94.0      Yes 23.51       No     No      66
    ## 1394 Female     201  39 151.0  94.0       No 31.48       No     No      74
    ## 1395 Female     180  39 124.0  83.0       No 22.91       No     No      77
    ## 1396 Female     248  46 115.0  82.0       No 28.92       No     No     100
    ## 1397   Male     240  52 120.0  77.0      Yes 32.27       No     No      62
    ## 1398 Female     304  48 102.0  66.5      Yes 28.90      Yes     No      66
    ## 1399   Male     271  52 121.0  73.0      Yes 21.85       No     No      86
    ## 1400   Male     250  57 127.5  80.0      Yes 29.38       No     No      80
    ## 1401 Female     200  44 111.0  79.0      Yes 27.29       No     No      74
    ## 1402 Female     272  57 112.5  70.0      Yes 23.08       No     No      58
    ## 1403 Female     211  49 103.5  66.5       No 24.17       No     No      87
    ## 1404 Female     230  42 142.5  97.5       No 29.94       No     No      75
    ## 1406 Female     240  54 114.0  80.0       No 30.78       No     No      71
    ## 1407   Male     225  42 110.0  73.0       No 27.67       No     No      65
    ## 1408   Male     185  64 114.0  73.0      Yes 34.53       No     No      97
    ## 1409 Female     260  49 123.0  80.0      Yes 23.10       No     No      65
    ## 1410 Female     292  59 114.0  84.0       No 27.39       No     No      72
    ## 1411   Male     246  37 124.0  83.0      Yes 30.93       No     No      85
    ## 1412   Male     165  33 141.5  95.0       No 26.74       No     No      77
    ## 1414   Male     198  38 119.0  73.0       No 30.27       No     No      70
    ## 1415   Male     204  40 115.0  83.0       No 25.05       No     No      76
    ## 1416   Male     253  66 163.0  86.0       No 24.35       No     No      91
    ## 1417   Male     210  64 120.0  68.0       No 24.77       No     No      77
    ## 1419   Male     217  47 145.0  89.0      Yes 28.88       No     No      75
    ## 1421 Female     190  39 137.0  81.0      Yes 19.57       No     No      85
    ## 1422   Male     205  55 127.0  76.0      Yes 22.24      Yes     No     325
    ## 1423 Female     235  41 144.0  88.0      Yes 24.16       No     No      82
    ## 1424 Female     259  56 138.0  87.0       No 30.73       No     No      75
    ## 1425 Female     289  65 149.5  78.0       No 25.37       No     No      86
    ## 1426 Female     272  53 146.0  89.0       No 25.50       No     No      67
    ## 1428 Female     249  42 101.0  70.0       No 21.74       No     No      60
    ## 1429   Male     194  45 133.0  83.0      Yes 20.41       No     No      76
    ## 1430 Female     340  55 140.0  83.0       No 26.18       No     No      83
    ## 1431   Male     314  53 128.0  77.0      Yes 24.74       No     No      78
    ## 1432   Male     226  62 106.0  67.5      Yes 23.88       No     No      87
    ## 1433   Male     225  54 113.5  74.0       No 25.63       No     No      80
    ## 1435 Female     232  43 138.0  88.0       No 22.53       No     No      96
    ## 1437   Male     278  49 135.0  76.0      Yes 29.03       No     No      75
    ## 1438 Female     305  66 138.0  86.0       No 20.74       No     No      62
    ## 1440 Female     317  64 182.5  88.0      Yes 20.52       No     No      79
    ## 1441 Female     292  43 109.0  73.0      Yes 22.87       No     No      93
    ## 1442 Female     175  42 132.0  86.0      Yes 20.53       No     No      88
    ## 1443   Male     295  46 172.5 116.5       No 27.18       No     No      77
    ## 1444   Male     262  43 121.0  79.0       No 24.01       No     No      87
    ## 1445 Female     234  53 126.0  73.0       No 27.60       No     No      90
    ## 1446   Male     240  50 112.0  82.0      Yes 24.39       No     No      62
    ## 1447 Female     251  46 121.0  81.0      Yes 23.05       No     No      84
    ## 1448 Female     358  51 134.0  87.0       No 29.36       No    Yes      87
    ## 1449 Female     239  67 154.0  90.0      Yes 28.56       No     No      90
    ## 1450   Male     185  51 125.0  85.0      Yes 29.43       No     No      72
    ## 1451 Female     217  45 109.0  72.0      Yes 33.65       No     No      68
    ## 1452   Male     225  65 117.0  79.0       No 27.17       No     No      80
    ## 1453   Male     220  44 105.0  70.0      Yes 21.01       No     No      68
    ## 1456   Male     272  44 115.0  76.0      Yes 24.16       No     No      77
    ## 1457 Female     211  48 149.0 100.0       No 30.91       No     No      62
    ## 1458   Male     155  36 123.0  78.0      Yes 22.05       No     No      78
    ## 1459   Male     252  63 149.0  90.0      Yes 33.49       No     No      83
    ## 1461   Male     172  36 122.5  82.5       No 28.53       No     No      75
    ## 1462 Female     244  65 162.0  98.0       No 24.50       No     No      82
    ## 1463 Female     246  48 129.0  86.0      Yes 25.04       No     No      87
    ## 1464 Female     276  63 144.0  90.0       No 21.35       No     No      78
    ## 1465   Male     229  40 137.0  85.0      Yes 35.20       No     No      55
    ## 1466 Female     204  40 122.0  94.0       No 29.86       No     No      82
    ## 1467 Female     170  42 113.0  79.0       No 21.31       No     No      65
    ## 1469 Female     200  66 175.0  67.5       No 27.95      Yes     No      93
    ## 1470 Female     192  63 143.0  87.0       No 23.64       No     No     100
    ## 1471   Male     234  36 133.0  88.0       No 26.78       No     No     100
    ## 1473 Female     161  46 100.0  64.0      Yes 20.66       No     No      60
    ## 1474   Male     251  39 115.0  77.0      Yes 24.01       No     No      98
    ## 1475 Female     255  39 142.0  85.5      Yes 24.89       No     No     108
    ## 1476 Female     345  51 142.0  88.0      Yes 19.05       No     No      73
    ## 1477   Male     180  68 132.0  70.0       No 28.36       No     No      68
    ## 1478   Male     212  48 139.0  86.0      Yes 20.27       No     No      62
    ## 1479 Female     210  57 131.0  85.0       No 26.59       No     No      77
    ## 1480   Male     245  41 146.0  86.0      Yes 24.50       No     No      72
    ## 1481   Male     231  42 123.0  87.0      Yes 21.48       No     No      44
    ## 1482 Female     187  38 118.0  78.0      Yes 30.06       No     No      63
    ## 1483 Female     190  41  95.0  57.0      Yes 20.00       No     No      77
    ## 1484 Female     174  43 158.5 100.5      Yes 35.99       No     No      88
    ## 1485 Female     243  55 142.0  92.0       No 30.24       No     No      85
    ## 1486 Female     160  38  95.0  65.0      Yes 21.99       No     No      77
    ## 1487 Female     250  46 115.0  74.0      Yes 22.70       No     No      69
    ## 1488   Male     263  45 115.0  76.0      Yes 24.94       No     No      78
    ## 1489 Female     221  53 131.0  89.0      Yes 24.09       No     No      95
    ## 1490   Male     271  44 136.0  90.0      Yes 25.24       No     No      64
    ## 1492   Male     260  51 123.0  72.0      Yes 26.83       No     No      65
    ## 1493   Male     213  39 130.0  72.0       No 22.32       No     No      78
    ## 1494 Female     239  51 127.5  77.5       No 26.65       No     No      79
    ## 1495 Female     252  55 130.0  82.0       No 29.17       No     No      85
    ## 1496   Male     201  43 121.0  82.0      Yes 23.84       No     No      91
    ## 1497 Female     261  67 135.0  80.0       No 21.80       No     No      93
    ## 1498 Female     329  55 145.0  82.0       No 23.43       No     No      95
    ## 1499 Female     229  39 125.0  80.0       No 24.10       No     No      58
    ## 1500   Male     254  54 136.5  83.0       No 20.55       No     No      95
    ## 1501   Male     293  41 115.0  77.5      Yes 26.26       No     No      57
    ## 1502 Female     250  48 117.0  81.0       No 27.04       No     No      65
    ## 1503 Female     262  45 116.0  66.0       No 21.56       No     No      76
    ## 1504   Male     256  61 165.0  80.0      Yes 24.12       No     No      97
    ## 1505   Male     239  38 141.0  98.0      Yes 27.80       No     No      85
    ## 1506 Female     186  49 120.0  74.0      Yes 19.39       No     No      69
    ## 1507 Female     232  47 110.0  70.0      Yes 25.86       No     No      82
    ## 1508 Female     260  59 162.5 105.0       No 24.39       No     No      72
    ## 1509 Female     391  60 114.0  64.0       No 24.57       No     No      83
    ## 1511   Male     275  37 127.0  80.0      Yes 27.22       No     No      93
    ## 1513 Female     231  48 115.0  69.0       No 25.48       No     No      77
    ## 1514   Male     219  54 110.0  72.0      Yes 26.05       No     No      86
    ## 1515   Male     215  39 112.0  65.0      Yes 23.60       No     No      78
    ## 1516   Male     279  46 118.0  82.0      Yes 22.78       No     No      74
    ## 1517 Female     249  58 151.0  97.5       No 22.88       No     No      84
    ## 1518   Male     198  46 109.0  81.0      Yes 23.28       No     No      85
    ## 1519 Female     246  65 119.0  76.0       No 19.83       No     No     156
    ## 1520 Female     246  59 135.0  70.0       No 18.43       No     No     107
    ## 1521   Male     310  52 135.0  89.0       No 29.51       No     No      74
    ## 1522 Female     219  40 100.0  60.0       No 19.78       No     No      60
    ## 1523   Male     215  54 117.5  70.5      Yes 26.77       No     No      75
    ## 1524 Female     214  51 139.0  93.0       No 29.80       No     No      82
    ## 1525 Female     318  60 132.0  75.5       No 18.87       No     No      82
    ## 1526   Male     335  51 125.5  94.0      Yes 27.77       No     No      67
    ## 1527 Female     206  44 121.0  81.0       No 24.12       No     No      77
    ## 1528 Female     325  47 160.0  95.0      Yes 32.07       No     No      87
    ## 1529   Male     237  59 111.0  80.0       No 29.77       No     No      72
    ## 1530   Male     234  40 127.0  79.0       No 26.56       No     No      92
    ## 1531   Male     340  50 134.0  95.0      Yes 30.46       No     No      86
    ## 1532 Female     258  54 148.0  93.0       No 20.51       No     No      95
    ## 1534 Female     260  51 107.5  70.0       No 23.53       No     No      67
    ## 1535   Male     176  61 125.0  82.0      Yes 29.82       No     No      75
    ## 1536   Male     206  41 124.0  90.0       No 30.69       No     No      83
    ## 1537 Female     278  47 156.0  96.0       No 27.86       No     No      70
    ## 1538 Female     214  46 128.0  71.0      Yes 21.82       No     No      66
    ## 1539 Female     334  60 132.0  94.0       No 25.38       No     No      98
    ## 1540 Female     290  49 137.5  92.0       No 24.46       No     No      74
    ## 1542 Female     172  45 137.0  92.5       No 30.35       No     No      83
    ## 1543 Female     194  39 115.0  70.0       No 25.73       No     No      54
    ## 1544   Male     295  61 131.5  85.0       No 27.33       No     No      83
    ## 1545 Female     185  40 117.5  72.0       No 21.12       No     No      79
    ## 1546   Male     255  55 129.0  78.0       No 26.36       No     No      79
    ## 1547   Male     179  61 112.0  66.5       No 24.38       No     No     100
    ## 1548   Male     258  52 132.0  80.0      Yes 27.52      Yes     No     268
    ## 1549   Male     206  54 141.0  92.0       No 35.85       No     No     120
    ## 1550 Female     215  50 130.5  98.0       No 24.94       No     No     100
    ## 1551   Male     254  37 122.5  82.5      Yes 23.87       No     No      83
    ## 1552 Female     230  55 109.0  70.0       No 28.27       No     No      73
    ## 1553 Female     326  61 141.0  75.0      Yes 26.11       No     No      72
    ## 1554 Female     232  63 126.0  66.0       No 22.62       No     No      79
    ## 1555 Female     306  57 135.0  88.0       No 28.36       No     No      70
    ## 1557 Female     242  58 123.0  69.0      Yes 23.38       No     No      72
    ## 1558 Female     160  57 105.0  70.0       No 27.01       No     No      61
    ## 1559 Female     304  40 121.0  88.0      Yes 22.52       No     No      80
    ## 1560   Male     193  60 123.0  68.0       No 27.47       No     No      85
    ## 1561   Male     214  55 132.5  85.5      Yes 29.25       No     No     103
    ## 1562 Female     232  43 122.0  70.0      Yes 23.09       No     No      77
    ## 1563 Female     165  46 127.5  87.0       No 23.29       No     No      75
    ## 1564   Male     262  52 100.0  68.0      Yes 18.65       No     No      74
    ## 1565 Female     234  52 126.0  85.0       No 26.36       No     No      96
    ## 1566   Male     217  64 129.0  61.0       No 21.85       No     No      81
    ## 1567   Male     275  62 111.0  63.0       No 22.68       No     No      78
    ## 1568 Female     230  42 142.5  79.0       No 25.15       No     No      99
    ## 1569 Female     241  47 122.5  77.0      Yes 22.18       No     No      78
    ## 1571   Male     265  52 143.0  94.5       No 26.00       No     No      75
    ## 1572   Male     220  42 112.5  80.0      Yes 29.00       No     No      60
    ## 1573   Male     217  39 107.0  73.0      Yes 23.98       No     No      67
    ## 1574 Female     320  51 142.5  93.5       No 33.66       No     No      80
    ## 1575 Female     255  38 125.0  85.0       No 23.05       No     No      73
    ## 1576 Female     222  63 146.0  78.0       No 16.92       No     No      74
    ## 1577 Female     235  41 143.5  90.0       No 26.22       No     No      83
    ## 1578   Male     220  58 129.0  82.0      Yes 26.33       No     No      80
    ## 1579 Female     252  47 132.5  85.0      Yes 20.05       No     No      80
    ## 1580   Male     219  36 121.0  66.0       No 20.86       No     No      76
    ## 1581 Female     248  56 112.5  60.0       No 22.69       No     No      92
    ## 1582 Female     253  47 129.0  81.0      Yes 22.18       No     No     122
    ## 1583   Male     273  49 142.0 108.0      Yes 23.19       No     No      72
    ## 1584 Female     306  41 199.0 106.0      Yes 38.75       No     No      75
    ## 1586   Male     269  51 134.0  92.0      Yes 30.39       No     No      81
    ## 1587 Female     190  40 112.0  80.0       No 26.13       No     No      78
    ## 1589 Female     169  44 179.0 107.0       No 44.55       No     No      77
    ## 1590 Female     220  36 125.0  85.0      Yes 21.34       No     No      82
    ## 1591 Female     210  61 179.0 100.0       No 21.64       No     No      95
    ## 1592 Female     326  51 101.0  67.0      Yes 22.73       No     No      87
    ## 1593 Female     195  43 104.0  57.0      Yes 20.86       No     No      78
    ## 1594 Female     205  48 112.0  71.0       No 17.11       No     No      75
    ## 1595 Female     169  37 104.0  66.0       No 20.84       No     No      72
    ## 1596   Male     264  44 138.0  92.0      Yes 28.32       No     No      67
    ## 1597 Female     218  54 130.0  85.0       No 20.55       No     No      85
    ## 1598   Male     211  54 120.0  72.0      Yes 25.13       No     No      60
    ## 1599 Female     239  59 124.0  72.0       No 19.34       No     No      70
    ## 1600 Female     185  44 133.0  69.0      Yes 22.34       No     No      76
    ## 1602 Female     220  45 108.0  81.0       No 25.68       No     No      70
    ## 1603 Female     183  53 129.0  80.0      Yes 26.51       No     No      80
    ## 1605   Male     214  46 118.0  82.0      Yes 29.41       No     No      94
    ## 1606 Female     293  63 186.5  97.0       No 30.47       No     No      96
    ## 1607   Male     238  52 131.0  99.0       No 31.19       No     No      86
    ## 1608   Male     209  57 134.0  75.0      Yes 24.72       No     No      87
    ## 1609   Male     410  42 116.0  83.0      Yes 21.68       No     No      83
    ## 1610 Female     193  48 138.5  87.5       No 25.10       No     No      90
    ## 1611   Male     212  48 120.0  72.0      Yes 22.01       No     No      77
    ## 1612 Female     182  42 138.0  91.0       No 20.02       No     No      74
    ## 1613 Female     276  59 127.5  85.5       No 22.91       No     No      60
    ## 1614   Male     213  40 145.0 100.5       No 27.34       No     No     117
    ## 1615 Female     213  45 150.0  90.0      Yes 22.35       No    Yes      72
    ## 1616   Male     280  69 121.0  71.0       No 29.23       No     No     161
    ## 1617   Male     274  38 120.0  80.0      Yes 25.17       No     No      68
    ## 1618 Female     182  46 117.0  78.0      Yes 22.15       No     No      59
    ## 1619   Male     285  67 155.0  90.0      Yes 30.42       No     No      77
    ## 1620 Female     297  53 164.0 102.0      Yes 24.50       No     No      95
    ## 1621 Female     209  38 122.5  76.5      Yes 24.51       No     No      73
    ## 1622 Female     213  43  96.0  62.0       No 19.38       No     No      80
    ## 1623   Male     294  59 170.0 103.0       No 31.60       No    Yes      70
    ## 1624 Female     250  61 173.0  89.0       No 29.25       No     No      87
    ## 1625   Male     167  52 134.0  80.0      Yes 29.77       No     No     102
    ## 1626 Female     258  53 186.0 101.0       No 28.90       No     No      70
    ## 1627   Male     222  67 154.0 106.0       No 26.71       No     No      74
    ## 1628 Female     190  37 102.0  65.0       No 20.68       No     No      87
    ## 1629 Female     302  54 160.0  94.0       No 29.40       No     No      75
    ## 1630   Male     179  54 103.0  73.0       No 21.03       No     No      84
    ## 1631   Male     219  55 115.0  84.0       No 26.08       No     No      93
    ## 1633 Female     312  62 204.0 118.0       No 24.83       No     No      86
    ## 1634   Male     233  52 114.0  78.0      Yes 22.81       No     No      73
    ## 1636 Female     187  47 127.5  90.0       No 24.63       No     No      72
    ## 1637   Male     219  47 116.0  82.0       No 24.05       No     No      87
    ## 1638 Female     323  39 131.5  85.0      Yes 24.79       No     No      93
    ## 1642 Female     271  46 158.0  94.0      Yes 25.17       No     No      71
    ## 1643 Female     226  51 130.0  80.0       No 23.24       No     No      63
    ## 1644 Female     255  38 123.5  84.5      Yes 25.33       No     No      88
    ## 1645   Male     204  61 120.0  80.0      Yes 25.71       No     No      83
    ## 1646   Male     286  40 111.0  95.5       No 29.42       No     No      89
    ## 1648 Female     272  52 112.5  75.5       No 22.69       No     No      83
    ## 1649   Male     242  36 118.5  84.5      Yes 24.04       No     No     103
    ## 1650   Male     203  40 112.5  70.0      Yes 22.71       No     No      72
    ## 1651 Female     234  63 140.0  93.0       No 28.69       No     No      87
    ## 1652 Female     246  55 139.0  90.0      Yes 29.00       No     No     100
    ## 1653   Male     188  67 124.0  79.0      Yes 21.35       No     No      75
    ## 1654 Female     212  42 110.0  65.0       No 23.64       No     No      63
    ## 1655   Male     259  55 217.0 112.0      Yes 29.60       No     No      77
    ## 1656 Female     235  39 196.0 116.0      Yes 29.70       No     No      87
    ## 1657 Female     309  57 130.0  75.5       No 25.99       No     No      75
    ## 1658 Female     170  60 146.0  89.0       No 32.41       No     No      81
    ## 1659   Male     157  53 123.0  83.0       No 19.94       No     No      88
    ## 1660 Female     260  64 175.0 110.0       No 37.70       No     No      74
    ## 1661 Female     186  39 114.0  77.0      Yes 21.01       No     No      85
    ## 1663 Female     249  62 176.0  89.0       No 24.49       No     No      81
    ## 1664   Male     270  45 140.0  94.0      Yes 30.39       No     No      80
    ## 1665   Male     229  62 138.0  82.0      Yes 24.65       No     No      73
    ## 1666 Female     168  40 111.0  78.0       No 20.82       No     No      74
    ## 1667 Female     289  58 156.5  85.0      Yes 25.46       No     No      86
    ## 1668 Female     272  36 113.0  66.5      Yes 20.69       No     No      59
    ## 1669   Male     235  66 142.0  76.0      Yes 26.37       No     No      93
    ## 1670 Female     274  62 167.0  94.0       No 28.18       No    Yes      80
    ## 1671   Male     235  61 113.0  70.0       No 20.57       No     No      77
    ## 1672   Male     246  63 193.0 104.0       No 23.08       No     No      73
    ## 1673 Female     200  61 187.0  95.5      Yes 21.57       No     No      64
    ## 1675 Female     249  41 107.5  75.0      Yes 23.69       No     No      78
    ## 1676   Male     271  47 147.0  97.0       No 24.99       No     No      78
    ## 1677   Male     286  44 112.5  72.5       No 24.72       No     No      76
    ## 1678 Female     237  46 196.0 120.0       No 31.64       No     No      60
    ## 1679 Female     234  45 189.0  87.0      Yes 23.10       No     No      90
    ## 1680 Female     356  61 168.0  98.0      Yes 27.30       No     No     106
    ## 1681 Female     255  43 130.0  85.0       No 29.56       No     No      78
    ## 1682 Female     190  46 126.5  85.0       No 19.03       No     No      75
    ## 1683   Male     219  52 136.0  76.0       No 24.49       No     No      73
    ## 1684 Female     254  64 196.0 119.0       No 35.22       No     No      79
    ## 1685 Female     170  49 112.0  79.0       No 21.00       No     No      80
    ## 1688   Male     180  50 116.0  82.0      Yes 30.11       No     No      92
    ## 1689   Male     178  59 170.0 118.0       No 33.45       No     No      81
    ## 1690 Female     338  56 190.0  97.0       No 26.10       No     No      83
    ## 1691 Female     241  49 120.0  70.5       No 23.29       No     No      84
    ## 1693   Male     221  46 125.0  88.0      Yes 24.81       No     No      87
    ## 1695 Female     233  49 158.0 102.0      Yes 25.31       No     No      72
    ## 1696   Male     247  41 104.0  73.0       No 23.70       No     No      68
    ## 1697 Female     244  46 107.5  70.0       No 22.72       No     No      88
    ## 1698 Female     229  50 114.0  68.0      Yes 23.20       No     No      70
    ## 1699 Female     197  44 125.5  82.5       No 25.13       No     No      50
    ## 1700   Male     200  57 117.5  80.0       No 25.41       No     No      80
    ## 1701   Male     264  43 131.0  85.0      Yes 25.49       No     No      67
    ## 1702 Female     312  58 161.0 100.5       No 21.51       No     No      64
    ## 1703   Male     313  41 114.0  79.0      Yes 25.63       No     No      93
    ## 1704 Female     240  43 110.5  66.0       No 24.09       No     No      80
    ## 1705   Male     204  36 132.5  82.5      Yes 21.27       No     No      84
    ## 1706   Male     246  52 122.0  81.0      Yes 27.61       No     No      98
    ## 1707   Male     211  44 130.0  85.5       No 26.98       No     No      82
    ## 1708 Female     200  36 121.5  72.5      Yes 23.09       No     No      75
    ## 1709   Male     265  38 129.0  85.0      Yes 31.61       No     No      68
    ## 1710 Female     298  54 151.5  96.5      Yes 30.29       No     No      77
    ## 1711   Male     325  39 122.0  82.0      Yes 26.04       No     No      93
    ## 1712 Female     332  60 185.0  86.0      Yes 25.73       No     No      73
    ## 1714 Female     271  48 110.0  74.0       No 21.79       No     No      67
    ## 1715 Female     229  39 113.5  72.5      Yes 22.33       No     No      82
    ## 1716 Female     159  41 119.0  76.0       No 27.49       No     No      70
    ## 1718 Female     216  51 154.0  98.0       No 32.35       No     No     103
    ## 1720 Female     268  60 123.0  82.0       No 29.47       No     No      85
    ## 1721 Female     231  43 155.5  99.5       No 34.95      Yes     No     274
    ## 1722   Male     231  42 122.0  84.0      Yes 27.63       No     No      72
    ## 1723 Female     246  51 111.5  74.0      Yes 25.09       No     No      78
    ## 1724 Female     195  39 119.0  84.0       No 24.60       No     No      73
    ## 1725   Male     148  39 101.0  62.0      Yes 24.47       No     No      81
    ## 1727 Female     222  46 131.0  80.0       No 25.46       No     No      72
    ## 1728 Female     174  40 130.0  86.0      Yes 25.05       No     No      83
    ## 1729   Male     180  41 127.0  86.0      Yes 20.72       No     No      63
    ## 1731 Female     248  47 143.5 109.0       No 32.43       No     No      66
    ## 1732   Male     222  56 159.0  91.5      Yes 27.12       No     No      80
    ## 1733 Female     173  44 121.5  69.0      Yes 23.72       No     No      77
    ## 1734 Female     209  54 139.0  75.0      Yes 25.82       No     No      95
    ## 1735 Female     320  63 155.0  81.0       No 31.71       No     No      80
    ## 1736 Female     246  53 115.0  61.0       No 25.96       No     No      60
    ## 1737   Male     258  50 124.0  78.0      Yes 24.33       No     No      83
    ## 1738 Female     240  44 109.0  71.0      Yes 23.75       No     No      83
    ## 1739 Female     285  56 145.0 100.0      Yes 30.14       No     No      86
    ## 1740   Male     202  43 104.0  69.0      Yes 25.82       No     No      63
    ## 1741 Female     305  47 128.0  92.5       No 27.64       No     No      62
    ## 1742   Male     272  37 114.5  80.0      Yes 27.60       No     No      57
    ## 1744   Male     233  53 130.0  94.0      Yes 30.63       No     No      75
    ## 1745   Male     192  41 122.0  82.0      Yes 25.03       No     No      66
    ## 1746   Male     214  56 115.0  80.0       No 25.09      Yes     No     292
    ## 1748 Female     229  57 126.5  90.0      Yes 26.00       No     No      58
    ## 1749 Female     187  40 144.0  90.0      Yes 22.17       No     No      93
    ## 1750   Male     270  66 120.0  76.0       No 19.09       No     No      98
    ## 1751 Female     175  34 117.5  73.5      Yes 22.15       No     No      75
    ## 1752 Female     268  62 143.5  90.0       No 29.64       No     No      83
    ## 1754 Female     194  49 103.0  65.0      Yes 20.41       No     No      75
    ## 1755 Female     165  36 115.0  71.0       No 21.27       No     No      86
    ## 1756 Female     164  39 112.0  63.0       No 22.01       No     No      85
    ## 1757   Male     240  54 146.0  91.0      Yes 24.41       No     No      70
    ## 1758 Female     225  53  92.0  69.0       No 24.17       No     No      68
    ## 1759   Male     215  35 133.0  86.0      Yes 28.70       No     No      84
    ## 1761 Female     280  66 193.0 132.0       No 32.02       No     No      73
    ## 1762 Female     155  40 121.0  86.0      Yes 23.16       No     No      59
    ## 1763 Female     169  41 119.0  72.0      Yes 19.78       No     No      74
    ## 1764 Female     263  44 141.5  93.0       No 42.98       No     No      74
    ## 1765 Female     259  43 120.0  87.0      Yes 19.88       No     No      71
    ## 1766 Female     372  64 169.0  85.0       No 26.01       No     No      79
    ## 1768 Female     170  43 134.0  90.0       No 32.93       No     No      73
    ## 1770 Female     225  66 110.0  76.0      Yes 18.23       No     No      85
    ## 1771 Female     226  51 131.0  87.0       No 24.36       No     No      73
    ## 1772 Female     250  57 117.5  71.0      Yes 23.84       No     No      75
    ## 1773 Female     229  62 140.0  93.0       No 33.97       No     No     111
    ## 1774 Female     201  51 147.5  95.0      Yes 22.34       No     No      67
    ## 1775 Female     170  35 110.0  69.0       No 23.48       No     No      83
    ## 1776 Female     246  42 120.0  70.0      Yes 19.42       No     No      78
    ## 1777 Female     246  47 120.0  78.0      Yes 24.71       No     No      75
    ## 1778   Male     265  40 133.5  81.5       No 23.78       No     No      84
    ## 1779 Female     223  58 146.5  77.5       No 21.47       No     No      85
    ## 1780 Female     203  43 110.0  71.5      Yes 24.56       No     No      65
    ## 1782 Female     234  48 120.0  81.0      Yes 23.22       No     No      83
    ## 1783 Female     206  49 107.0  74.0      Yes 20.23       No     No      83
    ## 1784 Female     269  55 130.0  85.0       No 21.05       No     No      74
    ## 1785 Female     366  57 146.5  80.0       No 24.19       No     No      73
    ## 1786 Female     262  51 102.0  64.0      Yes 28.06       No     No      66
    ## 1787 Female     208  50 126.5  84.0       No 22.01       No     No      79
    ## 1788   Male     227  44 146.5  97.0      Yes 26.92       No     No      67
    ## 1789   Male     257  39 118.0  76.0      Yes 22.92       No     No      76
    ## 1790 Female     277  47 138.5  99.0       No 39.64       No    Yes      81
    ## 1792   Male     201  65 166.0  93.0       No 28.16       No    Yes      91
    ## 1793 Female     270  43 145.5  82.0      Yes 21.10       No     No      87
    ## 1794   Male     268  48 116.5  82.0      Yes 21.34       No     No      82
    ## 1797   Male     260  49 142.0  54.0      Yes 25.40       No     No      95
    ## 1799 Female     247  39 122.0  70.0      Yes 18.70       No     No      65
    ## 1800   Male     222  52 110.0  71.0       No 29.82       No     No     104
    ## 1801 Female     286  49 119.0  85.5       No 22.29       No     No      72
    ## 1802 Female     281  55 153.0  75.0       No 26.59       No     No      78
    ## 1803 Female     239  46 166.5 107.0       No 19.27       No     No      70
    ## 1804   Male     248  63 135.0  80.0      Yes 23.06       No     No     118
    ## 1806   Male     307  56 124.0  75.0      Yes 22.97       No     No      84
    ## 1807   Male     241  58 151.0 102.0       No 26.00       No     No      90
    ## 1808 Female     221  52 124.0  69.0       No 23.37       No     No      81
    ## 1809 Female     191  39 119.0  78.0      Yes 20.93       No     No      73
    ## 1810 Female     190  39 106.0  72.0      Yes 25.64       No     No      75
    ## 1811 Female     226  54 148.0  89.0       No 34.13       No     No      92
    ## 1812   Male     279  38 124.0  87.0      Yes 26.68       No     No      75
    ## 1813 Female     293  38 124.0  78.0       No 23.66       No     No      76
    ## 1814 Female     289  55 141.0  83.5      Yes 24.99       No     No      76
    ## 1815 Female     294  59 122.0  70.0       No 23.76       No     No     100
    ## 1816   Male     150  45 105.5  57.5      Yes 23.21       No     No      87
    ## 1817   Male     333  54 127.0  74.0      Yes 27.97       No     No      62
    ## 1818 Female     255  48 120.0  77.5       No 28.60       No     No      75
    ## 1819 Female     225  52 159.0  95.0       No 30.18       No     No     114
    ## 1821   Male     296  56 123.0  86.0       No 25.59       No     No      63
    ## 1822 Female     192  38 130.0  80.0      Yes 27.51       No     No      90
    ## 1825   Male     208  50 166.5 106.5       No 29.13       No     No      84
    ## 1826   Male     280  56 123.0  75.0       No 27.82       No     No     112
    ## 1827   Male     197  37 104.0  65.0       No 22.90       No     No      97
    ## 1828 Female     215  47 202.0 132.0       No 20.49       No     No      77
    ## 1829   Male     285  60 131.5  82.0      Yes 27.87       No     No     103
    ## 1830 Female     250  43 110.0  70.0      Yes 21.14       No     No      85
    ## 1831 Female     164  36 100.0  64.0      Yes 19.87       No     No      65
    ## 1832   Male     234  40 116.0  79.5      Yes 24.77       No     No      87
    ## 1833 Female     212  44 132.0  82.0       No 28.72       No     No      73
    ## 1834 Female     195  45 111.0  79.0      Yes 23.22       No     No      85
    ## 1836   Male     274  57 173.0 102.0      Yes 27.26       No     No      75
    ## 1837   Male     250  48 177.0 124.0      Yes 26.40       No     No      69
    ## 1838 Female     207  57 161.0  97.0       No 36.46       No     No      67
    ## 1839 Female     283  39 159.0 105.0       No 30.06       No    Yes      76
    ## 1840 Female     327  51 117.0  70.0      Yes 18.52       No     No      76
    ## 1841 Female     270  60 154.0  82.0       No 27.82       No     No      90
    ## 1842   Male     260  61 148.0  74.0       No 26.84       No     No      91
    ## 1843   Male     320  39 123.0  90.0       No 24.44       No     No      69
    ## 1844 Female     285  51 136.5  86.5       No 24.81       No     No      83
    ## 1845 Female     281  63 135.0  83.0       No 24.91       No    Yes      68
    ## 1846 Female     249  55 109.0  66.5      Yes 24.79       No     No      85
    ## 1848   Male     193  57 104.0  64.0      Yes 26.00       No     No      87
    ## 1849   Male     180  55 170.0 105.0       No 26.79       No     No      90
    ## 1850   Male     235  48 150.5  98.0       No 32.40       No     No      92
    ## 1851 Female     295  58 138.0  78.0       No 23.45       No     No      81
    ## 1853   Male     158  37 129.0  87.0      Yes 24.66       No     No      67
    ## 1854 Female     253  48 124.0  77.0       No 20.27       No     No      72
    ## 1856   Male     226  43 103.5  60.0      Yes 26.11       No     No      67
    ## 1857 Female     230  51 134.0  84.0       No 23.54       No     No      78
    ## 1858 Female     235  41 158.0  93.0       No 22.18       No     No      70
    ## 1859 Female     245  65 171.0  89.0       No 23.07       No     No      93
    ## 1860   Male     157  62 134.0  84.0      Yes 25.95       No     No      76
    ## 1861   Male     204  42 110.0  73.0      Yes 23.72       No     No      75
    ## 1864 Female     184  45 147.0  88.0      Yes 23.48       No     No      76
    ## 1865 Female     215  58 119.5  73.0      Yes 29.86       No     No      93
    ## 1866 Female     264  46 150.0  99.0       No 26.67       No     No     102
    ## 1867   Male     176  51 140.0  88.0      Yes 27.56       No     No      83
    ## 1868   Male     170  61 132.0  94.0      Yes 22.16       No     No      82
    ## 1869   Male     209  40 123.5  83.0       No 28.06       No     No      63
    ## 1870   Male     189  56 120.0  70.0      Yes 21.41       No     No      70
    ## 1871 Female     255  45 129.0  80.0       No 31.62       No     No     100
    ## 1872   Male     202  53 119.0  80.0       No 23.98       No     No      78
    ## 1873   Male     211  42 138.0  90.0      Yes 25.49       No     No      73
    ## 1874 Female     240  59 155.5 100.5       No 33.54       No     No     116
    ## 1875 Female     188  40 105.0  65.0      Yes 21.15       No     No      70
    ## 1876   Male     275  60 132.0  82.0       No 23.53       No     No      62
    ## 1877   Male     222  63 119.0  83.0      Yes 18.92       No     No      95
    ## 1878   Male     205  56 129.0  83.5       No 28.82       No     No      64
    ## 1881   Male     160  40 123.0  79.0      Yes 21.19       No     No      93
    ## 1882   Male     224  43  97.0  64.0      Yes 23.05       No     No      68
    ## 1883 Female     285  55 158.0  98.0       No 30.23       No    Yes      88
    ## 1884 Female     275  49 155.0  89.0       No 25.01       No     No      85
    ## 1885   Male     238  59 112.5  75.5      Yes 26.42       No     No      67
    ## 1886 Female     234  60 105.5  67.5       No 23.35       No     No      73
    ## 1887 Female     236  40 135.0  83.0      Yes 23.48       No     No      90
    ## 1888 Female     226  51 140.0  91.0       No 25.01       No     No      77
    ## 1889 Female     212  38 107.0  67.5      Yes 20.40       No     No      87
    ## 1891 Female     195  55 135.0  80.0       No 32.91       No     No      75
    ## 1892   Male     238  45 141.0  87.0      Yes 26.46       No     No      68
    ## 1893   Male     271  48 132.5  85.0      Yes 27.57       No     No      68
    ## 1894   Male     180  43 131.0  92.0      Yes 27.18       No     No      85
    ## 1895   Male     256  57 135.0  88.0      Yes 29.51       No     No      75
    ## 1896   Male     225  38 125.0  79.0      Yes 26.23       No     No      74
    ## 1897 Female     180  40 136.5  89.0       No 17.61       No     No      77
    ## 1899   Male     224  53 120.0  70.0       No 26.14       No     No      76
    ## 1900 Female     206  56 158.0  70.0       No 28.34       No     No      88
    ## 1901 Female     342  64 128.0  71.0      Yes 20.52       No     No      62
    ## 1902   Male     234  67 130.0  83.0       No 27.78       No     No      74
    ## 1903 Female     233  43 100.0  70.0      Yes 22.90       No     No      78
    ## 1904 Female     280  55 144.0  79.0      Yes 19.50       No     No      75
    ## 1906   Male     188  48 120.0  82.5      Yes 31.67       No     No      68
    ## 1907   Male     244  51 129.0  72.0      Yes 22.67       No     No      80
    ## 1909   Male     245  57 122.0  69.0      Yes 24.17       No     No      92
    ## 1910 Female     200  55 141.0  92.0       No 23.48       No     No      84
    ## 1911 Female     167  38 102.5  60.0       No 22.58       No     No      57
    ## 1912   Male     196  39 124.5  99.0       No 28.41       No     No      73
    ## 1913 Female     246  63 163.0  82.0      Yes 24.38       No     No     108
    ## 1914 Female     265  46 115.0  79.5       No 30.00       No     No      61
    ## 1916 Female     203  49 128.0  72.0       No 24.87       No     No      82
    ## 1917 Female     235  61 154.0  94.0      Yes 22.91       No     No      66
    ## 1918 Female     256  42 129.0  85.5       No 31.84       No     No      74
    ## 1919 Female     210  58 102.0  60.0       No 26.98       No     No      90
    ## 1920   Male     206  64 180.0 100.0       No 26.22       No     No      74
    ## 1921 Female     233  41 115.0  70.0       No 25.03       No     No      97
    ## 1922 Female     183  61 150.0  86.0       No 25.05       No     No      70
    ## 1923 Female     344  62 154.0 110.0       No 35.11       No     No      76
    ## 1924   Male     201  44 142.5 104.5      Yes 34.59       No     No      67
    ## 1925 Female     340  61 121.0  78.0       No 23.33       No     No      73
    ## 1926 Female     227  56 140.0  88.0       No 24.59       No     No      78
    ## 1927   Male     256  41 107.0  73.0       No 26.38       No     No      65
    ## 1928   Male     252  49 156.0  91.0      Yes 25.35       No     No     114
    ## 1929   Male     194  56 127.0  83.0       No 26.05       No     No      73
    ## 1930   Male     196  35 107.5  66.5      Yes 22.64       No     No      79
    ## 1932   Male     190  40 121.5  74.5       No 29.35       No     No      93
    ## 1933   Male     232  62 120.0  67.5       No 23.24       No     No      62
    ## 1934 Female     289  47 126.5  84.0       No 21.71       No     No      70
    ## 1935   Male     252  48 104.0  73.5      Yes 23.03       No     No      77
    ## 1936 Female     196  42 145.0  87.0      Yes 29.65       No     No      60
    ## 1937 Female     211  42 134.0  82.0       No 20.93       No     No      80
    ## 1938 Female     229  50 105.0  72.5       No 26.25       No     No      79
    ## 1939 Female     214  39 123.0  78.0       No 38.06       No     No      62
    ## 1940 Female     294  58 195.0  90.0       No 27.73      Yes    Yes     127
    ## 1941   Male     182  66 151.0  88.0       No 25.22       No     No      80
    ## 1942   Male     204  61 120.0  67.0       No 24.84       No    Yes      75
    ## 1944   Male     241  60 119.5  79.0       No 24.66       No     No      78
    ## 1945   Male     251  44 167.0 114.0       No 28.20       No     No      88
    ## 1946 Female     326  67 157.5  78.0      Yes 24.63       No     No      77
    ## 1949   Male     283  42 137.0  91.0       No 25.41       No     No      67
    ## 1950 Female     235  50 121.0  78.0      Yes 23.01       No     No      78
    ## 1951   Male     283  48 122.0  78.0       No 25.68       No     No      87
    ## 1952   Male     227  50 114.0  87.0       No 33.10       No     No     120
    ## 1953 Female     249  43 155.0  93.5       No 21.99       No     No     115
    ## 1955   Male     235  57 124.0  77.0       No 24.19       No     No      86
    ## 1956 Female     344  56 119.0  82.0       No 26.82       No     No     105
    ## 1958 Female     365  47 127.0  76.0      Yes 24.44       No     No      80
    ## 1959 Female     310  56 128.5  82.0      Yes 25.36       No     No      85
    ## 1960 Female     213  36 126.0  80.0       No 19.53       No     No      79
    ## 1961   Male     233  63 130.0  81.0       No 25.82       No     No      81
    ## 1962 Female     298  62 137.0  85.0      Yes 26.73       No     No      87
    ## 1963 Female     229  44 119.0  75.5       No 25.09       No     No      88
    ## 1964 Female     200  41 124.0  76.0      Yes 24.20       No     No      86
    ## 1965 Female     301  56 127.5  91.5       No 25.39       No     No      74
    ## 1966   Male     217  56 200.0 120.0      Yes 33.71       No     No      72
    ## 1968 Female     287  50 147.5  87.5      Yes 30.36       No     No      72
    ## 1969 Female     310  55 135.0  76.5       No 26.31       No     No      74
    ## 1970 Female     251  38 126.0  76.0       No 29.19       No     No      66
    ## 1972 Female     221  47 144.0  91.0      Yes 35.78       No     No      66
    ## 1974   Male     219  54 113.0  70.0      Yes 20.41       No     No      53
    ## 1975 Female     290  40 125.0  90.0       No 32.81       No     No      87
    ## 1976   Male     225  63 146.0  82.0      Yes 27.17       No     No      85
    ## 1977 Female     288  65 146.0  94.5       No 26.54       No     No      74
    ## 1978 Female     225  50 132.0  81.0       No 23.62       No     No     103
    ## 1979   Male     215  54 120.0  85.0      Yes 29.93       No     No      75
    ## 1980 Female     280  39 152.0 104.0      Yes 24.22       No     No      82
    ## 1981   Male     228  58 127.0  90.0       No 29.10       No     No      76
    ## 1982   Male     235  46 136.5  92.0       No 22.92       No     No      89
    ## 1983 Female     239  53 112.5  67.0      Yes 25.63       No     No      74
    ## 1984   Male     230  48 135.5  90.0      Yes 25.34      Yes     No      91
    ## 1985 Female     300  51 128.0  78.0      Yes 26.69       No     No      97
    ## 1986   Male     283  64 113.0  77.0       No 22.73       No     No      67
    ## 1987   Male     230  63 127.0  82.0      Yes 19.97       No     No      67
    ## 1988 Female     238  55 137.5  87.0       No 26.89       No     No     107
    ## 1989 Female     225  39 112.0  74.0      Yes 27.26       No     No      85
    ## 1990 Female     261  59 141.0  78.0       No 25.32       No     No      76
    ## 1991 Female     161  52 180.0 114.0      Yes 32.52       No    Yes     104
    ## 1992 Female     186  48 107.0  76.0       No 26.39       No     No      90
    ## 1993   Male     232  64 113.5  70.0      Yes 21.03       No     No      58
    ## 1994 Female     229  64 145.0  85.0       No 29.67       No     No      74
    ## 1995 Female     286  54 110.0  74.0       No 26.28       No     No      90
    ## 1996   Male     245  55 133.0  78.0       No 29.05      Yes     No     115
    ## 1997 Female     210  47 113.0  65.0       No 21.33       No     No      62
    ## 1998   Male     220  42 119.0  73.5      Yes 23.31       No     No      63
    ## 1999 Female     175  55 107.5  65.0       No 20.17       No     No      79
    ## 2000 Female     219  42 126.0  73.0      Yes 22.65       No     No      65
    ## 2001 Female     228  47 118.0  84.0       No 18.67       No     No      90
    ## 2003   Male     300  62 108.0  73.0      Yes 20.87       No     No      80
    ## 2005   Male     261  37 111.5  72.5      Yes 26.48       No     No      74
    ## 2006   Male     184  58 127.5  70.0      Yes 25.62       No     No      80
    ## 2007   Male     215  51 147.0  96.0       No 28.59       No     No     100
    ## 2008 Female     188  48 170.0 110.0      Yes 26.03       No     No     118
    ## 2009   Male     256  51 130.0  75.0      Yes 28.76       No     No      83
    ## 2010 Female     223  52 132.0  82.0       No 26.06       No     No      63
    ## 2011 Female     292  66 143.0  95.0       No 31.11       No     No      77
    ## 2012   Male     211  37 116.5  77.5      Yes 24.50       No     No      78
    ## 2014   Male     167  61 105.0  67.5      Yes 27.28       No     No      86
    ## 2015   Male     235  58 146.0  85.0      Yes 19.51       No     No      71
    ## 2017 Female     286  42 133.5  80.0      Yes 26.25       No     No      65
    ## 2018 Female     193  42 129.0  91.5      Yes 27.78       No     No      74
    ## 2019   Male     222  40 112.0  82.0       No 23.71       No     No      85
    ## 2020 Female     288  59 158.0  90.0       No 32.84       No     No      87
    ## 2021 Female     240  58 150.0  80.0       No 26.45      Yes     No     255
    ## 2023   Male     224  41 114.0  68.0      Yes 21.42       No     No     107
    ## 2024 Female     246  51 135.0  82.0       No 24.67       No     No      77
    ## 2025   Male     202  62 111.0  79.5       No 27.91       No     No     100
    ## 2026 Female     307  49 112.5  70.0       No 23.86       No     No      72
    ## 2027   Male     192  53 173.0 112.0      Yes 24.50       No     No      67
    ## 2028 Female     410  59 142.0  79.0       No 25.58       No     No      90
    ## 2029 Female     243  44 129.0  88.0      Yes 30.85       No     No      83
    ## 2030   Male     296  64 142.0  84.0       No 27.01       No     No      83
    ## 2032 Female     237  45 118.0  84.0       No 22.53       No     No      78
    ## 2033   Male     193  65 134.0  72.0      Yes 25.31       No     No      83
    ## 2034 Female     237  63 155.0  92.5       No 31.50       No     No      83
    ## 2035   Male     227  38 108.0  65.5      Yes 25.45       No     No     110
    ## 2036 Female     220  39 137.5 101.5       No 22.85       No     No      88
    ## 2037   Male     161  45 122.0  82.0      Yes 26.09       No     No      91
    ## 2038 Female     232  55 119.0  81.0       No 30.00       No     No     100
    ## 2039 Female     213  63 182.0  92.0      Yes 26.87       No    Yes      63
    ## 2040   Male     250  41 150.0  92.0       No 33.29       No     No      85
    ## 2041 Female     295  46 145.0  90.0       No 25.87       No     No      79
    ## 2042 Female     232  58 145.0  94.0       No 26.38       No     No      80
    ## 2043   Male     270  64 157.0  71.0       No 27.35       No     No      93
    ## 2045   Male     210  59 134.0  84.0       No 25.64       No     No      77
    ## 2046   Male     222  46 120.5  73.5       No 27.23       No     No      60
    ## 2047 Female     216  51 128.0  83.5       No 24.41       No     No      75
    ## 2048 Female     190  42 121.0  85.5       No 22.19       No     No      85
    ## 2049 Female     262  62 175.0  85.0       No 30.91       No     No      72
    ## 2050   Male     260  47 111.0  70.0      Yes 23.46       No     No      71
    ## 2052 Female     229  48 111.0  85.0      Yes 24.10       No     No      74
    ## 2053 Female     289  57 142.0  83.0      Yes 35.17       No     No      72
    ## 2054 Female     192  43 119.5  69.5      Yes 24.67       No     No      83
    ## 2055   Male     229  64 105.0  72.5       No 26.76       No     No      77
    ## 2056   Male     347  46 125.0  83.0      Yes 30.71       No     No     120
    ## 2057   Male     170  45 145.5  99.0      Yes 26.74       No     No      85
    ## 2058 Female     238  49  97.0  67.0       No 23.17       No     No      77
    ## 2059 Female     244  45 119.0  73.0       No 26.51       No     No      80
    ## 2060   Male     226  53 139.0  80.0      Yes 23.62       No     No      69
    ## 2061   Male     250  45 126.0  89.5      Yes 28.68       No     No      92
    ## 2062 Female     249  52 112.0  75.0      Yes 22.54       No     No      71
    ## 2063 Female     266  61 171.0  87.0       No 32.77       No     No     123
    ## 2064   Male     276  54 139.0  88.0      Yes 24.93       No     No      79
    ## 2065 Female     212  46 132.0  82.0       No 34.52       No     No      72
    ## 2067   Male     176  48 116.0  86.0      Yes 21.45       No     No      76
    ## 2068   Male     236  48 112.5  75.0      Yes 30.43       No     No      67
    ## 2069   Male     283  66 151.0  95.0      Yes 26.09       No     No      75
    ## 2072   Male     224  41 126.0  79.0      Yes 22.94       No     No      68
    ## 2073 Female     330  64 108.0  82.0       No 23.09       No     No      80
    ## 2075   Male     238  53 126.0  72.0      Yes 23.40       No     No      73
    ## 2077   Male     319  42 116.0  74.0      Yes 26.62       No     No     136
    ## 2078 Female     162  38 105.0  70.0       No 21.35       No     No      71
    ## 2084 Female     174  46 115.5  65.0      Yes 29.84       No     No      80
    ## 2085   Male     238  61 232.0 136.0       No 24.83       No     No      79
    ## 2086 Female     279  48 128.5  73.5       No 27.49       No     No      77
    ## 2087   Male     198  57 119.0  80.0       No 30.05       No     No      79
    ## 2088 Female     144  40 122.5  80.0       No 27.46       No     No     123
    ## 2089 Female     193  49 134.0  88.0       No 25.77       No     No      76
    ## 2090   Male     234  43 173.0  96.0       No 27.99       No    Yes      76
    ## 2091 Female     241  50 132.0  85.0       No 23.81       No    Yes      84
    ## 2092 Female     224  51 111.5  77.0       No 24.08       No     No      74
    ## 2093 Female     260  58  85.5  51.0       No 20.76      Yes     No     206
    ## 2095 Female     195  40 132.0  81.0       No 24.26       No     No      86
    ## 2096 Female     228  46 110.0  80.0      Yes 19.74       No     No     127
    ## 2097   Male     260  44 134.0  80.0      Yes 22.25       No     No      81
    ## 2100 Female     192  44 112.0  62.0       No 30.47       No     No      82
    ## 2101   Male     300  36 102.0  66.5      Yes 25.68       No     No     100
    ## 2103 Female     266  62 162.5 100.0       No 28.59       No     No      73
    ## 2104   Male     266  64 191.0  81.0      Yes 25.33       No    Yes      78
    ## 2106 Female     261  47 133.0  77.0      Yes 27.96       No     No     105
    ## 2108 Female     227  54 168.0  94.0       No 22.70       No     No      70
    ## 2109 Female     279  52 148.0  99.0       No 26.64       No     No      85
    ## 2110   Male     190  52 118.0  80.0      Yes 24.47       No     No      88
    ## 2112 Female     150  38 108.0  70.5       No 20.42       No     No      88
    ## 2113   Male     280  45 128.0  82.0      Yes 29.17       No     No      62
    ## 2115   Male     225  42 110.0  80.0      Yes 22.51       No     No      77
    ## 2116   Male     260  47 126.0  91.0      Yes 27.01       No     No      84
    ## 2117 Female     240  43 141.5  93.0       No 38.43       No     No      77
    ## 2119   Male     280  46 142.0  91.0       No 28.09       No     No      65
    ## 2120 Female     233  42 132.0  82.0       No 26.81       No     No      71
    ## 2122 Female     273  62 150.5  97.0       No 22.01       No     No      74
    ## 2123 Female     390  62 184.5  83.0       No 18.99      Yes     No      47
    ## 2124 Female     207  57 175.0  80.0       No 20.86       No     No      75
    ## 2126 Female     226  39  95.0  59.0      Yes 22.88       No     No      83
    ## 2127   Male     245  36 122.0  74.5      Yes 24.27       No     No      85
    ## 2128 Female     282  38 135.0  80.0      Yes 29.14       No     No      89
    ## 2129   Male     287  67 145.0  92.0       No 36.04       No     No      77
    ## 2130 Female     173  64 144.0  82.0       No 22.54       No     No      77
    ## 2131 Female     196  38 100.0  80.0       No 22.90       No     No      74
    ## 2132   Male     244  38 118.5  88.0       No 28.68       No     No      65
    ## 2133 Female     165  43 113.0  79.0      Yes 28.96       No     No      72
    ## 2134   Male     221  38 130.0  87.0       No 26.43       No     No      55
    ## 2135 Female     321  54 150.0  93.0       No 22.50       No     No     131
    ## 2136   Male     226  67 139.0  82.5      Yes 27.99       No     No      76
    ## 2137 Female     192  41 103.0  72.5       No 22.72       No     No      70
    ## 2139   Male     262  39 126.0  91.0       No 31.38       No     No      84
    ## 2140 Female     216  50 100.0  70.0       No 23.88       No     No      73
    ## 2141   Male     184  68 157.0  97.0       No 33.16      Yes     No     148
    ## 2142   Male     222  64 156.0  95.0      Yes 28.46       No     No      75
    ## 2143   Male     225  54 134.5  92.5       No 30.62       No     No      69
    ## 2144 Female     224  49 130.0  87.0       No 29.01       No     No      73
    ## 2145 Female     237  40 112.5  77.5      Yes 23.58       No     No      84
    ## 2147   Male     229  38 116.0  70.0      Yes 26.65       No     No      71
    ## 2148 Female     253  48 105.0  59.0      Yes 19.42       No     No      83
    ## 2151   Male     234  50 122.0  81.5      Yes 29.64       No     No      76
    ## 2152 Female     281  50 107.5  74.0       No 23.44       No     No      79
    ## 2153   Male     273  44 120.0  80.0       No 29.89       No     No      87
    ## 2154   Male     163  44 105.0  72.5      Yes 21.02       No     No      62
    ## 2155 Female     229  45 130.0  80.0       No 21.72       No     No      77
    ## 2156 Female     210  52 146.5  82.0       No 32.27       No     No      72
    ## 2157   Male     163  51 126.0  78.0      Yes 26.56       No     No      78
    ## 2158   Male     251  57 114.0  75.0       No 26.35       No     No      72
    ## 2159   Male     305  60 112.5  75.0      Yes 22.70       No     No      98
    ## 2160   Male     165  40 117.0  77.0      Yes 21.71       No     No      66
    ## 2162   Male     257  57 138.5  90.0      Yes 25.14       No     No     101
    ## 2163   Male     203  64 120.0  73.0       No 27.08       No     No      81
    ## 2164 Female     220  56 122.0  74.0       No 25.66       No     No      93
    ## 2166   Male     255  54 149.0  86.0      Yes 20.12       No     No      65
    ## 2167   Male     274  42 120.0  83.5      Yes 23.95       No     No      79
    ## 2170   Male     210  36 117.0  84.0      Yes 28.35       No     No      96
    ## 2171 Female     165  35 117.5  72.5       No 27.86       No     No      67
    ## 2172 Female     288  54 124.0  77.0       No 29.88       No     No      92
    ## 2173   Male     245  62 127.0  79.0      Yes 23.50       No     No      73
    ## 2174 Female     216  39 116.0  72.0      Yes 24.25       No     No      71
    ## 2176 Female     166  37 112.0  73.5      Yes 21.64       No     No      93
    ## 2177   Male     287  64 121.0  70.0       No 23.64       No     No      75
    ## 2178 Female     220  59 145.0  86.0       No 25.90       No     No      77
    ## 2179 Female     182  40  95.5  64.0      Yes 25.21       No     No      72
    ## 2181   Male     266  52 117.0  78.0       No 27.39       No     No      88
    ## 2183 Female     237  40 130.0  72.0      Yes 23.54       No     No      80
    ## 2184 Female     308  62 152.0  98.0       No 35.42       No     No      76
    ## 2185 Female     270  48 134.0  87.0       No 24.91       No     No      77
    ## 2187 Female     263  55 155.0  84.0      Yes 27.87       No     No      60
    ## 2188   Male     213  57 141.0  90.0       No 30.77       No     No      77
    ## 2189   Male     259  50 171.0 120.0       No 29.38       No     No      85
    ## 2190   Male     236  50 102.0  69.0      Yes 21.98       No     No      73
    ## 2191 Female     232  42 113.0  69.0       No 21.29       No     No      70
    ## 2192   Male     202  65 186.0 109.0       No 30.66       No    Yes      75
    ## 2194 Female     185  42 123.0  74.0      Yes 24.41       No     No      92
    ## 2195   Male     228  66 188.0 128.0      Yes 29.58       No     No      67
    ## 2196   Male     254  62 161.0 104.0       No 27.43       No    Yes      82
    ## 2197 Female     177  41 107.0  76.0       No 22.37       No     No      65
    ## 2198   Male     265  57 124.0  81.0      Yes 28.18       No     No     100
    ## 2200   Male     206  43 162.5  93.5      Yes 29.01       No     No      77
    ## 2201 Female     209  39 104.0  72.0       No 23.96      Yes     No     103
    ## 2202 Female     243  50 157.0  98.0       No 23.82       No     No      78
    ## 2203 Female     305  50 114.0  80.0       No 24.33       No     No      80
    ## 2204 Female     200  33 119.0  74.0      Yes 23.80       No     No      74
    ## 2205 Female     264  36 123.0  94.5       No 28.59       No     No      70
    ## 2208 Female     309  53 130.0  86.0      Yes 22.37       No     No      80
    ## 2209   Male     250  67 120.0  83.0       No 22.36       No     No      78
    ## 2214 Female     275  55 144.5  88.5      Yes 27.05       No     No      79
    ## 2215   Male     260  56 126.5  76.0       No 25.14       No     No      70
    ## 2216   Male     217  39 119.0  85.0      Yes 23.32       No     No      96
    ## 2217   Male     217  62 107.5  75.0      Yes 26.21       No     No      66
    ## 2218   Male     258  57 146.0  76.0       No 24.94       No     No      87
    ## 2219 Female     180  39 113.0  73.0      Yes 17.65       No     No      73
    ## 2220   Male     177  56 127.0  79.0      Yes 20.12       No     No      88
    ## 2222   Male     258  54 146.0  98.5       No 26.05       No     No      68
    ## 2224 Female     195  39 137.0  93.0       No 26.39       No     No      75
    ## 2225 Female     160  37 104.0  61.0      Yes 20.22       No     No     100
    ## 2226 Female     266  60 135.0  87.0       No 26.72       No     No      82
    ## 2227 Female     254  62 167.5 102.5      Yes 27.15       No     No      83
    ## 2228   Male     200  58 136.0  88.0      Yes 26.25       No     No      73
    ## 2229 Female     195  51 154.0  96.0      Yes 28.38       No     No      75
    ## 2230   Male     234  48 141.0  98.0       No 21.06       No     No      82
    ## 2231 Female     225  58 146.0  77.0      Yes 24.60       No     No      53
    ## 2232   Male     225  54 133.0  83.0      Yes 22.18       No     No      65
    ## 2233 Female     276  63 134.0  85.0       No 23.64       No     No      86
    ## 2234   Male     179  37 125.0  82.0      Yes 19.53       No     No      70
    ## 2235 Female     216  59 205.0  92.5      Yes 25.86       No     No      84
    ## 2236   Male     298  55 169.5 104.5       No 27.51       No     No      78
    ## 2238 Female     143  40 125.5  80.0      Yes 21.99       No     No      95
    ## 2239 Female     290  56 185.0 107.5       No 26.45       No     No      84
    ## 2240 Female     195  41 120.5  76.0      Yes 22.91       No     No      70
    ## 2242 Female     164  37  96.5  67.0      Yes 24.99       No     No      67
    ## 2243 Female     259  50 152.0  97.0       No 33.68       No     No      76
    ## 2244 Female     206  41 130.0  88.0       No 22.25       No     No      79
    ## 2245   Male     298  41 132.0  85.5       No 31.06       No     No      90
    ## 2246 Female     283  49 127.0  86.0      Yes 23.68       No     No      78
    ## 2247   Male     234  51 129.0  94.0      Yes 23.90       No     No      96
    ## 2248   Male     202  62 149.5  85.0       No 25.42       No     No      80
    ## 2249   Male     258  42 110.0  69.0       No 26.25       No     No      73
    ## 2250 Female     293  59 124.0  74.0       No 25.56       No     No      77
    ## 2251 Female     245  52 148.0  92.0       No 22.70       No     No      60
    ## 2252 Female     154  37 106.0  59.5      Yes 22.71       No     No      50
    ## 2253   Male     310  36 126.0  87.0      Yes 28.09       No     No      78
    ## 2254   Male     186  36 121.0  79.5       No 27.08       No     No      65
    ## 2255 Female     170  39 110.5  69.0      Yes 22.19       No     No     103
    ## 2256 Female     271  55 146.5  80.0       No 20.69       No     No      89
    ## 2257 Female     200  40 116.0  69.0      Yes 23.90       No     No      78
    ## 2258   Male     305  38 130.0  95.0       No 23.10       No     No      78
    ## 2259 Female     265  44 110.0  78.0      Yes 20.88       No     No      68
    ## 2260   Male     260  54 116.0  77.0       No 28.56       No     No      61
    ## 2261   Male     162  40 129.0  76.5      Yes 24.12       No     No      73
    ## 2262 Female     280  66 188.0  92.0       No 24.49       No     No     101
    ## 2263 Female     212  42 115.0  72.0      Yes 23.72       No     No     100
    ## 2266 Female     212  66 220.0  96.0       No 44.71       No     No      95
    ## 2267 Female     274  54 116.0  79.0      Yes 24.77       No     No      65
    ## 2268 Female     274  57 116.0  83.0       No 25.77       No     No      83
    ## 2269   Male     328  38 124.0  82.0       No 29.08       No     No      69
    ## 2270 Female     165  42 146.0  90.0       No 28.78       No     No      74
    ## 2271 Female     240  42 127.0  90.0      Yes 38.54       No     No      69
    ## 2272 Female     255  45 111.0  72.0      Yes 17.32       No     No      65
    ## 2273   Male     210  48 134.0  83.0       No 25.08       No     No     107
    ## 2274   Male     249  44 182.0 111.0       No 39.04       No     No      67
    ## 2275 Female     205  62 118.0  79.5       No 30.21       No     No      75
    ## 2276 Female     203  41 124.0  86.0       No 28.25       No     No      80
    ## 2279   Male     247  49 150.0  88.0      Yes 27.92       No     No      74
    ## 2280 Female     205  54 107.5  67.5       No 20.26       No     No      88
    ## 2282 Female     229  43 124.0  79.0       No 30.28       No     No      77
    ## 2283 Female     273  61 210.0 120.0       No 25.11       No     No      83
    ## 2284   Male     206  65 106.0  64.0      Yes 24.74       No     No      81
    ## 2285   Male     229  63 109.0  75.0       No 38.42      Yes     No     120
    ## 2287   Male     300  38 120.0  74.0      Yes 28.74       No     No      78
    ## 2289 Female     279  53 132.0  81.0       No 26.18       No     No      77
    ## 2290 Female     250  45 130.0  80.0      Yes 20.24       No     No      86
    ## 2291 Female     293  51 151.0  92.0      Yes 30.67       No     No      77
    ## 2292   Male     238  52 184.0 102.0      Yes 28.88       No    Yes      94
    ## 2293   Male     297  40 105.0  73.0      Yes 27.15       No     No      88
    ## 2294   Male     232  42 130.0  91.0      Yes 25.77       No     No      70
    ## 2295   Male     165  37 134.5  91.0      Yes 27.97       No     No      80
    ## 2296 Female     232  52 115.0  80.0      Yes 28.79       No     No      68
    ## 2297   Male     209  44 127.0  86.0      Yes 26.41       No     No      88
    ## 2298   Male     339  44  97.0  62.0      Yes 22.19       No     No      85
    ## 2299   Male     272  42 128.0  83.0      Yes 33.26       No     No      63
    ## 2300   Male     288  45 124.0  81.0      Yes 27.94       No     No     118
    ## 2301   Male     190  39 114.0  70.0       No 27.80       No     No      80
    ## 2302 Female     293  63 193.0  63.0       No 30.00       No     No      76
    ## 2303   Male     190  42 121.5  79.0      Yes 24.20       No     No      77
    ## 2304   Male     260  55 120.0  80.0       No 28.89       No     No      68
    ## 2305 Female     230  46 154.0  98.0       No 28.23       No     No      90
    ## 2306   Male     258  45 128.0  79.0       No 32.03       No     No      75
    ## 2307 Female     209  37 115.0  69.0       No 24.66       No     No      77
    ## 2308 Female     232  50 120.0  75.0       No 23.88       No     No      74
    ## 2309   Male     260  52 158.0  96.0      Yes 21.17       No     No      82
    ## 2310   Male     260  53 120.0  80.0      Yes 29.35       No     No      73
    ## 2311   Male     405  46 181.5 102.5      Yes 26.33       No     No      97
    ## 2312 Female     197  35 109.0  73.0       No 26.38       No     No      95
    ## 2313 Female     359  42 115.0  71.0      Yes 24.46       No     No      68
    ## 2315 Female     234  52 130.0  75.0       No 28.35       No     No      72
    ## 2316   Male     280  63 185.0 105.0       No 28.41      Yes     No     177
    ## 2317   Male     230  41 150.0 101.0       No 28.54       No     No      65
    ## 2318 Female     174  48 154.0  84.0      Yes 31.76       No     No      90
    ## 2319   Male     295  62 139.0  80.0      Yes 21.28       No     No      97
    ## 2320 Female     243  53 188.5 106.5       No 29.82       No     No      70
    ## 2321   Male     260  60 127.5  72.5      Yes 25.06       No     No      75
    ## 2323   Male     195  39 106.0  80.0      Yes 23.57      Yes     No     132
    ## 2324   Male     226  43 132.5  85.0      Yes 26.64       No     No      58
    ## 2326   Male     182  42 120.0  83.0      Yes 27.26       No     No      87
    ## 2327   Male     288  60 122.5  68.5      Yes 22.13       No     No      88
    ## 2328 Female     322  65 165.0  95.0      Yes 22.84       No     No      81
    ## 2329 Female     230  45 116.0  75.0      Yes 21.35       No     No      77
    ## 2330 Female     309  38 113.0  68.0      Yes 21.35       No     No      75
    ## 2331   Male     240  38 122.5  80.0       No 23.97       No     No      43
    ## 2332 Female     205  59 142.5  74.5       No 25.63       No     No      83
    ## 2333 Female     232  55 170.0  92.0       No 26.09       No     No      74
    ## 2334   Male     193  62 132.5  80.0       No 27.20       No     No      78
    ## 2335 Female     167  37 118.0  76.0      Yes 19.61       No     No      67
    ## 2336   Male     238  38 119.0  83.5      Yes 24.18       No     No      71
    ## 2337   Male     285  52 110.0  79.0       No 23.41       No     No      65
    ## 2338 Female     239  45 112.0  70.0      Yes 23.48       No     No      95
    ## 2339 Female     235  41 107.5  68.0      Yes 21.10       No     No     113
    ## 2340   Male     219  52 125.0  82.0       No 24.06      Yes     No     173
    ## 2341   Male     207  59 111.5  67.0      Yes 20.12       No     No      84
    ## 2342 Female     260  53 139.0  80.0       No 20.31       No     No      76
    ## 2343   Male     284  63 182.5  97.5       No 23.98       No     No      87
    ## 2344 Female     246  48 165.0  84.0       No 27.60       No     No      73
    ## 2345   Male     288  54 145.0  92.5      Yes 26.20       No     No      98
    ## 2346   Male     255  39 108.0  75.0      Yes 23.90       No     No      70
    ## 2347   Male     270  46 138.0  97.0       No 33.79       No     No      65
    ## 2348   Male     239  45 116.0  68.0      Yes 26.00       No     No      79
    ## 2349 Female     202  50 138.0  72.0      Yes 25.03       No     No      98
    ## 2350   Male     233  46 108.0  74.0      Yes 23.97       No     No      82
    ## 2351 Female     200  36 108.0  62.0       No 20.79       No     No      69
    ## 2352 Female     243  53 164.0 111.0       No 39.53       No     No      70
    ## 2353   Male     250  50 148.0 108.0      Yes 24.00       No     No      86
    ## 2354 Female     235  63 131.5  91.0       No 24.69       No     No      68
    ## 2355 Female     178  38  96.0  67.0      Yes 20.40       No     No      82
    ## 2356 Female     207  56 120.0  60.0       No 22.89       No     No      71
    ## 2358 Female     220  38 114.0  73.5      Yes 27.06       No     No      67
    ## 2359   Male     205  48 148.0 103.0       No 28.31       No     No      60
    ## 2360   Male     197  56 140.0  86.0      Yes 25.16       No     No      71
    ## 2362   Male     213  41 120.0  78.0      Yes 28.78       No     No      70
    ## 2363   Male     232  48 130.5  90.0       No 29.07       No     No     118
    ## 2364   Male     238  57 128.5  87.5       No 25.38       No     No      89
    ## 2365   Male     258  44 116.0  67.0      Yes 28.30       No     No      79
    ## 2366   Male     208  40 120.0  80.0      Yes 25.98       No     No      66
    ## 2367 Female     268  45 100.0  70.0      Yes 23.45       No     No      87
    ## 2368   Male     314  57 109.5  72.0       No 25.62       No     No      71
    ## 2371   Male     260  42 132.0  86.0      Yes 29.76       No     No      66
    ## 2372   Male     148  38 115.0  75.0       No 26.49       No     No      74
    ## 2373   Male     199  41 139.0  80.0      Yes 25.51       No     No      62
    ## 2375 Female     275  63 148.0  75.0       No 28.87       No     No      83
    ## 2378   Male     346  52 133.0  96.0      Yes 25.95       No     No     126
    ## 2379 Female     237  48 124.5  66.5       No 33.29       No     No      91
    ## 2381   Male     193  44 134.0  88.0      Yes 23.77       No     No      77
    ## 2382 Female     350  49 174.0  90.0      Yes 18.44       No    Yes      78
    ## 2383   Male     166  42 110.0  70.0      Yes 19.97       No     No      69
    ## 2384 Female     219  42 120.0  70.0      Yes 24.10       No     No      73
    ## 2385 Female     281  57 192.0 105.0      Yes 27.04       No    Yes      75
    ## 2386   Male     272  57 109.0  77.0       No 25.57       No     No      89
    ## 2387 Female     243  54 138.0  79.0       No 21.93       No     No      95
    ## 2388   Male     208  39 146.0  92.0       No 25.91       No     No      74
    ## 2389 Female     240  50 176.5 115.0      Yes 27.71       No     No      83
    ## 2390 Female     240  49 110.0  71.0      Yes 22.02       No     No      84
    ## 2391   Male     206  49 101.0  69.0      Yes 28.40       No     No      82
    ## 2392   Male     251  48 148.0  91.0       No 25.79       No     No      83
    ## 2393 Female     270  39 110.0  78.0      Yes 22.00       No     No      68
    ## 2394 Female     292  54 137.0  84.0       No 24.96       No     No      94
    ## 2395 Female     178  68 154.5  83.0      Yes 18.88       No     No      86
    ## 2396   Male     239  48 127.5  62.5       No 24.39       No     No      79
    ## 2397 Female     235  51 110.5  69.0       No 21.80       No     No      86
    ## 2398   Male     212  54 153.0 100.0      Yes 27.03       No     No      82
    ## 2399 Female     200  58 177.0  97.0       No 28.40       No     No      73
    ## 2400 Female     238  45 127.5  83.5       No 27.73       No     No      78
    ## 2401   Male     227  60 119.0  76.0       No 24.80       No     No      92
    ## 2402   Male     256  56 140.0  81.0      Yes 27.10       No     No      92
    ## 2404 Female     254  40 132.0  92.0       No 23.45       No     No      62
    ## 2405 Female     213  64 133.0  77.5       No 35.62       No     No      74
    ## 2406 Female     190  47 160.0  96.0       No 28.02       No     No      64
    ## 2407   Male     195  54 113.5  72.0      Yes 21.78       No     No      67
    ## 2408   Male     230  60 149.0  95.0      Yes 26.68       No     No      92
    ## 2409   Male     336  48 183.0 108.0       No 28.11       No    Yes      96
    ## 2410 Female     204  41 129.0  86.0       No 20.72       No     No      70
    ## 2411 Female     180  36 116.0  85.5       No 26.32       No     No      81
    ## 2412   Male     230  37 123.0  92.0       No 28.61       No     No      80
    ## 2413 Female     175  39 114.0  68.0       No 24.20       No     No      93
    ## 2414   Male     235  63 199.0 114.0       No 29.76       No     No      84
    ## 2415   Male     199  66 159.0  92.0      Yes 26.35       No     No      74
    ## 2416 Female     202  48 171.0  97.0       No 32.67       No     No      78
    ## 2417 Female     274  37 197.5 125.0      Yes 43.48       No     No      94
    ## 2418 Female     246  64 190.0 100.0       No 32.68       No     No     102
    ## 2419 Female     275  46 116.0  88.5       No 28.04       No     No      70
    ## 2420   Male     203  38 116.0  81.0       No 30.19       No     No      80
    ## 2421 Female     245  43 112.5  80.0      Yes 23.43       No     No      77
    ## 2422 Female     205  44 120.0  83.5      Yes 24.30       No     No      77
    ## 2423   Male     213  44 116.0  77.5      Yes 26.09       No     No      83
    ## 2426 Female     233  46 174.0 100.0       No 26.72       No     No      73
    ## 2427 Female     200  46 133.0  84.0       No 26.15       No     No      83
    ## 2428   Male     195  60 107.5  72.0      Yes 26.53       No     No      80
    ## 2429   Male     312  61 110.0  66.0      Yes 26.28       No     No      96
    ## 2430 Female     192  51 133.0  87.0      Yes 20.72       No     No      77
    ## 2431   Male     176  63 116.0  83.0      Yes 27.80       No     No      75
    ## 2432   Male     157  61 195.0 108.0      Yes 25.08       No     No      78
    ## 2433   Male     183  36 140.0  85.0       No 24.39       No     No      84
    ## 2434   Male     249  43 145.0  85.0      Yes 28.77       No     No     100
    ## 2435   Male     233  60 135.0  75.0      Yes 22.17       No     No      60
    ## 2436   Male     259  46 101.0  71.0      Yes 20.10       No     No      73
    ## 2437 Female     252  48 143.0  81.0      Yes 24.00       No     No     101
    ## 2438   Male     273  66 145.0  88.0      Yes 25.41       No     No      74
    ## 2439   Male     193  40 122.0  78.0      Yes 28.40       No     No      93
    ## 2440 Female     199  42 111.0  71.0       No 27.23       No     No     107
    ## 2441 Female     229  44 146.0  78.0      Yes 25.25       No     No      75
    ## 2442   Male     181  46 122.0  87.5      Yes 29.15       No     No      83
    ## 2443 Female     168  50 120.0  80.0      Yes 25.26       No     No      60
    ## 2444 Female     267  40 150.0  93.0       No 31.77       No     No      74
    ## 2445   Male     261  50 180.0 100.0      Yes 25.90       No     No      66
    ## 2446 Female     235  49 109.0  70.0      Yes 28.66       No     No      73
    ## 2447 Female     281  42 115.5  79.0      Yes 22.90       No     No      71
    ## 2448 Female     201  39 121.0  78.5       No 24.26       No     No      90
    ## 2449   Male     215  49 132.0  85.0       No 33.14       No     No      75
    ## 2453 Female     218  41 102.0  67.0      Yes 19.23       No     No      78
    ## 2454   Male     266  37 110.0  72.5       No 26.09       No     No      73
    ## 2455   Male     217  43 115.0  80.0       No 28.82       No     No      70
    ## 2456 Female     241  64 136.5  85.0      Yes 26.42       No     No      77
    ## 2458 Female     226  40 118.0  72.0       No 24.75       No     No      79
    ## 2459   Male     216  48 144.0  79.0      Yes 25.09       No     No      85
    ## 2460 Female     380  63 175.0  78.0       No 20.15       No     No      95
    ## 2461 Female     233  51 120.0  81.0       No 28.25       No     No      75
    ## 2462 Female     270  64 146.0 106.0       No 27.96       No     No      77
    ## 2463 Female     281  58 150.0 101.0       No 36.91       No     No      97
    ## 2464 Female     299  48 132.0  81.0      Yes 24.35       No     No      70
    ## 2465   Male     243  37 125.0  75.0      Yes 29.19       No     No      78
    ## 2467   Male     300  44 139.0  92.5       No 25.66       No     No      59
    ## 2468 Female     293  55 165.0 106.0       No 24.71       No     No      79
    ## 2469 Female     270  42 126.0  86.0       No 23.98       No     No      77
    ## 2471   Male     262  43 105.0  72.5      Yes 28.36       No     No      70
    ## 2472 Female     212  44 128.0  81.5      Yes 27.51       No     No      87
    ## 2473   Male     232  37 129.0  74.0      Yes 24.46       No     No      88
    ## 2474 Female     254  50 132.0  96.0       No 29.35       No     No      67
    ## 2475   Male     200  58 128.0  83.0       No 29.63       No     No      80
    ## 2476   Male     355  42 113.0  81.0      Yes 26.17       No     No      71
    ## 2477 Female     250  55 131.0  96.0       No 31.34       No     No      80
    ## 2478   Male     210  48 134.5  84.0      Yes 27.97       No     No      87
    ## 2479 Female     242  32 111.0  70.0      Yes 29.84       No     No      88
    ## 2480 Female     240  39 118.0  84.0       No 25.57       No     No      75
    ## 2482 Female     390  54 150.0  94.0      Yes 27.34       No     No      71
    ## 2483   Male     223  53 142.5  87.0       No 25.56      Yes     No      75
    ## 2484 Female     190  50 148.0  92.5       No 22.99       No     No      90
    ## 2485 Female     212  60 146.0  80.0       No 28.53       No     No      72
    ## 2486 Female     315  54 176.0  87.0       No 29.23       No     No      72
    ## 2487 Female     214  58 158.5  94.5      Yes 29.14       No     No      97
    ## 2489 Female     230  46 115.0  75.0       No 22.90       No     No      79
    ## 2490 Female     189  50 144.0  88.0       No 39.08       No     No      87
    ## 2491   Male     218  45 133.0  87.0       No 31.90       No     No     115
    ## 2493   Male     264  60 137.0  80.0       No 29.99       No     No      58
    ## 2494 Female     279  52 135.0  86.0       No 27.02       No     No      72
    ## 2495   Male     287  58 145.0  72.5       No 26.27      Yes     No     206
    ## 2496 Female     186  38 166.0  96.0      Yes 33.14       No     No      75
    ## 2497   Male     235  59 136.0  96.0      Yes 28.61       No     No      85
    ## 2498 Female     156  41 121.0  88.0       No 23.57       No     No      86
    ## 2499   Male     275  46 137.0  88.0      Yes 29.28       No     No      88
    ## 2500 Female     310  53 146.0  91.0       No 29.30       No     No      72
    ## 2501   Male     259  50 132.5  92.0       No 30.41       No     No      75
    ## 2502 Female     240  57 142.0  85.0      Yes 22.55       No     No      77
    ## 2503   Male     332  40 114.0  76.0      Yes 26.51       No     No      90
    ## 2504 Female     228  41 108.0  62.0       No 23.92       No     No      99
    ## 2505 Female     280  41 129.0  89.0      Yes 39.69       No     No      65
    ## 2506   Male     197  63 112.0  80.0      Yes 25.42       No     No      76
    ## 2507 Female     254  56 113.0  61.0       No 25.50       No     No      87
    ## 2508 Female     205  49 137.0  95.0       No 25.29       No     No      74
    ## 2509   Male     265  47 137.5  88.5       No 23.75       No     No      90
    ## 2510   Male     227  64 102.0  61.0       No 23.29       No     No      81
    ## 2511 Female     198  57 208.0 104.0      Yes 24.91       No     No      82
    ## 2512 Female     342  59 137.0  83.5      Yes 25.18      Yes     No     140
    ## 2513   Male     302  40 122.0  95.0      Yes 23.81       No     No      82
    ## 2514 Female     152  45 118.5  76.0      Yes 25.38       No     No      68
    ## 2515 Female     193  46 110.0  68.0       No 22.19       No     No      65
    ## 2516 Female     196  42 136.0  88.0       No 34.55       No     No      78
    ## 2518 Female     326  41 124.0  78.0      Yes 25.27       No     No      67
    ## 2519 Female     135  35 105.0  69.0       No 22.88       No     No      76
    ## 2520 Female     156  37 120.0  87.0      Yes 21.80       No     No      89
    ## 2521   Male     230  55 130.0  85.0      Yes 30.16       No     No     108
    ## 2522   Male     303  56 136.5  97.0      Yes 26.64       No     No     106
    ## 2523   Male     220  51 125.0  82.0      Yes 24.10       No     No      73
    ## 2524 Female     252  49 123.0  69.0      Yes 21.45       No     No      89
    ## 2525 Female     233  62 130.0  87.0       No 21.34      Yes     No     386
    ## 2526 Female     323  49 123.5  78.0      Yes 22.86       No     No      63
    ## 2527 Female     288  44 150.0  89.0      Yes 21.11       No     No      97
    ## 2529   Male     190  39 118.0  70.0      Yes 23.57       No     No      69
    ## 2530 Female     255  51 102.5  64.5      Yes 24.14       No     No      71
    ## 2531   Male     218  48 113.0  79.0      Yes 27.33       No     No      73
    ## 2532   Male     179  43 129.0  81.0      Yes 19.05       No     No      77
    ## 2533 Female     212  65  94.0  62.0       No 25.83       No     No      88
    ## 2534 Female     216  65 163.0 102.0       No 30.12       No     No      73
    ## 2535   Male     205  56 210.0 130.0      Yes 25.49       No     No     127
    ## 2537   Male     188  45 132.0  91.0      Yes 28.04       No     No      77
    ## 2538 Female     260  61 112.5  75.0       No 21.97       No     No      74
    ## 2539 Female     180  53 184.5 110.5       No 27.49       No     No      74
    ## 2540 Female     221  65 141.0  82.5       No 29.48       No     No      93
    ## 2541 Female     159  34  92.5  70.0      Yes 22.15       No     No      68
    ## 2542   Male     216  46 126.0  85.0      Yes 30.10       No     No      70
    ## 2543   Male     193  49 172.0 105.0      Yes 19.70       No     No      77
    ## 2544 Female     229  63 120.0  82.0       No 25.58       No     No      73
    ## 2545   Male     268  55 128.0  87.0      Yes 28.57       No     No      87
    ## 2546   Male     205  42 110.0  73.0      Yes 22.40       No     No      66
    ## 2547 Female     214  44 102.0  68.0      Yes 32.82       No     No      80
    ## 2548   Male     193  49 115.0  79.0      Yes 21.86       No     No      88
    ## 2551   Male     189  49 132.0  76.0      Yes 20.26       No     No      60
    ## 2552 Female     248  52 128.0  74.0      Yes 21.84       No     No      79
    ## 2553 Female     185  41 125.0  88.0       No 24.88       No     No      75
    ## 2554   Male     155  36 127.5  77.5      Yes 30.20       No     No      61
    ## 2555 Female     161  50 145.0  89.0      Yes 20.30       No     No      81
    ## 2556   Male     266  59 165.0  99.0      Yes 29.07       No     No      82
    ## 2557   Male     246  66 144.0  98.0       No 35.05      Yes     No      95
    ## 2558   Male     275  44 127.0  82.5      Yes 27.43       No     No      76
    ## 2559 Female     224  49 126.0  82.5       No 27.49       No     No      67
    ## 2560   Male     223  51 117.5  77.5       No 27.78       No     No      97
    ## 2561 Female     325  56 160.0  97.5      Yes 23.40       No     No      86
    ## 2562   Male     278  38 136.5  87.0       No 29.84       No     No      75
    ## 2563   Male     124  66 138.0  91.0       No 32.33       No     No      96
    ## 2564 Female     253  54 127.0  68.5       No 22.54       No     No      76
    ## 2565 Female     250  48 103.0  76.5      Yes 23.25       No     No      66
    ## 2566 Female     213  63 172.0  95.0      Yes 27.68       No    Yes      67
    ## 2567   Male     264  63 129.5  82.5      Yes 24.91       No     No      83
    ## 2568 Female     214  50 176.0 113.0      Yes 22.17       No     No      71
    ## 2569   Male     205  51 125.0  94.0       No 29.11       No     No      90
    ## 2570   Male     243  48 116.0  82.0       No 26.09       No     No      85
    ## 2571 Female     295  60 162.0  85.0       No 32.76       No     No      71
    ## 2572   Male     187  48 102.0  69.0      Yes 24.24       No     No      72
    ## 2573   Male     189  58 136.0  86.0       No 23.97       No     No      60
    ## 2574 Female     242  44 124.0  72.5       No 23.07       No     No      83
    ## 2575 Female     236  39 127.0  78.0      Yes 17.51       No     No      76
    ## 2576 Female     185  37 115.0  76.0      Yes 23.55       No     No      80
    ## 2577   Male     270  56 109.0  75.0      Yes 24.58       No     No      64
    ## 2578   Male     285  59 128.0  91.0       No 28.23       No     No      80
    ## 2579 Female     185  53 202.5  85.0       No 24.65       No     No     104
    ## 2580   Male     299  59 120.5  78.0       No 25.45       No     No     100
    ## 2581   Male     266  65 140.0 100.0       No 29.36       No     No      77
    ## 2583   Male     328  54 136.5  85.0      Yes 33.71       No     No      85
    ## 2584 Female     290  58 145.0  88.0       No 36.65       No    Yes      86
    ## 2585 Female     261  40 112.0  67.0      Yes 21.83       No     No      61
    ## 2586   Male     237  37 102.0  72.0      Yes 19.68       No     No      83
    ## 2587   Male     189  64 156.0  69.0       No 21.68       No     No     100
    ## 2588   Male     208  62 140.0  97.5       No 27.27       No    Yes     121
    ## 2590   Male     276  46 155.0  85.0      Yes 28.46       No    Yes      64
    ## 2591   Male     259  48 135.0  90.0      Yes 20.72       No     No      81
    ## 2592   Male     212  37 120.0  72.0      Yes 23.51       No     No      80
    ## 2593   Male     202  52 136.0  94.0       No 29.93       No     No      67
    ## 2594   Male     227  51 139.0  74.0      Yes 29.29       No     No      67
    ## 2595 Female     261  66 154.0  97.0      Yes 32.60       No    Yes      81
    ## 2597 Female     239  54 142.5  83.5       No 30.47       No     No      84
    ## 2598   Male     165  55 166.0 101.0       No 24.80       No     No      90
    ## 2599   Male     198  57 128.0  85.0      Yes 28.18       No     No      99
    ## 2600   Male     232  39 122.5  78.5      Yes 26.11       No     No      73
    ## 2601   Male     252  59 146.0  92.0       No 27.88       No     No      80
    ## 2602 Female     260  55 136.5  87.5       No 25.41       No     No      60
    ## 2603   Male     230  42 112.0  66.0      Yes 25.31       No     No      69
    ## 2604 Female     260  45  98.0  74.0       No 19.16       No     No      76
    ## 2605 Female     271  48 130.0  84.0      Yes 21.97       No     No      85
    ## 2606 Female     220  50 122.0  80.0      Yes 24.22       No     No      72
    ## 2607 Female     200  37 112.5  68.0      Yes 25.87       No     No      67
    ## 2608 Female     273  57 131.0  84.0       No 22.99       No     No      76
    ## 2609   Male     235  57 150.0  95.5      Yes 27.56       No     No      73
    ## 2610 Female     201  42 141.0  84.5      Yes 26.58       No     No      97
    ## 2611 Female     371  67 166.0  85.0      Yes 25.35       No     No      86
    ## 2612   Male     201  44 120.0  81.0       No 26.49       No     No      78
    ## 2613 Female     309  59 141.0  77.5      Yes 25.97       No     No      79
    ## 2614   Male     190  47 120.0  81.0       No 24.54       No     No      73
    ## 2615 Female     312  59 175.0  82.0      Yes 39.82       No    Yes      85
    ## 2616   Male     300  46 146.0  98.5       No 30.41       No     No      79
    ## 2617   Male     275  61 164.0 104.0       No 23.37       No     No      83
    ## 2618 Female     227  39 138.0  89.0       No 26.74       No     No      60
    ## 2619   Male     215  45 104.0  72.0       No 30.34       No     No      79
    ## 2620   Male     235  53 130.0  80.0       No 28.15       No     No      78
    ## 2621 Female     270  67 137.5  72.5      Yes 35.01       No     No      73
    ## 2622   Male     300  57 121.0  74.0       No 28.09      Yes     No     155
    ## 2623 Female     232  46 115.0  70.0      Yes 25.18       No     No      59
    ## 2624   Male     277  45 140.0  84.0      Yes 28.74       No     No      74
    ## 2625   Male     270  60 145.0  81.0       No 29.37       No     No      73
    ## 2626   Male     219  53 141.0 105.0       No 26.86       No     No      60
    ## 2627   Male     309  43 124.0  85.0      Yes 26.91      Yes     No     215
    ## 2628 Female     199  45 124.0  78.0      Yes 21.94       No     No      78
    ## 2629 Female     258  68 158.0  94.0      Yes 31.64       No     No      84
    ## 2630 Female     199  57 117.0  83.0       No 24.76       No     No      82
    ## 2632 Female     199  54 159.0 102.0       No 22.91       No    Yes      93
    ## 2633 Female     223  44  96.0  59.0      Yes 23.82       No     No      87
    ## 2634 Female     212  50 169.0 117.0       No 27.08       No     No      68
    ## 2635 Female     206  39 102.5  65.0      Yes 19.80       No     No      85
    ## 2636 Female     172  43 149.0  82.0      Yes 22.35       No    Yes      64
    ## 2637   Male     226  42 119.0  80.0      Yes 25.29       No     No      98
    ## 2638 Female     211  36 100.0  61.5       No 22.19       No     No      73
    ## 2639   Male     334  54 133.5  80.0       No 23.40       No     No      77
    ## 2640 Female     197  41 113.0  70.0       No 23.78       No     No      90
    ## 2641   Male     280  60 114.0  82.0       No 23.96       No     No      84
    ## 2643   Male     181  41 125.0  79.0       No 19.09       No     No      70
    ## 2644   Male     224  52 128.0  82.0       No 23.96       No     No      68
    ## 2645   Male     216  53 110.0  79.0      Yes 24.76       No     No      74
    ## 2646 Female     222  37 110.0  71.0      Yes 18.30       No     No      67
    ## 2649   Male     113  38 120.0  83.5      Yes 30.34       No     No      85
    ## 2650 Female     212  40 110.0  70.0       No 22.98       No     No      85
    ## 2651   Male     309  37 107.5  80.0      Yes 25.23       No     No      84
    ## 2652   Male     238  65 122.0  81.0       No 23.95      Yes     No     150
    ## 2655 Female     180  39 112.5  85.0      Yes 25.31       No     No      58
    ## 2656 Female     317  41 149.5  93.0       No 35.42       No     No      87
    ## 2657   Male     201  48 117.5  80.0       No 23.68       No     No      85
    ## 2658   Male     167  50 159.0  95.0       No 25.20       No     No      87
    ## 2659 Female     302  67 121.0  83.0       No 30.12       No     No      64
    ## 2660   Male     228  41 113.0  82.5      Yes 25.67       No     No      70
    ## 2661   Male     226  46 140.0  86.0       No 31.93       No     No      72
    ## 2662 Female     232  43 126.0  73.0       No 20.43       No     No      75
    ## 2663 Female     215  55 122.0  86.0       No 30.61       No     No      87
    ## 2664   Male     214  55 110.0  71.0      Yes 24.24       No     No      72
    ## 2665 Female     266  48 155.0 100.0      Yes 27.86       No     No      84
    ## 2667 Female     234  48 144.0  90.0       No 29.34       No     No      70
    ## 2668 Female     155  38 122.0  81.0       No 27.14       No     No      70
    ## 2669 Female     240  48 123.0  80.0       No 24.57       No     No      74
    ## 2670 Female     178  40 119.0  78.5       No 23.28       No     No      75
    ## 2671 Female     254  57 182.5  97.0       No 27.38       No     No      72
    ## 2672 Female     274  40 132.5  85.5       No 24.87       No     No      70
    ## 2673 Female     266  64 166.0  90.0       No 23.33       No     No      87
    ## 2674 Female     237  53 126.0  84.0       No 27.04       No     No      88
    ## 2675   Male     260  46 137.5  94.0      Yes 32.37       No     No      78
    ## 2676 Female     190  51 131.5  89.0      Yes 23.66       No     No     100
    ## 2677 Female     188  40 123.0  76.5       No 28.00       No     No      76
    ## 2679   Male     200  52 141.0  81.5      Yes 26.56       No     No      85
    ## 2680 Female     167  41 147.5  87.5      Yes 32.52       No     No      80
    ## 2681 Female     248  50 154.5 104.0      Yes 19.88       No     No      87
    ## 2682   Male     176  36 125.0  85.0       No 25.88       No     No      83
    ## 2683 Female     226  51 105.0  71.0       No 27.73       No     No      79
    ## 2684   Male     233  40 122.0  87.0      Yes 24.91       No    Yes      69
    ## 2685 Female     253  45 135.0  85.0       No 24.35       No     No      84
    ## 2686 Female     170  40 142.5  95.0       No 36.79       No     No      62
    ## 2687 Female     273  43 119.0  72.0       No 24.59       No     No      75
    ## 2688   Male     207  54 146.0  98.0       No 23.63       No     No      91
    ## 2689 Female     227  51 160.0  90.0       No 23.48       No    Yes      57
    ## 2690   Male     150  38 109.0  75.0       No 21.43       No     No      67
    ## 2692 Female     204  36 112.0  78.0      Yes 28.74       No     No      82
    ## 2693 Female     215  52 180.0 108.0       No 37.02       No     No      89
    ## 2695 Female     226  39 115.0  80.0      Yes 25.19       No     No      74
    ## 2696   Male     333  58 139.0  96.0       No 28.38       No     No      78
    ## 2698   Male     304  49 147.0 102.0       No 31.67       No     No      77
    ## 2699 Female     350  52 175.0  85.0       No 30.99       No     No      63
    ## 2700   Male     275  40 112.5  85.0      Yes 28.04       No     No      71
    ## 2701 Female     315  55 123.0  77.5      Yes 26.21       No     No      84
    ## 2702   Male     290  58 116.5  83.0      Yes 23.27       No     No      85
    ## 2703   Male     251  48 112.0  77.0       No 30.28       No     No      88
    ## 2704 Female     248  54 131.0  71.0       No 24.89       No     No      73
    ## 2705   Male     231  60 122.0  78.0      Yes 28.09       No     No      67
    ## 2706 Female     180  41 115.0  84.0       No 20.11       No     No      64
    ## 2708 Female     256  43 129.0  86.0      Yes 25.89       No     No      72
    ## 2710 Female     168  49 165.0  99.0       No 27.10       No     No      86
    ## 2711   Male     207  53 102.5  72.5       No 26.50       No     No      95
    ## 2712   Male     246  44 142.0  92.0      Yes 23.85       No     No      65
    ## 2714 Female     295  44 114.0  85.0      Yes 23.10       No     No      84
    ## 2715 Female     216  42 119.0  75.0      Yes 27.01       No     No      73
    ## 2716 Female     249  37 108.0  74.0       No 23.28       No     No      66
    ## 2718   Male     221  45 105.0  70.0       No 23.95       No     No      83
    ## 2719 Female     212  53 177.0  75.0       No 32.29      Yes     No     285
    ## 2721 Female     214  36 107.5  66.0      Yes 22.02       No     No     103
    ## 2722   Male     214  64 155.0  99.0      Yes 22.46       No     No      82
    ## 2723 Female     213  54 144.0  82.0      Yes 29.45       No     No      72
    ## 2724   Male     240  55 135.0  95.0      Yes 28.69       No     No     108
    ## 2725   Male     230  59 182.0 102.0      Yes 25.91      Yes     No     147
    ## 2726   Male     252  55 129.0  84.0      Yes 23.55       No     No      72
    ## 2727 Female     354  60 130.0  82.5       No 26.76       No     No      79
    ## 2728 Female     320  57 125.0  72.5       No 25.27       No     No      82
    ## 2729 Female     193  42 137.5  82.5       No 30.55       No     No      69
    ## 2730   Male     243  59 162.0  91.0      Yes 33.00       No     No      81
    ## 2731 Female     268  58 170.0 104.0       No 24.51       No     No      75
    ## 2732 Female     382  57 133.0  77.0       No 24.27       No     No      81
    ## 2735   Male     231  51 110.0  72.0      Yes 21.19       No     No      81
    ## 2736   Male     288  49 128.0  89.0      Yes 35.96       No     No      75
    ## 2738 Female     211  52 129.0  73.0       No 29.09       No     No     117
    ## 2739   Male     240  41 166.0  71.0      Yes 35.53       No     No      57
    ## 2741   Male     210  44 103.0  78.0      Yes 34.89       No     No      77
    ## 2742   Male     287  40 141.0  86.0      Yes 27.42       No    Yes      94
    ## 2743 Female     205  39 139.5  87.0       No 20.70       No     No      85
    ## 2744 Female     195  38 116.0  72.0      Yes 24.45       No     No      90
    ## 2745   Male     188  64 191.0 106.0       No 37.38       No     No      84
    ## 2746   Male     240  44 134.0  88.0      Yes 31.62       No     No      61
    ## 2747   Male     241  53 113.0  84.0       No 28.27       No     No      73
    ## 2748 Female     308  43 110.0  70.0       No 24.83       No     No      83
    ## 2749 Female     213  53 133.0  89.0       No 17.89       No     No      73
    ## 2750 Female     232  60 152.5  85.0       No 23.03       No     No     123
    ## 2752 Female     225  48 132.0  80.0       No 21.35       No     No      74
    ## 2754 Female     312  50 165.0  90.0       No 30.47       No     No      83
    ## 2755 Female     280  45 116.0  69.0      Yes 26.45       No     No      92
    ## 2756 Female     288  56 123.0  70.0      Yes 18.62       No     No      96
    ## 2757   Male     222  61 139.0  67.0      Yes 22.01       No     No      72
    ## 2758   Male     267  44 131.0  79.0       No 30.32       No     No      79
    ## 2759 Female     202  43 121.5  86.5       No 20.82       No     No      77
    ## 2760 Female     345  59 182.5 103.0      Yes 31.52       No     No      76
    ## 2761 Female     200  54 134.0  87.0       No 26.27       No     No      86
    ## 2762 Female     295  41 118.0  83.0       No 29.01       No     No      75
    ## 2763 Female     173  42 105.0  70.0      Yes 21.98       No     No      79
    ## 2764   Male     215  61 114.0  72.5      Yes 25.86       No     No      61
    ## 2766 Female     260  43 127.0  79.5       No 26.01       No     No      77
    ## 2768 Female     243  47 121.0  61.0      Yes 20.32       No     No     110
    ## 2769   Male     246  55 147.0  86.5       No 24.80       No     No      76
    ## 2770   Male     205  38 125.0  85.0       No 28.03       No     No      84
    ## 2771   Male     219  46 118.0  79.0      Yes 24.17       No     No      90
    ## 2773 Female     262  64 160.0  82.0       No 21.11       No     No      90
    ## 2774 Female     155  40 120.0  73.5      Yes 21.32       No     No      78
    ## 2775 Female     254  60 114.0  78.0       No 33.88       No     No      84
    ## 2776 Female     150  63 152.0  88.0       No 36.54      Yes     No     170
    ## 2777 Female     255  64 147.5  93.0       No 25.16       No     No      84
    ## 2778   Male     279  45 143.5  90.0       No 21.07       No     No      87
    ## 2779   Male     210  47 131.0  74.0      Yes 18.88       No     No      86
    ## 2780   Male     226  54 122.5  77.5       No 27.93       No     No      71
    ## 2781 Female     171  42 111.0  81.0      Yes 29.47       No     No      77
    ## 2783 Female     274  56 147.0  98.0       No 19.20       No    Yes      82
    ## 2784   Male     220  46 119.0  80.0       No 25.86       No     No      95
    ## 2785 Female     208  55 190.0 130.0       No 56.80       No     No      86
    ## 2787   Male     211  44 145.0  88.0      Yes 23.39       No     No      79
    ## 2788   Male     248  52 119.0  82.0      Yes 26.55       No     No     100
    ## 2789 Female     344  65 120.0  75.0       No 25.41       No     No      98
    ## 2790 Female     285  51 151.5  79.0       No 26.91       No     No      76
    ## 2791 Female     158  42 127.5  90.0       No 31.75       No     No      92
    ## 2792   Male     272  42 125.5  80.5      Yes 25.35       No     No      71
    ## 2793   Male     175  58  83.5  58.0      Yes 29.66       No     No     115
    ## 2794 Female     209  40  96.0  67.0      Yes 23.46       No     No      70
    ## 2795   Male     292  41 152.0  73.0      Yes 25.21       No     No     112
    ## 2796 Female     318  57 161.0  90.0       No 38.11      Yes    Yes     112
    ## 2798 Female     230  65 159.0  87.0       No 22.91       No    Yes      65
    ## 2799   Male     276  68 144.0  74.0       No 27.99       No     No      77
    ## 2800 Female     364  59 142.0  84.0       No 26.24       No     No      70
    ## 2801   Male     339  51 137.5  81.0      Yes 24.22       No     No      85
    ## 2802   Male     242  36 115.0  75.0      Yes 25.64       No     No      83
    ## 2803   Male     230  63 152.5  88.0      Yes 24.10       No     No      74
    ## 2804   Male     179  37 131.5  81.0      Yes 24.99       No     No      68
    ## 2805   Male     215  51 115.0  69.0      Yes 25.70       No     No      77
    ## 2806   Male     210  52 128.0  87.0       No 26.25       No     No      61
    ## 2807 Female     250  47  98.0  73.0      Yes 24.39       No     No      88
    ## 2808 Female     307  56 118.0  80.0       No 29.67       No     No      70
    ## 2809 Female     253  39 127.0  68.0       No 23.29       No     No      75
    ## 2810 Female     242  68 146.5  71.0       No 25.21       No     No      96
    ## 2811 Female     291  52 150.0  94.0       No 28.68       No     No      74
    ## 2812 Female     295  65 210.0 135.0       No 24.73       No     No      93
    ## 2813   Male     173  47 117.5  77.5      Yes 20.23       No     No      75
    ## 2814 Female     248  54 131.5  83.0       No 26.84       No     No      63
    ## 2815   Male     194  67 171.0  84.0      Yes 21.64       No     No      88
    ## 2816   Male     204  63 158.0 109.0       No 28.40       No     No      71
    ## 2817 Female     220  58 121.0  75.0       No 25.33       No     No      76
    ## 2819 Female     246  63 116.0  69.0       No 23.44       No     No      78
    ## 2820   Male     207  58 110.0  80.0      Yes 23.55       No     No      78
    ## 2821 Female     177  36 115.0  63.5      Yes 22.54       No     No      73
    ## 2822 Female     208  46 122.0  76.0       No 18.40       No     No      70
    ## 2823   Male     259  64 115.0  81.0      Yes 20.09       No     No      78
    ## 2826 Female     262  56 150.0  89.0       No 40.21       No     No      86
    ## 2827 Female     322  68 148.0  95.0      Yes 20.98       No     No      72
    ## 2829 Female     219  63 124.0  75.0       No 28.57       No     No      76
    ## 2830 Female     288  49  99.0  60.0      Yes 22.19       No     No      76
    ## 2831   Male     335  40 136.5  84.0      Yes 23.60       No     No      76
    ## 2832 Female     192  46 122.0  83.0       No 23.99       No     No     106
    ## 2833 Female     198  51 128.0  90.0       No 26.32       No     No      75
    ## 2834   Male     274  44 150.0  95.0      Yes 25.46       No     No      70
    ## 2835 Female     228  67 144.0  85.0       No 27.59       No    Yes      75
    ## 2836 Female     215  38 106.5  75.0       No 23.82       No     No      67
    ## 2837   Male     263  43 132.0  90.0      Yes 28.85       No     No      71
    ## 2838   Male     285  46 121.0  82.0      Yes 27.62       No     No      79
    ## 2841   Male     309  41 113.0  74.0      Yes 25.40       No     No      64
    ## 2842   Male     242  65 176.0  95.0       No 20.34       No     No     103
    ## 2843 Female     341  54 155.0 105.0       No 27.01       No     No      62
    ## 2844   Male     265  41 130.0  80.0      Yes 26.58       No     No      63
    ## 2845 Female     195  48 121.0  78.0      Yes 26.27       No     No      80
    ## 2846   Male     240  62 152.0  95.0       No 25.37       No     No      70
    ## 2847 Female     282  58 143.0  76.0       No 27.87       No     No      92
    ## 2848   Male     216  67 148.0  70.0      Yes 32.63       No     No     103
    ## 2850 Female     260  41 101.0  68.0      Yes 22.49       No     No      77
    ## 2851   Male     205  37 142.0  80.0      Yes 27.93       No     No     103
    ## 2852 Female     286  62 123.0  77.0      Yes 20.56       No     No      86
    ## 2853   Male     225  44 130.0  77.5      Yes 21.19       No     No      82
    ## 2854 Female     275  48 177.0 101.0       No 25.22       No     No      82
    ## 2855 Female     231  44 103.0  73.0      Yes 22.02       No     No      83
    ## 2856 Female     212  42 116.0  68.0       No 21.49       No     No      79
    ## 2859   Male     238  46 140.5  92.5      Yes 26.97       No     No      83
    ## 2860 Female     211  46 120.0  84.0      Yes 22.53       No     No      87
    ## 2861   Male     300  56 165.0 112.0      Yes 23.68       No     No      78
    ## 2862 Female     232  40 140.0  92.0      Yes 26.56       No     No      73
    ## 2863   Male     190  53 141.0 115.0      Yes 21.01       No    Yes      86
    ## 2864 Female     180  63 170.5 100.5       No 27.69       No     No      63
    ## 2865   Male     236  41 127.0  84.0      Yes 31.12       No     No      83
    ## 2866   Male     211  55 142.0  82.5      Yes 19.81       No     No      70
    ## 2867   Male     212  45 119.0  72.0      Yes 22.33       No     No      91
    ## 2869 Female     249  43 108.0  76.0      Yes 27.49       No     No      76
    ## 2870 Female     200  45 117.5  83.5      Yes 23.68       No     No      73
    ## 2871 Female     159  38 132.0  85.0       No 28.21       No     No      76
    ## 2872   Male     228  40 131.0  80.0      Yes 26.77       No     No      74
    ## 2873 Female     187  40 107.0  74.0      Yes 22.37       No     No      67
    ## 2874   Male     275  43 111.0  76.0      Yes 23.95       No     No      65
    ## 2875   Male     264  58 156.0  86.0      Yes 26.05       No     No     103
    ## 2876 Female     210  45 121.0  82.0      Yes 23.08       No     No      71
    ## 2877   Male     232  59 108.0  73.0      Yes 24.92       No     No      81
    ## 2879 Female     305  64 126.5  67.0       No 25.77       No     No      66
    ## 2880   Male     229  44 132.0  94.0       No 34.39       No     No      80
    ## 2881   Male     285  39 120.0  77.5      Yes 28.33       No     No      79
    ## 2883 Female     205  47 107.0  68.0      Yes 23.29       No     No     108
    ## 2884 Female     175  45 114.0  74.0       No 30.53       No     No     105
    ## 2885   Male     238  45 175.0 114.0      Yes 31.14       No    Yes      70
    ## 2887   Male     293  50 150.0  85.0      Yes 26.09       No     No      68
    ## 2888 Female     208  54 138.0  78.0       No 30.47       No     No      73
    ## 2889 Female     170  41 113.5  65.5       No 31.71       No     No      93
    ## 2890 Female     244  44 132.5  87.0       No 24.17       No     No      86
    ## 2891 Female     205  44 109.0  73.0       No 17.48       No     No      57
    ## 2892   Male     234  45 126.0  80.0      Yes 23.14       No     No      88
    ## 2894   Male     197  41 134.0  87.0       No 25.48       No     No      75
    ## 2895 Female     209  41 107.0  65.0      Yes 27.27       No     No      87
    ## 2896 Female     231  58 165.0  94.5      Yes 27.02       No     No      80
    ## 2897   Male     253  65 111.0  60.0       No 24.12       No     No      69
    ## 2898 Female     226  44 123.0  82.0       No 24.67       No     No      70
    ## 2899 Female     228  41  93.0  71.0       No 31.57       No     No      85
    ## 2900 Female     261  60 123.5  79.0       No 29.28       No     No     103
    ## 2902 Female     240  40 115.0  72.0      Yes 18.82       No     No      68
    ## 2904 Female     280  56 101.0  71.0       No 28.22       No     No      67
    ## 2905 Female     260  52 171.0 118.0       No 28.33       No     No      80
    ## 2906   Male     235  54 132.0  87.0       No 26.13       No     No      75
    ## 2907 Female     278  42 118.5  72.0       No 21.99       No     No      78
    ## 2908   Male     238  42 133.0  88.0      Yes 28.84       No     No      58
    ## 2909   Male     210  43 181.0  97.5       No 21.83       No    Yes      55
    ## 2910 Female     239  49 143.0  93.0       No 28.38       No     No      87
    ## 2911   Male     344  47 102.5  71.0      Yes 27.73       No     No      80
    ## 2912   Male     193  53 142.5 100.0      Yes 24.15       No     No      69
    ## 2913 Female     281  62 126.0  85.5       No 23.29       No     No      63
    ## 2915 Female     189  41 113.5  61.0      Yes 23.08       No     No      73
    ## 2916 Female     179  38 116.5  72.5      Yes 21.49       No     No      76
    ## 2917 Female     243  50 128.0  76.0      Yes 25.65       No     No      80
    ## 2918   Male     230  44 138.0  98.0      Yes 25.12       No     No      93
    ## 2919   Male     274  55 148.0  94.0      Yes 28.55       No     No      65
    ## 2920 Female     268  53 156.0 100.0       No 26.73       No     No      69
    ## 2921 Female     220  54 142.5  83.5       No 25.09       No     No      84
    ## 2923   Male     286  47 148.0  98.0      Yes 29.98       No     No      93
    ## 2924   Male     133  51 138.0  78.0      Yes 16.98       No     No      65
    ## 2925 Female     295  60 149.0 103.5       No 22.69       No     No      52
    ## 2926   Male     280  49 166.0  98.0       No 23.03       No    Yes      72
    ## 2927   Male     367  43 141.0  82.5       No 25.62       No     No      90
    ## 2928   Male     207  57 111.0  80.0       No 37.15       No     No      70
    ## 2929 Female     224  61 155.0  71.0       No 25.98       No     No      86
    ## 2931 Female     160  47 197.0 109.0      Yes 34.91      Yes     No     320
    ## 2932 Female     313  62 164.0  94.0       No 27.40      Yes     No      81
    ## 2933   Male     191  62 156.0  91.0       No 31.20       No     No      75
    ## 2934 Female     219  35 116.0  74.0       No 21.63       No     No      68
    ## 2936   Male     212  51 146.0  89.0      Yes 24.49       No     No     132
    ## 2937   Male     241  42 118.0  85.0      Yes 30.03       No     No      74
    ## 2938 Female     215  40 110.0  70.0       No 19.64       No     No      87
    ## 2939 Female     312  62 162.5  93.5       No 25.33       No    Yes      83
    ## 2940 Female     188  44 122.5  67.5      Yes 23.13       No     No      85
    ## 2941 Female     214  46 144.0  79.0       No 23.62       No     No      80
    ## 2942 Female     193  40 105.0  60.0      Yes 22.54       No     No      85
    ## 2943   Male     201  40 121.5  76.5      Yes 29.38       No     No      94
    ## 2944   Male     214  48 114.0  80.0      Yes 26.27       No     No      99
    ## 2945 Female     208  43 119.0  69.0      Yes 19.48       No     No      73
    ## 2946   Male     194  48 102.5  69.0      Yes 18.55       No     No      83
    ## 2947 Female     230  43 116.0  86.0       No 27.78       No     No      78
    ## 2948   Male     265  38 123.0  86.0      Yes 28.88       No     No      80
    ## 2949   Male     263  48 114.0  80.0       No 25.14       No     No      77
    ## 2950   Male     279  58 181.0  74.0       No 22.49       No     No      65
    ## 2951 Female     304  51 132.0  80.0       No 25.15       No     No     106
    ## 2952   Male     218  42 129.0  90.0      Yes 23.22       No     No      73
    ## 2953   Male     272  46 145.0  88.0      Yes 28.12       No     No      65
    ## 2954   Male     285  48 116.0  81.0      Yes 23.10       No     No      58
    ## 2955 Female     240  55 107.5  70.0      Yes 18.06       No     No     140
    ## 2957 Female     232  41 117.5  77.5       No 20.62       No     No      75
    ## 2958 Female     252  57 139.0  82.0      Yes 26.36       No     No      70
    ## 2959   Male     237  51 135.0  87.5       No 24.87       No     No      63
    ## 2960   Male     220  42 131.0  90.0      Yes 24.21       No     No      79
    ## 2961 Female     177  41 108.0  73.0       No 24.73       No     No      72
    ## 2962   Male     231  67 138.0  76.0       No 18.99       No     No      46
    ## 2963 Female     230  45 138.0  85.0      Yes 22.54       No     No      67
    ## 2964 Female     219  61 102.0  59.0       No 18.14       No     No      62
    ## 2965 Female     185  43 125.0  84.0      Yes 23.18       No     No      55
    ## 2966   Male     230  54 126.0  93.0       No 25.36       No     No      84
    ## 2968 Female     297  59 198.0 108.0       No 28.06       No    Yes     109
    ## 2969   Male     203  40 123.0  82.0      Yes 24.74       No     No      78
    ## 2970   Male     173  60 135.0  97.5       No 30.71       No    Yes      86
    ## 2971   Male     282  49 134.0  79.0      Yes 26.87       No     No      78
    ## 2972 Female     220  52 163.0  87.0       No 25.12       No     No      94
    ## 2973   Male     264  51 119.0  82.0      Yes 28.55       No     No      63
    ## 2974 Female     257  49 115.0  75.0       No 25.61       No     No      74
    ## 2975 Female     208  52 141.0  85.0       No 24.85       No     No      68
    ## 2976   Male     184  42 114.5  76.0      Yes 29.89       No     No      59
    ## 2977   Male     237  44 124.0  83.0      Yes 27.17       No     No      88
    ## 2980   Male     192  62 163.0  83.0      Yes 20.99       No     No      97
    ## 2981 Female     262  54 113.5  80.0       No 26.42       No     No      67
    ## 2982   Male     172  48 131.0  79.0      Yes 35.12      Yes     No     108
    ## 2983   Male     240  52 135.0  82.0      Yes 27.73       No     No      87
    ## 2984 Female     263  42 150.0  88.0       No 23.68       No     No      78
    ## 2985 Female     280  50 126.0  76.0       No 21.82       No     No      74
    ## 2986   Male     165  46 107.0  82.0      Yes 24.67       No     No      77
    ## 2987   Male     188  59 111.0  72.0      Yes 21.48       No     No      44
    ## 2988   Male     314  59 135.0  77.5       No 22.17      Yes     No     170
    ## 2989   Male     320  50 190.0 110.0      Yes 25.45       No     No      71
    ## 2990 Female     153  39 107.0  82.0       No 20.99       No     No      82
    ## 2992   Male     255  58 136.0  85.5      Yes 25.33       No     No      62
    ## 2993 Female     205  47 106.0  75.0       No 23.44       No     No      97
    ## 2994 Female     163  57 112.0  71.0       No 22.66       No     No      78
    ## 2995 Female     281  63 125.0  80.0       No 21.35       No     No      99
    ## 2996 Female     194  43 122.5  82.0       No 35.16       No     No      80
    ## 2997   Male     226  56 175.5 113.0       No 30.77       No    Yes      70
    ## 2998   Male     163  50  96.5  72.5      Yes 20.72       No     No      77
    ## 2999   Male     214  47 132.0  95.0      Yes 24.23       No     No      77
    ## 3000 Female     215  48 114.0  64.0      Yes 21.51       No     No      64
    ## 3001   Male     212  56 130.0  79.0       No 27.73       No     No     100
    ## 3002   Male     152  40 119.0  86.0      Yes 23.35       No     No      66
    ## 3003   Male     253  63 172.0  82.0      Yes 24.19       No     No     137
    ## 3004 Female     280  52 118.0  71.0       No 21.84       No     No      62
    ## 3006 Female     240  59 149.0  88.0       No 27.48       No     No      82
    ## 3007 Female     228  42 122.0  85.0      Yes 23.99       No     No      68
    ## 3008   Male     264  41 126.5  82.0      Yes 23.96       No     No      78
    ## 3009 Female     247  48 142.5  86.5      Yes 24.87       No     No      68
    ## 3010   Male     291  49 160.0  99.0      Yes 29.91       No     No      88
    ## 3011   Male     175  43 103.0  64.5      Yes 17.81       No     No      78
    ## 3012 Female     280  53 135.0  87.0       No 27.17       No     No      80
    ## 3013   Male     211  49 128.0  89.0      Yes 31.07       No     No      76
    ## 3014 Female     260  67 140.0  88.0       No 24.67       No     No      65
    ## 3016 Female     234  59 181.0 107.0       No 39.40       No    Yes      90
    ## 3017   Male     238  46 163.5 102.0      Yes 28.65       No     No     100
    ## 3019 Female     205  47 122.0  78.0       No 23.78       No     No      83
    ## 3020 Female     201  36 124.0  73.0      Yes 25.25       No     No     100
    ## 3021 Female     292  51 166.0  89.0       No 27.51       No     No      88
    ## 3022 Female     189  49 102.0  66.5       No 23.88       No     No      70
    ## 3023 Female     241  54 103.0  66.0       No 22.73       No     No      66
    ## 3024   Male     248  49 130.5  82.0       No 27.29      Yes     No     254
    ## 3025   Male     144  67 113.0  62.0      Yes 22.90       No     No      69
    ## 3026   Male     346  62 102.5  66.5       No 17.17      Yes     No     394
    ## 3027 Female     224  38  90.0  70.0      Yes 18.18       No     No      57
    ## 3029 Female     250  50 109.0  70.0       No 20.74       No     No      77
    ## 3030 Female     290  45 124.0  72.5       No 24.24       No     No      87
    ## 3031   Male     178  48 132.5  73.0      Yes 24.67       No     No      70
    ## 3032   Male     231  58 171.0  95.0      Yes 26.11       No     No      77
    ## 3033 Female     203  36 112.5  73.5       No 24.47       No     No      73
    ## 3036 Female     174  57 120.0  62.0      Yes 25.13       No     No      77
    ## 3039   Male     206  43 133.5  89.5      Yes 21.73       No     No      89
    ## 3040 Female     179  63 149.0  89.0       No 40.81       No     No      77
    ## 3041 Female     229  59 127.5  76.0      Yes 23.65       No     No      70
    ## 3042 Female     303  67 204.0  96.0       No 27.86      Yes    Yes     394
    ## 3043 Female     223  42 129.0  79.0      Yes 28.04       No     No     100
    ## 3045 Female     252  39 125.0  86.0       No 20.74       No     No      93
    ## 3046 Female     183  42 120.0  76.0       No 21.12       No     No      72
    ## 3047 Female     230  43 118.0  70.5      Yes 26.24       No     No      67
    ## 3048 Female     263  52 132.5  87.0       No 30.42       No     No      64
    ## 3049 Female     210  45 120.0  72.0       No 22.01       No     No      93
    ## 3050   Male     243  61 142.0  89.0       No 27.30       No     No      67
    ## 3051 Female     192  56 122.0  82.5       No 28.61       No     No      58
    ## 3052   Male     296  43 137.0  90.0      Yes 23.97       No     No      97
    ## 3053   Male     228  56 151.5 103.0      Yes 23.58       No     No      66
    ## 3054   Male     345  55 134.0  89.0      Yes 27.38       No     No      60
    ## 3055   Male     273  46 132.5  69.0      Yes 26.83       No     No      88
    ## 3056   Male     252  43 112.0  78.0      Yes 24.25       No     No      65
    ## 3057   Male     206  43 107.5  73.5       No 24.17       No     No      71
    ## 3059 Female     258  50 127.0  80.0       No 24.10      Yes     No     124
    ## 3060   Male     230  48 140.5  89.0      Yes 23.34       No     No      80
    ## 3063 Female     283  57 207.5 118.0      Yes 38.61       No     No      83
    ## 3064   Male     218  51 115.0  75.0       No 25.62       No     No      62
    ## 3065 Female     280  40 115.0  81.0      Yes 21.32       No     No      84
    ## 3066   Male     240  51 127.0  80.0      Yes 28.85       No     No      67
    ## 3067   Male     276  39 129.5  85.5      Yes 24.74       No     No      85
    ## 3068   Male     163  41 102.0  59.5      Yes 22.27       No     No      82
    ## 3069 Female     254  59 181.0 101.0       No 24.67       No     No      68
    ## 3070 Female     295  61 191.0  97.0       No 31.27       No     No      90
    ## 3071 Female     276  43 138.0  88.0      Yes 20.40       No     No      70
    ## 3072 Female     225  65 145.0  91.0       No 29.80       No     No      83
    ## 3073 Female     176  42 136.0  87.0       No 30.91       No     No      73
    ## 3074   Male     217  41 110.0  79.0      Yes 21.27       No     No      95
    ## 3076   Male     220  56 113.0  83.0       No 25.12       No     No      86
    ## 3079 Female     287  56 150.5  94.0       No 26.01       No     No      99
    ## 3081 Female     279  59 132.0  88.0       No 26.48       No     No      88
    ## 3082   Male     217  54 116.0  72.0       No 22.49       No     No      77
    ## 3084 Female     254  45 104.0  73.0       No 19.46       No     No      77
    ## 3085 Female     328  53 121.0  70.0       No 30.75       No     No      77
    ## 3086 Female     306  54 195.0 110.0       No 32.10       No     No      68
    ## 3087 Female     281  59 130.5  86.0       No 25.11       No     No      83
    ## 3088   Male     260  41 137.5  80.0      Yes 26.89       No     No      55
    ## 3089 Female     237  47 110.0  77.0      Yes 25.62       No     No      83
    ## 3091 Female     160  44 131.0  81.5      Yes 25.71       No     No      70
    ## 3092   Male     180  54 147.5 100.0       No 25.11       No     No      70
    ## 3093 Female     180  55 165.0  88.0       No 22.57       No     No      77
    ## 3094   Male     260  53 142.0  88.0      Yes 23.65       No     No      77
    ## 3095   Male     258  52 177.0 111.0       No 30.38      Yes     No     270
    ## 3096   Male     274  51 108.0  75.0      Yes 23.60       No     No      68
    ## 3097 Female     228  49 124.0  72.0       No 21.18       No     No      88
    ## 3098   Male     193  68 145.0  67.0      Yes 23.13       No     No      72
    ## 3099   Male     202  39 108.0  74.0      Yes 25.51       No     No     104
    ## 3100   Male     248  54 155.0  92.5       No 29.86       No     No      66
    ## 3102   Male     237  46 105.0  72.5      Yes 33.49       No     No      79
    ## 3103 Female     262  49 147.5  97.5      Yes 24.96       No    Yes      67
    ## 3104 Female     229  46 127.5  80.0       No 22.34       No     No      77
    ## 3106 Female     352  51 136.5  87.0      Yes 25.79       No     No      67
    ## 3107   Male     280  55 134.0  85.5      Yes 29.86       No     No      75
    ## 3108 Female     270  61 174.5 101.5       No 29.87       No     No      76
    ## 3109   Male     220  38 126.0  88.0      Yes 24.46       No     No      74
    ## 3110   Male     252  36 126.0  77.0      Yes 25.23       No     No      63
    ## 3111   Male     225  50 133.5  86.0      Yes 29.38       No     No     106
    ## 3112   Male     234  42 125.0  86.0      Yes 28.09       No     No      70
    ## 3113 Female     258  45 148.5  88.0       No 23.46       No     No      90
    ## 3114 Female     280  48 105.0  85.0      Yes 25.50       No     No      79
    ## 3115   Male     233  50 158.0  88.0      Yes 28.26       No     No      94
    ## 3116 Female     265  52 130.0  84.0       No 27.09       No     No      69
    ## 3117   Male     258  50 198.0 106.0      Yes 26.73       No    Yes      73
    ## 3118   Male     241  58 152.5 105.0      Yes 25.18       No     No      85
    ## 3119 Female     432  57 153.0  85.0       No 26.13       No    Yes      75
    ## 3122   Male     189  44 130.0  90.0      Yes 22.33       No     No      64
    ## 3123 Female     230  40 107.5  75.0      Yes 26.38       No     No      76
    ## 3124 Female     240  58 150.0  84.0      Yes 26.85       No     No      94
    ## 3125 Female     252  45 160.0 105.0       No 31.72       No     No      83
    ## 3126 Female     172  38  98.0  53.0      Yes 22.18       No     No      82
    ## 3127   Male     239  61 122.0  83.0       No 28.85       No     No      94
    ## 3129 Female     317  62 126.0  75.0       No 23.29       No     No      86
    ## 3130 Female     229  50 121.0  85.5       No 23.09       No     No      75
    ## 3131   Male     200  36 103.0  67.5      Yes 30.82       No     No      72
    ## 3132   Male     234  44 135.0  88.0       No 25.65       No     No      83
    ## 3133   Male     260  44 127.0  80.0       No 31.69       No     No      73
    ## 3134   Male     232  44 125.0  86.0      Yes 23.23       No     No      80
    ## 3136 Female     216  55 125.0  80.0      Yes 27.18      Yes     No     244
    ## 3137   Male     258  48 150.0 105.0      Yes 25.94       No     No      60
    ## 3138   Male     185  64  94.0  62.0       No 26.11       No     No      68
    ## 3139 Female     181  53 163.5  87.0      Yes 34.69       No     No      71
    ## 3140 Female     217  37 126.0  80.0       No 25.91       No     No      62
    ## 3142   Male     210  40 132.0  86.0      Yes 31.57       No     No      80
    ## 3143   Male     264  45 118.5  81.0      Yes 26.35       No     No      90
    ## 3144   Male     188  39 105.0  65.0       No 22.85       No     No      76
    ## 3145 Female     204  55 131.0  87.0       No 24.68       No     No      70
    ## 3147   Male     191  42 118.0  80.0       No 24.98       No     No      63
    ## 3148   Male     200  52 131.0  94.0      Yes 26.77       No     No      82
    ## 3149 Female     165  35 106.0  64.0      Yes 19.14       No     No      70
    ## 3150 Female     216  39 126.0  86.0       No 22.72       No     No      70
    ## 3151 Female     195  54 122.0  88.5       No 24.11       No     No      67
    ## 3152 Female     240  56 125.0  79.0      Yes 27.38       No     No      82
    ## 3153 Female     273  49 147.0  89.0      Yes 24.26       No     No      62
    ## 3154   Male     273  39 148.0 103.0      Yes 28.62       No     No      74
    ## 3155   Male     276  53 130.0  86.0      Yes 24.21       No     No      82
    ## 3156 Female     250  59 120.5  80.5       No 27.59       No     No      75
    ## 3157 Female     300  52 137.0  79.0       No 25.40       No     No      73
    ## 3158   Male     266  45 109.0  70.0      Yes 23.72       No     No      71
    ## 3159   Male     239  53 122.0  82.0      Yes 29.40       No     No      84
    ## 3160 Female     273  66 197.0  91.0       No 23.22       No     No      80
    ## 3161 Female     261  41 147.5  97.0      Yes 31.65       No     No      73
    ## 3162 Female     282  45 115.0  71.0       No 22.54       No     No      96
    ## 3163 Female     290  46 131.0  84.0      Yes 18.28       No     No      68
    ## 3164 Female     232  41 125.0  80.0      Yes 25.38       No     No      79
    ## 3165 Female     313  59 186.5  99.0       No 25.65       No     No      84
    ## 3166 Female     310  44 141.0  76.0      Yes 20.52       No     No      90
    ## 3167   Male     174  64 119.0  80.0      Yes 20.03       No     No      78
    ## 3170   Male     225  53 108.5  73.0       No 26.55       No     No      73
    ## 3171 Female     226  34 112.5  77.5       No 24.99       No     No      72
    ## 3172   Male     226  38 142.0  93.0      Yes 29.76       No     No      63
    ## 3173   Male     172  53 131.0  81.0       No 26.30       No     No      84
    ## 3174   Male     195  44 114.0  79.0       No 25.01       No     No      76
    ## 3175 Female     271  45 121.0  77.0       No 23.96       No     No      82
    ## 3176   Male     287  57 149.0  86.0       No 26.33       No     No      65
    ## 3177   Male     257  61 124.5  80.0       No 23.84       No     No      76
    ## 3178 Female     177  38 126.0  80.0      Yes 23.84       No     No      79
    ## 3179   Male     235  57 123.0  78.0       No 28.53       No     No      74
    ## 3180 Female     266  47 107.0  77.0       No 30.61       No     No      64
    ## 3183   Male     239  58 147.0  98.0       No 26.51       No     No      68
    ## 3185 Female     270  60 130.0  72.5      Yes 20.84       No     No     102
    ## 3186 Female     241  63 153.0  89.0       No 32.57       No     No      75
    ## 3188   Male     243  56 128.5  75.0       No 27.68       No     No      89
    ## 3189 Female     278  54 127.5  83.0       No 23.08       No     No      78
    ## 3190 Female     315  57 193.0 109.0      Yes 27.99       No    Yes      74
    ## 3192 Female     351  44 106.0  65.0      Yes 25.34       No     No      69
    ## 3193   Male     174  60 164.0 113.0      Yes 18.64       No     No      93
    ## 3194   Male     196  41 154.0  99.0      Yes 23.46       No     No      56
    ## 3195 Female     275  50 123.0  83.0      Yes 24.29       No     No      64
    ## 3196 Female     254  49 131.0  92.0       No 29.22       No     No      65
    ## 3197   Male     306  48 164.0 120.0      Yes 26.44       No    Yes      74
    ## 3198   Male     239  50 148.5 100.0       No 27.83       No     No      90
    ## 3200 Female     280  52 127.5  81.5       No 33.90       No     No      70
    ## 3201   Male     165  33 136.0  75.0       No 24.95       No     No      90
    ## 3202   Male     328  46 130.0  86.5       No 27.17       No     No      82
    ## 3203   Male     243  52 116.0  64.0      Yes 23.78       No     No      70
    ## 3204   Male     217  45 130.0  82.0      Yes 26.95       No     No      85
    ## 3205 Female     196  46 114.0  75.0      Yes 21.01       No     No      69
    ## 3206   Male     256  63 138.5  77.5      Yes 28.85       No     No      84
    ## 3207 Female     192  41 123.0  72.0      Yes 19.16       No     No      90
    ## 3208 Female     193  41 127.0  83.0       No 21.49       No     No      80
    ## 3209   Male     214  52  98.0  67.0       No 23.43       No     No      90
    ## 3210   Male     297  46 133.0  92.0       No 25.98       No     No      64
    ## 3211 Female     184  63 196.0 101.0       No 28.27       No     No      82
    ## 3212   Male     257  37 105.0  70.0       No 31.47       No     No      80
    ## 3213 Female     269  43 139.0  96.0      Yes 24.38       No     No      71
    ## 3214 Female     223  67 111.0  73.0       No 27.89       No     No      63
    ## 3215 Female     245  42 117.5  90.0       No 22.65       No     No      96
    ## 3216 Female     177  51 141.0  92.0       No 29.64       No     No     130
    ## 3217 Female     289  63 158.0  80.0       No 32.66       No     No      84
    ## 3218   Male     163  36 111.0  73.0       No 30.18       No     No      90
    ## 3219   Male     260  51 127.5  92.0       No 32.98       No     No      93
    ## 3220   Male     177  43 147.5  92.0       No 24.81       No     No      90
    ## 3221 Female     220  40 131.5  82.5      Yes 24.35       No     No      78
    ## 3222   Male     278  49 152.0  93.0      Yes 29.76       No     No      63
    ## 3223   Male     167  57 127.0  69.0      Yes 28.13       No     No     107
    ## 3224   Male     216  39 122.5  77.0       No 24.06       No     No      67
    ## 3225 Female     274  58 159.0  90.0       No 28.40       No    Yes      81
    ## 3228   Male     204  53 152.0  74.0      Yes 24.80       No     No      89
    ## 3229 Female     298  60 133.0  89.0       No 25.09       No     No      81
    ## 3230   Male     192  40 141.5 108.5       No 36.01       No     No      78
    ## 3231 Female     214  49 172.0 111.0      Yes 40.51       No     No      70
    ## 3232 Female     254  52 114.0  80.0      Yes 16.59       No     No      74
    ## 3233   Male     249  57 139.0  83.0      Yes 23.81       No     No      76
    ## 3234 Female     303  46 115.0  78.0      Yes 22.50       No     No      67
    ## 3236 Female     230  40 127.0  72.0       No 29.62       No     No      70
    ## 3237 Female     252  54 146.0  87.0       No 33.11       No     No     100
    ## 3238 Female     262  51 150.5  95.0       No 30.79       No     No      77
    ## 3239   Male     177  54 162.5  99.5      Yes 22.97       No     No      93
    ## 3240 Female     236  51 115.0  71.0       No 23.48       No     No      84
    ## 3241 Female     337  39 120.0  79.5       No 23.07       No     No      77
    ## 3243 Female     321  47 132.0  88.0      Yes 28.14       No     No      74
    ## 3244   Male     248  63 130.0  73.5       No 23.50       No     No      83
    ## 3245   Male     285  39 121.0  82.0      Yes 27.62       No     No      65
    ## 3246 Female     256  40 124.0  79.0      Yes 24.23       No     No      81
    ## 3249   Male     210  47 163.5  97.0      Yes 28.24      Yes     No     183
    ## 3251   Male     275  38 117.5  85.0      Yes 28.94       No     No      74
    ## 3252   Male     238  44 132.0  86.0       No 27.22       No     No      85
    ## 3253 Female     195  45 129.5  83.0       No 25.74       No     No      74
    ## 3254 Female     268  53 163.0  97.0      Yes 27.88       No     No      65
    ## 3256   Male     260  54 116.0  86.0      Yes 24.05       No     No      85
    ## 3257 Female     248  48 114.5  70.0      Yes 27.69       No     No      75
    ## 3258 Female     267  38 179.5  97.0      Yes 20.44       No     No      67
    ## 3259 Female     240  54 125.0  87.0       No 28.76       No     No      76
    ## 3261 Female     198  51 142.5  80.0       No 23.86       No     No     100
    ## 3262   Male     221  39 126.0  80.0       No 23.90       No     No      80
    ## 3263   Male     182  46 120.0  78.0      Yes 20.23       No     No      85
    ## 3264 Female     213  44 126.0  76.0       No 19.05       No     No      84
    ## 3265   Male     274  48 134.0  84.0       No 31.78       No     No      80
    ## 3266 Female     260  42 118.5  74.5       No 22.19       No     No      75
    ## 3267   Male     222  52 125.0  76.0       No 25.23       No     No      77
    ## 3268 Female     326  62 160.0  90.0       No 33.70       No     No      84
    ## 3269   Male     279  54 155.0 105.0      Yes 26.47       No     No      66
    ## 3270 Female     280  64 133.0  82.0       No 28.92       No     No      65
    ## 3271   Male     272  65 137.0  81.0       No 25.74       No     No      97
    ## 3272 Female     286  49 144.0  91.0       No 29.35       No     No      67
    ## 3273 Female     254  60 177.0 101.0      Yes 23.27       No     No      79
    ## 3274 Female     231  70 136.0  84.0       No 31.78       No    Yes      95
    ## 3275   Male     273  44 114.0  83.0      Yes 27.33       No     No      65
    ## 3276 Female     226  40 138.0  99.0      Yes 35.02       No     No      73
    ## 3277 Female     221  47 140.0  94.0       No 28.84      Yes     No      85
    ## 3278 Female     243  39 116.0  79.0       No 22.44       No     No      82
    ## 3279 Female     257  57 133.0  94.0       No 17.71       No     No      75
    ## 3280 Female     275  58 140.0  78.0      Yes 19.18       No     No      74
    ## 3281 Female     295  67 170.0  89.0       No 35.35       No    Yes      63
    ## 3282 Female     310  42 124.0  72.5      Yes 22.32       No     No      74
    ## 3284 Female     246  52 136.0  84.0       No 20.15       No     No      86
    ## 3287 Female     305  50 151.0 106.5       No 25.38       No     No      75
    ## 3288   Male     293  41 132.0  86.0      Yes 24.62       No     No      75
    ## 3289 Female     193  38 107.0  73.0      Yes 20.73       No     No      72
    ## 3291 Female     219  57 125.0  83.0       No 23.22       No     No     115
    ## 3292   Male     200  35 131.0  87.0      Yes 23.93       No     No      70
    ## 3293   Male     270  42 122.0  82.0      Yes 30.31       No     No     142
    ## 3294   Male     188  37 123.5  77.0      Yes 26.62       No     No      80
    ## 3297   Male     234  35 122.5  76.5      Yes 25.16       No     No      85
    ## 3298   Male     245  40 139.0  74.0      Yes 25.48       No     No      93
    ## 3300 Female     240  62 154.0  92.0       No 29.49       No     No      67
    ## 3301   Male     234  42 103.0  71.0      Yes 21.88       No     No      82
    ## 3302 Female     316  60 212.0 116.0       No 27.39       No     No      76
    ## 3303   Male     213  60 123.0  73.0       No 23.65       No     No      79
    ## 3305   Male     292  39 153.0 100.0      Yes 28.09       No     No      69
    ## 3306   Male     363  44 140.0  87.0      Yes 26.44       No     No      79
    ## 3307 Female     226  59 108.0  72.0       No 24.87       No     No      86
    ## 3308   Male     217  60 167.0 109.0      Yes 24.86       No     No      72
    ## 3309 Female     265  53 132.0  80.0       No 26.25       No     No      76
    ## 3310 Female     189  40 115.0  81.0      Yes 22.73       No     No     103
    ## 3311   Male     253  53 115.0  81.0       No 28.09       No     No      84
    ## 3312   Male     201  41 122.0  76.5      Yes 23.81       No     No      73
    ## 3313 Female     279  60 140.5  89.0       No 22.43       No     No      69
    ## 3314 Female     209  36 107.0  73.5       No 21.59       No     No      73
    ## 3315 Female     311  45 117.5  76.0      Yes 26.27       No     No      67
    ## 3316 Female     256  50 136.5  81.0      Yes 23.07       No     No      78
    ## 3317   Male     195  43 113.0  81.0      Yes 34.32       No     No      87
    ## 3318   Male     167  60 162.0  90.0       No 26.36       No     No      82
    ## 3319   Male     248  45 121.0  72.0      Yes 27.88       No     No      88
    ## 3320 Female     278  49 131.0  93.0       No 31.40       No     No      66
    ## 3321   Male     223  42 119.5  87.0      Yes 23.56       No     No      73
    ## 3322 Female     331  55 128.0  84.0       No 21.18       No     No      76
    ## 3323   Male     271  45 164.0  98.0      Yes 26.05       No     No      81
    ## 3324 Female     267  53 152.0  89.0       No 28.49       No     No     103
    ## 3325 Female     188  39 113.0  81.0      Yes 26.44       No     No      87
    ## 3326 Female     226  38 117.5  72.0      Yes 20.71       No     No      73
    ## 3327   Male     252  42 151.5  95.0       No 25.74       No     No      73
    ## 3328 Female     235  48 120.0  81.0      Yes 23.36       No     No      80
    ## 3329   Male     164  62 126.0  79.0       No 27.89       No     No      69
    ## 3330   Male     190  63 148.0  90.0       No 27.13       No     No      86
    ## 3331 Female     165  49 120.0  66.5       No 21.45       No     No      74
    ## 3334 Female     207  43  95.5  70.0      Yes 19.78       No     No      79
    ## 3335 Female     155  45 116.0  81.5       No 25.56       No     No      82
    ## 3336 Female     220  38 105.0  69.0      Yes 24.69       No     No      87
    ## 3337   Male     187  46 140.0  94.0      Yes 24.30       No     No      67
    ## 3338   Male     195  48 124.0  80.0       No 23.96       No     No      73
    ## 3339   Male     218  45 130.0  87.0      Yes 24.10       No     No      64
    ## 3340 Female     225  59 132.0  80.0       No 28.41       No     No      88
    ## 3341   Male     237  43 104.0  74.0       No 23.02       No     No      77
    ## 3342 Female     297  45 142.0  91.0      Yes 35.02       No     No      86
    ## 3343   Male     195  46 111.0  72.0       No 25.72       No     No     108
    ## 3344 Female     277  49 120.0  80.0      Yes 19.72       No     No      75
    ## 3345   Male     165  38 128.0  80.0      Yes 25.62       No     No      85
    ## 3346 Female     214  47 144.0  92.0      Yes 22.73      Yes     No      57
    ## 3347   Male     275  45 105.0  86.0       No 32.92       No     No      92
    ## 3348 Female     224  49 140.0  88.0       No 23.79       No     No      86
    ## 3349   Male     253  37 122.0  81.0       No 27.98       No     No      93
    ## 3350   Male     208  49 115.0  83.0      Yes 24.85       No     No      73
    ## 3352 Female     184  37 137.5  88.5       No 16.48       No     No      68
    ## 3353   Male     216  45 127.5  82.5      Yes 27.15       No     No      67
    ## 3354   Male     206  57 124.0  80.0       No 23.04       No     No      60
    ## 3355 Female     313  46 170.0 107.0       No 27.78       No     No      62
    ## 3356   Male     208  48 138.5  85.5      Yes 23.85       No     No      64
    ## 3357   Male     210  46 154.0  91.0      Yes 26.25       No     No      82
    ## 3358 Female     403  49 170.0 107.5       No 38.43       No     No      93
    ## 3359 Female     194  57 199.5 107.0       No 26.84       No     No      69
    ## 3360 Female     226  55 141.0  84.0       No 27.60       No     No      73
    ## 3361 Female     219  39 161.0 103.0       No 29.06       No     No     106
    ## 3362 Female     196  40 126.0  85.0       No 22.81       No     No      96
    ## 3363 Female     225  42 111.0  71.0      Yes 23.43       No     No      85
    ## 3364 Female     182  65 170.0  81.0       No 21.28       No     No      82
    ## 3366 Female     279  65 152.0 102.0       No 30.43       No    Yes      78
    ## 3367 Female     195  36 109.0  69.0      Yes 23.24       No     No      70
    ## 3368 Female     200  44 128.0  82.0       No 23.24       No     No      73
    ## 3369 Female     258  50 123.0  70.0      Yes 19.72       No     No      71
    ## 3370   Male     168  69 110.0  68.0      Yes 19.15       No    Yes      54
    ## 3371 Female     202  43 124.0  92.0       No 21.26       No     No      74
    ## 3372 Female     236  50 152.0  92.0       No 24.47       No     No      67
    ## 3374 Female     300  43 120.0  78.0      Yes 28.18       No     No     106
    ## 3375 Female     172  41 108.0  73.0      Yes 22.50       No     No      73
    ## 3376 Female     290  55 132.0  81.0       No 27.86       No     No      92
    ## 3377   Male     187  58 136.5  84.0      Yes 25.08       No     No      68
    ## 3378 Female     315  46 165.0  85.0      Yes 32.89       No     No      91
    ## 3379 Female     201  45 128.5  73.0       No 24.56       No     No      82
    ## 3380 Female     263  66 193.0  95.0       No 27.07       No     No      94
    ## 3382 Female     204  51 119.0  80.0       No 24.03       No     No      87
    ## 3383 Female     255  46 117.0  77.0       No 26.35       No     No      76
    ## 3384 Female     200  43 115.0  80.0      Yes 26.66       No     No      81
    ## 3385   Male     327  42 134.0  93.0       No 25.14       No     No      72
    ## 3386 Female     272  64 131.0  85.0      Yes 21.82       No     No      80
    ## 3387   Male     230  39 134.0  92.0      Yes 28.44       No     No      56
    ## 3388 Female     285  45 132.5  97.5      Yes 24.74       No     No      77
    ## 3389   Male     281  52 133.0  93.0      Yes 32.27      Yes     No      80
    ## 3390   Male     258  50 129.0  84.0       No 24.56       No     No      75
    ## 3391   Male     235  60 146.0  92.0       No 26.85       No     No     108
    ## 3392   Male     265  62 135.0  80.0       No 27.94       No     No      80
    ## 3393   Male     230  51 115.5  85.0       No 30.75       No     No      85
    ## 3394   Male     212  51 122.5  73.0      Yes 22.34       No     No      74
    ## 3395 Female     272  40 123.0  75.0      Yes 23.08       No     No      63
    ## 3396 Female     275  39 131.0  80.0       No 25.79       No     No      84
    ## 3397   Male     281  47 119.0  81.0      Yes 23.72       No     No      65
    ## 3398 Female     259  48 129.0  81.0      Yes 21.08       No     No      65
    ## 3399   Male     205  40 131.0  81.0      Yes 23.74       No     No      87
    ## 3400 Female     202  43 114.0  78.0      Yes 26.61       No     No      87
    ## 3401   Male     342  51 110.0  70.0      Yes 28.86       No     No      87
    ## 3402   Male     239  59 168.0 100.0      Yes 27.87       No     No      96
    ## 3403 Female     254  44 145.0  85.0       No 21.27      Yes     No     137
    ## 3404 Female     239  64 143.0  84.0      Yes 20.06       No     No      73
    ## 3406 Female     165  37 108.0  75.0       No 21.84       No     No      83
    ## 3407   Male     265  54 121.0  82.0       No 23.52       No     No      67
    ## 3408   Male     227  60 122.0  80.0      Yes 25.64       No     No      80
    ## 3409 Female     286  69 117.0  73.0       No 20.92       No     No     103
    ## 3410   Male     226  36 124.5  84.5      Yes 21.63       No     No      74
    ## 3411 Female     261  51 127.0  81.0      Yes 20.24       No     No      96
    ## 3412 Female     224  40 106.0  72.0      Yes 23.59       No     No      71
    ## 3413 Female     258  59 166.0  88.0       No 30.66       No     No      78
    ## 3414   Male     225  58 105.5  74.0       No 25.68       No     No      93
    ## 3415 Female     250  55 196.0 103.0      Yes 27.53       No     No      75
    ## 3416   Male     316  44 135.0  95.0      Yes 25.48       No     No      68
    ## 3417   Male     216  44 113.5  77.0       No 28.23       No     No      77
    ## 3419   Male     197  42 102.0  70.5      Yes 24.68       No     No      45
    ## 3420 Female     305  48 118.5  73.0      Yes 20.99       No     No      84
    ## 3421   Male     185  46 121.0  85.0       No 31.31       No     No      97
    ## 3422 Female     194  44 117.0  78.0       No 24.94       No     No      79
    ## 3423 Female     167  47 115.0  70.0      Yes 22.71       No     No     100
    ## 3424   Male     270  46 116.5  87.0      Yes 28.53       No     No      71
    ## 3425 Female     159  41  99.0  62.0      Yes 19.09       No     No      67
    ## 3426   Male     260  61 115.0  79.0      Yes 23.65       No     No      78
    ## 3427 Female     218  42 122.0  75.0       No 18.11       No     No      85
    ## 3428   Male     232  69 151.0  74.0      Yes 24.14       No     No      62
    ## 3429 Female     307  53 142.0  94.0       No 25.70       No     No      60
    ## 3430 Female     241  66 112.0  66.0      Yes 23.36       No     No      74
    ## 3431   Male     229  45 117.0  78.0      Yes 22.79       No     No      86
    ## 3433   Male     200  38 123.0  85.0      Yes 25.63       No     No      84
    ## 3434   Male     200  58 112.0  74.5      Yes 23.37       No     No      62
    ## 3435 Female     194  39 112.5  77.5      Yes 21.51       No     No      84
    ## 3437 Female     290  56 164.5 102.0       No 30.33       No     No     105
    ## 3438 Female     175  44 130.0  80.0      Yes 19.18       No     No     117
    ## 3439 Female     239  59 127.0  88.5       No 27.20       No     No      78
    ## 3440 Female     242  63 142.5  85.0       No 28.25       No     No      73
    ## 3444 Female     205  41 120.0  80.0      Yes 20.67       No     No      64
    ## 3445   Male     283  55 137.0  82.0      Yes 28.49       No     No      85
    ## 3446 Female     193  46 124.5  72.0       No 26.84       No     No      80
    ## 3447 Female     250  46 112.5  60.0      Yes 22.72       No     No      74
    ## 3449 Female     215  48 128.0  86.5       No 22.72       No     No      85
    ## 3450 Female     267  46 119.0  65.0      Yes 29.15       No     No      75
    ## 3451 Female     253  48 120.0  77.5      Yes 24.53       No     No      98
    ## 3452   Male     213  57 136.5  87.0       No 25.51      Yes     No     119
    ## 3453 Female     165  40 101.0  59.0      Yes 23.06       No     No      76
    ## 3454   Male     262  61 129.0  85.0      Yes 21.77       No     No      74
    ## 3455   Male     259  47 139.0  79.0       No 29.34       No     No      71
    ## 3456 Female     229  47 127.0  76.5      Yes 23.48       No     No      64
    ## 3457 Female     200  42 119.0  75.0      Yes 22.91       No     No      69
    ## 3458   Male     193  56 114.0  80.0       No 28.41       No     No      88
    ## 3459 Female     225  41 110.0  60.0       No 25.54       No     No      58
    ## 3461   Male     235  58 127.5  76.0      Yes 21.02       No     No     135
    ## 3463   Male     177  65 119.0  82.5      Yes 21.18       No     No      88
    ## 3464 Female     259  61 134.5  87.0      Yes 22.91       No     No      91
    ## 3466 Female     220  47 124.0  75.5       No 24.71       No     No      68
    ## 3467 Female     243  41 159.0 100.0      Yes 27.78       No     No      71
    ## 3469   Male     235  60 108.5  73.5      Yes 21.76       No     No     102
    ## 3470 Female     254  37 119.0  62.5      Yes 28.78       No     No      69
    ## 3471   Male     215  36 110.0  67.0      Yes 23.10       No     No      84
    ## 3472 Female     264  52 181.0 112.5       No 24.80       No    Yes      77
    ## 3473 Female     265  58 143.5  85.0       No 21.68      Yes     No     107
    ## 3474   Male     225  57 177.5 120.0      Yes 25.48       No     No      96
    ## 3475 Female     198  45 119.0  80.0      Yes 22.18       No     No      79
    ## 3476 Female     226  58 125.0  75.0      Yes 24.00       No     No      73
    ## 3477   Male     243  62 157.0  96.0      Yes 28.83       No     No      71
    ## 3478 Female     235  61 102.0  74.5       No 19.22       No     No      83
    ## 3479   Male     241  44 111.0  78.0      Yes 27.78      Yes     No     167
    ## 3480   Male     219  65 128.0  77.5      Yes 19.76       No     No      73
    ## 3481   Male     175  48 131.0  52.0       No 24.95       No     No      77
    ## 3482   Male     194  65 122.0  68.0      Yes 26.08       No     No      73
    ## 3484 Female     186  39 124.0  80.0       No 29.70       No     No     102
    ## 3486 Female     173  40 125.0  75.0      Yes 25.67       No     No     102
    ## 3487   Male     228  63 141.0  82.0      Yes 25.84       No     No      81
    ## 3488 Female     200  39 141.0  87.0       No 24.04       No     No      70
    ## 3489   Male     185  62 111.0  67.0      Yes 21.84       No     No      77
    ## 3490 Female     240  50 107.5  70.0       No 26.76       No     No     102
    ## 3491   Male     253  46 159.0 100.0      Yes 27.31       No     No      70
    ## 3492 Female     198  47 143.0  87.0      Yes 20.86       No     No      79
    ## 3494 Female     285  64 160.0  87.0       No 29.97       No    Yes      87
    ## 3495 Female     254  59 154.5  93.0      Yes 21.82       No     No      89
    ## 3496 Female     297  63 133.5  92.0       No 25.09       No     No      74
    ## 3497 Female     326  51 112.0  83.0      Yes 20.82       No     No      70
    ## 3498   Male     190  47 147.5  92.5      Yes 31.31       No     No      82
    ## 3499   Male     150  38 123.0  76.0      Yes 20.39       No     No      82
    ## 3500 Female     288  45 135.0  88.0       No 25.23       No     No      79
    ## 3501   Male     195  61 113.0  77.0       No 25.51       No     No      96
    ## 3502   Male     290  51 168.0 103.0      Yes 29.11       No     No      64
    ## 3503   Male     257  45 117.5  82.5      Yes 24.62       No     No      86
    ## 3504 Female     197  54 138.0  80.0       No 31.82       No     No      73
    ## 3505 Female     259  64 195.0 110.0      Yes 20.09       No     No      63
    ## 3506 Female     232  49 125.0  83.0       No 22.17       No     No      77
    ## 3507 Female     223  59 159.0 102.0       No 28.30       No    Yes      83
    ## 3508   Male     240  48 126.0  84.0      Yes 25.82       No     No      60
    ## 3509   Male     306  48 140.0  94.0       No 23.95       No     No      70
    ## 3510   Male     291  63 134.5  80.0      Yes 22.49       No     No      70
    ## 3511   Male     275  55 115.0  78.0       No 28.52       No     No      79
    ## 3512 Female     231  42 110.0  80.0      Yes 19.12       No     No      70
    ## 3513   Male     274  59 144.5  91.5       No 26.08       No     No      63
    ## 3514 Female     272  64 149.5  86.0       No 26.91       No     No      62
    ## 3515 Female     261  51 142.0  90.0       No 28.58       No     No      84
    ## 3516   Male     250  54 114.0  82.0      Yes 24.39       No     No      62
    ## 3517 Female     312  66 150.0  74.0       No 25.59       No     No      90
    ## 3519   Male     233  36 128.0  86.5       No 35.31       No     No      76
    ## 3520 Female     320  54 165.0  84.0       No 38.31       No     No      99
    ## 3521 Female     231  47 102.5  66.0       No 23.40       No     No      78
    ## 3522 Female     181  41 112.0  75.5       No 25.14       No     No      95
    ## 3523   Male     260  62 127.5  75.0      Yes 25.35       No     No     103
    ## 3524 Female     274  65 143.5  93.0       No 26.82       No     No      63
    ## 3525 Female     234  50 143.0  90.0       No 29.36       No     No      86
    ## 3526   Male     215  55 121.5  81.5      Yes 28.45       No     No      98
    ## 3528   Male     210  37 120.0  72.0      Yes 23.80       No     No      97
    ## 3530   Male     253  60 141.0  92.0      Yes 24.10       No     No      76
    ## 3531 Female     205  47 133.0  93.0       No 27.82       No     No      60
    ## 3532   Male     184  38 102.0  73.0       No 25.69       No     No      78
    ## 3533 Female     221  47 116.5  81.0      Yes 25.85       No     No      75
    ## 3535   Male     227  51 158.0 105.0       No 27.22       No     No      96
    ## 3536 Female     254  61 168.0  92.0       No 31.63       No     No      90
    ## 3538   Male     245  41 139.0  84.0       No 28.76       No     No      68
    ## 3539   Male     206  47 125.0  72.0      Yes 21.14       No     No      74
    ## 3541   Male     249  37 112.0  70.0      Yes 22.79       No     No      76
    ## 3542   Male     163  59 138.0  80.0      Yes 31.08       No     No      70
    ## 3543   Male     265  50 110.0  65.0      Yes 24.45       No     No      70
    ## 3544 Female     178  44 113.0  78.0       No 31.93       No     No      74
    ## 3545 Female     217  61 182.0  86.0       No 26.98       No     No     113
    ## 3546 Female     214  43 121.0  84.0      Yes 24.68       No     No      74
    ## 3547   Male     157  62 130.0  87.0       No 28.73       No     No      47
    ## 3548 Female     382  57 140.0  94.0      Yes 21.20       No     No      70
    ## 3549   Male     228  39 122.5  87.0      Yes 31.60       No     No      73
    ## 3550   Male     179  46 111.0  80.0      Yes 20.87       No     No      76
    ## 3551   Male     270  49 160.5 106.5      Yes 30.33       No     No      65
    ## 3552 Female     239  44 103.0  67.0      Yes 26.58       No     No      73
    ## 3553 Female     288  62 118.5  71.0       No 26.18       No     No      87
    ## 3554   Male     281  55 165.0 108.0      Yes 24.14       No     No      66
    ## 3555   Male     229  58 140.0  89.5      Yes 25.96       No     No      83
    ## 3556   Male     228  53 132.5  55.0      Yes 19.97       No    Yes      83
    ## 3557   Male     250  46 123.0  76.0      Yes 21.66       No     No      78
    ## 3558 Female     258  62 162.0  97.5       No 30.53       No     No      87
    ## 3559   Male     225  52 132.0  88.0       No 23.35       No     No      69
    ## 3560   Male     212  41 112.0  63.5      Yes 25.20       No     No      76
    ## 3561   Male     158  44 150.5  87.0      Yes 21.44       No     No      98
    ## 3562   Male     232  66 157.0  91.0      Yes 31.20       No     No      67
    ## 3564 Female     241  58 153.0 106.0       No 26.94       No     No      84
    ## 3565 Female     240  51 112.0  83.0      Yes 24.10       No     No      77
    ## 3566   Male     188  63 113.0  82.0      Yes 29.01       No     No      93
    ## 3567 Female     212  65 192.5 110.0      Yes 23.48       No     No      71
    ## 3568 Female     204  42 108.0  70.5       No 27.71       No     No      65
    ## 3569   Male     229  44 142.0  92.0       No 25.21       No     No      76
    ## 3570 Female     286  57 153.0  87.0       No 27.64       No     No      83
    ## 3571 Female     346  49 130.0  80.0      Yes 22.54       No     No      77
    ## 3572 Female     219  53 108.0  65.0       No 22.19       No     No      76
    ## 3573   Male     250  54 123.0  75.0       No 25.91       No     No      71
    ## 3574 Female     361  63 167.0 100.0       No 27.31       No     No     103
    ## 3575 Female     179  65 228.0 130.0      Yes 19.74       No    Yes      60
    ## 3576   Male     197  63 113.0  70.0      Yes 23.72       No     No      67
    ## 3577 Female     238  42 125.0  77.5       No 22.90       No     No      77
    ## 3578   Male     232  41 129.0  83.0      Yes 28.08       No     No      68
    ## 3579 Female     260  66 167.0  96.0       No 29.04       No     No      82
    ## 3580   Male     220  52 151.0 102.0       No 25.82       No     No      85
    ## 3581 Female     281  48 127.0  79.5       No 22.83       No     No      80
    ## 3582 Female     256  53 128.0  90.0      Yes 23.65       No     No     102
    ## 3583 Female     333  45 122.0  80.0       No 24.39       No     No      80
    ## 3584   Male     208  40 119.0  66.0      Yes 28.09       No     No      66
    ## 3585   Male     215  39 145.5  92.5      Yes 28.35       No     No      94
    ## 3586 Female     230  51 168.5  97.0       No 26.36       No     No      77
    ## 3587   Male     232  45 165.0  86.0      Yes 23.86       No     No      68
    ## 3588   Male     210  44 133.0  85.5       No 25.14       No     No      90
    ## 3589 Female     213  55 163.0  91.0      Yes 28.66       No     No      66
    ## 3590 Female     266  56 114.0  72.0      Yes 22.64       No     No      83
    ## 3591   Male     255  51 115.0  67.0      Yes 26.97       No     No      58
    ## 3592   Male     240  42 132.0  89.5       No 29.35       No     No     103
    ## 3593   Male     285  52 135.0  86.0       No 27.78       No     No      93
    ## 3594   Male     203  36 101.5  67.0      Yes 24.43       No     No      74
    ## 3595   Male     236  42 124.0  80.0       No 21.50       No     No      60
    ## 3596   Male     170  48 122.0  70.0      Yes 23.62       No     No      73
    ## 3597 Female     271  50 112.5  60.0      Yes 23.29       No     No      61
    ## 3598 Female     180  36 118.0  80.0      Yes 29.59       No     No      84
    ## 3599   Male     263  62 112.0  61.0      Yes 24.46       No     No      95
    ## 3600 Female     194  61 148.0  89.0       No 23.48       No     No     101
    ## 3601 Female     222  36 147.0  94.0      Yes 26.79       No     No      71
    ## 3602 Female     271  58 146.0  92.0      Yes 23.07       No     No      83
    ## 3603   Male     310  57 147.5  90.0       No 32.09       No     No      73
    ## 3604 Female     221  49 136.0  90.0       No 28.30       No     No      80
    ## 3605 Female     241  48 129.0  86.0       No 20.41       No     No      67
    ## 3606 Female     241  56 174.0  97.0      Yes 29.22      Yes    Yes     135
    ## 3607 Female     231  61 128.0  87.0      Yes 26.30       No     No      93
    ## 3608 Female     242  51 108.0  77.0       No 22.91       No     No      80
    ## 3609 Female     195  42 126.0  81.0      Yes 22.26       No     No      77
    ## 3610 Female     206  58 159.5  93.5      Yes 18.53       No     No      58
    ## 3611   Male     254  47 138.0  96.0       No 29.73       No     No      69
    ## 3614 Female     216  52 125.0  72.0       No 24.98       No     No      95
    ## 3615   Male     245  57 132.0  77.0      Yes 23.01      Yes     No     207
    ## 3616 Female     219  61 120.0  72.5       No 22.35       No     No      92
    ## 3617 Female     312  62 156.0 105.0       No 22.35       No     No      82
    ## 3618   Male     214  56 142.0  66.0      Yes 24.47       No     No      94
    ## 3619   Male     178  57 123.0  79.0       No 26.36       No     No      78
    ## 3620   Male     248  60 135.0  85.0      Yes 23.06       No     No      61
    ## 3621   Male     259  50 108.0  81.0      Yes 22.81       No     No      72
    ## 3623   Male     248  52 128.0  83.5       No 25.88       No     No      75
    ## 3624 Female     230  54 135.0  85.0       No 19.18       No     No      89
    ## 3625 Female     253  44 118.0  68.0      Yes 22.72       No     No     110
    ## 3626 Female     249  59 138.0  72.0       No 25.02       No     No      70
    ## 3628 Female     180  62 163.0  89.0       No 20.21       No     No      83
    ## 3629 Female     230  45 128.5  87.5       No 23.75       No     No      62
    ## 3630   Male     258  40 112.0  78.0      Yes 28.57       No     No      70
    ## 3631   Male     245  43 144.5  95.0       No 27.15       No     No      45
    ## 3632   Male     453  42 158.0 108.0      Yes 28.89       No     No     110
    ## 3633 Female     156  40 110.0  74.0       No 20.79       No     No      76
    ## 3634   Male     165  48 115.0  80.0       No 26.79       No     No      78
    ## 3635 Female     232  43 127.0  79.0       No 30.79       No     No      54
    ## 3636   Male     188  62 127.5  75.5      Yes 24.38       No     No      72
    ## 3637   Male     232  59 151.5 110.0      Yes 26.89       No     No      69
    ## 3639   Male     222  63 159.0  90.0       No 21.90       No     No      95
    ## 3640 Female     238  50 158.0  74.0      Yes 35.68       No     No      98
    ## 3641   Male     256  38 123.0  92.0       No 25.42       No     No      82
    ## 3642 Female     194  35 100.5  69.0      Yes 17.92       No     No      73
    ## 3643   Male     254  56 166.0 107.0       No 21.97       No     No      83
    ## 3644   Male     240  68 161.0 105.0      Yes 24.29       No    Yes      92
    ## 3646 Female     185  40 111.0  79.0       No 22.90       No     No      81
    ## 3647 Female     168  35  83.5  55.0      Yes 16.71       No     No      63
    ## 3649 Female     298  62 248.0 130.0       No 37.10       No    Yes      77
    ## 3650 Female     233  61 126.0  73.0       No 23.16       No     No      62
    ## 3651   Male     173  54 121.0  79.0      Yes 26.21       No     No      68
    ## 3652   Male     281  38 111.0  72.5      Yes 27.22       No     No      80
    ## 3653 Female     273  63 152.0  70.0       No 19.69       No     No      79
    ## 3654 Female     259  60 155.0  90.0       No 27.94       No     No      95
    ## 3655 Female     185  49 108.0  70.0      Yes 20.13       No     No      58
    ## 3656   Male     325  53 172.5 112.5       No 28.47       No     No      74
    ## 3657   Male     204  50 122.0  81.0       No 32.22       No     No      73
    ## 3658   Male     237  62 163.0  94.0      Yes 25.62       No     No      84
    ## 3659 Female     235  55 123.0  81.0       No 31.44       No     No      78
    ## 3660 Female     221  46 108.0  73.0       No 20.06       No     No      85
    ## 3661 Female     242  40 112.5  62.5      Yes 27.65       No     No      70
    ## 3662   Male     183  38 107.5  71.0      Yes 23.74       No     No      74
    ## 3663 Female     256  55 143.0  82.5       No 23.81       No     No      90
    ## 3664 Female     231  40 129.0  87.0      Yes 23.29       No     No      99
    ## 3665   Male     149  58  98.0  60.0       No 24.73       No     No      71
    ## 3666 Female     210  49 116.0  76.0       No 28.86       No     No      71
    ## 3667   Male     214  50 114.0  72.0      Yes 22.93       No     No      83
    ## 3668   Male     252  45 124.0  89.0       No 22.82       No     No      71
    ## 3669   Male     306  42 196.0 109.0      Yes 27.72       No     No      87
    ## 3670   Male     251  51 160.0  98.0      Yes 24.63       No     No      85
    ## 3671 Female     267  56 122.5  85.0       No 24.22       No     No     100
    ## 3672 Female     294  53 156.0  95.0       No 26.05       No     No     115
    ## 3673 Female     240  38 126.5  75.5       No 24.38       No     No      64
    ## 3675   Male     205  48 108.0  75.0      Yes 17.50       No     No      70
    ## 3676   Male     211  40 122.0  81.0      Yes 30.55       No     No      91
    ## 3677 Female     186  62 176.5  92.0       No 22.53       No    Yes      60
    ## 3678   Male     197  36 115.0  65.0       No 20.42       No     No      77
    ## 3679   Male     218  43 121.0  69.0      Yes 24.21       No     No     103
    ## 3680 Female     246  56 131.0  79.0       No 27.69       No     No      65
    ## 3681   Male     190  59 127.0  77.0       No 28.47       No     No     100
    ## 3682 Female     258  43 161.5  96.0      Yes 38.96       No     No      84
    ## 3683   Male     285  43 100.5  66.0      Yes 22.05       No     No      75
    ## 3684 Female     195  50 131.5  83.0      Yes 24.61       No     No      78
    ## 3685 Female     217  46 131.5  77.0       No 20.25       No     No      76
    ## 3687 Female     264  57 135.5  77.0       No 24.50       No     No      93
    ## 3688   Male     175  43 125.0  76.0      Yes 24.92       No     No      95
    ## 3689   Male     217  61 131.0  83.0      Yes 26.46       No     No      81
    ## 3690   Male     224  60 143.5  77.5       No 26.13       No     No      81
    ## 3691 Female     309  48 136.0  90.0       No 26.83       No     No      75
    ## 3692 Female     201  38 123.5  78.0      Yes 27.14       No     No      77
    ## 3693   Male     352  44 164.0 119.0      Yes 28.92       No     No      72
    ## 3694   Male     226  64 145.0  87.0       No 27.02       No     No      77
    ## 3696   Male     223  52 133.0  82.0      Yes 21.18       No     No      77
    ## 3697 Female     260  60 139.0  81.0      Yes 24.68       No     No      70
    ## 3698   Male     185  62 124.0  84.0      Yes 24.53       No     No      70
    ## 3699   Male     320  58 139.0  81.5       No 23.65       No     No      82
    ## 3700 Female     246  48 113.0  87.0      Yes 18.01       No     No      63
    ## 3701 Female     302  44 116.0  77.0      Yes 22.67       No     No      98
    ## 3702   Male     237  47 182.0 110.0      Yes 28.88       No    Yes      83
    ## 3703 Female     216  44 113.0  74.0       No 22.19       No     No      76
    ## 3704 Female     275  46 126.0  71.0      Yes 24.91       No     No      71
    ## 3705 Female     192  37 112.0  67.0       No 24.61       No     No      58
    ## 3706   Male     242  57 130.0  74.0       No 28.90       No     No      53
    ## 3707   Male     202  55 104.0  76.0       No 31.68       No     No      74
    ## 3708 Female     241  54 106.0  77.0       No 27.64       No     No      74
    ## 3709 Female     309  55 177.5 110.0       No 22.89       No     No      73
    ## 3710 Female     262  46 121.0  78.0      Yes 24.24       No     No      72
    ## 3711 Female     197  39 134.0  78.0       No 30.36       No     No      83
    ## 3712 Female     209  37 110.0  77.5       No 17.93       No     No      80
    ## 3713 Female     235  43 128.5  80.0      Yes 18.83       No     No      70
    ## 3714   Male     217  46 115.0  72.0      Yes 21.34      Yes     No      73
    ## 3715 Female     241  68 154.0  96.0       No 30.12       No     No      70
    ## 3716 Female     280  46 202.0 124.0      Yes 28.06       No    Yes      63
    ## 3717   Male     260  41 151.0  85.0      Yes 33.08       No     No      91
    ## 3718 Female     265  55 154.0  87.0      Yes 20.92       No     No      66
    ## 3719   Male     186  39 162.0 109.0      Yes 29.72       No     No     129
    ## 3721 Female     210  46 117.0  78.5      Yes 22.54       No     No     115
    ## 3723 Female     222  53 123.0  82.0      Yes 25.52       No     No      67
    ## 3724 Female     292  58 110.0  83.0       No 26.98       No     No      73
    ## 3725   Male     221  57 117.0  72.0      Yes 25.51       No     No      83
    ## 3726   Male     161  40 122.0  85.0       No 30.80       No     No      85
    ## 3727 Female     273  49 154.0  80.0       No 20.26       No     No      63
    ## 3728 Female     291  43 106.0  65.0      Yes 23.83       No     No      82
    ## 3729   Male     212  37 114.5  77.0       No 31.22       No     No      68
    ## 3730 Female     276  43  99.0  62.0      Yes 22.17       No     No      80
    ## 3731 Female     165  46  99.5  66.0      Yes 21.67       No     No      66
    ## 3732 Female     276  62 185.0  95.0       No 26.21       No     No     110
    ## 3734 Female     310  61 108.0  70.0       No 30.23       No     No      65
    ## 3736 Female     184  49 121.0  82.0       No 21.14       No    Yes      78
    ## 3737   Male     186  56 116.0  67.0      Yes 24.62       No     No      83
    ## 3739 Female     207  54 158.0  89.0       No 31.44       No     No     112
    ## 3740 Female     190  39  85.0  70.0      Yes 22.43       No     No      60
    ## 3742 Female     217  65 169.0 111.0      Yes 32.54       No     No      78
    ## 3743 Female     279  59 116.0  66.0      Yes 21.83       No     No      79
    ## 3744 Female     209  40 130.0  84.5      Yes 39.94       No     No     104
    ## 3745 Female     199  41 117.0  86.0       No 34.54       No     No      84
    ## 3746 Female     202  36 105.5  67.0      Yes 22.66       No     No      63
    ## 3747   Male     235  48 135.0  88.0       No 27.61       No     No     137
    ## 3748 Female     193  52 146.0  89.0       No 25.37       No     No      84
    ## 3749   Male     207  63 165.0 100.0      Yes 21.33       No     No      77
    ## 3750   Male     206  61 143.0  96.0       No 27.04       No     No      87
    ## 3752   Male     246  53 176.0 109.0      Yes 32.32       No     No      59
    ## 3753 Female     140  44 118.0  74.0      Yes 26.51       No     No      82
    ## 3756   Male     246  44 136.0  94.0      Yes 24.56       No     No      86
    ## 3757 Female     244  48 117.0  81.5       No 28.96       No     No      78
    ## 3758 Female     268  45 130.0  94.0      Yes 34.27       No     No      93
    ## 3759 Female     252  58 135.0  84.0      Yes 28.24       No     No      79
    ## 3760   Male     200  40 122.5  75.0      Yes 20.25       No     No      67
    ## 3762   Male     260  65 135.0  97.0      Yes 25.55       No     No      81
    ## 3763 Female     175  47 107.0  69.0      Yes 23.64       No     No      70
    ## 3764 Female     270  44 167.5  92.5      Yes 21.28       No     No      77
    ## 3765 Female     226  47 122.5  80.0      Yes 24.62       No     No      68
    ## 3766 Female     238  57 133.0  72.0      Yes 18.09       No     No     115
    ## 3767 Female     275  60 138.0  87.0       No 29.64       No     No      86
    ## 3768   Male     375  65 151.0  91.0       No 27.72       No     No      67
    ## 3770 Female     206  42 101.0  75.0       No 18.73       No     No      84
    ## 3771   Male     250  58 109.0  78.5       No 25.26       No     No      83
    ## 3772   Male     189  45 132.0  78.0       No 28.40      Yes     No     177
    ## 3774 Female     238  63 123.0  66.0       No 28.60       No     No      98
    ## 3775   Male     167  49 119.0  67.0      Yes 25.83       No     No      83
    ## 3776 Female     220  51 137.0  79.0      Yes 21.66       No     No      74
    ## 3777   Male     193  40 130.0  89.0       No 28.32       No     No      84
    ## 3778 Female     308  49 128.0  78.0       No 24.82       No     No      70
    ## 3780 Female     212  55 161.0  85.0       No 28.26       No     No      63
    ## 3781   Male     270  38 120.0  75.0      Yes 23.76       No     No      96
    ## 3782 Female     262  54 230.0 110.0      Yes 24.76       No     No      97
    ## 3783 Female     180  40 118.0  86.5      Yes 22.69       No     No      67
    ## 3785   Male     235  54 143.0  75.0      Yes 32.99       No     No     102
    ## 3786   Male     189  66 140.0  71.0       No 27.56      Yes     No     119
    ## 3787   Male     178  63 155.0  79.0       No 25.90       No     No      61
    ## 3788   Male     170  52 153.0 102.5       No 29.35       No     No      67
    ## 3789 Female     202  48 146.5  79.0       No 22.19       No     No      95
    ## 3791   Male     208  54 137.5  82.5      Yes 25.58       No     No      63
    ## 3792 Female     182  38 138.0  72.0       No 21.67       No     No     108
    ## 3793   Male     215  65 147.5  95.0      Yes 29.08       No     No      88
    ## 3794   Male     196  53 129.0  80.0       No 26.32       No     No      89
    ## 3795 Female     244  44 101.0  66.0      Yes 25.38       No     No      76
    ## 3797 Female     190  47 162.0  85.0       No 30.59       No     No      80
    ## 3798 Female     285  57 197.0  72.0       No 23.41       No     No      78
    ## 3799 Female     231  44 133.0  89.0      Yes 29.29       No     No      83
    ## 3800 Female     286  43 164.0  89.0      Yes 24.44       No     No      87
    ## 3801   Male     202  50 189.0 121.0      Yes 33.81       No     No      72
    ## 3802   Male     214  57 132.0  74.0       No 27.98       No     No      87
    ## 3803 Female     226  41 125.5  82.0      Yes 23.80       No     No      75
    ## 3804 Female     297  45 134.0  93.0       No 28.81       No     No      74
    ## 3805 Female     193  46 118.0  92.0      Yes 21.14       No     No      78
    ## 3806 Female     270  50 142.5  85.0       No 21.86       No     No      75
    ## 3807 Female     216  42 124.0  82.0       No 28.74       No     No      67
    ## 3808 Female     240  58 126.0  52.0       No 25.66       No     No      63
    ## 3809   Male     240  48 136.0  95.5       No 26.36       No     No      73
    ## 3810 Female     188  36 112.0  78.0       No 22.54       No     No      73
    ## 3812   Male     157  68 106.0  48.0      Yes 26.73       No     No      65
    ## 3813 Female     275  37 118.0  71.0      Yes 23.10       No     No      95
    ## 3815 Female     204  50 147.0 100.0      Yes 39.94       No     No      90
    ## 3817   Male     211  62 128.0  78.0       No 27.99       No     No      83
    ## 3818 Female     178  40 142.0  84.0      Yes 34.46       No     No      77
    ## 3819 Female     268  41 140.0  92.5      Yes 24.71       No     No      90
    ## 3820 Female     164  38 113.0  68.0      Yes 25.75       No     No      75
    ## 3821   Male     270  53 151.0  89.0      Yes 26.76       No     No      75
    ## 3823 Female     207  46 144.0  88.0      Yes 23.65       No     No      86
    ## 3824 Female     190  51 153.0 102.5       No 39.22       No     No      69
    ## 3825 Female     226  39 146.0  86.0       No 24.41       No     No      85
    ## 3826 Female     257  37 141.0  93.0      Yes 41.29       No     No      58
    ## 3829 Female     200  46 112.5  71.0      Yes 18.68       No     No      77
    ## 3830   Male     199  50 113.0  79.0      Yes 24.11       No     No      64
    ## 3831 Female     227  48 100.0  76.0       No 29.45       No     No      67
    ## 3832   Male     273  48 140.0  90.0       No 27.48       No     No      78
    ## 3833   Male     220  49 114.0  82.0       No 24.68       No     No      60
    ## 3834   Male     351  58 130.0  77.5       No 31.84       No     No      80
    ## 3835   Male     226  60 155.0  92.5       No 30.85       No     No      87
    ## 3836   Male     155  49 128.0  82.0      Yes 23.58       No     No      77
    ## 3838 Female     243  41  97.0  63.0      Yes 22.53       No     No      64
    ## 3839 Female     175  38 112.0  73.0      Yes 19.49       No     No      71
    ## 3841   Male     208  56 167.0  92.0       No 24.66       No    Yes      75
    ## 3842 Female     240  54 113.0  73.0      Yes 24.21       No     No      77
    ## 3844 Female     203  69 166.0  90.0       No 25.40       No    Yes      80
    ## 3845   Male     223  67 214.0  94.0       No 25.86       No    Yes      87
    ## 3846   Male     192  43 143.0  88.0      Yes 27.94       No     No      79
    ## 3847   Male     194  48 134.5  90.5      Yes 25.68       No     No      72
    ## 3848   Male     263  48 132.0  91.0      Yes 40.08       No     No      91
    ## 3849 Female     230  42 124.0  80.0      Yes 24.87       No     No      77
    ## 3850 Female     211  47 159.5  82.5      Yes 34.08      Yes     No     250
    ## 3852 Female     282  62 125.0  75.0       No 29.88      Yes     No     136
    ## 3853   Male     260  57 123.0  73.0      Yes 27.51       No     No      83
    ## 3854   Male     163  38 117.5  75.0      Yes 28.30       No     No      70
    ## 3856 Female     315  51 119.0  75.0      Yes 25.79       No     No      55
    ## 3857   Male     186  40 131.0  81.0      Yes 22.14       No     No      87
    ## 3858 Female     250  43 123.0  74.0      Yes 26.01       No     No      90
    ## 3859   Male     155  34 117.5  72.5      Yes 23.51       No     No      65
    ## 3860 Female     185  55 140.0  84.0       No 25.94       No     No      90
    ## 3861 Female     205  44 116.0  70.0       No 21.99       No     No      85
    ## 3862   Male     252  50 114.0  75.0       No 30.89       No     No      69
    ## 3863   Male     225  52 119.0  65.0       No 26.89       No     No      74
    ## 3864 Female     235  63 125.0  79.0       No 24.38       No     No      83
    ## 3865   Male     240  43 147.5  88.0      Yes 25.60       No     No     113
    ## 3866   Male     242  41 124.5  86.5      Yes 28.80       No     No      67
    ## 3867 Female     180  34 111.0  56.0      Yes 21.51       No     No      78
    ## 3868   Male     195  37 141.0  84.0       No 25.66       No     No     117
    ## 3869 Female     210  50 105.0  77.0       No 23.96       No     No      86
    ## 3870 Female     212  36 127.0  83.0      Yes 26.82       No     No      75
    ## 3871 Female     254  59 126.5  79.0      Yes 25.92       No     No      77
    ## 3872   Male     248  40 136.5  78.0      Yes 24.60       No     No      99
    ## 3873   Male     283  42 145.0  95.0      Yes 27.53       No     No      84
    ## 3874   Male     186  42 134.0  86.5       No 25.71       No     No      92
    ## 3875   Male     212  43 135.0  86.0      Yes 30.22       No     No      75
    ## 3876   Male     309  49 145.0  92.0      Yes 32.13       No     No      73
    ## 3877 Female     308  47 138.0  94.0      Yes 24.33       No     No      90
    ## 3878   Male     225  60 149.0  96.0       No 27.73       No     No      60
    ## 3879 Female     275  46 170.0 118.0      Yes 36.12       No     No      84
    ## 3880   Male     154  46 141.0  90.0      Yes 22.76       No     No      65
    ## 3881   Male     246  49 141.0  92.0      Yes 27.92       No     No      76
    ## 3884   Male     198  53 142.5  82.0       No 23.84       No     No      78
    ## 3885 Female     270  57 120.0  79.0      Yes 24.83       No     No      81
    ## 3886 Female     208  44 142.0  88.0       No 31.29       No    Yes      77
    ## 3887 Female     286  60 172.5  85.0       No 22.00       No     No      71
    ## 3888 Female     236  47 128.0  81.0      Yes 27.42       No     No      93
    ## 3889 Female     214  40 109.5  69.0      Yes 20.32       No     No      81
    ## 3890   Male     175  44 104.0  78.0       No 26.26       No     No      82
    ## 3891 Female     271  65 144.5  88.0      Yes 32.41      Yes     No     116
    ## 3892   Male     208  49 118.0  73.0      Yes 24.16       No     No      75
    ## 3893 Female     265  56 150.0  84.0      Yes 28.66       No     No      73
    ## 3894   Male     240  41 158.0 102.0       No 23.68       No     No      81
    ## 3895 Female     212  46 148.0  91.5       No 26.22       No     No      75
    ## 3896 Female     260  56 120.0  84.0      Yes 36.18       No     No      76
    ## 3897 Female     330  59 155.0  90.0       No 23.94       No     No      96
    ## 3898   Male     185  45 108.0  70.0       No 20.50       No     No      87
    ## 3899   Male     161  63 196.0 102.0       No 28.43       No     No      88
    ## 3900 Female     194  45 115.0  85.0      Yes 27.77       No     No      82
    ## 3901   Male     259  59 116.0  83.0      Yes 27.93       No     No      83
    ## 3902   Male     220  60 167.5 110.0       No 30.41       No     No      84
    ## 3903 Female     340  48 143.0  93.0      Yes 23.08       No     No      83
    ## 3904 Female     165  50  96.5  62.5      Yes 23.48       No     No      78
    ## 3905 Female     220  45 133.5  85.5      Yes 25.38       No     No      73
    ## 3906 Female     257  55 132.0  81.0       No 27.49       No     No      86
    ## 3907   Male     237  55 124.0  84.0       No 28.18       No     No      97
    ## 3909 Female     254  57 146.5  81.0      Yes 41.61       No     No      85
    ## 3910 Female     193  65 153.0 100.0       No 32.00      Yes     No     107
    ## 3911   Male     230  36 122.0  78.0       No 26.53       No     No      78
    ## 3912   Male     203  60 140.0  95.0       No 28.04       No     No      83
    ## 3913 Female     284  44 143.0  92.0      Yes 21.19       No     No      88
    ## 3915 Female     244  55 145.5  94.0      Yes 28.86       No     No      72
    ## 3916 Female     240  61 163.0 112.5       No 26.80       No     No      82
    ## 3918   Male     307  55 122.0  86.0      Yes 25.04       No     No      74
    ## 3919 Female     233  59 149.0  85.0       No 24.67       No     No      72
    ## 3920 Female     275  60 134.0  74.0       No 19.91       No     No      85
    ## 3921 Female     223  54 110.0  67.5       No 21.22      Yes     No     294
    ## 3922 Female     232  45 163.0 101.0      Yes 23.26       No     No      74
    ## 3924 Female     225  56 113.0  75.5       No 27.51       No     No     104
    ## 3925 Female     190  35 115.0  77.0       No 23.95       No     No      80
    ## 3926   Male     185  34 100.5  66.0      Yes 24.42       No     No     115
    ## 3928 Female     274  46 158.0  97.0      Yes 22.83       No     No      78
    ## 3929   Male     212  45 118.0  71.0      Yes 21.03       No     No      83
    ## 3930 Female     261  50 110.0  76.0      Yes 23.31       No     No      85
    ## 3931 Female     262  55 122.5  84.0       No 28.68       No     No      76
    ## 3932   Male     231  53 130.0  71.5      Yes 25.23       No     No      75
    ## 3933 Female     246  60 160.0  92.0      Yes 26.38       No     No      73
    ## 3934   Male     318  39 118.0  75.5      Yes 27.52       No     No      70
    ## 3937 Female     176  40  99.0  59.0      Yes 22.13       No     No      78
    ## 3938   Male     239  52 125.0  88.0      Yes 23.86       No     No     123
    ## 3939   Male     227  35 106.0  73.0      Yes 29.27       No     No      79
    ## 3940   Male     262  49 115.0  83.0      Yes 22.86       No     No      57
    ## 3941   Male     201  43 126.0  82.0      Yes 29.48       No     No      92
    ## 3942   Male     290  40 120.0  89.0      Yes 27.99       No     No      66
    ## 3944   Male     215  39 102.0  64.5      Yes 24.50       No     No      62
    ## 3946 Female     295  51 131.0  87.0       No 24.41       No     No      77
    ## 3948   Male     228  37 113.0  83.0       No 24.81       No     No      73
    ## 3949 Female     286  61 141.0  81.0      Yes 23.61       No     No      52
    ## 3950 Female     186  47 139.0  85.0      Yes 27.90      Yes     No     125
    ## 3951 Female     281  67 101.0  59.0       No 23.10       No     No      85
    ## 3953 Female     204  55 100.5  62.0       No 29.44       No     No      88
    ## 3954 Female     252  41 133.0  90.0       No 26.83       No     No      64
    ## 3955 Female     247  60 131.0  81.0       No 22.19       No     No      94
    ## 3956 Female     220  42 129.0  81.0      Yes 19.74       No     No      61
    ## 3957 Female     214  38 115.0  90.0       No 25.69       No     No      65
    ## 3958   Male     277  47 137.0  86.0      Yes 26.25       No     No      85
    ## 3960   Male     238  61 119.0  78.0       No 25.36       No     No      75
    ## 3961   Male     232  52 154.0  92.0       No 26.51       No     No      74
    ## 3962 Female     246  56 136.0  87.0       No 26.21       No     No      72
    ## 3964 Female     227  34 102.0  68.0       No 24.96       No     No      75
    ## 3965   Male     241  58 164.0  97.0       No 32.18       No    Yes      54
    ## 3966 Female     325  54 170.0 107.0      Yes 25.07       No    Yes      64
    ## 3967 Female     224  48 152.5  90.0       No 29.80       No     No      85
    ## 3968 Female     254  59 116.0  71.0       No 25.48       No     No      98
    ## 3969   Male     234  53 113.0  68.0       No 24.80      Yes     No     108
    ## 3970   Male     210  43 127.5  82.5       No 27.94       No     No      80
    ## 3971   Male     267  54 141.5  91.0      Yes 25.36       No     No      87
    ## 3972 Female     274  48 110.0  84.0       No 22.51       No     No      78
    ## 3973 Female     262  48 125.0  80.0       No 25.23       No     No      72
    ## 3974 Female     205  46 115.0  75.0      Yes 19.48       No     No      78
    ## 3975   Male     188  43 137.0  87.0      Yes 26.28       No     No      73
    ## 3976 Female     263  53 173.0  89.0       No 23.03       No     No      82
    ## 3977   Male     227  51 162.5 104.0      Yes 34.97       No     No      65
    ## 3978 Female     410  52 105.0  67.5      Yes 27.33       No     No      90
    ## 3979 Female     223  43 100.0  70.0       No 22.73       No     No      68
    ## 3980   Male     238  65 146.0  86.0      Yes 29.47       No     No      66
    ## 3981   Male     232  50 148.5  94.0       No 25.78      Yes     No      88
    ## 3982 Female     350  49 135.0  86.5      Yes 25.56       No     No      83
    ## 3983   Male     223  35 128.0  82.0      Yes 19.99       No     No      67
    ## 3984 Female     252  61 119.0  77.0      Yes 23.20       No     No      65
    ## 3985 Female     257  44 129.0  93.0       No 27.56       No     No      76
    ## 3986   Male     289  41 109.0  74.0      Yes 25.80       No     No      86
    ## 3987   Male     266  60 115.5  82.5       No 23.68       No     No      83
    ## 3989   Male     296  47 141.0  93.0      Yes 28.50      Yes     No     332
    ## 3990 Female     235  41 129.0  94.0      Yes 23.71       No     No      81
    ## 3991 Female     230  53 170.0 113.0       No 29.55       No     No     115
    ## 3992   Male     231  48 121.0  81.0       No 27.56       No     No      85
    ## 3993   Male     206  39 124.0  78.0      Yes 19.98       No     No      80
    ## 3994 Female     185  38 100.0  72.0       No 22.15       No     No      83
    ## 3995 Female     223  65 158.0  90.0       No 27.26       No     No      93
    ## 3996 Female     192  39 109.0  61.0      Yes 23.36       No     No      84
    ## 3997 Female     185  37  99.0  59.0       No 22.52       No     No      69
    ## 3998 Female     181  37 120.5  80.0       No 26.29       No     No      67
    ## 3999 Female     190  39 133.0  87.0       No 32.09       No     No      80
    ## 4000 Female     188  52 130.0  71.0      Yes 23.88       No     No      89
    ## 4001   Male     232  55 132.5  87.0       No 20.72       No     No      71
    ## 4002   Male     252  62 156.5  93.0       No 28.65       No     No      97
    ## 4003   Male     288  46 131.0  81.0      Yes 25.94       No     No      80
    ## 4004   Male     222  52 108.0  70.0      Yes 23.09       No     No      61
    ## 4005   Male     276  68 127.0  66.5       No 25.78       No     No     104
    ## 4006   Male     207  45 111.0  72.0      Yes 21.49       No     No      87
    ## 4007 Female     266  62 124.0  69.0       No 22.90       No     No      82
    ## 4008 Female     279  53 167.0 101.0       No 30.85       No     No      87
    ## 4009 Female     215  42 111.0  72.0       No 25.38       No     No      55
    ## 4010 Female     342  65 168.5  98.5       No 23.87       No    Yes      60
    ## 4011   Male     260  63 155.5  98.0       No 30.08      Yes     No     109
    ## 4012 Female     267  49 107.0  74.0      Yes 27.80       No     No      75
    ## 4013 Female     209  52 111.0  79.0       No 23.63       No     No      70
    ## 4014 Female     219  45 123.5  80.5       No 24.45       No     No      83
    ## 4015   Male     218  53 120.0  80.0      Yes 29.87       No     No      73
    ## 4016 Female     358  62 215.0 110.0      Yes 37.62      Yes    Yes     368
    ## 4017   Male     190  40 110.0  70.0      Yes 24.63       No     No      72
    ## 4018   Male     321  57 192.5 113.0       No 25.94       No    Yes      90
    ## 4019 Female     171  41 135.0  82.5      Yes 24.35       No     No      82
    ## 4020   Male     170  37 111.0  74.0      Yes 26.00       No     No      67
    ## 4021 Female     233  46 106.0  60.0      Yes 20.84      Yes     No     348
    ## 4022   Male     266  38 149.0  91.0      Yes 31.08       No     No      61
    ## 4023   Male     228  48 110.0  70.0      Yes 21.47       No     No      76
    ## 4024 Female     242  62 137.0  75.0       No 30.51       No     No      78
    ## 4025 Female     260  55 172.5 100.5       No 32.27       No     No      72
    ## 4026   Male     237  68 130.0  62.0       No 33.52       No     No      82
    ## 4027   Male     213  45 130.0  80.0       No 27.25       No     No      75
    ## 4028 Female     210  63 172.0 108.0       No 33.02       No     No      93
    ## 4029   Male     290  58 124.0  76.0      Yes 21.65       No     No      81
    ## 4030 Female     199  46 102.0  56.0       No 21.96       No     No      84
    ## 4031 Female     163  43 104.5  65.0      Yes 17.84       No     No      71
    ## 4032 Female     204  43 133.0  86.5      Yes 26.01       No     No      79
    ## 4033 Female     254  38 120.5  76.0       No 22.73       No     No      78
    ## 4034 Female     200  64 151.0  90.0       No 21.99       No     No      87
    ## 4035 Female     245  67 169.0  82.0      Yes 26.05       No    Yes     122
    ## 4036   Male     255  45 111.0  80.0      Yes 30.92       No     No      80
    ## 4038   Male     347  45 157.0  98.0       No 26.63       No     No      80
    ## 4039 Female     315  53 140.0  86.5       No 25.21       No    Yes      81
    ## 4040   Male     259  48 147.5  87.5       No 25.10       No     No      73
    ## 4041   Male     228  61 123.0  88.0      Yes 26.88       No     No      67
    ## 4043 Female     279  64 172.0  87.0       No 24.01       No     No      70
    ## 4044 Female     277  57 148.0  89.5       No 31.82       No     No      89
    ## 4045   Male     256  45 144.0  96.0      Yes 26.60       No     No      83
    ## 4048   Male     248  38 110.0  61.0      Yes 22.17       No     No      55
    ## 4049 Female     214  56 147.0  65.0      Yes 17.68       No     No      87
    ## 4050 Female     244  59 160.0  85.0       No 29.41       No     No      85
    ## 4051 Female     205  48 140.0  99.0       No 24.41       No     No      92
    ## 4053   Male     334  40 120.0  77.0      Yes 22.66       No     No      92
    ## 4054   Male     246  41 111.5  67.0       No 18.76       No     No      60
    ## 4055 Female     285  45 120.0  80.0       No 27.45       No     No      93
    ## 4056   Male     162  59 117.0  73.0      Yes 25.01       No     No      87
    ## 4058   Male     270  47 152.5 108.0       No 26.09       No     No      88
    ## 4059 Female     181  39 112.0  71.0       No 21.80       No     No      67
    ## 4060   Male     205  65 106.0  73.0       No 23.14       No     No      87
    ## 4061   Male     189  37 127.0  85.0      Yes 32.35       No     No      84
    ## 4062   Male     192  40 128.0  81.0       No 25.41       No     No      76
    ## 4063 Female     272  53 165.0  95.0       No 28.92       No     No     100
    ## 4064 Female     167  38 105.0  70.0      Yes 19.76       No     No      80
    ## 4065 Female     236  57 164.0 100.0       No 25.45       No     No      67
    ## 4066   Male     143  42 119.0  73.0      Yes 22.04       No     No      70
    ## 4067   Male     228  54 143.5  92.0       No 26.65       No     No      74
    ## 4069   Male     267  61 121.0  81.0       No 26.58       No     No      95
    ## 4070 Female     200  39 110.0  68.0      Yes 20.24       No     No      62
    ## 4071 Female     304  65 139.0  75.0       No 25.12      Yes     No     116
    ## 4072   Male     247  58 156.0  90.0      Yes 24.13       No     No      94
    ## 4073 Female     263  46 110.0  65.0      Yes 27.27       No     No      73
    ## 4074   Male     230  58 128.0  88.0       No 25.98       No     No      58
    ## 4075   Male     258  68 128.0  79.0       No 29.54       No     No      86
    ## 4076   Male     240  40 150.0  98.0       No 40.38       No     No      74
    ## 4077   Male     210  38 116.5  74.0       No 21.19       No     No      89
    ## 4079   Male     232  52 155.0  80.5      Yes 29.60       No     No      67
    ## 4080   Male     204  56 119.0  75.5      Yes 29.30       No     No      79
    ## 4081   Male     178  66 128.0  65.0      Yes 28.74       No     No      95
    ## 4083 Female     260  54 173.0  88.0       No 27.20       No     No      87
    ## 4084   Male     212  44 130.0  94.0       No 26.97       No     No      64
    ## 4085 Female     230  47 137.0  79.0       No 27.13       No     No      76
    ## 4086 Female     205  42 112.0  65.0       No 27.35       No     No      60
    ## 4087   Male     206  59 187.0  97.0      Yes 26.00       No     No      69
    ## 4088   Male     243  55 116.0  72.0      Yes 24.73       No     No      72
    ## 4089 Female     202  40 139.0  85.0       No 22.01       No     No      64
    ## 4090   Male     200  43 141.0  89.0       No 26.59       No     No     105
    ## 4091   Male     184  43 127.5  81.0      Yes 28.31       No     No      75
    ## 4093 Female     233  56 173.0  98.5       No 21.88       No     No      76
    ## 4094 Female     203  52 129.0  97.0       No 37.58       No     No      77
    ## 4095   Male     373  62 138.5  85.0      Yes 23.35       No     No      67
    ## 4096 Female     222  44 130.0  86.0      Yes 27.42       No     No      84
    ## 4097 Female     316  63  92.5  66.5       No 23.10       No     No      76
    ## 4098   Male     155  57 143.0  96.0       No 31.01       No     No      88
    ## 4099   Male     305  38 114.0  80.0      Yes 28.61       No     No      71
    ## 4100 Female     164  38 110.0  76.0      Yes 23.85       No     No      83
    ## 4101   Male     155  38 131.5  77.5       No 25.94       No     No      84
    ## 4102   Male     248  59 151.5  71.0       No 27.14       No     No     110
    ## 4103 Female     292  52 157.0 112.0      Yes 29.56       No     No      84
    ## 4105 Female     185  35 131.5  84.0      Yes 20.32       No     No      76
    ## 4106 Female     225  61 194.0 111.0       No 51.28       No    Yes     103
    ## 4107   Male     261  67 170.0 100.0      Yes 22.71       No     No      79
    ## 4108 Female     273  56 130.5  82.0      Yes 25.48       No     No      91
    ## 4109 Female     237  55 153.0  80.0       No 28.90       No     No      72
    ## 4110   Male     249  48 132.0  95.0      Yes 29.79       No     No      53
    ## 4111 Female     275  51 150.0  99.0      Yes 23.17       No     No      65
    ## 4112 Female     214  43 133.0  86.0      Yes 22.72       No     No      77
    ## 4113   Male     191  40 130.5  63.0      Yes 23.92       No     No      66
    ## 4114 Female     212  46 118.0  79.0       No 26.83       No     No      72
    ## 4115   Male     208  62 144.0  80.0       No 31.42       No     No      66
    ## 4116 Female     175  35 121.5  74.5       No 20.86       No     No      93
    ## 4117 Female     227  49 150.0  91.0      Yes 24.30       No     No      83
    ## 4118   Male     239  49 116.0  84.0       No 31.47       No     No      76
    ## 4119 Female     260  47 126.5  81.0      Yes 26.58       No     No      82
    ## 4120 Female     202  45  93.5  58.0      Yes 21.25       No     No      60
    ## 4121 Female     197  39 126.5  76.5      Yes 19.71       No     No      63
    ## 4123 Female     346  51 152.0  96.5       No 25.29       No     No      79
    ## 4124 Female     250  55 132.0  88.0       No 22.62       No     No      70
    ## 4125 Female     250  57 152.5  92.5       No 32.31       No     No      94
    ## 4127 Female     179  41 111.0  79.0       No 25.87       No     No      82
    ## 4128   Male     260  48 124.5  75.5      Yes 20.42       No     No      85
    ## 4129 Female     202  48 111.5  72.0       No 26.18       No     No      77
    ## 4130 Female     200  47 126.0  86.0       No 26.32       No     No      92
    ## 4131 Female     291  53 137.0  83.0       No 38.94       No    Yes      73
    ## 4132 Female     174  44 174.0 130.0      Yes 33.99       No     No      63
    ## 4133 Female     212  62 173.0  67.0       No 27.88       No     No      97
    ## 4134   Male     217  63 110.0  68.0      Yes 21.99       No     No      68
    ## 4135   Male     313  59 150.0  82.0       No 27.27       No     No      94
    ## 4136 Female     240  41 118.0  81.0       No 25.48       No     No      80
    ## 4138   Male     291  39 177.5 100.0      Yes 25.74       No     No      91
    ## 4139   Male     211  40 117.0  76.0      Yes 25.74       No     No      79
    ## 4142 Female     304  66 161.0  90.0       No 23.48       No     No      57
    ## 4143 Female     292  52 125.0  87.0       No 31.92       No     No      67
    ## 4144 Female     213  39 125.0  87.0       No 16.73       No     No      75
    ## 4145 Female     223  56 144.0  87.0       No 21.75       No     No      92
    ## 4146   Male     263  64 102.5  63.5       No 28.82       No     No      99
    ## 4147   Male     270  42 129.5 100.0      Yes 30.13       No     No      88
    ## 4148 Female     172  42 120.5  80.0       No 19.93       No     No      61
    ## 4149 Female     205  54  96.0  66.0      Yes 23.26       No     No      75
    ## 4150 Female     250  43 112.5  76.5      Yes 25.23       No     No      63
    ## 4151   Male     195  64 176.0  78.0       No 24.90      Yes     No     370
    ## 4152 Female     306  56 120.0  87.0       No 25.38       No     No      84
    ## 4153 Female     309  62 146.0  82.5       No 28.55       No     No      71
    ## 4154 Female     262  64 147.0  90.0       No 26.51      Yes     No     173
    ## 4155   Male     240  43 126.0  79.0      Yes 21.38       No     No      40
    ## 4156 Female     287  56 169.0  91.0       No 26.36       No     No      83
    ## 4157 Female     202  46 157.0  94.0      Yes 19.37       No     No      65
    ## 4158   Male     180  61 110.0  80.0      Yes 20.72       No     No      86
    ## 4160 Female     273  55 125.0  80.0       No 23.05       No     No      67
    ## 4161 Female     235  61 207.0 122.5       No 31.64       No     No      72
    ## 4162 Female     238  54 136.5  85.0       No 24.96       No     No      71
    ## 4163 Female     188  41 145.0  99.0      Yes 28.60       No     No      74
    ## 4164   Male     180  62 134.0  88.0      Yes 20.56       No     No      78
    ## 4165   Male     199  51 113.0  68.0       No 22.46       No     No      68
    ## 4166 Female     185  35 100.0  66.5      Yes 24.08       No     No      75
    ## 4167 Female     203  46 129.0  77.0       No 29.29       No     No      97
    ## 4171 Female     266  55 131.0  76.0       No 26.45       No     No      84
    ## 4172   Male     163  39 129.0  84.0      Yes 21.01       No     No     120
    ## 4173 Female     179  43 101.0  68.5      Yes 19.83       No     No      76
    ## 4174 Female     190  40 122.0  78.0       No 28.18       No     No      87
    ## 4175 Female     246  54 153.0  80.0      Yes 37.30       No     No      74
    ## 4176   Male     200  41 118.0  87.0      Yes 21.28       No     No      71
    ## 4177 Female     325  46 119.0  86.0       No 35.13       No     No      64
    ## 4178 Female     205  44 132.5  82.5      Yes 30.98       No     No      66
    ## 4179 Female     205  37 111.0  60.5      Yes 21.80       No     No      82
    ## 4180   Male     298  54 133.0  84.0      Yes 25.59       No     No      94
    ## 4181   Male     200  46 110.0  72.0      Yes 28.61       No     No      75
    ## 4182   Male     237  60 122.5  75.0      Yes 24.87      Yes     No     193
    ## 4183 Female     385  58 165.0  95.0       No 41.66       No     No      91
    ## 4184 Female     277  46 122.5  77.5       No 27.42       No     No      77
    ## 4186   Male     186  39 126.0  67.0      Yes 22.04       No     No      72
    ## 4187 Female     309  52 142.0  87.0       No 24.22       No     No     110
    ## 4188   Male     242  40 115.0  74.0      Yes 23.09       No     No      80
    ## 4189   Male     254  46 135.0 100.0       No 27.86       No     No      75
    ## 4191   Male     240  64 141.0  76.0      Yes 24.94       No     No      60
    ## 4192   Male     245  48 118.0  73.0      Yes 21.84       No     No     102
    ## 4193   Male     254  44 123.0  82.0      Yes 24.56       No     No      68
    ## 4195   Male     271  39 118.0  74.0      Yes 22.66       No     No      76
    ## 4196 Female     268  55 140.0  88.0       No 26.99       No     No     117
    ## 4197 Female     159  38 108.0  72.0       No 27.68       No     No      84
    ## 4198   Male     264  62 129.0  85.0      Yes 26.15       No     No      63
    ## 4200   Male     217  63 125.0  71.0      Yes 25.91       No     No      74
    ## 4201   Male     310  43 140.0  92.0      Yes 26.00       No     No      88
    ## 4202 Female     266  55 107.0  70.0       No 24.51       No     No      77
    ## 4203   Male     232  53 136.0  73.0      Yes 22.26       No     No      73
    ## 4204   Male     250  60 157.0  94.0      Yes 29.89       No     No      68
    ## 4205   Male     233  45 147.0 101.0      Yes 24.32       No     No      99
    ## 4206   Male     215  48 133.0  90.0      Yes 27.24       No     No      50
    ## 4207 Female     332  46 162.5  92.5      Yes 26.13       No     No      67
    ## 4208   Male     242  55 130.0  85.0      Yes 26.79       No     No      93
    ## 4209   Male     153  36 143.0  87.0       No 28.30       No     No      82
    ## 4211   Male     232  60 173.0 106.0      Yes 28.63       No    Yes      64
    ## 4212 Female     196  34 108.5  68.0       No 25.67       No     No      82
    ## 4213   Male     193  40 130.0  86.5      Yes 23.48       No     No      78
    ## 4214   Male     325  40 112.0  67.0      Yes 25.09       No     No      94
    ## 4215 Female     233  49 112.5  80.0      Yes 27.87       No     No      80
    ## 4216 Female     213  40 130.0  80.0      Yes 19.98       No     No      76
    ## 4217 Female     205  52 159.0 110.0      Yes 28.18       No     No      83
    ## 4218 Female     166  58 185.5 115.5       No 27.97       No     No      85
    ## 4219 Female     259  46 129.0  83.0      Yes 22.91       No     No      84
    ## 4220 Female     181  39 103.0  62.0       No 20.68       No     No      69
    ## 4221   Male     225  52 126.0  75.0      Yes 22.18       No     No     100
    ## 4222   Male     226  66 213.0 133.0       No 25.29       No     No      67
    ## 4223 Female     199  56 160.0 105.0       No 25.71       No     No      83
    ## 4224   Male     208  44 175.0 101.0       No 27.93      Yes     No     193
    ## 4225   Male     296  56 111.5  74.0      Yes 23.38       No     No      71
    ## 4226 Female     195  39 129.5  93.5       No 34.84       No     No      85
    ## 4227   Male     177  56 124.0  77.0       No 27.81       No     No      88
    ## 4228   Male     214  46 110.0  73.0      Yes 18.10       No     No      60
    ## 4229   Male     315  37 118.0  79.5      Yes 22.52       No     No      70
    ## 4230 Female     345  59 148.0  95.0       No 23.72       No     No      60
    ## 4232   Male     262  50  97.5  62.5      Yes 21.55       No     No      84
    ## 4233   Male     206  41 124.0  89.0      Yes 27.63       No     No      70
    ## 4234   Male     253  39 104.0  64.0      Yes 27.48       No     No      74
    ## 4235 Female     160  44 107.0  69.0       No 18.63       No     No      78
    ## 4236 Female     193  48 127.0  81.0       No 25.85       No     No      70
    ## 4237   Male     311  44 115.0  80.0      Yes 25.43       No     No      90
    ## 4238   Male     254  44 130.0  80.0       No 28.15       No     No      74
    ## 4239 Female     272  57 157.0  80.0      Yes 25.15       No    Yes      95
    ## 4240 Female     224  48 192.5 115.0       No 23.03       No     No      85
    ## 4241   Male     239  50 135.0  80.0      Yes 25.56       No     No      67
    ## 4242   Male     198  55 176.0 109.0       No 28.45       No     No      85
    ## 4243   Male     197  38 121.0  82.0       No 24.57       No     No      74
    ## 4245 Female     210  61 112.0  72.0       No 24.69       No     No      77
    ## 4246 Female     246  46 142.5  95.0       No 24.28       No     No      99
    ## 4247   Male     155  41 107.0  81.0      Yes 20.96      Yes     No     191
    ## 4248 Female     232  53 116.0  71.5       No 21.31       No     No      83
    ## 4249   Male     119  47 117.0  78.5      Yes 26.40       No     No      78
    ## 4250   Male     209  42 130.0  86.0      Yes 24.01       No     No      55
    ## 4251   Male     240  42 133.5  97.5      Yes 28.94       No     No      73
    ## 4252 Female     264  51 135.0  83.0       No 26.68       No     No      74
    ## 4253 Female     202  40 158.0 103.0       No 28.35       No     No      80
    ## 4255 Female     280  52 205.5 123.5       No 35.03       No     No      93
    ## 4257 Female     296  52 140.0  93.0       No 26.81       No     No      74
    ## 4258   Male     275  53 109.0  79.0       No 27.75       No     No     104
    ## 4259   Male     241  51 118.0  79.0      Yes 20.09       No     No      56
    ## 4260 Female     220  40 136.0  80.5       No 22.72       No     No      70
    ## 4261   Male     203  63 192.5 125.0      Yes 26.18       No     No      83
    ## 4262   Male     176  54 131.0  84.0      Yes 25.80       No     No      72
    ## 4263 Female     265  61 200.0 125.0       No 29.50      Yes    Yes     256
    ## 4264   Male     222  43 115.0  72.5       No 25.46       No     No      69
    ## 4265   Male     214  42 120.0  81.0       No 28.47       No     No      77
    ## 4266   Male     260  62 104.0  69.0      Yes 24.02       No     No      93
    ## 4267   Male     240  58 124.0  70.0       No 25.56       No     No     122
    ## 4268 Female     195  36 116.0  76.0       No 22.16       No     No      77
    ## 4269 Female     253  63 128.0  81.0       No 26.92       No     No     102
    ## 4271 Female     247  46 115.0  71.0       No 27.72       No     No      69
    ## 4272 Female     326  54 187.0  95.0       No 29.94      Yes     No     235
    ## 4273 Female     239  58 121.0  81.0       No 35.19       No     No      70
    ## 4275 Female     262  56 126.0  74.0      Yes 27.35       No     No     115
    ## 4276   Male     290  47 143.0  88.0      Yes 27.40       No     No      77
    ## 4277 Female     232  64 149.5  84.0       No 20.49       No     No      96
    ## 4278   Male     220  51 151.0  87.5       No 22.01       No     No      86
    ## 4279   Male     198  43 116.0  74.0      Yes 23.99       No     No      78
    ## 4280 Female     221  40  93.0  62.5      Yes 18.84       No     No      73
    ## 4281 Female     229  46 125.0  80.0       No 27.27       No     No      80
    ## 4282   Male     180  64 121.0  79.0      Yes 20.45       No     No      92
    ## 4283   Male     213  60 140.5  83.0       No 28.59       No     No      69
    ## 4284 Female     273  56 136.0  80.0       No 27.73      Yes     No     210
    ## 4285   Male     245  54 152.0  82.0      Yes 23.71       No     No      75
    ## 4287 Female     271  61 122.0  67.5       No 22.02       No     No      73
    ## 4288   Male     245  69 123.0  77.0      Yes 26.58       No     No      81
    ## 4289   Male     266  59 133.0  97.0       No 28.23       No     No      70
    ## 4290 Female     217  56 134.0  75.0      Yes 29.59       No     No      92
    ## 4291   Male     256  50 130.0  84.0      Yes 28.67       No     No      67
    ## 4292 Female     238  51 123.0  80.0       No 22.19       No     No     100
    ## 4293 Female     175  38 142.0  86.0       No 22.01       No     No      73
    ## 4294   Male     171  43 110.0  74.5      Yes 25.09       No     No      85
    ## 4295 Female     282  46 176.0  98.0      Yes 33.02       No    Yes      78
    ## 4296 Female     199  43 104.0  79.0      Yes 20.12       No     No      64
    ## 4297 Female     260  56 158.0 102.5       No 26.89       No    Yes      88
    ## 4298   Male     299  54 146.5  92.0      Yes 26.38       No     No      71
    ## 4299   Male     220  45 139.0  74.0       No 27.34       No     No     100
    ## 4300 Female     257  61 141.0  80.0       No 33.90       No     No      60
    ## 4301   Male     230  49 108.0  68.0       No 26.17       No     No      82
    ## 4302 Female     270  64 142.0  68.0      Yes 21.32       No     No      80
    ## 4303 Female     280  65 115.0  73.0       No 19.76       No     No      58
    ## 4304 Female     250  63 190.0  88.0       No 24.16       No     No     118
    ## 4305   Male     213  40 149.0  83.0       No 26.68       No     No      94
    ## 4306   Male     283  47 146.5  97.5       No 26.25       No     No      73
    ## 4307 Female     160  37 137.0  82.0       No 21.03       No     No     113
    ## 4308   Male     195  51 122.0  72.0      Yes 21.51       No     No      64
    ## 4309 Female     255  60 140.0  77.0      Yes 36.29       No    Yes     100
    ## 4310   Male     169  40 127.0  81.0      Yes 25.82       No     No      83
    ## 4313   Male     268  51 206.0 116.0       No 26.35       No     No      70
    ## 4314   Male     201  57 108.5  70.5       No 22.90       No     No      84
    ## 4315 Female     232  64 164.0 102.0       No 33.37       No     No      73
    ## 4316   Male     216  38 124.0  84.0       No 28.12       No     No      75
    ## 4317 Female     287  48 155.0  80.0       No 28.54       No     No      75
    ## 4318   Male     228  53 114.0  80.0       No 20.54       No     No      75
    ## 4319   Male     176  57 147.5  87.5       No 24.15       No     No     100
    ## 4320 Female     291  45 125.0  82.0      Yes 21.26       No     No      72
    ## 4321   Male     210  57 120.0  77.5      Yes 27.14       No     No      71
    ## 4322 Female     259  57 170.0 101.0       No 38.17       No     No      75
    ## 4323 Female     195  39  97.5  60.0       No 20.62       No     No      68
    ## 4324   Male     247  46 134.0  96.0      Yes 32.47       No     No      72
    ## 4325 Female     179  42 115.0  78.0       No 25.75       No     No      77
    ## 4326 Female     235  42 128.0  86.0      Yes 24.05       No     No      70
    ## 4327 Female     238  42 118.0  80.0       No 33.19       No     No      76
    ## 4330 Female     186  56 126.0  73.0      Yes 25.81       No     No      65
    ## 4331   Male     305  49 135.0  82.0      Yes 26.29       No     No      65
    ## 4332   Male     206  48 118.0  81.0       No 28.13       No     No      87
    ## 4333   Male     234  47 162.0 110.0      Yes 27.51       No     No      85
    ## 4334   Male     231  55 105.0  82.0       No 27.73       No     No      66
    ## 4335   Male     237  59 131.5  84.0       No 24.17       No     No      94
    ## 4336   Male     239  55 144.0  96.5      Yes 28.82       No     No     102
    ## 4337   Male     235  42 103.0  70.0       No 21.48       No     No      73
    ## 4338 Female     202  48 128.0  74.0      Yes 25.11       No     No      75
    ## 4339   Male     226  42 108.0  78.0      Yes 22.54       No     No      68
    ## 4340 Female     198  55 136.0  93.0       No 22.54       No     No      83
    ## 4341 Female     207  62 127.5  75.0       No 22.91       No     No      80
    ## 4342 Female     257  36 103.0  72.5      Yes 27.86       No     No      65
    ## 4343 Female     212  57 147.0  85.0       No 35.19       No     No      85
    ## 4345   Male     269  52 157.5  83.0       No 26.60      Yes     No      80
    ## 4346   Male     227  47 126.0  84.0      Yes 19.14       No     No      74
    ## 4347 Female     227  38  99.0  62.0       No 27.16       No     No      90
    ## 4348 Female     325  61 125.5  85.5       No 24.40       No     No      70
    ## 4349 Female     220  69 143.0  81.0       No 26.27       No    Yes      77
    ## 4350   Male     241  62 135.0  97.5       No 24.88       No     No      96
    ## 4352 Female     224  37 109.0  72.0      Yes 22.81       No     No      93
    ## 4354   Male     236  65 118.5  77.5      Yes 24.30       No     No      65
    ## 4355   Male     287  39 136.0  86.0      Yes 19.00       No     No      83
    ## 4357 Female     230  44 130.0  81.5       No 25.74       No     No      77
    ## 4358   Male     194  57 133.0  78.0       No 29.02       No     No      92
    ## 4359 Female     235  45 106.0  58.0       No 26.79       No     No      79
    ## 4360 Female     325  62 180.0 108.0       No 35.16       No    Yes      81
    ## 4361 Female     214  38 101.0  70.0      Yes 21.83       No     No      77
    ## 4363   Male     249  41 125.0  87.0      Yes 27.13       No     No      81
    ## 4364   Male     223  58 154.5  88.0       No 27.91       No     No      81
    ## 4365 Female     302  54 210.0 127.5       No 31.98       No     No      79
    ## 4366   Male     267  51 129.5  80.0      Yes 25.98       No     No      79
    ## 4367   Male     215  51 135.0  94.0      Yes 30.61       No     No      78
    ## 4368   Male     266  52 107.0  75.0       No 25.64       No     No      98
    ## 4369 Female     199  41 111.0  70.0       No 21.99       No     No      82
    ## 4370 Female     218  52 121.5  57.0      Yes 20.78       No     No      85
    ## 4371   Male     210  41 132.5  85.0      Yes 28.62       No     No      70
    ## 4372 Female     207  54 137.5  89.0       No 25.43       No     No      72
    ## 4373   Male     285  43 129.0  95.0       No 26.64       No     No      74
    ## 4374   Male     235  38 121.0  83.0       No 25.85       No     No      75
    ## 4375   Male     196  65 157.0  86.0       No 26.36       No     No      80
    ## 4376 Female     252  60 182.0  99.0       No 40.23       No     No      60
    ## 4377 Female     208  37 118.5  70.0       No 25.09       No     No      85
    ## 4378   Male     232  44 137.5  87.5      Yes 30.03       No     No      70
    ## 4380 Female     215  49 106.0  63.0      Yes 19.22       No     No      66
    ## 4381 Female     199  42 143.0  79.0      Yes 18.68       No     No      76
    ## 4382 Female     180  44 110.0  70.0       No 23.98       No     No      67
    ## 4383   Male     261  54 117.0  74.0      Yes 20.88       No     No      77
    ## 4384 Female     229  41 150.0  89.0       No 36.07       No     No      92
    ## 4385 Female     173  37 101.0  69.0      Yes 20.02       No     No      73
    ## 4386   Male     241  45 129.0  80.0      Yes 27.11       No     No      65
    ## 4387 Female     306  63 195.0 105.0       No 27.96       No     No      87
    ## 4388   Male     246  65 179.0  96.0      Yes 19.34       No    Yes      76
    ## 4389 Female     207  40 124.0  78.0      Yes 22.90       No     No      66
    ## 4390 Female     222  46 116.0  80.0       No 24.62       No     No      87
    ## 4391 Female     257  63 170.0 105.0       No 25.49       No     No      87
    ## 4392 Female     240  52 157.5 105.0      Yes 29.64       No     No      80
    ## 4393 Female     271  61 133.0  83.0       No 25.31       No     No      60
    ## 4394 Female     169  38 115.0  60.0       No 26.87       No     No      60
    ## 4395   Male     214  62 130.0  80.0       No 24.35       No     No      77
    ## 4396 Female     236  42 160.5 100.0       No 27.01       No     No      74
    ## 4397   Male     240  63 146.0  84.0      Yes 30.48      Yes     No     120
    ## 4398 Female     216  38 142.0  72.0       No 22.01       No     No      75
    ## 4399 Female     213  46 136.0  77.0       No 31.02       No     No      73
    ## 4400 Female     190  43 103.0  67.5      Yes 24.08       No     No      69
    ## 4401 Female     261  61 124.0  76.5       No 23.06       No     No      83
    ## 4403   Male     286  65 135.0  80.0       No 28.06       No     No     116
    ## 4404   Male     282  50 126.5  88.0       No 27.30       No     No      87
    ## 4405   Male     245  42 105.0  70.0      Yes 22.41       No     No      69
    ## 4406   Male     209  48 144.0  88.0      Yes 29.11       No     No      60
    ## 4407   Male     193  48 141.0  95.0       No 27.89       No     No      84
    ## 4408 Female     230  56 123.0  78.5      Yes 24.71       No     No      87
    ## 4409 Female     236  63 155.0  82.0       No 39.17      Yes     No      79
    ## 4410 Female     251  59 125.0  80.0      Yes 22.18       No     No      70
    ## 4411 Female     266  63 167.0  94.0       No 25.23       No     No      94
    ## 4412 Female     265  52 137.5  84.5       No 26.91       No     No      86
    ## 4413   Male     204  48 125.0  84.5       No 22.37       No     No      75
    ## 4414   Male     252  60 128.0  82.0       No 21.18       No     No      70
    ## 4415   Male     260  50 119.0  74.0       No 21.85       No     No      72
    ## 4416   Male     289  53 188.0 110.0       No 26.70       No     No      63
    ## 4417   Male     287  56 149.0  98.0       No 21.68       No    Yes      75
    ## 4418   Male     198  47 120.0  80.0      Yes 25.23       No     No      76
    ## 4419   Male     216  45 137.5  85.0      Yes 24.24       No     No     105
    ## 4420   Male     217  43 122.0  84.0       No 32.18       No     No      88
    ## 4421   Male     233  58 125.5  84.0       No 26.05       No     No      76
    ## 4422   Male     187  43 129.5  88.0      Yes 25.62       No     No      75
    ## 4423 Female     260  50 190.0 130.0       No 43.67      Yes     No     260
    ## 4426   Male     187  58 141.0  81.0       No 24.96       No     No      81
    ## 4427   Male     176  68 168.0  97.0       No 23.14       No     No      79
    ## 4428   Male     313  50 179.0  92.0      Yes 25.97       No     No      86
    ## 4429   Male     207  51 126.5  80.0      Yes 19.71       No     No      68
    ## 4432 Female     269  52 133.5  83.0       No 21.47       No     No     107
    ## 4433   Male     185  40 141.0  98.0       No 25.60       No     No      72
    ## 4434 Female     196  39 133.0  86.0      Yes 20.91       No     No      80
    ##      EDUC
    ## 1       4
    ## 2       2
    ## 3       1
    ## 4       3
    ## 5       3
    ## 6       2
    ## 7       1
    ## 8       2
    ## 9       1
    ## 10      1
    ## 11      1
    ## 12      2
    ## 13      1
    ## 14      3
    ## 16      2
    ## 17      3
    ## 18      2
    ## 19      2
    ## 20      2
    ## 21      2
    ## 23      1
    ## 24      3
    ## 25      2
    ## 26      4
    ## 28      2
    ## 29      3
    ## 30      1
    ## 31      4
    ## 32      4
    ## 33      3
    ## 34      1
    ## 36      1
    ## 37      1
    ## 38      2
    ## 40      1
    ## 41      1
    ## 42      3
    ## 43      2
    ## 44      2
    ## 46      2
    ## 47      3
    ## 48      2
    ## 49      1
    ## 50      1
    ## 51      2
    ## 53      1
    ## 54      4
    ## 55      2
    ## 56      1
    ## 58      1
    ## 59      1
    ## 60      1
    ## 61      2
    ## 62      4
    ## 63      4
    ## 64      1
    ## 65      1
    ## 66      1
    ## 67      1
    ## 68      1
    ## 69      1
    ## 70      1
    ## 71      2
    ## 72      4
    ## 74      3
    ## 76      1
    ## 77      2
    ## 78      1
    ## 79      1
    ## 81      1
    ## 82      2
    ## 83      1
    ## 84      2
    ## 85      2
    ## 86      1
    ## 87      1
    ## 88      2
    ## 89      1
    ## 90      2
    ## 91      4
    ## 92      1
    ## 93      2
    ## 94      1
    ## 95      2
    ## 96      1
    ## 97      1
    ## 98      3
    ## 99      4
    ## 100     3
    ## 101     4
    ## 103     1
    ## 104     1
    ## 105     2
    ## 106     3
    ## 107     3
    ## 108     1
    ## 109     1
    ## 110     1
    ## 111     1
    ## 112     4
    ## 113     1
    ## 114     1
    ## 115     1
    ## 116     1
    ## 118     4
    ## 119     1
    ## 121     1
    ## 122     1
    ## 123     2
    ## 124     4
    ## 125     3
    ## 126     4
    ## 127     1
    ## 128     3
    ## 129     1
    ## 130     2
    ## 131     2
    ## 132     2
    ## 133     4
    ## 134     1
    ## 135     2
    ## 136     1
    ## 138     2
    ## 139     2
    ## 140     3
    ## 141     1
    ## 142     2
    ## 143     2
    ## 144     3
    ## 145     4
    ## 146     1
    ## 147     1
    ## 148     3
    ## 149     1
    ## 150     1
    ## 151     4
    ## 152     1
    ## 153     3
    ## 154     4
    ## 155     1
    ## 156     1
    ## 157     1
    ## 158     1
    ## 159     1
    ## 161     1
    ## 162     2
    ## 163     1
    ## 164     3
    ## 165     3
    ## 166     3
    ## 167     4
    ## 168     4
    ## 169     2
    ## 170     2
    ## 171     1
    ## 172     1
    ## 173     2
    ## 174     4
    ## 175     4
    ## 176     1
    ## 177     2
    ## 178     1
    ## 179     2
    ## 180     2
    ## 181     3
    ## 182     1
    ## 183     1
    ## 184     4
    ## 185     3
    ## 186     1
    ## 187     1
    ## 188     1
    ## 189     1
    ## 190     3
    ## 191     1
    ## 192     1
    ## 194     2
    ## 195     4
    ## 196     1
    ## 197     2
    ## 198     1
    ## 199     2
    ## 200     1
    ## 201     1
    ## 203     3
    ## 204     2
    ## 205     1
    ## 206     2
    ## 207     2
    ## 208     1
    ## 209     2
    ## 210     4
    ## 211     2
    ## 213     4
    ## 215     2
    ## 216     1
    ## 217     1
    ## 218     2
    ## 219     2
    ## 220     4
    ## 222     1
    ## 224     3
    ## 227     1
    ## 228     2
    ## 229     2
    ## 230     2
    ## 231     2
    ## 232     1
    ## 233     2
    ## 234     1
    ## 235     1
    ## 236     2
    ## 237     1
    ## 238     1
    ## 239     3
    ## 240     1
    ## 241     4
    ## 242     1
    ## 243     2
    ## 244     2
    ## 245     1
    ## 246     1
    ## 247     1
    ## 248     4
    ## 249     1
    ## 250     3
    ## 251     3
    ## 252     4
    ## 253     1
    ## 254     3
    ## 255     3
    ## 256     2
    ## 257     4
    ## 261     2
    ## 262     1
    ## 264     3
    ## 265     1
    ## 266     1
    ## 267     4
    ## 268     1
    ## 269     4
    ## 270     1
    ## 271     1
    ## 272     1
    ## 273     3
    ## 274     1
    ## 275     1
    ## 277     1
    ## 278     2
    ## 279     1
    ## 280     1
    ## 281     2
    ## 282     2
    ## 283     1
    ## 284     2
    ## 285     4
    ## 287     1
    ## 288     2
    ## 289     1
    ## 290     1
    ## 291     2
    ## 293     2
    ## 294     3
    ## 296     4
    ## 297     1
    ## 298     2
    ## 299     2
    ## 300     2
    ## 301     2
    ## 302     4
    ## 303     3
    ## 304     4
    ## 305     4
    ## 308     2
    ## 310     3
    ## 311     2
    ## 312     2
    ## 313     1
    ## 316     1
    ## 317     2
    ## 320     1
    ## 321     3
    ## 323     3
    ## 324     1
    ## 325     1
    ## 326     1
    ## 329     1
    ## 330     2
    ## 331     3
    ## 333     1
    ## 334     2
    ## 335     2
    ## 336     3
    ## 337     3
    ## 338     3
    ## 339     2
    ## 340     1
    ## 341     1
    ## 343     1
    ## 344     1
    ## 345     1
    ## 346     2
    ## 347     1
    ## 348     3
    ## 349     1
    ## 350     1
    ## 352     1
    ## 353     1
    ## 354     1
    ## 355     4
    ## 358     1
    ## 360     2
    ## 361     1
    ## 362     2
    ## 363     1
    ## 364     3
    ## 365     2
    ## 366     1
    ## 367     1
    ## 368     2
    ## 370     1
    ## 371     2
    ## 372     3
    ## 373     1
    ## 374     2
    ## 375     3
    ## 376     4
    ## 377     1
    ## 378     2
    ## 379     2
    ## 380     3
    ## 381     1
    ## 382     1
    ## 383     1
    ## 384     2
    ## 385     4
    ## 386     1
    ## 387     2
    ## 388     4
    ## 389     3
    ## 390     2
    ## 391     1
    ## 392     2
    ## 393     3
    ## 394     4
    ## 395     2
    ## 396     3
    ## 398     1
    ## 399     1
    ## 400     2
    ## 401     2
    ## 402     3
    ## 403     2
    ## 404     2
    ## 405     1
    ## 406     1
    ## 407     3
    ## 408     1
    ## 409     2
    ## 411     3
    ## 412     2
    ## 413     1
    ## 414     1
    ## 416     3
    ## 417     2
    ## 418     1
    ## 419     4
    ## 420     1
    ## 421     3
    ## 423     4
    ## 424     2
    ## 425     1
    ## 426     3
    ## 427     3
    ## 430     2
    ## 431     2
    ## 432     1
    ## 433     1
    ## 434     2
    ## 436     1
    ## 437     2
    ## 439     2
    ## 440     3
    ## 441     2
    ## 442     1
    ## 444     1
    ## 448     2
    ## 449     1
    ## 450     1
    ## 452     4
    ## 453     2
    ## 455     1
    ## 456     2
    ## 457     1
    ## 458     3
    ## 459     1
    ## 460     4
    ## 461     1
    ## 462     3
    ## 463     4
    ## 464     1
    ## 465     4
    ## 467     3
    ## 468     3
    ## 469     2
    ## 470     4
    ## 471     1
    ## 472     1
    ## 475     1
    ## 476     2
    ## 477     4
    ## 478     2
    ## 479     2
    ## 480     1
    ## 481     2
    ## 482     1
    ## 483     2
    ## 484     3
    ## 486     2
    ## 487     1
    ## 488     1
    ## 489     1
    ## 491     1
    ## 492     2
    ## 493     3
    ## 494     1
    ## 495     3
    ## 496     2
    ## 497     1
    ## 498     1
    ## 499     4
    ## 500     2
    ## 501     1
    ## 502     2
    ## 503     2
    ## 504     1
    ## 505     2
    ## 506     2
    ## 508     1
    ## 509     2
    ## 511     4
    ## 512     3
    ## 513     4
    ## 514     3
    ## 515     1
    ## 516     2
    ## 517     3
    ## 518     1
    ## 519     4
    ## 521     3
    ## 522     1
    ## 523     1
    ## 525     3
    ## 526     3
    ## 527     4
    ## 528     3
    ## 529     1
    ## 530     3
    ## 531     2
    ## 533     2
    ## 534     1
    ## 535     1
    ## 536     1
    ## 537     1
    ## 538     2
    ## 539     2
    ## 541     3
    ## 542     4
    ## 543     4
    ## 544     3
    ## 545     1
    ## 546     2
    ## 547     3
    ## 548     1
    ## 549     3
    ## 550     3
    ## 551     2
    ## 552     1
    ## 553     1
    ## 554     2
    ## 555     1
    ## 556     1
    ## 557     3
    ## 558     2
    ## 559     2
    ## 560     3
    ## 562     3
    ## 563     2
    ## 564     2
    ## 565     3
    ## 566     1
    ## 567     2
    ## 568     1
    ## 569     3
    ## 570     1
    ## 571     1
    ## 572     1
    ## 574     2
    ## 575     1
    ## 576     1
    ## 577     2
    ## 578     2
    ## 579     2
    ## 580     3
    ## 581     4
    ## 582     3
    ## 583     1
    ## 584     2
    ## 586     1
    ## 587     4
    ## 588     1
    ## 589     4
    ## 591     2
    ## 592     3
    ## 593     3
    ## 594     1
    ## 595     2
    ## 596     4
    ## 597     3
    ## 598     1
    ## 601     2
    ## 602     1
    ## 603     1
    ## 604     1
    ## 605     2
    ## 606     4
    ## 607     1
    ## 608     1
    ## 609     2
    ## 610     2
    ## 611     1
    ## 612     1
    ## 613     2
    ## 615     4
    ## 616     1
    ## 617     3
    ## 618     1
    ## 619     2
    ## 620     3
    ## 621     3
    ## 622     2
    ## 623     2
    ## 624     1
    ## 625     1
    ## 626     1
    ## 627     1
    ## 628     2
    ## 629     2
    ## 630     3
    ## 631     1
    ## 632     1
    ## 634     2
    ## 635     2
    ## 636     2
    ## 637     2
    ## 638     2
    ## 639     4
    ## 640     2
    ## 641     3
    ## 642     2
    ## 643     1
    ## 644     1
    ## 645     2
    ## 647     3
    ## 648     3
    ## 649     2
    ## 650     1
    ## 651     1
    ## 652     2
    ## 653     2
    ## 654     4
    ## 655     1
    ## 656     2
    ## 657     2
    ## 658     1
    ## 659     2
    ## 660     2
    ## 661     4
    ## 662     2
    ## 663     1
    ## 664     2
    ## 665     2
    ## 666     2
    ## 667     2
    ## 668     4
    ## 669     4
    ## 671     2
    ## 672     4
    ## 673     1
    ## 674     1
    ## 675     3
    ## 676     1
    ## 677     1
    ## 678     1
    ## 679     1
    ## 680     3
    ## 681     1
    ## 682     4
    ## 683     1
    ## 684     1
    ## 685     4
    ## 687     2
    ## 688     2
    ## 689     1
    ## 690     4
    ## 691     1
    ## 692     1
    ## 693     1
    ## 694     1
    ## 695     2
    ## 696     1
    ## 697     4
    ## 699     1
    ## 701     2
    ## 702     4
    ## 703     1
    ## 705     1
    ## 706     1
    ## 707     4
    ## 708     1
    ## 709     3
    ## 710     2
    ## 711     2
    ## 712     3
    ## 713     2
    ## 714     1
    ## 716     1
    ## 717     1
    ## 718     2
    ## 719     4
    ## 721     2
    ## 722     2
    ## 723     1
    ## 724     3
    ## 725     1
    ## 726     1
    ## 727     1
    ## 728     1
    ## 729     3
    ## 730     1
    ## 732     1
    ## 733     1
    ## 734     2
    ## 735     1
    ## 736     4
    ## 737     1
    ## 738     2
    ## 739     2
    ## 741     1
    ## 742     2
    ## 743     2
    ## 744     1
    ## 746     2
    ## 747     1
    ## 748     1
    ## 749     2
    ## 750     1
    ## 751     3
    ## 752     4
    ## 753     3
    ## 754     1
    ## 755     3
    ## 756     3
    ## 757     1
    ## 758     1
    ## 759     1
    ## 760     1
    ## 761     1
    ## 762     2
    ## 763     3
    ## 765     4
    ## 766     2
    ## 767     1
    ## 768     1
    ## 769     2
    ## 770     4
    ## 771     3
    ## 772     2
    ## 773     1
    ## 774     2
    ## 775     1
    ## 776     2
    ## 777     1
    ## 778     1
    ## 779     2
    ## 780     2
    ## 781     2
    ## 782     1
    ## 784     2
    ## 786     1
    ## 787     2
    ## 789     3
    ## 790     4
    ## 791     1
    ## 792     3
    ## 794     1
    ## 795     1
    ## 796     1
    ## 798     2
    ## 799     4
    ## 800     2
    ## 801     2
    ## 802     2
    ## 803     1
    ## 804     3
    ## 808     3
    ## 809     3
    ## 811     3
    ## 812     2
    ## 813     1
    ## 814     3
    ## 815     3
    ## 816     1
    ## 817     1
    ## 818     1
    ## 819     3
    ## 820     1
    ## 821     4
    ## 822     1
    ## 823     2
    ## 824     2
    ## 825     2
    ## 827     2
    ## 828     3
    ## 829     3
    ## 830     1
    ## 831     2
    ## 832     1
    ## 833     2
    ## 834     2
    ## 835     3
    ## 836     3
    ## 837     1
    ## 838     2
    ## 839     1
    ## 840     4
    ## 841     1
    ## 844     4
    ## 845     1
    ## 846     2
    ## 847     2
    ## 848     4
    ## 850     2
    ## 851     1
    ## 853     1
    ## 854     2
    ## 855     2
    ## 856     3
    ## 857     1
    ## 858     2
    ## 860     1
    ## 861     2
    ## 863     1
    ## 865     2
    ## 866     2
    ## 867     4
    ## 868     3
    ## 869     1
    ## 870     1
    ## 871     1
    ## 872     3
    ## 873     2
    ## 874     2
    ## 875     2
    ## 876     1
    ## 877     1
    ## 878     2
    ## 879     2
    ## 880     2
    ## 881     1
    ## 882     4
    ## 883     2
    ## 884     3
    ## 885     1
    ## 886     2
    ## 887     1
    ## 888     1
    ## 889     1
    ## 890     2
    ## 891     2
    ## 892     2
    ## 893     2
    ## 894     2
    ## 895     1
    ## 896     2
    ## 897     2
    ## 898     2
    ## 899     1
    ## 900     4
    ## 901     2
    ## 904     1
    ## 905     3
    ## 906     2
    ## 907     1
    ## 908     1
    ## 909     1
    ## 910     4
    ## 911     4
    ## 912     1
    ## 913     4
    ## 915     1
    ## 916     1
    ## 917     1
    ## 918     3
    ## 919     1
    ## 920     4
    ## 921     4
    ## 922     1
    ## 923     2
    ## 924     1
    ## 925     1
    ## 926     1
    ## 927     1
    ## 929     3
    ## 930     2
    ## 931     1
    ## 932     3
    ## 933     2
    ## 934     1
    ## 935     1
    ## 936     3
    ## 938     1
    ## 939     4
    ## 940     2
    ## 941     1
    ## 942     1
    ## 943     3
    ## 944     4
    ## 945     2
    ## 946     1
    ## 947     1
    ## 948     1
    ## 950     1
    ## 952     1
    ## 953     1
    ## 954     1
    ## 955     4
    ## 957     1
    ## 958     2
    ## 959     1
    ## 960     3
    ## 961     1
    ## 962     1
    ## 963     4
    ## 964     1
    ## 965     1
    ## 966     2
    ## 967     2
    ## 968     2
    ## 969     2
    ## 970     1
    ## 971     1
    ## 973     2
    ## 974     4
    ## 975     2
    ## 977     2
    ## 978     1
    ## 981     1
    ## 982     3
    ## 983     3
    ## 984     4
    ## 985     2
    ## 986     2
    ## 987     3
    ## 989     1
    ## 990     1
    ## 991     4
    ## 992     4
    ## 993     2
    ## 994     1
    ## 995     1
    ## 996     3
    ## 997     2
    ## 998     1
    ## 999     2
    ## 1000    2
    ## 1001    2
    ## 1002    1
    ## 1003    2
    ## 1004    3
    ## 1006    1
    ## 1007    1
    ## 1008    1
    ## 1009    1
    ## 1010    3
    ## 1012    1
    ## 1013    3
    ## 1014    1
    ## 1015    1
    ## 1016    1
    ## 1017    1
    ## 1018    1
    ## 1020    2
    ## 1021    1
    ## 1022    1
    ## 1023    3
    ## 1024    2
    ## 1025    1
    ## 1026    4
    ## 1027    1
    ## 1028    2
    ## 1029    2
    ## 1030    2
    ## 1031    1
    ## 1032    3
    ## 1033    2
    ## 1035    1
    ## 1036    3
    ## 1039    2
    ## 1040    3
    ## 1042    1
    ## 1043    1
    ## 1044    1
    ## 1045    4
    ## 1046    3
    ## 1047    1
    ## 1048    2
    ## 1049    1
    ## 1050    1
    ## 1051    2
    ## 1052    2
    ## 1053    4
    ## 1054    1
    ## 1055    3
    ## 1056    3
    ## 1057    2
    ## 1059    2
    ## 1060    3
    ## 1061    1
    ## 1062    1
    ## 1063    4
    ## 1064    2
    ## 1065    2
    ## 1068    3
    ## 1069    2
    ## 1070    1
    ## 1071    4
    ## 1072    2
    ## 1073    1
    ## 1074    4
    ## 1075    3
    ## 1076    1
    ## 1078    2
    ## 1079    1
    ## 1080    1
    ## 1081    2
    ## 1082    1
    ## 1084    1
    ## 1085    1
    ## 1086    2
    ## 1087    1
    ## 1088    4
    ## 1089    1
    ## 1090    2
    ## 1091    1
    ## 1092    2
    ## 1093    1
    ## 1094    2
    ## 1095    1
    ## 1096    1
    ## 1097    2
    ## 1098    4
    ## 1099    1
    ## 1100    2
    ## 1101    3
    ## 1102    2
    ## 1105    1
    ## 1106    1
    ## 1107    2
    ## 1108    1
    ## 1109    1
    ## 1110    1
    ## 1112    1
    ## 1113    1
    ## 1115    1
    ## 1116    1
    ## 1118    2
    ## 1119    3
    ## 1120    2
    ## 1121    3
    ## 1123    1
    ## 1124    1
    ## 1125    3
    ## 1126    2
    ## 1127    2
    ## 1128    2
    ## 1129    1
    ## 1130    2
    ## 1131    2
    ## 1132    4
    ## 1134    1
    ## 1135    1
    ## 1136    1
    ## 1137    2
    ## 1138    4
    ## 1139    1
    ## 1140    1
    ## 1141    1
    ## 1142    1
    ## 1143    4
    ## 1144    1
    ## 1145    4
    ## 1147    2
    ## 1148    4
    ## 1151    2
    ## 1152    1
    ## 1153    1
    ## 1154    1
    ## 1155    3
    ## 1156    2
    ## 1157    2
    ## 1158    3
    ## 1159    2
    ## 1160    2
    ## 1162    1
    ## 1163    1
    ## 1166    3
    ## 1167    1
    ## 1169    4
    ## 1170    2
    ## 1171    3
    ## 1172    1
    ## 1173    3
    ## 1174    3
    ## 1176    1
    ## 1177    2
    ## 1178    2
    ## 1179    2
    ## 1180    1
    ## 1181    2
    ## 1182    2
    ## 1183    4
    ## 1184    1
    ## 1186    2
    ## 1187    2
    ## 1188    1
    ## 1190    1
    ## 1191    3
    ## 1192    1
    ## 1193    3
    ## 1194    1
    ## 1195    1
    ## 1196    3
    ## 1197    2
    ## 1198    1
    ## 1199    2
    ## 1200    2
    ## 1201    3
    ## 1202    2
    ## 1204    1
    ## 1205    4
    ## 1206    3
    ## 1207    3
    ## 1208    1
    ## 1209    2
    ## 1211    1
    ## 1212    1
    ## 1213    4
    ## 1214    2
    ## 1215    2
    ## 1216    4
    ## 1217    3
    ## 1218    1
    ## 1219    2
    ## 1220    1
    ## 1221    3
    ## 1223    3
    ## 1224    2
    ## 1226    2
    ## 1227    2
    ## 1229    1
    ## 1231    4
    ## 1233    2
    ## 1234    1
    ## 1235    2
    ## 1236    1
    ## 1237    3
    ## 1238    1
    ## 1239    3
    ## 1240    1
    ## 1242    2
    ## 1243    1
    ## 1244    1
    ## 1245    1
    ## 1246    1
    ## 1247    1
    ## 1248    1
    ## 1249    4
    ## 1250    1
    ## 1251    1
    ## 1252    1
    ## 1253    1
    ## 1255    1
    ## 1256    3
    ## 1257    2
    ## 1259    4
    ## 1260    1
    ## 1261    3
    ## 1262    4
    ## 1263    1
    ## 1265    1
    ## 1266    1
    ## 1268    3
    ## 1269    2
    ## 1270    1
    ## 1271    1
    ## 1272    1
    ## 1273    1
    ## 1274    1
    ## 1275    3
    ## 1276    2
    ## 1277    2
    ## 1278    2
    ## 1279    1
    ## 1280    1
    ## 1281    2
    ## 1283    2
    ## 1284    2
    ## 1285    1
    ## 1286    4
    ## 1287    4
    ## 1288    1
    ## 1290    1
    ## 1291    1
    ## 1292    1
    ## 1293    2
    ## 1295    3
    ## 1296    1
    ## 1297    1
    ## 1298    2
    ## 1299    3
    ## 1300    2
    ## 1301    2
    ## 1302    2
    ## 1303    1
    ## 1304    3
    ## 1306    2
    ## 1307    3
    ## 1308    1
    ## 1310    1
    ## 1311    2
    ## 1312    3
    ## 1313    4
    ## 1314    1
    ## 1316    1
    ## 1317    1
    ## 1318    3
    ## 1319    1
    ## 1320    2
    ## 1321    1
    ## 1322    1
    ## 1323    3
    ## 1324    1
    ## 1325    1
    ## 1326    2
    ## 1327    1
    ## 1328    3
    ## 1329    1
    ## 1330    2
    ## 1332    1
    ## 1333    1
    ## 1334    1
    ## 1335    1
    ## 1336    2
    ## 1337    2
    ## 1338    2
    ## 1339    1
    ## 1341    2
    ## 1343    2
    ## 1345    1
    ## 1346    4
    ## 1347    3
    ## 1348    3
    ## 1349    1
    ## 1350    1
    ## 1351    2
    ## 1352    2
    ## 1354    1
    ## 1355    2
    ## 1356    1
    ## 1358    1
    ## 1359    1
    ## 1360    2
    ## 1361    1
    ## 1362    2
    ## 1363    2
    ## 1364    1
    ## 1365    2
    ## 1366    4
    ## 1368    1
    ## 1369    1
    ## 1370    1
    ## 1372    1
    ## 1373    1
    ## 1375    2
    ## 1376    2
    ## 1377    2
    ## 1378    1
    ## 1379    2
    ## 1380    1
    ## 1381    2
    ## 1382    1
    ## 1383    3
    ## 1384    2
    ## 1386    2
    ## 1387    1
    ## 1388    3
    ## 1389    4
    ## 1390    2
    ## 1391    1
    ## 1392    1
    ## 1394    2
    ## 1395    3
    ## 1396    1
    ## 1397    1
    ## 1398    1
    ## 1399    1
    ## 1400    1
    ## 1401    4
    ## 1402    3
    ## 1403    1
    ## 1404    2
    ## 1406    1
    ## 1407    4
    ## 1408    1
    ## 1409    3
    ## 1410    3
    ## 1411    2
    ## 1412    2
    ## 1414    1
    ## 1415    1
    ## 1416    1
    ## 1417    2
    ## 1419    2
    ## 1421    2
    ## 1422    1
    ## 1423    2
    ## 1424    1
    ## 1425    4
    ## 1426    2
    ## 1428    1
    ## 1429    2
    ## 1430    2
    ## 1431    2
    ## 1432    1
    ## 1433    1
    ## 1435    2
    ## 1437    2
    ## 1438    2
    ## 1440    1
    ## 1441    3
    ## 1442    3
    ## 1443    3
    ## 1444    4
    ## 1445    1
    ## 1446    3
    ## 1447    3
    ## 1448    2
    ## 1449    2
    ## 1450    4
    ## 1451    2
    ## 1452    1
    ## 1453    1
    ## 1456    3
    ## 1457    2
    ## 1458    4
    ## 1459    1
    ## 1461    4
    ## 1462    1
    ## 1463    3
    ## 1464    2
    ## 1465    4
    ## 1466    2
    ## 1467    1
    ## 1469    1
    ## 1470    2
    ## 1471    4
    ## 1473    2
    ## 1474    1
    ## 1475    1
    ## 1476    1
    ## 1477    1
    ## 1478    1
    ## 1479    3
    ## 1480    2
    ## 1481    3
    ## 1482    3
    ## 1483    1
    ## 1484    2
    ## 1485    1
    ## 1486    2
    ## 1487    2
    ## 1488    3
    ## 1489    3
    ## 1490    1
    ## 1492    3
    ## 1493    1
    ## 1494    2
    ## 1495    1
    ## 1496    1
    ## 1497    2
    ## 1498    3
    ## 1499    1
    ## 1500    1
    ## 1501    2
    ## 1502    1
    ## 1503    2
    ## 1504    1
    ## 1505    1
    ## 1506    3
    ## 1507    3
    ## 1508    1
    ## 1509    1
    ## 1511    1
    ## 1513    2
    ## 1514    1
    ## 1515    1
    ## 1516    1
    ## 1517    2
    ## 1518    1
    ## 1519    1
    ## 1520    1
    ## 1521    2
    ## 1522    2
    ## 1523    2
    ## 1524    3
    ## 1525    1
    ## 1526    1
    ## 1527    1
    ## 1528    1
    ## 1529    4
    ## 1530    2
    ## 1531    1
    ## 1532    3
    ## 1534    3
    ## 1535    1
    ## 1536    4
    ## 1537    1
    ## 1538    2
    ## 1539    2
    ## 1540    2
    ## 1542    1
    ## 1543    3
    ## 1544    1
    ## 1545    2
    ## 1546    1
    ## 1547    4
    ## 1548    2
    ## 1549    1
    ## 1550    2
    ## 1551    4
    ## 1552    1
    ## 1553    2
    ## 1554    4
    ## 1555    1
    ## 1557    2
    ## 1558    4
    ## 1559    3
    ## 1560    2
    ## 1561    1
    ## 1562    3
    ## 1563    2
    ## 1564    1
    ## 1565    3
    ## 1566    2
    ## 1567    1
    ## 1568    3
    ## 1569    4
    ## 1571    1
    ## 1572    1
    ## 1573    4
    ## 1574    2
    ## 1575    2
    ## 1576    1
    ## 1577    2
    ## 1578    2
    ## 1579    1
    ## 1580    4
    ## 1581    3
    ## 1582    2
    ## 1583    2
    ## 1584    1
    ## 1586    4
    ## 1587    2
    ## 1589    2
    ## 1590    3
    ## 1591    2
    ## 1592    3
    ## 1593    2
    ## 1594    1
    ## 1595    4
    ## 1596    2
    ## 1597    1
    ## 1598    3
    ## 1599    4
    ## 1600    1
    ## 1602    4
    ## 1603    1
    ## 1605    3
    ## 1606    2
    ## 1607    3
    ## 1608    1
    ## 1609    1
    ## 1610    1
    ## 1611    4
    ## 1612    4
    ## 1613    1
    ## 1614    4
    ## 1615    2
    ## 1616    1
    ## 1617    1
    ## 1618    1
    ## 1619    1
    ## 1620    4
    ## 1621    2
    ## 1622    2
    ## 1623    1
    ## 1624    1
    ## 1625    4
    ## 1626    2
    ## 1627    1
    ## 1628    3
    ## 1629    1
    ## 1630    1
    ## 1631    2
    ## 1633    3
    ## 1634    1
    ## 1636    2
    ## 1637    1
    ## 1638    2
    ## 1642    2
    ## 1643    1
    ## 1644    2
    ## 1645    2
    ## 1646    3
    ## 1648    3
    ## 1649    2
    ## 1650    2
    ## 1651    1
    ## 1652    2
    ## 1653    1
    ## 1654    4
    ## 1655    1
    ## 1656    3
    ## 1657    2
    ## 1658    1
    ## 1659    1
    ## 1660    1
    ## 1661    3
    ## 1663    3
    ## 1664    3
    ## 1665    1
    ## 1666    2
    ## 1667    1
    ## 1668    2
    ## 1669    4
    ## 1670    4
    ## 1671    4
    ## 1672    1
    ## 1673    1
    ## 1675    1
    ## 1676    1
    ## 1677    2
    ## 1678    3
    ## 1679    2
    ## 1680    1
    ## 1681    2
    ## 1682    2
    ## 1683    2
    ## 1684    2
    ## 1685    2
    ## 1688    2
    ## 1689    1
    ## 1690    2
    ## 1691    3
    ## 1693    1
    ## 1695    2
    ## 1696    2
    ## 1697    2
    ## 1698    1
    ## 1699    1
    ## 1700    1
    ## 1701    2
    ## 1702    2
    ## 1703    2
    ## 1704    4
    ## 1705    2
    ## 1706    1
    ## 1707    4
    ## 1708    3
    ## 1709    2
    ## 1710    1
    ## 1711    1
    ## 1712    1
    ## 1714    2
    ## 1715    4
    ## 1716    4
    ## 1718    1
    ## 1720    3
    ## 1721    1
    ## 1722    1
    ## 1723    1
    ## 1724    3
    ## 1725    3
    ## 1727    1
    ## 1728    2
    ## 1729    2
    ## 1731    3
    ## 1732    1
    ## 1733    1
    ## 1734    1
    ## 1735    1
    ## 1736    1
    ## 1737    2
    ## 1738    3
    ## 1739    1
    ## 1740    4
    ## 1741    1
    ## 1742    2
    ## 1744    3
    ## 1745    1
    ## 1746    3
    ## 1748    1
    ## 1749    3
    ## 1750    1
    ## 1751    1
    ## 1752    1
    ## 1754    4
    ## 1755    1
    ## 1756    1
    ## 1757    1
    ## 1758    1
    ## 1759    3
    ## 1761    3
    ## 1762    4
    ## 1763    2
    ## 1764    1
    ## 1765    1
    ## 1766    1
    ## 1768    2
    ## 1770    2
    ## 1771    1
    ## 1772    4
    ## 1773    2
    ## 1774    1
    ## 1775    1
    ## 1776    2
    ## 1777    3
    ## 1778    1
    ## 1779    1
    ## 1780    4
    ## 1782    3
    ## 1783    2
    ## 1784    3
    ## 1785    3
    ## 1786    1
    ## 1787    3
    ## 1788    1
    ## 1789    2
    ## 1790    1
    ## 1792    1
    ## 1793    3
    ## 1794    1
    ## 1797    3
    ## 1799    2
    ## 1800    1
    ## 1801    3
    ## 1802    2
    ## 1803    1
    ## 1804    4
    ## 1806    3
    ## 1807    2
    ## 1808    2
    ## 1809    4
    ## 1810    3
    ## 1811    1
    ## 1812    2
    ## 1813    1
    ## 1814    2
    ## 1815    1
    ## 1816    1
    ## 1817    1
    ## 1818    1
    ## 1819    1
    ## 1821    4
    ## 1822    3
    ## 1825    1
    ## 1826    3
    ## 1827    3
    ## 1828    2
    ## 1829    1
    ## 1830    3
    ## 1831    1
    ## 1832    2
    ## 1833    1
    ## 1834    3
    ## 1836    1
    ## 1837    2
    ## 1838    1
    ## 1839    2
    ## 1840    1
    ## 1841    3
    ## 1842    2
    ## 1843    4
    ## 1844    2
    ## 1845    1
    ## 1846    3
    ## 1848    1
    ## 1849    4
    ## 1850    3
    ## 1851    1
    ## 1853    2
    ## 1854    3
    ## 1856    4
    ## 1857    3
    ## 1858    1
    ## 1859    1
    ## 1860    1
    ## 1861    2
    ## 1864    3
    ## 1865    1
    ## 1866    2
    ## 1867    1
    ## 1868    1
    ## 1869    4
    ## 1870    2
    ## 1871    3
    ## 1872    2
    ## 1873    2
    ## 1874    1
    ## 1875    4
    ## 1876    1
    ## 1877    1
    ## 1878    3
    ## 1881    1
    ## 1882    3
    ## 1883    2
    ## 1884    2
    ## 1885    3
    ## 1886    4
    ## 1887    3
    ## 1888    1
    ## 1889    1
    ## 1891    3
    ## 1892    2
    ## 1893    2
    ## 1894    1
    ## 1895    4
    ## 1896    2
    ## 1897    3
    ## 1899    4
    ## 1900    1
    ## 1901    4
    ## 1902    1
    ## 1903    3
    ## 1904    1
    ## 1906    1
    ## 1907    2
    ## 1909    2
    ## 1910    1
    ## 1911    2
    ## 1912    2
    ## 1913    2
    ## 1914    1
    ## 1916    3
    ## 1917    2
    ## 1918    3
    ## 1919    1
    ## 1920    3
    ## 1921    3
    ## 1922    2
    ## 1923    2
    ## 1924    1
    ## 1925    3
    ## 1926    1
    ## 1927    3
    ## 1928    2
    ## 1929    1
    ## 1930    3
    ## 1932    2
    ## 1933    1
    ## 1934    1
    ## 1935    1
    ## 1936    2
    ## 1937    2
    ## 1938    3
    ## 1939    1
    ## 1940    1
    ## 1941    3
    ## 1942    4
    ## 1944    1
    ## 1945    3
    ## 1946    2
    ## 1949    3
    ## 1950    2
    ## 1951    3
    ## 1952    1
    ## 1953    1
    ## 1955    3
    ## 1956    3
    ## 1958    2
    ## 1959    2
    ## 1960    2
    ## 1961    2
    ## 1962    4
    ## 1963    1
    ## 1964    2
    ## 1965    1
    ## 1966    1
    ## 1968    1
    ## 1969    1
    ## 1970    2
    ## 1972    1
    ## 1974    2
    ## 1975    1
    ## 1976    1
    ## 1977    1
    ## 1978    2
    ## 1979    2
    ## 1980    1
    ## 1981    4
    ## 1982    1
    ## 1983    1
    ## 1984    4
    ## 1985    1
    ## 1986    2
    ## 1987    2
    ## 1988    3
    ## 1989    2
    ## 1990    1
    ## 1991    1
    ## 1992    1
    ## 1993    1
    ## 1994    1
    ## 1995    3
    ## 1996    4
    ## 1997    1
    ## 1998    3
    ## 1999    2
    ## 2000    2
    ## 2001    1
    ## 2003    3
    ## 2005    1
    ## 2006    1
    ## 2007    3
    ## 2008    1
    ## 2009    1
    ## 2010    1
    ## 2011    2
    ## 2012    2
    ## 2014    2
    ## 2015    1
    ## 2017    2
    ## 2018    3
    ## 2019    2
    ## 2020    1
    ## 2021    2
    ## 2023    2
    ## 2024    1
    ## 2025    3
    ## 2026    1
    ## 2027    3
    ## 2028    1
    ## 2029    2
    ## 2030    1
    ## 2032    3
    ## 2033    1
    ## 2034    1
    ## 2035    2
    ## 2036    3
    ## 2037    3
    ## 2038    1
    ## 2039    1
    ## 2040    3
    ## 2041    1
    ## 2042    3
    ## 2043    1
    ## 2045    4
    ## 2046    1
    ## 2047    3
    ## 2048    3
    ## 2049    1
    ## 2050    1
    ## 2052    1
    ## 2053    3
    ## 2054    3
    ## 2055    3
    ## 2056    1
    ## 2057    3
    ## 2058    3
    ## 2059    1
    ## 2060    4
    ## 2061    3
    ## 2062    3
    ## 2063    1
    ## 2064    1
    ## 2065    2
    ## 2067    2
    ## 2068    4
    ## 2069    4
    ## 2072    1
    ## 2073    3
    ## 2075    3
    ## 2077    1
    ## 2078    3
    ## 2084    2
    ## 2085    1
    ## 2086    1
    ## 2087    2
    ## 2088    1
    ## 2089    1
    ## 2090    2
    ## 2091    2
    ## 2092    2
    ## 2093    1
    ## 2095    2
    ## 2096    4
    ## 2097    1
    ## 2100    1
    ## 2101    1
    ## 2103    1
    ## 2104    1
    ## 2106    1
    ## 2108    3
    ## 2109    2
    ## 2110    4
    ## 2112    4
    ## 2113    1
    ## 2115    1
    ## 2116    2
    ## 2117    4
    ## 2119    2
    ## 2120    2
    ## 2122    1
    ## 2123    1
    ## 2124    1
    ## 2126    1
    ## 2127    3
    ## 2128    4
    ## 2129    2
    ## 2130    1
    ## 2131    1
    ## 2132    3
    ## 2133    2
    ## 2134    1
    ## 2135    2
    ## 2136    1
    ## 2137    3
    ## 2139    2
    ## 2140    2
    ## 2141    1
    ## 2142    1
    ## 2143    4
    ## 2144    4
    ## 2145    3
    ## 2147    2
    ## 2148    3
    ## 2151    1
    ## 2152    1
    ## 2153    1
    ## 2154    1
    ## 2155    2
    ## 2156    1
    ## 2157    3
    ## 2158    4
    ## 2159    2
    ## 2160    1
    ## 2162    3
    ## 2163    1
    ## 2164    3
    ## 2166    1
    ## 2167    4
    ## 2170    2
    ## 2171    2
    ## 2172    3
    ## 2173    3
    ## 2174    2
    ## 2176    1
    ## 2177    2
    ## 2178    1
    ## 2179    4
    ## 2181    1
    ## 2183    1
    ## 2184    1
    ## 2185    1
    ## 2187    4
    ## 2188    3
    ## 2189    4
    ## 2190    4
    ## 2191    2
    ## 2192    1
    ## 2194    2
    ## 2195    1
    ## 2196    1
    ## 2197    1
    ## 2198    1
    ## 2200    3
    ## 2201    3
    ## 2202    3
    ## 2203    4
    ## 2204    2
    ## 2205    3
    ## 2208    1
    ## 2209    3
    ## 2214    3
    ## 2215    1
    ## 2216    2
    ## 2217    1
    ## 2218    1
    ## 2219    4
    ## 2220    1
    ## 2222    3
    ## 2224    2
    ## 2225    4
    ## 2226    1
    ## 2227    4
    ## 2228    1
    ## 2229    1
    ## 2230    3
    ## 2231    1
    ## 2232    1
    ## 2233    2
    ## 2234    1
    ## 2235    3
    ## 2236    1
    ## 2238    3
    ## 2239    2
    ## 2240    3
    ## 2242    2
    ## 2243    2
    ## 2244    4
    ## 2245    2
    ## 2246    2
    ## 2247    2
    ## 2248    2
    ## 2249    1
    ## 2250    3
    ## 2251    1
    ## 2252    2
    ## 2253    2
    ## 2254    2
    ## 2255    4
    ## 2256    2
    ## 2257    4
    ## 2258    3
    ## 2259    4
    ## 2260    1
    ## 2261    2
    ## 2262    1
    ## 2263    3
    ## 2266    1
    ## 2267    1
    ## 2268    2
    ## 2269    4
    ## 2270    4
    ## 2271    1
    ## 2272    2
    ## 2273    1
    ## 2274    4
    ## 2275    1
    ## 2276    1
    ## 2279    1
    ## 2280    1
    ## 2282    1
    ## 2283    2
    ## 2284    4
    ## 2285    1
    ## 2287    2
    ## 2289    1
    ## 2290    2
    ## 2291    2
    ## 2292    1
    ## 2293    1
    ## 2294    1
    ## 2295    2
    ## 2296    1
    ## 2297    2
    ## 2298    2
    ## 2299    3
    ## 2300    2
    ## 2301    1
    ## 2302    1
    ## 2303    1
    ## 2304    4
    ## 2305    3
    ## 2306    1
    ## 2307    2
    ## 2308    2
    ## 2309    2
    ## 2310    1
    ## 2311    1
    ## 2312    2
    ## 2313    2
    ## 2315    1
    ## 2316    1
    ## 2317    3
    ## 2318    1
    ## 2319    1
    ## 2320    4
    ## 2321    1
    ## 2323    1
    ## 2324    2
    ## 2326    1
    ## 2327    1
    ## 2328    4
    ## 2329    1
    ## 2330    2
    ## 2331    4
    ## 2332    1
    ## 2333    1
    ## 2334    2
    ## 2335    2
    ## 2336    3
    ## 2337    4
    ## 2338    4
    ## 2339    2
    ## 2340    1
    ## 2341    1
    ## 2342    3
    ## 2343    2
    ## 2344    1
    ## 2345    3
    ## 2346    4
    ## 2347    2
    ## 2348    2
    ## 2349    2
    ## 2350    4
    ## 2351    1
    ## 2352    4
    ## 2353    4
    ## 2354    1
    ## 2355    4
    ## 2356    3
    ## 2358    1
    ## 2359    4
    ## 2360    1
    ## 2362    4
    ## 2363    1
    ## 2364    3
    ## 2365    1
    ## 2366    3
    ## 2367    3
    ## 2368    3
    ## 2371    4
    ## 2372    2
    ## 2373    2
    ## 2375    1
    ## 2378    4
    ## 2379    2
    ## 2381    1
    ## 2382    2
    ## 2383    4
    ## 2384    4
    ## 2385    1
    ## 2386    3
    ## 2387    2
    ## 2388    4
    ## 2389    3
    ## 2390    2
    ## 2391    1
    ## 2392    2
    ## 2393    2
    ## 2394    3
    ## 2395    1
    ## 2396    2
    ## 2397    2
    ## 2398    2
    ## 2399    2
    ## 2400    1
    ## 2401    4
    ## 2402    4
    ## 2404    2
    ## 2405    1
    ## 2406    1
    ## 2407    4
    ## 2408    1
    ## 2409    2
    ## 2410    4
    ## 2411    1
    ## 2412    4
    ## 2413    2
    ## 2414    1
    ## 2415    2
    ## 2416    1
    ## 2417    1
    ## 2418    2
    ## 2419    2
    ## 2420    1
    ## 2421    1
    ## 2422    1
    ## 2423    1
    ## 2426    1
    ## 2427    1
    ## 2428    1
    ## 2429    1
    ## 2430    2
    ## 2431    3
    ## 2432    1
    ## 2433    2
    ## 2434    1
    ## 2435    4
    ## 2436    4
    ## 2437    1
    ## 2438    4
    ## 2439    2
    ## 2440    1
    ## 2441    1
    ## 2442    2
    ## 2443    1
    ## 2444    2
    ## 2445    4
    ## 2446    1
    ## 2447    2
    ## 2448    3
    ## 2449    2
    ## 2453    4
    ## 2454    4
    ## 2455    4
    ## 2456    1
    ## 2458    3
    ## 2459    4
    ## 2460    3
    ## 2461    2
    ## 2462    3
    ## 2463    1
    ## 2464    4
    ## 2465    1
    ## 2467    4
    ## 2468    2
    ## 2469    3
    ## 2471    1
    ## 2472    1
    ## 2473    2
    ## 2474    3
    ## 2475    1
    ## 2476    2
    ## 2477    1
    ## 2478    2
    ## 2479    2
    ## 2480    3
    ## 2482    2
    ## 2483    3
    ## 2484    2
    ## 2485    1
    ## 2486    1
    ## 2487    1
    ## 2489    3
    ## 2490    1
    ## 2491    1
    ## 2493    1
    ## 2494    1
    ## 2495    1
    ## 2496    2
    ## 2497    4
    ## 2498    1
    ## 2499    2
    ## 2500    4
    ## 2501    4
    ## 2502    1
    ## 2503    1
    ## 2504    3
    ## 2505    2
    ## 2506    1
    ## 2507    3
    ## 2508    1
    ## 2509    1
    ## 2510    3
    ## 2511    1
    ## 2512    3
    ## 2513    3
    ## 2514    3
    ## 2515    3
    ## 2516    2
    ## 2518    1
    ## 2519    2
    ## 2520    4
    ## 2521    1
    ## 2522    4
    ## 2523    1
    ## 2524    1
    ## 2525    1
    ## 2526    1
    ## 2527    2
    ## 2529    2
    ## 2530    2
    ## 2531    2
    ## 2532    1
    ## 2533    1
    ## 2534    3
    ## 2535    2
    ## 2537    2
    ## 2538    3
    ## 2539    1
    ## 2540    1
    ## 2541    2
    ## 2542    1
    ## 2543    3
    ## 2544    1
    ## 2545    4
    ## 2546    2
    ## 2547    1
    ## 2548    2
    ## 2551    3
    ## 2552    3
    ## 2553    2
    ## 2554    4
    ## 2555    2
    ## 2556    2
    ## 2557    1
    ## 2558    4
    ## 2559    3
    ## 2560    2
    ## 2561    3
    ## 2562    2
    ## 2563    2
    ## 2564    2
    ## 2565    2
    ## 2566    2
    ## 2567    2
    ## 2568    4
    ## 2569    1
    ## 2570    3
    ## 2571    1
    ## 2572    3
    ## 2573    3
    ## 2574    1
    ## 2575    2
    ## 2576    2
    ## 2577    2
    ## 2578    3
    ## 2579    2
    ## 2580    2
    ## 2581    1
    ## 2583    2
    ## 2584    1
    ## 2585    2
    ## 2586    2
    ## 2587    1
    ## 2588    1
    ## 2590    3
    ## 2591    4
    ## 2592    4
    ## 2593    2
    ## 2594    2
    ## 2595    1
    ## 2597    3
    ## 2598    1
    ## 2599    1
    ## 2600    1
    ## 2601    1
    ## 2602    3
    ## 2603    3
    ## 2604    2
    ## 2605    2
    ## 2606    2
    ## 2607    1
    ## 2608    1
    ## 2609    1
    ## 2610    2
    ## 2611    2
    ## 2612    4
    ## 2613    2
    ## 2614    4
    ## 2615    1
    ## 2616    3
    ## 2617    1
    ## 2618    1
    ## 2619    2
    ## 2620    1
    ## 2621    1
    ## 2622    1
    ## 2623    2
    ## 2624    1
    ## 2625    1
    ## 2626    1
    ## 2627    1
    ## 2628    2
    ## 2629    3
    ## 2630    2
    ## 2632    2
    ## 2633    2
    ## 2634    1
    ## 2635    2
    ## 2636    2
    ## 2637    1
    ## 2638    1
    ## 2639    2
    ## 2640    2
    ## 2641    4
    ## 2643    1
    ## 2644    4
    ## 2645    1
    ## 2646    3
    ## 2649    1
    ## 2650    2
    ## 2651    2
    ## 2652    2
    ## 2655    2
    ## 2656    1
    ## 2657    1
    ## 2658    1
    ## 2659    1
    ## 2660    2
    ## 2661    2
    ## 2662    3
    ## 2663    1
    ## 2664    3
    ## 2665    3
    ## 2667    1
    ## 2668    3
    ## 2669    4
    ## 2670    2
    ## 2671    1
    ## 2672    3
    ## 2673    2
    ## 2674    1
    ## 2675    4
    ## 2676    2
    ## 2677    3
    ## 2679    1
    ## 2680    1
    ## 2681    1
    ## 2682    2
    ## 2683    1
    ## 2684    4
    ## 2685    4
    ## 2686    2
    ## 2687    2
    ## 2688    2
    ## 2689    1
    ## 2690    2
    ## 2692    2
    ## 2693    1
    ## 2695    4
    ## 2696    1
    ## 2698    1
    ## 2699    3
    ## 2700    2
    ## 2701    2
    ## 2702    4
    ## 2703    1
    ## 2704    3
    ## 2705    2
    ## 2706    2
    ## 2708    1
    ## 2710    4
    ## 2711    4
    ## 2712    2
    ## 2714    2
    ## 2715    2
    ## 2716    2
    ## 2718    4
    ## 2719    1
    ## 2721    3
    ## 2722    1
    ## 2723    3
    ## 2724    1
    ## 2725    1
    ## 2726    3
    ## 2727    2
    ## 2728    1
    ## 2729    1
    ## 2730    3
    ## 2731    1
    ## 2732    1
    ## 2735    4
    ## 2736    1
    ## 2738    1
    ## 2739    3
    ## 2741    4
    ## 2742    2
    ## 2743    1
    ## 2744    2
    ## 2745    1
    ## 2746    2
    ## 2747    4
    ## 2748    2
    ## 2749    2
    ## 2750    3
    ## 2752    1
    ## 2754    1
    ## 2755    2
    ## 2756    2
    ## 2757    1
    ## 2758    2
    ## 2759    2
    ## 2760    2
    ## 2761    2
    ## 2762    1
    ## 2763    2
    ## 2764    1
    ## 2766    1
    ## 2768    4
    ## 2769    1
    ## 2770    1
    ## 2771    1
    ## 2773    1
    ## 2774    1
    ## 2775    2
    ## 2776    1
    ## 2777    1
    ## 2778    2
    ## 2779    2
    ## 2780    4
    ## 2781    1
    ## 2783    3
    ## 2784    1
    ## 2785    1
    ## 2787    2
    ## 2788    3
    ## 2789    3
    ## 2790    1
    ## 2791    1
    ## 2792    2
    ## 2793    1
    ## 2794    1
    ## 2795    3
    ## 2796    2
    ## 2798    1
    ## 2799    1
    ## 2800    2
    ## 2801    1
    ## 2802    1
    ## 2803    1
    ## 2804    2
    ## 2805    4
    ## 2806    2
    ## 2807    3
    ## 2808    1
    ## 2809    2
    ## 2810    2
    ## 2811    1
    ## 2812    1
    ## 2813    3
    ## 2814    1
    ## 2815    1
    ## 2816    3
    ## 2817    2
    ## 2819    1
    ## 2820    2
    ## 2821    2
    ## 2822    1
    ## 2823    1
    ## 2826    1
    ## 2827    2
    ## 2829    1
    ## 2830    1
    ## 2831    2
    ## 2832    4
    ## 2833    1
    ## 2834    1
    ## 2835    1
    ## 2836    2
    ## 2837    1
    ## 2838    1
    ## 2841    1
    ## 2842    3
    ## 2843    3
    ## 2844    1
    ## 2845    1
    ## 2846    4
    ## 2847    1
    ## 2848    1
    ## 2850    2
    ## 2851    3
    ## 2852    2
    ## 2853    4
    ## 2854    1
    ## 2855    1
    ## 2856    3
    ## 2859    2
    ## 2860    2
    ## 2861    1
    ## 2862    2
    ## 2863    3
    ## 2864    1
    ## 2865    2
    ## 2866    1
    ## 2867    2
    ## 2869    1
    ## 2870    2
    ## 2871    2
    ## 2872    1
    ## 2873    1
    ## 2874    2
    ## 2875    1
    ## 2876    3
    ## 2877    2
    ## 2879    1
    ## 2880    1
    ## 2881    2
    ## 2883    4
    ## 2884    1
    ## 2885    3
    ## 2887    1
    ## 2888    1
    ## 2889    1
    ## 2890    3
    ## 2891    4
    ## 2892    1
    ## 2894    4
    ## 2895    1
    ## 2896    2
    ## 2897    4
    ## 2898    2
    ## 2899    1
    ## 2900    1
    ## 2902    2
    ## 2904    2
    ## 2905    2
    ## 2906    1
    ## 2907    4
    ## 2908    2
    ## 2909    1
    ## 2910    2
    ## 2911    2
    ## 2912    3
    ## 2913    4
    ## 2915    2
    ## 2916    2
    ## 2917    4
    ## 2918    4
    ## 2919    1
    ## 2920    1
    ## 2921    2
    ## 2923    1
    ## 2924    1
    ## 2925    1
    ## 2926    4
    ## 2927    4
    ## 2928    1
    ## 2929    2
    ## 2931    1
    ## 2932    2
    ## 2933    1
    ## 2934    3
    ## 2936    1
    ## 2937    1
    ## 2938    2
    ## 2939    1
    ## 2940    3
    ## 2941    2
    ## 2942    2
    ## 2943    1
    ## 2944    2
    ## 2945    2
    ## 2946    2
    ## 2947    1
    ## 2948    1
    ## 2949    2
    ## 2950    1
    ## 2951    3
    ## 2952    4
    ## 2953    2
    ## 2954    1
    ## 2955    1
    ## 2957    2
    ## 2958    1
    ## 2959    1
    ## 2960    2
    ## 2961    2
    ## 2962    1
    ## 2963    3
    ## 2964    2
    ## 2965    3
    ## 2966    1
    ## 2968    1
    ## 2969    3
    ## 2970    1
    ## 2971    3
    ## 2972    1
    ## 2973    1
    ## 2974    1
    ## 2975    1
    ## 2976    2
    ## 2977    3
    ## 2980    1
    ## 2981    1
    ## 2982    2
    ## 2983    2
    ## 2984    2
    ## 2985    1
    ## 2986    4
    ## 2987    1
    ## 2988    4
    ## 2989    3
    ## 2990    1
    ## 2992    2
    ## 2993    3
    ## 2994    2
    ## 2995    3
    ## 2996    2
    ## 2997    1
    ## 2998    2
    ## 2999    4
    ## 3000    3
    ## 3001    2
    ## 3002    2
    ## 3003    1
    ## 3004    1
    ## 3006    3
    ## 3007    1
    ## 3008    2
    ## 3009    4
    ## 3010    2
    ## 3011    1
    ## 3012    3
    ## 3013    2
    ## 3014    3
    ## 3016    1
    ## 3017    4
    ## 3019    1
    ## 3020    2
    ## 3021    1
    ## 3022    4
    ## 3023    3
    ## 3024    1
    ## 3025    1
    ## 3026    3
    ## 3027    1
    ## 3029    3
    ## 3030    1
    ## 3031    4
    ## 3032    1
    ## 3033    2
    ## 3036    1
    ## 3039    1
    ## 3040    1
    ## 3041    2
    ## 3042    2
    ## 3043    2
    ## 3045    2
    ## 3046    2
    ## 3047    2
    ## 3048    1
    ## 3049    1
    ## 3050    1
    ## 3051    3
    ## 3052    4
    ## 3053    4
    ## 3054    1
    ## 3055    1
    ## 3056    1
    ## 3057    2
    ## 3059    1
    ## 3060    2
    ## 3063    1
    ## 3064    3
    ## 3065    3
    ## 3066    1
    ## 3067    4
    ## 3068    1
    ## 3069    1
    ## 3070    1
    ## 3071    3
    ## 3072    1
    ## 3073    2
    ## 3074    3
    ## 3076    1
    ## 3079    3
    ## 3081    2
    ## 3082    1
    ## 3084    3
    ## 3085    1
    ## 3086    1
    ## 3087    1
    ## 3088    1
    ## 3089    3
    ## 3091    1
    ## 3092    4
    ## 3093    1
    ## 3094    1
    ## 3095    1
    ## 3096    2
    ## 3097    3
    ## 3098    1
    ## 3099    1
    ## 3100    3
    ## 3102    1
    ## 3103    4
    ## 3104    1
    ## 3106    1
    ## 3107    2
    ## 3108    2
    ## 3109    2
    ## 3110    1
    ## 3111    1
    ## 3112    1
    ## 3113    2
    ## 3114    4
    ## 3115    1
    ## 3116    3
    ## 3117    1
    ## 3118    1
    ## 3119    1
    ## 3122    2
    ## 3123    4
    ## 3124    1
    ## 3125    4
    ## 3126    2
    ## 3127    1
    ## 3129    2
    ## 3130    1
    ## 3131    4
    ## 3132    1
    ## 3133    1
    ## 3134    2
    ## 3136    3
    ## 3137    2
    ## 3138    1
    ## 3139    1
    ## 3140    2
    ## 3142    3
    ## 3143    1
    ## 3144    1
    ## 3145    1
    ## 3147    4
    ## 3148    4
    ## 3149    2
    ## 3150    4
    ## 3151    1
    ## 3152    1
    ## 3153    2
    ## 3154    2
    ## 3155    2
    ## 3156    2
    ## 3157    2
    ## 3158    3
    ## 3159    2
    ## 3160    1
    ## 3161    4
    ## 3162    1
    ## 3163    2
    ## 3164    3
    ## 3165    4
    ## 3166    1
    ## 3167    1
    ## 3170    4
    ## 3171    2
    ## 3172    4
    ## 3173    1
    ## 3174    1
    ## 3175    3
    ## 3176    1
    ## 3177    1
    ## 3178    4
    ## 3179    4
    ## 3180    1
    ## 3183    1
    ## 3185    3
    ## 3186    3
    ## 3188    1
    ## 3189    3
    ## 3190    2
    ## 3192    1
    ## 3193    1
    ## 3194    2
    ## 3195    2
    ## 3196    1
    ## 3197    2
    ## 3198    2
    ## 3200    1
    ## 3201    4
    ## 3202    1
    ## 3203    2
    ## 3204    2
    ## 3205    3
    ## 3206    1
    ## 3207    1
    ## 3208    1
    ## 3209    2
    ## 3210    3
    ## 3211    2
    ## 3212    3
    ## 3213    1
    ## 3214    1
    ## 3215    1
    ## 3216    1
    ## 3217    1
    ## 3218    2
    ## 3219    4
    ## 3220    2
    ## 3221    1
    ## 3222    2
    ## 3223    1
    ## 3224    2
    ## 3225    1
    ## 3228    3
    ## 3229    2
    ## 3230    4
    ## 3231    1
    ## 3232    4
    ## 3233    1
    ## 3234    1
    ## 3236    2
    ## 3237    2
    ## 3238    2
    ## 3239    1
    ## 3240    1
    ## 3241    2
    ## 3243    3
    ## 3244    1
    ## 3245    2
    ## 3246    2
    ## 3249    4
    ## 3251    4
    ## 3252    1
    ## 3253    1
    ## 3254    3
    ## 3256    3
    ## 3257    3
    ## 3258    4
    ## 3259    2
    ## 3261    3
    ## 3262    3
    ## 3263    1
    ## 3264    3
    ## 3265    2
    ## 3266    3
    ## 3267    1
    ## 3268    2
    ## 3269    2
    ## 3270    2
    ## 3271    4
    ## 3272    1
    ## 3273    1
    ## 3274    2
    ## 3275    1
    ## 3276    1
    ## 3277    2
    ## 3278    1
    ## 3279    1
    ## 3280    1
    ## 3281    2
    ## 3282    2
    ## 3284    2
    ## 3287    3
    ## 3288    2
    ## 3289    2
    ## 3291    1
    ## 3292    2
    ## 3293    1
    ## 3294    3
    ## 3297    2
    ## 3298    4
    ## 3300    1
    ## 3301    1
    ## 3302    4
    ## 3303    3
    ## 3305    2
    ## 3306    3
    ## 3307    3
    ## 3308    1
    ## 3309    1
    ## 3310    3
    ## 3311    3
    ## 3312    2
    ## 3313    1
    ## 3314    2
    ## 3315    2
    ## 3316    4
    ## 3317    2
    ## 3318    1
    ## 3319    2
    ## 3320    3
    ## 3321    2
    ## 3322    1
    ## 3323    4
    ## 3324    1
    ## 3325    2
    ## 3326    2
    ## 3327    3
    ## 3328    1
    ## 3329    1
    ## 3330    1
    ## 3331    2
    ## 3334    4
    ## 3335    2
    ## 3336    1
    ## 3337    1
    ## 3338    4
    ## 3339    1
    ## 3340    1
    ## 3341    3
    ## 3342    2
    ## 3343    2
    ## 3344    1
    ## 3345    3
    ## 3346    1
    ## 3347    1
    ## 3348    1
    ## 3349    4
    ## 3350    1
    ## 3352    1
    ## 3353    1
    ## 3354    1
    ## 3355    1
    ## 3356    2
    ## 3357    1
    ## 3358    4
    ## 3359    1
    ## 3360    1
    ## 3361    1
    ## 3362    2
    ## 3363    1
    ## 3364    1
    ## 3366    2
    ## 3367    1
    ## 3368    1
    ## 3369    2
    ## 3370    4
    ## 3371    1
    ## 3372    1
    ## 3374    2
    ## 3375    1
    ## 3376    1
    ## 3377    2
    ## 3378    3
    ## 3379    4
    ## 3380    1
    ## 3382    3
    ## 3383    4
    ## 3384    1
    ## 3385    1
    ## 3386    1
    ## 3387    2
    ## 3388    2
    ## 3389    3
    ## 3390    2
    ## 3391    1
    ## 3392    1
    ## 3393    4
    ## 3394    1
    ## 3395    1
    ## 3396    1
    ## 3397    3
    ## 3398    3
    ## 3399    4
    ## 3400    3
    ## 3401    4
    ## 3402    3
    ## 3403    4
    ## 3404    4
    ## 3406    4
    ## 3407    4
    ## 3408    3
    ## 3409    3
    ## 3410    2
    ## 3411    1
    ## 3412    2
    ## 3413    1
    ## 3414    2
    ## 3415    3
    ## 3416    1
    ## 3417    3
    ## 3419    1
    ## 3420    1
    ## 3421    1
    ## 3422    2
    ## 3423    3
    ## 3424    1
    ## 3425    2
    ## 3426    1
    ## 3427    2
    ## 3428    1
    ## 3429    1
    ## 3430    3
    ## 3431    2
    ## 3433    2
    ## 3434    3
    ## 3435    2
    ## 3437    1
    ## 3438    2
    ## 3439    2
    ## 3440    1
    ## 3444    3
    ## 3445    2
    ## 3446    1
    ## 3447    1
    ## 3449    1
    ## 3450    3
    ## 3451    1
    ## 3452    1
    ## 3453    2
    ## 3454    1
    ## 3455    1
    ## 3456    2
    ## 3457    1
    ## 3458    3
    ## 3459    1
    ## 3461    2
    ## 3463    1
    ## 3464    1
    ## 3466    1
    ## 3467    2
    ## 3469    1
    ## 3470    3
    ## 3471    2
    ## 3472    4
    ## 3473    2
    ## 3474    1
    ## 3475    3
    ## 3476    2
    ## 3477    1
    ## 3478    3
    ## 3479    4
    ## 3480    2
    ## 3481    3
    ## 3482    1
    ## 3484    1
    ## 3486    1
    ## 3487    1
    ## 3488    2
    ## 3489    1
    ## 3490    2
    ## 3491    1
    ## 3492    2
    ## 3494    1
    ## 3495    2
    ## 3496    3
    ## 3497    2
    ## 3498    1
    ## 3499    2
    ## 3500    1
    ## 3501    1
    ## 3502    2
    ## 3503    4
    ## 3504    1
    ## 3505    2
    ## 3506    4
    ## 3507    3
    ## 3508    1
    ## 3509    1
    ## 3510    3
    ## 3511    4
    ## 3512    1
    ## 3513    1
    ## 3514    2
    ## 3515    2
    ## 3516    3
    ## 3517    1
    ## 3519    1
    ## 3520    2
    ## 3521    3
    ## 3522    1
    ## 3523    4
    ## 3524    1
    ## 3525    1
    ## 3526    3
    ## 3528    4
    ## 3530    1
    ## 3531    1
    ## 3532    1
    ## 3533    3
    ## 3535    2
    ## 3536    1
    ## 3538    2
    ## 3539    1
    ## 3541    2
    ## 3542    1
    ## 3543    1
    ## 3544    1
    ## 3545    1
    ## 3546    3
    ## 3547    1
    ## 3548    4
    ## 3549    1
    ## 3550    2
    ## 3551    1
    ## 3552    2
    ## 3553    4
    ## 3554    1
    ## 3555    4
    ## 3556    1
    ## 3557    2
    ## 3558    1
    ## 3559    4
    ## 3560    2
    ## 3561    1
    ## 3562    1
    ## 3564    1
    ## 3565    4
    ## 3566    1
    ## 3567    2
    ## 3568    3
    ## 3569    1
    ## 3570    2
    ## 3571    2
    ## 3572    4
    ## 3573    1
    ## 3574    1
    ## 3575    1
    ## 3576    1
    ## 3577    4
    ## 3578    2
    ## 3579    1
    ## 3580    1
    ## 3581    3
    ## 3582    2
    ## 3583    4
    ## 3584    1
    ## 3585    1
    ## 3586    1
    ## 3587    1
    ## 3588    2
    ## 3589    1
    ## 3590    4
    ## 3591    1
    ## 3592    2
    ## 3593    1
    ## 3594    1
    ## 3595    4
    ## 3596    1
    ## 3597    3
    ## 3598    1
    ## 3599    1
    ## 3600    2
    ## 3601    1
    ## 3602    1
    ## 3603    4
    ## 3604    2
    ## 3605    4
    ## 3606    1
    ## 3607    2
    ## 3608    3
    ## 3609    2
    ## 3610    2
    ## 3611    4
    ## 3614    4
    ## 3615    3
    ## 3616    2
    ## 3617    1
    ## 3618    2
    ## 3619    1
    ## 3620    1
    ## 3621    1
    ## 3623    1
    ## 3624    3
    ## 3625    4
    ## 3626    1
    ## 3628    3
    ## 3629    3
    ## 3630    4
    ## 3631    2
    ## 3632    2
    ## 3633    1
    ## 3634    2
    ## 3635    2
    ## 3636    1
    ## 3637    4
    ## 3639    1
    ## 3640    1
    ## 3641    2
    ## 3642    2
    ## 3643    1
    ## 3644    1
    ## 3646    2
    ## 3647    2
    ## 3649    2
    ## 3650    4
    ## 3651    1
    ## 3652    1
    ## 3653    1
    ## 3654    3
    ## 3655    1
    ## 3656    2
    ## 3657    1
    ## 3658    1
    ## 3659    1
    ## 3660    2
    ## 3661    2
    ## 3662    1
    ## 3663    2
    ## 3664    1
    ## 3665    2
    ## 3666    3
    ## 3667    2
    ## 3668    1
    ## 3669    1
    ## 3670    4
    ## 3671    3
    ## 3672    3
    ## 3673    3
    ## 3675    1
    ## 3676    3
    ## 3677    1
    ## 3678    2
    ## 3679    3
    ## 3680    3
    ## 3681    1
    ## 3682    1
    ## 3683    3
    ## 3684    1
    ## 3685    1
    ## 3687    1
    ## 3688    1
    ## 3689    4
    ## 3690    3
    ## 3691    3
    ## 3692    1
    ## 3693    3
    ## 3694    3
    ## 3696    1
    ## 3697    2
    ## 3698    1
    ## 3699    4
    ## 3700    2
    ## 3701    1
    ## 3702    1
    ## 3703    4
    ## 3704    2
    ## 3705    1
    ## 3706    1
    ## 3707    2
    ## 3708    1
    ## 3709    2
    ## 3710    1
    ## 3711    2
    ## 3712    1
    ## 3713    4
    ## 3714    1
    ## 3715    3
    ## 3716    1
    ## 3717    1
    ## 3718    2
    ## 3719    2
    ## 3721    1
    ## 3723    1
    ## 3724    1
    ## 3725    2
    ## 3726    4
    ## 3727    3
    ## 3728    2
    ## 3729    2
    ## 3730    2
    ## 3731    3
    ## 3732    2
    ## 3734    2
    ## 3736    1
    ## 3737    2
    ## 3739    1
    ## 3740    3
    ## 3742    1
    ## 3743    3
    ## 3744    1
    ## 3745    1
    ## 3746    3
    ## 3747    1
    ## 3748    2
    ## 3749    1
    ## 3750    1
    ## 3752    1
    ## 3753    1
    ## 3756    1
    ## 3757    1
    ## 3758    1
    ## 3759    1
    ## 3760    2
    ## 3762    1
    ## 3763    2
    ## 3764    2
    ## 3765    1
    ## 3766    2
    ## 3767    1
    ## 3768    1
    ## 3770    2
    ## 3771    3
    ## 3772    4
    ## 3774    1
    ## 3775    3
    ## 3776    2
    ## 3777    2
    ## 3778    2
    ## 3780    1
    ## 3781    3
    ## 3782    2
    ## 3783    3
    ## 3785    1
    ## 3786    2
    ## 3787    1
    ## 3788    1
    ## 3789    3
    ## 3791    3
    ## 3792    1
    ## 3793    4
    ## 3794    2
    ## 3795    2
    ## 3797    1
    ## 3798    1
    ## 3799    4
    ## 3800    2
    ## 3801    1
    ## 3802    4
    ## 3803    2
    ## 3804    3
    ## 3805    1
    ## 3806    2
    ## 3807    1
    ## 3808    1
    ## 3809    4
    ## 3810    3
    ## 3812    1
    ## 3813    1
    ## 3815    1
    ## 3817    3
    ## 3818    2
    ## 3819    3
    ## 3820    1
    ## 3821    1
    ## 3823    2
    ## 3824    1
    ## 3825    2
    ## 3826    2
    ## 3829    3
    ## 3830    4
    ## 3831    1
    ## 3832    4
    ## 3833    2
    ## 3834    2
    ## 3835    1
    ## 3836    1
    ## 3838    4
    ## 3839    3
    ## 3841    1
    ## 3842    1
    ## 3844    3
    ## 3845    3
    ## 3846    2
    ## 3847    4
    ## 3848    1
    ## 3849    2
    ## 3850    2
    ## 3852    1
    ## 3853    1
    ## 3854    1
    ## 3856    2
    ## 3857    2
    ## 3858    2
    ## 3859    3
    ## 3860    1
    ## 3861    4
    ## 3862    2
    ## 3863    3
    ## 3864    1
    ## 3865    3
    ## 3866    2
    ## 3867    2
    ## 3868    4
    ## 3869    2
    ## 3870    4
    ## 3871    1
    ## 3872    3
    ## 3873    2
    ## 3874    1
    ## 3875    1
    ## 3876    4
    ## 3877    2
    ## 3878    1
    ## 3879    1
    ## 3880    3
    ## 3881    2
    ## 3884    2
    ## 3885    2
    ## 3886    2
    ## 3887    3
    ## 3888    2
    ## 3889    1
    ## 3890    4
    ## 3891    1
    ## 3892    3
    ## 3893    1
    ## 3894    3
    ## 3895    1
    ## 3896    3
    ## 3897    3
    ## 3898    3
    ## 3899    3
    ## 3900    1
    ## 3901    1
    ## 3902    1
    ## 3903    4
    ## 3904    4
    ## 3905    2
    ## 3906    1
    ## 3907    1
    ## 3909    1
    ## 3910    1
    ## 3911    2
    ## 3912    3
    ## 3913    3
    ## 3915    1
    ## 3916    1
    ## 3918    1
    ## 3919    2
    ## 3920    1
    ## 3921    1
    ## 3922    2
    ## 3924    2
    ## 3925    3
    ## 3926    1
    ## 3928    3
    ## 3929    1
    ## 3930    4
    ## 3931    3
    ## 3932    1
    ## 3933    1
    ## 3934    2
    ## 3937    2
    ## 3938    1
    ## 3939    2
    ## 3940    2
    ## 3941    2
    ## 3942    1
    ## 3944    2
    ## 3946    2
    ## 3948    2
    ## 3949    1
    ## 3950    1
    ## 3951    4
    ## 3953    1
    ## 3954    1
    ## 3955    4
    ## 3956    2
    ## 3957    1
    ## 3958    1
    ## 3960    1
    ## 3961    3
    ## 3962    4
    ## 3964    3
    ## 3965    3
    ## 3966    2
    ## 3967    1
    ## 3968    1
    ## 3969    2
    ## 3970    2
    ## 3971    3
    ## 3972    1
    ## 3973    3
    ## 3974    3
    ## 3975    2
    ## 3976    4
    ## 3977    1
    ## 3978    2
    ## 3979    4
    ## 3980    1
    ## 3981    3
    ## 3982    1
    ## 3983    2
    ## 3984    1
    ## 3985    3
    ## 3986    4
    ## 3987    1
    ## 3989    2
    ## 3990    2
    ## 3991    2
    ## 3992    1
    ## 3993    2
    ## 3994    3
    ## 3995    1
    ## 3996    3
    ## 3997    2
    ## 3998    1
    ## 3999    3
    ## 4000    1
    ## 4001    4
    ## 4002    1
    ## 4003    1
    ## 4004    2
    ## 4005    1
    ## 4006    3
    ## 4007    1
    ## 4008    2
    ## 4009    3
    ## 4010    1
    ## 4011    1
    ## 4012    3
    ## 4013    2
    ## 4014    1
    ## 4015    1
    ## 4016    3
    ## 4017    3
    ## 4018    1
    ## 4019    1
    ## 4020    4
    ## 4021    2
    ## 4022    2
    ## 4023    2
    ## 4024    1
    ## 4025    4
    ## 4026    1
    ## 4027    1
    ## 4028    2
    ## 4029    1
    ## 4030    1
    ## 4031    2
    ## 4032    1
    ## 4033    4
    ## 4034    3
    ## 4035    3
    ## 4036    4
    ## 4038    3
    ## 4039    3
    ## 4040    1
    ## 4041    3
    ## 4043    2
    ## 4044    1
    ## 4045    1
    ## 4048    4
    ## 4049    1
    ## 4050    1
    ## 4051    3
    ## 4053    4
    ## 4054    2
    ## 4055    1
    ## 4056    1
    ## 4058    4
    ## 4059    4
    ## 4060    4
    ## 4061    1
    ## 4062    2
    ## 4063    1
    ## 4064    3
    ## 4065    4
    ## 4066    1
    ## 4067    2
    ## 4069    1
    ## 4070    4
    ## 4071    1
    ## 4072    1
    ## 4073    1
    ## 4074    2
    ## 4075    1
    ## 4076    4
    ## 4077    4
    ## 4079    1
    ## 4080    1
    ## 4081    1
    ## 4083    1
    ## 4084    1
    ## 4085    3
    ## 4086    2
    ## 4087    3
    ## 4088    1
    ## 4089    2
    ## 4090    2
    ## 4091    4
    ## 4093    2
    ## 4094    1
    ## 4095    1
    ## 4096    2
    ## 4097    1
    ## 4098    1
    ## 4099    1
    ## 4100    3
    ## 4101    2
    ## 4102    1
    ## 4103    2
    ## 4105    2
    ## 4106    1
    ## 4107    2
    ## 4108    2
    ## 4109    3
    ## 4110    2
    ## 4111    3
    ## 4112    4
    ## 4113    1
    ## 4114    1
    ## 4115    1
    ## 4116    2
    ## 4117    3
    ## 4118    2
    ## 4119    1
    ## 4120    4
    ## 4121    3
    ## 4123    2
    ## 4124    1
    ## 4125    2
    ## 4127    1
    ## 4128    2
    ## 4129    1
    ## 4130    1
    ## 4131    1
    ## 4132    2
    ## 4133    1
    ## 4134    4
    ## 4135    1
    ## 4136    4
    ## 4138    1
    ## 4139    2
    ## 4142    1
    ## 4143    1
    ## 4144    3
    ## 4145    2
    ## 4146    1
    ## 4147    2
    ## 4148    2
    ## 4149    1
    ## 4150    2
    ## 4151    3
    ## 4152    2
    ## 4153    1
    ## 4154    1
    ## 4155    2
    ## 4156    4
    ## 4157    1
    ## 4158    3
    ## 4160    2
    ## 4161    4
    ## 4162    3
    ## 4163    2
    ## 4164    1
    ## 4165    4
    ## 4166    4
    ## 4167    2
    ## 4171    2
    ## 4172    4
    ## 4173    1
    ## 4174    1
    ## 4175    2
    ## 4176    2
    ## 4177    1
    ## 4178    1
    ## 4179    1
    ## 4180    3
    ## 4181    4
    ## 4182    1
    ## 4183    1
    ## 4184    3
    ## 4186    1
    ## 4187    4
    ## 4188    4
    ## 4189    2
    ## 4191    1
    ## 4192    4
    ## 4193    1
    ## 4195    4
    ## 4196    1
    ## 4197    3
    ## 4198    2
    ## 4200    2
    ## 4201    1
    ## 4202    1
    ## 4203    1
    ## 4204    3
    ## 4205    3
    ## 4206    1
    ## 4207    3
    ## 4208    1
    ## 4209    2
    ## 4211    1
    ## 4212    3
    ## 4213    2
    ## 4214    1
    ## 4215    1
    ## 4216    2
    ## 4217    2
    ## 4218    2
    ## 4219    3
    ## 4220    3
    ## 4221    1
    ## 4222    3
    ## 4223    2
    ## 4224    1
    ## 4225    1
    ## 4226    1
    ## 4227    3
    ## 4228    2
    ## 4229    1
    ## 4230    2
    ## 4232    2
    ## 4233    2
    ## 4234    1
    ## 4235    1
    ## 4236    3
    ## 4237    2
    ## 4238    2
    ## 4239    2
    ## 4240    2
    ## 4241    1
    ## 4242    2
    ## 4243    4
    ## 4245    3
    ## 4246    1
    ## 4247    2
    ## 4248    4
    ## 4249    2
    ## 4250    2
    ## 4251    3
    ## 4252    3
    ## 4253    3
    ## 4255    1
    ## 4257    2
    ## 4258    2
    ## 4259    1
    ## 4260    3
    ## 4261    4
    ## 4262    1
    ## 4263    1
    ## 4264    2
    ## 4265    1
    ## 4266    1
    ## 4267    4
    ## 4268    1
    ## 4269    3
    ## 4271    2
    ## 4272    1
    ## 4273    1
    ## 4275    1
    ## 4276    1
    ## 4277    1
    ## 4278    1
    ## 4279    4
    ## 4280    2
    ## 4281    2
    ## 4282    1
    ## 4283    2
    ## 4284    1
    ## 4285    2
    ## 4287    2
    ## 4288    1
    ## 4289    2
    ## 4290    1
    ## 4291    2
    ## 4292    2
    ## 4293    2
    ## 4294    4
    ## 4295    2
    ## 4296    2
    ## 4297    2
    ## 4298    1
    ## 4299    1
    ## 4300    3
    ## 4301    3
    ## 4302    3
    ## 4303    4
    ## 4304    2
    ## 4305    3
    ## 4306    1
    ## 4307    2
    ## 4308    1
    ## 4309    1
    ## 4310    1
    ## 4313    3
    ## 4314    1
    ## 4315    1
    ## 4316    3
    ## 4317    1
    ## 4318    3
    ## 4319    1
    ## 4320    1
    ## 4321    2
    ## 4322    1
    ## 4323    3
    ## 4324    1
    ## 4325    1
    ## 4326    2
    ## 4327    2
    ## 4330    4
    ## 4331    3
    ## 4332    4
    ## 4333    2
    ## 4334    1
    ## 4335    2
    ## 4336    1
    ## 4337    3
    ## 4338    1
    ## 4339    1
    ## 4340    2
    ## 4341    1
    ## 4342    2
    ## 4343    1
    ## 4345    2
    ## 4346    4
    ## 4347    4
    ## 4348    4
    ## 4349    2
    ## 4350    4
    ## 4352    3
    ## 4354    4
    ## 4355    2
    ## 4357    1
    ## 4358    4
    ## 4359    1
    ## 4360    1
    ## 4361    1
    ## 4363    2
    ## 4364    4
    ## 4365    2
    ## 4366    2
    ## 4367    1
    ## 4368    1
    ## 4369    4
    ## 4370    2
    ## 4371    3
    ## 4372    1
    ## 4373    2
    ## 4374    2
    ## 4375    1
    ## 4376    1
    ## 4377    1
    ## 4378    1
    ## 4380    3
    ## 4381    2
    ## 4382    3
    ## 4383    2
    ## 4384    1
    ## 4385    3
    ## 4386    1
    ## 4387    1
    ## 4388    1
    ## 4389    2
    ## 4390    1
    ## 4391    1
    ## 4392    1
    ## 4393    1
    ## 4394    2
    ## 4395    2
    ## 4396    2
    ## 4397    1
    ## 4398    2
    ## 4399    1
    ## 4400    2
    ## 4401    4
    ## 4403    2
    ## 4404    1
    ## 4405    2
    ## 4406    2
    ## 4407    2
    ## 4408    3
    ## 4409    1
    ## 4410    3
    ## 4411    1
    ## 4412    1
    ## 4413    1
    ## 4414    1
    ## 4415    1
    ## 4416    3
    ## 4417    4
    ## 4418    2
    ## 4419    4
    ## 4420    1
    ## 4421    1
    ## 4422    4
    ## 4423    1
    ## 4426    3
    ## 4427    1
    ## 4428    1
    ## 4429    3
    ## 4432    2
    ## 4433    3
    ## 4434    3

You can select a range of columns using the `:` operator, or use the
`contains()` function to find all columns containing a specific word or
combination of letters.

``` r
fhsBP <- select(fhs, contains("BP"))
fhsBP
```

    ##      SYSBP DIABP BPMEDS
    ## 1    106.0  70.0     No
    ## 2    121.0  81.0     No
    ## 3    127.5  80.0     No
    ## 4    150.0  95.0     No
    ## 5    130.0  84.0     No
    ## 6    180.0 110.0     No
    ## 7    138.0  71.0     No
    ## 8    100.0  71.0     No
    ## 9    141.5  89.0     No
    ## 10   162.0 107.0     No
    ## 11   133.0  76.0     No
    ## 12   131.0  88.0     No
    ## 13   142.0  94.0     No
    ## 14   124.0  88.0    Yes
    ## 16   140.0  90.0     No
    ## 17   138.0  90.0     No
    ## 18   112.0  78.0     No
    ## 19   122.0  84.5     No
    ## 20   139.0  88.0     No
    ## 21   108.0  70.5     No
    ## 23   148.0  78.0     No
    ## 24   132.0  82.0     No
    ## 25   137.5  90.0     No
    ## 26   102.0  68.0     No
    ## 28   132.0  91.0     No
    ## 29   182.0 121.0     No
    ## 30   130.0  88.0     No
    ## 31   102.0  68.0     No
    ## 32   115.0  85.5     No
    ## 33   125.0  90.0     No
    ## 34   150.0  85.0     No
    ## 36   125.0  74.5     No
    ## 37   147.0  74.0     No
    ## 38   124.5  92.5     No
    ## 40   160.0  98.0     No
    ## 41   153.0 101.0     No
    ## 42   111.0  73.0     No
    ## 43   116.5  80.0     No
    ## 44   122.0  78.0     No
    ## 46   132.0  83.5     No
    ## 47   206.0  92.0    Yes
    ## 48    96.0  63.0     No
    ## 49   179.5 114.0     No
    ## 50   119.0  77.5     No
    ## 51   116.0  69.0     No
    ## 53   156.5  92.5     No
    ## 54   112.0  66.0     No
    ## 55   130.0  78.0     No
    ## 56   145.0  82.5     No
    ## 58   116.0  71.0     No
    ## 59   114.0  76.0     No
    ## 60   143.5  81.0     No
    ## 61   115.0  69.0     No
    ## 62   158.0 102.0     No
    ## 63   121.0  79.0     No
    ## 64   157.0  89.0     No
    ## 65   123.5  75.0     No
    ## 66   126.5  80.0     No
    ## 67   136.0  84.0     No
    ## 68   154.0  87.0     No
    ## 69   190.0  99.0     No
    ## 70   107.0  73.0     No
    ## 71   112.5  60.0     No
    ## 72   110.0  67.5     No
    ## 74   116.0  72.5     No
    ## 76   150.0 106.0     No
    ## 77   110.0  76.0     No
    ## 78   138.5  85.0     No
    ## 79   155.0  85.0     No
    ## 81   151.0 101.0     No
    ## 82   152.0  90.0     No
    ## 83   179.0  94.0     No
    ## 84   155.0 110.0     No
    ## 85   138.0  86.5     No
    ## 86   124.0  78.0     No
    ## 87   182.0 101.0     No
    ## 88   113.0  72.5     No
    ## 89   138.0  82.0     No
    ## 90   200.0 104.0     No
    ## 91   124.0  86.0     No
    ## 92   117.5  80.0     No
    ## 93   121.0  61.5     No
    ## 94   138.0  89.0     No
    ## 95   132.5  80.0     No
    ## 96   102.0  71.5     No
    ## 97   154.0  80.0     No
    ## 98   136.0  96.0     No
    ## 99   126.0  79.0     No
    ## 100  123.0  76.5     No
    ## 101  134.0  80.0     No
    ## 103  180.0  90.0     No
    ## 104  121.0  75.0     No
    ## 105  132.5  87.0     No
    ## 106  141.0  84.0     No
    ## 107  110.0  64.0     No
    ## 108  141.0  82.0     No
    ## 109  145.0  77.0     No
    ## 110  100.0  63.0     No
    ## 111  135.0  82.0     No
    ## 112  115.0  79.0     No
    ## 113  138.0  90.0     No
    ## 114  187.0  88.0     No
    ## 115  123.0  69.0     No
    ## 116  139.0  80.0     No
    ## 118  138.0  88.5     No
    ## 119  158.0 105.0     No
    ## 121  127.0  82.0     No
    ## 122  160.5  96.0     No
    ## 123  100.0  68.0     No
    ## 124  112.0  85.5     No
    ## 125  157.0  88.0     No
    ## 126  105.0  75.0     No
    ## 127  112.0  70.0     No
    ## 128  109.0  71.0     No
    ## 129  135.0  97.0     No
    ## 130  158.0  86.0     No
    ## 131  112.5  80.0     No
    ## 132  128.0  77.0     No
    ## 133  118.0  77.0     No
    ## 134  115.0  77.5     No
    ## 135  112.0  73.0     No
    ## 136  106.0  67.5     No
    ## 138  134.0  89.0     No
    ## 139  128.0  64.0     No
    ## 140  121.0  84.0     No
    ## 141  131.0  81.0     No
    ## 142  131.0  81.0     No
    ## 143  124.0  70.0     No
    ## 144  154.0 100.0     No
    ## 145  127.5  81.5     No
    ## 146  138.0  82.0     No
    ## 147  100.0  72.5     No
    ## 148  117.5  72.5     No
    ## 149  151.0  95.0     No
    ## 150  136.0  94.0     No
    ## 151  149.0 100.0     No
    ## 152  141.0  86.0     No
    ## 153  153.0  98.0     No
    ## 154  123.0  81.0     No
    ## 155  180.5 106.5     No
    ## 156  130.0  77.0     No
    ## 157  142.0  93.0     No
    ## 158  136.5  76.0     No
    ## 159  212.0 104.0     No
    ## 161  134.0  82.0     No
    ## 162  153.0  80.5     No
    ## 163  150.0  77.5     No
    ## 164  191.0 124.5    Yes
    ## 165  108.0  81.0     No
    ## 166  121.5  73.0     No
    ## 167  139.0  86.0    Yes
    ## 168  180.0  97.5     No
    ## 169  118.0  86.5     No
    ## 170  102.0  61.0     No
    ## 171  117.5  80.0     No
    ## 172  173.0  89.0     No
    ## 173  109.0  78.0     No
    ## 174  110.0  78.0     No
    ## 175  134.0  83.0     No
    ## 176  109.0  77.0     No
    ## 177  144.0  78.0     No
    ## 178  129.5  83.0     No
    ## 179  122.0  80.5     No
    ## 180  110.0  74.0     No
    ## 181  111.0  64.0     No
    ## 182  124.0  76.0     No
    ## 183  158.0  86.5     No
    ## 184  151.0 101.0     No
    ## 185  182.0 106.0     No
    ## 186  122.0  83.0     No
    ## 187  117.5  71.0     No
    ## 188  110.0  67.0     No
    ## 189  112.0  60.0     No
    ## 190  128.0  82.0     No
    ## 191  118.0  76.0     No
    ## 192  128.0  77.0     No
    ## 194  158.0  98.0     No
    ## 195  117.0  74.5     No
    ## 196  138.0  97.0     No
    ## 197  125.0  86.5     No
    ## 198  144.5  83.5     No
    ## 199  133.0  92.0     No
    ## 200  109.0  70.0     No
    ## 201  127.0  81.0     No
    ## 203  128.0  87.0     No
    ## 204  106.0  71.0     No
    ## 205  117.0  76.0     No
    ## 206  157.0  82.5     No
    ## 207  170.0  86.0     No
    ## 208  137.0  84.0     No
    ## 209  111.0  68.0     No
    ## 210  110.0  67.5     No
    ## 211   94.0  66.5     No
    ## 213  135.0  91.0     No
    ## 215  144.0  78.0     No
    ## 216  155.0  86.0     No
    ## 217  117.5  72.5     No
    ## 218  128.0  68.0     No
    ## 219  125.0  65.0     No
    ## 220  124.0  84.0     No
    ## 222  166.0  98.0     No
    ## 224  126.5  93.0     No
    ## 227  114.0  85.5     No
    ## 228  123.0  72.0     No
    ## 229   96.0  70.0     No
    ## 230  117.5  77.5     No
    ## 231  124.0  66.0     No
    ## 232  130.0  70.0     No
    ## 233  106.0  82.0     No
    ## 234  160.0  85.0     No
    ## 235  134.0  60.0     No
    ## 236  142.0  89.0     No
    ## 237  135.0  92.0     No
    ## 238  174.0  92.5    Yes
    ## 239  177.5  98.0     No
    ## 240  131.0  93.0     No
    ## 241  136.5  99.5     No
    ## 242  128.0  91.0     No
    ## 243  129.0  70.0     No
    ## 244  112.5  79.0     No
    ## 245  159.0  90.0     No
    ## 246  128.0  74.0     No
    ## 247  108.0  74.0     No
    ## 248  155.0  85.0     No
    ## 249  124.0  85.0     No
    ## 250  122.0  84.0     No
    ## 251  129.0  88.0     No
    ## 252  139.0  88.0     No
    ## 253  148.0  90.0     No
    ## 254  134.0  81.0     No
    ## 255  178.0 127.0     No
    ## 256  124.0  76.0     No
    ## 257  130.5  82.0     No
    ## 261  107.5  73.0     No
    ## 262  200.0 122.5     No
    ## 264  144.0  84.0     No
    ## 265  109.0  70.0     No
    ## 266  124.0  74.0     No
    ## 267  125.0  80.0     No
    ## 268  189.0 110.0    Yes
    ## 269  168.0  98.0     No
    ## 270  139.0  90.0     No
    ## 271  108.0  72.5     No
    ## 272  197.5 105.0     No
    ## 273  146.0  89.0     No
    ## 274  137.0  91.0     No
    ## 275  150.0  81.0     No
    ## 277  134.0  97.0     No
    ## 278  107.5  72.5     No
    ## 279  156.5  86.0     No
    ## 280  140.0  93.0     No
    ## 281  134.0  98.0     No
    ## 282  125.0  85.0     No
    ## 283  139.0  98.0     No
    ## 284  174.0  84.5     No
    ## 285  129.0  83.0     No
    ## 287  136.0  87.0     No
    ## 288  127.5  86.0     No
    ## 289  129.0  80.0     No
    ## 290  147.0  92.0     No
    ## 291  127.5  91.0     No
    ## 293  132.0  78.0     No
    ## 294  123.0  79.0     No
    ## 296  100.0  74.5     No
    ## 297  142.0  85.0     No
    ## 298   98.0  64.0     No
    ## 299  105.0  57.0     No
    ## 300  105.0  69.0     No
    ## 301  111.0  68.0     No
    ## 302  126.0  72.0     No
    ## 303  131.5  76.0     No
    ## 304  116.0  72.0     No
    ## 305  125.0  76.0     No
    ## 308  117.5  80.0     No
    ## 310  101.0  63.0     No
    ## 311  130.0  80.0     No
    ## 312  130.0  80.0    Yes
    ## 313  121.0  86.0     No
    ## 316  158.5  81.0     No
    ## 317   97.0  65.0     No
    ## 320  101.0  61.0     No
    ## 321   97.0  64.0     No
    ## 323  114.0  77.0     No
    ## 324  151.5  88.0     No
    ## 325  110.0  65.0     No
    ## 326  140.0  94.0     No
    ## 329  168.0  98.0     No
    ## 330  157.0 101.0     No
    ## 331  125.0  78.0     No
    ## 333  118.0  81.0     No
    ## 334  148.0  85.0     No
    ## 335  131.0  82.0     No
    ## 336  112.0  74.0     No
    ## 337  133.0  78.0     No
    ## 338  127.0  70.0     No
    ## 339   97.5  57.5     No
    ## 340  125.0  79.0     No
    ## 341  189.0 111.0     No
    ## 343  120.0  72.0     No
    ## 344  101.0  67.0     No
    ## 345  110.0  78.5     No
    ## 346  204.0  94.0    Yes
    ## 347  145.0  99.0     No
    ## 348  126.0  87.0     No
    ## 349  126.5  78.5     No
    ## 350  143.5  89.0     No
    ## 352  160.0  82.0     No
    ## 353  158.0 101.0     No
    ## 354  105.0  75.0     No
    ## 355  111.0  60.0     No
    ## 358  136.0  66.0     No
    ## 360  144.0  98.0    Yes
    ## 361  128.0  86.0     No
    ## 362  130.0  86.0     No
    ## 363  110.0  69.0     No
    ## 364  157.5 104.5     No
    ## 365  110.0  67.0     No
    ## 366  106.0  71.0     No
    ## 367  121.0  83.0     No
    ## 368  137.0  89.5     No
    ## 370  140.5  90.0     No
    ## 371  127.5  83.0     No
    ## 372  118.0  79.0     No
    ## 373  171.0 112.0     No
    ## 374  146.0  95.0     No
    ## 375  100.0  70.0     No
    ## 376  118.0  78.0     No
    ## 377   98.0  57.0     No
    ## 378  106.0  66.0     No
    ## 379  130.0  89.0     No
    ## 380  133.0  80.0     No
    ## 381  114.0  67.0     No
    ## 382  128.0  85.0     No
    ## 383  215.0 105.0     No
    ## 384  122.0  85.0     No
    ## 385  159.0 102.0     No
    ## 386  145.0  89.0     No
    ## 387   95.0  55.0     No
    ## 388  156.0  93.0    Yes
    ## 389  155.0  82.5     No
    ## 390  116.0  87.0     No
    ## 391  122.0  74.0     No
    ## 392  109.0  73.0     No
    ## 393  170.0  98.0     No
    ## 394  123.0  89.0     No
    ## 395  125.0  80.0     No
    ## 396  154.0 100.0    Yes
    ## 398  135.0  85.0     No
    ## 399  122.5  66.5     No
    ## 400  120.0  80.0     No
    ## 401  178.0  88.5     No
    ## 402  124.0  87.0     No
    ## 403  127.5  86.5     No
    ## 404  111.0  67.0     No
    ## 405  146.5  83.0     No
    ## 406  117.0  76.0     No
    ## 407  124.0  74.0     No
    ## 408  134.0  86.5     No
    ## 409  149.0  96.0     No
    ## 411  136.5  79.0     No
    ## 412  122.5  77.5     No
    ## 413  132.5  66.0     No
    ## 414  159.0  82.0     No
    ## 416  113.0  68.0     No
    ## 417  113.5  75.0     No
    ## 418  116.0  67.0     No
    ## 419  119.0  75.0     No
    ## 420  138.0  84.0     No
    ## 421  141.0  91.0     No
    ## 423  176.0 119.0     No
    ## 424  137.0  86.0     No
    ## 425  177.5 120.0     No
    ## 426  131.0  86.0     No
    ## 427  156.0 100.0     No
    ## 430  197.0 118.0     No
    ## 431  221.0 118.0    Yes
    ## 432  130.0  80.0     No
    ## 433  111.0  72.0     No
    ## 434  148.0 106.0     No
    ## 436  119.0  82.0     No
    ## 437   90.0  60.0     No
    ## 439  112.5  70.0     No
    ## 440  160.0 120.0     No
    ## 441  110.0  80.0     No
    ## 442  147.0  90.0     No
    ## 444  135.0  74.0     No
    ## 448  109.5  72.0     No
    ## 449   96.0  64.0     No
    ## 450  117.5  80.0     No
    ## 452  120.0  80.0     No
    ## 453  165.0 105.0     No
    ## 455  141.0  89.0     No
    ## 456  124.0  76.0     No
    ## 457  113.0  80.0     No
    ## 458  116.0  81.0     No
    ## 459  140.0  87.0     No
    ## 460  107.0  68.0     No
    ## 461   95.5  59.0     No
    ## 462  146.0  91.0     No
    ## 463  130.0  70.0     No
    ## 464  209.0 133.0     No
    ## 465  112.5  65.0     No
    ## 467  115.0  75.0     No
    ## 468  112.0  76.0     No
    ## 469  132.0  82.0     No
    ## 470  107.5  80.0     No
    ## 471  127.5  83.5     No
    ## 472  136.5  77.5     No
    ## 475  129.0  81.0     No
    ## 476  150.0  85.0    Yes
    ## 477  109.0  71.0     No
    ## 478  124.0  83.0     No
    ## 479  110.0  70.0     No
    ## 480  135.0  89.0     No
    ## 481  146.0 104.0     No
    ## 482  152.0  99.0     No
    ## 483  154.0  95.5     No
    ## 484  125.0  74.0     No
    ## 486  110.0  64.0     No
    ## 487  162.5  92.5     No
    ## 488  134.0  86.0     No
    ## 489  177.5  95.0    Yes
    ## 491  107.0  66.0     No
    ## 492  150.0  96.5     No
    ## 493  116.0  74.5     No
    ## 494  127.0  74.0     No
    ## 495  125.0  80.0     No
    ## 496  138.0  90.0     No
    ## 497  145.0 100.0     No
    ## 498  103.0  70.0     No
    ## 499  145.0  85.0     No
    ## 500  110.0  75.0     No
    ## 501  295.0 135.0     No
    ## 502  158.0  90.0     No
    ## 503  125.0  80.0     No
    ## 504  126.0  81.0     No
    ## 505  145.0  90.0     No
    ## 506  120.0  80.0     No
    ## 508  150.0 101.0     No
    ## 509   98.0  64.5     No
    ## 511  130.0  72.0     No
    ## 512  134.0  78.0     No
    ## 513  135.0  93.0     No
    ## 514  134.0  89.5     No
    ## 515  114.0  81.0     No
    ## 516  137.0  76.0     No
    ## 517  151.0  95.0     No
    ## 518  121.0  76.0     No
    ## 519  111.0  68.5     No
    ## 521  103.0  71.0     No
    ## 522  155.0  93.0     No
    ## 523  132.0  82.0     No
    ## 525  144.0  88.0     No
    ## 526  108.0  71.0     No
    ## 527  160.5  98.5     No
    ## 528  139.0  81.0     No
    ## 529  134.5  89.0     No
    ## 530  115.5  62.0     No
    ## 531  132.0  84.0     No
    ## 533  155.0  84.0     No
    ## 534  137.5  89.5     No
    ## 535  133.0  84.0     No
    ## 536  135.0  95.0     No
    ## 537  130.0  78.0     No
    ## 538  155.0 100.0     No
    ## 539  122.0  72.5     No
    ## 541  117.5  65.0     No
    ## 542  129.0  80.0     No
    ## 543  148.0  96.0     No
    ## 544  132.0  84.5     No
    ## 545  173.0 117.0     No
    ## 546  119.0  80.0     No
    ## 547  130.0  82.0     No
    ## 548  122.5  76.0     No
    ## 549  105.0  74.0     No
    ## 550  120.0  73.0     No
    ## 551  125.0  83.0     No
    ## 552  189.0 121.0     No
    ## 553  113.0  75.0     No
    ## 554  123.0  83.0     No
    ## 555  145.0 102.0     No
    ## 556  105.0  59.5     No
    ## 557  130.0  85.0     No
    ## 558  171.0 101.0     No
    ## 559  114.0  79.0     No
    ## 560  125.0  83.0     No
    ## 562  131.0  88.0    Yes
    ## 563  115.0  72.5     No
    ## 564  145.0  74.0     No
    ## 565  110.0  71.0     No
    ## 566  132.0  92.0     No
    ## 567  117.0  65.0     No
    ## 568  112.0  83.0     No
    ## 569  121.5  78.0     No
    ## 570  144.0  96.5     No
    ## 571  114.0  80.0     No
    ## 572  128.0  94.0     No
    ## 574  142.5  90.0     No
    ## 575  174.5 103.0     No
    ## 576  140.0  82.0     No
    ## 577  110.0  75.0     No
    ## 578  129.0  76.0     No
    ## 579  130.0  85.0     No
    ## 580  136.5  83.5     No
    ## 581  120.0  79.0     No
    ## 582  122.5  75.0     No
    ## 583  113.5  75.5     No
    ## 584  128.0  93.0     No
    ## 586  164.0 102.5     No
    ## 587  116.0  72.0     No
    ## 588  125.0  80.0     No
    ## 589  165.0  77.0     No
    ## 591  108.0  66.0     No
    ## 592  135.0  89.0     No
    ## 593  121.0  82.0     No
    ## 594  133.0  91.0     No
    ## 595  120.0  74.0     No
    ## 596  150.0  96.5     No
    ## 597  148.0 100.0     No
    ## 598  122.0  76.0     No
    ## 601  121.5  82.5     No
    ## 602  163.0  78.0     No
    ## 603  139.0  92.0     No
    ## 604  117.0  80.0     No
    ## 605  146.0  94.0     No
    ## 606  125.0  74.0     No
    ## 607  118.5  85.5     No
    ## 608  142.0  87.0     No
    ## 609   98.0  65.0     No
    ## 610  116.0  73.0     No
    ## 611  185.0 114.0     No
    ## 612  126.0  82.0     No
    ## 613  220.0 118.0    Yes
    ## 615  115.0  72.5     No
    ## 616  123.0  88.0     No
    ## 617  125.0  73.5     No
    ## 618  119.0  82.0     No
    ## 619  121.5  74.0     No
    ## 620  144.0  78.0     No
    ## 621  120.0  70.0     No
    ## 622  123.0  90.0     No
    ## 623  113.0  75.0     No
    ## 624  103.0  73.0     No
    ## 625  132.0  91.0     No
    ## 626  140.0  73.0     No
    ## 627  164.0  88.0     No
    ## 628  153.0 108.0     No
    ## 629  145.0  79.0     No
    ## 630  120.5  78.0     No
    ## 631  145.0  91.0     No
    ## 632  116.0  86.0     No
    ## 634  155.0  92.5     No
    ## 635   98.5  69.5     No
    ## 636  107.0  61.0     No
    ## 637  108.0  70.0     No
    ## 638  127.5  72.5     No
    ## 639  120.5  75.0     No
    ## 640  117.5  77.5     No
    ## 641  161.0  96.0     No
    ## 642  148.0  98.0     No
    ## 643  139.5  89.0     No
    ## 644  127.5  80.0     No
    ## 645  127.5  87.5     No
    ## 647  105.0  68.0     No
    ## 648  180.0 101.0     No
    ## 649  123.0  92.5     No
    ## 650  115.0  71.0     No
    ## 651   96.0  67.0     No
    ## 652  130.5  85.5     No
    ## 653  117.0  75.0     No
    ## 654  168.5 108.0     No
    ## 655  103.5  70.0     No
    ## 656  129.5  93.0     No
    ## 657  148.0 108.0     No
    ## 658  134.5  78.5     No
    ## 659  117.5  67.5     No
    ## 660  176.0  99.0     No
    ## 661  111.0  73.0     No
    ## 662  109.0  71.0     No
    ## 663  125.0  80.0     No
    ## 664  141.0  93.0     No
    ## 665  131.0  85.5     No
    ## 666  119.0  77.0     No
    ## 667  123.5  81.0     No
    ## 668  163.5 103.0     No
    ## 669  117.0  90.0     No
    ## 671  127.5  81.0     No
    ## 672  128.5  83.0     No
    ## 673  128.5  74.5     No
    ## 674  134.0  84.0     No
    ## 675  134.0  81.0     No
    ## 676  111.0  70.0     No
    ## 677  128.5  82.0     No
    ## 678  148.0  91.0     No
    ## 679  107.0  68.5     No
    ## 680  129.0  73.0     No
    ## 681  130.0  86.5     No
    ## 682  121.0  75.0     No
    ## 683  112.0  76.0     No
    ## 684  127.0  78.0     No
    ## 685  126.5  92.5    Yes
    ## 687  138.0  91.0     No
    ## 688  167.0  89.5     No
    ## 689  205.5 104.5     No
    ## 690  150.0  92.0    Yes
    ## 691  136.0  92.0     No
    ## 692  118.5  73.0     No
    ## 693  112.0  68.5     No
    ## 694  145.0  81.0     No
    ## 695  138.0  79.0     No
    ## 696  115.0  70.0     No
    ## 697  110.0  78.0     No
    ## 699  119.5  85.5     No
    ## 701  119.0  80.0     No
    ## 702  167.5  95.0     No
    ## 703  130.0  78.0     No
    ## 705  136.5  85.5     No
    ## 706  148.0  92.5     No
    ## 707  166.0 102.0     No
    ## 708  144.0  80.0     No
    ## 709  137.0  79.0     No
    ## 710  152.5  82.0     No
    ## 711  127.0  81.0     No
    ## 712  119.0  83.0     No
    ## 713  125.0  74.0     No
    ## 714  172.5  75.0     No
    ## 716  126.0  84.0     No
    ## 717  176.0  98.0     No
    ## 718  145.0  95.0     No
    ## 719  146.0  93.5     No
    ## 721  132.5  91.0     No
    ## 722  105.0  70.0     No
    ## 723  117.0  68.0     No
    ## 724  186.0 102.0     No
    ## 725  183.0  93.0    Yes
    ## 726  118.0  73.0     No
    ## 727  110.0  67.5     No
    ## 728  136.0  86.0     No
    ## 729  109.0  72.5     No
    ## 730  135.0  85.5     No
    ## 732  140.5  84.0     No
    ## 733  114.0  79.0     No
    ## 734  148.0  84.0     No
    ## 735  120.0  80.0     No
    ## 736  134.0  80.0     No
    ## 737  134.0  90.5     No
    ## 738  107.5  67.5     No
    ## 739  126.0  78.0     No
    ## 741  127.5  80.0     No
    ## 742  141.0  99.0     No
    ## 743  124.0  84.0     No
    ## 744  141.0  79.0     No
    ## 746  120.0  80.0     No
    ## 747  153.5 105.0     No
    ## 748  134.0  79.0     No
    ## 749  152.0  99.0    Yes
    ## 750  163.5  93.0     No
    ## 751  135.0  87.0     No
    ## 752  130.0  77.5     No
    ## 753  103.0  68.0     No
    ## 754  117.0  74.0     No
    ## 755  107.0  73.0     No
    ## 756  122.0  78.0     No
    ## 757  155.0  92.5     No
    ## 758  117.0  73.0     No
    ## 759  120.0  75.0     No
    ## 760  131.0  74.0     No
    ## 761  144.0  84.0     No
    ## 762  147.5  92.5     No
    ## 763  212.0 103.0     No
    ## 765  183.0 114.5    Yes
    ## 766  120.0  74.0     No
    ## 767  178.0  91.0     No
    ## 768  128.0  76.0     No
    ## 769  126.0  80.0     No
    ## 770  132.0  78.0     No
    ## 771  127.0  76.5     No
    ## 772  143.5  92.0     No
    ## 773  175.0  94.0     No
    ## 774  160.0 108.0     No
    ## 775  119.0  82.0     No
    ## 776  142.5  90.0     No
    ## 777  100.0  73.0     No
    ## 778  192.0 102.0    Yes
    ## 779  127.0  75.0     No
    ## 780  110.0  62.5     No
    ## 781  129.0  85.0     No
    ## 782  117.5  75.0     No
    ## 784  134.5  83.0     No
    ## 786  131.0  87.0     No
    ## 787  160.0  92.0     No
    ## 789   96.5  71.0     No
    ## 790  107.0  66.5     No
    ## 791  157.0  94.0     No
    ## 792  125.0  89.0     No
    ## 794  180.0 108.0     No
    ## 795  185.0 100.0     No
    ## 796  147.0  94.0     No
    ## 798  123.0  69.0     No
    ## 799  148.0  85.5     No
    ## 800  141.0 102.0     No
    ## 801  133.0  83.0     No
    ## 802  118.0  84.0     No
    ## 803  160.0 100.0     No
    ## 804  114.0  66.5     No
    ## 808  118.0  82.0     No
    ## 809  112.0  62.5     No
    ## 811  150.0 101.0     No
    ## 812  137.0  81.0     No
    ## 813  118.5  82.5     No
    ## 814  126.0  82.0    Yes
    ## 815  121.5  74.0     No
    ## 816  133.0  76.0     No
    ## 817  140.5  94.0     No
    ## 818  159.5  91.0     No
    ## 819  118.5  69.0     No
    ## 820  164.0 107.0     No
    ## 821  118.5  82.0     No
    ## 822  129.0  90.0     No
    ## 823  130.0  92.0     No
    ## 824  177.0 110.0    Yes
    ## 825  156.0 104.0     No
    ## 827  116.0  78.5     No
    ## 828  152.0  82.0    Yes
    ## 829  117.0  80.0     No
    ## 830  127.0  79.0     No
    ## 831  150.0  90.0     No
    ## 832  162.0  87.5     No
    ## 833  124.0  86.5     No
    ## 834  132.0  91.0     No
    ## 835  130.0  82.0     No
    ## 836  120.0  70.0     No
    ## 837  142.0  93.0     No
    ## 838  125.0  72.5     No
    ## 839  164.0  94.5     No
    ## 840  112.5  82.5     No
    ## 841  108.0  66.0     No
    ## 844  142.0  82.0     No
    ## 845  135.0  79.0     No
    ## 846  111.0  70.0     No
    ## 847  128.0  82.0     No
    ## 848  124.0  79.0     No
    ## 850  129.0  86.0     No
    ## 851  168.0  82.0     No
    ## 853  102.5  65.0     No
    ## 854  110.0  78.5     No
    ## 855  134.0  87.0     No
    ## 856  127.5  85.0     No
    ## 857  122.0  70.0     No
    ## 858  126.0  85.0     No
    ## 860  105.0  84.0     No
    ## 861  147.0  90.0     No
    ## 863  200.0 140.0     No
    ## 865  127.0  75.0     No
    ## 866  128.0  91.0     No
    ## 867  137.0  96.5     No
    ## 868  138.0  83.0     No
    ## 869  133.0  80.0     No
    ## 870  120.0  73.0     No
    ## 871  122.0  82.0     No
    ## 872  121.5  78.0     No
    ## 873  116.0  79.0     No
    ## 874  132.0  84.0     No
    ## 875  131.5  78.5     No
    ## 876  118.0  66.0     No
    ## 877  118.5  80.5     No
    ## 878  114.0  70.5     No
    ## 879  102.0  69.0     No
    ## 880  131.0  82.0     No
    ## 881  119.0  85.0     No
    ## 882  107.0  76.5     No
    ## 883  122.5  67.5     No
    ## 884  131.0  80.0     No
    ## 885  114.0  83.0     No
    ## 886  111.0  62.0     No
    ## 887  118.0  86.0     No
    ## 888  146.0  88.0     No
    ## 889  141.0  81.0     No
    ## 890  107.0  73.0     No
    ## 891  128.0  83.0     No
    ## 892  163.0  85.0     No
    ## 893  136.0  99.0     No
    ## 894  244.0 124.0    Yes
    ## 895  134.0  92.0     No
    ## 896  113.0  77.0     No
    ## 897  113.5  74.0     No
    ## 898  104.0  76.0     No
    ## 899  135.0  82.0     No
    ## 900  114.0  81.0     No
    ## 901  131.0  79.0     No
    ## 904  127.5  75.0     No
    ## 905  122.0  73.0     No
    ## 906  128.0  86.5     No
    ## 907  127.5  83.5     No
    ## 908  153.0  75.0     No
    ## 909  168.0 102.0     No
    ## 910  127.5  81.0     No
    ## 911  134.0  91.0     No
    ## 912  113.5  73.0     No
    ## 913  156.0  98.0     No
    ## 915  112.0  69.0     No
    ## 916  164.0  85.5     No
    ## 917  176.0  97.0     No
    ## 918  121.0  80.0     No
    ## 919  141.0  83.0     No
    ## 920  148.0  93.0     No
    ## 921  123.0  73.5     No
    ## 922  152.5  90.0    Yes
    ## 923  108.0  72.0     No
    ## 924  111.0  61.0     No
    ## 925  126.0  85.0     No
    ## 926  141.0  92.0     No
    ## 927  130.0  82.5     No
    ## 929  108.0  77.0     No
    ## 930  168.0 102.0     No
    ## 931  146.0  89.0     No
    ## 932  128.0  77.0     No
    ## 933  134.5  87.0     No
    ## 934  124.5  82.0     No
    ## 935  213.0  94.5     No
    ## 936  115.0  81.0     No
    ## 938  127.0  76.0     No
    ## 939  129.0  92.0     No
    ## 940  126.0  82.0     No
    ## 941  167.5  91.5     No
    ## 942  108.0  66.0     No
    ## 943  125.0  75.0     No
    ## 944  153.0  94.0     No
    ## 945  117.0  71.0     No
    ## 946  152.0 103.0     No
    ## 947  118.0  74.0     No
    ## 948  154.0  88.5     No
    ## 950  114.0  84.0     No
    ## 952  120.0  90.0     No
    ## 953  111.0  67.5     No
    ## 954  116.0  82.0     No
    ## 955  115.0  80.0     No
    ## 957  159.0 115.0     No
    ## 958  134.0  84.0     No
    ## 959  206.0  98.0     No
    ## 960  123.0  71.0     No
    ## 961  157.0  99.0     No
    ## 962  142.0  74.0     No
    ## 963  106.0  70.5     No
    ## 964  139.0  88.0     No
    ## 965  242.0 141.0    Yes
    ## 966  111.0  70.0     No
    ## 967  151.0  96.0     No
    ## 968  199.0  83.0     No
    ## 969  184.0 109.0     No
    ## 970  129.5  87.0     No
    ## 971  127.5  80.0     No
    ## 973  171.5  83.0     No
    ## 974  129.0  90.0     No
    ## 975  133.0  86.0     No
    ## 977  167.5 102.5     No
    ## 978  160.0  85.0     No
    ## 981  127.5  75.0     No
    ## 982  115.5  68.0     No
    ## 983  128.0  85.0     No
    ## 984  118.0  77.0     No
    ## 985  198.0 107.0     No
    ## 986  117.5  67.5     No
    ## 987  127.0  86.0     No
    ## 989  114.5  62.5     No
    ## 990  137.5  77.5     No
    ## 991  151.0  98.0     No
    ## 992  172.5  88.0     No
    ## 993  151.0  89.5    Yes
    ## 994  158.5 107.0     No
    ## 995  118.0  84.0     No
    ## 996  160.0  97.0     No
    ## 997  137.0  82.5     No
    ## 998  131.5  92.0     No
    ## 999  123.5  84.0     No
    ## 1000 140.0  82.0     No
    ## 1001  98.0  66.0     No
    ## 1002 109.5  67.0     No
    ## 1003 115.5  65.5     No
    ## 1004 106.0  79.0     No
    ## 1006 122.0  81.0     No
    ## 1007 137.0  88.0     No
    ## 1008 160.0 109.0     No
    ## 1009 134.0  98.0     No
    ## 1010 125.5  86.0     No
    ## 1012 120.5  67.5     No
    ## 1013 131.0  74.0     No
    ## 1014 159.0 109.0     No
    ## 1015 170.0 100.0    Yes
    ## 1016 124.0  90.0     No
    ## 1017 154.0 100.0     No
    ## 1018 132.0  75.0     No
    ## 1020 111.0  72.0     No
    ## 1021 110.0  77.0     No
    ## 1022 149.0  81.0     No
    ## 1023 135.0  82.5     No
    ## 1024 127.0  89.0     No
    ## 1025 121.0  67.0     No
    ## 1026 115.0  82.0     No
    ## 1027 139.0  88.0     No
    ## 1028 122.0  80.0     No
    ## 1029 113.0  80.0     No
    ## 1030 124.0  86.0     No
    ## 1031 139.0  92.0     No
    ## 1032 118.0  77.5     No
    ## 1033 122.0  76.0     No
    ## 1035 107.5  71.0     No
    ## 1036 125.0  85.0     No
    ## 1039 102.0  64.0     No
    ## 1040 124.0  85.0     No
    ## 1042 206.0 104.0     No
    ## 1043 117.0  67.0     No
    ## 1044 145.0  94.0     No
    ## 1045 107.0  73.5     No
    ## 1046 142.0  94.0     No
    ## 1047 167.0  96.5     No
    ## 1048 142.0  82.5     No
    ## 1049 111.5  70.0     No
    ## 1050 148.0  81.0     No
    ## 1051 133.0  88.0     No
    ## 1052 100.0  71.0     No
    ## 1053 132.0  86.0     No
    ## 1054 139.0  84.0     No
    ## 1055 119.5  74.0     No
    ## 1056 127.0  83.0     No
    ## 1057 147.0  71.5     No
    ## 1059 136.5  95.0     No
    ## 1060 149.0  90.0     No
    ## 1061 155.0  74.0     No
    ## 1062 133.0  84.0     No
    ## 1063 131.5  84.0     No
    ## 1064 117.5  77.5     No
    ## 1065 131.0  88.0     No
    ## 1068 105.5  67.0     No
    ## 1069 138.0  86.0     No
    ## 1070 126.0  84.0     No
    ## 1071 145.0  90.0     No
    ## 1072 122.0  83.0     No
    ## 1073 141.5  91.0     No
    ## 1074 111.0  71.0     No
    ## 1075 142.0  90.0     No
    ## 1076 131.0  88.0     No
    ## 1078 117.0  80.0     No
    ## 1079 143.0  98.0     No
    ## 1080 165.0 115.0    Yes
    ## 1081 161.5  95.0     No
    ## 1082 122.0  66.5     No
    ## 1084 164.5  93.5     No
    ## 1085 120.0  80.0     No
    ## 1086 141.5  98.0     No
    ## 1087 124.0  82.0     No
    ## 1088 112.0  77.0     No
    ## 1089 137.0  97.0     No
    ## 1090 163.0 105.0     No
    ## 1091 123.0  76.5     No
    ## 1092 116.0  62.0     No
    ## 1093 171.5 105.5     No
    ## 1094 150.0  97.0     No
    ## 1095 126.0  68.0     No
    ## 1096 132.5  90.0     No
    ## 1097 146.0  88.0     No
    ## 1098 177.0 103.5     No
    ## 1099 105.0  74.0     No
    ## 1100 108.5  63.5     No
    ## 1101 110.0  60.0     No
    ## 1102 102.0  72.0     No
    ## 1105 129.5  92.5     No
    ## 1106 152.5  91.0     No
    ## 1107 132.0  85.0     No
    ## 1108 131.0  80.0     No
    ## 1109 146.5  82.0     No
    ## 1110 114.0  77.0     No
    ## 1112 124.0  70.0     No
    ## 1113 160.0  93.0    Yes
    ## 1115 153.0  93.0     No
    ## 1116 109.5  69.0     No
    ## 1118 132.0  67.0     No
    ## 1119 131.0  79.0     No
    ## 1120 128.0  80.5     No
    ## 1121 125.0  83.0    Yes
    ## 1123 201.0  93.0    Yes
    ## 1124 120.5  84.0     No
    ## 1125 145.0  90.0     No
    ## 1126 144.0  86.0     No
    ## 1127 127.0  70.0     No
    ## 1128 122.5  82.5     No
    ## 1129 131.0  81.0     No
    ## 1130 118.0  79.0     No
    ## 1131 112.5  76.5     No
    ## 1132 135.0  91.0     No
    ## 1134 118.0  78.0     No
    ## 1135 123.0  75.0     No
    ## 1136 112.0  66.0     No
    ## 1137 137.0  87.0     No
    ## 1138 114.0  81.0     No
    ## 1139 137.0  76.0     No
    ## 1140 125.0  71.0     No
    ## 1141 113.5  66.5     No
    ## 1142 113.0  76.5     No
    ## 1143 114.0  78.0     No
    ## 1144 110.0  72.5     No
    ## 1145 125.0  87.0     No
    ## 1147 131.0  87.0     No
    ## 1148 149.0  73.0     No
    ## 1151 137.0  99.0     No
    ## 1152 213.0 103.0    Yes
    ## 1153 105.0  65.0     No
    ## 1154 130.0  77.5     No
    ## 1155 145.0  90.0     No
    ## 1156 148.5  90.0    Yes
    ## 1157 159.5  94.0     No
    ## 1158 143.5  89.0     No
    ## 1159 140.0  80.0     No
    ## 1160 119.0  75.0     No
    ## 1162 105.5  64.5     No
    ## 1163 127.5  79.5     No
    ## 1166 114.0  71.0     No
    ## 1167 101.0  72.0     No
    ## 1169 133.0  72.5     No
    ## 1170 138.0  92.0     No
    ## 1171 173.0 102.0     No
    ## 1172 112.5  64.5     No
    ## 1173 109.0  74.0     No
    ## 1174 189.0 103.0     No
    ## 1176 141.5  88.5     No
    ## 1177 159.0  64.0     No
    ## 1178 146.0  77.0     No
    ## 1179 140.0  84.0     No
    ## 1180 122.0  87.5     No
    ## 1181 130.0  80.0     No
    ## 1182 125.0  80.5     No
    ## 1183 132.0  79.0     No
    ## 1184 159.0  89.0     No
    ## 1186 144.0  99.0     No
    ## 1187 133.0  88.0     No
    ## 1188 147.0  96.0     No
    ## 1190 141.0  73.5     No
    ## 1191 156.0  90.0     No
    ## 1192 174.0 112.0     No
    ## 1193 150.0  84.0     No
    ## 1194 140.0  88.0     No
    ## 1195 129.0  82.0     No
    ## 1196 130.0  84.0     No
    ## 1197 100.0  60.0     No
    ## 1198 113.0  68.0     No
    ## 1199 128.0  84.5     No
    ## 1200 145.0  77.0     No
    ## 1201 172.0  84.0     No
    ## 1202 167.0 107.0     No
    ## 1204 116.0  69.0     No
    ## 1205 122.0  73.0     No
    ## 1206 107.5  65.0     No
    ## 1207 144.0  91.0     No
    ## 1208 188.0 107.0    Yes
    ## 1209 119.0  73.0     No
    ## 1211 100.0  65.0     No
    ## 1212 130.0  80.0     No
    ## 1213 116.0  79.0     No
    ## 1214 147.0  86.0     No
    ## 1215 140.0  90.0     No
    ## 1216 114.0  81.0     No
    ## 1217 153.5 103.5     No
    ## 1218 127.5  67.5     No
    ## 1219 126.0  86.0     No
    ## 1220 138.0  76.0     No
    ## 1221 175.0 107.5     No
    ## 1223 123.0  69.0     No
    ## 1224 165.0  91.0     No
    ## 1226 131.0  83.0     No
    ## 1227 135.0  94.5     No
    ## 1229 131.0  84.0     No
    ## 1231 145.0  87.5     No
    ## 1233 130.0  94.0     No
    ## 1234 184.0 106.0     No
    ## 1235 135.0  80.0     No
    ## 1236 118.0  76.5     No
    ## 1237 160.5 109.0     No
    ## 1238 114.0  64.0     No
    ## 1239 177.0  96.0     No
    ## 1240 243.0 142.5     No
    ## 1242 178.0 103.0     No
    ## 1243 139.0  80.0     No
    ## 1244 127.0  77.0     No
    ## 1245 184.0  90.0     No
    ## 1246 145.5  87.5     No
    ## 1247 133.0  87.0     No
    ## 1248 187.5  85.5    Yes
    ## 1249 119.0  76.0     No
    ## 1250 125.0  86.0     No
    ## 1251 165.0  96.0     No
    ## 1252 138.0  92.0     No
    ## 1253 126.0  84.0     No
    ## 1255 178.0 106.0     No
    ## 1256 138.0 105.0     No
    ## 1257 123.0  85.5     No
    ## 1259 110.0  69.0     No
    ## 1260 122.0  82.0     No
    ## 1261 119.0  76.0     No
    ## 1262 118.0  72.0     No
    ## 1263 110.0  78.0     No
    ## 1265 155.0  85.0     No
    ## 1266 111.0  81.0     No
    ## 1268 150.0  88.0     No
    ## 1269 116.0  75.0     No
    ## 1270 117.0  76.0     No
    ## 1271 173.0  84.0     No
    ## 1272  99.0  67.0     No
    ## 1273 115.0  77.0     No
    ## 1274 181.0  90.0     No
    ## 1275 115.0  75.0     No
    ## 1276 174.0 110.0     No
    ## 1277  96.0  61.0     No
    ## 1278 122.0  72.0     No
    ## 1279 120.5  85.5     No
    ## 1280 127.0  82.0     No
    ## 1281 102.5  67.5     No
    ## 1283 130.0  80.0     No
    ## 1284 124.0  81.0     No
    ## 1285 180.0 109.5     No
    ## 1286 112.0  74.5     No
    ## 1287 133.5  76.0     No
    ## 1288 100.5  66.0     No
    ## 1290 120.0  86.0     No
    ## 1291 124.0  70.0     No
    ## 1292 156.0  92.0     No
    ## 1293 156.5 105.0     No
    ## 1295 140.0  88.0     No
    ## 1296 128.0  69.0     No
    ## 1297 119.0  84.0     No
    ## 1298 130.0  68.0     No
    ## 1299 118.0  78.0     No
    ## 1300 156.0  90.0     No
    ## 1301 148.0 100.0     No
    ## 1302 113.0  77.0     No
    ## 1303 159.0  88.0     No
    ## 1304 127.0  93.0     No
    ## 1306 143.0 101.0     No
    ## 1307 125.0  92.0     No
    ## 1308 144.0  86.0     No
    ## 1310 122.0  70.0     No
    ## 1311 129.0  89.0     No
    ## 1312 108.5  71.5     No
    ## 1313 162.0  99.0     No
    ## 1314 135.5  80.0     No
    ## 1316 144.0  85.0     No
    ## 1317 123.0  71.0     No
    ## 1318 128.0  76.0    Yes
    ## 1319 126.5  67.5     No
    ## 1320 131.5  89.0     No
    ## 1321 120.0  74.0     No
    ## 1322 143.0 104.0     No
    ## 1323 128.0  74.0     No
    ## 1324 160.0  96.0     No
    ## 1325 122.0  74.0     No
    ## 1326  97.0  63.0     No
    ## 1327 146.0  95.0     No
    ## 1328 155.0 102.0     No
    ## 1329 128.0  90.0     No
    ## 1330 114.0  85.0     No
    ## 1332 175.0  95.0     No
    ## 1333 125.0  94.0     No
    ## 1334 121.0  85.5     No
    ## 1335 116.0  77.0     No
    ## 1336 117.5  80.0     No
    ## 1337 108.0  64.0     No
    ## 1338 132.5  87.5     No
    ## 1339 108.0  70.0     No
    ## 1341 108.5  75.0     No
    ## 1343 110.0  70.0     No
    ## 1345 140.0  85.0     No
    ## 1346 110.0  73.0     No
    ## 1347 103.0  61.0     No
    ## 1348 122.5  80.0     No
    ## 1349 115.0  81.0     No
    ## 1350 163.0 105.0     No
    ## 1351  95.0  58.0     No
    ## 1352 134.0  90.0     No
    ## 1354 173.0  59.0     No
    ## 1355 148.0  99.0     No
    ## 1356 131.0  87.0     No
    ## 1358 172.5  98.0     No
    ## 1359 174.0 100.0     No
    ## 1360 125.0  89.0     No
    ## 1361 126.0  82.0     No
    ## 1362 120.0  72.5     No
    ## 1363 137.0  75.0     No
    ## 1364 164.0  81.0     No
    ## 1365 121.0  74.0     No
    ## 1366 132.0  87.0     No
    ## 1368 115.0  70.0     No
    ## 1369 128.0  87.0     No
    ## 1370 108.0  67.0     No
    ## 1372 114.0  70.0     No
    ## 1373 136.0  70.0     No
    ## 1375 127.5  85.0     No
    ## 1376 106.0  78.0     No
    ## 1377 127.0  81.0     No
    ## 1378 126.0  70.0     No
    ## 1379 158.0  97.0     No
    ## 1380 154.0  98.0     No
    ## 1381 130.0  86.0     No
    ## 1382 124.0  78.0     No
    ## 1383 112.5  73.5     No
    ## 1384 135.0  90.0     No
    ## 1386 140.0 108.0     No
    ## 1387 133.5  87.5     No
    ## 1388 122.0  76.0     No
    ## 1389 120.0  84.0     No
    ## 1390 113.0  80.0     No
    ## 1391 114.0  72.0     No
    ## 1392 128.0  94.0     No
    ## 1394 151.0  94.0     No
    ## 1395 124.0  83.0     No
    ## 1396 115.0  82.0     No
    ## 1397 120.0  77.0     No
    ## 1398 102.0  66.5     No
    ## 1399 121.0  73.0     No
    ## 1400 127.5  80.0     No
    ## 1401 111.0  79.0     No
    ## 1402 112.5  70.0     No
    ## 1403 103.5  66.5     No
    ## 1404 142.5  97.5     No
    ## 1406 114.0  80.0     No
    ## 1407 110.0  73.0     No
    ## 1408 114.0  73.0     No
    ## 1409 123.0  80.0     No
    ## 1410 114.0  84.0     No
    ## 1411 124.0  83.0     No
    ## 1412 141.5  95.0     No
    ## 1414 119.0  73.0     No
    ## 1415 115.0  83.0     No
    ## 1416 163.0  86.0     No
    ## 1417 120.0  68.0     No
    ## 1419 145.0  89.0     No
    ## 1421 137.0  81.0     No
    ## 1422 127.0  76.0     No
    ## 1423 144.0  88.0     No
    ## 1424 138.0  87.0     No
    ## 1425 149.5  78.0     No
    ## 1426 146.0  89.0     No
    ## 1428 101.0  70.0     No
    ## 1429 133.0  83.0     No
    ## 1430 140.0  83.0     No
    ## 1431 128.0  77.0     No
    ## 1432 106.0  67.5     No
    ## 1433 113.5  74.0     No
    ## 1435 138.0  88.0     No
    ## 1437 135.0  76.0     No
    ## 1438 138.0  86.0     No
    ## 1440 182.5  88.0     No
    ## 1441 109.0  73.0     No
    ## 1442 132.0  86.0     No
    ## 1443 172.5 116.5     No
    ## 1444 121.0  79.0     No
    ## 1445 126.0  73.0     No
    ## 1446 112.0  82.0     No
    ## 1447 121.0  81.0     No
    ## 1448 134.0  87.0    Yes
    ## 1449 154.0  90.0     No
    ## 1450 125.0  85.0     No
    ## 1451 109.0  72.0     No
    ## 1452 117.0  79.0     No
    ## 1453 105.0  70.0     No
    ## 1456 115.0  76.0     No
    ## 1457 149.0 100.0     No
    ## 1458 123.0  78.0     No
    ## 1459 149.0  90.0     No
    ## 1461 122.5  82.5     No
    ## 1462 162.0  98.0     No
    ## 1463 129.0  86.0     No
    ## 1464 144.0  90.0     No
    ## 1465 137.0  85.0     No
    ## 1466 122.0  94.0     No
    ## 1467 113.0  79.0     No
    ## 1469 175.0  67.5     No
    ## 1470 143.0  87.0     No
    ## 1471 133.0  88.0     No
    ## 1473 100.0  64.0     No
    ## 1474 115.0  77.0     No
    ## 1475 142.0  85.5     No
    ## 1476 142.0  88.0     No
    ## 1477 132.0  70.0     No
    ## 1478 139.0  86.0     No
    ## 1479 131.0  85.0     No
    ## 1480 146.0  86.0     No
    ## 1481 123.0  87.0     No
    ## 1482 118.0  78.0     No
    ## 1483  95.0  57.0     No
    ## 1484 158.5 100.5     No
    ## 1485 142.0  92.0     No
    ## 1486  95.0  65.0     No
    ## 1487 115.0  74.0     No
    ## 1488 115.0  76.0     No
    ## 1489 131.0  89.0     No
    ## 1490 136.0  90.0     No
    ## 1492 123.0  72.0     No
    ## 1493 130.0  72.0     No
    ## 1494 127.5  77.5     No
    ## 1495 130.0  82.0     No
    ## 1496 121.0  82.0     No
    ## 1497 135.0  80.0     No
    ## 1498 145.0  82.0     No
    ## 1499 125.0  80.0     No
    ## 1500 136.5  83.0     No
    ## 1501 115.0  77.5     No
    ## 1502 117.0  81.0     No
    ## 1503 116.0  66.0     No
    ## 1504 165.0  80.0     No
    ## 1505 141.0  98.0     No
    ## 1506 120.0  74.0     No
    ## 1507 110.0  70.0     No
    ## 1508 162.5 105.0     No
    ## 1509 114.0  64.0     No
    ## 1511 127.0  80.0     No
    ## 1513 115.0  69.0     No
    ## 1514 110.0  72.0     No
    ## 1515 112.0  65.0     No
    ## 1516 118.0  82.0     No
    ## 1517 151.0  97.5     No
    ## 1518 109.0  81.0     No
    ## 1519 119.0  76.0     No
    ## 1520 135.0  70.0     No
    ## 1521 135.0  89.0     No
    ## 1522 100.0  60.0     No
    ## 1523 117.5  70.5     No
    ## 1524 139.0  93.0     No
    ## 1525 132.0  75.5     No
    ## 1526 125.5  94.0     No
    ## 1527 121.0  81.0     No
    ## 1528 160.0  95.0     No
    ## 1529 111.0  80.0     No
    ## 1530 127.0  79.0     No
    ## 1531 134.0  95.0     No
    ## 1532 148.0  93.0     No
    ## 1534 107.5  70.0     No
    ## 1535 125.0  82.0     No
    ## 1536 124.0  90.0     No
    ## 1537 156.0  96.0     No
    ## 1538 128.0  71.0     No
    ## 1539 132.0  94.0     No
    ## 1540 137.5  92.0     No
    ## 1542 137.0  92.5     No
    ## 1543 115.0  70.0     No
    ## 1544 131.5  85.0     No
    ## 1545 117.5  72.0     No
    ## 1546 129.0  78.0     No
    ## 1547 112.0  66.5     No
    ## 1548 132.0  80.0     No
    ## 1549 141.0  92.0     No
    ## 1550 130.5  98.0     No
    ## 1551 122.5  82.5     No
    ## 1552 109.0  70.0     No
    ## 1553 141.0  75.0     No
    ## 1554 126.0  66.0     No
    ## 1555 135.0  88.0     No
    ## 1557 123.0  69.0     No
    ## 1558 105.0  70.0     No
    ## 1559 121.0  88.0     No
    ## 1560 123.0  68.0     No
    ## 1561 132.5  85.5     No
    ## 1562 122.0  70.0     No
    ## 1563 127.5  87.0     No
    ## 1564 100.0  68.0     No
    ## 1565 126.0  85.0     No
    ## 1566 129.0  61.0     No
    ## 1567 111.0  63.0     No
    ## 1568 142.5  79.0     No
    ## 1569 122.5  77.0     No
    ## 1571 143.0  94.5     No
    ## 1572 112.5  80.0     No
    ## 1573 107.0  73.0     No
    ## 1574 142.5  93.5     No
    ## 1575 125.0  85.0     No
    ## 1576 146.0  78.0     No
    ## 1577 143.5  90.0     No
    ## 1578 129.0  82.0     No
    ## 1579 132.5  85.0     No
    ## 1580 121.0  66.0     No
    ## 1581 112.5  60.0     No
    ## 1582 129.0  81.0     No
    ## 1583 142.0 108.0     No
    ## 1584 199.0 106.0     No
    ## 1586 134.0  92.0     No
    ## 1587 112.0  80.0     No
    ## 1589 179.0 107.0     No
    ## 1590 125.0  85.0     No
    ## 1591 179.0 100.0     No
    ## 1592 101.0  67.0     No
    ## 1593 104.0  57.0     No
    ## 1594 112.0  71.0     No
    ## 1595 104.0  66.0     No
    ## 1596 138.0  92.0     No
    ## 1597 130.0  85.0     No
    ## 1598 120.0  72.0     No
    ## 1599 124.0  72.0     No
    ## 1600 133.0  69.0     No
    ## 1602 108.0  81.0     No
    ## 1603 129.0  80.0     No
    ## 1605 118.0  82.0     No
    ## 1606 186.5  97.0     No
    ## 1607 131.0  99.0     No
    ## 1608 134.0  75.0     No
    ## 1609 116.0  83.0     No
    ## 1610 138.5  87.5     No
    ## 1611 120.0  72.0     No
    ## 1612 138.0  91.0     No
    ## 1613 127.5  85.5     No
    ## 1614 145.0 100.5     No
    ## 1615 150.0  90.0    Yes
    ## 1616 121.0  71.0     No
    ## 1617 120.0  80.0     No
    ## 1618 117.0  78.0     No
    ## 1619 155.0  90.0     No
    ## 1620 164.0 102.0     No
    ## 1621 122.5  76.5     No
    ## 1622  96.0  62.0     No
    ## 1623 170.0 103.0    Yes
    ## 1624 173.0  89.0     No
    ## 1625 134.0  80.0     No
    ## 1626 186.0 101.0     No
    ## 1627 154.0 106.0     No
    ## 1628 102.0  65.0     No
    ## 1629 160.0  94.0     No
    ## 1630 103.0  73.0     No
    ## 1631 115.0  84.0     No
    ## 1633 204.0 118.0     No
    ## 1634 114.0  78.0     No
    ## 1636 127.5  90.0     No
    ## 1637 116.0  82.0     No
    ## 1638 131.5  85.0     No
    ## 1642 158.0  94.0     No
    ## 1643 130.0  80.0     No
    ## 1644 123.5  84.5     No
    ## 1645 120.0  80.0     No
    ## 1646 111.0  95.5     No
    ## 1648 112.5  75.5     No
    ## 1649 118.5  84.5     No
    ## 1650 112.5  70.0     No
    ## 1651 140.0  93.0     No
    ## 1652 139.0  90.0     No
    ## 1653 124.0  79.0     No
    ## 1654 110.0  65.0     No
    ## 1655 217.0 112.0     No
    ## 1656 196.0 116.0     No
    ## 1657 130.0  75.5     No
    ## 1658 146.0  89.0     No
    ## 1659 123.0  83.0     No
    ## 1660 175.0 110.0     No
    ## 1661 114.0  77.0     No
    ## 1663 176.0  89.0     No
    ## 1664 140.0  94.0     No
    ## 1665 138.0  82.0     No
    ## 1666 111.0  78.0     No
    ## 1667 156.5  85.0     No
    ## 1668 113.0  66.5     No
    ## 1669 142.0  76.0     No
    ## 1670 167.0  94.0    Yes
    ## 1671 113.0  70.0     No
    ## 1672 193.0 104.0     No
    ## 1673 187.0  95.5     No
    ## 1675 107.5  75.0     No
    ## 1676 147.0  97.0     No
    ## 1677 112.5  72.5     No
    ## 1678 196.0 120.0     No
    ## 1679 189.0  87.0     No
    ## 1680 168.0  98.0     No
    ## 1681 130.0  85.0     No
    ## 1682 126.5  85.0     No
    ## 1683 136.0  76.0     No
    ## 1684 196.0 119.0     No
    ## 1685 112.0  79.0     No
    ## 1688 116.0  82.0     No
    ## 1689 170.0 118.0     No
    ## 1690 190.0  97.0     No
    ## 1691 120.0  70.5     No
    ## 1693 125.0  88.0     No
    ## 1695 158.0 102.0     No
    ## 1696 104.0  73.0     No
    ## 1697 107.5  70.0     No
    ## 1698 114.0  68.0     No
    ## 1699 125.5  82.5     No
    ## 1700 117.5  80.0     No
    ## 1701 131.0  85.0     No
    ## 1702 161.0 100.5     No
    ## 1703 114.0  79.0     No
    ## 1704 110.5  66.0     No
    ## 1705 132.5  82.5     No
    ## 1706 122.0  81.0     No
    ## 1707 130.0  85.5     No
    ## 1708 121.5  72.5     No
    ## 1709 129.0  85.0     No
    ## 1710 151.5  96.5     No
    ## 1711 122.0  82.0     No
    ## 1712 185.0  86.0     No
    ## 1714 110.0  74.0     No
    ## 1715 113.5  72.5     No
    ## 1716 119.0  76.0     No
    ## 1718 154.0  98.0     No
    ## 1720 123.0  82.0     No
    ## 1721 155.5  99.5     No
    ## 1722 122.0  84.0     No
    ## 1723 111.5  74.0     No
    ## 1724 119.0  84.0     No
    ## 1725 101.0  62.0     No
    ## 1727 131.0  80.0     No
    ## 1728 130.0  86.0     No
    ## 1729 127.0  86.0     No
    ## 1731 143.5 109.0     No
    ## 1732 159.0  91.5     No
    ## 1733 121.5  69.0     No
    ## 1734 139.0  75.0     No
    ## 1735 155.0  81.0     No
    ## 1736 115.0  61.0     No
    ## 1737 124.0  78.0     No
    ## 1738 109.0  71.0     No
    ## 1739 145.0 100.0     No
    ## 1740 104.0  69.0     No
    ## 1741 128.0  92.5     No
    ## 1742 114.5  80.0     No
    ## 1744 130.0  94.0     No
    ## 1745 122.0  82.0     No
    ## 1746 115.0  80.0     No
    ## 1748 126.5  90.0     No
    ## 1749 144.0  90.0     No
    ## 1750 120.0  76.0     No
    ## 1751 117.5  73.5     No
    ## 1752 143.5  90.0     No
    ## 1754 103.0  65.0     No
    ## 1755 115.0  71.0     No
    ## 1756 112.0  63.0     No
    ## 1757 146.0  91.0     No
    ## 1758  92.0  69.0     No
    ## 1759 133.0  86.0     No
    ## 1761 193.0 132.0     No
    ## 1762 121.0  86.0     No
    ## 1763 119.0  72.0     No
    ## 1764 141.5  93.0     No
    ## 1765 120.0  87.0     No
    ## 1766 169.0  85.0     No
    ## 1768 134.0  90.0     No
    ## 1770 110.0  76.0     No
    ## 1771 131.0  87.0     No
    ## 1772 117.5  71.0     No
    ## 1773 140.0  93.0     No
    ## 1774 147.5  95.0     No
    ## 1775 110.0  69.0     No
    ## 1776 120.0  70.0     No
    ## 1777 120.0  78.0     No
    ## 1778 133.5  81.5     No
    ## 1779 146.5  77.5     No
    ## 1780 110.0  71.5     No
    ## 1782 120.0  81.0     No
    ## 1783 107.0  74.0     No
    ## 1784 130.0  85.0     No
    ## 1785 146.5  80.0     No
    ## 1786 102.0  64.0     No
    ## 1787 126.5  84.0     No
    ## 1788 146.5  97.0     No
    ## 1789 118.0  76.0     No
    ## 1790 138.5  99.0    Yes
    ## 1792 166.0  93.0    Yes
    ## 1793 145.5  82.0     No
    ## 1794 116.5  82.0     No
    ## 1797 142.0  54.0     No
    ## 1799 122.0  70.0     No
    ## 1800 110.0  71.0     No
    ## 1801 119.0  85.5     No
    ## 1802 153.0  75.0     No
    ## 1803 166.5 107.0     No
    ## 1804 135.0  80.0     No
    ## 1806 124.0  75.0     No
    ## 1807 151.0 102.0     No
    ## 1808 124.0  69.0     No
    ## 1809 119.0  78.0     No
    ## 1810 106.0  72.0     No
    ## 1811 148.0  89.0     No
    ## 1812 124.0  87.0     No
    ## 1813 124.0  78.0     No
    ## 1814 141.0  83.5     No
    ## 1815 122.0  70.0     No
    ## 1816 105.5  57.5     No
    ## 1817 127.0  74.0     No
    ## 1818 120.0  77.5     No
    ## 1819 159.0  95.0     No
    ## 1821 123.0  86.0     No
    ## 1822 130.0  80.0     No
    ## 1825 166.5 106.5     No
    ## 1826 123.0  75.0     No
    ## 1827 104.0  65.0     No
    ## 1828 202.0 132.0     No
    ## 1829 131.5  82.0     No
    ## 1830 110.0  70.0     No
    ## 1831 100.0  64.0     No
    ## 1832 116.0  79.5     No
    ## 1833 132.0  82.0     No
    ## 1834 111.0  79.0     No
    ## 1836 173.0 102.0     No
    ## 1837 177.0 124.0     No
    ## 1838 161.0  97.0     No
    ## 1839 159.0 105.0    Yes
    ## 1840 117.0  70.0     No
    ## 1841 154.0  82.0     No
    ## 1842 148.0  74.0     No
    ## 1843 123.0  90.0     No
    ## 1844 136.5  86.5     No
    ## 1845 135.0  83.0    Yes
    ## 1846 109.0  66.5     No
    ## 1848 104.0  64.0     No
    ## 1849 170.0 105.0     No
    ## 1850 150.5  98.0     No
    ## 1851 138.0  78.0     No
    ## 1853 129.0  87.0     No
    ## 1854 124.0  77.0     No
    ## 1856 103.5  60.0     No
    ## 1857 134.0  84.0     No
    ## 1858 158.0  93.0     No
    ## 1859 171.0  89.0     No
    ## 1860 134.0  84.0     No
    ## 1861 110.0  73.0     No
    ## 1864 147.0  88.0     No
    ## 1865 119.5  73.0     No
    ## 1866 150.0  99.0     No
    ## 1867 140.0  88.0     No
    ## 1868 132.0  94.0     No
    ## 1869 123.5  83.0     No
    ## 1870 120.0  70.0     No
    ## 1871 129.0  80.0     No
    ## 1872 119.0  80.0     No
    ## 1873 138.0  90.0     No
    ## 1874 155.5 100.5     No
    ## 1875 105.0  65.0     No
    ## 1876 132.0  82.0     No
    ## 1877 119.0  83.0     No
    ## 1878 129.0  83.5     No
    ## 1881 123.0  79.0     No
    ## 1882  97.0  64.0     No
    ## 1883 158.0  98.0    Yes
    ## 1884 155.0  89.0     No
    ## 1885 112.5  75.5     No
    ## 1886 105.5  67.5     No
    ## 1887 135.0  83.0     No
    ## 1888 140.0  91.0     No
    ## 1889 107.0  67.5     No
    ## 1891 135.0  80.0     No
    ## 1892 141.0  87.0     No
    ## 1893 132.5  85.0     No
    ## 1894 131.0  92.0     No
    ## 1895 135.0  88.0     No
    ## 1896 125.0  79.0     No
    ## 1897 136.5  89.0     No
    ## 1899 120.0  70.0     No
    ## 1900 158.0  70.0     No
    ## 1901 128.0  71.0     No
    ## 1902 130.0  83.0     No
    ## 1903 100.0  70.0     No
    ## 1904 144.0  79.0     No
    ## 1906 120.0  82.5     No
    ## 1907 129.0  72.0     No
    ## 1909 122.0  69.0     No
    ## 1910 141.0  92.0     No
    ## 1911 102.5  60.0     No
    ## 1912 124.5  99.0     No
    ## 1913 163.0  82.0     No
    ## 1914 115.0  79.5     No
    ## 1916 128.0  72.0     No
    ## 1917 154.0  94.0     No
    ## 1918 129.0  85.5     No
    ## 1919 102.0  60.0     No
    ## 1920 180.0 100.0     No
    ## 1921 115.0  70.0     No
    ## 1922 150.0  86.0     No
    ## 1923 154.0 110.0     No
    ## 1924 142.5 104.5     No
    ## 1925 121.0  78.0     No
    ## 1926 140.0  88.0     No
    ## 1927 107.0  73.0     No
    ## 1928 156.0  91.0     No
    ## 1929 127.0  83.0     No
    ## 1930 107.5  66.5     No
    ## 1932 121.5  74.5     No
    ## 1933 120.0  67.5     No
    ## 1934 126.5  84.0     No
    ## 1935 104.0  73.5     No
    ## 1936 145.0  87.0     No
    ## 1937 134.0  82.0     No
    ## 1938 105.0  72.5     No
    ## 1939 123.0  78.0     No
    ## 1940 195.0  90.0    Yes
    ## 1941 151.0  88.0     No
    ## 1942 120.0  67.0    Yes
    ## 1944 119.5  79.0     No
    ## 1945 167.0 114.0     No
    ## 1946 157.5  78.0     No
    ## 1949 137.0  91.0     No
    ## 1950 121.0  78.0     No
    ## 1951 122.0  78.0     No
    ## 1952 114.0  87.0     No
    ## 1953 155.0  93.5     No
    ## 1955 124.0  77.0     No
    ## 1956 119.0  82.0     No
    ## 1958 127.0  76.0     No
    ## 1959 128.5  82.0     No
    ## 1960 126.0  80.0     No
    ## 1961 130.0  81.0     No
    ## 1962 137.0  85.0     No
    ## 1963 119.0  75.5     No
    ## 1964 124.0  76.0     No
    ## 1965 127.5  91.5     No
    ## 1966 200.0 120.0     No
    ## 1968 147.5  87.5     No
    ## 1969 135.0  76.5     No
    ## 1970 126.0  76.0     No
    ## 1972 144.0  91.0     No
    ## 1974 113.0  70.0     No
    ## 1975 125.0  90.0     No
    ## 1976 146.0  82.0     No
    ## 1977 146.0  94.5     No
    ## 1978 132.0  81.0     No
    ## 1979 120.0  85.0     No
    ## 1980 152.0 104.0     No
    ## 1981 127.0  90.0     No
    ## 1982 136.5  92.0     No
    ## 1983 112.5  67.0     No
    ## 1984 135.5  90.0     No
    ## 1985 128.0  78.0     No
    ## 1986 113.0  77.0     No
    ## 1987 127.0  82.0     No
    ## 1988 137.5  87.0     No
    ## 1989 112.0  74.0     No
    ## 1990 141.0  78.0     No
    ## 1991 180.0 114.0    Yes
    ## 1992 107.0  76.0     No
    ## 1993 113.5  70.0     No
    ## 1994 145.0  85.0     No
    ## 1995 110.0  74.0     No
    ## 1996 133.0  78.0     No
    ## 1997 113.0  65.0     No
    ## 1998 119.0  73.5     No
    ## 1999 107.5  65.0     No
    ## 2000 126.0  73.0     No
    ## 2001 118.0  84.0     No
    ## 2003 108.0  73.0     No
    ## 2005 111.5  72.5     No
    ## 2006 127.5  70.0     No
    ## 2007 147.0  96.0     No
    ## 2008 170.0 110.0     No
    ## 2009 130.0  75.0     No
    ## 2010 132.0  82.0     No
    ## 2011 143.0  95.0     No
    ## 2012 116.5  77.5     No
    ## 2014 105.0  67.5     No
    ## 2015 146.0  85.0     No
    ## 2017 133.5  80.0     No
    ## 2018 129.0  91.5     No
    ## 2019 112.0  82.0     No
    ## 2020 158.0  90.0     No
    ## 2021 150.0  80.0     No
    ## 2023 114.0  68.0     No
    ## 2024 135.0  82.0     No
    ## 2025 111.0  79.5     No
    ## 2026 112.5  70.0     No
    ## 2027 173.0 112.0     No
    ## 2028 142.0  79.0     No
    ## 2029 129.0  88.0     No
    ## 2030 142.0  84.0     No
    ## 2032 118.0  84.0     No
    ## 2033 134.0  72.0     No
    ## 2034 155.0  92.5     No
    ## 2035 108.0  65.5     No
    ## 2036 137.5 101.5     No
    ## 2037 122.0  82.0     No
    ## 2038 119.0  81.0     No
    ## 2039 182.0  92.0    Yes
    ## 2040 150.0  92.0     No
    ## 2041 145.0  90.0     No
    ## 2042 145.0  94.0     No
    ## 2043 157.0  71.0     No
    ## 2045 134.0  84.0     No
    ## 2046 120.5  73.5     No
    ## 2047 128.0  83.5     No
    ## 2048 121.0  85.5     No
    ## 2049 175.0  85.0     No
    ## 2050 111.0  70.0     No
    ## 2052 111.0  85.0     No
    ## 2053 142.0  83.0     No
    ## 2054 119.5  69.5     No
    ## 2055 105.0  72.5     No
    ## 2056 125.0  83.0     No
    ## 2057 145.5  99.0     No
    ## 2058  97.0  67.0     No
    ## 2059 119.0  73.0     No
    ## 2060 139.0  80.0     No
    ## 2061 126.0  89.5     No
    ## 2062 112.0  75.0     No
    ## 2063 171.0  87.0     No
    ## 2064 139.0  88.0     No
    ## 2065 132.0  82.0     No
    ## 2067 116.0  86.0     No
    ## 2068 112.5  75.0     No
    ## 2069 151.0  95.0     No
    ## 2072 126.0  79.0     No
    ## 2073 108.0  82.0     No
    ## 2075 126.0  72.0     No
    ## 2077 116.0  74.0     No
    ## 2078 105.0  70.0     No
    ## 2084 115.5  65.0     No
    ## 2085 232.0 136.0     No
    ## 2086 128.5  73.5     No
    ## 2087 119.0  80.0     No
    ## 2088 122.5  80.0     No
    ## 2089 134.0  88.0     No
    ## 2090 173.0  96.0    Yes
    ## 2091 132.0  85.0    Yes
    ## 2092 111.5  77.0     No
    ## 2093  85.5  51.0     No
    ## 2095 132.0  81.0     No
    ## 2096 110.0  80.0     No
    ## 2097 134.0  80.0     No
    ## 2100 112.0  62.0     No
    ## 2101 102.0  66.5     No
    ## 2103 162.5 100.0     No
    ## 2104 191.0  81.0    Yes
    ## 2106 133.0  77.0     No
    ## 2108 168.0  94.0     No
    ## 2109 148.0  99.0     No
    ## 2110 118.0  80.0     No
    ## 2112 108.0  70.5     No
    ## 2113 128.0  82.0     No
    ## 2115 110.0  80.0     No
    ## 2116 126.0  91.0     No
    ## 2117 141.5  93.0     No
    ## 2119 142.0  91.0     No
    ## 2120 132.0  82.0     No
    ## 2122 150.5  97.0     No
    ## 2123 184.5  83.0     No
    ## 2124 175.0  80.0     No
    ## 2126  95.0  59.0     No
    ## 2127 122.0  74.5     No
    ## 2128 135.0  80.0     No
    ## 2129 145.0  92.0     No
    ## 2130 144.0  82.0     No
    ## 2131 100.0  80.0     No
    ## 2132 118.5  88.0     No
    ## 2133 113.0  79.0     No
    ## 2134 130.0  87.0     No
    ## 2135 150.0  93.0     No
    ## 2136 139.0  82.5     No
    ## 2137 103.0  72.5     No
    ## 2139 126.0  91.0     No
    ## 2140 100.0  70.0     No
    ## 2141 157.0  97.0     No
    ## 2142 156.0  95.0     No
    ## 2143 134.5  92.5     No
    ## 2144 130.0  87.0     No
    ## 2145 112.5  77.5     No
    ## 2147 116.0  70.0     No
    ## 2148 105.0  59.0     No
    ## 2151 122.0  81.5     No
    ## 2152 107.5  74.0     No
    ## 2153 120.0  80.0     No
    ## 2154 105.0  72.5     No
    ## 2155 130.0  80.0     No
    ## 2156 146.5  82.0     No
    ## 2157 126.0  78.0     No
    ## 2158 114.0  75.0     No
    ## 2159 112.5  75.0     No
    ## 2160 117.0  77.0     No
    ## 2162 138.5  90.0     No
    ## 2163 120.0  73.0     No
    ## 2164 122.0  74.0     No
    ## 2166 149.0  86.0     No
    ## 2167 120.0  83.5     No
    ## 2170 117.0  84.0     No
    ## 2171 117.5  72.5     No
    ## 2172 124.0  77.0     No
    ## 2173 127.0  79.0     No
    ## 2174 116.0  72.0     No
    ## 2176 112.0  73.5     No
    ## 2177 121.0  70.0     No
    ## 2178 145.0  86.0     No
    ## 2179  95.5  64.0     No
    ## 2181 117.0  78.0     No
    ## 2183 130.0  72.0     No
    ## 2184 152.0  98.0     No
    ## 2185 134.0  87.0     No
    ## 2187 155.0  84.0     No
    ## 2188 141.0  90.0     No
    ## 2189 171.0 120.0     No
    ## 2190 102.0  69.0     No
    ## 2191 113.0  69.0     No
    ## 2192 186.0 109.0    Yes
    ## 2194 123.0  74.0     No
    ## 2195 188.0 128.0     No
    ## 2196 161.0 104.0    Yes
    ## 2197 107.0  76.0     No
    ## 2198 124.0  81.0     No
    ## 2200 162.5  93.5     No
    ## 2201 104.0  72.0     No
    ## 2202 157.0  98.0     No
    ## 2203 114.0  80.0     No
    ## 2204 119.0  74.0     No
    ## 2205 123.0  94.5     No
    ## 2208 130.0  86.0     No
    ## 2209 120.0  83.0     No
    ## 2214 144.5  88.5     No
    ## 2215 126.5  76.0     No
    ## 2216 119.0  85.0     No
    ## 2217 107.5  75.0     No
    ## 2218 146.0  76.0     No
    ## 2219 113.0  73.0     No
    ## 2220 127.0  79.0     No
    ## 2222 146.0  98.5     No
    ## 2224 137.0  93.0     No
    ## 2225 104.0  61.0     No
    ## 2226 135.0  87.0     No
    ## 2227 167.5 102.5     No
    ## 2228 136.0  88.0     No
    ## 2229 154.0  96.0     No
    ## 2230 141.0  98.0     No
    ## 2231 146.0  77.0     No
    ## 2232 133.0  83.0     No
    ## 2233 134.0  85.0     No
    ## 2234 125.0  82.0     No
    ## 2235 205.0  92.5     No
    ## 2236 169.5 104.5     No
    ## 2238 125.5  80.0     No
    ## 2239 185.0 107.5     No
    ## 2240 120.5  76.0     No
    ## 2242  96.5  67.0     No
    ## 2243 152.0  97.0     No
    ## 2244 130.0  88.0     No
    ## 2245 132.0  85.5     No
    ## 2246 127.0  86.0     No
    ## 2247 129.0  94.0     No
    ## 2248 149.5  85.0     No
    ## 2249 110.0  69.0     No
    ## 2250 124.0  74.0     No
    ## 2251 148.0  92.0     No
    ## 2252 106.0  59.5     No
    ## 2253 126.0  87.0     No
    ## 2254 121.0  79.5     No
    ## 2255 110.5  69.0     No
    ## 2256 146.5  80.0     No
    ## 2257 116.0  69.0     No
    ## 2258 130.0  95.0     No
    ## 2259 110.0  78.0     No
    ## 2260 116.0  77.0     No
    ## 2261 129.0  76.5     No
    ## 2262 188.0  92.0     No
    ## 2263 115.0  72.0     No
    ## 2266 220.0  96.0     No
    ## 2267 116.0  79.0     No
    ## 2268 116.0  83.0     No
    ## 2269 124.0  82.0     No
    ## 2270 146.0  90.0     No
    ## 2271 127.0  90.0     No
    ## 2272 111.0  72.0     No
    ## 2273 134.0  83.0     No
    ## 2274 182.0 111.0     No
    ## 2275 118.0  79.5     No
    ## 2276 124.0  86.0     No
    ## 2279 150.0  88.0     No
    ## 2280 107.5  67.5     No
    ## 2282 124.0  79.0     No
    ## 2283 210.0 120.0     No
    ## 2284 106.0  64.0     No
    ## 2285 109.0  75.0     No
    ## 2287 120.0  74.0     No
    ## 2289 132.0  81.0     No
    ## 2290 130.0  80.0     No
    ## 2291 151.0  92.0     No
    ## 2292 184.0 102.0    Yes
    ## 2293 105.0  73.0     No
    ## 2294 130.0  91.0     No
    ## 2295 134.5  91.0     No
    ## 2296 115.0  80.0     No
    ## 2297 127.0  86.0     No
    ## 2298  97.0  62.0     No
    ## 2299 128.0  83.0     No
    ## 2300 124.0  81.0     No
    ## 2301 114.0  70.0     No
    ## 2302 193.0  63.0     No
    ## 2303 121.5  79.0     No
    ## 2304 120.0  80.0     No
    ## 2305 154.0  98.0     No
    ## 2306 128.0  79.0     No
    ## 2307 115.0  69.0     No
    ## 2308 120.0  75.0     No
    ## 2309 158.0  96.0     No
    ## 2310 120.0  80.0     No
    ## 2311 181.5 102.5     No
    ## 2312 109.0  73.0     No
    ## 2313 115.0  71.0     No
    ## 2315 130.0  75.0     No
    ## 2316 185.0 105.0     No
    ## 2317 150.0 101.0     No
    ## 2318 154.0  84.0     No
    ## 2319 139.0  80.0     No
    ## 2320 188.5 106.5     No
    ## 2321 127.5  72.5     No
    ## 2323 106.0  80.0     No
    ## 2324 132.5  85.0     No
    ## 2326 120.0  83.0     No
    ## 2327 122.5  68.5     No
    ## 2328 165.0  95.0     No
    ## 2329 116.0  75.0     No
    ## 2330 113.0  68.0     No
    ## 2331 122.5  80.0     No
    ## 2332 142.5  74.5     No
    ## 2333 170.0  92.0     No
    ## 2334 132.5  80.0     No
    ## 2335 118.0  76.0     No
    ## 2336 119.0  83.5     No
    ## 2337 110.0  79.0     No
    ## 2338 112.0  70.0     No
    ## 2339 107.5  68.0     No
    ## 2340 125.0  82.0     No
    ## 2341 111.5  67.0     No
    ## 2342 139.0  80.0     No
    ## 2343 182.5  97.5     No
    ## 2344 165.0  84.0     No
    ## 2345 145.0  92.5     No
    ## 2346 108.0  75.0     No
    ## 2347 138.0  97.0     No
    ## 2348 116.0  68.0     No
    ## 2349 138.0  72.0     No
    ## 2350 108.0  74.0     No
    ## 2351 108.0  62.0     No
    ## 2352 164.0 111.0     No
    ## 2353 148.0 108.0     No
    ## 2354 131.5  91.0     No
    ## 2355  96.0  67.0     No
    ## 2356 120.0  60.0     No
    ## 2358 114.0  73.5     No
    ## 2359 148.0 103.0     No
    ## 2360 140.0  86.0     No
    ## 2362 120.0  78.0     No
    ## 2363 130.5  90.0     No
    ## 2364 128.5  87.5     No
    ## 2365 116.0  67.0     No
    ## 2366 120.0  80.0     No
    ## 2367 100.0  70.0     No
    ## 2368 109.5  72.0     No
    ## 2371 132.0  86.0     No
    ## 2372 115.0  75.0     No
    ## 2373 139.0  80.0     No
    ## 2375 148.0  75.0     No
    ## 2378 133.0  96.0     No
    ## 2379 124.5  66.5     No
    ## 2381 134.0  88.0     No
    ## 2382 174.0  90.0    Yes
    ## 2383 110.0  70.0     No
    ## 2384 120.0  70.0     No
    ## 2385 192.0 105.0    Yes
    ## 2386 109.0  77.0     No
    ## 2387 138.0  79.0     No
    ## 2388 146.0  92.0     No
    ## 2389 176.5 115.0     No
    ## 2390 110.0  71.0     No
    ## 2391 101.0  69.0     No
    ## 2392 148.0  91.0     No
    ## 2393 110.0  78.0     No
    ## 2394 137.0  84.0     No
    ## 2395 154.5  83.0     No
    ## 2396 127.5  62.5     No
    ## 2397 110.5  69.0     No
    ## 2398 153.0 100.0     No
    ## 2399 177.0  97.0     No
    ## 2400 127.5  83.5     No
    ## 2401 119.0  76.0     No
    ## 2402 140.0  81.0     No
    ## 2404 132.0  92.0     No
    ## 2405 133.0  77.5     No
    ## 2406 160.0  96.0     No
    ## 2407 113.5  72.0     No
    ## 2408 149.0  95.0     No
    ## 2409 183.0 108.0    Yes
    ## 2410 129.0  86.0     No
    ## 2411 116.0  85.5     No
    ## 2412 123.0  92.0     No
    ## 2413 114.0  68.0     No
    ## 2414 199.0 114.0     No
    ## 2415 159.0  92.0     No
    ## 2416 171.0  97.0     No
    ## 2417 197.5 125.0     No
    ## 2418 190.0 100.0     No
    ## 2419 116.0  88.5     No
    ## 2420 116.0  81.0     No
    ## 2421 112.5  80.0     No
    ## 2422 120.0  83.5     No
    ## 2423 116.0  77.5     No
    ## 2426 174.0 100.0     No
    ## 2427 133.0  84.0     No
    ## 2428 107.5  72.0     No
    ## 2429 110.0  66.0     No
    ## 2430 133.0  87.0     No
    ## 2431 116.0  83.0     No
    ## 2432 195.0 108.0     No
    ## 2433 140.0  85.0     No
    ## 2434 145.0  85.0     No
    ## 2435 135.0  75.0     No
    ## 2436 101.0  71.0     No
    ## 2437 143.0  81.0     No
    ## 2438 145.0  88.0     No
    ## 2439 122.0  78.0     No
    ## 2440 111.0  71.0     No
    ## 2441 146.0  78.0     No
    ## 2442 122.0  87.5     No
    ## 2443 120.0  80.0     No
    ## 2444 150.0  93.0     No
    ## 2445 180.0 100.0     No
    ## 2446 109.0  70.0     No
    ## 2447 115.5  79.0     No
    ## 2448 121.0  78.5     No
    ## 2449 132.0  85.0     No
    ## 2453 102.0  67.0     No
    ## 2454 110.0  72.5     No
    ## 2455 115.0  80.0     No
    ## 2456 136.5  85.0     No
    ## 2458 118.0  72.0     No
    ## 2459 144.0  79.0     No
    ## 2460 175.0  78.0     No
    ## 2461 120.0  81.0     No
    ## 2462 146.0 106.0     No
    ## 2463 150.0 101.0     No
    ## 2464 132.0  81.0     No
    ## 2465 125.0  75.0     No
    ## 2467 139.0  92.5     No
    ## 2468 165.0 106.0     No
    ## 2469 126.0  86.0     No
    ## 2471 105.0  72.5     No
    ## 2472 128.0  81.5     No
    ## 2473 129.0  74.0     No
    ## 2474 132.0  96.0     No
    ## 2475 128.0  83.0     No
    ## 2476 113.0  81.0     No
    ## 2477 131.0  96.0     No
    ## 2478 134.5  84.0     No
    ## 2479 111.0  70.0     No
    ## 2480 118.0  84.0     No
    ## 2482 150.0  94.0     No
    ## 2483 142.5  87.0     No
    ## 2484 148.0  92.5     No
    ## 2485 146.0  80.0     No
    ## 2486 176.0  87.0     No
    ## 2487 158.5  94.5     No
    ## 2489 115.0  75.0     No
    ## 2490 144.0  88.0     No
    ## 2491 133.0  87.0     No
    ## 2493 137.0  80.0     No
    ## 2494 135.0  86.0     No
    ## 2495 145.0  72.5     No
    ## 2496 166.0  96.0     No
    ## 2497 136.0  96.0     No
    ## 2498 121.0  88.0     No
    ## 2499 137.0  88.0     No
    ## 2500 146.0  91.0     No
    ## 2501 132.5  92.0     No
    ## 2502 142.0  85.0     No
    ## 2503 114.0  76.0     No
    ## 2504 108.0  62.0     No
    ## 2505 129.0  89.0     No
    ## 2506 112.0  80.0     No
    ## 2507 113.0  61.0     No
    ## 2508 137.0  95.0     No
    ## 2509 137.5  88.5     No
    ## 2510 102.0  61.0     No
    ## 2511 208.0 104.0     No
    ## 2512 137.0  83.5     No
    ## 2513 122.0  95.0     No
    ## 2514 118.5  76.0     No
    ## 2515 110.0  68.0     No
    ## 2516 136.0  88.0     No
    ## 2518 124.0  78.0     No
    ## 2519 105.0  69.0     No
    ## 2520 120.0  87.0     No
    ## 2521 130.0  85.0     No
    ## 2522 136.5  97.0     No
    ## 2523 125.0  82.0     No
    ## 2524 123.0  69.0     No
    ## 2525 130.0  87.0     No
    ## 2526 123.5  78.0     No
    ## 2527 150.0  89.0     No
    ## 2529 118.0  70.0     No
    ## 2530 102.5  64.5     No
    ## 2531 113.0  79.0     No
    ## 2532 129.0  81.0     No
    ## 2533  94.0  62.0     No
    ## 2534 163.0 102.0     No
    ## 2535 210.0 130.0     No
    ## 2537 132.0  91.0     No
    ## 2538 112.5  75.0     No
    ## 2539 184.5 110.5     No
    ## 2540 141.0  82.5     No
    ## 2541  92.5  70.0     No
    ## 2542 126.0  85.0     No
    ## 2543 172.0 105.0     No
    ## 2544 120.0  82.0     No
    ## 2545 128.0  87.0     No
    ## 2546 110.0  73.0     No
    ## 2547 102.0  68.0     No
    ## 2548 115.0  79.0     No
    ## 2551 132.0  76.0     No
    ## 2552 128.0  74.0     No
    ## 2553 125.0  88.0     No
    ## 2554 127.5  77.5     No
    ## 2555 145.0  89.0     No
    ## 2556 165.0  99.0     No
    ## 2557 144.0  98.0     No
    ## 2558 127.0  82.5     No
    ## 2559 126.0  82.5     No
    ## 2560 117.5  77.5     No
    ## 2561 160.0  97.5     No
    ## 2562 136.5  87.0     No
    ## 2563 138.0  91.0     No
    ## 2564 127.0  68.5     No
    ## 2565 103.0  76.5     No
    ## 2566 172.0  95.0    Yes
    ## 2567 129.5  82.5     No
    ## 2568 176.0 113.0     No
    ## 2569 125.0  94.0     No
    ## 2570 116.0  82.0     No
    ## 2571 162.0  85.0     No
    ## 2572 102.0  69.0     No
    ## 2573 136.0  86.0     No
    ## 2574 124.0  72.5     No
    ## 2575 127.0  78.0     No
    ## 2576 115.0  76.0     No
    ## 2577 109.0  75.0     No
    ## 2578 128.0  91.0     No
    ## 2579 202.5  85.0     No
    ## 2580 120.5  78.0     No
    ## 2581 140.0 100.0     No
    ## 2583 136.5  85.0     No
    ## 2584 145.0  88.0    Yes
    ## 2585 112.0  67.0     No
    ## 2586 102.0  72.0     No
    ## 2587 156.0  69.0     No
    ## 2588 140.0  97.5    Yes
    ## 2590 155.0  85.0    Yes
    ## 2591 135.0  90.0     No
    ## 2592 120.0  72.0     No
    ## 2593 136.0  94.0     No
    ## 2594 139.0  74.0     No
    ## 2595 154.0  97.0    Yes
    ## 2597 142.5  83.5     No
    ## 2598 166.0 101.0     No
    ## 2599 128.0  85.0     No
    ## 2600 122.5  78.5     No
    ## 2601 146.0  92.0     No
    ## 2602 136.5  87.5     No
    ## 2603 112.0  66.0     No
    ## 2604  98.0  74.0     No
    ## 2605 130.0  84.0     No
    ## 2606 122.0  80.0     No
    ## 2607 112.5  68.0     No
    ## 2608 131.0  84.0     No
    ## 2609 150.0  95.5     No
    ## 2610 141.0  84.5     No
    ## 2611 166.0  85.0     No
    ## 2612 120.0  81.0     No
    ## 2613 141.0  77.5     No
    ## 2614 120.0  81.0     No
    ## 2615 175.0  82.0    Yes
    ## 2616 146.0  98.5     No
    ## 2617 164.0 104.0     No
    ## 2618 138.0  89.0     No
    ## 2619 104.0  72.0     No
    ## 2620 130.0  80.0     No
    ## 2621 137.5  72.5     No
    ## 2622 121.0  74.0     No
    ## 2623 115.0  70.0     No
    ## 2624 140.0  84.0     No
    ## 2625 145.0  81.0     No
    ## 2626 141.0 105.0     No
    ## 2627 124.0  85.0     No
    ## 2628 124.0  78.0     No
    ## 2629 158.0  94.0     No
    ## 2630 117.0  83.0     No
    ## 2632 159.0 102.0    Yes
    ## 2633  96.0  59.0     No
    ## 2634 169.0 117.0     No
    ## 2635 102.5  65.0     No
    ## 2636 149.0  82.0    Yes
    ## 2637 119.0  80.0     No
    ## 2638 100.0  61.5     No
    ## 2639 133.5  80.0     No
    ## 2640 113.0  70.0     No
    ## 2641 114.0  82.0     No
    ## 2643 125.0  79.0     No
    ## 2644 128.0  82.0     No
    ## 2645 110.0  79.0     No
    ## 2646 110.0  71.0     No
    ## 2649 120.0  83.5     No
    ## 2650 110.0  70.0     No
    ## 2651 107.5  80.0     No
    ## 2652 122.0  81.0     No
    ## 2655 112.5  85.0     No
    ## 2656 149.5  93.0     No
    ## 2657 117.5  80.0     No
    ## 2658 159.0  95.0     No
    ## 2659 121.0  83.0     No
    ## 2660 113.0  82.5     No
    ## 2661 140.0  86.0     No
    ## 2662 126.0  73.0     No
    ## 2663 122.0  86.0     No
    ## 2664 110.0  71.0     No
    ## 2665 155.0 100.0     No
    ## 2667 144.0  90.0     No
    ## 2668 122.0  81.0     No
    ## 2669 123.0  80.0     No
    ## 2670 119.0  78.5     No
    ## 2671 182.5  97.0     No
    ## 2672 132.5  85.5     No
    ## 2673 166.0  90.0     No
    ## 2674 126.0  84.0     No
    ## 2675 137.5  94.0     No
    ## 2676 131.5  89.0     No
    ## 2677 123.0  76.5     No
    ## 2679 141.0  81.5     No
    ## 2680 147.5  87.5     No
    ## 2681 154.5 104.0     No
    ## 2682 125.0  85.0     No
    ## 2683 105.0  71.0     No
    ## 2684 122.0  87.0    Yes
    ## 2685 135.0  85.0     No
    ## 2686 142.5  95.0     No
    ## 2687 119.0  72.0     No
    ## 2688 146.0  98.0     No
    ## 2689 160.0  90.0    Yes
    ## 2690 109.0  75.0     No
    ## 2692 112.0  78.0     No
    ## 2693 180.0 108.0     No
    ## 2695 115.0  80.0     No
    ## 2696 139.0  96.0     No
    ## 2698 147.0 102.0     No
    ## 2699 175.0  85.0     No
    ## 2700 112.5  85.0     No
    ## 2701 123.0  77.5     No
    ## 2702 116.5  83.0     No
    ## 2703 112.0  77.0     No
    ## 2704 131.0  71.0     No
    ## 2705 122.0  78.0     No
    ## 2706 115.0  84.0     No
    ## 2708 129.0  86.0     No
    ## 2710 165.0  99.0     No
    ## 2711 102.5  72.5     No
    ## 2712 142.0  92.0     No
    ## 2714 114.0  85.0     No
    ## 2715 119.0  75.0     No
    ## 2716 108.0  74.0     No
    ## 2718 105.0  70.0     No
    ## 2719 177.0  75.0     No
    ## 2721 107.5  66.0     No
    ## 2722 155.0  99.0     No
    ## 2723 144.0  82.0     No
    ## 2724 135.0  95.0     No
    ## 2725 182.0 102.0     No
    ## 2726 129.0  84.0     No
    ## 2727 130.0  82.5     No
    ## 2728 125.0  72.5     No
    ## 2729 137.5  82.5     No
    ## 2730 162.0  91.0     No
    ## 2731 170.0 104.0     No
    ## 2732 133.0  77.0     No
    ## 2735 110.0  72.0     No
    ## 2736 128.0  89.0     No
    ## 2738 129.0  73.0     No
    ## 2739 166.0  71.0     No
    ## 2741 103.0  78.0     No
    ## 2742 141.0  86.0    Yes
    ## 2743 139.5  87.0     No
    ## 2744 116.0  72.0     No
    ## 2745 191.0 106.0     No
    ## 2746 134.0  88.0     No
    ## 2747 113.0  84.0     No
    ## 2748 110.0  70.0     No
    ## 2749 133.0  89.0     No
    ## 2750 152.5  85.0     No
    ## 2752 132.0  80.0     No
    ## 2754 165.0  90.0     No
    ## 2755 116.0  69.0     No
    ## 2756 123.0  70.0     No
    ## 2757 139.0  67.0     No
    ## 2758 131.0  79.0     No
    ## 2759 121.5  86.5     No
    ## 2760 182.5 103.0     No
    ## 2761 134.0  87.0     No
    ## 2762 118.0  83.0     No
    ## 2763 105.0  70.0     No
    ## 2764 114.0  72.5     No
    ## 2766 127.0  79.5     No
    ## 2768 121.0  61.0     No
    ## 2769 147.0  86.5     No
    ## 2770 125.0  85.0     No
    ## 2771 118.0  79.0     No
    ## 2773 160.0  82.0     No
    ## 2774 120.0  73.5     No
    ## 2775 114.0  78.0     No
    ## 2776 152.0  88.0     No
    ## 2777 147.5  93.0     No
    ## 2778 143.5  90.0     No
    ## 2779 131.0  74.0     No
    ## 2780 122.5  77.5     No
    ## 2781 111.0  81.0     No
    ## 2783 147.0  98.0    Yes
    ## 2784 119.0  80.0     No
    ## 2785 190.0 130.0     No
    ## 2787 145.0  88.0     No
    ## 2788 119.0  82.0     No
    ## 2789 120.0  75.0     No
    ## 2790 151.5  79.0     No
    ## 2791 127.5  90.0     No
    ## 2792 125.5  80.5     No
    ## 2793  83.5  58.0     No
    ## 2794  96.0  67.0     No
    ## 2795 152.0  73.0     No
    ## 2796 161.0  90.0    Yes
    ## 2798 159.0  87.0    Yes
    ## 2799 144.0  74.0     No
    ## 2800 142.0  84.0     No
    ## 2801 137.5  81.0     No
    ## 2802 115.0  75.0     No
    ## 2803 152.5  88.0     No
    ## 2804 131.5  81.0     No
    ## 2805 115.0  69.0     No
    ## 2806 128.0  87.0     No
    ## 2807  98.0  73.0     No
    ## 2808 118.0  80.0     No
    ## 2809 127.0  68.0     No
    ## 2810 146.5  71.0     No
    ## 2811 150.0  94.0     No
    ## 2812 210.0 135.0     No
    ## 2813 117.5  77.5     No
    ## 2814 131.5  83.0     No
    ## 2815 171.0  84.0     No
    ## 2816 158.0 109.0     No
    ## 2817 121.0  75.0     No
    ## 2819 116.0  69.0     No
    ## 2820 110.0  80.0     No
    ## 2821 115.0  63.5     No
    ## 2822 122.0  76.0     No
    ## 2823 115.0  81.0     No
    ## 2826 150.0  89.0     No
    ## 2827 148.0  95.0     No
    ## 2829 124.0  75.0     No
    ## 2830  99.0  60.0     No
    ## 2831 136.5  84.0     No
    ## 2832 122.0  83.0     No
    ## 2833 128.0  90.0     No
    ## 2834 150.0  95.0     No
    ## 2835 144.0  85.0    Yes
    ## 2836 106.5  75.0     No
    ## 2837 132.0  90.0     No
    ## 2838 121.0  82.0     No
    ## 2841 113.0  74.0     No
    ## 2842 176.0  95.0     No
    ## 2843 155.0 105.0     No
    ## 2844 130.0  80.0     No
    ## 2845 121.0  78.0     No
    ## 2846 152.0  95.0     No
    ## 2847 143.0  76.0     No
    ## 2848 148.0  70.0     No
    ## 2850 101.0  68.0     No
    ## 2851 142.0  80.0     No
    ## 2852 123.0  77.0     No
    ## 2853 130.0  77.5     No
    ## 2854 177.0 101.0     No
    ## 2855 103.0  73.0     No
    ## 2856 116.0  68.0     No
    ## 2859 140.5  92.5     No
    ## 2860 120.0  84.0     No
    ## 2861 165.0 112.0     No
    ## 2862 140.0  92.0     No
    ## 2863 141.0 115.0    Yes
    ## 2864 170.5 100.5     No
    ## 2865 127.0  84.0     No
    ## 2866 142.0  82.5     No
    ## 2867 119.0  72.0     No
    ## 2869 108.0  76.0     No
    ## 2870 117.5  83.5     No
    ## 2871 132.0  85.0     No
    ## 2872 131.0  80.0     No
    ## 2873 107.0  74.0     No
    ## 2874 111.0  76.0     No
    ## 2875 156.0  86.0     No
    ## 2876 121.0  82.0     No
    ## 2877 108.0  73.0     No
    ## 2879 126.5  67.0     No
    ## 2880 132.0  94.0     No
    ## 2881 120.0  77.5     No
    ## 2883 107.0  68.0     No
    ## 2884 114.0  74.0     No
    ## 2885 175.0 114.0    Yes
    ## 2887 150.0  85.0     No
    ## 2888 138.0  78.0     No
    ## 2889 113.5  65.5     No
    ## 2890 132.5  87.0     No
    ## 2891 109.0  73.0     No
    ## 2892 126.0  80.0     No
    ## 2894 134.0  87.0     No
    ## 2895 107.0  65.0     No
    ## 2896 165.0  94.5     No
    ## 2897 111.0  60.0     No
    ## 2898 123.0  82.0     No
    ## 2899  93.0  71.0     No
    ## 2900 123.5  79.0     No
    ## 2902 115.0  72.0     No
    ## 2904 101.0  71.0     No
    ## 2905 171.0 118.0     No
    ## 2906 132.0  87.0     No
    ## 2907 118.5  72.0     No
    ## 2908 133.0  88.0     No
    ## 2909 181.0  97.5    Yes
    ## 2910 143.0  93.0     No
    ## 2911 102.5  71.0     No
    ## 2912 142.5 100.0     No
    ## 2913 126.0  85.5     No
    ## 2915 113.5  61.0     No
    ## 2916 116.5  72.5     No
    ## 2917 128.0  76.0     No
    ## 2918 138.0  98.0     No
    ## 2919 148.0  94.0     No
    ## 2920 156.0 100.0     No
    ## 2921 142.5  83.5     No
    ## 2923 148.0  98.0     No
    ## 2924 138.0  78.0     No
    ## 2925 149.0 103.5     No
    ## 2926 166.0  98.0    Yes
    ## 2927 141.0  82.5     No
    ## 2928 111.0  80.0     No
    ## 2929 155.0  71.0     No
    ## 2931 197.0 109.0     No
    ## 2932 164.0  94.0     No
    ## 2933 156.0  91.0     No
    ## 2934 116.0  74.0     No
    ## 2936 146.0  89.0     No
    ## 2937 118.0  85.0     No
    ## 2938 110.0  70.0     No
    ## 2939 162.5  93.5    Yes
    ## 2940 122.5  67.5     No
    ## 2941 144.0  79.0     No
    ## 2942 105.0  60.0     No
    ## 2943 121.5  76.5     No
    ## 2944 114.0  80.0     No
    ## 2945 119.0  69.0     No
    ## 2946 102.5  69.0     No
    ## 2947 116.0  86.0     No
    ## 2948 123.0  86.0     No
    ## 2949 114.0  80.0     No
    ## 2950 181.0  74.0     No
    ## 2951 132.0  80.0     No
    ## 2952 129.0  90.0     No
    ## 2953 145.0  88.0     No
    ## 2954 116.0  81.0     No
    ## 2955 107.5  70.0     No
    ## 2957 117.5  77.5     No
    ## 2958 139.0  82.0     No
    ## 2959 135.0  87.5     No
    ## 2960 131.0  90.0     No
    ## 2961 108.0  73.0     No
    ## 2962 138.0  76.0     No
    ## 2963 138.0  85.0     No
    ## 2964 102.0  59.0     No
    ## 2965 125.0  84.0     No
    ## 2966 126.0  93.0     No
    ## 2968 198.0 108.0    Yes
    ## 2969 123.0  82.0     No
    ## 2970 135.0  97.5    Yes
    ## 2971 134.0  79.0     No
    ## 2972 163.0  87.0     No
    ## 2973 119.0  82.0     No
    ## 2974 115.0  75.0     No
    ## 2975 141.0  85.0     No
    ## 2976 114.5  76.0     No
    ## 2977 124.0  83.0     No
    ## 2980 163.0  83.0     No
    ## 2981 113.5  80.0     No
    ## 2982 131.0  79.0     No
    ## 2983 135.0  82.0     No
    ## 2984 150.0  88.0     No
    ## 2985 126.0  76.0     No
    ## 2986 107.0  82.0     No
    ## 2987 111.0  72.0     No
    ## 2988 135.0  77.5     No
    ## 2989 190.0 110.0     No
    ## 2990 107.0  82.0     No
    ## 2992 136.0  85.5     No
    ## 2993 106.0  75.0     No
    ## 2994 112.0  71.0     No
    ## 2995 125.0  80.0     No
    ## 2996 122.5  82.0     No
    ## 2997 175.5 113.0    Yes
    ## 2998  96.5  72.5     No
    ## 2999 132.0  95.0     No
    ## 3000 114.0  64.0     No
    ## 3001 130.0  79.0     No
    ## 3002 119.0  86.0     No
    ## 3003 172.0  82.0     No
    ## 3004 118.0  71.0     No
    ## 3006 149.0  88.0     No
    ## 3007 122.0  85.0     No
    ## 3008 126.5  82.0     No
    ## 3009 142.5  86.5     No
    ## 3010 160.0  99.0     No
    ## 3011 103.0  64.5     No
    ## 3012 135.0  87.0     No
    ## 3013 128.0  89.0     No
    ## 3014 140.0  88.0     No
    ## 3016 181.0 107.0    Yes
    ## 3017 163.5 102.0     No
    ## 3019 122.0  78.0     No
    ## 3020 124.0  73.0     No
    ## 3021 166.0  89.0     No
    ## 3022 102.0  66.5     No
    ## 3023 103.0  66.0     No
    ## 3024 130.5  82.0     No
    ## 3025 113.0  62.0     No
    ## 3026 102.5  66.5     No
    ## 3027  90.0  70.0     No
    ## 3029 109.0  70.0     No
    ## 3030 124.0  72.5     No
    ## 3031 132.5  73.0     No
    ## 3032 171.0  95.0     No
    ## 3033 112.5  73.5     No
    ## 3036 120.0  62.0     No
    ## 3039 133.5  89.5     No
    ## 3040 149.0  89.0     No
    ## 3041 127.5  76.0     No
    ## 3042 204.0  96.0    Yes
    ## 3043 129.0  79.0     No
    ## 3045 125.0  86.0     No
    ## 3046 120.0  76.0     No
    ## 3047 118.0  70.5     No
    ## 3048 132.5  87.0     No
    ## 3049 120.0  72.0     No
    ## 3050 142.0  89.0     No
    ## 3051 122.0  82.5     No
    ## 3052 137.0  90.0     No
    ## 3053 151.5 103.0     No
    ## 3054 134.0  89.0     No
    ## 3055 132.5  69.0     No
    ## 3056 112.0  78.0     No
    ## 3057 107.5  73.5     No
    ## 3059 127.0  80.0     No
    ## 3060 140.5  89.0     No
    ## 3063 207.5 118.0     No
    ## 3064 115.0  75.0     No
    ## 3065 115.0  81.0     No
    ## 3066 127.0  80.0     No
    ## 3067 129.5  85.5     No
    ## 3068 102.0  59.5     No
    ## 3069 181.0 101.0     No
    ## 3070 191.0  97.0     No
    ## 3071 138.0  88.0     No
    ## 3072 145.0  91.0     No
    ## 3073 136.0  87.0     No
    ## 3074 110.0  79.0     No
    ## 3076 113.0  83.0     No
    ## 3079 150.5  94.0     No
    ## 3081 132.0  88.0     No
    ## 3082 116.0  72.0     No
    ## 3084 104.0  73.0     No
    ## 3085 121.0  70.0     No
    ## 3086 195.0 110.0     No
    ## 3087 130.5  86.0     No
    ## 3088 137.5  80.0     No
    ## 3089 110.0  77.0     No
    ## 3091 131.0  81.5     No
    ## 3092 147.5 100.0     No
    ## 3093 165.0  88.0     No
    ## 3094 142.0  88.0     No
    ## 3095 177.0 111.0     No
    ## 3096 108.0  75.0     No
    ## 3097 124.0  72.0     No
    ## 3098 145.0  67.0     No
    ## 3099 108.0  74.0     No
    ## 3100 155.0  92.5     No
    ## 3102 105.0  72.5     No
    ## 3103 147.5  97.5    Yes
    ## 3104 127.5  80.0     No
    ## 3106 136.5  87.0     No
    ## 3107 134.0  85.5     No
    ## 3108 174.5 101.5     No
    ## 3109 126.0  88.0     No
    ## 3110 126.0  77.0     No
    ## 3111 133.5  86.0     No
    ## 3112 125.0  86.0     No
    ## 3113 148.5  88.0     No
    ## 3114 105.0  85.0     No
    ## 3115 158.0  88.0     No
    ## 3116 130.0  84.0     No
    ## 3117 198.0 106.0    Yes
    ## 3118 152.5 105.0     No
    ## 3119 153.0  85.0    Yes
    ## 3122 130.0  90.0     No
    ## 3123 107.5  75.0     No
    ## 3124 150.0  84.0     No
    ## 3125 160.0 105.0     No
    ## 3126  98.0  53.0     No
    ## 3127 122.0  83.0     No
    ## 3129 126.0  75.0     No
    ## 3130 121.0  85.5     No
    ## 3131 103.0  67.5     No
    ## 3132 135.0  88.0     No
    ## 3133 127.0  80.0     No
    ## 3134 125.0  86.0     No
    ## 3136 125.0  80.0     No
    ## 3137 150.0 105.0     No
    ## 3138  94.0  62.0     No
    ## 3139 163.5  87.0     No
    ## 3140 126.0  80.0     No
    ## 3142 132.0  86.0     No
    ## 3143 118.5  81.0     No
    ## 3144 105.0  65.0     No
    ## 3145 131.0  87.0     No
    ## 3147 118.0  80.0     No
    ## 3148 131.0  94.0     No
    ## 3149 106.0  64.0     No
    ## 3150 126.0  86.0     No
    ## 3151 122.0  88.5     No
    ## 3152 125.0  79.0     No
    ## 3153 147.0  89.0     No
    ## 3154 148.0 103.0     No
    ## 3155 130.0  86.0     No
    ## 3156 120.5  80.5     No
    ## 3157 137.0  79.0     No
    ## 3158 109.0  70.0     No
    ## 3159 122.0  82.0     No
    ## 3160 197.0  91.0     No
    ## 3161 147.5  97.0     No
    ## 3162 115.0  71.0     No
    ## 3163 131.0  84.0     No
    ## 3164 125.0  80.0     No
    ## 3165 186.5  99.0     No
    ## 3166 141.0  76.0     No
    ## 3167 119.0  80.0     No
    ## 3170 108.5  73.0     No
    ## 3171 112.5  77.5     No
    ## 3172 142.0  93.0     No
    ## 3173 131.0  81.0     No
    ## 3174 114.0  79.0     No
    ## 3175 121.0  77.0     No
    ## 3176 149.0  86.0     No
    ## 3177 124.5  80.0     No
    ## 3178 126.0  80.0     No
    ## 3179 123.0  78.0     No
    ## 3180 107.0  77.0     No
    ## 3183 147.0  98.0     No
    ## 3185 130.0  72.5     No
    ## 3186 153.0  89.0     No
    ## 3188 128.5  75.0     No
    ## 3189 127.5  83.0     No
    ## 3190 193.0 109.0    Yes
    ## 3192 106.0  65.0     No
    ## 3193 164.0 113.0     No
    ## 3194 154.0  99.0     No
    ## 3195 123.0  83.0     No
    ## 3196 131.0  92.0     No
    ## 3197 164.0 120.0    Yes
    ## 3198 148.5 100.0     No
    ## 3200 127.5  81.5     No
    ## 3201 136.0  75.0     No
    ## 3202 130.0  86.5     No
    ## 3203 116.0  64.0     No
    ## 3204 130.0  82.0     No
    ## 3205 114.0  75.0     No
    ## 3206 138.5  77.5     No
    ## 3207 123.0  72.0     No
    ## 3208 127.0  83.0     No
    ## 3209  98.0  67.0     No
    ## 3210 133.0  92.0     No
    ## 3211 196.0 101.0     No
    ## 3212 105.0  70.0     No
    ## 3213 139.0  96.0     No
    ## 3214 111.0  73.0     No
    ## 3215 117.5  90.0     No
    ## 3216 141.0  92.0     No
    ## 3217 158.0  80.0     No
    ## 3218 111.0  73.0     No
    ## 3219 127.5  92.0     No
    ## 3220 147.5  92.0     No
    ## 3221 131.5  82.5     No
    ## 3222 152.0  93.0     No
    ## 3223 127.0  69.0     No
    ## 3224 122.5  77.0     No
    ## 3225 159.0  90.0    Yes
    ## 3228 152.0  74.0     No
    ## 3229 133.0  89.0     No
    ## 3230 141.5 108.5     No
    ## 3231 172.0 111.0     No
    ## 3232 114.0  80.0     No
    ## 3233 139.0  83.0     No
    ## 3234 115.0  78.0     No
    ## 3236 127.0  72.0     No
    ## 3237 146.0  87.0     No
    ## 3238 150.5  95.0     No
    ## 3239 162.5  99.5     No
    ## 3240 115.0  71.0     No
    ## 3241 120.0  79.5     No
    ## 3243 132.0  88.0     No
    ## 3244 130.0  73.5     No
    ## 3245 121.0  82.0     No
    ## 3246 124.0  79.0     No
    ## 3249 163.5  97.0     No
    ## 3251 117.5  85.0     No
    ## 3252 132.0  86.0     No
    ## 3253 129.5  83.0     No
    ## 3254 163.0  97.0     No
    ## 3256 116.0  86.0     No
    ## 3257 114.5  70.0     No
    ## 3258 179.5  97.0     No
    ## 3259 125.0  87.0     No
    ## 3261 142.5  80.0     No
    ## 3262 126.0  80.0     No
    ## 3263 120.0  78.0     No
    ## 3264 126.0  76.0     No
    ## 3265 134.0  84.0     No
    ## 3266 118.5  74.5     No
    ## 3267 125.0  76.0     No
    ## 3268 160.0  90.0     No
    ## 3269 155.0 105.0     No
    ## 3270 133.0  82.0     No
    ## 3271 137.0  81.0     No
    ## 3272 144.0  91.0     No
    ## 3273 177.0 101.0     No
    ## 3274 136.0  84.0    Yes
    ## 3275 114.0  83.0     No
    ## 3276 138.0  99.0     No
    ## 3277 140.0  94.0     No
    ## 3278 116.0  79.0     No
    ## 3279 133.0  94.0     No
    ## 3280 140.0  78.0     No
    ## 3281 170.0  89.0    Yes
    ## 3282 124.0  72.5     No
    ## 3284 136.0  84.0     No
    ## 3287 151.0 106.5     No
    ## 3288 132.0  86.0     No
    ## 3289 107.0  73.0     No
    ## 3291 125.0  83.0     No
    ## 3292 131.0  87.0     No
    ## 3293 122.0  82.0     No
    ## 3294 123.5  77.0     No
    ## 3297 122.5  76.5     No
    ## 3298 139.0  74.0     No
    ## 3300 154.0  92.0     No
    ## 3301 103.0  71.0     No
    ## 3302 212.0 116.0     No
    ## 3303 123.0  73.0     No
    ## 3305 153.0 100.0     No
    ## 3306 140.0  87.0     No
    ## 3307 108.0  72.0     No
    ## 3308 167.0 109.0     No
    ## 3309 132.0  80.0     No
    ## 3310 115.0  81.0     No
    ## 3311 115.0  81.0     No
    ## 3312 122.0  76.5     No
    ## 3313 140.5  89.0     No
    ## 3314 107.0  73.5     No
    ## 3315 117.5  76.0     No
    ## 3316 136.5  81.0     No
    ## 3317 113.0  81.0     No
    ## 3318 162.0  90.0     No
    ## 3319 121.0  72.0     No
    ## 3320 131.0  93.0     No
    ## 3321 119.5  87.0     No
    ## 3322 128.0  84.0     No
    ## 3323 164.0  98.0     No
    ## 3324 152.0  89.0     No
    ## 3325 113.0  81.0     No
    ## 3326 117.5  72.0     No
    ## 3327 151.5  95.0     No
    ## 3328 120.0  81.0     No
    ## 3329 126.0  79.0     No
    ## 3330 148.0  90.0     No
    ## 3331 120.0  66.5     No
    ## 3334  95.5  70.0     No
    ## 3335 116.0  81.5     No
    ## 3336 105.0  69.0     No
    ## 3337 140.0  94.0     No
    ## 3338 124.0  80.0     No
    ## 3339 130.0  87.0     No
    ## 3340 132.0  80.0     No
    ## 3341 104.0  74.0     No
    ## 3342 142.0  91.0     No
    ## 3343 111.0  72.0     No
    ## 3344 120.0  80.0     No
    ## 3345 128.0  80.0     No
    ## 3346 144.0  92.0     No
    ## 3347 105.0  86.0     No
    ## 3348 140.0  88.0     No
    ## 3349 122.0  81.0     No
    ## 3350 115.0  83.0     No
    ## 3352 137.5  88.5     No
    ## 3353 127.5  82.5     No
    ## 3354 124.0  80.0     No
    ## 3355 170.0 107.0     No
    ## 3356 138.5  85.5     No
    ## 3357 154.0  91.0     No
    ## 3358 170.0 107.5     No
    ## 3359 199.5 107.0     No
    ## 3360 141.0  84.0     No
    ## 3361 161.0 103.0     No
    ## 3362 126.0  85.0     No
    ## 3363 111.0  71.0     No
    ## 3364 170.0  81.0     No
    ## 3366 152.0 102.0    Yes
    ## 3367 109.0  69.0     No
    ## 3368 128.0  82.0     No
    ## 3369 123.0  70.0     No
    ## 3370 110.0  68.0    Yes
    ## 3371 124.0  92.0     No
    ## 3372 152.0  92.0     No
    ## 3374 120.0  78.0     No
    ## 3375 108.0  73.0     No
    ## 3376 132.0  81.0     No
    ## 3377 136.5  84.0     No
    ## 3378 165.0  85.0     No
    ## 3379 128.5  73.0     No
    ## 3380 193.0  95.0     No
    ## 3382 119.0  80.0     No
    ## 3383 117.0  77.0     No
    ## 3384 115.0  80.0     No
    ## 3385 134.0  93.0     No
    ## 3386 131.0  85.0     No
    ## 3387 134.0  92.0     No
    ## 3388 132.5  97.5     No
    ## 3389 133.0  93.0     No
    ## 3390 129.0  84.0     No
    ## 3391 146.0  92.0     No
    ## 3392 135.0  80.0     No
    ## 3393 115.5  85.0     No
    ## 3394 122.5  73.0     No
    ## 3395 123.0  75.0     No
    ## 3396 131.0  80.0     No
    ## 3397 119.0  81.0     No
    ## 3398 129.0  81.0     No
    ## 3399 131.0  81.0     No
    ## 3400 114.0  78.0     No
    ## 3401 110.0  70.0     No
    ## 3402 168.0 100.0     No
    ## 3403 145.0  85.0     No
    ## 3404 143.0  84.0     No
    ## 3406 108.0  75.0     No
    ## 3407 121.0  82.0     No
    ## 3408 122.0  80.0     No
    ## 3409 117.0  73.0     No
    ## 3410 124.5  84.5     No
    ## 3411 127.0  81.0     No
    ## 3412 106.0  72.0     No
    ## 3413 166.0  88.0     No
    ## 3414 105.5  74.0     No
    ## 3415 196.0 103.0     No
    ## 3416 135.0  95.0     No
    ## 3417 113.5  77.0     No
    ## 3419 102.0  70.5     No
    ## 3420 118.5  73.0     No
    ## 3421 121.0  85.0     No
    ## 3422 117.0  78.0     No
    ## 3423 115.0  70.0     No
    ## 3424 116.5  87.0     No
    ## 3425  99.0  62.0     No
    ## 3426 115.0  79.0     No
    ## 3427 122.0  75.0     No
    ## 3428 151.0  74.0     No
    ## 3429 142.0  94.0     No
    ## 3430 112.0  66.0     No
    ## 3431 117.0  78.0     No
    ## 3433 123.0  85.0     No
    ## 3434 112.0  74.5     No
    ## 3435 112.5  77.5     No
    ## 3437 164.5 102.0     No
    ## 3438 130.0  80.0     No
    ## 3439 127.0  88.5     No
    ## 3440 142.5  85.0     No
    ## 3444 120.0  80.0     No
    ## 3445 137.0  82.0     No
    ## 3446 124.5  72.0     No
    ## 3447 112.5  60.0     No
    ## 3449 128.0  86.5     No
    ## 3450 119.0  65.0     No
    ## 3451 120.0  77.5     No
    ## 3452 136.5  87.0     No
    ## 3453 101.0  59.0     No
    ## 3454 129.0  85.0     No
    ## 3455 139.0  79.0     No
    ## 3456 127.0  76.5     No
    ## 3457 119.0  75.0     No
    ## 3458 114.0  80.0     No
    ## 3459 110.0  60.0     No
    ## 3461 127.5  76.0     No
    ## 3463 119.0  82.5     No
    ## 3464 134.5  87.0     No
    ## 3466 124.0  75.5     No
    ## 3467 159.0 100.0     No
    ## 3469 108.5  73.5     No
    ## 3470 119.0  62.5     No
    ## 3471 110.0  67.0     No
    ## 3472 181.0 112.5    Yes
    ## 3473 143.5  85.0     No
    ## 3474 177.5 120.0     No
    ## 3475 119.0  80.0     No
    ## 3476 125.0  75.0     No
    ## 3477 157.0  96.0     No
    ## 3478 102.0  74.5     No
    ## 3479 111.0  78.0     No
    ## 3480 128.0  77.5     No
    ## 3481 131.0  52.0     No
    ## 3482 122.0  68.0     No
    ## 3484 124.0  80.0     No
    ## 3486 125.0  75.0     No
    ## 3487 141.0  82.0     No
    ## 3488 141.0  87.0     No
    ## 3489 111.0  67.0     No
    ## 3490 107.5  70.0     No
    ## 3491 159.0 100.0     No
    ## 3492 143.0  87.0     No
    ## 3494 160.0  87.0    Yes
    ## 3495 154.5  93.0     No
    ## 3496 133.5  92.0     No
    ## 3497 112.0  83.0     No
    ## 3498 147.5  92.5     No
    ## 3499 123.0  76.0     No
    ## 3500 135.0  88.0     No
    ## 3501 113.0  77.0     No
    ## 3502 168.0 103.0     No
    ## 3503 117.5  82.5     No
    ## 3504 138.0  80.0     No
    ## 3505 195.0 110.0     No
    ## 3506 125.0  83.0     No
    ## 3507 159.0 102.0    Yes
    ## 3508 126.0  84.0     No
    ## 3509 140.0  94.0     No
    ## 3510 134.5  80.0     No
    ## 3511 115.0  78.0     No
    ## 3512 110.0  80.0     No
    ## 3513 144.5  91.5     No
    ## 3514 149.5  86.0     No
    ## 3515 142.0  90.0     No
    ## 3516 114.0  82.0     No
    ## 3517 150.0  74.0     No
    ## 3519 128.0  86.5     No
    ## 3520 165.0  84.0     No
    ## 3521 102.5  66.0     No
    ## 3522 112.0  75.5     No
    ## 3523 127.5  75.0     No
    ## 3524 143.5  93.0     No
    ## 3525 143.0  90.0     No
    ## 3526 121.5  81.5     No
    ## 3528 120.0  72.0     No
    ## 3530 141.0  92.0     No
    ## 3531 133.0  93.0     No
    ## 3532 102.0  73.0     No
    ## 3533 116.5  81.0     No
    ## 3535 158.0 105.0     No
    ## 3536 168.0  92.0     No
    ## 3538 139.0  84.0     No
    ## 3539 125.0  72.0     No
    ## 3541 112.0  70.0     No
    ## 3542 138.0  80.0     No
    ## 3543 110.0  65.0     No
    ## 3544 113.0  78.0     No
    ## 3545 182.0  86.0     No
    ## 3546 121.0  84.0     No
    ## 3547 130.0  87.0     No
    ## 3548 140.0  94.0     No
    ## 3549 122.5  87.0     No
    ## 3550 111.0  80.0     No
    ## 3551 160.5 106.5     No
    ## 3552 103.0  67.0     No
    ## 3553 118.5  71.0     No
    ## 3554 165.0 108.0     No
    ## 3555 140.0  89.5     No
    ## 3556 132.5  55.0    Yes
    ## 3557 123.0  76.0     No
    ## 3558 162.0  97.5     No
    ## 3559 132.0  88.0     No
    ## 3560 112.0  63.5     No
    ## 3561 150.5  87.0     No
    ## 3562 157.0  91.0     No
    ## 3564 153.0 106.0     No
    ## 3565 112.0  83.0     No
    ## 3566 113.0  82.0     No
    ## 3567 192.5 110.0     No
    ## 3568 108.0  70.5     No
    ## 3569 142.0  92.0     No
    ## 3570 153.0  87.0     No
    ## 3571 130.0  80.0     No
    ## 3572 108.0  65.0     No
    ## 3573 123.0  75.0     No
    ## 3574 167.0 100.0     No
    ## 3575 228.0 130.0    Yes
    ## 3576 113.0  70.0     No
    ## 3577 125.0  77.5     No
    ## 3578 129.0  83.0     No
    ## 3579 167.0  96.0     No
    ## 3580 151.0 102.0     No
    ## 3581 127.0  79.5     No
    ## 3582 128.0  90.0     No
    ## 3583 122.0  80.0     No
    ## 3584 119.0  66.0     No
    ## 3585 145.5  92.5     No
    ## 3586 168.5  97.0     No
    ## 3587 165.0  86.0     No
    ## 3588 133.0  85.5     No
    ## 3589 163.0  91.0     No
    ## 3590 114.0  72.0     No
    ## 3591 115.0  67.0     No
    ## 3592 132.0  89.5     No
    ## 3593 135.0  86.0     No
    ## 3594 101.5  67.0     No
    ## 3595 124.0  80.0     No
    ## 3596 122.0  70.0     No
    ## 3597 112.5  60.0     No
    ## 3598 118.0  80.0     No
    ## 3599 112.0  61.0     No
    ## 3600 148.0  89.0     No
    ## 3601 147.0  94.0     No
    ## 3602 146.0  92.0     No
    ## 3603 147.5  90.0     No
    ## 3604 136.0  90.0     No
    ## 3605 129.0  86.0     No
    ## 3606 174.0  97.0    Yes
    ## 3607 128.0  87.0     No
    ## 3608 108.0  77.0     No
    ## 3609 126.0  81.0     No
    ## 3610 159.5  93.5     No
    ## 3611 138.0  96.0     No
    ## 3614 125.0  72.0     No
    ## 3615 132.0  77.0     No
    ## 3616 120.0  72.5     No
    ## 3617 156.0 105.0     No
    ## 3618 142.0  66.0     No
    ## 3619 123.0  79.0     No
    ## 3620 135.0  85.0     No
    ## 3621 108.0  81.0     No
    ## 3623 128.0  83.5     No
    ## 3624 135.0  85.0     No
    ## 3625 118.0  68.0     No
    ## 3626 138.0  72.0     No
    ## 3628 163.0  89.0     No
    ## 3629 128.5  87.5     No
    ## 3630 112.0  78.0     No
    ## 3631 144.5  95.0     No
    ## 3632 158.0 108.0     No
    ## 3633 110.0  74.0     No
    ## 3634 115.0  80.0     No
    ## 3635 127.0  79.0     No
    ## 3636 127.5  75.5     No
    ## 3637 151.5 110.0     No
    ## 3639 159.0  90.0     No
    ## 3640 158.0  74.0     No
    ## 3641 123.0  92.0     No
    ## 3642 100.5  69.0     No
    ## 3643 166.0 107.0     No
    ## 3644 161.0 105.0    Yes
    ## 3646 111.0  79.0     No
    ## 3647  83.5  55.0     No
    ## 3649 248.0 130.0    Yes
    ## 3650 126.0  73.0     No
    ## 3651 121.0  79.0     No
    ## 3652 111.0  72.5     No
    ## 3653 152.0  70.0     No
    ## 3654 155.0  90.0     No
    ## 3655 108.0  70.0     No
    ## 3656 172.5 112.5     No
    ## 3657 122.0  81.0     No
    ## 3658 163.0  94.0     No
    ## 3659 123.0  81.0     No
    ## 3660 108.0  73.0     No
    ## 3661 112.5  62.5     No
    ## 3662 107.5  71.0     No
    ## 3663 143.0  82.5     No
    ## 3664 129.0  87.0     No
    ## 3665  98.0  60.0     No
    ## 3666 116.0  76.0     No
    ## 3667 114.0  72.0     No
    ## 3668 124.0  89.0     No
    ## 3669 196.0 109.0     No
    ## 3670 160.0  98.0     No
    ## 3671 122.5  85.0     No
    ## 3672 156.0  95.0     No
    ## 3673 126.5  75.5     No
    ## 3675 108.0  75.0     No
    ## 3676 122.0  81.0     No
    ## 3677 176.5  92.0    Yes
    ## 3678 115.0  65.0     No
    ## 3679 121.0  69.0     No
    ## 3680 131.0  79.0     No
    ## 3681 127.0  77.0     No
    ## 3682 161.5  96.0     No
    ## 3683 100.5  66.0     No
    ## 3684 131.5  83.0     No
    ## 3685 131.5  77.0     No
    ## 3687 135.5  77.0     No
    ## 3688 125.0  76.0     No
    ## 3689 131.0  83.0     No
    ## 3690 143.5  77.5     No
    ## 3691 136.0  90.0     No
    ## 3692 123.5  78.0     No
    ## 3693 164.0 119.0     No
    ## 3694 145.0  87.0     No
    ## 3696 133.0  82.0     No
    ## 3697 139.0  81.0     No
    ## 3698 124.0  84.0     No
    ## 3699 139.0  81.5     No
    ## 3700 113.0  87.0     No
    ## 3701 116.0  77.0     No
    ## 3702 182.0 110.0    Yes
    ## 3703 113.0  74.0     No
    ## 3704 126.0  71.0     No
    ## 3705 112.0  67.0     No
    ## 3706 130.0  74.0     No
    ## 3707 104.0  76.0     No
    ## 3708 106.0  77.0     No
    ## 3709 177.5 110.0     No
    ## 3710 121.0  78.0     No
    ## 3711 134.0  78.0     No
    ## 3712 110.0  77.5     No
    ## 3713 128.5  80.0     No
    ## 3714 115.0  72.0     No
    ## 3715 154.0  96.0     No
    ## 3716 202.0 124.0    Yes
    ## 3717 151.0  85.0     No
    ## 3718 154.0  87.0     No
    ## 3719 162.0 109.0     No
    ## 3721 117.0  78.5     No
    ## 3723 123.0  82.0     No
    ## 3724 110.0  83.0     No
    ## 3725 117.0  72.0     No
    ## 3726 122.0  85.0     No
    ## 3727 154.0  80.0     No
    ## 3728 106.0  65.0     No
    ## 3729 114.5  77.0     No
    ## 3730  99.0  62.0     No
    ## 3731  99.5  66.0     No
    ## 3732 185.0  95.0     No
    ## 3734 108.0  70.0     No
    ## 3736 121.0  82.0    Yes
    ## 3737 116.0  67.0     No
    ## 3739 158.0  89.0     No
    ## 3740  85.0  70.0     No
    ## 3742 169.0 111.0     No
    ## 3743 116.0  66.0     No
    ## 3744 130.0  84.5     No
    ## 3745 117.0  86.0     No
    ## 3746 105.5  67.0     No
    ## 3747 135.0  88.0     No
    ## 3748 146.0  89.0     No
    ## 3749 165.0 100.0     No
    ## 3750 143.0  96.0     No
    ## 3752 176.0 109.0     No
    ## 3753 118.0  74.0     No
    ## 3756 136.0  94.0     No
    ## 3757 117.0  81.5     No
    ## 3758 130.0  94.0     No
    ## 3759 135.0  84.0     No
    ## 3760 122.5  75.0     No
    ## 3762 135.0  97.0     No
    ## 3763 107.0  69.0     No
    ## 3764 167.5  92.5     No
    ## 3765 122.5  80.0     No
    ## 3766 133.0  72.0     No
    ## 3767 138.0  87.0     No
    ## 3768 151.0  91.0     No
    ## 3770 101.0  75.0     No
    ## 3771 109.0  78.5     No
    ## 3772 132.0  78.0     No
    ## 3774 123.0  66.0     No
    ## 3775 119.0  67.0     No
    ## 3776 137.0  79.0     No
    ## 3777 130.0  89.0     No
    ## 3778 128.0  78.0     No
    ## 3780 161.0  85.0     No
    ## 3781 120.0  75.0     No
    ## 3782 230.0 110.0     No
    ## 3783 118.0  86.5     No
    ## 3785 143.0  75.0     No
    ## 3786 140.0  71.0     No
    ## 3787 155.0  79.0     No
    ## 3788 153.0 102.5     No
    ## 3789 146.5  79.0     No
    ## 3791 137.5  82.5     No
    ## 3792 138.0  72.0     No
    ## 3793 147.5  95.0     No
    ## 3794 129.0  80.0     No
    ## 3795 101.0  66.0     No
    ## 3797 162.0  85.0     No
    ## 3798 197.0  72.0     No
    ## 3799 133.0  89.0     No
    ## 3800 164.0  89.0     No
    ## 3801 189.0 121.0     No
    ## 3802 132.0  74.0     No
    ## 3803 125.5  82.0     No
    ## 3804 134.0  93.0     No
    ## 3805 118.0  92.0     No
    ## 3806 142.5  85.0     No
    ## 3807 124.0  82.0     No
    ## 3808 126.0  52.0     No
    ## 3809 136.0  95.5     No
    ## 3810 112.0  78.0     No
    ## 3812 106.0  48.0     No
    ## 3813 118.0  71.0     No
    ## 3815 147.0 100.0     No
    ## 3817 128.0  78.0     No
    ## 3818 142.0  84.0     No
    ## 3819 140.0  92.5     No
    ## 3820 113.0  68.0     No
    ## 3821 151.0  89.0     No
    ## 3823 144.0  88.0     No
    ## 3824 153.0 102.5     No
    ## 3825 146.0  86.0     No
    ## 3826 141.0  93.0     No
    ## 3829 112.5  71.0     No
    ## 3830 113.0  79.0     No
    ## 3831 100.0  76.0     No
    ## 3832 140.0  90.0     No
    ## 3833 114.0  82.0     No
    ## 3834 130.0  77.5     No
    ## 3835 155.0  92.5     No
    ## 3836 128.0  82.0     No
    ## 3838  97.0  63.0     No
    ## 3839 112.0  73.0     No
    ## 3841 167.0  92.0    Yes
    ## 3842 113.0  73.0     No
    ## 3844 166.0  90.0    Yes
    ## 3845 214.0  94.0    Yes
    ## 3846 143.0  88.0     No
    ## 3847 134.5  90.5     No
    ## 3848 132.0  91.0     No
    ## 3849 124.0  80.0     No
    ## 3850 159.5  82.5     No
    ## 3852 125.0  75.0     No
    ## 3853 123.0  73.0     No
    ## 3854 117.5  75.0     No
    ## 3856 119.0  75.0     No
    ## 3857 131.0  81.0     No
    ## 3858 123.0  74.0     No
    ## 3859 117.5  72.5     No
    ## 3860 140.0  84.0     No
    ## 3861 116.0  70.0     No
    ## 3862 114.0  75.0     No
    ## 3863 119.0  65.0     No
    ## 3864 125.0  79.0     No
    ## 3865 147.5  88.0     No
    ## 3866 124.5  86.5     No
    ## 3867 111.0  56.0     No
    ## 3868 141.0  84.0     No
    ## 3869 105.0  77.0     No
    ## 3870 127.0  83.0     No
    ## 3871 126.5  79.0     No
    ## 3872 136.5  78.0     No
    ## 3873 145.0  95.0     No
    ## 3874 134.0  86.5     No
    ## 3875 135.0  86.0     No
    ## 3876 145.0  92.0     No
    ## 3877 138.0  94.0     No
    ## 3878 149.0  96.0     No
    ## 3879 170.0 118.0     No
    ## 3880 141.0  90.0     No
    ## 3881 141.0  92.0     No
    ## 3884 142.5  82.0     No
    ## 3885 120.0  79.0     No
    ## 3886 142.0  88.0    Yes
    ## 3887 172.5  85.0     No
    ## 3888 128.0  81.0     No
    ## 3889 109.5  69.0     No
    ## 3890 104.0  78.0     No
    ## 3891 144.5  88.0     No
    ## 3892 118.0  73.0     No
    ## 3893 150.0  84.0     No
    ## 3894 158.0 102.0     No
    ## 3895 148.0  91.5     No
    ## 3896 120.0  84.0     No
    ## 3897 155.0  90.0     No
    ## 3898 108.0  70.0     No
    ## 3899 196.0 102.0     No
    ## 3900 115.0  85.0     No
    ## 3901 116.0  83.0     No
    ## 3902 167.5 110.0     No
    ## 3903 143.0  93.0     No
    ## 3904  96.5  62.5     No
    ## 3905 133.5  85.5     No
    ## 3906 132.0  81.0     No
    ## 3907 124.0  84.0     No
    ## 3909 146.5  81.0     No
    ## 3910 153.0 100.0     No
    ## 3911 122.0  78.0     No
    ## 3912 140.0  95.0     No
    ## 3913 143.0  92.0     No
    ## 3915 145.5  94.0     No
    ## 3916 163.0 112.5     No
    ## 3918 122.0  86.0     No
    ## 3919 149.0  85.0     No
    ## 3920 134.0  74.0     No
    ## 3921 110.0  67.5     No
    ## 3922 163.0 101.0     No
    ## 3924 113.0  75.5     No
    ## 3925 115.0  77.0     No
    ## 3926 100.5  66.0     No
    ## 3928 158.0  97.0     No
    ## 3929 118.0  71.0     No
    ## 3930 110.0  76.0     No
    ## 3931 122.5  84.0     No
    ## 3932 130.0  71.5     No
    ## 3933 160.0  92.0     No
    ## 3934 118.0  75.5     No
    ## 3937  99.0  59.0     No
    ## 3938 125.0  88.0     No
    ## 3939 106.0  73.0     No
    ## 3940 115.0  83.0     No
    ## 3941 126.0  82.0     No
    ## 3942 120.0  89.0     No
    ## 3944 102.0  64.5     No
    ## 3946 131.0  87.0     No
    ## 3948 113.0  83.0     No
    ## 3949 141.0  81.0     No
    ## 3950 139.0  85.0     No
    ## 3951 101.0  59.0     No
    ## 3953 100.5  62.0     No
    ## 3954 133.0  90.0     No
    ## 3955 131.0  81.0     No
    ## 3956 129.0  81.0     No
    ## 3957 115.0  90.0     No
    ## 3958 137.0  86.0     No
    ## 3960 119.0  78.0     No
    ## 3961 154.0  92.0     No
    ## 3962 136.0  87.0     No
    ## 3964 102.0  68.0     No
    ## 3965 164.0  97.0    Yes
    ## 3966 170.0 107.0    Yes
    ## 3967 152.5  90.0     No
    ## 3968 116.0  71.0     No
    ## 3969 113.0  68.0     No
    ## 3970 127.5  82.5     No
    ## 3971 141.5  91.0     No
    ## 3972 110.0  84.0     No
    ## 3973 125.0  80.0     No
    ## 3974 115.0  75.0     No
    ## 3975 137.0  87.0     No
    ## 3976 173.0  89.0     No
    ## 3977 162.5 104.0     No
    ## 3978 105.0  67.5     No
    ## 3979 100.0  70.0     No
    ## 3980 146.0  86.0     No
    ## 3981 148.5  94.0     No
    ## 3982 135.0  86.5     No
    ## 3983 128.0  82.0     No
    ## 3984 119.0  77.0     No
    ## 3985 129.0  93.0     No
    ## 3986 109.0  74.0     No
    ## 3987 115.5  82.5     No
    ## 3989 141.0  93.0     No
    ## 3990 129.0  94.0     No
    ## 3991 170.0 113.0     No
    ## 3992 121.0  81.0     No
    ## 3993 124.0  78.0     No
    ## 3994 100.0  72.0     No
    ## 3995 158.0  90.0     No
    ## 3996 109.0  61.0     No
    ## 3997  99.0  59.0     No
    ## 3998 120.5  80.0     No
    ## 3999 133.0  87.0     No
    ## 4000 130.0  71.0     No
    ## 4001 132.5  87.0     No
    ## 4002 156.5  93.0     No
    ## 4003 131.0  81.0     No
    ## 4004 108.0  70.0     No
    ## 4005 127.0  66.5     No
    ## 4006 111.0  72.0     No
    ## 4007 124.0  69.0     No
    ## 4008 167.0 101.0     No
    ## 4009 111.0  72.0     No
    ## 4010 168.5  98.5    Yes
    ## 4011 155.5  98.0     No
    ## 4012 107.0  74.0     No
    ## 4013 111.0  79.0     No
    ## 4014 123.5  80.5     No
    ## 4015 120.0  80.0     No
    ## 4016 215.0 110.0    Yes
    ## 4017 110.0  70.0     No
    ## 4018 192.5 113.0    Yes
    ## 4019 135.0  82.5     No
    ## 4020 111.0  74.0     No
    ## 4021 106.0  60.0     No
    ## 4022 149.0  91.0     No
    ## 4023 110.0  70.0     No
    ## 4024 137.0  75.0     No
    ## 4025 172.5 100.5     No
    ## 4026 130.0  62.0     No
    ## 4027 130.0  80.0     No
    ## 4028 172.0 108.0     No
    ## 4029 124.0  76.0     No
    ## 4030 102.0  56.0     No
    ## 4031 104.5  65.0     No
    ## 4032 133.0  86.5     No
    ## 4033 120.5  76.0     No
    ## 4034 151.0  90.0     No
    ## 4035 169.0  82.0    Yes
    ## 4036 111.0  80.0     No
    ## 4038 157.0  98.0     No
    ## 4039 140.0  86.5    Yes
    ## 4040 147.5  87.5     No
    ## 4041 123.0  88.0     No
    ## 4043 172.0  87.0     No
    ## 4044 148.0  89.5     No
    ## 4045 144.0  96.0     No
    ## 4048 110.0  61.0     No
    ## 4049 147.0  65.0     No
    ## 4050 160.0  85.0     No
    ## 4051 140.0  99.0     No
    ## 4053 120.0  77.0     No
    ## 4054 111.5  67.0     No
    ## 4055 120.0  80.0     No
    ## 4056 117.0  73.0     No
    ## 4058 152.5 108.0     No
    ## 4059 112.0  71.0     No
    ## 4060 106.0  73.0     No
    ## 4061 127.0  85.0     No
    ## 4062 128.0  81.0     No
    ## 4063 165.0  95.0     No
    ## 4064 105.0  70.0     No
    ## 4065 164.0 100.0     No
    ## 4066 119.0  73.0     No
    ## 4067 143.5  92.0     No
    ## 4069 121.0  81.0     No
    ## 4070 110.0  68.0     No
    ## 4071 139.0  75.0     No
    ## 4072 156.0  90.0     No
    ## 4073 110.0  65.0     No
    ## 4074 128.0  88.0     No
    ## 4075 128.0  79.0     No
    ## 4076 150.0  98.0     No
    ## 4077 116.5  74.0     No
    ## 4079 155.0  80.5     No
    ## 4080 119.0  75.5     No
    ## 4081 128.0  65.0     No
    ## 4083 173.0  88.0     No
    ## 4084 130.0  94.0     No
    ## 4085 137.0  79.0     No
    ## 4086 112.0  65.0     No
    ## 4087 187.0  97.0     No
    ## 4088 116.0  72.0     No
    ## 4089 139.0  85.0     No
    ## 4090 141.0  89.0     No
    ## 4091 127.5  81.0     No
    ## 4093 173.0  98.5     No
    ## 4094 129.0  97.0     No
    ## 4095 138.5  85.0     No
    ## 4096 130.0  86.0     No
    ## 4097  92.5  66.5     No
    ## 4098 143.0  96.0     No
    ## 4099 114.0  80.0     No
    ## 4100 110.0  76.0     No
    ## 4101 131.5  77.5     No
    ## 4102 151.5  71.0     No
    ## 4103 157.0 112.0     No
    ## 4105 131.5  84.0     No
    ## 4106 194.0 111.0    Yes
    ## 4107 170.0 100.0     No
    ## 4108 130.5  82.0     No
    ## 4109 153.0  80.0     No
    ## 4110 132.0  95.0     No
    ## 4111 150.0  99.0     No
    ## 4112 133.0  86.0     No
    ## 4113 130.5  63.0     No
    ## 4114 118.0  79.0     No
    ## 4115 144.0  80.0     No
    ## 4116 121.5  74.5     No
    ## 4117 150.0  91.0     No
    ## 4118 116.0  84.0     No
    ## 4119 126.5  81.0     No
    ## 4120  93.5  58.0     No
    ## 4121 126.5  76.5     No
    ## 4123 152.0  96.5     No
    ## 4124 132.0  88.0     No
    ## 4125 152.5  92.5     No
    ## 4127 111.0  79.0     No
    ## 4128 124.5  75.5     No
    ## 4129 111.5  72.0     No
    ## 4130 126.0  86.0     No
    ## 4131 137.0  83.0    Yes
    ## 4132 174.0 130.0     No
    ## 4133 173.0  67.0     No
    ## 4134 110.0  68.0     No
    ## 4135 150.0  82.0     No
    ## 4136 118.0  81.0     No
    ## 4138 177.5 100.0     No
    ## 4139 117.0  76.0     No
    ## 4142 161.0  90.0     No
    ## 4143 125.0  87.0     No
    ## 4144 125.0  87.0     No
    ## 4145 144.0  87.0     No
    ## 4146 102.5  63.5     No
    ## 4147 129.5 100.0     No
    ## 4148 120.5  80.0     No
    ## 4149  96.0  66.0     No
    ## 4150 112.5  76.5     No
    ## 4151 176.0  78.0     No
    ## 4152 120.0  87.0     No
    ## 4153 146.0  82.5     No
    ## 4154 147.0  90.0     No
    ## 4155 126.0  79.0     No
    ## 4156 169.0  91.0     No
    ## 4157 157.0  94.0     No
    ## 4158 110.0  80.0     No
    ## 4160 125.0  80.0     No
    ## 4161 207.0 122.5     No
    ## 4162 136.5  85.0     No
    ## 4163 145.0  99.0     No
    ## 4164 134.0  88.0     No
    ## 4165 113.0  68.0     No
    ## 4166 100.0  66.5     No
    ## 4167 129.0  77.0     No
    ## 4171 131.0  76.0     No
    ## 4172 129.0  84.0     No
    ## 4173 101.0  68.5     No
    ## 4174 122.0  78.0     No
    ## 4175 153.0  80.0     No
    ## 4176 118.0  87.0     No
    ## 4177 119.0  86.0     No
    ## 4178 132.5  82.5     No
    ## 4179 111.0  60.5     No
    ## 4180 133.0  84.0     No
    ## 4181 110.0  72.0     No
    ## 4182 122.5  75.0     No
    ## 4183 165.0  95.0     No
    ## 4184 122.5  77.5     No
    ## 4186 126.0  67.0     No
    ## 4187 142.0  87.0     No
    ## 4188 115.0  74.0     No
    ## 4189 135.0 100.0     No
    ## 4191 141.0  76.0     No
    ## 4192 118.0  73.0     No
    ## 4193 123.0  82.0     No
    ## 4195 118.0  74.0     No
    ## 4196 140.0  88.0     No
    ## 4197 108.0  72.0     No
    ## 4198 129.0  85.0     No
    ## 4200 125.0  71.0     No
    ## 4201 140.0  92.0     No
    ## 4202 107.0  70.0     No
    ## 4203 136.0  73.0     No
    ## 4204 157.0  94.0     No
    ## 4205 147.0 101.0     No
    ## 4206 133.0  90.0     No
    ## 4207 162.5  92.5     No
    ## 4208 130.0  85.0     No
    ## 4209 143.0  87.0     No
    ## 4211 173.0 106.0    Yes
    ## 4212 108.5  68.0     No
    ## 4213 130.0  86.5     No
    ## 4214 112.0  67.0     No
    ## 4215 112.5  80.0     No
    ## 4216 130.0  80.0     No
    ## 4217 159.0 110.0     No
    ## 4218 185.5 115.5     No
    ## 4219 129.0  83.0     No
    ## 4220 103.0  62.0     No
    ## 4221 126.0  75.0     No
    ## 4222 213.0 133.0     No
    ## 4223 160.0 105.0     No
    ## 4224 175.0 101.0     No
    ## 4225 111.5  74.0     No
    ## 4226 129.5  93.5     No
    ## 4227 124.0  77.0     No
    ## 4228 110.0  73.0     No
    ## 4229 118.0  79.5     No
    ## 4230 148.0  95.0     No
    ## 4232  97.5  62.5     No
    ## 4233 124.0  89.0     No
    ## 4234 104.0  64.0     No
    ## 4235 107.0  69.0     No
    ## 4236 127.0  81.0     No
    ## 4237 115.0  80.0     No
    ## 4238 130.0  80.0     No
    ## 4239 157.0  80.0    Yes
    ## 4240 192.5 115.0     No
    ## 4241 135.0  80.0     No
    ## 4242 176.0 109.0     No
    ## 4243 121.0  82.0     No
    ## 4245 112.0  72.0     No
    ## 4246 142.5  95.0     No
    ## 4247 107.0  81.0     No
    ## 4248 116.0  71.5     No
    ## 4249 117.0  78.5     No
    ## 4250 130.0  86.0     No
    ## 4251 133.5  97.5     No
    ## 4252 135.0  83.0     No
    ## 4253 158.0 103.0     No
    ## 4255 205.5 123.5     No
    ## 4257 140.0  93.0     No
    ## 4258 109.0  79.0     No
    ## 4259 118.0  79.0     No
    ## 4260 136.0  80.5     No
    ## 4261 192.5 125.0     No
    ## 4262 131.0  84.0     No
    ## 4263 200.0 125.0    Yes
    ## 4264 115.0  72.5     No
    ## 4265 120.0  81.0     No
    ## 4266 104.0  69.0     No
    ## 4267 124.0  70.0     No
    ## 4268 116.0  76.0     No
    ## 4269 128.0  81.0     No
    ## 4271 115.0  71.0     No
    ## 4272 187.0  95.0     No
    ## 4273 121.0  81.0     No
    ## 4275 126.0  74.0     No
    ## 4276 143.0  88.0     No
    ## 4277 149.5  84.0     No
    ## 4278 151.0  87.5     No
    ## 4279 116.0  74.0     No
    ## 4280  93.0  62.5     No
    ## 4281 125.0  80.0     No
    ## 4282 121.0  79.0     No
    ## 4283 140.5  83.0     No
    ## 4284 136.0  80.0     No
    ## 4285 152.0  82.0     No
    ## 4287 122.0  67.5     No
    ## 4288 123.0  77.0     No
    ## 4289 133.0  97.0     No
    ## 4290 134.0  75.0     No
    ## 4291 130.0  84.0     No
    ## 4292 123.0  80.0     No
    ## 4293 142.0  86.0     No
    ## 4294 110.0  74.5     No
    ## 4295 176.0  98.0    Yes
    ## 4296 104.0  79.0     No
    ## 4297 158.0 102.5    Yes
    ## 4298 146.5  92.0     No
    ## 4299 139.0  74.0     No
    ## 4300 141.0  80.0     No
    ## 4301 108.0  68.0     No
    ## 4302 142.0  68.0     No
    ## 4303 115.0  73.0     No
    ## 4304 190.0  88.0     No
    ## 4305 149.0  83.0     No
    ## 4306 146.5  97.5     No
    ## 4307 137.0  82.0     No
    ## 4308 122.0  72.0     No
    ## 4309 140.0  77.0    Yes
    ## 4310 127.0  81.0     No
    ## 4313 206.0 116.0     No
    ## 4314 108.5  70.5     No
    ## 4315 164.0 102.0     No
    ## 4316 124.0  84.0     No
    ## 4317 155.0  80.0     No
    ## 4318 114.0  80.0     No
    ## 4319 147.5  87.5     No
    ## 4320 125.0  82.0     No
    ## 4321 120.0  77.5     No
    ## 4322 170.0 101.0     No
    ## 4323  97.5  60.0     No
    ## 4324 134.0  96.0     No
    ## 4325 115.0  78.0     No
    ## 4326 128.0  86.0     No
    ## 4327 118.0  80.0     No
    ## 4330 126.0  73.0     No
    ## 4331 135.0  82.0     No
    ## 4332 118.0  81.0     No
    ## 4333 162.0 110.0     No
    ## 4334 105.0  82.0     No
    ## 4335 131.5  84.0     No
    ## 4336 144.0  96.5     No
    ## 4337 103.0  70.0     No
    ## 4338 128.0  74.0     No
    ## 4339 108.0  78.0     No
    ## 4340 136.0  93.0     No
    ## 4341 127.5  75.0     No
    ## 4342 103.0  72.5     No
    ## 4343 147.0  85.0     No
    ## 4345 157.5  83.0     No
    ## 4346 126.0  84.0     No
    ## 4347  99.0  62.0     No
    ## 4348 125.5  85.5     No
    ## 4349 143.0  81.0    Yes
    ## 4350 135.0  97.5     No
    ## 4352 109.0  72.0     No
    ## 4354 118.5  77.5     No
    ## 4355 136.0  86.0     No
    ## 4357 130.0  81.5     No
    ## 4358 133.0  78.0     No
    ## 4359 106.0  58.0     No
    ## 4360 180.0 108.0    Yes
    ## 4361 101.0  70.0     No
    ## 4363 125.0  87.0     No
    ## 4364 154.5  88.0     No
    ## 4365 210.0 127.5     No
    ## 4366 129.5  80.0     No
    ## 4367 135.0  94.0     No
    ## 4368 107.0  75.0     No
    ## 4369 111.0  70.0     No
    ## 4370 121.5  57.0     No
    ## 4371 132.5  85.0     No
    ## 4372 137.5  89.0     No
    ## 4373 129.0  95.0     No
    ## 4374 121.0  83.0     No
    ## 4375 157.0  86.0     No
    ## 4376 182.0  99.0     No
    ## 4377 118.5  70.0     No
    ## 4378 137.5  87.5     No
    ## 4380 106.0  63.0     No
    ## 4381 143.0  79.0     No
    ## 4382 110.0  70.0     No
    ## 4383 117.0  74.0     No
    ## 4384 150.0  89.0     No
    ## 4385 101.0  69.0     No
    ## 4386 129.0  80.0     No
    ## 4387 195.0 105.0     No
    ## 4388 179.0  96.0    Yes
    ## 4389 124.0  78.0     No
    ## 4390 116.0  80.0     No
    ## 4391 170.0 105.0     No
    ## 4392 157.5 105.0     No
    ## 4393 133.0  83.0     No
    ## 4394 115.0  60.0     No
    ## 4395 130.0  80.0     No
    ## 4396 160.5 100.0     No
    ## 4397 146.0  84.0     No
    ## 4398 142.0  72.0     No
    ## 4399 136.0  77.0     No
    ## 4400 103.0  67.5     No
    ## 4401 124.0  76.5     No
    ## 4403 135.0  80.0     No
    ## 4404 126.5  88.0     No
    ## 4405 105.0  70.0     No
    ## 4406 144.0  88.0     No
    ## 4407 141.0  95.0     No
    ## 4408 123.0  78.5     No
    ## 4409 155.0  82.0     No
    ## 4410 125.0  80.0     No
    ## 4411 167.0  94.0     No
    ## 4412 137.5  84.5     No
    ## 4413 125.0  84.5     No
    ## 4414 128.0  82.0     No
    ## 4415 119.0  74.0     No
    ## 4416 188.0 110.0     No
    ## 4417 149.0  98.0    Yes
    ## 4418 120.0  80.0     No
    ## 4419 137.5  85.0     No
    ## 4420 122.0  84.0     No
    ## 4421 125.5  84.0     No
    ## 4422 129.5  88.0     No
    ## 4423 190.0 130.0     No
    ## 4426 141.0  81.0     No
    ## 4427 168.0  97.0     No
    ## 4428 179.0  92.0     No
    ## 4429 126.5  80.0     No
    ## 4432 133.5  83.0     No
    ## 4433 141.0  98.0     No
    ## 4434 133.0  86.0     No

This can come in very useful when subsetting data frames. Other useful
functions include `starts_with()` and `ends_with()`, for column name
prefixes and suffixes respectively.

-----

In a similar manner to subsetting columns, `filter()` can be used to
filter observations or rows.

``` r
fhsBMI40 <- filter(fhs, BMI >= 40)
fhsBMI40
```

    ##       SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1    Male     178  52 160.0  98.0       No 40.11      Yes     No     225
    ## 2  Female     183  45 151.0 101.0       No 45.80       No     No      63
    ## 3  Female     278  66 187.0  88.0       No 40.52       No     No      84
    ## 4  Female     266  62 173.0  89.0       No 42.00       No     No      75
    ## 5  Female     180  60 200.0 122.5      Yes 44.27      Yes     No     150
    ## 6  Female     251  67 192.0 102.0       No 44.09       No    Yes      62
    ## 7  Female     248  53 200.0 140.0       No 43.30      Yes     No     130
    ## 8  Female     199  42 141.0  92.0      Yes 43.69       No     No      60
    ## 9  Female     256  61 160.0 109.0       No 42.53       No     No      79
    ## 10 Female     169  44 179.0 107.0       No 44.55       No     No      77
    ## 11 Female     263  44 141.5  93.0       No 42.98       No     No      74
    ## 12 Female     212  66 220.0  96.0       No 44.71       No     No      95
    ## 13 Female     274  37 197.5 125.0      Yes 43.48       No     No      94
    ## 14 Female     208  55 190.0 130.0       No 56.80       No     No      86
    ## 15 Female     262  56 150.0  89.0       No 40.21       No     No      86
    ## 16 Female     179  63 149.0  89.0       No 40.81       No     No      77
    ## 17 Female     214  49 172.0 111.0      Yes 40.51       No     No      70
    ## 18 Female     257  37 141.0  93.0      Yes 41.29       No     No      58
    ## 19   Male     263  48 132.0  91.0      Yes 40.08       No     No      91
    ## 20 Female     254  57 146.5  81.0      Yes 41.61       No     No      85
    ## 21   Male     240  40 150.0  98.0       No 40.38       No     No      74
    ## 22 Female     225  61 194.0 111.0       No 51.28       No    Yes     103
    ## 23 Female     385  58 165.0  95.0       No 41.66       No     No      91
    ## 24 Female     252  60 182.0  99.0       No 40.23       No     No      60
    ## 25 Female     260  50 190.0 130.0       No 43.67      Yes     No     260
    ##    EDUC  MI
    ## 1     1  No
    ## 2     1  No
    ## 3     1  No
    ## 4     1 Yes
    ## 5     1  No
    ## 6     1  No
    ## 7     1  No
    ## 8     1  No
    ## 9     1  No
    ## 10    2  No
    ## 11    1  No
    ## 12    1  No
    ## 13    1 Yes
    ## 14    1  No
    ## 15    1  No
    ## 16    1  No
    ## 17    1  No
    ## 18    2  No
    ## 19    1  No
    ## 20    1 Yes
    ## 21    4  No
    ## 22    1 Yes
    ## 23    1  No
    ## 24    1  No
    ## 25    1 Yes

Here, we can quickly see the 25 individuals who are obese in our data
set. A `,` can be used as a logical AND between filters. Perhaps we are
interested in girls who are current smokers.

``` r
fhsSmokeF <- filter(fhs, SEX == 'Female', CURSMOKE == 'Yes')
fhsSmokeF
```

    ##        SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1   Female     225  61 150.0  95.0      Yes 28.58       No     No     103
    ## 2   Female     285  46 130.0  84.0      Yes 23.10       No     No      85
    ## 3   Female     313  45 100.0  71.0      Yes 21.68       No     No      78
    ## 4   Female     221  38 140.0  90.0      Yes 21.35       No     No      70
    ## 5   Female     291  46 112.0  78.0      Yes 23.38       No     No      89
    ## 6   Female     195  38 122.0  84.5      Yes 23.24       No     No      78
    ## 7   Female     190  42 108.0  70.5      Yes 21.59       No     No      85
    ## 8   Female     215  52 132.0  82.0      Yes 25.11       No     No      75
    ## 9   Female     233  42 153.0 101.0      Yes 28.93       No     No      90
    ## 10  Female     243  43 116.5  80.0      Yes 26.87       No     No      78
    ## 11  Female     237  41 122.0  78.0      Yes 23.28       No     No      74
    ## 12  Female     179  63 116.0  69.0      Yes 22.15       No     No      75
    ## 13  Female     267  63 156.5  92.5      Yes 27.10       No     No      79
    ## 14  Female     237  47 130.0  78.0      Yes 19.66       No     No      75
    ## 15  Female     250  46 116.0  71.0      Yes 20.35       No     No      94
    ## 16  Female     266  54 114.0  76.0      Yes 17.61      Yes     No      55
    ## 17  Female     205  40 158.0 102.0      Yes 25.45       No     No      87
    ## 18  Female     235  57 126.5  80.0      Yes 24.88       No     No      72
    ## 19  Female     300  47 112.5  60.0      Yes 20.13       No     No      83
    ## 20  Female     189  41 150.0 106.0      Yes 33.80       No     No      75
    ## 21  Female     221  44 110.0  76.0      Yes 22.16       No     No      83
    ## 22  Female     170  52 124.0  78.0      Yes 26.03       No     No      82
    ## 23  Female     197  36 113.0  72.5      Yes 22.73       No     No      65
    ## 24  Female     326  61 200.0 104.0      Yes 38.46       No     No      78
    ## 25  Female     239  63 134.0  80.0      Yes 26.64      Yes     No     126
    ## 26  Female     269  56 121.0  75.0      Yes 22.36       No     No      66
    ## 27  Female     220  47 132.5  87.0      Yes 27.98       No     No      75
    ## 28  Female     268  45 110.0  64.0      Yes 20.68       No     No      71
    ## 29  Female     173  42 100.0  63.0      Yes 23.25       No     No      99
    ## 30  Female     185  37 100.0  68.0      Yes 18.38       No     No      72
    ## 31  Female     237  46 112.0  70.0      Yes 20.20       No     No      62
    ## 32  Female     246  56 128.0  64.0      Yes 25.54       No     No      92
    ## 33  Female     150  50 121.0  84.0      Yes 28.69       No     No      88
    ## 34  Female     187  41 154.0 100.0      Yes 20.50       No     No      78
    ## 35  Female     279  56 136.0  94.0      Yes 32.99       No     No     102
    ## 36  Female     259  59 141.0  86.0      Yes 25.97       No     No      86
    ## 37  Female     159  36 121.5  73.0      Yes 20.41       No     No      75
    ## 38  Female     205  37 110.0  78.0      Yes 24.43       No     No      75
    ## 39  Female     197  46 144.0  78.0      Yes 22.51       No     No      60
    ## 40  Female     265  50 110.0  74.0      Yes 25.26       No     No      88
    ## 41  Female     200  39 111.0  64.0      Yes 19.24       No     No      60
    ## 42  Female     308  55 124.0  76.0      Yes 27.23       No     No      68
    ## 43  Female     229  61 122.0  83.0      Yes 25.45       No     No      61
    ## 44  Female     236  39 117.5  71.0      Yes 27.27       No     No      74
    ## 45  Female     214  42 110.0  67.0      Yes 22.54       No     No      75
    ## 46  Female     225  53 128.0  77.0      Yes 23.95       No     No      78
    ## 47  Female     224  45 117.0  74.5      Yes 16.75       No     No      87
    ## 48  Female     226  49 106.0  71.0      Yes 22.89       No     No      57
    ## 49  Female     248  55 157.0  82.5      Yes 22.91       No     No      83
    ## 50  Female     215  58 170.0  86.0      Yes 29.06       No     No      98
    ## 51  Female     186  37 135.0  91.0      Yes 21.48       No     No      84
    ## 52  Female     212  56 117.5  72.5      Yes 27.30       No     No      75
    ## 53  Female     150  36 117.5  77.5      Yes 23.71       No     No      74
    ## 54  Female     180  38 124.0  66.0      Yes 29.29       No     No      68
    ## 55  Female     135  36 108.0  74.0      Yes 22.53       No     No      75
    ## 56  Female     273  55 122.0  84.0      Yes 27.15       No     No      97
    ## 57  Female     310  62 178.0 127.0      Yes 31.77       No     No      79
    ## 58  Female     197  40 124.0  76.0      Yes 18.06       No     No      69
    ## 59  Female     180  60 200.0 122.5      Yes 44.27      Yes     No     150
    ## 60  Female     254  57 174.0  84.5      Yes 24.22       No     No      76
    ## 61  Female     262  53 127.5  86.0      Yes 24.11       No     No      73
    ## 62  Female     228  36 111.0  68.0      Yes 23.86       No     No      68
    ## 63  Female     242  54 125.0  76.0      Yes 22.16       No     No      87
    ## 64  Female     206  34 101.0  63.0      Yes 21.50       No     No      66
    ## 65  Female     188  46  97.0  65.0      Yes 21.17       No     No      60
    ## 66  Female     269  39  97.0  64.0      Yes 23.09       No     No      67
    ## 67  Female     260  63 168.0  98.0      Yes 21.05       No     No      73
    ## 68  Female     240  53 131.0  82.0      Yes 24.22       No     No      80
    ## 69  Female     220  51 112.0  74.0      Yes 31.23       No     No      66
    ## 70  Female     186  43 120.0  72.0      Yes 24.33       No     No      86
    ## 71  Female     197  48 101.0  67.0      Yes 21.35       No     No     100
    ## 72  Female     218  42 126.0  87.0      Yes 22.50       No     No      73
    ## 73  Female     237  42 105.0  75.0      Yes 23.85       No     No      87
    ## 74  Female     207  41 111.0  60.0      Yes 18.48       No     No      76
    ## 75  Female     238  63 136.0  66.0      Yes 20.20       No     No      92
    ## 76  Female     276  53 130.0  86.0      Yes 27.09       No     No      56
    ## 77  Female     231  54 127.5  83.0      Yes 21.31      Yes     No     115
    ## 78  Female     210  40 118.0  79.0      Yes 21.21       No     No      84
    ## 79  Female     203  38 100.0  70.0      Yes 22.73       No     No      80
    ## 80  Female     244  46  98.0  57.0      Yes 24.01       No     No      95
    ## 81  Female     200  42  95.0  55.0      Yes 23.68       No     No      83
    ## 82  Female     290  58 155.0  82.5      Yes 29.50       No     No     113
    ## 83  Female     285  45 116.0  87.0      Yes 23.85       No     No      55
    ## 84  Female     185  39 111.0  67.0      Yes 23.87       No     No      87
    ## 85  Female     155  47 122.5  77.5      Yes 21.34       No     No      78
    ## 86  Female     323  55 197.0 118.0      Yes 27.51       No     No     112
    ## 87  Female     206  61 130.0  80.0      Yes 21.93       No     No      82
    ## 88  Female     235  52 119.0  82.0      Yes 24.25       No     No      79
    ## 89  Female     177  42 112.5  70.0      Yes 20.62       No     No      83
    ## 90  Female     176  38 110.0  80.0      Yes 24.03       No     No     113
    ## 91  Female     268  48 117.5  80.0      Yes 36.11       No     No      67
    ## 92  Female     199  33 116.0  81.0      Yes 21.61       No     No      93
    ## 93  Female     265  43 107.0  68.0      Yes 21.08       No     No      95
    ## 94  Female     266  48 115.0  75.0      Yes 31.16       No     No      90
    ## 95  Female     292  50 132.0  82.0      Yes 22.54      Yes     No     110
    ## 96  Female     186  47 150.0  85.0      Yes 22.53       No    Yes      93
    ## 97  Female     195  48 109.0  71.0      Yes 21.10       No     No      65
    ## 98  Female     284  40 124.0  83.0      Yes 27.90       No     No      71
    ## 99  Female     177  35 110.0  70.0      Yes 25.71       No     No      84
    ## 100 Female     250  57 125.0  74.0      Yes 21.08       No     No      72
    ## 101 Female     152  44 110.0  64.0      Yes 25.71       No     No      83
    ## 102 Female     247  46 125.0  80.0      Yes 21.51       No     No      80
    ## 103 Female     250  37 138.0  90.0      Yes 19.56       No     No      74
    ## 104 Female     320  51 145.0  85.0      Yes 24.03       No     No      98
    ## 105 Female     181  44 150.0 101.0      Yes 23.74       No     No      86
    ## 106 Female     269  46 134.0  78.0      Yes 26.80       No     No     104
    ## 107 Female     210  40 103.0  71.0      Yes 24.40       No     No      68
    ## 108 Female     201  41 108.0  71.0      Yes 20.47       No     No      75
    ## 109 Female     248  51 139.0  81.0      Yes 31.16       No     No      95
    ## 110 Female     218  46 115.5  62.0      Yes 23.48       No     No      77
    ## 111 Female     210  53 132.0  84.5      Yes 27.08       No     No      84
    ## 112 Female     241  46 130.0  82.0      Yes 34.84       No     No      93
    ## 113 Female     176  41 113.0  75.0      Yes 22.29       No     No      55
    ## 114 Female     247  40 125.0  83.0      Yes 22.55       No     No      80
    ## 115 Female     233  55 128.0  94.0      Yes 36.62       No     No      95
    ## 116 Female     213  50 140.0  82.0      Yes 22.18       No     No      72
    ## 117 Female     250  44 136.5  83.5      Yes 21.33       No     No      95
    ## 118 Female     240  37 120.0  79.0      Yes 23.09       No     No      80
    ## 119 Female     242  44 135.0  89.0      Yes 23.29       No     No      77
    ## 120 Female     247  49 121.0  82.0      Yes 29.07       No     No      69
    ## 121 Female     277  39 148.0 100.0      Yes 24.12       No     No      72
    ## 122 Female     176  51 146.0  94.0      Yes 27.42       No     No      85
    ## 123 Female     240  46 125.0  74.0      Yes 22.89       No     No      76
    ## 124 Female     210  52 185.0 114.0      Yes 27.01       No     No      83
    ## 125 Female     294  60 220.0 118.0      Yes 24.22       No    Yes      59
    ## 126 Female     220  46 123.0  88.0      Yes 32.49       No     No      79
    ## 127 Female     205  40 125.0  73.5      Yes 20.68       No     No      99
    ## 128 Female     216  42 120.0  70.0      Yes 21.93       No     No      88
    ## 129 Female     246  47 113.0  75.0      Yes 21.66       No     No      68
    ## 130 Female     245  61 140.0  73.0      Yes 30.74       No     No      91
    ## 131 Female     215  53 153.0 108.0      Yes 27.31       No     No      70
    ## 132 Female     250  64 145.0  79.0      Yes 25.16       No     No      86
    ## 133 Female     244  52 127.5  72.5      Yes 24.29       No     No     118
    ## 134 Female     398  51 161.0  96.0      Yes 23.63       No     No      83
    ## 135 Female     264  51 139.5  89.0      Yes 29.38       No     No      76
    ## 136 Female     272  47 127.5  87.5      Yes 22.35       No     No      72
    ## 137 Female     187  42  96.0  67.0      Yes 24.23       No     No      84
    ## 138 Female     190  52 117.0  75.0      Yes 21.48       No     No      67
    ## 139 Female     195  41 148.0 108.0      Yes 18.21       No     No      69
    ## 140 Female     266  56 134.5  78.5      Yes 30.78       No     No      84
    ## 141 Female     229  38 117.5  67.5      Yes 23.47       No     No      80
    ## 142 Female     273  42 111.0  73.0      Yes 19.27       No     No      89
    ## 143 Female     314  60 141.0  93.0      Yes 25.23       No     No      94
    ## 144 Female     160  39 128.5  74.5      Yes 20.56       No     No      83
    ## 145 Female     281  58 134.0  81.0      Yes 22.54       No     No      74
    ## 146 Female     292  56 111.0  70.0      Yes 23.17       No     No      74
    ## 147 Female     203  50 128.5  82.0      Yes 18.99       No     No      84
    ## 148 Female     240  41 107.0  68.5      Yes 23.47       No     No      83
    ## 149 Female     206  59 167.0  89.5      Yes 25.83       No     No      75
    ## 150 Female     289  56 150.0  92.0      Yes 25.68       No    Yes      84
    ## 151 Female     199  38 112.0  68.5      Yes 23.88       No     No      67
    ## 152 Female     248  49 137.0  79.0      Yes 21.60       No     No      74
    ## 153 Female     156  45 119.0  83.0      Yes 22.02       No     No      78
    ## 154 Female     272  54 132.5  91.0      Yes 23.09       No     No      78
    ## 155 Female     203  46 117.0  68.0      Yes 21.50       No     No      85
    ## 156 Female     212  60 186.0 102.0      Yes 23.06       No     No      60
    ## 157 Female     262  54 136.0  86.0      Yes 23.28       No     No      69
    ## 158 Female     245  35 148.0  84.0      Yes 23.74       No     No      73
    ## 159 Female     200  38 124.0  84.0      Yes 20.67       No     No      75
    ## 160 Female     240  39 120.0  80.0      Yes 24.79       No     No      75
    ## 161 Female     168  37 117.0  74.0      Yes 21.51       No     No      77
    ## 162 Female     296  48 117.0  73.0      Yes 24.59       No     No      78
    ## 163 Female     290  61 178.0  91.0      Yes 28.87       No     No      80
    ## 164 Female     249  48 132.0  78.0      Yes 23.10       No     No     137
    ## 165 Female     210  45 127.0  76.5      Yes 21.67       No     No      72
    ## 166 Female     193  41 134.5  83.0      Yes 22.28       No     No     127
    ## 167 Female     278  53 131.0  87.0      Yes 33.38       No     No      74
    ## 168 Female     192  42  96.5  71.0      Yes 26.03       No     No      68
    ## 169 Female     253  46 118.0  82.0      Yes 19.70       No     No      70
    ## 170 Female     199  43 137.0  81.0      Yes 21.85       No     No      72
    ## 171 Female     230  64 177.0 110.0      Yes 28.91       No    Yes     113
    ## 172 Female     167  59 156.0 104.0      Yes 15.96       No     No      45
    ## 173 Female     242  53 127.0  79.0      Yes 19.64       No     No      74
    ## 174 Female     231  35 150.0  90.0      Yes 23.09       No     No      72
    ## 175 Female     191  47 125.0  72.5      Yes 23.81       No     No      85
    ## 176 Female     239  60 164.0  94.5      Yes 25.01       No     No      89
    ## 177 Female     178  36 102.5  65.0      Yes 20.87       No     No      94
    ## 178 Female     217  37 110.0  78.5      Yes 32.26       No     No      84
    ## 179 Female     328  47 134.0  87.0      Yes 22.34       No     No      99
    ## 180 Female     234  47 128.0  91.0      Yes 25.59       No     No      93
    ## 181 Female     240  59 122.5  67.5      Yes 25.40       No     No      81
    ## 182 Female     283  64 163.0  85.0      Yes 21.17       No     No      68
    ## 183 Female     261  53 136.0  99.0      Yes 21.02       No     No      94
    ## 184 Female     234  58 113.0  77.0      Yes 20.68       No     No      67
    ## 185 Female     202  40 104.0  76.0      Yes 19.93       No     No      62
    ## 186 Female     221  50 112.0  69.0      Yes 24.07       No     No      79
    ## 187 Female     233  39 126.0  85.0      Yes 22.89       No     No      87
    ## 188 Female     199  42 141.0  92.0      Yes 43.69       No     No      60
    ## 189 Female     200  57 108.0  77.0      Yes 18.55       No     No      87
    ## 190 Female     244  62 168.0 102.0      Yes 26.39       No     No     105
    ## 191 Female     217  44 124.5  82.0      Yes 22.36       No     No      68
    ## 192 Female     201  43 129.0  92.0      Yes 24.54       No     No      63
    ## 193 Female     154  35 125.0  75.0      Yes 23.10       No     No      75
    ## 194 Female     253  46 118.0  74.0      Yes 26.42       No     No      64
    ## 195 Female     194  42 111.0  67.5      Yes 21.34       No     No      47
    ## 196 Female     353  60 116.0  82.0      Yes 22.66       No     No      71
    ## 197 Female     196  45 123.0  71.0      Yes 20.56       No     No      76
    ## 198 Female     247  53 139.0  88.0      Yes 23.71       No     No      53
    ## 199 Female     206  62 242.0 141.0      Yes 39.86       No    Yes      94
    ## 200 Female     229  37 111.0  70.0      Yes 20.24       No     No      70
    ## 201 Female     232  47 133.0  86.0      Yes 20.15       No     No      74
    ## 202 Female     247  47 160.0  85.0      Yes 27.05       No     No      77
    ## 203 Female     215  47 128.0  85.0      Yes 20.89       No     No      90
    ## 204 Female     220  34 117.5  67.5      Yes 20.79       No     No      86
    ## 205 Female     170  39 137.5  77.5      Yes 27.35       No     No      70
    ## 206 Female     160  36  98.0  66.0      Yes 25.07       No     No      73
    ## 207 Female     218  42 109.5  67.0      Yes 23.48       No     No      71
    ## 208 Female     284  39 115.5  65.5      Yes 20.39       No     No      78
    ## 209 Female     205  37 120.5  67.5      Yes 22.89       No     No     113
    ## 210 Female     258  45 111.0  72.0      Yes 26.24       No     No      65
    ## 211 Female     192  48 135.0  82.5      Yes 32.67       No     No      69
    ## 212 Female     300  53 127.0  89.0      Yes 25.46       No     No      70
    ## 213 Female     177  50 121.0  67.0      Yes 22.02       No     No      77
    ## 214 Female     310  61 118.0  77.5      Yes 24.03       No     No      70
    ## 215 Female     220  51 142.0  82.5      Yes 21.02       No     No      78
    ## 216 Female     187  54 133.0  88.0      Yes 31.82       No     No      77
    ## 217 Female     264  62 142.0  90.0      Yes 31.78       No     No      97
    ## 218 Female     285  56 165.0 115.0      Yes 24.25       No    Yes     116
    ## 219 Female     280  49 120.0  80.0      Yes 22.33       No     No      75
    ## 220 Female     236  59 123.0  76.5      Yes 30.67       No     No     107
    ## 221 Female     252  55 108.5  63.5      Yes 25.23       No     No     121
    ## 222 Female     140  41 110.0  60.0      Yes 23.38       No     No      82
    ## 223 Female     170  45 109.5  69.0      Yes 17.38       No     No      66
    ## 224 Female     225  54 131.0  79.0      Yes 25.91       No     No      62
    ## 225 Female     265  49 144.0  86.0      Yes 25.57       No     No      68
    ## 226 Female     346  55 131.0  81.0      Yes 22.69       No     No      77
    ## 227 Female     216  40 112.5  76.5      Yes 27.22       No     No      77
    ## 228 Female     239  40 118.0  78.0      Yes 23.48       No     No      75
    ## 229 Female     246  52 113.5  66.5      Yes 19.47       No     No      60
    ## 230 Female     279  37 110.0  72.5      Yes 24.89       No     No      70
    ## 231 Female     352  60 149.0  73.0      Yes 25.96       No     No      79
    ## 232 Female     209  42 105.0  65.0      Yes 23.80       No     No      64
    ## 233 Female     175  39 105.5  64.5      Yes 25.83       No     No      72
    ## 234 Female     253  42 109.0  74.0      Yes 24.38       No     No      60
    ## 235 Female     236  63 189.0 103.0      Yes 27.91       No     No      74
    ## 236 Female     245  44 125.0  80.5      Yes 24.58       No     No      80
    ## 237 Female     192  45 132.0  79.0      Yes 24.53       No     No     112
    ## 238 Female     392  46 113.0  68.0      Yes 23.35       No     No      63
    ## 239 Female     242  67 172.0  84.0      Yes 19.81       No     No     111
    ## 240 Female     242  50 116.0  69.0      Yes 21.65       No     No      73
    ## 241 Female     291  46 107.5  65.0      Yes 24.10       No     No      78
    ## 242 Female     223  41 119.0  73.0      Yes 24.22       No     No      77
    ## 243 Female     334  52 147.0  86.0      Yes 29.01      Yes     No      63
    ## 244 Female     315  64 135.0  80.0      Yes 25.23       No     No      89
    ## 245 Female     242  41 139.0  80.0      Yes 19.68       No     No      60
    ## 246 Female     216  66 133.0  87.0      Yes 30.06       No     No      91
    ## 247 Female     214  36 119.0  76.0      Yes 21.67       No     No      75
    ## 248 Female     199  40 122.0  82.0      Yes 22.16       No     No      77
    ## 249 Female     304  42 119.0  76.0      Yes 32.52       No     No      80
    ## 250 Female     234  52 111.0  81.0      Yes 22.35       No     No      77
    ## 251 Female     155  38  96.0  61.0      Yes 24.19       No     No      68
    ## 252 Female     149  39 122.0  72.0      Yes 21.30       No     No      75
    ## 253 Female     160  38 102.5  67.5      Yes 21.16       No     No      68
    ## 254 Female     226  41 130.0  80.0      Yes 25.25       No     No      73
    ## 255 Female     197  42 124.0  81.0      Yes 21.50       No     No      63
    ## 256 Female     231  47 133.5  76.0      Yes 25.77       No     No      73
    ## 257 Female     261  45 140.0  88.0      Yes 21.44       No     No      78
    ## 258 Female     298  50 156.0  90.0      Yes 24.24       No     No     100
    ## 259 Female     213  43 113.0  77.0      Yes 29.34       No     No      73
    ## 260 Female     266  49 159.0  88.0      Yes 20.66       No     No      84
    ## 261 Female     285  47 122.0  70.0      Yes 23.48       No     No      82
    ## 262 Female     248  46 128.0  76.0      Yes 28.87       No    Yes      77
    ## 263 Female     184  59 122.0  74.0      Yes 24.66       No     No      67
    ## 264 Female     270  59 175.0  95.0      Yes 29.69       No     No      76
    ## 265 Female     244  40 110.0  73.0      Yes 21.84       No     No      67
    ## 266 Female     222  55 103.0  61.0      Yes 23.18       No     No      75
    ## 267 Female     222  45  95.0  58.0      Yes 21.68       No     No      77
    ## 268 Female     186  35 106.0  78.0      Yes 24.73       No     No      70
    ## 269 Female     235  53 154.0  98.0      Yes 26.91       No     No      65
    ## 270 Female     152  42 122.0  76.0      Yes 21.26       No     No      78
    ## 271 Female     304  48 102.0  66.5      Yes 28.90      Yes     No      66
    ## 272 Female     200  44 111.0  79.0      Yes 27.29       No     No      74
    ## 273 Female     272  57 112.5  70.0      Yes 23.08       No     No      58
    ## 274 Female     260  49 123.0  80.0      Yes 23.10       No     No      65
    ## 275 Female     190  39 137.0  81.0      Yes 19.57       No     No      85
    ## 276 Female     235  41 144.0  88.0      Yes 24.16       No     No      82
    ## 277 Female     317  64 182.5  88.0      Yes 20.52       No     No      79
    ## 278 Female     292  43 109.0  73.0      Yes 22.87       No     No      93
    ## 279 Female     175  42 132.0  86.0      Yes 20.53       No     No      88
    ## 280 Female     251  46 121.0  81.0      Yes 23.05       No     No      84
    ## 281 Female     239  67 154.0  90.0      Yes 28.56       No     No      90
    ## 282 Female     217  45 109.0  72.0      Yes 33.65       No     No      68
    ## 283 Female     246  48 129.0  86.0      Yes 25.04       No     No      87
    ## 284 Female     161  46 100.0  64.0      Yes 20.66       No     No      60
    ## 285 Female     255  39 142.0  85.5      Yes 24.89       No     No     108
    ## 286 Female     345  51 142.0  88.0      Yes 19.05       No     No      73
    ## 287 Female     187  38 118.0  78.0      Yes 30.06       No     No      63
    ## 288 Female     190  41  95.0  57.0      Yes 20.00       No     No      77
    ## 289 Female     174  43 158.5 100.5      Yes 35.99       No     No      88
    ## 290 Female     160  38  95.0  65.0      Yes 21.99       No     No      77
    ## 291 Female     250  46 115.0  74.0      Yes 22.70       No     No      69
    ## 292 Female     221  53 131.0  89.0      Yes 24.09       No     No      95
    ## 293 Female     186  49 120.0  74.0      Yes 19.39       No     No      69
    ## 294 Female     232  47 110.0  70.0      Yes 25.86       No     No      82
    ## 295 Female     325  47 160.0  95.0      Yes 32.07       No     No      87
    ## 296 Female     214  46 128.0  71.0      Yes 21.82       No     No      66
    ## 297 Female     326  61 141.0  75.0      Yes 26.11       No     No      72
    ## 298 Female     242  58 123.0  69.0      Yes 23.38       No     No      72
    ## 299 Female     304  40 121.0  88.0      Yes 22.52       No     No      80
    ## 300 Female     232  43 122.0  70.0      Yes 23.09       No     No      77
    ## 301 Female     241  47 122.5  77.0      Yes 22.18       No     No      78
    ## 302 Female     252  47 132.5  85.0      Yes 20.05       No     No      80
    ## 303 Female     253  47 129.0  81.0      Yes 22.18       No     No     122
    ## 304 Female     306  41 199.0 106.0      Yes 38.75       No     No      75
    ## 305 Female     220  36 125.0  85.0      Yes 21.34       No     No      82
    ## 306 Female     326  51 101.0  67.0      Yes 22.73       No     No      87
    ## 307 Female     195  43 104.0  57.0      Yes 20.86       No     No      78
    ## 308 Female     185  44 133.0  69.0      Yes 22.34       No     No      76
    ## 309 Female     183  53 129.0  80.0      Yes 26.51       No     No      80
    ## 310 Female     213  45 150.0  90.0      Yes 22.35       No    Yes      72
    ## 311 Female     182  46 117.0  78.0      Yes 22.15       No     No      59
    ## 312 Female     297  53 164.0 102.0      Yes 24.50       No     No      95
    ## 313 Female     209  38 122.5  76.5      Yes 24.51       No     No      73
    ## 314 Female     323  39 131.5  85.0      Yes 24.79       No     No      93
    ## 315 Female     271  46 158.0  94.0      Yes 25.17       No     No      71
    ## 316 Female     255  38 123.5  84.5      Yes 25.33       No     No      88
    ## 317 Female     246  55 139.0  90.0      Yes 29.00       No     No     100
    ## 318 Female     235  39 196.0 116.0      Yes 29.70       No     No      87
    ## 319 Female     186  39 114.0  77.0      Yes 21.01       No     No      85
    ## 320 Female     289  58 156.5  85.0      Yes 25.46       No     No      86
    ## 321 Female     272  36 113.0  66.5      Yes 20.69       No     No      59
    ## 322 Female     200  61 187.0  95.5      Yes 21.57       No     No      64
    ## 323 Female     249  41 107.5  75.0      Yes 23.69       No     No      78
    ## 324 Female     234  45 189.0  87.0      Yes 23.10       No     No      90
    ## 325 Female     356  61 168.0  98.0      Yes 27.30       No     No     106
    ## 326 Female     233  49 158.0 102.0      Yes 25.31       No     No      72
    ## 327 Female     229  50 114.0  68.0      Yes 23.20       No     No      70
    ## 328 Female     200  36 121.5  72.5      Yes 23.09       No     No      75
    ## 329 Female     298  54 151.5  96.5      Yes 30.29       No     No      77
    ## 330 Female     332  60 185.0  86.0      Yes 25.73       No     No      73
    ## 331 Female     229  39 113.5  72.5      Yes 22.33       No     No      82
    ## 332 Female     246  51 111.5  74.0      Yes 25.09       No     No      78
    ## 333 Female     174  40 130.0  86.0      Yes 25.05       No     No      83
    ## 334 Female     173  44 121.5  69.0      Yes 23.72       No     No      77
    ## 335 Female     209  54 139.0  75.0      Yes 25.82       No     No      95
    ## 336 Female     240  44 109.0  71.0      Yes 23.75       No     No      83
    ## 337 Female     285  56 145.0 100.0      Yes 30.14       No     No      86
    ## 338 Female     229  57 126.5  90.0      Yes 26.00       No     No      58
    ## 339 Female     187  40 144.0  90.0      Yes 22.17       No     No      93
    ## 340 Female     175  34 117.5  73.5      Yes 22.15       No     No      75
    ## 341 Female     194  49 103.0  65.0      Yes 20.41       No     No      75
    ## 342 Female     155  40 121.0  86.0      Yes 23.16       No     No      59
    ## 343 Female     169  41 119.0  72.0      Yes 19.78       No     No      74
    ## 344 Female     259  43 120.0  87.0      Yes 19.88       No     No      71
    ## 345 Female     225  66 110.0  76.0      Yes 18.23       No     No      85
    ## 346 Female     250  57 117.5  71.0      Yes 23.84       No     No      75
    ## 347 Female     201  51 147.5  95.0      Yes 22.34       No     No      67
    ## 348 Female     246  42 120.0  70.0      Yes 19.42       No     No      78
    ## 349 Female     246  47 120.0  78.0      Yes 24.71       No     No      75
    ## 350 Female     203  43 110.0  71.5      Yes 24.56       No     No      65
    ## 351 Female     234  48 120.0  81.0      Yes 23.22       No     No      83
    ## 352 Female     206  49 107.0  74.0      Yes 20.23       No     No      83
    ## 353 Female     262  51 102.0  64.0      Yes 28.06       No     No      66
    ## 354 Female     270  43 145.5  82.0      Yes 21.10       No     No      87
    ## 355 Female     247  39 122.0  70.0      Yes 18.70       No     No      65
    ## 356 Female     191  39 119.0  78.0      Yes 20.93       No     No      73
    ## 357 Female     190  39 106.0  72.0      Yes 25.64       No     No      75
    ## 358 Female     289  55 141.0  83.5      Yes 24.99       No     No      76
    ## 359 Female     192  38 130.0  80.0      Yes 27.51       No     No      90
    ## 360 Female     250  43 110.0  70.0      Yes 21.14       No     No      85
    ## 361 Female     164  36 100.0  64.0      Yes 19.87       No     No      65
    ## 362 Female     195  45 111.0  79.0      Yes 23.22       No     No      85
    ## 363 Female     327  51 117.0  70.0      Yes 18.52       No     No      76
    ## 364 Female     249  55 109.0  66.5      Yes 24.79       No     No      85
    ## 365 Female     184  45 147.0  88.0      Yes 23.48       No     No      76
    ## 366 Female     215  58 119.5  73.0      Yes 29.86       No     No      93
    ## 367 Female     188  40 105.0  65.0      Yes 21.15       No     No      70
    ## 368 Female     236  40 135.0  83.0      Yes 23.48       No     No      90
    ## 369 Female     212  38 107.0  67.5      Yes 20.40       No     No      87
    ## 370 Female     342  64 128.0  71.0      Yes 20.52       No     No      62
    ## 371 Female     233  43 100.0  70.0      Yes 22.90       No     No      78
    ## 372 Female     280  55 144.0  79.0      Yes 19.50       No     No      75
    ## 373 Female     246  63 163.0  82.0      Yes 24.38       No     No     108
    ## 374 Female     235  61 154.0  94.0      Yes 22.91       No     No      66
    ## 375 Female     196  42 145.0  87.0      Yes 29.65       No     No      60
    ## 376 Female     326  67 157.5  78.0      Yes 24.63       No     No      77
    ## 377 Female     235  50 121.0  78.0      Yes 23.01       No     No      78
    ## 378 Female     365  47 127.0  76.0      Yes 24.44       No     No      80
    ## 379 Female     310  56 128.5  82.0      Yes 25.36       No     No      85
    ## 380 Female     298  62 137.0  85.0      Yes 26.73       No     No      87
    ## 381 Female     200  41 124.0  76.0      Yes 24.20       No     No      86
    ## 382 Female     287  50 147.5  87.5      Yes 30.36       No     No      72
    ## 383 Female     221  47 144.0  91.0      Yes 35.78       No     No      66
    ## 384 Female     280  39 152.0 104.0      Yes 24.22       No     No      82
    ## 385 Female     239  53 112.5  67.0      Yes 25.63       No     No      74
    ## 386 Female     300  51 128.0  78.0      Yes 26.69       No     No      97
    ## 387 Female     225  39 112.0  74.0      Yes 27.26       No     No      85
    ## 388 Female     161  52 180.0 114.0      Yes 32.52       No    Yes     104
    ## 389 Female     219  42 126.0  73.0      Yes 22.65       No     No      65
    ## 390 Female     188  48 170.0 110.0      Yes 26.03       No     No     118
    ## 391 Female     286  42 133.5  80.0      Yes 26.25       No     No      65
    ## 392 Female     193  42 129.0  91.5      Yes 27.78       No     No      74
    ## 393 Female     243  44 129.0  88.0      Yes 30.85       No     No      83
    ## 394 Female     213  63 182.0  92.0      Yes 26.87       No    Yes      63
    ## 395 Female     229  48 111.0  85.0      Yes 24.10       No     No      74
    ## 396 Female     289  57 142.0  83.0      Yes 35.17       No     No      72
    ## 397 Female     192  43 119.5  69.5      Yes 24.67       No     No      83
    ## 398 Female     249  52 112.0  75.0      Yes 22.54       No     No      71
    ## 399 Female     174  46 115.5  65.0      Yes 29.84       No     No      80
    ## 400 Female     228  46 110.0  80.0      Yes 19.74       No     No     127
    ## 401 Female     261  47 133.0  77.0      Yes 27.96       No     No     105
    ## 402 Female     226  39  95.0  59.0      Yes 22.88       No     No      83
    ## 403 Female     282  38 135.0  80.0      Yes 29.14       No     No      89
    ## 404 Female     165  43 113.0  79.0      Yes 28.96       No     No      72
    ## 405 Female     237  40 112.5  77.5      Yes 23.58       No     No      84
    ## 406 Female     253  48 105.0  59.0      Yes 19.42       No     No      83
    ## 407 Female     216  39 116.0  72.0      Yes 24.25       No     No      71
    ## 408 Female     166  37 112.0  73.5      Yes 21.64       No     No      93
    ## 409 Female     182  40  95.5  64.0      Yes 25.21       No     No      72
    ## 410 Female     237  40 130.0  72.0      Yes 23.54       No     No      80
    ## 411 Female     263  55 155.0  84.0      Yes 27.87       No     No      60
    ## 412 Female     185  42 123.0  74.0      Yes 24.41       No     No      92
    ## 413 Female     200  33 119.0  74.0      Yes 23.80       No     No      74
    ## 414 Female     309  53 130.0  86.0      Yes 22.37       No     No      80
    ## 415 Female     275  55 144.5  88.5      Yes 27.05       No     No      79
    ## 416 Female     180  39 113.0  73.0      Yes 17.65       No     No      73
    ## 417 Female     160  37 104.0  61.0      Yes 20.22       No     No     100
    ## 418 Female     254  62 167.5 102.5      Yes 27.15       No     No      83
    ## 419 Female     195  51 154.0  96.0      Yes 28.38       No     No      75
    ## 420 Female     225  58 146.0  77.0      Yes 24.60       No     No      53
    ## 421 Female     216  59 205.0  92.5      Yes 25.86       No     No      84
    ## 422 Female     143  40 125.5  80.0      Yes 21.99       No     No      95
    ## 423 Female     195  41 120.5  76.0      Yes 22.91       No     No      70
    ## 424 Female     164  37  96.5  67.0      Yes 24.99       No     No      67
    ## 425 Female     283  49 127.0  86.0      Yes 23.68       No     No      78
    ## 426 Female     154  37 106.0  59.5      Yes 22.71       No     No      50
    ## 427 Female     170  39 110.5  69.0      Yes 22.19       No     No     103
    ## 428 Female     200  40 116.0  69.0      Yes 23.90       No     No      78
    ## 429 Female     265  44 110.0  78.0      Yes 20.88       No     No      68
    ## 430 Female     212  42 115.0  72.0      Yes 23.72       No     No     100
    ## 431 Female     274  54 116.0  79.0      Yes 24.77       No     No      65
    ## 432 Female     240  42 127.0  90.0      Yes 38.54       No     No      69
    ## 433 Female     255  45 111.0  72.0      Yes 17.32       No     No      65
    ## 434 Female     250  45 130.0  80.0      Yes 20.24       No     No      86
    ## 435 Female     293  51 151.0  92.0      Yes 30.67       No     No      77
    ## 436 Female     232  52 115.0  80.0      Yes 28.79       No     No      68
    ## 437 Female     359  42 115.0  71.0      Yes 24.46       No     No      68
    ## 438 Female     174  48 154.0  84.0      Yes 31.76       No     No      90
    ## 439 Female     322  65 165.0  95.0      Yes 22.84       No     No      81
    ## 440 Female     230  45 116.0  75.0      Yes 21.35       No     No      77
    ## 441 Female     309  38 113.0  68.0      Yes 21.35       No     No      75
    ## 442 Female     167  37 118.0  76.0      Yes 19.61       No     No      67
    ## 443 Female     239  45 112.0  70.0      Yes 23.48       No     No      95
    ## 444 Female     235  41 107.5  68.0      Yes 21.10       No     No     113
    ## 445 Female     202  50 138.0  72.0      Yes 25.03       No     No      98
    ## 446 Female     178  38  96.0  67.0      Yes 20.40       No     No      82
    ## 447 Female     220  38 114.0  73.5      Yes 27.06       No     No      67
    ## 448 Female     268  45 100.0  70.0      Yes 23.45       No     No      87
    ## 449 Female     350  49 174.0  90.0      Yes 18.44       No    Yes      78
    ## 450 Female     219  42 120.0  70.0      Yes 24.10       No     No      73
    ## 451 Female     281  57 192.0 105.0      Yes 27.04       No    Yes      75
    ## 452 Female     240  50 176.5 115.0      Yes 27.71       No     No      83
    ## 453 Female     240  49 110.0  71.0      Yes 22.02       No     No      84
    ## 454 Female     270  39 110.0  78.0      Yes 22.00       No     No      68
    ## 455 Female     178  68 154.5  83.0      Yes 18.88       No     No      86
    ## 456 Female     274  37 197.5 125.0      Yes 43.48       No     No      94
    ## 457 Female     245  43 112.5  80.0      Yes 23.43       No     No      77
    ## 458 Female     205  44 120.0  83.5      Yes 24.30       No     No      77
    ## 459 Female     192  51 133.0  87.0      Yes 20.72       No     No      77
    ## 460 Female     252  48 143.0  81.0      Yes 24.00       No     No     101
    ## 461 Female     229  44 146.0  78.0      Yes 25.25       No     No      75
    ## 462 Female     168  50 120.0  80.0      Yes 25.26       No     No      60
    ## 463 Female     235  49 109.0  70.0      Yes 28.66       No     No      73
    ## 464 Female     281  42 115.5  79.0      Yes 22.90       No     No      71
    ## 465 Female     218  41 102.0  67.0      Yes 19.23       No     No      78
    ## 466 Female     241  64 136.5  85.0      Yes 26.42       No     No      77
    ## 467 Female     299  48 132.0  81.0      Yes 24.35       No     No      70
    ## 468 Female     212  44 128.0  81.5      Yes 27.51       No     No      87
    ## 469 Female     242  32 111.0  70.0      Yes 29.84       No     No      88
    ## 470 Female     390  54 150.0  94.0      Yes 27.34       No     No      71
    ## 471 Female     214  58 158.5  94.5      Yes 29.14       No     No      97
    ## 472 Female     186  38 166.0  96.0      Yes 33.14       No     No      75
    ## 473 Female     240  57 142.0  85.0      Yes 22.55       No     No      77
    ## 474 Female     280  41 129.0  89.0      Yes 39.69       No     No      65
    ## 475 Female     198  57 208.0 104.0      Yes 24.91       No     No      82
    ## 476 Female     342  59 137.0  83.5      Yes 25.18      Yes     No     140
    ## 477 Female     152  45 118.5  76.0      Yes 25.38       No     No      68
    ## 478 Female     326  41 124.0  78.0      Yes 25.27       No     No      67
    ## 479 Female     156  37 120.0  87.0      Yes 21.80       No     No      89
    ## 480 Female     252  49 123.0  69.0      Yes 21.45       No     No      89
    ## 481 Female     323  49 123.5  78.0      Yes 22.86       No     No      63
    ## 482 Female     288  44 150.0  89.0      Yes 21.11       No     No      97
    ## 483 Female     255  51 102.5  64.5      Yes 24.14       No     No      71
    ## 484 Female     159  34  92.5  70.0      Yes 22.15       No     No      68
    ## 485 Female     214  44 102.0  68.0      Yes 32.82       No     No      80
    ## 486 Female     248  52 128.0  74.0      Yes 21.84       No     No      79
    ## 487 Female     161  50 145.0  89.0      Yes 20.30       No     No      81
    ## 488 Female     325  56 160.0  97.5      Yes 23.40       No     No      86
    ## 489 Female     250  48 103.0  76.5      Yes 23.25       No     No      66
    ## 490 Female     213  63 172.0  95.0      Yes 27.68       No    Yes      67
    ## 491 Female     214  50 176.0 113.0      Yes 22.17       No     No      71
    ## 492 Female     236  39 127.0  78.0      Yes 17.51       No     No      76
    ## 493 Female     185  37 115.0  76.0      Yes 23.55       No     No      80
    ## 494 Female     261  40 112.0  67.0      Yes 21.83       No     No      61
    ## 495 Female     261  66 154.0  97.0      Yes 32.60       No    Yes      81
    ## 496 Female     271  48 130.0  84.0      Yes 21.97       No     No      85
    ## 497 Female     220  50 122.0  80.0      Yes 24.22       No     No      72
    ## 498 Female     200  37 112.5  68.0      Yes 25.87       No     No      67
    ## 499 Female     201  42 141.0  84.5      Yes 26.58       No     No      97
    ## 500 Female     371  67 166.0  85.0      Yes 25.35       No     No      86
    ## 501 Female     309  59 141.0  77.5      Yes 25.97       No     No      79
    ## 502 Female     312  59 175.0  82.0      Yes 39.82       No    Yes      85
    ## 503 Female     270  67 137.5  72.5      Yes 35.01       No     No      73
    ## 504 Female     232  46 115.0  70.0      Yes 25.18       No     No      59
    ## 505 Female     199  45 124.0  78.0      Yes 21.94       No     No      78
    ## 506 Female     258  68 158.0  94.0      Yes 31.64       No     No      84
    ## 507 Female     223  44  96.0  59.0      Yes 23.82       No     No      87
    ## 508 Female     206  39 102.5  65.0      Yes 19.80       No     No      85
    ## 509 Female     172  43 149.0  82.0      Yes 22.35       No    Yes      64
    ## 510 Female     222  37 110.0  71.0      Yes 18.30       No     No      67
    ## 511 Female     180  39 112.5  85.0      Yes 25.31       No     No      58
    ## 512 Female     266  48 155.0 100.0      Yes 27.86       No     No      84
    ## 513 Female     190  51 131.5  89.0      Yes 23.66       No     No     100
    ## 514 Female     167  41 147.5  87.5      Yes 32.52       No     No      80
    ## 515 Female     248  50 154.5 104.0      Yes 19.88       No     No      87
    ## 516 Female     204  36 112.0  78.0      Yes 28.74       No     No      82
    ## 517 Female     226  39 115.0  80.0      Yes 25.19       No     No      74
    ## 518 Female     315  55 123.0  77.5      Yes 26.21       No     No      84
    ## 519 Female     256  43 129.0  86.0      Yes 25.89       No     No      72
    ## 520 Female     295  44 114.0  85.0      Yes 23.10       No     No      84
    ## 521 Female     216  42 119.0  75.0      Yes 27.01       No     No      73
    ## 522 Female     214  36 107.5  66.0      Yes 22.02       No     No     103
    ## 523 Female     213  54 144.0  82.0      Yes 29.45       No     No      72
    ## 524 Female     195  38 116.0  72.0      Yes 24.45       No     No      90
    ## 525 Female     280  45 116.0  69.0      Yes 26.45       No     No      92
    ## 526 Female     288  56 123.0  70.0      Yes 18.62       No     No      96
    ## 527 Female     345  59 182.5 103.0      Yes 31.52       No     No      76
    ## 528 Female     173  42 105.0  70.0      Yes 21.98       No     No      79
    ## 529 Female     243  47 121.0  61.0      Yes 20.32       No     No     110
    ## 530 Female     155  40 120.0  73.5      Yes 21.32       No     No      78
    ## 531 Female     171  42 111.0  81.0      Yes 29.47       No     No      77
    ## 532 Female     209  40  96.0  67.0      Yes 23.46       No     No      70
    ## 533 Female     250  47  98.0  73.0      Yes 24.39       No     No      88
    ## 534 Female     177  36 115.0  63.5      Yes 22.54       No     No      73
    ## 535 Female     322  68 148.0  95.0      Yes 20.98       No     No      72
    ## 536 Female     288  49  99.0  60.0      Yes 22.19       No     No      76
    ## 537 Female     195  48 121.0  78.0      Yes 26.27       No     No      80
    ## 538 Female     260  41 101.0  68.0      Yes 22.49       No     No      77
    ## 539 Female     286  62 123.0  77.0      Yes 20.56       No     No      86
    ## 540 Female     231  44 103.0  73.0      Yes 22.02       No     No      83
    ## 541 Female     211  46 120.0  84.0      Yes 22.53       No     No      87
    ## 542 Female     232  40 140.0  92.0      Yes 26.56       No     No      73
    ## 543 Female     249  43 108.0  76.0      Yes 27.49       No     No      76
    ## 544 Female     200  45 117.5  83.5      Yes 23.68       No     No      73
    ## 545 Female     187  40 107.0  74.0      Yes 22.37       No     No      67
    ## 546 Female     210  45 121.0  82.0      Yes 23.08       No     No      71
    ## 547 Female     205  47 107.0  68.0      Yes 23.29       No     No     108
    ## 548 Female     209  41 107.0  65.0      Yes 27.27       No     No      87
    ## 549 Female     231  58 165.0  94.5      Yes 27.02       No     No      80
    ## 550 Female     240  40 115.0  72.0      Yes 18.82       No     No      68
    ## 551 Female     189  41 113.5  61.0      Yes 23.08       No     No      73
    ## 552 Female     179  38 116.5  72.5      Yes 21.49       No     No      76
    ## 553 Female     243  50 128.0  76.0      Yes 25.65       No     No      80
    ## 554 Female     160  47 197.0 109.0      Yes 34.91      Yes     No     320
    ## 555 Female     188  44 122.5  67.5      Yes 23.13       No     No      85
    ## 556 Female     193  40 105.0  60.0      Yes 22.54       No     No      85
    ## 557 Female     208  43 119.0  69.0      Yes 19.48       No     No      73
    ## 558 Female     240  55 107.5  70.0      Yes 18.06       No     No     140
    ## 559 Female     252  57 139.0  82.0      Yes 26.36       No     No      70
    ## 560 Female     230  45 138.0  85.0      Yes 22.54       No     No      67
    ## 561 Female     185  43 125.0  84.0      Yes 23.18       No     No      55
    ## 562 Female     215  48 114.0  64.0      Yes 21.51       No     No      64
    ## 563 Female     228  42 122.0  85.0      Yes 23.99       No     No      68
    ## 564 Female     247  48 142.5  86.5      Yes 24.87       No     No      68
    ## 565 Female     201  36 124.0  73.0      Yes 25.25       No     No     100
    ## 566 Female     224  38  90.0  70.0      Yes 18.18       No     No      57
    ## 567 Female     174  57 120.0  62.0      Yes 25.13       No     No      77
    ## 568 Female     229  59 127.5  76.0      Yes 23.65       No     No      70
    ## 569 Female     223  42 129.0  79.0      Yes 28.04       No     No     100
    ## 570 Female     230  43 118.0  70.5      Yes 26.24       No     No      67
    ## 571 Female     283  57 207.5 118.0      Yes 38.61       No     No      83
    ## 572 Female     280  40 115.0  81.0      Yes 21.32       No     No      84
    ## 573 Female     276  43 138.0  88.0      Yes 20.40       No     No      70
    ## 574 Female     237  47 110.0  77.0      Yes 25.62       No     No      83
    ## 575 Female     160  44 131.0  81.5      Yes 25.71       No     No      70
    ## 576 Female     262  49 147.5  97.5      Yes 24.96       No    Yes      67
    ## 577 Female     352  51 136.5  87.0      Yes 25.79       No     No      67
    ## 578 Female     280  48 105.0  85.0      Yes 25.50       No     No      79
    ## 579 Female     230  40 107.5  75.0      Yes 26.38       No     No      76
    ## 580 Female     240  58 150.0  84.0      Yes 26.85       No     No      94
    ## 581 Female     172  38  98.0  53.0      Yes 22.18       No     No      82
    ## 582 Female     216  55 125.0  80.0      Yes 27.18      Yes     No     244
    ## 583 Female     181  53 163.5  87.0      Yes 34.69       No     No      71
    ## 584 Female     165  35 106.0  64.0      Yes 19.14       No     No      70
    ## 585 Female     240  56 125.0  79.0      Yes 27.38       No     No      82
    ## 586 Female     273  49 147.0  89.0      Yes 24.26       No     No      62
    ## 587 Female     261  41 147.5  97.0      Yes 31.65       No     No      73
    ## 588 Female     290  46 131.0  84.0      Yes 18.28       No     No      68
    ## 589 Female     232  41 125.0  80.0      Yes 25.38       No     No      79
    ## 590 Female     310  44 141.0  76.0      Yes 20.52       No     No      90
    ## 591 Female     177  38 126.0  80.0      Yes 23.84       No     No      79
    ## 592 Female     270  60 130.0  72.5      Yes 20.84       No     No     102
    ## 593 Female     315  57 193.0 109.0      Yes 27.99       No    Yes      74
    ## 594 Female     351  44 106.0  65.0      Yes 25.34       No     No      69
    ## 595 Female     275  50 123.0  83.0      Yes 24.29       No     No      64
    ## 596 Female     196  46 114.0  75.0      Yes 21.01       No     No      69
    ## 597 Female     192  41 123.0  72.0      Yes 19.16       No     No      90
    ## 598 Female     269  43 139.0  96.0      Yes 24.38       No     No      71
    ## 599 Female     220  40 131.5  82.5      Yes 24.35       No     No      78
    ## 600 Female     214  49 172.0 111.0      Yes 40.51       No     No      70
    ## 601 Female     254  52 114.0  80.0      Yes 16.59       No     No      74
    ## 602 Female     303  46 115.0  78.0      Yes 22.50       No     No      67
    ## 603 Female     321  47 132.0  88.0      Yes 28.14       No     No      74
    ## 604 Female     256  40 124.0  79.0      Yes 24.23       No     No      81
    ## 605 Female     268  53 163.0  97.0      Yes 27.88       No     No      65
    ## 606 Female     248  48 114.5  70.0      Yes 27.69       No     No      75
    ## 607 Female     267  38 179.5  97.0      Yes 20.44       No     No      67
    ## 608 Female     254  60 177.0 101.0      Yes 23.27       No     No      79
    ## 609 Female     226  40 138.0  99.0      Yes 35.02       No     No      73
    ## 610 Female     275  58 140.0  78.0      Yes 19.18       No     No      74
    ## 611 Female     310  42 124.0  72.5      Yes 22.32       No     No      74
    ## 612 Female     193  38 107.0  73.0      Yes 20.73       No     No      72
    ## 613 Female     189  40 115.0  81.0      Yes 22.73       No     No     103
    ## 614 Female     311  45 117.5  76.0      Yes 26.27       No     No      67
    ## 615 Female     256  50 136.5  81.0      Yes 23.07       No     No      78
    ## 616 Female     188  39 113.0  81.0      Yes 26.44       No     No      87
    ## 617 Female     226  38 117.5  72.0      Yes 20.71       No     No      73
    ## 618 Female     235  48 120.0  81.0      Yes 23.36       No     No      80
    ## 619 Female     207  43  95.5  70.0      Yes 19.78       No     No      79
    ## 620 Female     220  38 105.0  69.0      Yes 24.69       No     No      87
    ## 621 Female     297  45 142.0  91.0      Yes 35.02       No     No      86
    ## 622 Female     277  49 120.0  80.0      Yes 19.72       No     No      75
    ## 623 Female     214  47 144.0  92.0      Yes 22.73      Yes     No      57
    ## 624 Female     225  42 111.0  71.0      Yes 23.43       No     No      85
    ## 625 Female     195  36 109.0  69.0      Yes 23.24       No     No      70
    ## 626 Female     258  50 123.0  70.0      Yes 19.72       No     No      71
    ## 627 Female     300  43 120.0  78.0      Yes 28.18       No     No     106
    ## 628 Female     172  41 108.0  73.0      Yes 22.50       No     No      73
    ## 629 Female     315  46 165.0  85.0      Yes 32.89       No     No      91
    ## 630 Female     200  43 115.0  80.0      Yes 26.66       No     No      81
    ## 631 Female     272  64 131.0  85.0      Yes 21.82       No     No      80
    ## 632 Female     285  45 132.5  97.5      Yes 24.74       No     No      77
    ## 633 Female     272  40 123.0  75.0      Yes 23.08       No     No      63
    ## 634 Female     259  48 129.0  81.0      Yes 21.08       No     No      65
    ## 635 Female     202  43 114.0  78.0      Yes 26.61       No     No      87
    ## 636 Female     239  64 143.0  84.0      Yes 20.06       No     No      73
    ## 637 Female     261  51 127.0  81.0      Yes 20.24       No     No      96
    ## 638 Female     224  40 106.0  72.0      Yes 23.59       No     No      71
    ## 639 Female     250  55 196.0 103.0      Yes 27.53       No     No      75
    ## 640 Female     305  48 118.5  73.0      Yes 20.99       No     No      84
    ## 641 Female     167  47 115.0  70.0      Yes 22.71       No     No     100
    ## 642 Female     159  41  99.0  62.0      Yes 19.09       No     No      67
    ## 643 Female     241  66 112.0  66.0      Yes 23.36       No     No      74
    ## 644 Female     194  39 112.5  77.5      Yes 21.51       No     No      84
    ## 645 Female     175  44 130.0  80.0      Yes 19.18       No     No     117
    ## 646 Female     205  41 120.0  80.0      Yes 20.67       No     No      64
    ## 647 Female     250  46 112.5  60.0      Yes 22.72       No     No      74
    ## 648 Female     267  46 119.0  65.0      Yes 29.15       No     No      75
    ## 649 Female     253  48 120.0  77.5      Yes 24.53       No     No      98
    ## 650 Female     165  40 101.0  59.0      Yes 23.06       No     No      76
    ## 651 Female     229  47 127.0  76.5      Yes 23.48       No     No      64
    ## 652 Female     200  42 119.0  75.0      Yes 22.91       No     No      69
    ## 653 Female     259  61 134.5  87.0      Yes 22.91       No     No      91
    ## 654 Female     243  41 159.0 100.0      Yes 27.78       No     No      71
    ## 655 Female     254  37 119.0  62.5      Yes 28.78       No     No      69
    ## 656 Female     198  45 119.0  80.0      Yes 22.18       No     No      79
    ## 657 Female     226  58 125.0  75.0      Yes 24.00       No     No      73
    ## 658 Female     173  40 125.0  75.0      Yes 25.67       No     No     102
    ## 659 Female     198  47 143.0  87.0      Yes 20.86       No     No      79
    ## 660 Female     254  59 154.5  93.0      Yes 21.82       No     No      89
    ## 661 Female     326  51 112.0  83.0      Yes 20.82       No     No      70
    ## 662 Female     259  64 195.0 110.0      Yes 20.09       No     No      63
    ## 663 Female     231  42 110.0  80.0      Yes 19.12       No     No      70
    ## 664 Female     221  47 116.5  81.0      Yes 25.85       No     No      75
    ## 665 Female     214  43 121.0  84.0      Yes 24.68       No     No      74
    ## 666 Female     382  57 140.0  94.0      Yes 21.20       No     No      70
    ## 667 Female     239  44 103.0  67.0      Yes 26.58       No     No      73
    ## 668 Female     240  51 112.0  83.0      Yes 24.10       No     No      77
    ## 669 Female     212  65 192.5 110.0      Yes 23.48       No     No      71
    ## 670 Female     346  49 130.0  80.0      Yes 22.54       No     No      77
    ## 671 Female     179  65 228.0 130.0      Yes 19.74       No    Yes      60
    ## 672 Female     256  53 128.0  90.0      Yes 23.65       No     No     102
    ## 673 Female     213  55 163.0  91.0      Yes 28.66       No     No      66
    ## 674 Female     266  56 114.0  72.0      Yes 22.64       No     No      83
    ## 675 Female     271  50 112.5  60.0      Yes 23.29       No     No      61
    ## 676 Female     180  36 118.0  80.0      Yes 29.59       No     No      84
    ## 677 Female     222  36 147.0  94.0      Yes 26.79       No     No      71
    ## 678 Female     271  58 146.0  92.0      Yes 23.07       No     No      83
    ## 679 Female     241  56 174.0  97.0      Yes 29.22      Yes    Yes     135
    ## 680 Female     231  61 128.0  87.0      Yes 26.30       No     No      93
    ## 681 Female     195  42 126.0  81.0      Yes 22.26       No     No      77
    ## 682 Female     206  58 159.5  93.5      Yes 18.53       No     No      58
    ## 683 Female     253  44 118.0  68.0      Yes 22.72       No     No     110
    ## 684 Female     238  50 158.0  74.0      Yes 35.68       No     No      98
    ## 685 Female     194  35 100.5  69.0      Yes 17.92       No     No      73
    ## 686 Female     168  35  83.5  55.0      Yes 16.71       No     No      63
    ## 687 Female     185  49 108.0  70.0      Yes 20.13       No     No      58
    ## 688 Female     242  40 112.5  62.5      Yes 27.65       No     No      70
    ## 689 Female     231  40 129.0  87.0      Yes 23.29       No     No      99
    ## 690 Female     258  43 161.5  96.0      Yes 38.96       No     No      84
    ## 691 Female     195  50 131.5  83.0      Yes 24.61       No     No      78
    ## 692 Female     201  38 123.5  78.0      Yes 27.14       No     No      77
    ## 693 Female     260  60 139.0  81.0      Yes 24.68       No     No      70
    ## 694 Female     246  48 113.0  87.0      Yes 18.01       No     No      63
    ## 695 Female     302  44 116.0  77.0      Yes 22.67       No     No      98
    ## 696 Female     275  46 126.0  71.0      Yes 24.91       No     No      71
    ## 697 Female     262  46 121.0  78.0      Yes 24.24       No     No      72
    ## 698 Female     235  43 128.5  80.0      Yes 18.83       No     No      70
    ## 699 Female     280  46 202.0 124.0      Yes 28.06       No    Yes      63
    ## 700 Female     265  55 154.0  87.0      Yes 20.92       No     No      66
    ## 701 Female     210  46 117.0  78.5      Yes 22.54       No     No     115
    ## 702 Female     222  53 123.0  82.0      Yes 25.52       No     No      67
    ## 703 Female     291  43 106.0  65.0      Yes 23.83       No     No      82
    ## 704 Female     276  43  99.0  62.0      Yes 22.17       No     No      80
    ## 705 Female     165  46  99.5  66.0      Yes 21.67       No     No      66
    ## 706 Female     190  39  85.0  70.0      Yes 22.43       No     No      60
    ## 707 Female     217  65 169.0 111.0      Yes 32.54       No     No      78
    ## 708 Female     279  59 116.0  66.0      Yes 21.83       No     No      79
    ## 709 Female     209  40 130.0  84.5      Yes 39.94       No     No     104
    ## 710 Female     202  36 105.5  67.0      Yes 22.66       No     No      63
    ## 711 Female     140  44 118.0  74.0      Yes 26.51       No     No      82
    ## 712 Female     268  45 130.0  94.0      Yes 34.27       No     No      93
    ## 713 Female     252  58 135.0  84.0      Yes 28.24       No     No      79
    ## 714 Female     175  47 107.0  69.0      Yes 23.64       No     No      70
    ## 715 Female     270  44 167.5  92.5      Yes 21.28       No     No      77
    ## 716 Female     226  47 122.5  80.0      Yes 24.62       No     No      68
    ## 717 Female     238  57 133.0  72.0      Yes 18.09       No     No     115
    ## 718 Female     220  51 137.0  79.0      Yes 21.66       No     No      74
    ## 719 Female     262  54 230.0 110.0      Yes 24.76       No     No      97
    ## 720 Female     180  40 118.0  86.5      Yes 22.69       No     No      67
    ## 721 Female     244  44 101.0  66.0      Yes 25.38       No     No      76
    ## 722 Female     231  44 133.0  89.0      Yes 29.29       No     No      83
    ## 723 Female     286  43 164.0  89.0      Yes 24.44       No     No      87
    ## 724 Female     226  41 125.5  82.0      Yes 23.80       No     No      75
    ## 725 Female     193  46 118.0  92.0      Yes 21.14       No     No      78
    ## 726 Female     275  37 118.0  71.0      Yes 23.10       No     No      95
    ## 727 Female     204  50 147.0 100.0      Yes 39.94       No     No      90
    ## 728 Female     178  40 142.0  84.0      Yes 34.46       No     No      77
    ## 729 Female     268  41 140.0  92.5      Yes 24.71       No     No      90
    ## 730 Female     164  38 113.0  68.0      Yes 25.75       No     No      75
    ## 731 Female     207  46 144.0  88.0      Yes 23.65       No     No      86
    ## 732 Female     257  37 141.0  93.0      Yes 41.29       No     No      58
    ## 733 Female     200  46 112.5  71.0      Yes 18.68       No     No      77
    ## 734 Female     243  41  97.0  63.0      Yes 22.53       No     No      64
    ## 735 Female     175  38 112.0  73.0      Yes 19.49       No     No      71
    ## 736 Female     240  54 113.0  73.0      Yes 24.21       No     No      77
    ## 737 Female     230  42 124.0  80.0      Yes 24.87       No     No      77
    ## 738 Female     211  47 159.5  82.5      Yes 34.08      Yes     No     250
    ## 739 Female     315  51 119.0  75.0      Yes 25.79       No     No      55
    ## 740 Female     250  43 123.0  74.0      Yes 26.01       No     No      90
    ## 741 Female     180  34 111.0  56.0      Yes 21.51       No     No      78
    ## 742 Female     212  36 127.0  83.0      Yes 26.82       No     No      75
    ## 743 Female     254  59 126.5  79.0      Yes 25.92       No     No      77
    ## 744 Female     308  47 138.0  94.0      Yes 24.33       No     No      90
    ## 745 Female     275  46 170.0 118.0      Yes 36.12       No     No      84
    ## 746 Female     270  57 120.0  79.0      Yes 24.83       No     No      81
    ## 747 Female     236  47 128.0  81.0      Yes 27.42       No     No      93
    ## 748 Female     214  40 109.5  69.0      Yes 20.32       No     No      81
    ## 749 Female     271  65 144.5  88.0      Yes 32.41      Yes     No     116
    ## 750 Female     265  56 150.0  84.0      Yes 28.66       No     No      73
    ## 751 Female     260  56 120.0  84.0      Yes 36.18       No     No      76
    ## 752 Female     194  45 115.0  85.0      Yes 27.77       No     No      82
    ## 753 Female     340  48 143.0  93.0      Yes 23.08       No     No      83
    ## 754 Female     165  50  96.5  62.5      Yes 23.48       No     No      78
    ## 755 Female     220  45 133.5  85.5      Yes 25.38       No     No      73
    ## 756 Female     254  57 146.5  81.0      Yes 41.61       No     No      85
    ## 757 Female     284  44 143.0  92.0      Yes 21.19       No     No      88
    ## 758 Female     244  55 145.5  94.0      Yes 28.86       No     No      72
    ## 759 Female     232  45 163.0 101.0      Yes 23.26       No     No      74
    ## 760 Female     274  46 158.0  97.0      Yes 22.83       No     No      78
    ## 761 Female     261  50 110.0  76.0      Yes 23.31       No     No      85
    ## 762 Female     246  60 160.0  92.0      Yes 26.38       No     No      73
    ## 763 Female     176  40  99.0  59.0      Yes 22.13       No     No      78
    ## 764 Female     286  61 141.0  81.0      Yes 23.61       No     No      52
    ## 765 Female     186  47 139.0  85.0      Yes 27.90      Yes     No     125
    ## 766 Female     220  42 129.0  81.0      Yes 19.74       No     No      61
    ## 767 Female     325  54 170.0 107.0      Yes 25.07       No    Yes      64
    ## 768 Female     205  46 115.0  75.0      Yes 19.48       No     No      78
    ## 769 Female     410  52 105.0  67.5      Yes 27.33       No     No      90
    ## 770 Female     350  49 135.0  86.5      Yes 25.56       No     No      83
    ## 771 Female     252  61 119.0  77.0      Yes 23.20       No     No      65
    ## 772 Female     235  41 129.0  94.0      Yes 23.71       No     No      81
    ## 773 Female     192  39 109.0  61.0      Yes 23.36       No     No      84
    ## 774 Female     188  52 130.0  71.0      Yes 23.88       No     No      89
    ## 775 Female     267  49 107.0  74.0      Yes 27.80       No     No      75
    ## 776 Female     358  62 215.0 110.0      Yes 37.62      Yes    Yes     368
    ## 777 Female     171  41 135.0  82.5      Yes 24.35       No     No      82
    ## 778 Female     233  46 106.0  60.0      Yes 20.84      Yes     No     348
    ## 779 Female     163  43 104.5  65.0      Yes 17.84       No     No      71
    ## 780 Female     204  43 133.0  86.5      Yes 26.01       No     No      79
    ## 781 Female     245  67 169.0  82.0      Yes 26.05       No    Yes     122
    ## 782 Female     214  56 147.0  65.0      Yes 17.68       No     No      87
    ## 783 Female     167  38 105.0  70.0      Yes 19.76       No     No      80
    ## 784 Female     200  39 110.0  68.0      Yes 20.24       No     No      62
    ## 785 Female     263  46 110.0  65.0      Yes 27.27       No     No      73
    ## 786 Female     222  44 130.0  86.0      Yes 27.42       No     No      84
    ## 787 Female     164  38 110.0  76.0      Yes 23.85       No     No      83
    ## 788 Female     292  52 157.0 112.0      Yes 29.56       No     No      84
    ## 789 Female     185  35 131.5  84.0      Yes 20.32       No     No      76
    ## 790 Female     273  56 130.5  82.0      Yes 25.48       No     No      91
    ## 791 Female     275  51 150.0  99.0      Yes 23.17       No     No      65
    ## 792 Female     214  43 133.0  86.0      Yes 22.72       No     No      77
    ## 793 Female     227  49 150.0  91.0      Yes 24.30       No     No      83
    ## 794 Female     260  47 126.5  81.0      Yes 26.58       No     No      82
    ## 795 Female     202  45  93.5  58.0      Yes 21.25       No     No      60
    ## 796 Female     197  39 126.5  76.5      Yes 19.71       No     No      63
    ## 797 Female     174  44 174.0 130.0      Yes 33.99       No     No      63
    ## 798 Female     205  54  96.0  66.0      Yes 23.26       No     No      75
    ## 799 Female     250  43 112.5  76.5      Yes 25.23       No     No      63
    ## 800 Female     202  46 157.0  94.0      Yes 19.37       No     No      65
    ## 801 Female     188  41 145.0  99.0      Yes 28.60       No     No      74
    ## 802 Female     185  35 100.0  66.5      Yes 24.08       No     No      75
    ## 803 Female     179  43 101.0  68.5      Yes 19.83       No     No      76
    ## 804 Female     246  54 153.0  80.0      Yes 37.30       No     No      74
    ## 805 Female     205  44 132.5  82.5      Yes 30.98       No     No      66
    ## 806 Female     205  37 111.0  60.5      Yes 21.80       No     No      82
    ## 807 Female     332  46 162.5  92.5      Yes 26.13       No     No      67
    ## 808 Female     233  49 112.5  80.0      Yes 27.87       No     No      80
    ## 809 Female     213  40 130.0  80.0      Yes 19.98       No     No      76
    ## 810 Female     205  52 159.0 110.0      Yes 28.18       No     No      83
    ## 811 Female     259  46 129.0  83.0      Yes 22.91       No     No      84
    ## 812 Female     272  57 157.0  80.0      Yes 25.15       No    Yes      95
    ## 813 Female     262  56 126.0  74.0      Yes 27.35       No     No     115
    ## 814 Female     221  40  93.0  62.5      Yes 18.84       No     No      73
    ## 815 Female     217  56 134.0  75.0      Yes 29.59       No     No      92
    ## 816 Female     282  46 176.0  98.0      Yes 33.02       No    Yes      78
    ## 817 Female     199  43 104.0  79.0      Yes 20.12       No     No      64
    ## 818 Female     270  64 142.0  68.0      Yes 21.32       No     No      80
    ## 819 Female     255  60 140.0  77.0      Yes 36.29       No    Yes     100
    ## 820 Female     291  45 125.0  82.0      Yes 21.26       No     No      72
    ## 821 Female     235  42 128.0  86.0      Yes 24.05       No     No      70
    ## 822 Female     186  56 126.0  73.0      Yes 25.81       No     No      65
    ## 823 Female     202  48 128.0  74.0      Yes 25.11       No     No      75
    ## 824 Female     257  36 103.0  72.5      Yes 27.86       No     No      65
    ## 825 Female     224  37 109.0  72.0      Yes 22.81       No     No      93
    ## 826 Female     214  38 101.0  70.0      Yes 21.83       No     No      77
    ## 827 Female     218  52 121.5  57.0      Yes 20.78       No     No      85
    ## 828 Female     215  49 106.0  63.0      Yes 19.22       No     No      66
    ## 829 Female     199  42 143.0  79.0      Yes 18.68       No     No      76
    ## 830 Female     173  37 101.0  69.0      Yes 20.02       No     No      73
    ## 831 Female     207  40 124.0  78.0      Yes 22.90       No     No      66
    ## 832 Female     240  52 157.5 105.0      Yes 29.64       No     No      80
    ## 833 Female     190  43 103.0  67.5      Yes 24.08       No     No      69
    ## 834 Female     230  56 123.0  78.5      Yes 24.71       No     No      87
    ## 835 Female     251  59 125.0  80.0      Yes 22.18       No     No      70
    ## 836 Female     196  39 133.0  86.0      Yes 20.91       No     No      80
    ##     EDUC  MI
    ## 1      3  No
    ## 2      3  No
    ## 3      2  No
    ## 4      2  No
    ## 5      2  No
    ## 6      2  No
    ## 7      2  No
    ## 8      3  No
    ## 9      1  No
    ## 10     2  No
    ## 11     2  No
    ## 12     2  No
    ## 13     1  No
    ## 14     2  No
    ## 15     1  No
    ## 16     1  No
    ## 17     4  No
    ## 18     1 Yes
    ## 19     2  No
    ## 20     1  No
    ## 21     2  No
    ## 22     1  No
    ## 23     2  No
    ## 24     2  No
    ## 25     4  No
    ## 26     1  No
    ## 27     2  No
    ## 28     3  No
    ## 29     1  No
    ## 30     2  No
    ## 31     1  No
    ## 32     2  No
    ## 33     3  No
    ## 34     3  No
    ## 35     1  No
    ## 36     1 Yes
    ## 37     3  No
    ## 38     4  No
    ## 39     2  No
    ## 40     2  No
    ## 41     3  No
    ## 42     1  No
    ## 43     1  No
    ## 44     1  No
    ## 45     1 Yes
    ## 46     1  No
    ## 47     4  No
    ## 48     2  No
    ## 49     2  No
    ## 50     2  No
    ## 51     4  No
    ## 52     1  No
    ## 53     2  No
    ## 54     2  No
    ## 55     1  No
    ## 56     3  No
    ## 57     3  No
    ## 58     2  No
    ## 59     1  No
    ## 60     2  No
    ## 61     2  No
    ## 62     2  No
    ## 63     4  No
    ## 64     3  No
    ## 65     2  No
    ## 66     3  No
    ## 67     1  No
    ## 68     2  No
    ## 69     3  No
    ## 70     1  No
    ## 71     1  No
    ## 72     3 Yes
    ## 73     1  No
    ## 74     4  No
    ## 75     1 Yes
    ## 76     2  No
    ## 77     2  No
    ## 78     3  No
    ## 79     3  No
    ## 80     1  No
    ## 81     2  No
    ## 82     3  No
    ## 83     2  No
    ## 84     2  No
    ## 85     2  No
    ## 86     2  No
    ## 87     1  No
    ## 88     1  No
    ## 89     2  No
    ## 90     2  No
    ## 91     1  No
    ## 92     3  No
    ## 93     4  No
    ## 94     3  No
    ## 95     2 Yes
    ## 96     2  No
    ## 97     4  No
    ## 98     2  No
    ## 99     2  No
    ## 100    3  No
    ## 101    2  No
    ## 102    3  No
    ## 103    2  No
    ## 104    4  No
    ## 105    1  No
    ## 106    3  No
    ## 107    3  No
    ## 108    3  No
    ## 109    3  No
    ## 110    3  No
    ## 111    3  No
    ## 112    3  No
    ## 113    1  No
    ## 114    3  No
    ## 115    1  No
    ## 116    1  No
    ## 117    3  No
    ## 118    4  No
    ## 119    3 Yes
    ## 120    3  No
    ## 121    3  No
    ## 122    2  No
    ## 123    4  No
    ## 124    1  No
    ## 125    2  No
    ## 126    1  No
    ## 127    3  No
    ## 128    3 Yes
    ## 129    2  No
    ## 130    1 Yes
    ## 131    2 Yes
    ## 132    2 Yes
    ## 133    2  No
    ## 134    3  No
    ## 135    1 Yes
    ## 136    2  No
    ## 137    1  No
    ## 138    2  No
    ## 139    2  No
    ## 140    1  No
    ## 141    2  No
    ## 142    4  No
    ## 143    2  No
    ## 144    1  No
    ## 145    3  No
    ## 146    1  No
    ## 147    1  No
    ## 148    1  No
    ## 149    2 Yes
    ## 150    4 Yes
    ## 151    1  No
    ## 152    3  No
    ## 153    3  No
    ## 154    2  No
    ## 155    1 Yes
    ## 156    3  No
    ## 157    1 Yes
    ## 158    2  No
    ## 159    2  No
    ## 160    2  No
    ## 161    1  No
    ## 162    1  No
    ## 163    1  No
    ## 164    4  No
    ## 165    3  No
    ## 166    2  No
    ## 167    1  No
    ## 168    3  No
    ## 169    3  No
    ## 170    2  No
    ## 171    2  No
    ## 172    2 Yes
    ## 173    1  No
    ## 174    2  No
    ## 175    2  No
    ## 176    1  No
    ## 177    1  No
    ## 178    2  No
    ## 179    2  No
    ## 180    2  No
    ## 181    2  No
    ## 182    2  No
    ## 183    2  No
    ## 184    2  No
    ## 185    2  No
    ## 186    1  No
    ## 187    1 Yes
    ## 188    1  No
    ## 189    3  No
    ## 190    2  No
    ## 191    1  No
    ## 192    4  No
    ## 193    3  No
    ## 194    1  No
    ## 195    1  No
    ## 196    1  No
    ## 197    3  No
    ## 198    1  No
    ## 199    1  No
    ## 200    2  No
    ## 201    2  No
    ## 202    1  No
    ## 203    3  No
    ## 204    2  No
    ## 205    1 Yes
    ## 206    2  No
    ## 207    1  No
    ## 208    2  No
    ## 209    1  No
    ## 210    2  No
    ## 211    3  No
    ## 212    2  No
    ## 213    1  No
    ## 214    3  No
    ## 215    2  No
    ## 216    2  No
    ## 217    3 Yes
    ## 218    1 Yes
    ## 219    1  No
    ## 220    1  No
    ## 221    2  No
    ## 222    3  No
    ## 223    1  No
    ## 224    3  No
    ## 225    2  No
    ## 226    1 Yes
    ## 227    2  No
    ## 228    1  No
    ## 229    1  No
    ## 230    1 Yes
    ## 231    4  No
    ## 232    1  No
    ## 233    1  No
    ## 234    3  No
    ## 235    3  No
    ## 236    2 Yes
    ## 237    4  No
    ## 238    1  No
    ## 239    3  No
    ## 240    1  No
    ## 241    3  No
    ## 242    2  No
    ## 243    2 Yes
    ## 244    2 Yes
    ## 245    1  No
    ## 246    1  No
    ## 247    4  No
    ## 248    1  No
    ## 249    3  No
    ## 250    1  No
    ## 251    2  No
    ## 252    2  No
    ## 253    2  No
    ## 254    2  No
    ## 255    2  No
    ## 256    4  No
    ## 257    3  No
    ## 258    2  No
    ## 259    2  No
    ## 260    1  No
    ## 261    1  No
    ## 262    3  No
    ## 263    1  No
    ## 264    1  No
    ## 265    4  No
    ## 266    3  No
    ## 267    2  No
    ## 268    2  No
    ## 269    1  No
    ## 270    3  No
    ## 271    1 Yes
    ## 272    4  No
    ## 273    3  No
    ## 274    3  No
    ## 275    2  No
    ## 276    2  No
    ## 277    1 Yes
    ## 278    3  No
    ## 279    3  No
    ## 280    3  No
    ## 281    2  No
    ## 282    2  No
    ## 283    3  No
    ## 284    2  No
    ## 285    1  No
    ## 286    1  No
    ## 287    3  No
    ## 288    1  No
    ## 289    2  No
    ## 290    2  No
    ## 291    2  No
    ## 292    3  No
    ## 293    3  No
    ## 294    3  No
    ## 295    1  No
    ## 296    2  No
    ## 297    2 Yes
    ## 298    2  No
    ## 299    3  No
    ## 300    3  No
    ## 301    4  No
    ## 302    1  No
    ## 303    2  No
    ## 304    1  No
    ## 305    3  No
    ## 306    3  No
    ## 307    2  No
    ## 308    1  No
    ## 309    1  No
    ## 310    2  No
    ## 311    1  No
    ## 312    4  No
    ## 313    2  No
    ## 314    2 Yes
    ## 315    2  No
    ## 316    2  No
    ## 317    2  No
    ## 318    3 Yes
    ## 319    3  No
    ## 320    1  No
    ## 321    2  No
    ## 322    1 Yes
    ## 323    1  No
    ## 324    2  No
    ## 325    1  No
    ## 326    2  No
    ## 327    1  No
    ## 328    3  No
    ## 329    1  No
    ## 330    1 Yes
    ## 331    4  No
    ## 332    1  No
    ## 333    2  No
    ## 334    1  No
    ## 335    1  No
    ## 336    3  No
    ## 337    1  No
    ## 338    1  No
    ## 339    3  No
    ## 340    1  No
    ## 341    4 Yes
    ## 342    4  No
    ## 343    2  No
    ## 344    1  No
    ## 345    2  No
    ## 346    4  No
    ## 347    1  No
    ## 348    2  No
    ## 349    3  No
    ## 350    4  No
    ## 351    3  No
    ## 352    2  No
    ## 353    1  No
    ## 354    3  No
    ## 355    2  No
    ## 356    4  No
    ## 357    3  No
    ## 358    2  No
    ## 359    3  No
    ## 360    3  No
    ## 361    1  No
    ## 362    3  No
    ## 363    1  No
    ## 364    3  No
    ## 365    3  No
    ## 366    1 Yes
    ## 367    4  No
    ## 368    3  No
    ## 369    1  No
    ## 370    4 Yes
    ## 371    3  No
    ## 372    1  No
    ## 373    2  No
    ## 374    2  No
    ## 375    2  No
    ## 376    2  No
    ## 377    2  No
    ## 378    2  No
    ## 379    2  No
    ## 380    4  No
    ## 381    2  No
    ## 382    1 Yes
    ## 383    1  No
    ## 384    1  No
    ## 385    1  No
    ## 386    1  No
    ## 387    2  No
    ## 388    1 Yes
    ## 389    2  No
    ## 390    1  No
    ## 391    2  No
    ## 392    3  No
    ## 393    2  No
    ## 394    1  No
    ## 395    1  No
    ## 396    3  No
    ## 397    3  No
    ## 398    3  No
    ## 399    2  No
    ## 400    4  No
    ## 401    1  No
    ## 402    1  No
    ## 403    4  No
    ## 404    2  No
    ## 405    3  No
    ## 406    3  No
    ## 407    2  No
    ## 408    1  No
    ## 409    4  No
    ## 410    1  No
    ## 411    4  No
    ## 412    2 Yes
    ## 413    2  No
    ## 414    1  No
    ## 415    3  No
    ## 416    4  No
    ## 417    4  No
    ## 418    4  No
    ## 419    1  No
    ## 420    1  No
    ## 421    3  No
    ## 422    3  No
    ## 423    3  No
    ## 424    2  No
    ## 425    2  No
    ## 426    2  No
    ## 427    4  No
    ## 428    4  No
    ## 429    4  No
    ## 430    3  No
    ## 431    1  No
    ## 432    1  No
    ## 433    2  No
    ## 434    2 Yes
    ## 435    2  No
    ## 436    1  No
    ## 437    2  No
    ## 438    1 Yes
    ## 439    4  No
    ## 440    1  No
    ## 441    2  No
    ## 442    2  No
    ## 443    4  No
    ## 444    2  No
    ## 445    2  No
    ## 446    4  No
    ## 447    1  No
    ## 448    3  No
    ## 449    2  No
    ## 450    4  No
    ## 451    1 Yes
    ## 452    3  No
    ## 453    2  No
    ## 454    2  No
    ## 455    1  No
    ## 456    1 Yes
    ## 457    1  No
    ## 458    1  No
    ## 459    2  No
    ## 460    1  No
    ## 461    1  No
    ## 462    1  No
    ## 463    1  No
    ## 464    2  No
    ## 465    4  No
    ## 466    1  No
    ## 467    4  No
    ## 468    1  No
    ## 469    2  No
    ## 470    2 Yes
    ## 471    1 Yes
    ## 472    2  No
    ## 473    1  No
    ## 474    2  No
    ## 475    1  No
    ## 476    3 Yes
    ## 477    3  No
    ## 478    1 Yes
    ## 479    4  No
    ## 480    1  No
    ## 481    1  No
    ## 482    2  No
    ## 483    2  No
    ## 484    2  No
    ## 485    1  No
    ## 486    3  No
    ## 487    2  No
    ## 488    3 Yes
    ## 489    2  No
    ## 490    2  No
    ## 491    4  No
    ## 492    2  No
    ## 493    2  No
    ## 494    2  No
    ## 495    1  No
    ## 496    2  No
    ## 497    2  No
    ## 498    1  No
    ## 499    2  No
    ## 500    2  No
    ## 501    2  No
    ## 502    1  No
    ## 503    1  No
    ## 504    2  No
    ## 505    2  No
    ## 506    3  No
    ## 507    2  No
    ## 508    2  No
    ## 509    2  No
    ## 510    3  No
    ## 511    2  No
    ## 512    3  No
    ## 513    2  No
    ## 514    1  No
    ## 515    1 Yes
    ## 516    2  No
    ## 517    4  No
    ## 518    2  No
    ## 519    1  No
    ## 520    2  No
    ## 521    2  No
    ## 522    3  No
    ## 523    3  No
    ## 524    2  No
    ## 525    2 Yes
    ## 526    2  No
    ## 527    2 Yes
    ## 528    2  No
    ## 529    4  No
    ## 530    1  No
    ## 531    1  No
    ## 532    1  No
    ## 533    3  No
    ## 534    2  No
    ## 535    2  No
    ## 536    1  No
    ## 537    1  No
    ## 538    2  No
    ## 539    2  No
    ## 540    1  No
    ## 541    2  No
    ## 542    2  No
    ## 543    1  No
    ## 544    2 Yes
    ## 545    1  No
    ## 546    3  No
    ## 547    4 Yes
    ## 548    1  No
    ## 549    2 Yes
    ## 550    2  No
    ## 551    2  No
    ## 552    2  No
    ## 553    4  No
    ## 554    1  No
    ## 555    3  No
    ## 556    2  No
    ## 557    2  No
    ## 558    1  No
    ## 559    1  No
    ## 560    3  No
    ## 561    3  No
    ## 562    3  No
    ## 563    1  No
    ## 564    4  No
    ## 565    2  No
    ## 566    1  No
    ## 567    1  No
    ## 568    2  No
    ## 569    2  No
    ## 570    2  No
    ## 571    1  No
    ## 572    3  No
    ## 573    3 Yes
    ## 574    3  No
    ## 575    1  No
    ## 576    4  No
    ## 577    1  No
    ## 578    4  No
    ## 579    4  No
    ## 580    1  No
    ## 581    2  No
    ## 582    3  No
    ## 583    1  No
    ## 584    2  No
    ## 585    1  No
    ## 586    2 Yes
    ## 587    4 Yes
    ## 588    2  No
    ## 589    3  No
    ## 590    1  No
    ## 591    4  No
    ## 592    3  No
    ## 593    2  No
    ## 594    1  No
    ## 595    2  No
    ## 596    3  No
    ## 597    1  No
    ## 598    1  No
    ## 599    1  No
    ## 600    1  No
    ## 601    4  No
    ## 602    1  No
    ## 603    3  No
    ## 604    2 Yes
    ## 605    3  No
    ## 606    3  No
    ## 607    4  No
    ## 608    1  No
    ## 609    1 Yes
    ## 610    1  No
    ## 611    2  No
    ## 612    2  No
    ## 613    3  No
    ## 614    2 Yes
    ## 615    4  No
    ## 616    2  No
    ## 617    2  No
    ## 618    1  No
    ## 619    4  No
    ## 620    1  No
    ## 621    2  No
    ## 622    1  No
    ## 623    1  No
    ## 624    1  No
    ## 625    1  No
    ## 626    2  No
    ## 627    2 Yes
    ## 628    1  No
    ## 629    3  No
    ## 630    1  No
    ## 631    1  No
    ## 632    2 Yes
    ## 633    1  No
    ## 634    3  No
    ## 635    3  No
    ## 636    4  No
    ## 637    1  No
    ## 638    2  No
    ## 639    3  No
    ## 640    1  No
    ## 641    3  No
    ## 642    2  No
    ## 643    3  No
    ## 644    2  No
    ## 645    2  No
    ## 646    3  No
    ## 647    1  No
    ## 648    3 Yes
    ## 649    1  No
    ## 650    2  No
    ## 651    2  No
    ## 652    1  No
    ## 653    1  No
    ## 654    2  No
    ## 655    3  No
    ## 656    3  No
    ## 657    2  No
    ## 658    1  No
    ## 659    2  No
    ## 660    2  No
    ## 661    2 Yes
    ## 662    2 Yes
    ## 663    1  No
    ## 664    3  No
    ## 665    3  No
    ## 666    4  No
    ## 667    2  No
    ## 668    4  No
    ## 669    2 Yes
    ## 670    2  No
    ## 671    1  No
    ## 672    2  No
    ## 673    1  No
    ## 674    4  No
    ## 675    3  No
    ## 676    1  No
    ## 677    1  No
    ## 678    1  No
    ## 679    1  No
    ## 680    2  No
    ## 681    2  No
    ## 682    2  No
    ## 683    4  No
    ## 684    1 Yes
    ## 685    2  No
    ## 686    2  No
    ## 687    1  No
    ## 688    2  No
    ## 689    1  No
    ## 690    1 Yes
    ## 691    1  No
    ## 692    1  No
    ## 693    2  No
    ## 694    2 Yes
    ## 695    1  No
    ## 696    2  No
    ## 697    1  No
    ## 698    4  No
    ## 699    1  No
    ## 700    2  No
    ## 701    1  No
    ## 702    1  No
    ## 703    2  No
    ## 704    2  No
    ## 705    3  No
    ## 706    3  No
    ## 707    1  No
    ## 708    3  No
    ## 709    1  No
    ## 710    3  No
    ## 711    1  No
    ## 712    1  No
    ## 713    1  No
    ## 714    2  No
    ## 715    2 Yes
    ## 716    1  No
    ## 717    2  No
    ## 718    2  No
    ## 719    2  No
    ## 720    3  No
    ## 721    2  No
    ## 722    4  No
    ## 723    2  No
    ## 724    2  No
    ## 725    1  No
    ## 726    1  No
    ## 727    1  No
    ## 728    2  No
    ## 729    3  No
    ## 730    1  No
    ## 731    2  No
    ## 732    2  No
    ## 733    3  No
    ## 734    4  No
    ## 735    3  No
    ## 736    1  No
    ## 737    2  No
    ## 738    2  No
    ## 739    2  No
    ## 740    2  No
    ## 741    2  No
    ## 742    4  No
    ## 743    1  No
    ## 744    2  No
    ## 745    1 Yes
    ## 746    2  No
    ## 747    2  No
    ## 748    1  No
    ## 749    1  No
    ## 750    1  No
    ## 751    3  No
    ## 752    1  No
    ## 753    4  No
    ## 754    4  No
    ## 755    2  No
    ## 756    1 Yes
    ## 757    3  No
    ## 758    1  No
    ## 759    2 Yes
    ## 760    3  No
    ## 761    4  No
    ## 762    1  No
    ## 763    2  No
    ## 764    1  No
    ## 765    1  No
    ## 766    2  No
    ## 767    2  No
    ## 768    3  No
    ## 769    2 Yes
    ## 770    1 Yes
    ## 771    1  No
    ## 772    2  No
    ## 773    3  No
    ## 774    1  No
    ## 775    3  No
    ## 776    3  No
    ## 777    1  No
    ## 778    2  No
    ## 779    2  No
    ## 780    1 Yes
    ## 781    3  No
    ## 782    1  No
    ## 783    3  No
    ## 784    4  No
    ## 785    1  No
    ## 786    2  No
    ## 787    3  No
    ## 788    2  No
    ## 789    2  No
    ## 790    2  No
    ## 791    3  No
    ## 792    4  No
    ## 793    3  No
    ## 794    1  No
    ## 795    4  No
    ## 796    3  No
    ## 797    2 Yes
    ## 798    1  No
    ## 799    2  No
    ## 800    1  No
    ## 801    2  No
    ## 802    4  No
    ## 803    1  No
    ## 804    2  No
    ## 805    1  No
    ## 806    1  No
    ## 807    3 Yes
    ## 808    1  No
    ## 809    2  No
    ## 810    2  No
    ## 811    3  No
    ## 812    2 Yes
    ## 813    1  No
    ## 814    2  No
    ## 815    1  No
    ## 816    2  No
    ## 817    2  No
    ## 818    3 Yes
    ## 819    1  No
    ## 820    1  No
    ## 821    2  No
    ## 822    4  No
    ## 823    1  No
    ## 824    2  No
    ## 825    3  No
    ## 826    1  No
    ## 827    2  No
    ## 828    3  No
    ## 829    2  No
    ## 830    3  No
    ## 831    2  No
    ## 832    1  No
    ## 833    2  No
    ## 834    3  No
    ## 835    3  No
    ## 836    3  No

-----

The pipe operator `%>%` is commonly used by `dplyr` and pipes output
from one function to another. Instead of nesting functions, which can be
difficult to read at a glance, `dplyr` instead pipes them.

Here, we pipe our data frame to the `filter()` function, and then show
the result using `head`:

``` r
fhs %>%
  filter(CURSMOKE == 'Yes', BMI < 20) %>%
  head
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1 Female     237  47 130.0    78      Yes 19.66       No     No      75
    ## 2 Female     266  54 114.0    76      Yes 17.61      Yes     No      55
    ## 3   Male     220  53 123.5    75      Yes 19.64       No     No      73
    ## 4   Male     188  53 138.0    89      Yes 18.23       No     No      75
    ## 5 Female     185  37 100.0    68      Yes 18.38       No     No      72
    ## 6 Female     200  39 111.0    64      Yes 19.24       No     No      60
    ##   EDUC MI
    ## 1    2 No
    ## 2    1 No
    ## 3    1 No
    ## 4    1 No
    ## 5    2 No
    ## 6    3 No

-----

To re-order rows by a particular columns values, you can use the
`arrange()` function. Perhaps you want to have a look at the youngest
individuals in your data frame.

``` r
fhs %>%
  arrange(AGE) %>%
  head
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1 Female     242  32 111.0    70      Yes 29.84       No     No      88
    ## 2 Female     199  33 116.0    81      Yes 21.61       No     No      93
    ## 3 Female     158  33 108.0    67       No 19.84       No     No      69
    ## 4   Male     165  33 141.5    95       No 26.74       No     No      77
    ## 5 Female     200  33 119.0    74      Yes 23.80       No     No      74
    ## 6   Male     165  33 136.0    75       No 24.95       No     No      90
    ##   EDUC MI
    ## 1    2 No
    ## 2    3 No
    ## 3    1 No
    ## 4    2 No
    ## 5    2 No
    ## 6    4 No

You can look at the oldest observations by using the `desc()` function
inside `arrange()`.

``` r
fhs %>%
  arrange(desc(AGE)) %>%
  head
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1 Female     231  70   136    84       No 31.78       No    Yes      95
    ## 2   Male     280  69   121    71       No 29.23       No     No     161
    ## 3   Male     168  69   110    68      Yes 19.15       No    Yes      54
    ## 4 Female     286  69   117    73       No 20.92       No     No     103
    ## 5   Male     232  69   151    74      Yes 24.14       No     No      62
    ## 6 Female     203  69   166    90       No 25.40       No    Yes      80
    ##   EDUC  MI
    ## 1    2 Yes
    ## 2    1 Yes
    ## 3    4 Yes
    ## 4    3  No
    ## 5    1  No
    ## 6    3  No

-----

Sometimes, you would like to create new columns as a function of the
existing ones. In this instance, `mutate()` comes in handy.

To create a column for pulse pressure, we use the difference between
systolic and diastolic blood pressure.

``` r
fhs %>%
  mutate(PP = SYSBP - DIABP) %>%
  select(PP, SYSBP, DIABP) %>%
  arrange(desc(PP)) %>%
  head
```

    ##    PP SYSBP DIABP
    ## 1 160   295   135
    ## 2 130   193    63
    ## 3 125   197    72
    ## 4 124   220    96
    ## 5 120   244   124
    ## 6 120   230   110

Combining `mutate()` with `ifelse()` can help create helpful categorical
variables.

``` r
fhs %>%
  filter(DIABETES == 'Yes') %>%
  mutate(HYPO = ifelse(GLUCOSE <= 70,1,0)) %>%
  select(DIABETES, GLUCOSE, HYPO) %>%
  head
```

    ##   DIABETES GLUCOSE HYPO
    ## 1      Yes     225    0
    ## 2      Yes     215    0
    ## 3      Yes      55    1
    ## 4      Yes     202    0
    ## 5      Yes     195    0
    ## 6      Yes     126    0

Hopefully, you can start to see how using `dplyr` can help you get to
know your data.

-----

As well as manipulating data, you may also want to `summarise()` it.
This can make exploring the distribution of variables very intuitive.

``` r
fhs %>%
  summarise(TC_AVG = mean(TOTCHOL),
            TC_MIN = min(TOTCHOL),
            TC_MAX = max(TOTCHOL),
            TOTAL_N = n())
```

    ##     TC_AVG TC_MIN TC_MAX TOTAL_N
    ## 1 237.1641    113    600    3851

Multiple summary statistics are available, including `sd()`, `median()`,
and `sum()`.

It can be really useful to combine functions with `group_by()`, which
splits the data frame, performs the function, and then combines the
output. Here, we show the above output grouped by eduction level.

``` r
fhs %>%
  group_by(EDUC) %>%
    summarise(TC_AVG = mean(TOTCHOL),
            TC_MIN = min(TOTCHOL),
            TC_MAX = max(TOTCHOL),
            TOTAL_N = n())
```

    ## # A tibble: 4 x 5
    ##    EDUC TC_AVG TC_MIN TC_MAX TOTAL_N
    ##   <int>  <dbl>  <int>  <int>   <int>
    ## 1     1   238.    113    432    1630
    ## 2     2   235.    119    600    1135
    ## 3     3   238.    140    464     636
    ## 4     4   236.    150    403     450

-----

#### Question 4:

  - What is the sex and smoking status of the youngest diabetic?
  - What is the average systolic blood pressure smokers who currently
    take blood pressure medication?
  - For diabetics who experienced MI, what is their mean age, BMI, and
    glucose level?
  - How many overweight smokers are in each education category?

-----

## See you tomorrow\!

Hopefully, this introduction to the basics of R was informative and
helpful. Tomorrow, we will start by diving into data visualisation and
performing some statistical tests.

-----
