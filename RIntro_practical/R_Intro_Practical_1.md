These introductory practicals are designed to teach you the basics of R.
At the end, we expect you to be able to follow along with the other
practicals of this course. Please make sure that you turn your answers
in to the blackboard at the end.

-----

## Part 1: RStudio and R Markdown

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

Using the `up arrow` and `down arrow` in the console, you can scroll
through historical commands. These can also be viewed in the **History**
tab of the *Environment* section (see Figure).

-----

You can also type your code in the *R-script* (see Figure) and run a
line or a selection using `ctrl + enter`. This is preferable to typing
in the *Console* directly, since you can save the scripts for later use
or share them with collaborators.

When installing packages for this FOS course, we created a new library
for R.

Set your library path to this destination by typing the following code
and running it, and check that it worked. You will need to do this at
the start of all FOS practicals in R.

``` r
.libPaths("C:/fos_2019/library")
.libPaths()
```

-----

We advise you to use `R Markdown` documents during this course. These
will not only let you save your script for later use, but also easily
`knit` the results to a HTML and share those.

Open an `R Markdown` document by clicking the new document icon in the
top left, and selecting *R Markdown…* from the drop-down menu. Give your
document a name and keep HTML as the default output.

-----

You will now see several different chunks, similar to this:

![Figure
1](https://github.com/molepi/Molecular-Data-Science/blob/master/RIntro_practical/markdown.PNG)

-----

Save your `R Markdown` document to a directory that you would like to
work in, using the **Save** button (red).

The **Knit** button (blue) at the top can be used to render a HTML
document of your work. Try it now.

At the top of your document, between the `---` markers (yellow) is some
information about your document. This includes the title you chose and
the output type.

Any plain text (orange) you write here can be used to describe or
comment your code. You can specify different levels of header using `#`,
bold with `**`, and italic with `*`.

**Code chunks** (green) are initialized using the Insert drop down menu
(top right) and selecting R.

If you want to run a single code chunk, you can do this using the **Run
Current Chunk** button (pink).

-----

Delete all the chunks **with the exception of the header** (yellow).

Write a brief description of what this document will contain along with
a title. If you want to check how it will look, you can `Knit` your
document at any time.

`R Markdown` is a great way to ensure reproducibility in research, and
using it is a good habit to pick up early.

Copy the following equations in your `R Markdown` document in a code
chunk, place your cursor on the first line, and run the equations
line-by-line using `ctrl + enter`. You can also select multiple lines
and run them together, or run the entire chunk.

    5 - 2         # subtraction
    2 * 2         # multiplication
    6 / 2         # division
    3 ^ 2         # exponentiation
    sqrt(9)       # square root
    abs(-5)       # absolute value

-----

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
be installed together. It contains packages like `ggplot2`, mentioned in
the lecture, for creating graphics, and `dplyr`, which we will use today
for data manipulation.

A list of currently installed packages is shown in the **Packages** tab.
The functions contained within these can be loaded into R with the
`library()` function.

-----

Now, let’s load `tidyverse`:

``` r
library(tidyverse)
```

We will use this package later in the practical. After loading a package
into R, the functions inside become available for use.

Lastly, the **Help** tab can be used to view documentation on any loaded
function. `Tab` can be used as autocomplete in R, but can also be used
to browse options for completion. Find help on the function, `mean`.

    ?mean

-----

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

For character class objects, `nchar()` can be used to count the number
of letters.

``` r
nchar("apple")
```

    ## [1] 5

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

#### Question 2: Obtain the following vectors using the variables you have created. There are multiple ways to do this, but you only need to give one of them in your answers.

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

Create a data frame, using the `data.frame()` function.

``` r
df <- data.frame(y, z)  
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

You can add row names and column names to a data frame using
`rownames()` and `colnames()`.

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
df$col3 <- w
df
```

    ##      col1 col2  col3
    ## row1    5    a  TRUE
    ## row2    6    b  TRUE
    ## row3    7    c FALSE

Adding extra rows requires the `rbind()` function, and for values of the
new row to be of the correct class.

``` r
df2 <- data.frame(col1 = 8, col2 = "d", col3 = F, row.names="row4") 
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
  - Make a new data frame `df4`, which is a subset of data frame `df`
    but only contains columns `col1` and `col2`, and rows `row2` and
    `row3`.
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

First, we want to explore the data.

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

``` r
str(fhs)
```

    ## 'data.frame':    4434 obs. of  12 variables:
    ##  $ SEX     : Factor w/ 2 levels "Female","Male": 2 1 2 1 1 1 1 1 2 2 ...
    ##  $ TOTCHOL : int  195 250 245 225 285 228 205 313 260 225 ...
    ##  $ AGE     : int  39 46 48 61 46 43 63 45 52 43 ...
    ##  $ SYSBP   : num  106 121 128 150 130 ...
    ##  $ DIABP   : num  70 81 80 95 84 110 71 71 89 107 ...
    ##  $ CURSMOKE: Factor w/ 2 levels "No","Yes": 1 1 2 2 2 1 1 2 1 2 ...
    ##  $ BMI     : num  27 28.7 25.3 28.6 23.1 ...
    ##  $ DIABETES: Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ BPMEDS  : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ GLUCOSE : int  77 76 70 103 85 99 85 78 79 88 ...
    ##  $ EDUC    : int  4 2 1 3 3 2 1 2 1 1 ...
    ##  $ MI      : Factor w/ 2 levels "No","Yes": 2 1 1 1 1 2 1 1 1 1 ...

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

-----

## Part 5: Data Manipulation

As mentioned above, `dplyr` is a `tidyverse` package that is frequently
used for data summary and transformation. `dplyr` code is often more
intuitive to read and use than base R packages for data manipulation.

This package contains a number of functions, which we will introduce one
by one and then combine.

-----

The `select()` function allows you to subset specific columns out of a
dataset, similar to subsetting introduced above.

``` r
fhsBP <- select(fhs, SYSBP, DIABP, BPMEDS)
head(fhsBP)
```

    ##   SYSBP DIABP BPMEDS
    ## 1 106.0    70     No
    ## 2 121.0    81     No
    ## 3 127.5    80     No
    ## 4 150.0    95     No
    ## 5 130.0    84     No
    ## 6 180.0   110     No

You can also perform negative indexing, using the `-` operator.

``` r
fhsCOV <- select(fhs, -MI)
head(fhsCOV)
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1   Male     195  39 106.0    70       No 26.97       No     No      77
    ## 2 Female     250  46 121.0    81       No 28.73       No     No      76
    ## 3   Male     245  48 127.5    80      Yes 25.34       No     No      70
    ## 4 Female     225  61 150.0    95      Yes 28.58       No     No     103
    ## 5 Female     285  46 130.0    84      Yes 23.10       No     No      85
    ## 6 Female     228  43 180.0   110       No 30.30       No     No      99
    ##   EDUC
    ## 1    4
    ## 2    2
    ## 3    1
    ## 4    3
    ## 5    3
    ## 6    2

You can select a range of columns using the `:` operator, or use the
`contains()` function to find all columns containing a specific word or
combination of letters.

``` r
fhsBP <- select(fhs, contains("BP"))
head(fhsBP)
```

    ##   SYSBP DIABP BPMEDS
    ## 1 106.0    70     No
    ## 2 121.0    81     No
    ## 3 127.5    80     No
    ## 4 150.0    95     No
    ## 5 130.0    84     No
    ## 6 180.0   110     No

This can come in very useful when subsetting data frames. Other useful
functions include `starts_with()` and `ends_with()`, for column name
prefixes and suffixes respectively.

-----

In a similar manner to subsetting columns, `filter()` can be used to
filter observations or rows.

``` r
fhsBMI40 <- filter(fhs, BMI >= 40)
head(fhsBMI40)
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1   Male     178  52   160  98.0       No 40.11      Yes     No     225
    ## 2 Female     183  45   151 101.0       No 45.80       No     No      63
    ## 3 Female     278  66   187  88.0       No 40.52       No     No      84
    ## 4 Female     266  62   173  89.0       No 42.00       No     No      75
    ## 5 Female     180  60   200 122.5      Yes 44.27      Yes     No     150
    ## 6 Female     251  67   192 102.0       No 44.09       No    Yes      62
    ##   EDUC  MI
    ## 1    1  No
    ## 2    1  No
    ## 3    1  No
    ## 4    1 Yes
    ## 5    1  No
    ## 6    1  No

Here, we can quickly see the individuals who are morbidly obese in our
data set. A `,` can be used as a logical AND between filters. Perhaps we
are interested in girls who are current smokers.

``` r
fhsSmokeF <- filter(fhs, SEX == 'Female', CURSMOKE == 'Yes')
head(fhsSmokeF)
```

    ##      SEX TOTCHOL AGE SYSBP DIABP CURSMOKE   BMI DIABETES BPMEDS GLUCOSE
    ## 1 Female     225  61   150  95.0      Yes 28.58       No     No     103
    ## 2 Female     285  46   130  84.0      Yes 23.10       No     No      85
    ## 3 Female     313  45   100  71.0      Yes 21.68       No     No      78
    ## 4 Female     221  38   140  90.0      Yes 21.35       No     No      70
    ## 5 Female     291  46   112  78.0      Yes 23.38       No     No      89
    ## 6 Female     195  38   122  84.5      Yes 23.24       No     No      78
    ##   EDUC MI
    ## 1    3 No
    ## 2    3 No
    ## 3    2 No
    ## 4    2 No
    ## 5    2 No
    ## 6    2 No

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

Combining `mutate()` with `ifelse()` can help create helpful binary
variables.

The syntax of the `ifelse()` statement is `ifelse` (`condition`, `value
if TRUE`, `value if FALSE`). So, below we create a new variable `HYPO`
that takes the value `1` if `GLUCOSE` is 70 or less, and `0` otherwise.

If you want to use multiple conditions within `ifelse()`, you can use
the `&` and `|` operators.

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

#### Question 5:

  - What is the sex and smoking status of the youngest diabetic?
  - What is the average systolic blood pressure of smokers who currently
    take blood pressure medication?
  - For diabetics who experienced MI, what is their mean age, BMI, and
    glucose level?
  - How many obese smokers are in each education category? Do you notice
    a trend?
  - Now look at the total number of people in each education category
    alongside the percentage who are obese smokers. Is there a trend
    now?

-----

## See you tomorrow\!

Hopefully, this introduction to the basics of R was informative and
helpful. Tomorrow, we will start by diving into data visualisation and
performing some statistical tests.

-----
