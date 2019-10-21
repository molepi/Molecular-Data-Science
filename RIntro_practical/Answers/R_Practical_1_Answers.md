-----

## Part 1: RStudio and R Markdown

-----

#### Answers to Question 1:

  - `na.rm` indicates whether missing (NA) values should be stripped
    before the `mean` is calculated
  - `sqrt` accepts numeric or complex vectors and arrays

-----

## Part 2: Vectors

Define numeric variable `x`

``` r
x <- 1 + 1
x
```

    ## [1] 2

Load `tidyverse`

``` r
library(tidyverse)
```

Define numeric vector `y`

``` r
y <- c(5, 6, 7)
y
```

    ## [1] 5 6 7

Define character vector `z`

``` r
z <- c("a", "b", "c")
z
```

    ## [1] "a" "b" "c"

Define logical vector `w`

``` r
w <- c(T, T, F)
w
```

    ## [1]  TRUE  TRUE FALSE

-----

#### Answers to Question 2:

``` r
y[3]
```

    ## [1] 7

``` r
z[c(1,3)]
```

    ## [1] "a" "c"

``` r
y + 0:x
```

    ## [1] 5 7 9

``` r
rev(z)
```

    ## [1] "c" "b" "a"

-----

## Part 3: Data Frames

Combine `y` and `z` into data frame `df`

``` r
df <- data.frame(col1 = y, col2 = z, row.names = c("row1", "row2", "row3"))
df
```

    ##      col1 col2
    ## row1    5    a
    ## row2    6    b
    ## row3    7    c

Add `col3` to `df`

``` r
df$col3 <- w
df
```

    ##      col1 col2  col3
    ## row1    5    a  TRUE
    ## row2    6    b  TRUE
    ## row3    7    c FALSE

Add `row4` to `df` to make `df3`

``` r
df2 <- data.frame(col1 = 8, col2 = "d", col3 = F, row.names = "row4")
df3 <- rbind(df, df2)
df3
```

    ##      col1 col2  col3
    ## row1    5    a  TRUE
    ## row2    6    b  TRUE
    ## row3    7    c FALSE
    ## row4    8    d FALSE

#### Answers to Question 3:

Add a new column with the name `col4` to data frame `df` with the values
`1, 1, 1`

``` r
df$col4 <- c(1,1,1)
df
```

    ##      col1 col2  col3 col4
    ## row1    5    a  TRUE    1
    ## row2    6    b  TRUE    1
    ## row3    7    c FALSE    1

Make a new data frame `df4`, which is a subset of data frame `df` but
only contains columns `col1` and `col2`, and rows `row2` and `row3`

``` r
df4 <- df[c("row2", "row3"), c("col1", "col2")]
df4
```

    ##      col1 col2
    ## row2    6    b
    ## row3    7    c

Calculate the sum of `col1` of data frame `df3`

``` r
sum(df3["col1"])
```

    ## [1] 26

Change the row names of data frame `df4` to `rowX` and `rowY`

``` r
rownames(df4) <- c("rowX", "rowY")
df4
```

    ##      col1 col2
    ## rowX    6    b
    ## rowY    7    c

-----

## Part 4: Real Data

Read in the FHS data

``` r
fhs <- read.csv(url("https://raw.githubusercontent.com/molepi/Molecular-Data-Science/master/RIntro_practical/data.csv"))
```

Remove observations with missing values

``` r
fhs <- fhs[complete.cases(fhs), ]
sum(is.na(fhs))
```

    ## [1] 0

Define formula `f`

``` r
f <- y ~ x + z
f
```

    ## y ~ x + z

#### Question 4:

Make a table of current smoking status

``` r
xtabs(~ CURSMOKE, fhs)
```

    ## CURSMOKE
    ##   No  Yes 
    ## 1962 1889

Make a proportion table of BP medication

``` r
prop.table(xtabs(~ BPMEDS, fhs))
```

    ## BPMEDS
    ##         No        Yes 
    ## 0.96598286 0.03401714

Make a table showing the number of smokers for each education level

``` r
xtabs(~ CURSMOKE + EDUC, fhs)
```

    ##         EDUC
    ## CURSMOKE   1   2   3   4
    ##      No  871 529 340 222
    ##      Yes 759 606 296 228

-----

## Part 5: Data Manipulation

Create a subset of FHS containing BP variables

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

Create a subset of FHS containing all variables but MI

``` r
fhsCOV <- select(fhs, -MI)
```

Create a subset of FHS containing BP variables using the `contains()`
function

``` r
fhsBP <- select(fhs, contains("BP"))
```

Show only observations where `BMI` is greater or equal to 40

``` r
fhsBMI40 <- filter(fhs, BMI >= 40)
```

Show only observations for female smokers

``` r
fhsSmokeF <- filter(fhs, SEX == 'Female', CURSMOKE == 'Yes')
```

-----

#### Answers to Question 5:

What is the sex and smoking status of the youngest diabetic? **Female
non-smoker**

``` r
fhs %>%
  select(SEX, CURSMOKE, AGE, DIABETES) %>%
  filter(DIABETES == "Yes") %>%
  arrange(AGE) %>%
  head(n = 1)
```

    ##      SEX CURSMOKE AGE DIABETES
    ## 1 Female       No  36      Yes

What is the average systolic blood pressure of smokers who currently
take BP medication? **168 mmHg**

``` r
fhs %>%
  filter(CURSMOKE == "Yes", BPMEDS == "Yes") %>%
  summarise(AVG_SYSBP = mean(SYSBP))
```

    ##   AVG_SYSBP
    ## 1  168.3444

For diabetics who experienced MI, what is their mean age, BMI, and
glucose level? **55.4 years, 27.9 kg/m^2, and 165 mg/dL**

``` r
fhs %>%
  filter(DIABETES == "Yes", MI == "Yes") %>%
  summarise(AVG_AGE = mean(AGE),
            AVG_BMI = mean(BMI),
            AVG_GLU = mean(GLUCOSE))
```

    ##    AVG_AGE  AVG_BMI  AVG_GLU
    ## 1 55.44444 27.92067 165.8444

How many overweight smokers are in each education category? Do you
notice a trend? **Yes, there are less smokers in higher education
groups**

``` r
fhs %>%
  filter(BMI >= 30, CURSMOKE == "Yes") %>%
  group_by(EDUC) %>%
  summarise(N_EDUC = n())
```

    ## # A tibble: 4 x 2
    ##    EDUC N_EDUC
    ##   <int>  <int>
    ## 1     1     93
    ## 2     2     37
    ## 3     3     24
    ## 4     4     15

Now look at the total number of people in each education category
alongside the percentage who are overweight smokers. Is there a trend
now? **No, when looking at the percentage there doesn’t appear to be a
trend**

``` r
fhs %>%
  mutate(OVRSM = ifelse(BMI >= 30 & CURSMOKE == "Yes", 1, 0)) %>%
  group_by(EDUC) %>%
  summarise(EDUC_N = n(),
            OVRSM_N = sum(OVRSM),
            PERC = OVRSM_N * 100 / EDUC_N)
```

    ## # A tibble: 4 x 4
    ##    EDUC EDUC_N OVRSM_N  PERC
    ##   <int>  <int>   <dbl> <dbl>
    ## 1     1   1630      93  5.71
    ## 2     2   1135      37  3.26
    ## 3     3    636      24  3.77
    ## 4     4    450      15  3.33

-----
