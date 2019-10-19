We continue from yesterday by exploring the relationship between
variables, both with visualisations and statistical tests. Don’t forget
to turn your R scripts into the Blackboard.

-----

## Part 1: Introduction to ggplot2

First, we must load in the dataset again, keeping only complete cases.

``` r
fhs <- read.csv(url("https://raw.githubusercontent.com/molepi/Molecular-Data-Science/master/RIntro_practical/data.csv"))
fhs <- fhs[complete.cases(fhs), ]
```

We should also set our library location and load in the `tidyverse`
package.

``` r
.libPaths("C:/fos_2019/library")
```

``` r
library(tidyverse)
```

-----

The package `ggplot2` is commonly used to produce visualisations of
data. The plots will appear in the **Plots** pane, you can save them
using the Export button if you would like.

The `ggplot()` function alone just initializes a blank plot. You can
specify your x and y axis using `aes()`, but no lines or points will be
drawn. You can pipe your data to `ggplot()` using `%>%` or specify it as
the first argument.

``` r
fhs %>%
  ggplot(aes(BMI, GLUCOSE))
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

-----

Titles can be specified by adding `ggtitle()` to the `ggplot()` function
with the `+` operator. X and Y axis labels can also be implemented with
`xlab()` and `ylab()`, respectively.

``` r
fhs %>%
  ggplot(aes(BMI, GLUCOSE)) +
  ggtitle("Graph of BMI against Glucose levels") +
  xlab("BMI (kg/m^2)") +
  ylab("Glucose level (mg/dL")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

-----

In order to add data visualisation to our graph, we need to add `geom`
elements. One of these, `geom_bar()` is used to create bar charts for
categorical variables.

The `aes()` function can be moved inside the `geom_bar()`, which is
useful if you want to overlap multiple `geom` elements showcasing
different aspects of your data.

Additionally, `alpha`, `colour`, and `fill` options can specify the
appearance of the bars. Note, you can spell colour either the American
or English way - both are accepted as arguments. Changing these options
can be done for other `geom` elements as well.

``` r
fhs %>%
  ggplot() +
  geom_bar(aes(EDUC), alpha=0.7, colour='grey30', fill='cadetblue2') +
  ggtitle("Bar plot of Education levels") + 
  xlab("Education level") 
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

-----

Bars can be coloured by other variables to see the distribution of one
category within another. This is done by moving the `fill` option inside
the `aes()` function and specifying a categorical variable.

``` r
fhs %>%
  ggplot() +
  geom_bar(aes(EDUC, fill=SEX), alpha=0.7, color='grey30') +
  ggtitle("Bar plot of sex within Education levels") + 
  xlab("Education level") 
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

-----

Subsetting, as shown yesterday with `dplyr` can also be piped into
`ggplot()`. We can show the above graph only for those who experienced
MI, for example.

``` r
fhs %>%
  filter(MI == 'Yes') %>%
  ggplot() +
  geom_bar(aes(EDUC, fill=SEX), alpha=0.7, color='grey30') +
  ggtitle("Bar plot of sex within Education levels") + 
  xlab("Education level") 
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

-----

To view bars side-by-side instead of stacked, set `position` to `dodge`.

``` r
fhs %>%
  filter(MI == 'Yes') %>%
  ggplot() +
  geom_bar(aes(EDUC, fill=SEX), alpha=0.7, color='grey30', position="dodge") +
  ggtitle("Bar plot of sex within Education levels") + 
  xlab("Education level") 
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

-----

For continuous variables, we divide the x axis into bins and visualize
using `geom_histogram()`.

``` r
fhs %>%
  ggplot() +
  geom_histogram(aes(AGE), binwidth=2, alpha=0.7, colour='grey30', fill='cadetblue3') +
  ggtitle("Histogram of age distribution") +
  xlab("Age (years)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

-----

To show the distribution of a continuous variable within categories, you
can use `facet_wrap()` which takes a `formula` class as an argument.
Separating plots by group using `facet_wrap()` is possible for most
`ggplot` graphics.

``` r
fhs %>%
  ggplot() +
  geom_histogram(aes(AGE), binwidth=2, alpha=0.7, colour='grey30', fill='cadetblue3') +
  facet_wrap(~EDUC) +
  ggtitle("Histogram of age distribution within education levels") +
  xlab("Age (years)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

-----

In order to show a smoothed version of a histogram, a density plot, we
use `geom_density`.

``` r
fhs %>%
  ggplot() +
  geom_density(aes(SYSBP), alpha=0.7, colour='grey30', fill='cadetblue2') +
  ggtitle("Density plot of systolic blood pressure") +
  xlab("Systolic BP (mmHg)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

-----

If you want to look at a variable’s distribution between groups, use the
`fill` option within `aes()`.

``` r
fhs %>%
  ggplot() +
  geom_density(aes(SYSBP, fill=BPMEDS), alpha=0.7, colour='grey30') +
  ggtitle("Density plot of systolic blood pressure within BP medication status") +
  xlab("Systolic BP (mmHg)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

-----

This will show the distributions separately, but you can also stack them
with a `position` argument.

``` r
fhs %>%
  ggplot() +
  geom_density(aes(SYSBP, fill=BPMEDS), alpha=0.7, colour='grey30', position="stack") +
  ggtitle("Density plot of systolic blood pressure within BP medication status") +
  xlab("Systolic BP (mmHg)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

-----

#### Question 1: Make the following:

  - A bar plot of MI, coloured by smoking status
  - A histogram of total cholesterol in overweight individuals
  - Four density plots of age for each education level
  - Two bar plots of sex in current and non-smokers for those who
    experienced MI

-----

## Part 2: Comparative Plots

A scatter plot can be a good way to compare two continuous variables. In
`ggplot` this is implemented using `geom_point()`.

``` r
fhs %>%
  ggplot() +
  geom_point(aes(SYSBP, DIABP), alpha=0.7) +
  ggtitle("Scatter plot of Systolic BP against Diastolic BP") +
  xlab("Systolic BP (mmHg)") +
  ylab("Diastolic BP (mmHg)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

-----

Reference lines can be added to plots using `geom_vline()`,
`geom_hline()`, and `geom_abline()`. These are very useful for
annotating plots.

``` r
fhs %>%
  ggplot() +
  geom_point(aes(SYSBP, DIABP), alpha=0.7) +
  geom_vline(aes(xintercept=160), color='red', size=1.5, alpha=0.7) +
  ggtitle("Scatter plot of Systolic BP against Diastolic BP") +
  xlab("Systolic BP (mmHg)") +
  ylab("Diastolic BP (mmHg)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

-----

You can also use the `cut()` function alongside `seq()` to create
categorical variables, to help with visualisations.

``` r
fhs %>%
  ggplot() +
  geom_point(aes(SYSBP, DIABP), alpha=0.7) +
  geom_vline(aes(xintercept=160), color='red', size=1.5, alpha=0.7) +
  ggtitle("Scatter plot of Systolic BP against Diastolic BP") +
  xlab("Systolic BP (mmHg)") +
  ylab("Diastolic BP (mmHg)") +
  facet_wrap(. ~ cut(AGE,seq(30,70,10)))
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

-----

To show the distribution of a continuous variable across categories, you
can use `geom_boxplot()`. This visualizes the median, upper and lower
quartiles, and two whiskers, which extend to 1.5 IQRs. All outlying
points are shown individually.

``` r
fhs %>%
  ggplot() +
  geom_boxplot(aes(CURSMOKE, BMI), color='grey30', fill='cadetblue2') +
  ggtitle("Distribution of BMI within smoking status") +
  xlab("Current smoking status") +
  ylab("BMI (kg/m^2)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

-----

If you want to view a horizontal boxplot, you can use `coord_fip()`.
Outliers can also be highlighted by changing their colour and shape.

``` r
fhs %>%
  ggplot() +
  geom_boxplot(aes(CURSMOKE, BMI), color='grey30', fill='cadetblue2', outlier.color='palegreen3', outlier.shape=1) + 
  coord_flip() +
  ggtitle("Distribution of BMI within smoking status") +
  xlab("Current smoking status") +
  ylab("BMI (kg/m^2)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

-----

Finally, box plots can also be split by a further category using the
`fill` or `color` options within `aes()`.

``` r
fhs %>%
  ggplot() +
  geom_boxplot(aes(CURSMOKE, BMI, fill=SEX), color='grey30', outlier.color='palegreen3', outlier.shape=1) + 
  coord_flip() +
  ggtitle("Distribution of BMI within smoking status") +
  xlab("Current smoking status") +
  ylab("BMI (kg/m^2)")
```

![](R_Intro_Practical_2_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

-----

#### Question 2: Make the following:

  - A scatter plot of glucose levels against BMI, coloured by diabetes
    status.
  - A horizontal box plot of systolic BP by sex for people over the age
    of 60
  - A scatter plot of total cholesterol against glucose for categories
    of BMI

-----

## Part 3: Statistical Tests

Next we will perform some statistical tests. Sometimes, it can be
helpful to save tables so that they can be used in other functions.

``` r
SMOKE_MI <- xtabs(~CURSMOKE+MI, fhs)
SMOKE_MI
```

    ##         MI
    ## CURSMOKE   No  Yes
    ##      No  1671  291
    ##      Yes 1537  352

-----

Then, we can use this table in a chi-square test, to see if there is a
significant association between smoking status and MI.

``` r
chisq <- chisq.test(SMOKE_MI)
chisq
```

    ## 
    ##  Pearson's Chi-squared test with Yates' continuity correction
    ## 
    ## data:  SMOKE_MI
    ## X-squared = 9.7325, df = 1, p-value = 0.00181

In this instance, we have sufficient evidence (X2=9.73, p=0.002) that
smoking and MI are significantly associated. Pearson residuals can be
inspected to explore which groups are contributing to this evidence.

``` r
round(chisq$residuals, 3)
```

    ##         MI
    ## CURSMOKE     No    Yes
    ##      No   0.905 -2.022
    ##      Yes -0.923  2.061

-----

Visualizations can often inform future statistical investigations.
Perhaps, the scatter plot above suggested a positive correlation between
systolic and diastolic blood pressure, but we want to investigate this
further.

We can use `cor.test` to test if the correlation between two variables
is significant.

``` r
cor <- cor.test(~ SYSBP + DIABP, fhs) 
cor
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  SYSBP and DIABP
    ## t = 78.845, df = 3849, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.7734942 0.7976652
    ## sample estimates:
    ##       cor 
    ## 0.7858797

This shows that we have significant evidence from our data that systolic
and diastolic blood pressure are strongly correlated (r=0.79, p\<0.001).

-----

The box plot shown above might suggest that older individuals are more
likely to have experienced MI. We can formally test this using a t-test.

``` r
tt <- t.test(AGE ~ MI, fhs) 
tt
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  AGE by MI
    ## t = -9.0315, df = 918.42, p-value < 2.2e-16
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -4.072404 -2.618479
    ## sample estimates:
    ##  mean in group No mean in group Yes 
    ##          49.34352          52.68896

This test concludes that we have sufficient evidence to suggest a
significant difference in age between those who had and did not have MI
(t=-9.03, p\<0.001).

-----

To run logistic regression in R, you use the `glm()` function,
specifying the `family` as binomial. So explore the fit, use the
`summary()` function.

``` r
glmfit <- glm(MI ~ BMI + AGE + SEX, fhs, family = "binomial") 
summary(glmfit)
```

    ## 
    ## Call:
    ## glm(formula = MI ~ BMI + AGE + SEX, family = "binomial", data = fhs)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.2636  -0.6548  -0.4721  -0.3237   2.4734  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept) -6.14563    0.39969 -15.376  < 2e-16 ***
    ## BMI          0.05965    0.01089   5.476 4.36e-08 ***
    ## AGE          0.04537    0.00522   8.691  < 2e-16 ***
    ## SEXMale      1.19593    0.09413  12.705  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 3474  on 3850  degrees of freedom
    ## Residual deviance: 3191  on 3847  degrees of freedom
    ## AIC: 3199
    ## 
    ## Number of Fisher Scoring iterations: 5

This returns the estimate, standard errors, z-score, and p-values for
each of the coefficients. Here, we have strong evidence (p\<0.001) that
higher BMI is associated with an increased risk of MI, after adjusting
for age and sex.

To calculate the OR, we take the exponent of the coefficients.

``` r
exp(coef(glmfit))
```

    ## (Intercept)         BMI         AGE     SEXMale 
    ## 0.002142822 1.061466118 1.046414460 3.306636063

Here, we can see that the odds of MI for males is 3.3 times higher than
for females, after adjusting for BMI and age.

-----

#### Question 3:

  - Is age correlated with total cholesterol?
  - Is MI associated with glucose levels?
  - Is MI associated with glucose levels in individuals with diabetes?

-----

## Part 4: Loops

One of the advantages of using R is to use loops to perform multiple
analyses at once. As an example, the following loop calculates the
association of MI with all the other variables and stores them in a data
frame.

``` r
results <- data.frame()       
for (i in colnames(fhs)[-ncol(fhs)]) {  
  fit <- glm(MI ~ get(i), fhs, family = "binomial") 
  results <- rbind(results, summary(fit)$coefficients[2, , drop = FALSE]) 
}
rownames(results) <- colnames(fhs)[-ncol(fhs)] 
results$OR <- exp(results$Estimate) 
results[order(results$"Pr(>|z|)"), ] 
```

    ##              Estimate   Std. Error   z value     Pr(>|z|)        OR
    ## SEX       1.153381422 0.0921407197 12.517608 5.981047e-36 3.1688902
    ## SYSBP     0.018909221 0.0017797581 10.624602 2.289869e-26 1.0190891
    ## AGE       0.044548141 0.0050507477  8.820108 1.143491e-18 1.0455553
    ## TOTCHOL   0.008106879 0.0009512993  8.521902 1.569545e-17 1.0081398
    ## DIABP     0.028722602 0.0034015327  8.444018 3.066132e-17 1.0291391
    ## DIABETES  1.307439829 0.1995849643  6.550793 5.723227e-11 3.6966974
    ## BMI       0.065226891 0.0099964058  6.525034 6.798595e-11 1.0674012
    ## GLUCOSE   0.008529646 0.0014836425  5.749125 8.970668e-09 1.0085661
    ## BPMEDS    0.782504278 0.1961813118  3.988679 6.644226e-05 2.1869421
    ## CURSMOKE  0.273897693 0.0867548074  3.157147 1.593209e-03 1.3150803
    ## EDUC     -0.059809137 0.0427093732 -1.400375 1.614011e-01 0.9419443

This is quite a complicated loop, but we can dissect the different parts
to get an idea of what loops in R consist of.

-----

An `if` statement is only run if the condition is true.

``` r
if (1 < 2) { 
  print("1 is smaller than 2") 
}
```

    ## [1] "1 is smaller than 2"

-----

An `else` statement is only run if the condition of the if-statement is
false.

``` r
if (1 > 2) {
  print("1 is greater than 2")
} else {
  print("1 is not greater than 2")
}
```

    ## [1] "1 is not greater than 2"

-----

With `for` loops you assign a range of values to a variable and perform
the actions in the loop for all these values in order.

``` r
for (i in 1:10) { # loop over values 1 to 10
  if (i < 5) {
    print(paste(i, "is smaller than 5"))
  } else {
    print(paste(i, "is not smaller than 5"))
  }
}
```

    ## [1] "1 is smaller than 5"
    ## [1] "2 is smaller than 5"
    ## [1] "3 is smaller than 5"
    ## [1] "4 is smaller than 5"
    ## [1] "5 is not smaller than 5"
    ## [1] "6 is not smaller than 5"
    ## [1] "7 is not smaller than 5"
    ## [1] "8 is not smaller than 5"
    ## [1] "9 is not smaller than 5"
    ## [1] "10 is not smaller than 5"

-----

You can also utilise loops and make your own functions in R.

``` r
greater <- function(x, y) { 
  if (x > y) {
    paste(x, "is greater than", y)
  } else {
    paste(x, "is not greater than", y)
  }
}
greater(1, 2)
```

    ## [1] "1 is not greater than 2"

``` r
greater(2, 1)
```

    ## [1] "2 is greater than 1"

-----

``` r
greater <- function(x, y, z) { 
  paste(x*y*z)
}
greater(1, 2, 3)
```

    ## [1] "6"

``` r
greater(2, 1, 4)
```

    ## [1] "8"

#### Question 6:

  - Make a function that requires as arguments x, y and z and returns
    their product.
  - Make a loop that returns the association of all columns of data
    frame fhs with age.

-----

You have just written your own simple function, but hopefully you can
see how complicated they could become. You can do so many things in R,
and new functions are being written all the time, but if a function to
do something doesn’t exist, you can always try to create your own.

See you after lunch\!

-----
