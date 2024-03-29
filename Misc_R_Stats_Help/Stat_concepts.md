---
title: "Statistics Concepts" 
author: Nikita Telkar 
output: 
  html_document:
    keep_md: true
theme: cosmo 
---   

***  


```
## Warning: package 'rmarkdown' was built under R version 3.6.3
```

```
## Warning: package 'knitr' was built under R version 3.6.3
```


**Poisson Distribution**  

+ Where the *mean* and *variance* of the data are equal (λ).  
+ Lower bound at 0 but no upper bound (discrete).  

Drawback: May not accurately describe the variability of the counts.  

![](https://i2.wp.com/www.theanalysisfactor.com/wp-content/uploads/2015/03/poisson.gif?resize=642%2C481)<!-- -->


***

**Normal Distribution** - [Article Link](https://www.theanalysisfactor.com/differences-between-normal-and-poisson-distributions/#:~:text=One%20difference%20is%20that%20in,t%20quite%20behave%20the%20same.)  

+ No bounds.  
+ Any value from -∞ to ∞ is possible  
+ Mean is not equal to variance.  

***

**Negative Binomial Distribution** -  [Article Link](https://towardsdatascience.com/negative-binomial-regression-f99031bb25b4)   

+ Regression model taht does not make the equi-dispersion assumption i.e.not assume that variance = mean.  
+ Variance greater than the mean = over-dispersion, variance less than the mean = under-dispersion.  
+ Extension of a Poisson distribution to account for *overdispersed* data.









