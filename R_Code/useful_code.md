---
title: "Useful R Code"
author: Nikita Telkar | Rob Lab
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 4
    theme: journal
    keep_md: true
---


***
**Data Types:** https://swcarpentry.github.io/r-novice-inflammation/13-supp-data-structures/

*** 

#### **Installing Packages**
If `install.packages()` doesn't work:


```
library(devtools)
install_github("gagolews/stringi")

OR

install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))

```

*** 

#### **Colour Change + other formatting options**
``` 
<span style="color:green"> Hello World! </span>  
Or  
<span style="color:#A763FF"> Hello World! </span>
Or
<span style='font-weight:bold;font-family:"Palatino Linotype";'><span style='color:#689107;'>Hello World!</span></span>
```

***
#### **Convert df to tibble**
```
library(tibble)
as_tibble()
```

***
#### **Convert df to dropdown table**  
```
library(DT)
datatable()
```
***

#### **Transpose Data: Interchange rows and columns**  
```
df <- t(df)  #Imp: Changes data type
```
***


#### **Printing**
```
print(raw_miRNA[2:28, 1], n = 27)
```  
  Or  
```
#alt way to print --> convert into separate variable and print
var <- raw_miRNA[2:28, 1]       
print(as_tibble(var), n = 27)
```  
Or  
```
paged_table() #best table display when datatable() is too heavy to be processed
```


***

#### **Return unique elements**  
```
unique(t_raw3[,1]) #data type is character  
                      #include [] inside of ()
```
***

#### **Deleting/Selecting columns/rows**  
```
sub_ti <- t_mi [c(1:4)]   #select first 4 columns
sub_ti <- head(sub_ti, 10)   #select first 10 rows
sub_ti <- tail(sub_ti, -10)   #delete first 10 rows
sub_ti <- head(sub_ti, -10)   #delete last 10 rows

```
***

#### **Finding Values from one in another**  

```
meta_dat$caseID[meta_dat$caseID %in% samp$Case_ID]
``` 

Negation:
```
meta_dat$caseID[!meta_dat$caseID %in% samp$Case_ID]

```
***
