---
title: "Progress Report"
author: "The GINNS"
date: "3/13/2020"
output: 
  html_document: 
    keep_md: yes
---

## What has changed based on the final proposal? (2 pt.)

+ **Did your dataset change? If so, why?** Our dataset hasn't changed, we are still working on the analysis of [GSE82218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82218).


+ **Have you decided to do a different analysis than what was mentioned in your proposal? If so, Why?** So far, our initial proposal for an analysis remains the same. If anything, it is possible that we will not have enough time to pursue our idea of analysing the X chromosome separately, but this was more of a side-goal from the beginning. Elaborating some more on this aspect that was mentioned in our project proposal, the idea would be to look at X-linked methylation and gene expression of female controls vs. disease females.

+ **Are there any changes in task assignments of group members?** No, Sierra & Naila are still working on the gene expression aspect, whereas Iciar & Nikita are working on methylation, and we will all come together for the integrative analysis aspect.

## What is the progress of the analyses? (5 pts.)

+ Briefly and concisely explain your **methodology and progress** for the aims you have investigated so far. Which parts were modified and which parts remained the same?

  + **Methylation analysis**: Thus far, we have completed QC and normalization, and ran some diagnostics / exploratory analysis of our data. Detailed explanations for all of our steps are provided in each of the R markdown files relevant to our analysis, found in the [02_methylation](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Code/02_methylation) folder within the [Code](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Code) folder of this repo. Briefly, our QC and filtering steps for the dataset include removal of XY probes, cross-hybridizing probes, SNPs, and polymorphic probes. Normalization presented somewhat of a challenge due to the limitation of not having IDATs, which greatly reduces the preprocessing methods that can be used. In addition, the lack of technical replicates is also a hindrance in evaluating the performance of different normalization methods. We ultimately chose BMIQ, which is widely recommended for all datasets and yielded good results. In terms of EDA, we have initially found sex to have a greater effect size than disease state. While this is certainly surprising, we have a few ideas as to what may be going on. Flaws in experiment design, such as running male and female samples in separate arrays, is a possibility. We cannot know for sure because there is no Sentrix ID nor plate or array information in the metadata. If due to true biological differences, we could think that given the large sex bias in autoimmune conditions, it may actually be representative of our dataset. However, it is likely a result of both. We need to further explore these results with PCA. In addition, we found a couple samples that appear to be sex-mismatched during our EDA, and we are considering whether we should re-annotate them or drop them altogether seeing as we cannot be certain that their associated metadata is correct either. 

  + **Expression Analysis**: Diagnostic checks were done on the raw data after filtering using the pvalues associated with signal (as given by the authors) and by probe quality for this Illumina platform. This exploratory analysis can be found [here](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/exp_proc_analysis.md). Based on the results of that and the authors' analysis, the expression data was Quantile normalized using [this script](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/Quantile-Normalization.md) and this normalized dataset is [here](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Data/gene-expression_GSE81622/Qunatile.normalized.txt). The diagnostic and quality control checks were [redone](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/diagnostics.md) on this dataset. Next, [differential expression](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/multivar_DE.md) using multiple covariates was performed.
  
  
  
  
+ What **R packages or other tools** are you using for your analyses? You do not need to provide your scripts in your report.

  + **Methylation analysis**: We are taking a cross-package approach to the analysis of our data, and we have so far used: minfi, missMethyl, wateRmelon, and methylumi in our progress. For the differential methylation analysis, which we will be completing this coming week, we will be using limma and DMRcate. 
  
  + **Expression analysis**: Throughout all analyses, `dplyr` was used to manipulate dataframes. Data exploration involved using base r for histograms and `ggplot2` for other graphs except in the case where pheatmap was used to do clustering and generate heatmaps. The gene names, locations, and probe qualities were taken from the Illumina database in the package `illuminaHumanv4.db`. The library `genefilter` was used to filter genes. The multivariate linear regression was performed using limma. For [univariate]("https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/Univariate_DE.md") analysis, we aimed at reporducing the results in the paper using the same methodology. We  incorporated `findMarkers` function from the `scran` library to preform t-tests.
We tested `fgsea` package as well for Gene Set Enrichment Analysis (GSEA). We downloaded [GO Gene sets]("https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp") from the Broad Institute MSigDB Collections(C5: GO gene sets). 

+ Provide the links to any markdown reports within your repo to refer to the relevant analysis. Provide references.

  + All the **methylation analysis** scripts (available in .Rmd, .md and html) can be found in the [methylation 02](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Code/02_methylation) folder. Within each report, we have referenced the package vignettes, any source code attributed to other authors, and relevant publications that we have used as a reference for our chosen methods.
  
  + The **expression analysis** scripts can be found in the [01_gene-expression](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Code/01_gene-expression) and are linked more specifically in the first bullet above.

## Results (3 pts.) 

+ What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results. List some challenges that you encountered? How will you address them?

  + **Methylation analysis**: Unfortunately, we haven't reached the stage of exploring our hypothesis; it has been quite time-consuming to preprocess the data, and to fully understand the steps that we are taking. However, we plan to complete the differential methylation analysis this week now that our data is ready, and to start the integrative analysis. Linear modelling will be performed using limma to find differentially expressed methylation probes between cases and controls, within cases, and within females, and interaction of covariates will also be assessed. We also intend on using the DMRcate package for our analysis.
  


Below is a summary of some results obtained so far.

![figure-1](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/01qc%2Bnormalization_methylation_files/figure-html/densityplots.jpg)

  We compared several normalization methods and decided on the BMIQ normalization.
  
  Then, we compared the MDS plots of samples before and after QC+normalization. We found quite interesting - and unexpected - that samples clustered by sex rather than by status (disease/normal) both before and after, as explained earlier in the report.
  
![figure-2](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/02exploratory_methylation_files/figure-html/mds_preQC_all.jpg)
  
![figure-3](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/02exploratory_methylation_files/figure-html/mds_postQC_all.jpg)
  
  + **Expression analysis**: Filtering by p-value associated with raw signal data and by Illumina probe quality gave data that clustered seemingly randomly with pheatmap. It also gave concerning results in preliminary t-tests in terms of p-value distributions violating the null hypothesis when testing DE between disease status SLE and SLE-LN. Thus, we decided to do further analyses with the quantile normalized data. This dataset is also what the authors used. Overall, PCA analysis didn't show clustering of the patients/samples based on disease status.
  The DE analysis with limma showed many genes differentially expressed between the normal and SLE individuals. Interestingly, one of the top differentially expressed probes is not associated with a gene according to Illumina. The SLE to SLE-LN comparison is as expected based on the data diagnostics and does not have differentially expressed genes. Further analyses will be done to compare Males to Females.
  In terms of univariate analysis, after preforming t-test for the first comparison( normal vs SLE), we selected top significant genes ( p<0.05) and ranked them based on logFC to run GSEA. We didn't obtain a significantly enriched pathway( padj<0.05) in this comparison which is alarming. We are in the process of diagnosing our code and reconciling those results.
  See this [summary](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/analysis_summary.md) for more details of the analyses.


