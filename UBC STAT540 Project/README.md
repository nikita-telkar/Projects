# Analysing DNA methylation & gene expression patterns in systemic lupus erythematosus (SLE)  
## A project by the GINNS for STAT540

***  

| Group Member     | Degree                      | Affiliation  |   Job Assignment|
| ---------------- | --------------------------- | ------------ | ----------------------------------------|
| Icíar Fernández  | Genome Science & Technology | UBC & BCCHR  | Project background and motivation, DNAm analysis |
| Naila Adam       | Genome Science & Technology | UBC & BCCRC  | Aims and Methodology, RNA expression analysis|
| Nikita Telkar    | Medical Genetics            | BCCHR & BCCRC | Project Aims and Methodology, DNAm analysis |
| Sierra Gillis    | Genome Science & Technology | UBC & BCCRC  |Dataset, RNA expression analysis |  

***
**Repo directory**

* The [**code**](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Code) folder is the location for all the scripts that we will be working on. Within the folder, scripts are divided by methylation and gene expression analysis, and there will eventually be a third folder for integrative analysis.

* The [**data**](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Data) folder contains the methylation and gene expression data & metadata. All data was downloaded from GEO, accession number [GSE82221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82221) for the SuperSeries, with [GSE81622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81622) containing gene expression, and [GSE82218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82218) containing methylation.

* The [**docs**](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Docs) folder holds other documents relevant to the project, such as our [project proposal](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Docs/project-proposal.md) and eventually, our progress report and other deliverables.

* The [**publication**](https://arthritis-research.biomedcentral.com/articles/10.1186/s13075-016-1050-x) where the data was initially published.

***

**Data**

+ [**gene-expression_GSE81622**](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Data/gene-expression_GSE81622) holds data files related to the gene expression dataset, with GEO Accession # [GSE81622](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81622).

+ [**methylation_GSE82218**](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/tree/master/Data/methylation_GSE82218) contains data files related to the methylation dataset, with GEO Accession # [GSE82218](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82218).

Both datasets are part of the [GSE82221 SuperSeries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82221). 

***  

**Analysis and results**

* Gene Expression - we performed [univariate](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/Univariate_DE.md) and [multivariate](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/multivar_DE.md) linear regression as well as [enrichment analysis](https://github.com/STAT540-UBC/Repo_team-The-GINNS_W2020/blob/master/Code/01_gene-expression/MultiVarModel_Enrichment.md) in order to find differentially expressed genes and their potential involvement in the disease.  
  + **Results** - Even though DE genes did correspond with those found in the publication, we found that no pathways were significantly downregulated, but there was upregulated in pathways for viral/bacterial infections, interferon type I, leukocyte, cytokines and granulocytes response.  
  
  
  
* Methylation Expression - We analysed methyldation data in order to reproduce the [differential methylation expression](https://github.com/nikita-telkar/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/03_differential_DNAm_analysis.md) results by way of linear modelling, and also [pathway analysis](https://github.com/nikita-telkar/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/04_DMPs_analysis.md). A separate insight into analysing only the [X Chromosome probes](https://github.com/nikita-telkar/Repo_team-The-GINNS_W2020/blob/master/Code/02_methylation/05xchr_analysis.Rmd) was also undertaken.
    + **Results** - Only 5 hypomethylated probes were reproducible by our analysis at the statistical threshold selection provided in the paper. None of the probes showed any significant pathway enrichment at p<0.5.  
