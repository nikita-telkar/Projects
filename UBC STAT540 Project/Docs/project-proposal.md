# STAT540 Group Project Proposal  
## Authors: The GINNS  
## Due Date: 13 February 2020  

***

## 1. Motivation and Background Work (1 pt)

Systemic lupus erythematosus (SLE) [OMIM 152700](https://www.omim.org/entry/152700) is a complex autoimmune disease characterized by tissue deposition of autoantibodies and immune complexes, which leads to chronic inflammation of connective tissue [1,2]. With the lack of a cure, treatment is limited to symptom relief [3]. As is the case with most autoimmune diseases, prevalence of SLE is higher in women, and approximately 20-70 per 10000 people are affected globally [3].

SLE is considered a disease of multifactorial etiology, attributable to both genetic and environmental factors [3,4]. GWAS studies have identified over 60 risk loci for SLE susceptibility across populations to date, most of which are regulatory loci (e.g. HLA) [4]. Epigenetic mechanisms play an important role in gene expression regulation, and they are directly influenced by the environment. Given the multifactorial etiology of SLE, the integration of gene expression and DNA methylation (DNAm) patterns rather than their individual analysis must be considered to expand our knowledge of the mechanisms underlying this complex disease.

The **primary aim** of our project is to investigate the interaction between DNAm and gene expression in SLE patients using the matched DNAm and gene expression dataset [GSE82221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82221), and validate findings of the [original study](https://arthritis-research.biomedcentral.com/articles/10.1186/s13075-016-1050-x) with this same dataset [5]. While the literature is rich in studies that investigate SLE-relevant DNAm and gene expression patterns separately, it lacks a comprehensive analysis of the interaction between the two data types. The detail provided in the methods section of this publication are somewhat scarce in terms of the steps that were taken for DNAm and gene expression analysis; no information is provided regarding QC, and whether the X and Y chromosome probes were removed to eliminate the impact of sex is not clarified. More information regarding our intended approach can be found in the *Dataset* section of this proposal, but we will be using raw data and performing the QC, normalization and analysis from scratch to validate their findings and create a reproducible pipeline.

As briefly mentioned earlier, SLE is more frequent in females than males. However, the mechanisms underlying the genetic basis of this sex-specific prevalence that characterizes many autoimmune diseases remains largely unclear. Past studies have identified abnormal X chromosome inactivation in SLE [6], and susceptibility loci have also been identified on this chromosome [7]. Following this rationale, our **second aim** is to separately analyse interactions between DNAm and gene expression patterns in the X chromosome of female SLE patients and female controls, both of which constitute the majority of the samples.

***

## 2. Division of labour (1 pt)

> Who is going to do what? State assignment of tasks and projected contributions for each group members.  

* **Sierra & Naila** - Differential expression analysis  

* **Icíar & Nikita** - Methylation differential analysis

* **Everyone** - Integrative Analysis

Which genes are found as significant when: examining each data type alone? Modeling disease by both expression and methylation? How do these sets compare?   
Enrichment analysis of genes found in the combined model.

> Please provide a table of group members with their background, degree, affiliations and job assignments.

| Group Member     | Degree                      | Affiliation  |   Job Assignment                                 |
| ---------------- | --------------------------- | ------------ | ------------------------------------------------ |
| Icíar Fernández  | Genome Science & Technology | UBC & BCCHR  | Project background and motivation, DNAm analysis |
| Naila Adam       | Genome Science & Technology | UBC & BCCRC  | Aims and Methodology, RNA expression analysis|
| Nikita Telkar    | Medical Genetics            | BCCHR & BCCRC | Project Aims and Methodology, DNAm analysis |
| Sierra Gillis    | Genome Science & Technology | UBC & BCCRC  |Dataset, RNA expression analysis |

***

## 3. Dataset (1 pt) - Sierra

* What kind of data are you working with? What is the general description and characteristics of the data?

**Whole genome DNA methylation data and RNA expression data** was taken from Peripheral Blood Mononuclear Cells collected from patients via blood sample. The cohort includes 25 normal controls of individuals without SLE, 15 individuals with SLE that are negative for Lupus Nephritis (LN-), and 15 patients with SLE and Lupus Nephritis (LN+). The controls were matched with the SLE patients for age, sex, and ethnicity, and were found to be not significant when the researchers tested for difference between patients and controls in terms of these factors. Between SLE LN- and SLE LN+ groups, disease duration varies by a couple months, but this factor was also said to be not significant. One factor that was found to be significant was the disease activity (SLEDAI score) and level of Complement Component 3. We suppose this is due to the LN+ group being a progression of the SLE disease, so a difference in disease activity is expected between the LN- and LN+ groups. Analysis the blood tests further, the authors found that expected factors are in the normal range for SLE patients. Also, these patients did not have lymphopenia or cytopenia. It is not explicitly said that these can be confounders, but it is likely why researchers take these other diseases into account.

* State the technology used to generate the data.

To generate the DNA methylation data, the researchers used the **Illumina HumanMethylation 450k BeadChip array** to find the status of 485,000 CpG sites (99% of RefSeq genes and 96% of CpG islands). The beadchip was scanned and the tool GenomeStudio Methylation module software was used to measure the intensities. The score for methylation of each CpG is the fluorescent intensity ratio known as a beta value, where 0 means the site is non-methylated, and 1 means the site is completely methylated. The authors further extracted, filtered and normalized the values using the R package RnBeads, but the GEO dataset is labelled as RAW so we will be conducting our own normalization and filtering based on diagnostics.

RNA expression data, more specifically mRNA transcription levels, were collected using **Illumina HumanHT-12 v4.0 Expression Beadchips,** which can capture 47323 different transcripts. Again, the authors performed extraction and normalization using GenomeStudio Module, but we have access to the raw data. As this is expression data from a beadchip, it is also a ratio of fluorescent intensity instead of the count data obtained from RNA sequencing in order to measure expression.

The authors also collected cytokine protein levels from urine for many proteins. To not add too many variables into our model and risk over-explanation, our group has decided to focus on DNA methylation and RNA expression data.

Using GEOquery to explore the data structures:    
GSE82221 refers to the entire experiment. GSE82218 is a Biobase ExpressionSet object for the methylation data. GSE81622 is the RNA expression data. This can be seen by the GPL section on GEO. It may be easier to look at the data by looking at each GSE separately.    
So for GSE81622, the expression data can be found as follows:

> gEset <- getGEO("GSE81622", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]       
exprs(gEset) #this will print it

The sample variables can be obtained by:

> pData(gEset) # not all of these columns will go into the MetaData

The experimental design is as follows (when shorten to relevant columns):   
Note: "##" is replaced with a sample specific number.    
55 x 5   

|                 | title                                                                                     | type  | description.1 | age:ch1 | Sex:ch1            |
|:-----------------|:-------------------------------------------------------------------------------------------|:-------|:---------------|:--------|:--------------------:|
| Row Names (GSM21598##) | 3 options: "normal control-##",  "SLE patient-##",  "SLE patient with lupus nephritis-##" | "RNA" | "replicate #" | "##y"   | "Female" or "Male" |

Expression Values:   
47323 x 55

|           | GSM21598##        |
|:-----------|:-------------------:|
| ProbeName | Real number value |

Methylation Data:   
Note: The 'metadata' is different for this data for the columns of interest to us, so we will extract the needed info and adjust the dataframe to be the same for each dataset where possible.    
gEset <- getGEO("GSE82218", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]    
485577 x 55

|           | GSM21598##        |
|:-----------|:-------------------:|
| ProbeName? cg######## | 0 < value < 1 |

**Warning**: This matrix contains NA values

***

## 4. Aims and Methodology (2 pt) - Naila & Nikita

> What are the specific questions/aims that you will address in order to achieve the final goals of the project?

1. Are there any genes that are differentially expressed within SLE cases (LN- and LN+), and with controls?

2. Is there any correlation between the gene expression and the methylation status (DNAm) within cases and with controls?

3. What genes are differentially expressed between female cases and controls, and is there any interaction between their gene expression and DNAm status?

> Computational & statistical approaches:

  * We will explore normalization methods for microarray data and understand the Quality Control metrics based on which the data will be filtered for downstream analysis. for RNA expression data, the authors performed Quantile normalization using GenomeStudio. They also filtered the data further by detection p value and centered median using Genespring software. We will aim at building our workflow entirely in R and avoid using such external software for this project.

 * As for DNA methylation, `RnBeads` is a package that was cited for extracting, filtering and normalizing the data in the publication. This is a start-to-finish pipeline for methylation analysis. Instead, we will be following a framework that involves separate R packages for each step.

* To explore the data, PCA will be used for visualization and batch effect diagnosis. Moreover, hierarchical clustering of the samples will be performed.

* For comparative analysis, one-way ANOVA/two-sample *t* test and ANCOVAR for multivariable adjustments were used in the publication. We will further investigate differential expression using  `limma`. Additionally, we will investigate multiple testing correction methods which are used interchangeably in the paper (e.g. FDR/ Bonferroni correction).

* To infer biological pathways activity for the differentially expressed genes, IPA( Ingenuity Pathway Analysis) software was used by the authors. We will alternatively use `fgsea` package from R to carry out Gene Set Enrichment Analysis (GSEA) using GO terms curated lists available from the Broad Institute GSEA database[8].

***

### **References:**

1. Online Mendelian Inheritance in Man, OMIM®. McKusick-Nathans Institute of Genetic Medicine, Johns Hopkins University (Baltimore, MD). World Wide Web URL: https://omim.org/
2. Systemic lupus erythematosus. Genetics Home Reference. https://ghr.nlm.nih.gov/condition/systemic-lupus-erythematosus. Published 2020.
3. Danchenko N, Satia J, Anthony M. Epidemiology of systemic lupus erythematosus: a comparison of worldwide disease burden. Lupus. 2006;15(5):308-318. doi:10.1191/0961203306lu2305xx
4. Teruel M, Alarcón-Riquelme M. The genetic basis of systemic lupus erythematosus: What are the risk factors and what have we learned. J Autoimmun. 2016;74:161-175. doi:10.1016/j.jaut.2016.08.001
5. Zhu, H., Mi, W., Luo, H. et al. Whole-genome transcription and DNA methylation analysis of peripheral blood mononuclear cells identified aberrant gene regulation pathways in systemic lupus erythematosus. Arthritis Res Ther 18, 162 (2016). https://doi.org/10.1186/s13075-016-1050-x
6. Syrett CM, Paneru B, Sandoval-Heglund D, et al. Altered X-chromosome inactivation in T cells may promote sex-biased autoimmune diseases. JCI Insight. 2019;4(7):e126751. Published 2019 Apr 4. doi:10.1172/jci.insight.126751
7. Zhu, Z., Liang, Z., Liany, H. et al. Discovery of a novel genetic susceptibility locus on X chromosome for systemic lupus erythematosus. Arthritis Res Ther 17, 349 (2015). https://doi.org/10.1186/s13075-015-0857-1
8. Sergushichev A (2016). “An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation.” bioRxiv. doi: 10.1101/060012, http://biorxiv.org/content/early/2016/06/20/060012.
