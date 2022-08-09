# easyVAF
Somatic sequence variants are associated with a cancer diagnosis, prognostic stratification, and treatment response. Variant allele frequency (VAF) is the percentage of sequence reads with a specific DNA variant over the read depth at that locus. VAFs on targeted loci under different (experimental) conditions are often compared. We present our R package  ‘ esayVAF’ for parametric and non-parametric comparison of VAFs among multiple treatment groups. 


To install and use the package, you may download the directory or follow the instructions below.

# Install package
if (!require("devtools")) install.packages("devtools")

devtools::install_github("zhuangyh/easyVAF")

# Load package
library(easyVAF)
In the vignettes folder, users can find a documentation that illustrates how to implement easyVAF with an example data set. The data file is included under the data folder. Details on all the functions included in the package are documented in the package manual under the package folder. Users may directly download the package tarball for all functions and example data, which is also under the package folder.

Please report any issues at https://github.com/zhaungyh/easyVAF/issues.

# References

> Junxiao Hu,Vida Alami, Yonghua Zhuang, Dexiang Gao.  “easyVAF, a R package for VAF comparison among groups”. Journal of Open Source Software, 2022. (*Submitted*)
