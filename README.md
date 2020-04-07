# Cardiovascular Signature Temporal Clustering Platform

During cardiovascular disease progression, multiple cellular systems (e.g., cardiac proteome) undergo diverse molecular changes. The temporal patterns of individual biomolecules depict their unique responses towards pathological drivers and contribute to underlying biological processes disrupted by disease progression. Advances in high-throughput omics technology have enabled cost-effective temporal profiling of targeted systems in animal models of human disease.

```CV.Signature.TCP``` is a open source R package to identifying temporal molecular signatures from time-course omics data. It is consisted of 3 modules, “preprocessing”, “clustering”, and “evaluation”.  First, two independent methods cubic splines and reduced rank PCA are provided for data preprocessing (i.e., missing data imputation and denoising). Next, two independent classification methods, k-means or hierarchical clustering are available to classify temporal patterns of biological variables. Finally, a jackstraw test is employed for evaluation, identifying biological variables strongly related to the representative pattern of temporal clusters.

We found that our platform produced biological meaningful clusters, enabling further functional delineation. Its flexible parameter settings and analytical routes allow a broader adaptation to other time-course omics data.

# Installation

To install this package from GitHub, please use ```devtools``` and also set repositories to both CRAN and Bioconductor:

```R
install.packages("devtools")
library("devtools")

setRepositories(ind=1:2)
install_github("UCLA-BD2K/CV.Signature.TCP")
```

#### Bioconductor dependencies

Some of Bioconductor dependencies may fail to be automatically installed. This is a known issue at the moment:
https://github.com/r-lib/devtools/issues/700#issuecomment-235127291

Therefore, if you get an error message about "package or namespace load failed", install them manually:

```R
source("https://bioconductor.org/biocLite.R")
biocLite(c('multtest'))
```

# Data

This package includes a partial dataset from the following publication:

J Wang, H Choi, NC Chung, Q Cao, DCM Ng, B Mirza, SB Scruggs, D Wang, AO Garlid, P Ping (2018). Integrated dissection of the cysteine oxidative post-translational modification proteome during cardiac hypertrophy. Journal of Proteome Research <https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.8b00372>