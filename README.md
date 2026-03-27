# Liability Threshold Burden Test (LTBT)
LTBT implements a novel rare-variant association test for binary traits that jointly models the contribution of rare variants and PRS under a liability threshold model. We define the total liability as the additive sum of rare-variant effect (γ), PRS, and unexplained risk factors (**Figure A**). Under a liability threshold model, the disease outcome is determined by whether an individual’s total liability is over the threshold determined by the prevalence of disease in population. By modeling the conditional distribution of PRS given the rare variant effect size γ (**Figure B**), we compare the likelihood of observed data under the null (γ=0) and alternative (γ≠0) hypotheses. 

![LTBT overview figure](LTBT-GitHub-Fig.png)

For inquiries on LTBT, please contact Sung Chun (sung.g.chun@gmail.com). 

## Installation
LTBT is implemented as an R module. Please download and install [the LTBT R module](https://github.com/sgchun/ltbt/blob/main/ltbt_1.0.tar.gz) as below:
   ```bash
   R CMD install ltbt_1.0.tar.gz
   ```
## How to Run LTBT
We provide a deidentified summary-level sample test data [RTEL1.tsv](https://github.com/sgchun/ltbt/blob/main/RTEL1.tsv) and test code [test.r](https://github.com/sgchun/ltbt/blob/main/test.r). 
   ```R
   library(ltbt)

   gamma.seq <- seq(-5, 5, 0.05)

   exdata <- read.delim("RTEL1.tsv", sep="\t")

   res <- run.ltbt(gamma.seq, exdata$gt, exdata$prs, exdata$outcome, 0.005)

   res$h2Lx
   # 0.1266585
   res$pvalue
   # 4.824949e-05
   res$gamma.mle
   # 1.45

   plot(res$gamma, res$logL, ty="l", xlab="gamma", ylab="Log Likelihood Ratio")
   abline(v=res$gamma.mle, col="red", lty=2)
   ```
![LTBT Log Likelihood Ratio Curve](RTEL1.llr.png)

## Citation
> Sung Chun*, Ahmad Samiei, Lauren Flynn, Heidi Makrynioti, Matthew Moll, Anna L. Peljto,
> David A. Schwartz, Michael H Cho, Shamil R Sunyaev, Ivan O Rosas, Gary M. Hunninghake,
> and Benjamin A Raby*. A new liability-based rare-variant burden test identifies a novel
> genetic association between HTRA3 and idiopathic pulmonary fibrosis. 
> (Preprint will be available soon)> * Joint contribution.

