library(data.table)

########################################################################
## gnomAD reference exome data (lifted over to GRCh38)
########################################################################

gnomad <- NULL

for (I in 1:22) {
    file <- paste("~/data/gnomad/2.1.1.liftover/",
                  "gnomad.exomes.r2.1.1.sites.", I,
                  ".liftover_grch38.probdam.popmax.tsv.gz", sep='')
    gnomad0 <- read.delim(file ,sep=" ", stringsAsFactors=FALSE, header=FALSE)

    gnomad <- rbind(gnomad, gnomad0)
}

colnames(gnomad)[1:6] <- c("CHRBP", "SNPID", "REF", "ALT", "FILTER", "POPMAX")

key <- paste(gsub(":", "_", gnomad$CHRBP, fixed=TRUE),
             gnomad$REF, gnomad$ALT,
             sep="_")

gnomad <- cbind(gnomad, key, stringsAsFactors=FALSE)
gnomad$POPMAX <- as.numeric(gnomad$POPMAX)

gnomad <- gnomad[!is.na(gnomad$POPMAX), ]

########################################################################
## Exome data from UKB (controls)
########################################################################

## Read HCLOF variants from UKB exome VCF into a data frame as below: 
# 
# locus.contig locus.position    alleles            gene           AF
#       <char>          <int>     <char>          <char>        <num>
#        chr10      100042575 ['T', 'G'] ENSG00000120054 1.064212e-06
# gene_symbol  qual           het homalt 
#      <char> <num>        <char> <char> 
#        CPN1    38 100000,123456   <NA> 
# ... 
ukb <- READ.HCLOF.FROM.UKB.EXOME.VCF()

## Parse ukb$alleles into A1 and A2
alleles <- gsub(" ", "", gsub("'", "", ukb$alleles))
alleles <- sapply(alleles, function (alstr) {
       strsplit(substr(alstr, 2, nchar(alstr) - 1), ",")
   })
n.alleles <- sapply(alleles, function (al) {
       length(al)
   }, simplify=TRUE)
stopifnot(all(n.alleles == 2))
A1 <- sapply(alleles, function (al) {
       al[1]
   }, simplify=TRUE)
A2 <- sapply(alleles, function (al) {
       al[2]
   }, simplify=TRUE)

## Create a unique variant identifier for each row 
keys <- paste(ukb$locus.contig, ukb$locus.position, A1, A2, sep="_")
ukb <- cbind(ukb, key=keys, stringsAsFactors=FALSE)
ukb <- as.data.table(ukb)
setkey(ukb, key)

## Read QC-failed variants from UKB exome VCF
## See variant_qc.py
ukb.qcfail.key <- READ.UKB.QCFAIL.VARIDS()

########################################################################
## Phenotype data from UKB (controls)
########################################################################

## Read death registry + summary diagnosis as below: 
#
#       eid f.40001.0.1 f.40001.0.2 ...
#    100001        L409        L408
icd10 <- READ.UKB.ICD10(c("40001", "40002", "41270"))

## IPF: J84.1
has.ipf.icd10 <-
    apply(as.matrix(icd10[, 2:ncol(icd10)]), 1,
          function (v) (length(grep("J841", v, fixed=TRUE)) > 0))
                               
cases.eid <- icd10$eid[has.ipf.icd10]

## UKB QC fields:
## See https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100313
## Sample QC files needs to be prepared separately.
##
## eid - subject ID
## sex - from 31
## genetic sex - from 22001
## het.missing.outliers - from Field 22027
## putative.sex.chromosome.aneuploidy - from Field 22019
## excluded.from.kinship.inference - from Field 22021
## in.white.British.ancestry.subset - from Field 21000
qc <- READ.UKB.QC(c("22019", "22027", "31", "22001", "22021", "21000"))

qc.ok.eid <- qc$eid[is.na(qc$sex.chromosome.aneuploidy) & 
                    is.na(qc$outliers.for.heterozygosity.or.missing.rate) &
                    (qc$sex == qc$genetic.sex) &
                    # not excluded from kinship inference
                    (qc$genetic.kinship.to.other.participants != -1) &
                    # White British
                    !is.na(qc$genetic.ethnic.grouping) &
                    (qc$genetic.ethnic.grouping == 1)]

## Filter cases and controls
cases.eid <- intersect(cases.eid, qc.ok.eid)
controls.eid <- setdiff(qc.ok.eid, cases.eid)

## ---------------------------------------------------------------------
## Load UKB kinship table and remove related individuals
#
#    ID1    ID2
# 100001 100002
rel <- READ.UKB.KINSHIP.TABLE()

## Make sure we do not have related individuals by andomly removing one 
all.eid <- union(cases.eid, controls.eid)
rel.sel <- rel[rel$ID1 %in% all.eid & rel$ID2 %in% all.eid, ]

filter.ID1 <- runif(nrow(rel.sel)) < 0.5

cases.eid <- setdiff(cases.eid, 
                     union(rel.sel$ID1[filter.ID1], rel.sel$ID2[!filter.ID1]))

controls.eid <- setdiff(controls.eid, 
                        union(rel.sel$ID1[filter.ID1], rel.sel$ID2[!filter.ID1]))

## ---------------------------------------------------------------------
## Read subject-level data
subj <- READ.UKB.SUBJECT.DATA()

## age >= 65 only 
## The phenotype data is from 2021, so we use 2021 as the reference year for age calculation.
age <- 2021 - subj$birthyear
subj$age <- age

subj <- subj[subj$age >= 65, ]
dim(subj)

subj <- subj[subj$eid %in% c(cases.eid, controls.eid), ]
subj <- subj[order(subj$eid), ]

## ---------------------------------------------------------------------
## Read 14-SNP IPF PRS from UKB cohort as:
# 
#    eid    PRS
# 100001 0.1234
#
prs <- READ.UKB.14SNP.IPF.PRS()

prs <- prs[prs$eid %in% subj$eid, ]
prs <- prs[order(prs$eid), ]

ukb.prs.mean <- mean(prs$PRS)
ukb.prs.mean
ukb.prs.sd <- sd(prs$PRS)
ukb.prs.sd

## standardize including both cases and controls
prs$PRS <- (prs$PRS - ukb.prs.mean) / ukb.prs.sd

## ---------------------------------------------------------------------
## Identify all subjects with Interstitial Lung Disease (ILD): J84

has.ild.icd10 <-
    apply(as.matrix(icd10[, 2:ncol(icd10)]), 1,
          function (v) (length(grep("J84", v, fixed=TRUE)) > 0))

icd10.ild.eid <- icd10$eid[has.ild.icd10]

## self-reported non-cancer illness table
#
#      eid f.20002.0.0 f.20002.0.1 f.20002.0.2	...
#   100001	      1485        1473       1353  ...
#
self.report <- READ.UKB.SELF.REPORTED.ILLNESS()

## self-reported ILD Codes: 
## abestosis: 1120, pulmonary fibrosis: 1121, fibrosing alveolitis: 1122
self.report.ila <-
    apply(as.matrix(self.report[, 2:ncol(self.report)]), 1,
          function (v) any(v[!is.na(v)] %in% c(1120:1122)))
self.report.ila.eid <- unique(self.report[self.report.ila, 1])
length(self.report.ila.eid)

ila.eid <- unique(c(icd10.ild.eid, self.report.ila.eid))

# Remove subjects with ILD from cases and controls
subj <- subj[!(subj$eid %in% ila.eid), ]
dim(subj)

prs <- prs[!(prs$eid %in% ila.eid), ]
dim(prs)

stopifnot(all(prs$eid == subj$eid))

########################################################################
## Exome data from CGS-PF (case)
########################################################################

## Read HCLOF variants from CGS-PF exome VCF into a data frame as below: 
# 
# locus.contig locus.position    alleles            gene           AF
#       <char>          <int>     <char>          <char>        <num>
#        chr10      100042575 ['T', 'G'] ENSG00000120054 1.064212e-06
# gene_symbol  qual           het homalt 
#      <char> <num>  <char> <char> 
#        CPN1    38 100,123   <NA> 
# ... 
ipf <- READ.HCLOF.FROM.CGSPF.EXOME.VCF()

alleles <- gsub(" ", "", gsub("'", "", ipf$alleles))
alleles <- sapply(alleles, function (alstr) {
       strsplit(substr(alstr, 2, nchar(alstr) - 1), ",")
   })
n.alleles <- sapply(alleles, function (al) {
       length(al)
   }, simplify=TRUE)
stopifnot(all(n.alleles == 2))
A1 <- sapply(alleles, function (al) {
       al[1]
   }, simplify=TRUE)
A2 <- sapply(alleles, function (al) {
       al[2]
   }, simplify=TRUE)

keys <- paste(ipf$locus.contig, ipf$locus.position, A1, A2, sep="_")
ipf <- cbind(ipf, key=keys, stringsAsFactors=FALSE)

## Read QC-failed variants from CGS-PF exome VCF
## See variant_qc.py
ipf.qcfail.key <- READ.CGSPF.QCFAIL.VARIDS()

## ---------------------------------------------------------------------
## Read 14-SNP IPF PRS from CGS-PF cohort as:
# 
#  id    PRS
# 100 0.1234
#
ipf.prs <- READ.CGSPF.14SNP.IPF.PRS()

## Standarize per UKB White population
ipf.prs$PRS <- (ipf.prs$PRS - ukb.prs.mean) / ukb.prs.sd

######################################################################
## Harmonization
########################################################################

## Exclude QC-failed sites from both CGS-PF and UKB datasets
ipf.exclude <- which(ipf$key %in% ukb.qcfail.key)
ipf <- ipf[ -ipf.exclude, ]

ukb.exclude <- which(ukb$key %in% ipf.qcfail.key)
ukb <- ukb[ -ukb.exclude, ]

# Since all IPF cases are unrelated, Allele Count (AC) over 2 is likely 
# due to sequencing errors, so we can filter out those variants. 
ipf <- ipf[ipf$AC > 0 & ipf$AC <= 2, ]

## ---------------------------------------------------------------------
## Apply the AF cutoff
af.co <- 0.001

# Calculate the combined AF across cases and controls for each variant, 
# and filter out variants with combined AF above the cutoff.
ipfaf <- sapply(which(ukb$key %in% ipf$key), function(I) {
    ipf$AF[ipf$key == ukb$key[I]]
  }, simplify=TRUE)

ukb <- cbind(ukb, ipfaf=0)
ukb$ipfaf[which(ukb$key %in% ipf$key)] <- ipfaf

cmbaf <- (ukb$AF * N.CTRL * 2 + ukb$ipfaf * N.CASE * 2) /
    (N.CTRL * 2 + N.CASE * 2)

ukb <- ukb[cmbaf < af.co, ]

ukbaf <- sapply(which(ipf$key %in% ukb$key), function(I) {
    ukb[ipf$key[I], AF]
  }, simplify=TRUE)

ipf <- cbind(ipf, ukbaf=0)
ipf$ukbaf[which(ipf$key %in% ukb$key)] <- ukbaf

cmbaf <- (ipf$AF * N.CASE * 2 + ipf$ukbaf * N.CTRL * 2) /
    (N.CTRL * 2 + N.CASE * 2)

ipf <- ipf[cmbaf < af.co, ]

## ---------------------------------------------------------------------
## Apply gnomad popmax AF cutoff
popmax.af <- 0.01

ipf <- ipf[!(ipf$key %in% gnomad$key[gnomad$POPMAX > popmax.af]), ]
ukb <- ukb[!(ukb$key %in% gnomad$key[gnomad$POPMAX > popmax.af]), ]
