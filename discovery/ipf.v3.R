## For questions, please contact:
## Sung Chun (SungGook.Chun@childrens.harvard.edu)

## Idiopathic Pulmonary Fibrosis

icd10 <- 
    read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.icd10.tsv", header=TRUE, stringsAsFactors=FALSE,
               sep="\t")
dim(icd10)

## death registry + summary diagnosis
icd10 <- icd10[, c(1, grep("f.40001.", colnames(icd10), fixed=TRUE),
                   grep("f.40002.", colnames(icd10), fixed=TRUE),
                   grep("f.41270.", colnames(icd10), fixed=TRUE)
                   ) ]
dim(icd10)


## Has ICD10 codes for IPF
## Changed J841 from J84
has.ipf.icd10 <-
    apply(as.matrix(icd10[, 2:ncol(icd10)]), 1,
          function (v) (length(grep("J841", v, fixed=TRUE)) > 0))
                               
cases.eid <- icd10$eid[has.ipf.icd10]

length(cases.eid)

## Controls
icd10 <- 
    read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.icd10.tsv", header=TRUE, stringsAsFactors=FALSE,
               sep="\t")
dim(icd10)

controls.eid <- icd10$eid
length(controls.eid)

controls.eid <- setdiff(controls.eid, cases.eid)
length(controls.eid)

all.eid <- union(cases.eid, controls.eid)
mean(all.eid %in% cases.eid)

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
## excess.relatives - from Field 22021
## in.white.British.ancestry.subset - from Field 21000

## Field 21000 (REDUNDANT)
## ancestry <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.ancestry.tsv",
##                        header=TRUE, sep="\t", stringsAsFactors=FALSE)
## ancestry <-
##     ancestry[, c(1, grep("f.21000.", colnames(ancestry), fixed=TRUE))]

## dim(ancestry)

qc <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.qc.tsv",
                 header=TRUE, sep="\t", stringsAsFactors=FALSE)
dim(qc)

qc.ok.eid <- qc$eid[is.na(qc$sex.chromosome.aneuploidy) & 
                    is.na(qc$outliers.for.heterozygosity.or.missing.rate) &
                    (qc$sex == qc$genetic.sex) &
                    # not excluded from kinship inference
                    (qc$genetic.kinship.to.other.participants != -1) &
                    # White British
                    !is.na(qc$genetic.ethnic.grouping) &
                    (qc$genetic.ethnic.grouping == 1)]
length(qc.ok.eid)

## White/British ethnicity (REDUNDANT)
## is.british <-
##     apply(as.matrix(ancestry[, 2:ncol(ancestry)]), 1,
##           function (v) any(v[!is.na(v)] == 1001))

## british.eid <- ancestry$eid[is.british]
## length(british.eid)

## qc.ok.eid <- intersect(british.eid, qc.ok.eid)
## length(qc.ok.eid)

## Filter cases and controls
cases.eid <- intersect(cases.eid, qc.ok.eid)
length(cases.eid)

controls.eid <- intersect(qc.ok.eid, controls.eid)
length(controls.eid)

all.eid <- union(cases.eid, controls.eid)
mean(all.eid %in% cases.eid)

######################################################################

## Kinship filtering
## UKB provides kinship inference table
rel <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/ukb_rel_a62142_s488221.dat", 
       header=TRUE, sep=" ", stringsAsFactors=FALSE)
dim(rel)

## Make sure we do not have related individuals
## Randomly remove one 
rel.sel <- rel[rel$ID1 %in% all.eid & rel$ID2 %in% all.eid, ]
nrow(rel.sel)

filter.ID1 <- runif(nrow(rel.sel)) < 0.5

sum(filter.ID1)

cases.eid <- setdiff(cases.eid, union(rel.sel$ID1[filter.ID1], 
                                      rel.sel$ID2[!filter.ID1]))
length(cases.eid)

controls.eid <- setdiff(controls.eid, union(rel.sel$ID1[filter.ID1], 
					    rel.sel$ID2[!filter.ID1]))
length(controls.eid)

all.eid <- union(cases.eid, controls.eid)
stopifnot(sum(rel$ID1 %in% all.eid & rel$ID2 %in% all.eid) == 0)

mean(all.eid %in% cases.eid)

## Save sample ids
write.table(data.frame(eid=cases.eid),
            file="ipf_cases.v3.eid",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(eid=controls.eid),
            file="ipf_controls.v3.eid",
            quote=FALSE, row.names=FALSE, col.names=FALSE)


length(cases.eid)
length(controls.eid)
