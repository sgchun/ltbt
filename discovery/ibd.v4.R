## For questions, please contact:
## Sung Chun (SungGook.Chun@childrens.harvard.edu)
## Shamil Sunyaev (ssunyaev@rics.bwh.harvard.edu)

## UKB phenotype file has to be parsed separately in the following steps:
## 1. Parse UKB fields 40001-* and 41202-*.
## 2. Prepare icd10_ibd.tsv as an example below:
##     eid	41202-0.1	41202-0.2	41202-0.1	...
##     100001     L409	L408
##     100002                                     E063
##     ...
## 

# ICD10 codes of major autoimmune diseases EXCLUDING IBD
autoimmune.icd10 <-
    c(
        ## "K50", "K500", "K501", "K508", "K509", 
        ## "K51", "K510", "K511", "K512", "K513", "K514", "K515", "K518",
        ## "K519",
        "D510", "D591", "D693", "E050", "E063", "E271", "G610",
        "K743", "K900", 
        "L809", "M313", "M315", "M316", "M353",
        "G700", "M321", "M329", "M350", "M459",
        "E10", "E100", "E101", "E102", "E103", "E104", "E105", "E106", "E107", 
        "E108", "E109",
        "G35",
        "H20", "H200", "H201", "H202", "H208", "H209",
        "K73", "K730", "K731", "K732", "K738", "K739",
        "L10", "L100", "L102", "L108", "L109",
        "L12", "L120", "L121", "L128", "L129",
        "L40", "L400", "L401", "L403", "L404", "L405", "L408", "L409",
        "L63", "L630", "L638", "L639",
        "M05", "M050", "M051", "M052", "M053", "M058", "M059",
        "M06", "M060", "M061", "M062", "M063", "M064", "M068", "M069",
        "M08", "M080", "M081", "M082", "M083", "M084", "M089",
        "M33", "M331", "M332", "M339",
        "M34", "M340", "M341", "M342", "M348", "M349")

icd10 <- 
    read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.icd10.tsv", header=TRUE, stringsAsFactors=FALSE,
               sep="\t")
dim(icd10)

icd10 <- icd10[, c(1, grep("f.40001.", colnames(icd10), fixed=TRUE),
                   grep("f.41202.", colnames(icd10), fixed=TRUE)) ]
dim(icd10)


## The ICD10 codes for IBD are K50 & K51
has.IBD <-
    apply(as.matrix(icd10[, 2:ncol(icd10)]), 1,
          function (v) (length(union(grep("K50", v, fixed=TRUE),
                                     grep("K51", v, fixed=TRUE))) > 0))

cases.eid <- icd10$eid[has.IBD]

## We exclude individuals with any common autoimmune disorders from controls
has.other.autoimmune <-
    apply(as.matrix(icd10[, 2:ncol(icd10)]), 1,
          function (v) any(v %in% autoimmune.icd10))

other.autoimmune.eid <- icd10$eid[has.other.autoimmune]

## UKB phenotype file has to be parsed separately in the following steps:
## 1. Parse UKB fields 20002-* (self-reported noncancer illness codes)
## 2. Extract all individuals with any values
## 3. Prepare self_reported_non_cancer.tsv as below:
##     eid	20002-0.0	20002-0.1	20002-0.2	20002-0.3	20002-0.4	20002-0.5	20002-0.6	20002-0.7	20002-0.8	20002-0.9	20002-0.10	20002-0.11	20002-0.12	20002-0.13	20002-0.14	20002-0.15	20002-0.16	20002-0.17	20002-0.18	20002-0.19	20002-0.20	20002-0.21	20002-0.22	20002-0.23	20002-0.24	20002-0.25	20002-0.26	20002-0.27	20002-0.28	20002-1.0	20002-1.1	20002-1.2	...
##     1000001	1485	1473	1353
##     1000002	1065	1111
##     1000003
##     1000004	1456
##     1000005	1113
##     1000006	1065	1286	1485
##     ...
## If you know which columns to extract, you can prepare the files by:
##   sed -E 's/"([^"]*)",/\1\t/g' ukbXXXXX.csv | cut -f <your_column_positions>

self.report <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/processed_phenotypes/ukb47653.self_reported_illness.tsv",
                          header=TRUE, sep="\t", stringsAsFactors=FALSE)
dim(self.report)


## self-reported illness codes of common autoimmune diseases, 
## excluding IBD, Crohns, and UC: 
# Pernicious anemia: 1331
# Thyrotoxicosis: 1225
# type 1 diabetes: 1222
# adrenocortical insufficiency: 1234
# multiple sclerosis: 1261
# Guillain-Barre 1256
# primary biliary cirrhosis: 1506
# celiac: 1456
# Pemphigus/Pemphigoid: 1345
# Psoriasis: 1453
# Alopecia: 1667
# Vitiligo: 1661
# rheumatoid arthritis: 1464
# wegners granulmatosis: 1378
# Dermatopolymyositis: 1383
# polymyalgia rheumatica: 1377
# myasthenia gravis: 1437
# systemic sclerosis: 1384
# systemic lupus erythematosis/sle: 1381
# sjogren's syndrome: 1382
# ankylosing spondylitis: 1313
autoimmune.codes <- c(1331, 1225, 1222, 1234, 1261, 1256, 1506, 1456, 1345,
                      1453, 1667, 1661, 1464, 1378, 1383, 1377, 1437,
                      1384, 1381, 1382, 1313)

## self-reported illness codes of IBD:
# IBD 1461
# Crohns 1462
# UC 1463
self.report.ibd <-
    apply(as.matrix(self.report[, 2:ncol(self.report)]), 1,
          function (v) any(v[!is.na(v)] %in% 1461:1463))

## Individuals with self-reported autoimmune diseases are excluded from
## controls
self.report.other.autoimmune <-
    apply(as.matrix(self.report[, 2:ncol(self.report)]), 1,
          function (v) any(v[!is.na(v)] %in% autoimmune.codes))

self.report.ibd.eid <-
    unique(self.report[self.report.ibd, 1])

self.report.other.autoimmune.eid <-
    unique(self.report[self.report.other.autoimmune, 1])

length(self.report.ibd.eid)
length(self.report.other.autoimmune.eid)

## Combine ICD10 diagnosis codes and self-reported illness
cases.eid <- union(cases.eid, self.report.ibd.eid)

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

## sex
female <-
    qc$eid[!is.na(qc$genetic.sex) & (qc$genetic.sex == qc$sex) &
                      (qc$genetic.sex == 0)]
male <-
    qc$eid[!is.na(qc$genetic.sex) & (qc$genetic.sex == qc$sex) &
                      (qc$genetic.sex == 1)]

length(female)
length(male)

## Filter cases and controls
cases.eid <- intersect(cases.eid, qc.ok.eid)
length(cases.eid)

## Individuals with any autoimmune diseases are excluded from controls
controls.eid <-
    setdiff(qc.ok.eid, union(union(cases.eid, other.autoimmune.eid), 
                             self.report.other.autoimmune.eid))

length(controls.eid)

######################################################################
## Cross-check
ibd1 <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/proj/nps/chun2020/ukb31063/ibdtr.phen.updated.phen",
                  header=T, sep="\t")

ibd2 <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/proj/nps/chun2020/ukb31063/ibdval_v3.phen.updated.phen",
                  header=T, sep="\t")

length(cases.eid)
sum(ibd1$Outcome == 1) + sum(ibd2$Outcome == 1)
length(intersect(cases.eid, union(ibd1$IID[ibd1$Outcome == 1],
                                  ibd2$IID[ibd2$Outcome == 1])))

######################################################################

## Kinship filtering
## UKB provides kinship inference table
rel <- read.delim("/lab-share/Pulmonary-Chun-e2/Public/data/ukb/ukb_rel_a62142_s488221.dat", 
       header=TRUE, sep=" ", stringsAsFactors=FALSE)
dim(rel)

## Make sure we do not have related individuals among cases
## should be 0
rel.cases <- rel[rel$ID1 %in% cases.eid & rel$ID2 %in% cases.eid, ]
nrow(rel.cases)

filter.ID1 <- runif(nrow(rel.cases)) < 0.5

sum(filter.ID1)

cases.eid <- setdiff(cases.eid, union(rel.cases$ID1[filter.ID1], 
                                      rel.cases$ID2[!filter.ID1]))
length(cases.eid)

stopifnot(sum(rel$ID1 %in% cases.eid & rel$ID2 %in% cases.eid) == 0)


## Exclude related individuals between cases and controls
## Keep cases and drop controls.
rel.ctrls1 <- rel[rel$ID2 %in% cases.eid & rel$ID1 %in% controls.eid, ]
rel.ctrls2 <- rel[rel$ID1 %in% cases.eid & rel$ID2 %in% controls.eid, ]

dim(rel.ctrls1)
dim(rel.ctrls2)

controls.eid <- setdiff(controls.eid, union(rel.ctrls1$ID1, rel.ctrls2$ID2))

stopifnot(sum(rel$ID2 %in% cases.eid & rel$ID1 %in% controls.eid) == 0)
stopifnot(sum(rel$ID1 %in% cases.eid & rel$ID2 %in% controls.eid) == 0)

## Exclude related controls
## Randomly remove one 
rel.ctrls <- rel[rel$ID1 %in% controls.eid & rel$ID2 %in% controls.eid, ]
nrow(rel.ctrls)

filter.ID1 <- runif(nrow(rel.ctrls)) < 0.5

sum(filter.ID1)

controls.eid <- setdiff(controls.eid, union(rel.ctrls$ID1[filter.ID1], 
					    rel.ctrls$ID2[!filter.ID1]))
length(controls.eid)

stopifnot(sum(rel$ID1 %in% controls.eid & rel$ID2 %in% controls.eid) == 0)

## Save sample ids
write.table(data.frame(eid=cases.eid),
            file="ibd_cases.eid",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(data.frame(eid=controls.eid),
            file="ibd_controls.eid",
            quote=FALSE, row.names=FALSE, col.names=FALSE)


