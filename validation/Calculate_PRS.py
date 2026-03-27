import os
import hail as hl

hl.init(default_reference='GRCh38')

# See the accompanied Allen2019.13SNPs.b38.qctool.txt
# rs2077551 was excluded because it is not present in All of Us genotype data.
modelfile = MODEL_FILE_PATH

scores = hl.import_table(modelfile, delimiter=' ',
                        types={'additive_beta': hl.tfloat32,'position': hl.tint})

scores = scores.drop('SNPID', 'rsid', 'heterozygote_beta', 'risk_score_identifier')
scores = scores.annotate(rsid=scores.chromosome + ":" + hl.str(scores.position) + ":" + scores.alleleA + ":" + scores.alleleB )
scores = scores.annotate(chr='chr' + scores.chromosome)
scores = scores.rename({'additive_beta':'score'})

scores = scores.annotate(chrpos=hl.locus(scores.chr, scores.position))
scores = scores.annotate(alleleAB = [scores.alleleA, scores.alleleB])
scores = scores.key_by("chrpos", "alleleAB" )
scores.describe()

gcBgenDir = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/bgen/"

## Calculate PRS for each chromosome and export results to tsv files.
## Excluding chr17 for missing rs2077551
for CHR in ["3", "4", "5", "6", "7", "11", "13", "15", "19", "20"]:
    mt = hl.import_bgen(gcBgenDir + "chr" + CHR + ".bgen",
                        entry_fields = ['GT'],
                        sample_file = gcBgenDir + "chr" + CHR + ".sample",
                        variants=scores.filter(scores.chr == 'chr' + CHR))
    mt = mt.filter_rows(hl.is_indel(mt.alleles[0], mt.alleles[1]), keep=False)
    mt.describe()
    print(mt.count())
    mt = hl.variant_qc(mt)
    mt.show(10)

    scores = scores.key_by("chrpos")

    mt = mt.annotate_rows(**scores[mt.locus])
    mt = mt.filter_rows(((mt.alleleB == mt.alleles[0]) & (mt.alleleA == mt.alleles[1])) | 
                        ((mt.alleleA == mt.alleles[0]) & (mt.alleleB == mt.alleles[1])), keep=True)
    flip = hl.case().when(mt.alleleB == mt.alleles[0], True).when(mt.alleleB == mt.alleles[1], False).or_missing()
    mt = mt.annotate_rows(flip=flip)
    mt.show()
    mt.flip.show()

    mt = mt.annotate_rows(prior=2 * hl.if_else(mt.flip, mt.variant_qc.AF[0], mt.variant_qc.AF[1]))
    mt = mt.annotate_cols(prs=hl.agg.sum( mt.score * hl.coalesce( hl.if_else(mt.flip, 2 - mt.GT.n_alt_alleles(), mt.GT.n_alt_alleles()), mt.prior)))

    mt.cols().export(OUTPUT_FILENAME + "." + CHR + ".tsv")


## After calculating PRS for each chromosome, we can sum them up to get the final PRS for each sample.
import pandas as pd
import statistics as st

CHR = "3"
df = pd.read_table(OUTPUT_FILENAME + "." + CHR + ".tsv")

rest_of_chroms = ["4", "5", "6", "7", "11", "13", "15", "19", "20"]

for CHR in rest_of_chroms: 
    print(CHR)
    df2 = pd.read_table(OUTPUT_FILENAME + "." + CHR + ".tsv")
    if not all(x == y for x, y in zip(df['s'], df2['s'])):
        print("ERROR: samples not matching")
        break
    df['prs'] = df['prs'] + df2['prs']

df.to_csv(OUTPUT_FILENAME + '.PRS.tsv', sep='\t', index=False)

