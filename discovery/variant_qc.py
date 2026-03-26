# Initialize
import pyspark
sc = pyspark.SparkContext()

import hail as hl
hl.init(sc=sc)

# DNA Nexus control
import dxpy

# Load vcf file 
mt0 = hl.import_vcf(INPUT_VCF_FILE, reference_genome="GRCh38", 
                    force_bgz=True, array_elements_required=False)

mt0.describe()
print(mt0.count_rows())
print(mt0.count_cols())

# For UK Biobank, QC was done with White British subjects only
with open("file:///mnt/project/Work/white_british.eid") as file: 
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

print(len(lines))
british = hl.literal(set(lines))
mt = mt0.filter_cols(british.contains(mt0['s']))

# For sample QC
mt = hl.sample_qc(mt)

# Filter genotype entries by allele balance check 
AB = mt.AD[1] / hl.sum(mt.AD)

filter_condition_AB = (
    hl.case()
    .when(mt.GT.is_hom_ref(), AB < 0.1)
    .when(mt.GT.is_het(), (AB > 0.2) & (AB < 0.8))
    .default(AB > 0.9) # hom-var
)

mt = mt.filter_entries(filter_condition_AB)

# For variant QC
# FILTER: by HWE P and missing rate 
mt = hl.variant_qc(mt)

mt.variant_qc.describe()

mt = mt.filter_rows((mt.variant_qc.p_value_hwe > 1e-15) & (mt.variant_qc.call_rate > 0.9))

# FILTER: >= 90% samples with DP >= 10 
mt = mt.filter_rows(hl.agg.mean(mt.DP >= 10) >= 0.9)

# write to a hail table
dbid = "YOUR DATABASE ID"

mt.write("dnax://" + dbid + "/" + OUTPUT_FILE + ".mt", overwrite=True)

# Identify QC failed variants by anti-joining with the original vcf file
qcfails = mt0.anti_join_rows(mt.rows())
qcfails = qcfails.rows()
qcfails_df = qcfails.select(qcfails.rsid, qcfails.qual).to_pandas()
qcfails_df.to_csv(QCFAIL_OUTPUT_FILE)

