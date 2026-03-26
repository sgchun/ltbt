Sign into https://ukbiobank.dnanexus.com/ 
Select “Tools” on top 
Select “Jupyter Lab” 
Select “New Jupyter Lab” on the top right corner
Add the “Feature” of “HAIL-0.2.61” 
Launch and Open
Select “Python 3” Notebook



# Initialize
import pyspark
sc = pyspark.SparkContext()

import hail as hl
hl.init(sc=sc)

# DNA Nexus control
import dxpy

projectid = "project-G48jpz8J7ZfgQB5jJ35VJKfF"

# File
# indir="file:///mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/"
indir="file:///mnt/project/Work/"

# filepre="ukb23156_c22_b0_v1"
filepre="ukb23156_c22_pilot"


# Load vcf file 
mt = hl.import_vcf(indir+filepre+".vcf.gz", reference_genome="GRCh38", force_bgz=True, array_elements_required=False)

mt.describe()
# print(mt.count_rows())
print(mt.count_cols())

# mt.filter_rows(mt.filters.length() > 0).filters.show(10)

# British Whites only
with open("file:///mnt/project/Work/ukb47653.british.eid") as file: 
    lines = file.readlines()
    lines = [line.rstrip() for line in lines]

print(len(lines))
british = hl.literal(set(lines))
mt = mt.filter_cols(british.contains(mt['s']))

# print(mt.count_rows())
# print(mt.count_cols())


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

# mt.show()

# print(mt.count_rows())


# FILTER: >= 90% samples with DP >= 10 
mt = mt.filter_rows(hl.agg.mean(mt.DP >= 10) >= 0.9)

# print(mt.count_rows())
# print(mt.count_cols())

# mt.DP.show(n_rows=5, n_cols=100)


# mt.AD.show(n_rows=5, n_cols=100)
# mt.GT.show(n_rows=5, n_cols=100)

# write back to vcf
# vcfoutfile = "/opt/notebooks/" + filepre +".varqc.vcf.gz"
# hl.export_vcf(mt, "file://" + vcfoutfile)

# write back as a hail table

# The following line create the spark db named “QC” and “test” (I have already created, no need to repeat)
# spark = pyspark.sql.SparkSession(sc)
# spark.sql("CREATE DATABASE QC LOCATION  'dnax://'")
# spark.sql("CREATE DATABASE test LOCATION  'dnax://'")

# db ID for QC 
# dbid =  "database-G5P0jzQJ7Zfqv2QV8ByqKXX5"

# db ID for test 
dbid = "database-G5P0QPjJ7ZfYvG6421P0P3fb"

mt.write("dnax://" + dbid + "/" + filepre + ".mt", overwrite=True)

# Write sample QC info
qcoutfile = "/opt/notebooks/" + filepre +".sample_qc.tsv"
mt.sample_qc.export("file://" + qcoutfile)

callratefile = "/opt/notebooks/" + filepre +".sample_call_rate.tsv"
mt.sample_qc.call_rate.export("file://" + callratefile)

# Upload
dxpy.upload_local_file(qcoutfile)
dxpy.upload_local_file(callratefile)










