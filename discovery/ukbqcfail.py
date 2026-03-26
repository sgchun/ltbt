start = 11 

# Initialize
import pyspark
sc = pyspark.SparkContext()

import dxpy

import hail as hl
hl.init(sc=sc)

################# ukb
dbid = "database-G5jjVB0J7Zfgy2j726Q5g86X"

import glob, os
inputFolder="//mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format/"
fileList = glob.glob(inputFolder+"ukb23156_c*.vcf.gz")
print(len(fileList))

# Read from Spark database
for row in range(start, start + 20):
    print(row)
    cFile = fileList[row]
    filepre = os.path.basename(cFile.split(".")[0])
    print(filepre)
    if "c23" in filepre:
        continue
    if "c24" in filepre:
        continue
    if os.path.exists("//mnt/project/Work/" + filepre + ".qcfail.csv"):
        continue
    ukbqc = hl.read_matrix_table("dnax://" + dbid + "/" + filepre + ".mt")
    ukb0 = hl.import_vcf("file:/"+cFile, reference_genome="GRCh38", force_bgz=True, array_elements_required=False)
    
    qcfails = ukb0.anti_join_rows(ukbqc.rows())
    qcfails = qcfails.rows()
    qcfailsdf = qcfails.select(qcfails.rsid, qcfails.qual).to_pandas()
    print(qcfailsdf.shape)
    print(qcfailsdf.head())
    
    qcfailsdf.to_csv(filepre + ".qcfail.csv")
    dxpy.upload_local_file(filepre + ".qcfail.csv", folder="/Work/", parents=True)
