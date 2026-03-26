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
row = 0
for row in range(0, len(fileList)):
    print(row)
    cFile = fileList[row]
    filepre = os.path.basename(cFile.split(".")[0])
    print(filepre)
    if "c23" in filepre:
        continue
    if "c24" in filepre:
        continue
    ukb = hl.read_matrix_table("dnax://" + dbid + "/" + filepre + ".mt")

    ## Synonymous
    ukb2 = ukb.filter_rows((ukb.vep.most_severe_consequence == 'synonymous_variant') &
                           (ukb.variant_qc.AF[1] > 0)) 

    vcfoutfile = "/opt/notebooks/" + filepre +".QC.syn.vcf.gz"
    hl.export_vcf(ukb2, "file://" + vcfoutfile)
    dxpy.upload_local_file(vcfoutfile, folder="/Work/", parents=True)

    ukb2syn = ukb2.annotate_rows(
        vep=ukb2.vep.annotate(
            transcript_consequences=ukb2.vep.transcript_consequences.filter(
                lambda csq: csq.consequence_terms.contains('synonymous_variant')))
    )

    ukb2syn = ukb2syn.annotate_rows(gene=ukb2syn.vep.transcript_consequences[0].gene_id,
                                    AF=ukb2syn.variant_qc.AF[1])
    ukb2syn = ukb2syn.rows()
    syndf = ukb2syn.select(ukb2syn.gene, ukb2syn.AF).to_pandas()
    # syndf.head()

    syndf.to_csv(filepre + ".syn_all.csv")
    dxpy.upload_local_file(filepre + ".syn_all.csv", folder="/Work/", parents=True)

    ## Damaging
    ukb2 = ukb.filter_rows((ukb.vep.most_severe_consequence == 'missense_variant') &
                           (ukb.variant_qc.AF[1] > 0))

    ukb2dam = ukb2.annotate_rows(
                vep=ukb2.vep.annotate(
                  transcript_consequences=ukb2.vep.transcript_consequences.filter(
                    lambda csq: csq.polyphen_prediction.contains('probably_damaging')))
        )

    ukb2dam = ukb2dam.filter_rows(ukb2dam.vep.transcript_consequences.length() > 0)
    
    vcfoutfile = "/opt/notebooks/" + filepre +".QC.probdam.vcf.gz"
    hl.export_vcf(ukb2dam, "file://" + vcfoutfile)
    dxpy.upload_local_file(vcfoutfile, folder="/Work/", parents=True)
    
    ukb2dam = ukb2dam.annotate_rows(gene=ukb2dam.vep.transcript_consequences[0].gene_id,
                                      AF=ukb2dam.variant_qc.AF[1])
    ukb2dam = ukb2dam.rows()
    damdf = ukb2dam.select(ukb2dam.gene, ukb2dam.AF).to_pandas()
    # damdf.head()

    damdf.to_csv(filepre + ".probdam_all.csv")
    dxpy.upload_local_file(filepre + ".probdam_all.csv", folder="/Work/", parents=True)

    ## LoF
    ukb2 = ukb.filter_rows((ukb.variant_qc.AF[1] > 0))

    ukb2lof = ukb2.annotate_rows(
        vep=ukb2.vep.annotate(
            transcript_consequences=ukb2.vep.transcript_consequences.filter(
                lambda csq: csq.lof.contains('HC')))
    )

    ukb2lof = ukb2lof.filter_rows(ukb2lof.vep.transcript_consequences.length() > 0)
    
    vcfoutfile = "/opt/notebooks/" + filepre +".QC.hclof.vcf.gz"
    hl.export_vcf(ukb2lof, "file://" + vcfoutfile)
    dxpy.upload_local_file(vcfoutfile, folder="/Work/", parents=True)
    
    ukb2lof = ukb2lof.annotate_rows(gene=ukb2lof.vep.transcript_consequences[0].gene_id,
                                      AF=ukb2lof.variant_qc.AF[1],
                                      gene_symbol=ukb2lof.vep.transcript_consequences[0].gene_symbol
                                      )
    ukb2lof = ukb2lof.rows()
    lofdf = ukb2lof.select(ukb2lof.gene, ukb2lof.AF, ukb2lof.gene_symbol).to_pandas()
    lofdf.head()

    lofdf.to_csv(filepre + ".hclof_all.csv")
    dxpy.upload_local_file(filepre + ".hclof_all.csv", folder="/Work/", parents=True)

