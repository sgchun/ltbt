Sign into https://ukbiobank.dnanexus.com/ 
Select “Tools” on top 
Select “Jupyter Lab” 
Select “New Jupyter Lab” on the top right corner
Select the “Feature” of “HAIL-0.2.61-VEP-1.0.3” 
Launch and Open
Select “Python 3” Notebook



# Initialize
import pyspark
sc = pyspark.SparkContext()

import hail as hl
hl.init(sc=sc)

# spark = pyspark.sql.SparkSession(sc)
# spark.sql("CREATE DATABASE VEP LOCATION  'dnax://'")

# db ID for VEP output
# dbid =  "database-G5PZ0v8J7Zfy39z97vbY39p6"

# db ID for test only 
dbid = "database-G5P0QPjJ7ZfYvG6421P0P3fb"

# File
# filepre="ukb23156_c22_b0_v1"
filepre="ukb23156_c22_pilot"

# Read from Spark database
mt2 = hl.read_matrix_table("dnax://" + dbid + "/" + filepre + ".mt")

# VEP annotation
mt2 = hl.methods.vep(mt2, config="file:///mnt/project/Work/vep.json")

# Upload back

mt.write("dnax://" + dbid + "/" + filepre + ".mt", overwrite=True)
