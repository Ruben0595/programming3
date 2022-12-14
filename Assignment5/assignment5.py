from pyspark.sql import SQLContext
from pyspark import SparkContext
from pyspark.sql import functions as f
import pandas as pd
import os
import csv
import os
  
if not os.path.exists("./output/"):
    os.makedirs("./output/")

sc = SparkContext('local[16]')
path= "/data/dataprocessing/interproscan/all_bacilli.tsv"
df = SQLContext(sc).read.csv(path, sep=r'\t', header=False, inferSchema= True)

#i got the column name list from github
df = df.withColumnRenamed('_c0', 'Protein_accession')\
    .withColumnRenamed('_c1', 'MD5')\
    .withColumnRenamed('_c2', 'Seq_len')\
    .withColumnRenamed('_c3', 'Analysis')\
    .withColumnRenamed('_c4', 'Signature_accession')\
    .withColumnRenamed('_c5', 'Signature_description')\
    .withColumnRenamed('_c6', 'Start')\
    .withColumnRenamed('_c7', 'Stop')\
    .withColumnRenamed('_c8', 'Score')\
    .withColumnRenamed('_c9', 'Status')\
    .withColumnRenamed('_c10', 'Date')\
    .withColumnRenamed('_c11', 'InterPro_accession')\
    .withColumnRenamed('_c12', 'InterPro_description')\
    .withColumnRenamed('_c13', 'GO_annotations')\
    .withColumnRenamed('_c14', 'Pathways')

#1. How many distinct protein annotations are found in the dataset? I.e. how many distinc InterPRO numbers are there?
Q1 = df.filter(df.InterPro_accession!="-").select("InterPro_accession").distinct()
E1 = Q1._sc._jvm.PythonSQLUtils.explainString(Q1._jdf.queryExecution(),"simple")
A1 = Q1.count()
print('q1: ', A1)


#2. How many annotations does a protein have on average?
Q2 = df.filter(df.InterPro_accession!="-").groupBy('Protein_accession')
Q2 =Q2.agg(f.count('InterPro_accession')).agg(f.mean('count(InterPro_accession)'))
E2 = Q2._sc._jvm.PythonSQLUtils.explainString(Q2._jdf.queryExecution(),"simple")
A2 = Q2.collect()[0][0]
print('q2: ', A2)

#3. What is the most common GO Term found?
Q3 = df.filter(df.GO_annotations!="-")
Q3 = Q3.withColumn("go_a", f.explode(f.split(df.GO_annotations, "[|]")))
Q3 = Q3.groupBy('go_a').count().orderBy('count', ascending=False)
E3 = Q3._sc._jvm.PythonSQLUtils.explainString(Q3._jdf.queryExecution(),"simple")
A3 = Q3.collect()[0][0]
print('q3: ', A3)


#4. What is the average size of an InterPRO feature found in the dataset?
Q4 = df.withColumn('size', df['Stop'] - df['Start'])
Q4 = Q4.groupBy().avg('size')
E4 = Q4._sc._jvm.PythonSQLUtils.explainString(Q4._jdf.queryExecution(),"simple")
A4 = Q4.collect()[0][0]
print('q4: ', A4)

#5. What is the top 10 most common InterPRO features?
Q5 = df.filter(df.InterPro_accession!="-")
Q5 = Q5.groupBy('InterPro_accession').count().orderBy('count', ascending=False)
E5 = Q5._sc._jvm.PythonSQLUtils.explainString(Q5._jdf.queryExecution(),"simple")
A5 = [i[0] for i in Q5.take(10)]
print('q5: ', A5)

#6. If you select InterPRO features that are almost the same size (within 90-100%) as the protein itself, what is the top10 then?
Q6 = df.filter(df.InterPro_accession!="-")
Q6 = Q6.withColumn('size_ipf', Q6['Stop'] - Q6['Start'])
Q6 = Q6.withColumn('size', Q6['size_ipf']/Q6['Seq_len'])
Q6 = Q6.filter(Q6['size'] > 0.9).orderBy('size', ascending=False)
E6 = Q6._sc._jvm.PythonSQLUtils.explainString(Q6._jdf.queryExecution(),"simple")
A6 = [i[11] for i in Q6.take(10)]
print('q6: ', A6)

#7. If you look at those features which also have textual annotation, what is the top 10 most common word found in that annotation?
Q7 = df.filter(df.InterPro_description!="-")
Q7 = Q7.withColumn("desc", f.explode(f.split(Q7.InterPro_description, " ")))
Q7 = Q7.groupBy('desc').count().orderBy('count', ascending=False)
E7 = Q7._sc._jvm.PythonSQLUtils.explainString(Q7._jdf.queryExecution(),"simple")
A7 = [i[0] for i in Q7.take(10)]
print('q7: ', A7)

#8. And the top 10 least common?
Q8 = df.filter(df.InterPro_description!="-")
Q8 = Q8.withColumn("desc", f.explode(f.split(Q8.InterPro_description, " ")))
Q8 = Q8.groupBy('desc').count().orderBy('count', ascending=True)
E8 = Q8._sc._jvm.PythonSQLUtils.explainString(Q7._jdf.queryExecution(),"simple")
A8 = [i[0] for i in Q8.take(10)]
print('q8: ', A8)

#9. Combining your answers for Q6 and Q7, what are the 10 most commons words found for the largest InterPRO features?
Q9 = df.filter(df.InterPro_accession!="-")
Q9 = Q9.withColumn('size_ipf', Q9['Stop'] - Q9['Start'])
Q9 = Q9.withColumn('size', Q9['size_ipf']/Q9['Seq_len'])
Q9 = Q9.filter(Q9['size'] > 0.9).orderBy('size', ascending=False)
Q9 = Q9.withColumn("desc", f.explode(f.split(Q9.InterPro_description, " ")))
Q9 = Q9.groupBy('desc').count().orderBy('count', ascending=False)
E9 = Q9._sc._jvm.PythonSQLUtils.explainString(Q9._jdf.queryExecution(),"simple")
A9 = [i[0] for i in Q9.take(10)]
print('q9: ', A9)

#10 What is the coefficient of correlation (R**2) between the size of the protein and the number of features found?
Q10 = df.filter(df.InterPro_accession!="-")
Q10 = Q10.groupby(Q10.Protein_accession, 'Seq_len').count()
E10 = Q10._sc._jvm.PythonSQLUtils.explainString(Q10._jdf.queryExecution(),"simple")
A10 = Q10.corr('Seq_len', 'count')**2
print('q10: ', A10)


with open('./output/output.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerow(['1', A1, E1])
    writer.writerow(['2', A2, E2])
    writer.writerow(['3', A3, E3])
    writer.writerow(['4', A4, E4])
    writer.writerow(['5', A5, E5])
    writer.writerow(['6', A6, E6])
    writer.writerow(['7', A7, E7])
    writer.writerow(['8', A8, E8])
    writer.writerow(['9', A9, E9])
    writer.writerow(['10', A10, E10])