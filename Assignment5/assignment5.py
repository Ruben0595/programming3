from pyspark.sql.types import StructType, StructField, IntegerType, StringType,FloatType
from pyspark.sql.functions import *
from pyspark.sql import SparkSession
import pandas as pd
import os

# special thxs to kai for the explain function. :)
def explain(data):
    return data._sc._jvm.PythonSQLUtils.explainString(data._jdf.queryExecution(), 'simple')

def spark_df():
    path = "/data/dataprocessing/interproscan/all_bacilli.tsv"
    schema = StructType([
    StructField("Protein_accession", StringType(), True),
    StructField("Sequence_MD5_digest", StringType(), True),
    StructField("Sequence_length", IntegerType(), True),
    StructField("Analysis", StringType(), True),
    StructField("Signature_accession", StringType(), True),
    StructField("Signature_description", StringType(), True),
    StructField("Start_location", IntegerType(), True),
    StructField("Stop_location", IntegerType(), True),
    StructField("Score", FloatType(), True),
    StructField("Status", StringType(), True),
    StructField("Date", StringType(), True),
    StructField("InterPro_annotations_accession", StringType(), True),
    StructField("InterPro_annotations_description", StringType(), True),
    StructField("GO_annotations", StringType(), True),
    StructField("Pathways_annotations", StringType(), True)])
    spark = SparkSession.builder.master("local[16]").appName("InterPro").getOrCreate()
    return spark.read.option("sep", "\t").option("header", "False").csv(path, schema=schema)

#1 How many distinct protein annotations are found in the dataset? I.e. how many distinc InterPRO numbers are there?
def q1(df):
    print("Awsering question 1")
    data1 = df.select('InterPro_annotations_accession')\
        .filter(df.InterPro_annotations_accession != "-")\
        .distinct()
    exp1 = explain(data1)
    data1 = data1.count()

    return data1, exp1

#2 How many annotations does a protein have on average?
def q2(df):
    print("Awsering question 2")
    data2 = df.select("Protein_accession",'InterPro_annotations_accession')\
                .filter(df.InterPro_annotations_accession != "-")\
                .groupBy("Protein_accession")\
                .count()\
                .select(mean("count"))
    exp2 = explain(data2)      
    data2 = data2.collect()[0].__getitem__(0)

    return data2, exp2

#3 What is the most common GO Term found?
def q3(df):
    print("Awsering question 3")
    data3 = df.select(df.GO_annotations, explode(split(col("GO_annotations"),"\|"))\
                        .alias("Split_col"))
    data3 = data3.filter(data3.Split_col != "-")\
                .select("Split_col")\
                .groupby("Split_col")\
                .count()\
                .sort("count",ascending=False)
    exp3 = explain(data3) 
    data3 = [data[0] for data in data3.take(1)]
    data3 = data3[0]

    return data3, exp3

#4 What is the average size of an InterPRO feature found in the dataset?
def q4(df):
    print("Awsering question 4")
    data4 = df.withColumn('Sub', (df['Stop_location'] - df['Start_location'])).summary("mean")
    exp4 = explain(data4)  
    data4 = data4.collect()[0].__getitem__(-1)

    return data4, exp4

#5 What is the top 10 most common InterPRO features?
def q5(df):
    print("Awsering question 5")
    data5 = df.select('InterPro_annotations_accession')\
        .filter(df.InterPro_annotations_accession != "-")\
        .groupBy('InterPro_annotations_accession')\
        .count()\
        .sort("count",ascending=False)\
        .select("InterPro_annotations_accession")
    exp5 = explain(data5) 
    data5 = [data[0] for data in data5.take(10)]

    return data5, exp5

#6 If you select InterPRO features that are almost the same size (within 90-100%) as the protein itself, what is the top10 then?
def q6(df):
    print("Awsering question 6")
    data6 = df.select('InterPro_annotations_accession',"Sequence_length",'Stop_location','Start_location')\
                .filter((df['Stop_location'] - df['Start_location'])/df["Sequence_length"]>=0.9)\
                .filter(df.InterPro_annotations_accession != "-")\
                .groupBy('InterPro_annotations_accession')\
                .count()\
                .sort("count",ascending=False)\
                .select("InterPro_annotations_accession")
    exp6 = explain(data6)
    data6 = [data[0] for data in data6.take(10)]
 
    return data6, exp6

#7 If you look at those features which also have textual annotation, what is the top 10 most common word found in that annotation?
def q7(df):
    print("Awsering question 7")
    data7 = df.select(df.InterPro_annotations_description,explode(split(col("InterPro_annotations_description")," |,"))\
                .alias("Split_col"))
    data7 = data7.select("Split_col")\
                .filter(data7.Split_col != "")\
                .filter(data7.Split_col != "-")\
                .groupby("Split_col")\
                .count()\
                .sort("count",ascending=False)\
                .select("Split_col")
    exp7 = explain(data7)
    data7 = [data[0] for data in data7.take(10)]

    return data7,exp7

#8 top 10 least common?
def q8(df):
    print("Awsering question 8")
    data8 = df.select(df.InterPro_annotations_description,explode(split(col("InterPro_annotations_description")," |,"))\
                .alias("Split_col"))
    data8 = data8.select("Split_col")\
                .filter(data8.Split_col != "")\
                .filter(data8.Split_col != "-")\
                .groupby("Split_col")\
                .count()\
                .sort("count",ascending=True)\
                .select("Split_col")
    exp8 = explain(data8) 
    data8 = [data[0] for data in data8.take(10)]

    return data8, exp8

#9 Combining your answers for Q6 and Q7, what are the 10 most commons words found for the largest InterPRO features?
def q9(df):
    print("Awsering question 9")
    data9 = df.select(df.InterPro_annotations_accession,df.InterPro_annotations_description)\
                .filter(df.InterPro_annotations_accession.isin(data6))\
                .distinct()
    data9 = data9.select(data9.InterPro_annotations_description,explode(split(col("InterPro_annotations_description")," |,")))\
                    .groupby("col")\
                    .count()
    data9 = data9.select(data9["col"], data9["count"])\
                    .filter(data9["col"] != "")\
                    .sort("count",ascending=False)
    exp9 = explain(data9)
    data9 = [data[0] for data in data9.take(10)]
    
    return data9, exp9

#10 What is the coefficient of correlation (R^2) between the size of the protein and the number of features found?
def q10(df):
    print("Awsering question 10")
    data10=df.select(df.Protein_accession,df.InterPro_annotations_accession,df.Sequence_length)\
                .filter(df.InterPro_annotations_accession != "-")\
                .groupby(df.Protein_accession,"Sequence_length")\
                .count()
    exp10 = explain(data10)
    data10 = data10.corr('Sequence_length', 'count')**2

    return data10, exp10


if __name__ == "__main__":
    df = spark_df()
    # assign the questions to the functions and the explain data
    data1, exp1 =  q1(df)
    data2, exp2 =  q2(df)
    data3, exp3 =  q3(df)
    data4, exp4 =  q4(df)
    data5, exp5 =  q5(df)
    data6, exp6 =  q6(df)
    data7, exp7 =  q7(df)
    data8, exp8 =  q8(df)
    data9, exp9 =  q9(df)
    data10,exp10 = q10(df)
    # create the columns with the awsers, explain data and the question number. 
    column1 = list(range(1,11))
    column2 = [data1,data2,data3,data4,data5,data6,data7,data8,data9,data10]
    column3 = [exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9,exp10]
    # write columns to csv
    try:
        d = {'Question': column1, 'Answer': column2,"Explain":column3}
        df = pd.DataFrame(data = d)
        if not os.path.exists("output"):
            os.makedirs("output")
        df.to_csv("output/output.csv",index=False)
    except Exception:
        print("Error while writing to csv")
