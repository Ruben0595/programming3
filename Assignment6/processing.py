import dask.dataframe as dd
from dask_ml import preprocessing

def file_cleaner(df):
    columns = ["Protein_accession", "Sequence_MD5_digest", "Sequence_length", "Analysis", "Signature_accession", "Signature_description", "Start_location", "Stop_location", "Score", "Status", "Date", "InterPro_annotations_accession", "InterPro_annotations_description", "GO_annotations", "Pathways_annotations"]
    df.columns = columns
    data = df.drop(["Analysis", "Sequence_MD5_digest", "Signature_accession", "Signature_description", "Score", "Status", "Date", "InterPro_annotations_description", "GO_annotations", "Pathways_annotations"], axis=1)
    data = data.dropna()
    '''drop every row containing '-' in the column 'InterPro_annotations_accession'''
    data = data[data['InterPro_annotations_accession'] != '-']
    return data

def find_long_short(data):
    '''create column, which substracts columns "Start_location" and "Stop_location"'''
    data['Length %'] = (data['Stop_location'] - data['Start_location']) / data['Sequence_length']

    '''if length % is greater than 0.9, greather_than column is 1, else 0'''
    data['greater_than'] = data['Length %'] >= 0.9

    '''make protein_90% df with column "Length %" > 0.9'''
    protein_90 = data[data['Length %'] >= 0.9]
    protein_0 = data[data['Length %'] < 0.9]

    protein_90 = protein_90.drop(['Sequence_length', 'Start_location', 'Stop_location', 'greater_than', 'Length %'], axis=1)
    protein_0 = protein_0.drop(['Sequence_length', 'Start_location', 'Stop_location', 'greater_than', 'Length %'], axis=1)

    protein_90 = protein_90.rename(columns = {'InterPro_annotations_accession': 'Long'})
    protein_0 = protein_0.rename(columns = {'InterPro_annotations_accession': 'Short'})

    
    return protein_0, protein_90

def create_matrix(data):
    data['count'] = 1
    data = data.drop('Long', axis = 1)
    data = data.categorize(columns=['Short'])
    Matrix = dd.reshape.pivot_table(data, index="Protein_accession",
                                    columns='Short', values="count",
                                    aggfunc='sum')
    return Matrix

def create_feats_lbls(matrix):
    x = matrix.drop(matrix.columns[[0, -1]], axis=1)
    y = matrix[["Long"]]

    y = y.categorize()
    y = preprocessing.OneHotEncoder().fit_transform(y)
    return x, y