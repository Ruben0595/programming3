'''imports'''
import dask.dataframe as dd
from dask_ml import model_selection
from sklearn import metrics, ensemble
import processing

'''open sm-all_bacilli.tsv using dask'''
df = dd.read_csv('sm-all_bacilli.tsv', sep='\t', header=None, dtype={8: 'object'})

'''extract data we need'''
data = processing.file_cleaner(df)

'''find the long and short proteins'''
short, long = processing.find_long_short(data)

'''merging the long and short protein accessions, so further processing is eassier'''
merged = long.merge(short, on='Protein_accession', how='outer')
Long = merged.drop('Short', axis = 1)

'''create a matrix, which contains the information which short a...'''
matrix = processing.create_matrix(merged)

matrix = matrix.merge(Long, how="left", on=["Protein_accession"])
matrix = matrix.dropna()
print(matrix.head())
'''extract features and labels from the matrix'''
x,y = processing.create_feats_lbls(matrix)

'''split the features and labels in a train/test data set, ratio 70/30 respectively'''
X_train, X_test, y_train, y_test = model_selection.train_test_split(x, y, test_size=0.3)

'''train the model'''
clf = ensemble.RandomForestClassifier(n_estimators=200, n_jobs=-1)
clf.fit(X_train, y_train)

'''predict new values'''
y_pred = clf.predict(X_test)

print(metrics.accuracy_score(y_test, y_pred))