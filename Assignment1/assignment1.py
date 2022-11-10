
from Bio import Entrez
import xml.etree.ElementTree as ET
from multiprocessing import Pool
import sys
import os

path = 'output'

# Check whether the specified path exists or not
isExist = os.path.exists(path)

if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(path)


Entrez.api_key = '553ac3e1d61512611b1e05ed4675b9366d08'
Entrez.email = 'rubenotter@gmail.com'

def get_papers(pmid):
	handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text")
	with open(f'./output/{pmid}.xml', 'wb') as file:
		file.write(handle.read())

if __name__ == '__main__':
	pmid = sys.argv[1]
	results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                   db="pmc",
                                   LinkName="pubmed_pmc_refs",
                                   id=pmid,
                                   api_key='553ac3e1d61512611b1e05ed4675b9366d08'))
	references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]

	print(references[:10])
	with Pool(8) as p:
		p.map(get_papers, references[:10])