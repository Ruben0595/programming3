import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue
import server_script as sv
import client_script as cl
import argparse
from Bio import Entrez
import pickle

Entrez.api_key = '553ac3e1d61512611b1e05ed4675b9366d08'
Entrez.email = 'rubenotter@gmail.com'

path = 'output'
# Check whether the specified path exists or not
isExist = os.path.exists(path)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(path)

parser = argparse.ArgumentParser()
parser.add_argument("-n", help="number of peons", type=int, required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-c", help="client mode", action='store_true')
group.add_argument("-s", help="server mode", action='store_true')
parser.add_argument("--port", help="port number", type=int, required=True)
parser.add_argument("--host", help="host number", type=str, required=True)
parser.add_argument("-a", help="number of articles", type=int, required=True)
parser.add_argument("id", help="PubMedID", type=int)
args = parser.parse_args()


def capitalize(word):
    """Capitalizes the word you pass in and returns it"""
    return word.upper()

data = ["Always", "look", "on", "the", "bright", "side", "of", "life!"]

def get_references(pmid):
    results = Entrez.read(Entrez.elink(dbfrom="pubmed",
                                   db="pmc",
                                   LinkName="pubmed_pmc_refs",
                                   id=pmid,
                                   api_key='553ac3e1d61512611b1e05ed4675b9366d08'))
    references = [f'{link["Id"]}' for link in results[0]["LinkSetDb"][0]["Link"]]
    return references


def get_authors(pubid):
        """
        Get the author of a provided pubmed id.
        """
        sep = os.path.sep

        try:
            handle = Entrez.esummary(db="pubmed", id=pubid, retmode="xml")
            records = Entrez.parse(handle)
            for record in records:
                pickle.dump(tuple(record['AuthorList']), open(f'output{sep}{pubid}.authors.pickle', 'wb'))
        except RuntimeError:
            pickle.dump(f"Authors for record {pubid} not found in database", open(f'output{sep}{pubid}.authors.pickle', 'wb'))


def get_papers(pmid):
    get_authors(pmid)
    handle = Entrez.efetch(db="pmc", id=pmid, rettype="XML", retmode="text")
    with open(f'output/{pmid}.xml', 'wb') as file:
        file.write(handle.read())

POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
AUTHKEY = b'whathasitgotinitspocketsesss?'

id_list = get_references(args.id)

if args.s:
    server_side = sv.Server(args.port, POISONPILL)
    server = mp.Process(target=server_side.runserver, args=(get_papers, id_list[:10], args.port, args.host))
    server.start()
    time.sleep(1)
    server.join()

if args.c:
    client_side = cl.Client(args.host, args.port, AUTHKEY, POISONPILL, ERROR)
    client = mp.Process(target=client_side.runclient, args=(4,))
    client.start()
    client.join()

