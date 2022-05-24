import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue
import server_script as sv
import client_script as cl
import argparse


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

POISONPILL = "MEMENTOMORI"
ERROR = "DOH"
#Have to log into this nuc for it to work
#use same host for server and client
IP = args.host
# PORTNUM = 5345
AUTHKEY = b'whathasitgotinitspocketsesss?'
PORTNUM = args.port


if args.s:
    server = sv.Server(PORTNUM, POISONPILL)
    server.runserver(downloader.download_files, ids[:args.a])

if args.c:
    client = cl.Client(IP, PORTNUM, AUTHKEY, POISONPILL, ERROR)
    client.runclient(args.n)