import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue

POISONPILL = "MEMENTOMORI"
ERROR = "DOH"

class Client():
    def __init__(self, ip, port, authkey, poisonpill, error):
        self.ip = ip
        self.port = port
        self.authkey = authkey
        self.poisonpill = poisonpill
        self.error = error

    def make_client_manager(self, ip, port, authkey):
        """ Create a manager for a client. This manager connects to a server on the
            given address and exposes the get_job_q and get_result_q methods for
            accessing the shared queues from the server.
            Return a manager object.
        """
        class ServerQueueManager(BaseManager):
            pass

        ServerQueueManager.register('get_job_q')
        ServerQueueManager.register('get_result_q')

        manager = ServerQueueManager(address=(ip, port), authkey=authkey)
        manager.connect()

        print('Client connected to %s:%s' % (ip, port))
        return manager



    def runclient(self, num_processes):
        manager = self.make_client_manager(self.ip, self.port, self.authkey)
        job_q = manager.get_job_q()
        result_q = manager.get_result_q()
        self.run_workers(job_q, result_q, num_processes)
        
    def run_workers(self, job_q, result_q, num_processes):
        processes = []
        for p in range(num_processes):
            temP = mp.Process(target=self.peon, args=(job_q, result_q))
            processes.append(temP)
            temP.start()
        print("Started %s workers!" % len(processes))
        for temP in processes:
            temP.join()

    def peon(self, job_q, result_q):
        my_name = mp.current_process().name
        while True:
            try:
                job = job_q.get_nowait()
                if job == POISONPILL:
                    job_q.put(POISONPILL)
                    print("Aaaaaaargh", my_name)
                    return
                else:
                    try:
                        result = job['fn'](job['arg'])
                        print("Peon %s Workwork on %s!" % (my_name, job['arg']))
                        result_q.put({'job': job, 'result' : result})
                    except NameError:
                        print("Can't find yer fun Bob!")
                        result_q.put({'job': job, 'result' : ERROR})

            except queue.Empty:
                print("sleepytime for", my_name)
                time.sleep(1)