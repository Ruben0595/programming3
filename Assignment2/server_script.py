import multiprocessing as mp
from multiprocessing.managers import BaseManager, SyncManager
import os, sys, time, queue

POISONPILL = "MEMENTOMORI"

class Server():
    def make_server_manager(port, authkey):
        """ Create a manager for the server, listening on the given port.
            Return a manager object with get_job_q and get_result_q methods.
        """
        job_q = queue.Queue()
        result_q = queue.Queue()

        # This is based on the examples in the official docs of multiprocessing.
        # get_{job|result}_q return synchronized proxies for the actual Queue
        # objects.
        class QueueManager(BaseManager):
            pass

        QueueManager.register('get_job_q', callable=lambda: job_q)
        QueueManager.register('get_result_q', callable=lambda: result_q)

        manager = QueueManager(address=('', port), authkey=authkey)
        manager.start()
        print('Server started at port %s' % port)
        return manager


    def runserver(fn, data, portnumber):
        # Start a shared manager server and access its queues
        manager = make_server_manager(portnumber, b'whathasitgotinitspocketsesss?')
        shared_job_q = manager.get_job_q()
        shared_result_q = manager.get_result_q()
        
        if not data:
            print("Gimme something to do here!")
            return
        
        print("Sending data!")
        for d in data:
            shared_job_q.put({'fn' : fn, 'arg' : d})
        
        time.sleep(2)  
        
        results = []
        while True:
            try:
                result = shared_result_q.get_nowait()
                results.append(result)
                print("Got result!", result)
                if len(results) == len(data):
                    print("Got all results!")
                    break
            except queue.Empty:
                time.sleep(1)
                continue
        # Tell the client process no more data will be forthcoming
        print("Time to kill some peons!")
        shared_job_q.put(POISONPILL)
        # Sleep a bit before shutting down the server - to give clients time to
        # realize the job queue is empty and exit in an orderly way.
        time.sleep(5)
        print("Aaaaaand we're done for the server!")
        manager.shutdown()
        print(results)
