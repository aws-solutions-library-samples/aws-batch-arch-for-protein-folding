import os
import platform
import subprocess
from collections import defaultdict
import time
import psutil
import threading
import hashlib
import json
import numpy as np

_run_monitor_signal = True



class InfoCollector(object):
    """An object for collecting information about a job to use in reporting and record keeping.
    It is mostly a container for information, but it can also be used to generate a report.
    
    This one will be written specifically for collecting info about an AF2 job.
    """

    _protocol_version = '1.0'

    def __init__(self, tag=None, dropoff_dir="/opt/superfold/info_collection", collect_data_file="/opt/superfold/info_collection/config.json") -> None:
        """Initialize an InfoCollector object.
        JSON fields:
        username: STRING
        timestamp: INT
        runtime: FLOAT
        compute-node: STRING
        num-cpu: INT
        cpu-ids: LIST of INT
        cpu-types: LIST of STRING
        gpu-types: LIST of STRING
        gpu-ids: LIST of STRING
        script: STRING
        interpreter: STRING
        tag: STRING
        max_system_memory: INT
        max_gpu_memory: INT
        af2-version: STRING
        sequence: STRING
        padded-length: INT
        seed: INT
        input-pdb: STRING
        source: STRING
        pLDDT: LIST of FLOAT
        LDDT: FLOAT
        model-num: INT
        num-recycles: INT
        pae-matrix: LIST of LIST of FLOAT
        rmsd: FLOAT
        TMscore: FLOAT
        pTMscore: FLOAT
        iptm: FLOAT
        used-msa: BOOL
        used-initial-guess: BOOL
        used-templates: BOOL
        output-number: INT

        """

        self._dropoff_dir = dropoff_dir
        self._collect_data_file = collect_data_file
        self._check_and_create_collect_data_file()

        collection_configuration = json.load(open(self._collect_data_file))
        self._run_collection_flag = collection_configuration['run_collection']
        if not self._run_collection_flag:
            return

        self._time_start = time.time()
        self._process = psutil.Process(os.getpid())

        process_info = self._process.as_dict()
        # for k,v in process_info.items():
        #     print(k, v, sep=': ')
        #     print("="*50)
        try:
            interpreter = process_info['environ']['SINGULARITY_CONTAINER']
        except KeyError:
            interpreter = process_info['exe']
        script = process_info['cmdline'][1]
        #get the absolute path to the script
        script = os.path.abspath(script)

        self._info = {}
        self._info['protocol-version'] = self._protocol_version
        self._info['username'] = os.environ['USER']
        self._info['timestamp'] = None
        self._info['runtime'] = None
        self._info['compute-node'] = platform.node()
        self._info['num-cpu'] = len(os.sched_getaffinity(0))
        self._info['cpu-ids'] = sorted(list(os.sched_getaffinity(0)))

        cpu_types = []
        my_processors = self._info['cpu-ids']
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).decode().strip()
        cpu_types_by_id = {}
        for processor_lines in all_info.split('\n\n'):
            processor_dict = {}
            for line in processor_lines.split('\n'):
                if line == '':
                    continue
                key, value = line.split(':')
                processor_dict[key.strip()] = value.strip()
            cpu_types_by_id[int(processor_dict['processor'])] = processor_dict['model name']

        for processor_id in my_processors:
            cpu_types.append(cpu_types_by_id[processor_id])
        self._info['cpu-types'] = cpu_types

        self._info['gpu-types'] = [g['name'] for g in self._get_gpu_info()]
        self._info['gpu-ids'] = [g['uuid'] for g in self._get_gpu_info()]
        self._info['script'] = script
        self._info['interpreter'] = interpreter
        self._info['tag'] = tag

        self._info['max_system_memory'] = -1
        self._info['max_gpu_memory'] = -1



        #the remainder get initized to none and will have to be implemented in the af2 script
        self._info['af2-version'] = None
        self._info['sequence'] = None
        self._info['padded-length'] = None
        self._info['seed'] = None
        self._info['input-pdb'] = None
        self._info['source'] = None
        self._info['pLDDT'] = None
        self._info['LDDT'] = None
        self._info['model-num'] = None
        self._info['num-recycles'] = None
        self._info['pae-matrix'] = None
        self._info['rmsd'] = None
        self._info['mmalign_rmsd'] = None
        self._info['TMscore'] = None
        self._info['pTMscore'] = None
        self._info['iptm'] = None
        self._info['used-msa'] = None
        self._info['used-initial-guess'] = None
        self._info['used-templates'] = None
        self._info['output-number'] = None
    
    def _check_and_create_collect_data_file(self):
        if not os.path.exists(self._collect_data_file):
            with open(self._collect_data_file, 'w') as file:
                json.dump({"run_collection": True}, file)
                
    def _get_system_memory_usage(self) -> int:
        """Get the current system memory usage in bytes."""
        process = self._process
        try:
            return process.memory_info().rss
        except psutil.NoSuchProcess:
            return None
        except psutil.AccessDenied:
            return None

    def _get_gpu_info(self) -> dict:
        """Get information about the GPUs available to this process."""
        try:
            # Run the nvidia-smi command
            result = subprocess.check_output(["nvidia-smi", "--query-gpu=name,uuid,memory.free,memory.total", "--format=csv,noheader,nounits"])

            # Split the output into lines
            lines = result.decode('utf-8').strip().split('\n')

            # Parse the output and store it in a list of dictionaries
            gpu_info = []
            for line in lines:
                parts = line.strip().split(', ')
                if len(parts) == 4:
                    name, uuid, free_memory, total_memory = parts
                    gpu_info.append({
                        'name': name,
                        'uuid': uuid,
                        'free_memory': int(free_memory),
                        'total_memory': int(total_memory),
                        'used_memory': int(total_memory) - int(free_memory)
                    })

            return gpu_info

        except Exception as e:
            #print("Error: ", e)
            return []

    def _get_gpu_memory_usage(self) -> int:
        """Get the current GPU memory usage in MiB."""
        gpu_info = self._get_gpu_info()
        if gpu_info:
            return gpu_info[0]['used_memory']
        else:
            return None 

    def _start_memory_monitor(self) -> None:
        """Start monitoring GPU & CPU usage.
        """
        # Create a lock to protect the dict
        self._monitor_lock = threading.Lock()
        # shared dict to store data
        self._monitor_info = {}

        def check_memory():
            #while self._run_monitor:
            global _run_monitor_signal
            while _run_monitor_signal:
                print('checking memory')
                with self._monitor_lock:
                    current_sysmem = self._get_system_memory_usage()
                    current_gpumem = self._get_gpu_memory_usage()
                    if current_sysmem:
                        self._info['max_system_memory'] = max(self._info['max_system_memory'], current_sysmem)

                    if current_gpumem:
                        self._info['max_gpu_memory'] = max(self._info['max_gpu_memory'], current_gpumem)

                time.sleep(1)
            print('stopped checking memory')

        # check that a monitor thread is not already running
        if not hasattr(self, '_monitor_thread'):
            # Create a thread to monitor memory usage
            #self._run_monitor = True
            global _run_monitor_signal
            _run_monitor_signal = True
            self._monitor_thread = threading.Thread(target=check_memory)
            self._monitor_thread.start()

    def _stop_memory_monitor(self) -> None:
        """Stop monitoring GPU & CPU usage.

        Stops the subprocess that is monitoring GPU and CPU usage.
        """

        # Stop the thread by setting the control variable
        # self._run_monitor = False
        global _run_monitor_signal
        _run_monitor_signal = False
       
        #stop the thread
        if hasattr(self, '_monitor_thread'):           
            # Wait for the thread to finish with a timeout
            self._monitor_thread.join(timeout=1)

        #take a final measurement
        current_sysmem = self._get_system_memory_usage()
        current_gpumem = self._get_gpu_memory_usage()
        if current_sysmem:
            self._info['max_system_memory'] = max(self._info['max_system_memory'], current_sysmem)
        if current_gpumem:
            self._info['max_gpu_memory'] = max(self._info['max_gpu_memory'], current_gpumem)

    def __getitem__(self, key):
        """Get an item from the info dict."""
        return self._info[key]
    
    def __setitem__(self, key, value):
        """Set an item in the info dict."""
        if key in self._info:
            self._info[key] = value
        else :
            raise KeyError(f"Key {key} not found in info dict and cannot be set.")

    def report(self) -> None:
        """Print a report of the information collected."""
        if not self._run_collection_flag:
            return
        
        self._stop_memory_monitor()
        stop_time = time.time()
        self._info['runtime'] = stop_time - self._time_start
        #record timestamp as unix epoch time
        self._info['timestamp'] = int(time.time())

        # cast devicearray to serializable type
        for key in self._info:
            self._info[key] = np.array(self._info[key]).tolist()

        file_contents = json.dumps(self._info)

        unique_filename = hashlib.md5(file_contents.encode()).hexdigest() + ".json"

        with open(self._dropoff_dir + "/" + unique_filename, 'w') as f:
            f.write(file_contents)

def main():
    """Main function for testing."""
    info_collector = InfoCollector()
    time.sleep(0.5)
    info_collector.report()


if __name__ == '__main__':
    main()