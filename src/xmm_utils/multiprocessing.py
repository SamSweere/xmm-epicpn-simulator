import multiprocessing as mp
from datetime import datetime
from multiprocessing.pool import Pool
from typing import Dict

from psutil import virtual_memory


def get_num_processes(mp_conf: Dict[str, float]):
    # num_processes : the number of processes, if None the calculation is done based on the gb_per_process limit
    # gb_per_process : the expected amount of ram one process needs, if None there is no limit
    num_processes = mp_conf["num_processes"]
    gb_per_process = mp_conf["gb_per_process"]

    cpu_count = mp.cpu_count()
    mem = virtual_memory().total * 1e-9  # System memory in gb
    print(f"CPU count: {cpu_count}\nSystem memory: {mem}")

    if num_processes:
        print(f"Num processes defined by user: {num_processes}")
    else:
        if isinstance(gb_per_process, float) and gb_per_process > 0:
            max_processes_mem = max(int(mem / float(gb_per_process)), 1)
            if max_processes_mem < cpu_count:
                print(f"Num processes limited by memory: {max_processes_mem}")
                num_processes = max_processes_mem
            else:
                print(f"Num processes limited by cpu count: {cpu_count}")
                num_processes = cpu_count
        else:
            print(f"Num processes limited by cpu count: {cpu_count}")
            num_processes = cpu_count

    return num_processes


def mp_run(func, argument_list, mp_conf: Dict[str, float]) -> None:
    print(f"Running {func.__qualname__}. This may take some time...")
    start = datetime.now()
    print(f"STARTED AT: {start}")
    with Pool(get_num_processes(mp_conf)) as pool:
        jobs = [pool.apply_async(func, args=args) for args in argument_list]
        for job in jobs:
            job.get()
    end = datetime.now()
    print(f"FINISHED AT: {end}")
    print(f"DURATION: {end - start}")
