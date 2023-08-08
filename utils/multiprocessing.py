import multiprocessing as mp
from datetime import datetime
from multiprocessing.pool import Pool
from typing import Dict

from psutil import virtual_memory


def _get_pool(mp_conf: Dict[str, float]):
    # num_processes : the number of processes, if None the calculation is done based on the gb_per_process limit
    # gb_per_process : the expected amount of ram one process needs, if None there is no limit
    num_processes = mp_conf["num_processes"]
    gb_per_process = mp_conf["gb_per_process"]

    cpu_count = mp.cpu_count()
    print("Number of cpu : ", cpu_count)
    mem = virtual_memory().total * 1e-9  # System memory in gb
    print("System memory: ", mem)

    if num_processes:
        print('Num processes defined by user')
        num_processes = num_processes
    else:
        if gb_per_process:
            max_precesses_mem = max(int(mem / float(gb_per_process)), 1)
            if max_precesses_mem < cpu_count:
                print('Num processes limited by memory')
                num_processes = max_precesses_mem
            else:
                print('Num processes limited by cpu count')
                num_processes = cpu_count
        else:
            print('Num processes limited by cpu count')
            num_processes = cpu_count

    print("Num processes:", num_processes)

    pool = Pool(processes=num_processes)

    return pool


def mp_run(func, argument_list, mp_conf: Dict[str, float]) -> None:
    print(f"Running {func.__qualname__} in {mp_conf} processes. This may take some time...")
    start = datetime.now()
    print(f"STARTED AT: {start}")
    with _get_pool(mp_conf) as pool:
        pool.starmap(func=func, iterable=argument_list)
    end = datetime.now()
    print(f"FINISHED AT: {end}")
    print(f"DURATION: {start - end}")
