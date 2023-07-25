import multiprocessing as mp
from typing import Iterable

from psutil import virtual_memory


def get_multiprocessing_pool(num_processes=None, gb_per_process=None):
    # num_processes : the number of processes, if None the calculation is done based on the gb_per_process limit
    # gb_per_process : the expected amount of ram one process needs, if None there is no limit

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

    pool = mp.Pool(processes=num_processes)

    return pool


def run_apply_async_multiprocessing(func, argument_list, num_processes=None, gb_per_process=None):
    # From: https://leimao.github.io/blog/Python-tqdm-Multiprocessing/
    #
    # pool = get_multiprocessing_pool(num_processes=num_processes, gb_per_process=gb_per_process)
    #
    # jobs = [pool.apply_async(func=func, args=(*argument,)) if isinstance(argument, tuple)
    #         else pool.apply_async(func=func, args=(argument,)) for argument in argument_list]
    #
    # pool.close()
    # result_list_tqdm = []
    # for job in tqdm(jobs):
    #     result_list_tqdm.append(job.get())
    #
    # return result_list_tqdm
    #
    results = []
    with mp.Pool(processes=num_processes) as pool:
        jobs = [pool.apply_async(func=func, args=(*argument,)) if isinstance(argument, tuple)
                else pool.apply_async(func=func, args=(argument,)) for argument in argument_list]
        for job in jobs:
            res = job.get()
            if isinstance(res, Iterable):
                results.extend(res)
            else:
                results.append(res)
        pool.close()
        pool.join()
    return results
