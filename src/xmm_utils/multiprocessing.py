from collections.abc import Generator
from datetime import datetime, timedelta
from functools import partial
from multiprocessing.pool import Pool

from loguru import logger


def handle_error(error):
    logger.exception(error)


def mp_run(func, kwds_generator: Generator, num_processes: int, debug: bool = False) -> tuple[list, timedelta]:
    start = datetime.now()
    with Pool(processes=num_processes) as pool:
        mp_apply = pool.apply if debug else partial(pool.apply_async, error_callback=handle_error)
        apply_results = [mp_apply(func, kwds=kwd) for kwd in kwds_generator]
        pool.close()
        pool.join()
    end = datetime.now()
    results = []
    print(f"apply_results: {apply_results}")
    flattened_apply_results = [
        item for sublist in apply_results for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    print(f"flattened_apply_results: {flattened_apply_results}")
    for apply_result in flattened_apply_results:
        if apply_result is not None:
            results.append(apply_result)
    return results, end - start
