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
    for apply_result in apply_results:
        result = apply_result if debug else apply_result.get()
        if result is not None:
            if isinstance(result, list):
                results.extend(result)
            else:
                results.append(result)
    return results, end - start
