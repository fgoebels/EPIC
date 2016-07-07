#!/usr/bin/env python

import multiprocessing as mp
import numpy as np
import time

def foo_pool(x):
    time.sleep(2)
    return x*x

def apply_async_with_callback():
    arr = mp.Array([[1,2], [3,4]])
    pool = mp.Pool()
    result_list = []
    tasks = range(10)
    results = pool.map_async(foo_pool, tasks)
    results.wait()
    print(results.get())

if __name__ == '__main__':
    apply_async_with_callback()
