import threading


def multiply_split_row(a, b, result, start, end):
    rows_a, cols_a = len(a), len(a[0])
    rows_b, cols_b = len(b), len(b[0])
    for i in range(start, end):
        for j in range(cols_b):
            result[i][j] = sum(a[i][k] * b[k][j] for k in range(len(b)))


def multiply_multithreading_row(a, b, num_threads):
    rows_a, cols_a = len(a), len(a[0])
    rows_b, cols_b = len(b), len(b[0])

    if cols_a != rows_b:
        return "The matrices are not compatible for multiplication."

    result = [[0 for _ in range(cols_b)] for _ in range(rows_a)]
    threads = []

    rows_per_thread = rows_a // num_threads
    for i in range(num_threads):
        start = i * rows_per_thread
        # The last thread should get the remaining rows if the number of rows is not divisible by num_threads
        end = rows_a if i == num_threads - 1 else start + rows_per_thread
        thread = threading.Thread(target=multiply_split_row, args=(a, b, result, start, end))
        threads.append(thread)
        thread.start()

    for thread in threads:
        thread.join()

    return result


import math


def multiply_block(a, b, result, start_row, end_row, start_col, end_col):
    for i in range(start_row, end_row):
        for j in range(start_col, end_col):
            result[i][j] = sum(a[i][k] * b[k][j] for k in range(len(b)))


def multiply_multithreading_block(a, b, num_threads):
    rows_a, cols_a = len(a), len(a[0])
    rows_b, cols_b = len(b), len(b[0])

    if cols_a != rows_b:
        return "The matrices are not compatible for multiplication."

    result = [[0 for _ in range(cols_b)] for _ in range(rows_a)]
    threads = []

    block_size = int(math.sqrt(num_threads))
    rows_per_block = rows_a // block_size
    cols_per_block = cols_b // block_size

    for i in range(block_size):
        for j in range(block_size):
            start_row = i * rows_per_block
            end_row = rows_a if i == block_size - 1 else start_row + rows_per_block
            start_col = j * cols_per_block
            end_col = cols_b if j == block_size - 1 else start_col + cols_per_block
            thread = threading.Thread(target=multiply_block,
                                      args=(a, b, result, start_row, end_row, start_col, end_col))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()

    return result
