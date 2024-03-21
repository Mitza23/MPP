def multiply_sequential(a, b):
    result = []
    for row in a:
        result_row = []
        for col in zip(*b):
            result_row.append(sum(x * y for x, y in zip(row, col)))
        result.append(result_row)
    return result
