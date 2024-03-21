import random
import pickle


def generate_matrix(n, m):
    return [[random.randint(-10000, 10000) for _ in range(m)] for _ in range(n)]

def read_matrix_from_file(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)
def write_matrix_to_file(matrix, filename):
    with open(filename, 'wb') as f:
        pickle.dump(matrix, f)

def generate_matrix_to_file(n, m, filename):
    matrix = generate_matrix(n, m)
    write_matrix_to_file(matrix, filename)

def generate_data(n, m, k, filename_a, filename_b):
    generate_matrix_to_file(n, k, filename_a)
    generate_matrix_to_file(k, m, filename_b)