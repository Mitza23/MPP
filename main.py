from setup.data_util import *
from solutions.l1 import *
from solutions.l2 import *


def test_read():
    generate_data(3, 3, 3, 'a.pkl', 'b.pkl')
    print(read_matrix_from_file('a.pkl'))
    print(read_matrix_from_file('b.pkl'))


def generate_test_data():
    a = [[1, 2, 1], [2, 1, 2], [1, 2, 1]]
    b = [[2, 1, 2], [1, 2, 1], [2, 1, 2]]
    write_matrix_to_file(a, 'a.pkl')
    write_matrix_to_file(b, 'b.pkl')


def test_multiply_sequential():
    generate_test_data()
    a = read_matrix_from_file('a.pkl')
    b = read_matrix_from_file('b.pkl')
    print(multiply_sequential(a, b))


def test_multiply_multithreading_row():
    generate_data(100, 100, 64, 'a.pkl', 'b.pkl')
    a = read_matrix_from_file('a.pkl')
    b = read_matrix_from_file('b.pkl')
    expected = multiply_sequential(a, b)
    c = multiply_multithreading_row(a, b, 16)
    print(c == expected)

def test_multiply_multithreading_block():
    generate_data(100, 100, 64, 'a.pkl', 'b.pkl')
    a = read_matrix_from_file('a.pkl')
    b = read_matrix_from_file('b.pkl')
    expected = multiply_sequential(a, b)
    c = multiply_multithreading_block(a, b, 10)
    print(c == expected)


def solve(filename_a, filename_b, filename_result, multiplication_method):
    a = read_matrix_from_file(filename_a)
    b = read_matrix_from_file(filename_b)
    write_matrix_to_file(multiplication_method(a, b), filename_result)


if __name__ == '__main__':
    test_multiply_multithreading_block()
