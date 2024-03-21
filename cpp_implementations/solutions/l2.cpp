#include <vector>
#include <thread>
#include <cmath>

using namespace std;

void multiply_rows(const vector<vector<int>>& a, const vector<vector<int>>& b, vector<vector<int>>& result, int start, int end) {
    for (int i = start; i < end; ++i) {
        for (size_t j = 0; j < b[0].size(); ++j) {
            for (size_t k = 0; k < b.size(); ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

vector<vector<int>> multiplyMultithreadingRows(const vector<vector<int>>& a, const vector<vector<int>>& b, int num_threads) {
    vector<vector<int>> result(a.size(), vector<int>(b[0].size(), 0));
    vector<thread> threads;

    int rows_per_thread = a.size() / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        int start = i * rows_per_thread;
        int end = (i == num_threads - 1) ? a.size() : start + rows_per_thread;
        threads.emplace_back(multiply_rows, cref(a), cref(b), ref(result), start, end);
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return result;
}


void multiply_block(const vector<vector<int>>& a, const vector<vector<int>>& b, vector<vector<int>>& result, int start_row, int end_row, int start_col, int end_col) {
    for (int i = start_row; i < end_row; ++i) {
        for (int j = start_col; j < end_col; ++j) {
            for (size_t k = 0; k < b.size(); ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

vector<vector<int>> multiplyMultithreadingBlocks(const vector<vector<int>>& a, const vector<vector<int>>& b, int num_threads) {
    vector<vector<int>> result(a.size(), vector<int>(b[0].size(), 0));
    vector<thread> threads;

    int block_size = sqrt(num_threads);
    int rows_per_block = a.size() / block_size;
    int cols_per_block = b[0].size() / block_size;

    for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < block_size; ++j) {
            int start_row = i * rows_per_block;
            int end_row = (i == block_size - 1) ? a.size() : start_row + rows_per_block;
            int start_col = j * cols_per_block;
            int end_col = (j == block_size - 1) ? b[0].size() : start_col + cols_per_block;
            threads.emplace_back(multiply_block, cref(a), cref(b), ref(result), start_row, end_row, start_col, end_col);
        }
    }

    for (auto& thread : threads) {
        thread.join();
    }

    return result;
}