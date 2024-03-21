#include <random>
#include <fstream>
#include <vector>

using namespace std;

vector<vector<int>> generateMatrix(int n, int m) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(-10000, 10000);

    vector<vector<int>> matrix(n, vector<int>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            matrix[i][j] = dis(gen);
        }
    }
    return matrix;
}

#include <fstream>
#include <vector>

void writeMatrixToFile(const vector<vector<int>>& matrix, const string& filename) {
    ofstream ofs(filename, ios::binary);
    for (const auto& row : matrix) {
        ofs.write(reinterpret_cast<const char*>(row.data()), row.size()*sizeof(int));
    }
}

vector<vector<int>> readMatrixFromFile(const string& filename, int n, int m) {
    ifstream ifs(filename, ios::binary);
    vector<vector<int>> matrix(n, vector<int>(m));

    for (auto& row : matrix) {
        ifs.read(reinterpret_cast<char*>(row.data()), row.size()*sizeof(int));
    }

    return matrix;
}

void generateMatrixToFile(int n, int m, const string &filename) {
    vector<vector<int>> matrix = generateMatrix(n, m);
    writeMatrixToFile(matrix, filename);
}

void generateData(int n, int m, int k, const string &filename_a, const string &filename_b) {
    generateMatrixToFile(n, k, filename_a);
    generateMatrixToFile(k, m, filename_b);
}