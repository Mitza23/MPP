#include <vector>

using namespace std;

vector<vector<int>> multiplySequential(const vector<vector<int>>& a, const vector<vector<int>>& b) {
    vector<vector<int>> result(a.size(), vector<int>(b[0].size(), 0));

    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b[0].size(); ++j) {
            for (size_t k = 0; k < b.size(); ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return result;
}