#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <chrono>
#include "data/dataUtil.cpp"
#include "solutions/l1.cpp"
#include "solutions/l2.cpp"

using namespace std;


void generate_test_data() {
    vector<vector<int>> a = {{1, 2, 1},
                             {2, 1, 2},
                             {1, 2, 1}};
    vector<vector<int>> b = {{2, 1, 2},
                             {1, 2, 1},
                             {2, 1, 2}};
    writeMatrixToFile(a, "a.pkl");
    writeMatrixToFile(b, "b.pkl");
}

void test_multiply_sequential() {
    generate_test_data();
    vector<vector<int>> a = readMatrixFromFile("a.pkl", 3, 3);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 3, 3);
    vector<vector<int>> result = multiplySequential(a, b);

    for (const auto &row: result) {
        for (const auto &elem: row) {
            cout << elem << ' ';
        }
        cout << '\n';
    }
}


void test_multiply_multithreading_row(const int nrThreads) {
    auto start = chrono::high_resolution_clock::now();  // Start timing here

    generateData(1024, 1024, 512, "a.pkl", "b.pkl");
    vector<vector<int>> a = readMatrixFromFile("a.pkl", 1024, 512);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 512, 1024);
    vector<vector<int>> expected = multiplySequential(a, b);
    vector<vector<int>> c = multiplyMultithreadingRows(a, b, nrThreads);

    assert(c == expected);

    auto end = chrono::high_resolution_clock::now();  // End timing here
    chrono::duration<double> elapsed = end - start;

    cout << "Test passed for " << nrThreads << "threads in " << elapsed.count() << " seconds.\n";
}

void test_multiply_multithreading_block() {
    generate_test_data();
    vector<vector<int>> a = readMatrixFromFile("a.pkl", 3, 3);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 3, 3);
    vector<vector<int>> expected = multiplySequential(a, b);
    vector<vector<int>> c = multiplyMultithreadingBlocks(a, b, 2);

    assert(c == expected);
    cout << "Test passed.\n";
}

//////////////////////////////////////////////////

double benchmarkMultiplyMultithreadingRow(const int nrThreads) {
//    auto start = chrono::high_resolution_clock::now();

    vector<vector<int>> a = readMatrixFromFile("a.pkl", 1024, 512);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 512, 1024);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<int>> c = multiplyMultithreadingRows(a, b, nrThreads);
    auto end = chrono::high_resolution_clock::now();

    writeMatrixToFile(c, "c.pkl");

//    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    return elapsed.count();
}

double benchmarkMultiplyMultithreadingBlock(const int nrThreads) {
//    auto start = chrono::high_resolution_clock::now();

    vector<vector<int>> a = readMatrixFromFile("a.pkl", 1024, 512);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 512, 1024);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<int>> c = multiplyMultithreadingBlocks(a, b, nrThreads);
    auto end = chrono::high_resolution_clock::now();

    writeMatrixToFile(c, "c.pkl");

//    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    return elapsed.count();
}

double runBenchmark(double (*benchmarkFunction)(const int), const int nrThreads) {
    double totalTime = 0;
    for (int i = 0; i < 10; i++) {
        cout << ". . . ";
        totalTime += benchmarkFunction(nrThreads);
    }
    cout << endl;
    return totalTime / 10;
}


void assessPerformance(double (*benchmarkFunction)(const int)) {
    generateData(1024, 1024, 512, "a.pkl", "b.pkl");

    cout << "Starting performance assessment for 4 threads.\n";
    double time4 = runBenchmark(benchmarkFunction, 4);
    cout << "Performance assessment for 4 threads: " << time4 << " seconds.\n\n";

    cout << "Starting performance assessment for 8 threads.\n";
    double time8 = runBenchmark(benchmarkFunction, 8);
    cout << "Performance assessment for 8 threads: " << time8 << " seconds.\n\n";

    cout << "Starting performance assessment for 16 threads.\n";
    double time16 = runBenchmark(benchmarkFunction, 16);
    cout << "Performance assessment for 16 threads: " << time16 << " seconds.\n\n";

    cout << "Starting performance assessment for 32 threads.\n";
    double time32 = runBenchmark(benchmarkFunction, 32);
    cout << "Performance assessment for 32 threads: " << time32 << " seconds.\n\n";

}

double benchmarkMultiplySequential() {
//    auto start = chrono::high_resolution_clock::now();

    vector<vector<int>> a = readMatrixFromFile("a.pkl", 1024, 512);
    vector<vector<int>> b = readMatrixFromFile("b.pkl", 512, 1024);

    auto start = chrono::high_resolution_clock::now();
    vector<vector<int>> c = multiplySequential(a, b);
    auto end = chrono::high_resolution_clock::now();

    writeMatrixToFile(c, "c.pkl");

//    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    return elapsed.count();
}

void assessSequentialPerformance() {
    generateData(1024, 1024, 512, "a.pkl", "b.pkl");

    double totalTime = 0;
    for (int i = 0; i < 10; i++) {
        cout << ". . . ";
        totalTime += benchmarkMultiplySequential();
    }
    cout << endl;

    cout << "Performance assessment for sequential multiplication: " << totalTime / 10 << " seconds.\n\n";
}

int main() {
    cout << "Testing MultiplySequential \n\n";
    assessSequentialPerformance();

    cout << "\n\n\n\n";

    cout << "Testing MultiplyMultithreadingRow \n\n";
    assessPerformance(benchmarkMultiplyMultithreadingRow);

    cout << "\n\n\n\n";

    cout << "Testing MultiplyMultithreadingBlock \n\n";
    assessPerformance(benchmarkMultiplyMultithreadingBlock);
    return 0;

}