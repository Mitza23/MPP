#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

// Define constants
const double Ra = 500.0;          // Rayleigh number
const int N1 = 200;               // Number of nodes in x-direction
const int N2 = 200;               // Number of nodes in y-direction
const int N3 = 100;               // Number of nodes in z-direction
const double Ax = 2.0;            // Length of the cavity in x-direction
const double Ay = 2.0;            // Length of the cavity in y-direction
const double Az = 1.0;            // Length of the cavity in z-direction
const int maxIterations = 10000; // Maximum number of iterations
const int step = 500;             // Number of steps to output data
const double tolerance = 1e-9;    // Convergence tolerance

void recordParameters(const vector<vector<vector<double>>> &u,
                      const vector<vector<vector<double>>> &v,
                      const vector<vector<vector<double>>> &t,
                      int N1,
                      int N2,
                      int N3,
                      int step) {
    ofstream fileOut("iteration_" + to_string(step) + ".txt");
    fileOut << N1 << " " << N2 << " " << N3 << "\n";
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            for (int k = 0; k < N3; ++k) {
                fileOut << i << " " << j << " " << k << " " << u[i][j][k] << " " << v[i][j][k] << " "
                        << t[i][j][k]
                        << "\n";
            }
        }
    }
    fileOut.close();
    cout << "recorded progress at iteration " << step << endl;
}

double computeJacobiForVelocityU(const vector<vector<vector<double>>> &velocityVectorU,
                                 const vector<vector<vector<double>>> &temperatureVector,
                                 const double cellSize,
                                 int i,
                                 int j,
                                 int k) {
    return Ra * cellSize / 12.0 * (temperatureVector[i][j + 1][k] - temperatureVector[i][j - 1][k]) +
           1.0 / 6.0 * (velocityVectorU[i - 1][j][k] + velocityVectorU[i + 1][j][k] +
                        velocityVectorU[i][j - 1][k] + velocityVectorU[i][j + 1][k] +
                        velocityVectorU[i][j][k - 1] + velocityVectorU[i][j][k + 1]);
}


double computeJacobiForVelocityV(const vector<vector<vector<double>>> &velocityVectorV,
                                 const vector<vector<vector<double>>> &temperatureVector,
                                 const double cellSize,
                                 int i,
                                 int j,
                                 int k) {
    return Ra * cellSize / 12.0 * (temperatureVector[i + 1][j][k] - temperatureVector[i - 1][j][k]) +
           1.0 / 6.0 * (velocityVectorV[i - 1][j][k] + velocityVectorV[i + 1][j][k] +
                        velocityVectorV[i][j - 1][k] + velocityVectorV[i][j + 1][k] +
                        velocityVectorV[i][j][k - 1] + velocityVectorV[i][j][k + 1]);
}

double computeJacobiForTemperature(const vector<vector<vector<double>>> &temperatureVector,
                                   const vector<vector<vector<double>>> &velocityVectorU,
                                   const vector<vector<vector<double>>> &velocityVectorV,
                                   const double cellSize,
                                   int i,
                                   int j,
                                   int k) {
    return 1.0 / 6.0 *
           (temperatureVector[i - 1][j][k] + temperatureVector[i + 1][j][k] + temperatureVector[i][j - 1][k] +
            temperatureVector[i][j + 1][k] +
            temperatureVector[i][j][k - 1] + temperatureVector[i][j][k + 1] + cellSize * cellSize -
            1.0 / 4.0 * (((velocityVectorU[i][j][k + 1] - velocityVectorU[i][j][k - 1]) *
                          (temperatureVector[i][j + 1][k] - temperatureVector[i][j - 1][k]) -
                          (velocityVectorU[i][j + 1][k] - velocityVectorU[i][j - 1][k]) *
                          (temperatureVector[i][j][k + 1] - temperatureVector[i][j][k - 1])) +
                         ((velocityVectorV[i + 1][j][k] - velocityVectorV[i - 1][j][k]) *
                          (temperatureVector[i][j][k + 1] - temperatureVector[i][j][k - 1]) -
                          (velocityVectorV[i][j][k + 1] - velocityVectorV[i][j][k - 1]) *
                          (temperatureVector[i + 1][j][k] - temperatureVector[i - 1][j][k]))));
}

int main() {
    double h = Ax / (N1 - 1); // Cubic grid cell size

    // Initialize arrays
    vector<vector<vector<double>>> u(N1, vector<vector<double>>(N2, vector<double>(N3, 0.0)));
    vector<vector<vector<double>>> u_new = u;
    vector<vector<vector<double>>> v(N1, vector<vector<double>>(N2, vector<double>(N3, 0.0)));
    vector<vector<vector<double>>> v_new = v;
    vector<vector<vector<double>>> t(N1, vector<vector<double>>(N2, vector<double>(N3, 0.0)));
    vector<vector<vector<double>>> t_new = t;

    int iteration = 0;
    bool done = false;

    while (!done && iteration < maxIterations) {
        iteration++;
        double err_u = 0.0, err_v = 0.0, err_T = 0.0;

        for (int i = 1; i < N1 - 1; ++i) {
            for (int j = 1; j < N2 - 1; ++j) {
                for (int k = 1; k < N3 - 1; ++k) {
                    u_new[i][j][k] = computeJacobiForVelocityU(u, t, h, i, j, k);
                    v_new[i][j][k] = computeJacobiForVelocityV(v, t, h, i, j, k);
                    t_new[i][j][k] = computeJacobiForTemperature(t, u, v, h, i, j, k);

                    // Error evaluation
                    err_u = max(err_u, abs(u_new[i][j][k] - u[i][j][k]));
                    err_v = max(err_v, abs(v_new[i][j][k] - v[i][j][k]));
                    err_T = max(err_T, abs(t_new[i][j][k] - t[i][j][k]));
                }
            }
        }

        // Apply adiabatic condition
        for (int i = 0; i < N1; ++i) {
            for (int j = 0; j < N2; ++j) {
                t_new[i][j][0] = t_new[i][j][1];
                t_new[i][j][N3 - 1] = t_new[i][j][N3 - 2];
            }
        }

        if (err_u < tolerance && err_v < tolerance && err_T < tolerance) {
            done = true;
            recordParameters(u, v, t, N1, N2, N3, iteration);
        } else {
            u = u_new;
            v = v_new;
            t = t_new;
            // Record parameters for visualization
            if (iteration % step == 0) {
                recordParameters(u, v, t, N1, N2, N3, iteration);
            }
        }
    }

    if (iteration >= maxIterations) {
        cout << "Maximum iterations limit reached without desired convergence" << endl;
    } else {
        cout << "Converged in " << iteration << " iterations" << endl;
    }

    return 0;
}
