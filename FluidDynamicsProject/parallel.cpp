#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mpi.h>

using namespace std;

// Define constants
// Define constants
const double Ra = 100.0;          // Rayleigh number
const int N = 10;                 // Axis length
const double Ax = 1.0;            // Length of the cavity in x-direction
const double Ay = 1.0;            // Length of the cavity in y-direction
const double Az = 1.0;            // Length of the cavity in z-direction
const int maxIterations = 1000; // Maximum number of iterations
const int step = 50;             // Number of steps to output data
const double tolerance = 1e-5;    // Convergence tolerance

void recordParameters(const std::vector<double> &u,
                      const std::vector<double> &v,
                      const std::vector<double> &t,
                      int N1,
                      int N2,
                      int N3,
                      int step) {
    std::ofstream fileOut("parallel_data_step_" + std::to_string(step) + ".txt");
    fileOut << N1 << " " << N2 << " " << N3 << "\n";
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            for (int k = 0; k < N3; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                fileOut << i << " " << j << " " << k << " " << u[index] << " " << v[index] << " " << t[index] << "\n";
            }
        }
    }
    fileOut.close();
    cout << "recorded progress at iteration " << step << endl;
}

void computeVelocityVectorV(int start,
                            int end,
                            int N2,
                            int N3,
                            double cellSize,
                            double Ra,
                            const std::vector<double> &velocityVectorV,
                            const std::vector<double> &temperatureVector,
                            std::vector<double> &v_new,
                            double &err_v) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                v_new[index] = Ra * cellSize / 12.0 * (temperatureVector[(i + 1) * N2 * N3 + j * N3 + k] -
                                                       temperatureVector[(i - 1) * N2 * N3 + j * N3 + k]) +
                               1.0 / 6.0 * (velocityVectorV[(i - 1) * N2 * N3 + j * N3 + k] +
                                            velocityVectorV[(i + 1) * N2 * N3 + j * N3 + k] +
                                            velocityVectorV[i * N2 * N3 + (j - 1) * N3 + k] +
                                            velocityVectorV[i * N2 * N3 + (j + 1) * N3 + k] +
                                            velocityVectorV[i * N2 * N3 + j * N3 + (k - 1)] +
                                            velocityVectorV[i * N2 * N3 + j * N3 + (k + 1)]);
                err_v = std::max(err_v, std::abs(v_new[index] - velocityVectorV[index]));
            }
        }
    }
}

void computeVelocityVectorU(int start,
                            int end,
                            int N2,
                            int N3,
                            double cellSize,
                            double Ra,
                            const std::vector<double> &velocityVectorU,
                            const std::vector<double> &temperatureVector,
                            std::vector<double> &u_new,
                            double &err_u) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                u_new[index] = Ra * cellSize / 12.0 * (temperatureVector[i * N2 * N3 + (j + 1) * N3 + k] -
                                                       temperatureVector[i * N2 * N3 + (j - 1) * N3 + k]) +
                               1.0 / 6.0 * (velocityVectorU[(i - 1) * N2 * N3 + j * N3 + k] +
                                            velocityVectorU[(i + 1) * N2 * N3 + j * N3 + k] +
                                            velocityVectorU[i * N2 * N3 + (j - 1) * N3 + k] +
                                            velocityVectorU[i * N2 * N3 + (j + 1) * N3 + k] +
                                            velocityVectorU[i * N2 * N3 + j * N3 + (k - 1)] +
                                            velocityVectorU[i * N2 * N3 + j * N3 + (k + 1)]);
                err_u = std::max(err_u, std::abs(u_new[index] - velocityVectorU[index]));
            }
        }
    }
}

void computeTemperatureVector(int start,
                              int end, int N2, int N3, double h,
                              const std::vector<double> &velocityVectorU,
                              const std::vector<double> &velocityVectorV,
                              const std::vector<double> &temperatureVector,
                              std::vector<double> &t_new,
                              double &err_T) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                t_new[index] = 1.0 / 6.0 * (temperatureVector[(i - 1) * N2 * N3 + j * N3 + k] +
                                            temperatureVector[(i + 1) * N2 * N3 + j * N3 + k] +
                                            temperatureVector[i * N2 * N3 + (j - 1) * N3 + k] +
                                            temperatureVector[i * N2 * N3 + (j + 1) * N3 + k] +
                                            temperatureVector[i * N2 * N3 + j * N3 + (k - 1)] +
                                            temperatureVector[i * N2 * N3 + j * N3 + (k + 1)] +
                                            h * h - 1.0 / 4.0 * (
                        (
                                (velocityVectorU[i * N2 * N3 + j * N3 + (k + 1)] -
                                 velocityVectorU[i * N2 * N3 + j * N3 + (k - 1)]) *
                                (temperatureVector[i * N2 * N3 + (j + 1) * N3 + k] -
                                 temperatureVector[i * N2 * N3 + (j - 1) * N3 + k]) -
                                (velocityVectorU[i * N2 * N3 + (j + 1) * N3 + k] -
                                 velocityVectorU[i * N2 * N3 + (j - 1) * N3 + k]) *
                                (temperatureVector[i * N2 * N3 + j * N3 + (k + 1)] -
                                 temperatureVector[i * N2 * N3 + j * N3 + (k - 1)])
                        ) +
                        (
                                (velocityVectorV[(i + 1) * N2 * N3 + j * N3 + k] -
                                 velocityVectorV[(i - 1) * N2 * N3 + j * N3 + k]) *
                                (temperatureVector[i * N2 * N3 + j * N3 + (k + 1)] -
                                 temperatureVector[i * N2 * N3 + j * N3 + (k - 1)]) -
                                (velocityVectorV[i * N2 * N3 + j * N3 + (k + 1)] -
                                 velocityVectorV[i * N2 * N3 + j * N3 + (k - 1)]) *
                                (temperatureVector[(i + 1) * N2 * N3 + j * N3 + k] -
                                 temperatureVector[(i - 1) * N2 * N3 + j * N3 + k])
                        )));
                err_T = std::max(err_T, std::abs(t_new[index] - temperatureVector[index]));
            }
        }
    }
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int N1 = N * Ax;          // Number of nodes in x-direction
    int N2 = N * Ay;          // Number of nodes in y-direction
    int N3 = N;               // Number of nodes in z-direction
    double h = Ax / (N1 - 1); // Cubic grid cell size

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Compute chunk indexes
    int chunk_size = (N1 + size - 1) / size;
    int start = rank * chunk_size;
    int end = min(start + chunk_size, N1);

    // Initialize arrays
    vector<double> u(N1 * N2 * N3, 0.0);
    vector<double> u_new = u;
    vector<double> v(N1 * N2 * N3, 0.0);
    vector<double> v_new = v;
    vector<double> t(N1 * N2 * N3, 0.0);
    vector<double> t_new = t;

    int iteration = 0;
    bool stop = false;
    cout << "Process " << rank << " starting" << endl;
    while (!stop && iteration < maxIterations) {
        iteration++;
        double err_u = 0.0, err_v = 0.0, err_T = 0.0;

        thread threadU(computeVelocityVectorU, start, end, N2, N3, h, Ra, ref(u), ref(t), ref(u_new), ref(err_u));
        thread threadV(computeVelocityVectorV, start, end, N2, N3, h, Ra, ref(v), ref(t), ref(v_new), ref(err_v));
        thread threadT(computeTemperatureVector, start, end, N2, N3, h, ref(u), ref(v), ref(t), ref(t_new), ref(err_T));

        threadU.join();
        threadV.join();
        threadT.join();

        // Parent collecting data from all child processes
        if (rank == 0) {
            for (int p = 1; p < size; ++p) {
                int p_start = p * chunk_size;
                int p_end = min(p_start + chunk_size, N1);

                MPI_Recv(&u_new[p_start * N2 * N3], (p_end - p_start) * N2 * N3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&v_new[p_start * N2 * N3], (p_end - p_start) * N2 * N3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&t_new[p_start * N2 * N3], (p_end - p_start) * N2 * N3, MPI_DOUBLE, p, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            }

            // Apply adiabatic condition
            for (int i = 0; i < N1; ++i) {
                for (int j = 0; j < N2; ++j) {
                    t_new[i * N2 * N3 + j * N3 + 0] = t_new[i * N2 * N3 + j * N3 + 1];
                    t_new[i * N2 * N3 + j * N3 + (N3 - 1)] = t_new[i * N2 * N3 + j * N3 + (N3 - 2)];
                }
            }
        } else {
            // Children return
            MPI_Send(&u_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&v_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&t_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        // Parent broadcasting result to all children for next iteration
        MPI_Bcast(&u_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&v_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&t_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (iteration % step == 0 && rank == 0) {
            recordParameters(u, v, t, N1, N2, N3, iteration);
        }

        if (err_u < tolerance && err_v < tolerance && err_T < tolerance) {
            stop = true;
            if (rank == 0)
                recordParameters(u, v, t, N1, N2, N3, iteration);
        }

        u = u_new;
        v = v_new;
        t = t_new;

        MPI_Allreduce(MPI_IN_PLACE, &stop, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        if (iteration >= maxIterations) {
            cout << "Maximum iteration limit reached without desired convergence" << endl;
        } else {
            cout << "Converged in " << iteration << " iteration" << endl;
        }
    }

    MPI_Finalize();
    return 0;
}