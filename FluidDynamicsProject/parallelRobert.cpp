#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mpi.h>

// Define constants
const double Ra = 100.0;          // Rayleigh number
const int N = 10;                 // Axis length
const double Ax = 1.0;            // Length of the cavity in x-direction
const double Ay = 1.0;            // Length of the cavity in y-direction
const double Az = 1.0;            // Length of the cavity in z-direction
const int maxIterations = 10000; // Maximum number of iterations
const int step = 100;             // Number of steps to output data

void writeData(const std::vector<double> &u,
               const std::vector<double> &v,
               const std::vector<double> &t,
               int N1,
               int N2,
               int N3,
               int step) {
    std::ofstream file("parallel_data_step_" + std::to_string(step) + ".txt");
    file << N1 << " " << N2 << " " << N3 << "\n";
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            for (int k = 0; k < N3; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                file << i << " " << j << " " << k << " " << u[index] << " " << v[index] << " " << t[index] << "\n";
            }
        }
    }
    file.close();
    std::cout << "Wrote data at step " << step << std::endl;
}

void computeVelocityVectorU(int start,
                            int end,
                            int N2,
                            int N3,
                            double h,
                            double Ra,
                            const std::vector<double> &u,
                            const std::vector<double> &t,
                            std::vector<double> &u_new,
                            double &err_u) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                u_new[index] = Ra * h / 12.0 * (t[i * N2 * N3 + (j + 1) * N3 + k] - t[i * N2 * N3 + (j - 1) * N3 + k]) +
                               1.0 / 6.0 * (u[(i - 1) * N2 * N3 + j * N3 + k] + u[(i + 1) * N2 * N3 + j * N3 + k] +
                                            u[i * N2 * N3 + (j - 1) * N3 + k] + u[i * N2 * N3 + (j + 1) * N3 + k] +
                                            u[i * N2 * N3 + j * N3 + (k - 1)] + u[i * N2 * N3 + j * N3 + (k + 1)]);
                err_u = std::max(err_u, std::abs(u_new[index] - u[index]));
            }
        }
    }
}

void computeV(int start,
              int end,
              int N2,
              int N3,
              double h,
              double Ra,
              const std::vector<double> &v,
              const std::vector<double> &t,
              std::vector<double> &v_new,
              double &err_v) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                v_new[index] = Ra * h / 12.0 * (t[(i + 1) * N2 * N3 + j * N3 + k] - t[(i - 1) * N2 * N3 + j * N3 + k]) +
                               1.0 / 6.0 * (v[(i - 1) * N2 * N3 + j * N3 + k] + v[(i + 1) * N2 * N3 + j * N3 + k] +
                                            v[i * N2 * N3 + (j - 1) * N3 + k] + v[i * N2 * N3 + (j + 1) * N3 + k] +
                                            v[i * N2 * N3 + j * N3 + (k - 1)] + v[i * N2 * N3 + j * N3 + (k + 1)]);
                err_v = std::max(err_v, std::abs(v_new[index] - v[index]));
            }
        }
    }
}

void computeTemperatureVector(int start,
                              int end, int N2, int N3, double h,
                              const std::vector<double> &u,
                              const std::vector<double> &v,
                              const std::vector<double> &t,
                              std::vector<double> &t_new,
                              double &err_T) {
    for (int i = start + (start == 0 ? 1 : 0); i < end - (end == N2 ? 1 : 0); ++i) {
        for (int j = 1; j < N2 - 1; ++j) {
            for (int k = 1; k < N3 - 1; ++k) {
                int index = i * N2 * N3 + j * N3 + k;
                t_new[index] = 1.0 / 6.0 * (t[(i - 1) * N2 * N3 + j * N3 + k] + t[(i + 1) * N2 * N3 + j * N3 + k] +
                                            t[i * N2 * N3 + (j - 1) * N3 + k] + t[i * N2 * N3 + (j + 1) * N3 + k] +
                                            t[i * N2 * N3 + j * N3 + (k - 1)] + t[i * N2 * N3 + j * N3 + (k + 1)] +
                                            h * h - 1.0 / 4.0 * (((u[i * N2 * N3 + j * N3 + (k + 1)] -
                                                                   u[i * N2 * N3 + j * N3 + (k - 1)]) *
                                                                  (t[i * N2 * N3 + (j + 1) * N3 + k] -
                                                                   t[i * N2 * N3 + (j - 1) * N3 + k]) -
                                                                  (u[i * N2 * N3 + (j + 1) * N3 + k] -
                                                                   u[i * N2 * N3 + (j - 1) * N3 + k]) *
                                                                  (t[i * N2 * N3 + j * N3 + (k + 1)] -
                                                                   t[i * N2 * N3 + j * N3 + (k - 1)])) +
                                                                 ((v[(i + 1) * N2 * N3 + j * N3 + k] -
                                                                   v[(i - 1) * N2 * N3 + j * N3 + k]) *
                                                                  (t[i * N2 * N3 + j * N3 + (k + 1)] -
                                                                   t[i * N2 * N3 + j * N3 + (k - 1)]) -
                                                                  (v[i * N2 * N3 + j * N3 + (k + 1)] -
                                                                   v[i * N2 * N3 + j * N3 + (k - 1)]) *
                                                                  (t[(i + 1) * N2 * N3 + j * N3 + k] -
                                                                   t[(i - 1) * N2 * N3 + j * N3 + k]))));
                err_T = std::max(err_T, std::abs(t_new[index] - t[index]));
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
    int end = std::min(start + chunk_size, N1);

    // Initialize arrays
    std::vector<double> u(N1 * N2 * N3, 0.0);
    std::vector<double> u_new = u;
    std::vector<double> v(N1 * N2 * N3, 0.0);
    std::vector<double> v_new = v;
    std::vector<double> t(N1 * N2 * N3, 0.0);
    std::vector<double> t_new = t;

    int nr_it = 0;
    bool stop = false;
    std::cout << "Process " << rank << " starting" << std::endl;
    while (!stop && nr_it < maxIterations) {
        nr_it++;
        double err_u = 0.0, err_v = 0.0, err_T = 0.0;

        // Create threads for u, v, and t computation
        std::thread threadU(computeVelocityVectorU, start, end, N2, N3, h, Ra, std::ref(u), std::ref(t), std::ref(u_new),
                            std::ref(err_u));
        std::thread threadV(computeV, start, end, N2, N3, h, Ra, std::ref(v), std::ref(t), std::ref(v_new),
                            std::ref(err_v));
        std::thread threadT(computeTemperatureVector, start, end, N2, N3, h, std::ref(u), std::ref(v), std::ref(t), std::ref(t_new),
                            std::ref(err_T));

        // Wait for threads to complete
        threadU.join();
        threadV.join();
        threadT.join();

        if (rank == 0) {
            // Root process receives all chunks and updates the arrays
            for (int p = 1; p < size; ++p) {
                int p_start = p * chunk_size;
                int p_end = std::min(p_start + chunk_size, N1);

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
            // Other processes send their chunks to the root process
            MPI_Send(&u_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&v_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&t_new[start * N2 * N3], (end - start) * N2 * N3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        // Root process sends the updated arrays back to all processes
        MPI_Bcast(&u_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&v_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&t_new[0], N1 * N2 * N3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (nr_it % step == 0 && rank == 0) {
            writeData(u, v, t, N1, N2, N3, nr_it);
        }

        if (err_u < 1e-9 && err_v < 1e-9 && err_T < 1e-9) {
            stop = true;
            if (rank == 0)
                writeData(u, v, t, N1, N2, N3, nr_it);
        }

        u = u_new;
        v = v_new;
        t = t_new;

        MPI_Allreduce(MPI_IN_PLACE, &stop, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        if (nr_it >= maxIterations) {
            std::cout << "Reached maximum iterations without converging" << std::endl;
        } else {
            std::cout << "Converged in " << nr_it << " iterations" << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}
