#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

// Define constants
const double Ra = 100.0;          // Rayleigh number
const int N = 50;                 // Axis length
const double Ax = 1.0;            // Length of the cavity in x-direction
const double Ay = 1.0;            // Length of the cavity in y-direction
const double Az = 1.0;            // Length of the cavity in z-direction
const int maxIterations = 10000; // Maximum number of iterations
const int step = 100;             // Number of steps to output data

void writeData(const std::vector<std::vector<std::vector<double>>> &u,
               const std::vector<std::vector<std::vector<double>>> &v,
               const std::vector<std::vector<std::vector<double>>> &t,
               int N1, 
               int N2, 
               int N3, 
               int step)
{
    std::ofstream file("sequential_data_step_" + std::to_string(step) + ".txt");
    file << N1 << " " << N2 << " " << N3 << "\n";
    for (int i = 0; i < N1; ++i)
    {
        for (int j = 0; j < N2; ++j)
        {
            for (int k = 0; k < N3; ++k)
            {
                file << i << " " << j << " " << k << " " << u[i][j][k] << " " << v[i][j][k] << " " << t[i][j][k] << "\n";
            }
        }
    }
    file.close();
    std::cout << "Wrote data at step " << step << std::endl;
}

int main()
{
    int N1 = N * Ax;          // Number of nodes in x-direction
    int N2 = N * Ay;          // Number of nodes in y-direction
    int N3 = N;               // Number of nodes in z-direction
    double h = Ax / (N1 - 1); // Cubic grid cell size

    // Initialize arrays
    std::vector<std::vector<std::vector<double>>> u(N1, std::vector<std::vector<double>>(N2, std::vector<double>(N3, 0.0)));
    std::vector<std::vector<std::vector<double>>> u_new = u;
    std::vector<std::vector<std::vector<double>>> v(N1, std::vector<std::vector<double>>(N2, std::vector<double>(N3, 0.0)));
    std::vector<std::vector<std::vector<double>>> v_new = v;
    std::vector<std::vector<std::vector<double>>> t(N1, std::vector<std::vector<double>>(N2, std::vector<double>(N3, 0.0)));
    std::vector<std::vector<std::vector<double>>> t_new = t;

    int nr_it = 0;
    bool stop = false;
    std::cout << "Starting" << std::endl;
    while (!stop && nr_it < maxIterations)
    {
        nr_it++;
        double err_u = 0.0, err_v = 0.0, err_T = 0.0;

        for (int i = 1; i < N1 - 1; ++i)
        {
            for (int j = 1; j < N2 - 1; ++j)
            {
                for (int k = 1; k < N3 - 1; ++k)
                {
                    u_new[i][j][k] = Ra * h / 12.0 * (t[i][j + 1][k] - t[i][j - 1][k]) +
                                     1.0 / 6.0 * (u[i - 1][j][k] + u[i + 1][j][k] + u[i][j - 1][k] + u[i][j + 1][k] + u[i][j][k - 1] + u[i][j][k + 1]);

                    v_new[i][j][k] = Ra * h / 12.0 * (t[i + 1][j][k] - t[i - 1][j][k]) +
                                     1.0 / 6.0 * (v[i - 1][j][k] + v[i + 1][j][k] + v[i][j - 1][k] + v[i][j + 1][k] + v[i][j][k - 1] + v[i][j][k + 1]);

                    t_new[i][j][k] = 1.0 / 6.0 * (t[i - 1][j][k] + t[i + 1][j][k] + t[i][j - 1][k] + t[i][j + 1][k] + t[i][j][k - 1] + t[i][j][k + 1] + h * h - 1.0 / 4.0 * (((u[i][j][k + 1] - u[i][j][k - 1]) * (t[i][j + 1][k] - t[i][j - 1][k]) - (u[i][j + 1][k] - u[i][j - 1][k]) * (t[i][j][k + 1] - t[i][j][k - 1])) + ((v[i + 1][j][k] - v[i - 1][j][k]) * (t[i][j][k + 1] - t[i][j][k - 1]) - (v[i][j][k + 1] - v[i][j][k - 1]) * (t[i + 1][j][k] - t[i - 1][j][k]))));

                    // Error evaluation
                    err_u = std::max(err_u, std::abs(u_new[i][j][k] - u[i][j][k]));
                    err_v = std::max(err_v, std::abs(v_new[i][j][k] - v[i][j][k]));
                    err_T = std::max(err_T, std::abs(t_new[i][j][k] - t[i][j][k]));
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

        if (nr_it % step == 0)
        {
            writeData(u, v, t, N1, N2, N3, nr_it);
        }

        if (err_u < 1e-9 && err_v < 1e-9 && err_T < 1e-9)
        {
            stop = true;
            writeData(u, v, t, N1, N2, N3, nr_it);
        }
        else
        {
            u = u_new;
            v = v_new;
            t = t_new;
        }
    }

    if (nr_it >= maxIterations)
    {
        std::cout << "Reached maximum iterations without converging" << std::endl;
    }
    else
    {
        std::cout << "Converged in " << nr_it << " iterations" << std::endl;
    }

    return 0;
}
