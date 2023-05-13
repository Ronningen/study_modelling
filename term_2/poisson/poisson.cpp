#include <fstream>
#include <vector>
#include <ctime>
#include <iostream>

using namespace std;

#define threads 1

struct Current
{
    double l, r, t, b;
    double integral()
    {
        // return l + r + b + t;
        return l - r + b - t;
    }
};

struct Cell
{
    double sigma, potential;
    Current current;
};

struct Mesh
{
    double dx;
    size_t _n, _m;
    vector<vector<Cell>> mesh;
    vector<double> potential_bc;

    /// @brief create rectangulare mesh m x n
    /// @param n - number of rows
    /// @param m - number of columns
    Mesh(double dx, size_t n, size_t m, initializer_list<double> potential_bc,
         double (*sigma)(size_t, size_t),
         double (*initial_potential)(size_t, size_t)) : dx(dx), _n(n + 2), _m(m + 2), potential_bc(potential_bc)
    {
        mesh = vector<vector<Cell>>(_n);
        for (size_t row = 0; row < _n; row += _n - 1)
        {
            mesh[row] = vector<Cell>(_m);
            // for (size_t col = 0; col < _m; col++)
            //     bottom-top bc
        }
        for (size_t row = 1; row < _n - 1; row++)
        {
            mesh[row] = vector<Cell>(_n);
            for (size_t col = 0; col < 2; col++)
            {
                // mesh[row][col * (_m - 1)].potential = *(potential_bc.begin() + col);
                mesh[row][col * (_m - 1)].sigma = sigma(row - 1, col * (_m - 3));
            }
            for (size_t col = 1; col < _m - 1; col++)
            {
                mesh[row][col].sigma = sigma(row - 1, col - 1);
                mesh[row][col].potential = initial_potential(row - 1, col - 1);
            }
        }
    }

    void traversal(double dt)
    {
        // tuning boundary conditions
        // #pragma omp parallel for num_threads(threads) collapse(2)
        for (size_t row = 1; row < _n - 1; row++)
        {
            for (size_t col = 0; col < 2; col++)
            {
                double potential = *(potential_bc.begin() + col);
                mesh[row][col * (_m - 1)].potential = 2 * potential - mesh[row][col * (_m - 3) + 1].potential;
            }
        }
        // calculationg currents
        #pragma omp parallel for num_threads(threads) collapse(2)
        for (size_t row = 1; row < _n - 1; row++)
        {
            for (size_t col = 1; col < _m - 1; col++)
            {
                mesh[row][col].current.l = (mesh[row][col - 1].potential - mesh[row][col].potential) / dx * 2 / (1 / mesh[row][col - 1].sigma + 1 / mesh[row][col].sigma);
                mesh[row][col].current.r = (mesh[row][col].potential - mesh[row][col + 1].potential) / dx * 2 / (1 / mesh[row][col].sigma + 1 / mesh[row][col + 1].sigma);
                mesh[row][col].current.b = (mesh[row - 1][col].potential - mesh[row][col].potential) / dx * 2 / (1 / mesh[row - 1][col].sigma + 1 / mesh[row][col].sigma);
                mesh[row][col].current.t = (mesh[row][col].potential - mesh[row + 1][col].potential) / dx * 2 / (1 / mesh[row][col].sigma + 1 / mesh[row + 1][col].sigma);
            }
        }
        // evolution
        #pragma omp parallel for num_threads(threads) collapse(2)
        for (size_t row = 1; row < _n - 1; row++)
        {
            for (size_t col = 1; col < _m - 1; col++)
            {
                mesh[row][col].potential += mesh[row][col].current.integral() / dx * dt;
            }
        }
    }
};

void print(Mesh mesh, ofstream &f, bool border = false)
{
    for (size_t row = 1 - border; row < mesh._n - 1 + border; row++)
    {
        for (size_t col = 1 - border; col < mesh._m - 1 + border; col++)
        {
            f << mesh.mesh[row][col].potential << ',';
        }
        f << endl;
    }
}

int main(int argc, char **argv)
{   
    auto start = clock();

    long long N = 1000;
    if (argc > 1)
        N = atoll(argv[1]);

    size_t n = 10, m = 10, scale = 300;
    double dx = 1, dt = 0.1;

    Mesh mesh(
        dx / scale, n * scale, m * scale, {200, 000},
        [](size_t r, size_t c)
        { return double(c > 500 ? 1 : 5); },
        [](size_t r, size_t c)
        { return double(0); });

    for (long long i = 0; i < N; i++)
    {
        mesh.traversal(dt);
    }

    auto end = clock();
    cout << (end - start) << endl;

    ofstream file;
    file.open("data.csv");
    print(mesh, file, false);
}