#include <iostream>
#include <vector>
#include <map>

using namespace std;

double next(double x, double a)
{
    return a*x*(1-x);
}

int main(int argc, char** argv)
{
    double x = atof(argv[1]);
    map<double, vector<double>> result;
    for (double a = 3; a <= 4; a += 0.001)
    {
        result[a] = vector<double>({x});
        for (int i = 0; i < 500; i++)
        {
            double t = next(result[a].back(), a);
            result[a].push_back(t);
        }
    }
    for (auto r : result)
    {
        cout << r.first << ": ";
        for (auto e : r.second)
            cout << e << " ";
        cout << "\n";
    }
}