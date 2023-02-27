#include <iostream>
#include <random>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
    vector<int> steps{-1,1};

    default_random_engine e(random_device{}());
    uniform_int_distribution<unsigned long> d{0, steps.size()-1};
    auto step = [&](){return steps[d(e)];};

    double N = atoi(argv[1]), M = atoi(argv[2]);

    vector<double> offsets(N);
    vector<double> offsets_abs(N);
    vector<double> offsets_2(N);

    for (int i = 0; i < M; i++)
    {
        vector<int> sample({0});
        for (int j = 0; j < N; j++)
        {
            auto next = sample[j] + step();
            sample.push_back(next);
            offsets[j] += next/M;
            offsets_abs[j] += abs(next)/M;
            offsets_2[j] += next*next/M;
        }
    }

    for (int i = 0; i < N; i++)
        cout << offsets[i] << ",";
    cout << "|";

    for (int i = 0; i < N; i++)
        cout << offsets_abs[i] << ",";
    cout << "|";

    for (int i = 0; i < N; i++)
        cout << offsets_2[i] << ",";
    cout << "|";
}