#include <iostream>
#include <random>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
    default_random_engine e(random_device{}());
    // vector<int> steps{-1,1};
    // uniform_int_distribution<unsigned long> d{0, steps.size()-1};
    // auto step = [&](){return steps[d(e)];};
    normal_distribution<> d{0, 2};
    auto step = [&](){return d(e);};

    double N = atoi(argv[1]), M = atoi(argv[2]);

    vector<double> offsets(N);
    vector<double> offsets_l(N);
    vector<double> offsets_abs(N);
    vector<double> offsets_2(N);

    for (int i = 0; i < M; i++)
    {
        vector<double> sample({0});
        for (int j = 0; j < N; j++)
        {
            auto next = sample[j] + step();
            if (abs(sample[j])>10) continue;
            sample.push_back(next);
            offsets[j] += next;
            offsets_abs[j] += abs(next);
            offsets_2[j] += next*next;
            offsets_l[j]++;
        }
        double s = 0;
        // for (int j = 0; j < N; j++)
        // {
        //     s += step();
        // }
        // cout << s << " ";
    }

    for (int i = 0; i < N; i++)
        cout << offsets[i]/offsets_l[i] << ",";
    cout << "|";

    for (int i = 0; i < N; i++)
        cout << offsets_abs[i]/offsets_l[i] << ",";
    cout << "|";

    for (int i = 0; i < N; i++)
        cout << offsets_2[i]/offsets_l[i] << ",";
    cout << "|";
}