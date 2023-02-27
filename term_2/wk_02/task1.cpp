#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <random>
#include <unordered_set>
#include <cmath>

using namespace std;

int main(int argc, char **argv)
{
    int K = atoi(argv[1]);
    // int N = atoi(argv[2]);
    
    for (int N = 1; pow(K, N) < 1000000; N++)
    {
        vector<int> row(N);
        int all = 0;
        int valid = 0;
        while (row[N - 1] < K && all < pow(K, N))
        {
            all++;

            unordered_set<int> s(row.begin(), row.end());
            valid += int(s.size() == K);

            for (int i = 0; i < N - 1; i++)
                if (++row[i] >= K)
                    row[i] = 0;
                else
                    break;
        }
        cout << double(valid) / double(all) << " ";
    }

    random_device r;
    auto choas = std::bind(uniform_int_distribution<int>(0, K-1), default_random_engine(r()));
    cout << "| ";

    int limit = 50;
    for (int N = 1; N < limit; N++)
    {
        vector<int> row(N);
        int all = 0;
        int valid = 0;
        for (int i = 0; i < 100000/(limit-N); i++)
        {
            all++;

            for (auto p = row.begin(); p != row.end(); p++)
                *p = choas();

            unordered_set<int> s(row.begin(), row.end());
            valid += int(s.size() == K);
        }
        cout << double(valid) / double(all) << " ";
    }
}