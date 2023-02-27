#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <random>

using namespace std;

int main(int argc, char **argv)
{
    int N = atoi(argv[1]);
    // int N = 10;
    vector<char> dict{'1', '2'};//, '3'};
    size_t len = 2;

    random_device r;
    auto monkey = std::bind(uniform_int_distribution<int>(0, dict.size() - 1), default_random_engine(r()));
    vector<pair<string, int>> master_res;

    // for (int m = 0; m < len; m++)
        for (int l = 0; l < len; l++)
            for (int k = 0; k < len; k++)
            {
                string phrase({/*dict.at(m),*/dict.at(l),dict.at(k)});
                cout << phrase << " ";
                vector<int> res;
                for (int j = 0; j < N; j++)
                {
                    string paper("");
                    for (int i = 0; i < len; i++)
                        paper.append({dict.at(monkey())});
                    while (paper.find(phrase) == string::npos)                        
                        paper.append({dict.at(monkey())});
                    res.push_back(paper.length());
                }
                auto av = accumulate(res.begin(), res.end(), 0.0) / (double)N;
                cout << av << endl;
                master_res.push_back({phrase, av});
            }
    // for (auto res : master_res)
    //     cout << res.first << " " << res.second << endl;
}