#include <iostream>
#include <string>
#include <map>
#include <random>
#include <cmath>

template <typename DIST, typename URNG>
std::vector<double> run(DIST dist, URNG engine, int N)
{
    std::vector<double> res;
    auto bind = std::bind(dist, engine);
    double buf;
    for (int i = 0; i < N; i++)
    {
        buf = bind();
        std::cout << buf << " ";
        res.push_back(buf);
    }
    return res;
}

int main(int argc, char **argv)
{
    float N, K, S;
    if (argc <= 2)
    {
        std::cout << "print N: ";
        std::cin >> N;
        std::cout << "print K: ";
        std::cin >> K;
        std::cout << "print Sigma: ";
        std::cin >> S;
    }

    N = atoi(argv[2]);
    K = atof(argv[3]);

    std::default_random_engine e(0);
    std::vector<double> sample;
    switch (atoi(argv[1]))
    {
    case 1:
        sample = run(std::uniform_int_distribution<int>(1, K), e, N);
        break;
    case 2:
        sample = run(std::uniform_real_distribution<double>(1, K), e, N);
        break;
    case 3:
        sample = run(std::normal_distribution<double>(K, atof(argv[4])), e, N);
        break;
    default:
        std::cout << "fail";
        return -1;
    }

    std::vector<double> bins({-10, -5, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 5, 10});
    std::vector<int> weights(bins.size() - 1);
    for (auto s : sample)
        if (s >= bins[0])
        {
            for (int i = 1; i < bins.size(); i++)
                if (s < bins[i])
                {
                    weights[i-1]++;
                    break;
                }
            if (s == bins.back())
                weights[weights.size()-1]++;
        }
    std::cout << std::endl << "bins: ";
    for (auto bin : bins)
        std::cout << bin << " ";
    
    std::cout << std::endl << "weights: ";
    for (auto weight : weights)
        std::cout << weight << " ";
}