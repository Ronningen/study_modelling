#include <iostream>
#include <random>
#include <algorithm>

template <typename DIST, typename URNG>
std::vector<double> run(DIST dist, URNG engine, int N)
{
    std::vector<double> res;
    double buf;
    for (int i = 0; i < N; i++)
    {
        buf = dist(engine);
        std::cout << buf << " ";
        res.push_back(buf);
    }
    return res;
}

std::vector<int> do_hist(std::vector<double> sample, std::vector<double> bins)
{
    std::vector<int> weights(bins.size() - 1);
    for (auto s : sample)
        if (s >= bins[0])
        {
            for (int i = 1; i < bins.size(); i++)
                if (s < bins[i])
                {
                    weights[i - 1]++;
                    break;
                }
            if (s == bins.back())
                weights[weights.size() - 1]++;
        }
    return weights;
}

template <typename T = double>
class any_dist
{
    T (*F)(T);
    std::uniform_real_distribution<T> base{0, 1};
    T y0;

    T const STEP = 10;
    T const ERR = 0.01;

public:
    // where F is Cumalative Distribution Function
    any_dist(T (*F)(T)) : F(F), y0(F(0)) {}

    // mutate x ~ uniform[0,1] to y ~ F solving x = F(y) with binary search
    template <class URNG>
    T operator()(URNG &engine)
    {
        T y = base(engine);
        T x = 0;
        if (y == y0)
            return y0;
        T d = STEP * ((y > y0)-0.5)*2;
        while (y * d > F(x) * d)
            x += d;

        T l = std::min(x, x - d);
        T r = std::max(x, x - d);

        while (r - l > ERR)
        {
            T c = (r + l)/2;
            if (F(c) > y)
                r = c;
            else
                l = c;
        }

        return l;
    }
};

template <typename T = double>
T triangle_F(T x)
{
    T a = -5, b = 5;
    if (x < a) return 0;
    if (x > b) return 1;
    T c = (a + b)/2, na = 0.5/(a-c)/(a-c), nb = 0.5/(b-c)/(b-c);
    if (x < c) return (x-a)*(x-a)*na;
    if (x > c) return 1-(x-b)*(x-b)*nb;
    return 0.5;
}

int main(int argc, char **argv)
{
    float N, K, S;
    int mod;
    if (argc <= 1)
    {
        mod = 4, N = 2;
        // std::cout << "print N: ";
        // std::cin >> N;
        // std::cout << "print K: ";
        // std::cin >> K;
        // std::cout << "print Sigma: ";
        // std::cin >> S;
    }
    else
    {
        mod = atoi(argv[1]);
        N = atoi(argv[2]);
        K = atof(argv[3]);
    }

    std::default_random_engine e(0);
    std::vector<double> sample;
    switch (mod)
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
    case 4:
        sample = run(any_dist<double>(triangle_F<double>), e, N);
        // sample = run(any_dist<double>([](double x){return x;}), e, N);
        break;
    default:
        std::cout << "fail";
        return -1;
    }

    std::vector<double> bins({-6,-5,-4,-3,-2,-1,-0.6,-0.3,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6});
    std::vector<int> weights = do_hist(sample, bins);

    std::cout << std::endl
              << "sample: ";
    for (auto samp : sample)
        std::cout << samp << " ";

    std::cout << std::endl
              << "bins: ";
    for (auto bin : bins)
        std::cout << bin << " ";

    std::cout << std::endl
              << "weights: ";
    for (auto weight : weights)
        std::cout << weight << " ";
}