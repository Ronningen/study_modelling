#include <iostream>
#include <fstream>
#include "json.hpp"
#include <valarray>
#include <vector>
#include <ctime>
#include <memory>

#pragma region Vectors

template <typename type>
using NaiveVector = std::valarray<type>; // Vector for naive ariphmetics
template <typename type>
type norm2(const NaiveVector<type> &v)
{
    return (v * v).sum();
}

template <typename type>
struct KahanVector final // Vector with Kahan summation algorithm
{
    // Store data in serias, which is summed only on +=, -=, [] operators being invoked

    KahanVector(std::initializer_list<type> ls) // ls.size() > 0
        : value({std::valarray<type>(ls)}), error(std::valarray<type>(ls.size()))
    {
    }

    KahanVector(size_t n) // n > 0
        : value({std::valarray<type>(n)}), error(std::valarray<type>(n))
    {
    }

    KahanVector<type> operator+(const KahanVector<type> &r) const &
    {
        KahanVector<type> tmp(*this);
        tmp.value.insert(tmp.value.end(), r.value.begin(), r.value.end());
        return tmp;
    }

    KahanVector<type> &operator+(const KahanVector<type> &r) &&
    {
        this->value.insert(this->value.end(), r.value.begin(), r.value.end());
        return *this;
    }

    KahanVector<type> operator-(const KahanVector<type> &r) const &
    {
        return *this + -r;
    }

    KahanVector<type> &operator-(const KahanVector<type> &r) &&
    {
        return std::move(*this) + -r;
    }

    KahanVector<type> &operator+=(const KahanVector<type> &r)
    {
        value.insert(value.end(), r.value.begin(), r.value.end());
        collapse();
        return *this;
    }

    KahanVector<type> &operator-=(const KahanVector<type> &r)
    {
        return *this += -r;
    }

    KahanVector<type> operator*(type r) const &
    {
        return KahanVector<type>(*this) *= r;
    }

    KahanVector<type> &operator*(type r) &&
    {
        return *this *= r;
    }

    KahanVector<type> operator/(type r) const &
    {
        return KahanVector<type>(*this) /= r;
    }

    KahanVector<type> &operator/(type r) &&
    {
        return *this /= r;
    }

    KahanVector<type> &operator*=(type r)
    {
        for (auto v = value.begin(); v != value.end(); v++)
            *v *= r;
        error *= r;
        return *this;
    }

    KahanVector<type> &operator/=(type r)
    {
        return *this *= (1 / r);
    }

    KahanVector<type> operator-() const &
    {
        return *this * -1;
    }

    KahanVector<type> &operator-() &&
    {
        return *this * -1;
    }

    template <typename index>
    auto operator[](index i)
    {
        collapse();
        return value.front()[i];
    }

    template <typename index>
    const auto operator[](index i) const
    {
        collapse();
        return value.front()[i];
    }

    size_t size() const
    {
        return error.size();
    }

    template <typename type_>
    friend type_ norm2(const KahanVector<type_> &v);

private:
    void collapse()
    {
        if (value.size() == 0)
        {
            std::cerr << "emty vector";
            throw std::logic_error("emty vector");
        }
        for (auto v = value.begin() + 1; v != value.end(); v++)
        {
            auto y = *v - error;
            auto t = value.front() + y;
            error = (t - value.front()) - y;
            value.front() = t;
        }
        value.resize(1);
    }
    void collapse() const
    {
        const_cast<KahanVector<type> *>(this)->collapse();
    }

    std::vector<std::valarray<type>> value;
    std::valarray<type> error;
};

template <typename type_, typename type>
KahanVector<type> operator*(type_ r, const KahanVector<type> &v)
{
    return v * r;
}

template <typename type_, typename type>
KahanVector<type> &operator*(type_ r, KahanVector<type> &&v)
{
    return std::move(v) * r;
}

template <typename type>
type norm2(const KahanVector<type> &v)
{
    v.collapse();
    return norm2<type>(v.value.front());
}

#pragma endregion // TODO: NeumaierVector, KleinVector, etc.

#pragma region Problems

template <template <typename type> class vector, typename type>
struct Problem // Generelized Cauchy problem as system of first-order ODE y' = f(x, y)
{
    Problem(const vector<type> &y0) : y0(y0) {}
    virtual vector<type> operator()(const type x, const vector<type> &y) const & = 0;

    const vector<type> y0;
};
template <template <typename type> class vector, typename type>
struct IAnalyticalProblem // Problem with khown analytical solution
{
    virtual vector<type> AnalyticalValue(const type x) const & = 0;
};
template <template <typename type> class vector, typename type>
struct IHaveInvariantProblem // Problem with an invariant with respect to x (integral of motion)
{
    virtual type Invariant(const vector<type> &y) const & = 0;
};

template <template <typename type> class vector, typename type>
struct SimplestOscillator final : Problem<vector, type>, IAnalyticalProblem<vector, type>, IHaveInvariantProblem<vector, type> // 1-dimensional harmonic oscillator x = y[0], v = y[1]
{
    SimplestOscillator(const vector<type> &y0, type w) : Problem<vector, type>(y0),
                                                         w(w), w2(w * w),
                                                         A(sqrt(y0[1] * y0[1] / w2 + y0[0] * y0[0])),
                                                         initial_phase(atan(y0[1] / y0[0] / w)) {}

    vector<type> operator()(const type x, const vector<type> &y) const & override
    {
        return {{y[1], -w2 * y[0]}};
    }

    vector<type> AnalyticalValue(const type x) const & override
    {
        return {{A * cos(w * x + initial_phase),
                 -A * w * sin(w * x + initial_phase)}};
    }

    type Invariant(const vector<type> &y) const & override
    {
        type x_ = y[0], v_ = y[1];
        return (v_ * v_ + w2 * x_ * x_) / 2;
    }

    const type w, w2, A, initial_phase;
};

template <template <typename type> class vector, typename type>
struct IdialPhysicalPendulum final : Problem<vector, type> // 1-dimensional harmonic oscillator phi = y[0], omega = y[1]
{
    IdialPhysicalPendulum(const vector<type> &y0, type w) : Problem<vector, type>(y0), w(w) {}

    vector<type> operator()(const type x, const vector<type> &y) const & override
    {
        return {{y[1], -w * sin(y[0])}};
    }

    const type w;
};

template <template <typename type> class vector, typename type>
struct RealPhysicalPendulum final : Problem<vector, type> // 1-dimensional harmonic oscillator phi = y[0], omega = y[1]
{
    RealPhysicalPendulum(const vector<type> &y0, type w = 1, type gamma = 0.1) : Problem<vector, type>(y0), w(w), gamma(gamma) {}

    vector<type> operator()(const type x, const vector<type> &y) const & override
    {
        return {{y[1], -2 * gamma * y[1] - w * sin(y[0])}};
    }

    const type w, gamma;
};

template <template <typename type> class vector, typename type>
struct LimitHollowEarth final : Problem<vector, type>, IHaveInvariantProblem<vector, type> // 1-dimensional hollow Earth with thin surface problem x = y[0], v = y[1]
{
    LimitHollowEarth(const vector<type> &y0, type GM, type R) : Problem<vector, type>(y0), GM(GM), R(R), U0(-GM / R) {}

    vector<type> operator()(const type x, const vector<type> &y) const & override
    {
        type x_ = y[0], v_ = y[1];

        if (x_ > R)
            return {{v_, -GM / x_ / x_}};
        else if (x_ < -R)
            return {{v_, GM / x_ / x_}};
        else
            return {{v_, 0}};
    }

    type Invariant(const vector<type> &y) const & override
    {
        type x_ = y[0], v_ = y[1];

        if (abs(x_) < R)
            return v_ * v_ / 2 + U0;
        else
            return v_ * v_ / 2 - GM / abs(x_);
    }

    const type GM, R, U0;
};

template <template <typename type> class vector, typename type>
struct HollowEarth final : Problem<vector, type>, IHaveInvariantProblem<vector, type> // 1-dimensional hollow Earth with bold surface problem x = y[0], v = y[1]
{
    HollowEarth(const vector<type> &y0, type GM, type R, type r) : Problem<vector, type>(y0), GM(GM), R(R), r(r), dR3(R * R * R - r * r * r), r3(r * r * r) {}

    vector<type> operator()(const type x, const vector<type> &y) const & override
    {
        type x_ = y[0], v_ = y[1];

        if (x_ > R)
            return {{v_, -GM / x_ / x_}};
        else if (x_ > r)
            return {{v_, -GM * (x_ * x_ * x_ - r3) / x_ / x_ / dR3}};

        else if (x_ < -R)
            return {{v_, GM / x_ / x_}};
        else if (x_ < -r)
            return {{v_, GM * (-x_ * x_ * x_ - r3) / x_ / x_ / dR3}};

        else
            return {{y[1], 0}};
    }

    type Invariant(const vector<type> &y) const & override
    {
        type x_ = y[0], v_ = y[1];

        if (abs(x_) > R)
            return v_ * v_ / 2 - GM / abs(x_);
        else if (abs(x_) > r)
            return v_ * v_ / 2 - GM * (abs(x_ * x_ * x_) - r3) / abs(x_) / dR3; // do not know how to calculate correctly
        else
            return 0;
    }

    const type GM, R, r, dR3, r3;
};

// parse:

template <template <typename type> class vector, typename type>
const std::shared_ptr<Problem<vector, type>> parse_problem(const nlohmann::json &run)
{
    auto problem_j = run["problem"];
    std::string problem_s = problem_j["type"];
    if (problem_s == "simplest_oscillator")
        return std::shared_ptr<Problem<vector, type>>{new SimplestOscillator<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["w"])};
    else if (problem_s == "ideal_physical_pendulum")
        return std::shared_ptr<Problem<vector, type>>{new IdialPhysicalPendulum<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["w"])};
    else if (problem_s == "physical_pendulum")
        return std::shared_ptr<Problem<vector, type>>{new RealPhysicalPendulum<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["w"] /*, problem_j["gamma"]*/)};
    else if (problem_s == "hollow_earth")
        return std::shared_ptr<Problem<vector, type>>{new HollowEarth<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["GM"], problem_j["R"], problem_j["r"])};
    else if (problem_s == "limit_hollow_earth")
        return std::shared_ptr<Problem<vector, type>>{new LimitHollowEarth<vector, type>({problem_j["x0"], problem_j["v0"]}, problem_j["GM"], problem_j["R"])};

    throw std::runtime_error("invalid configuration json");
}

#pragma endregion

#pragma region Constraints

template <template <typename type> class vector, typename type>
struct IConstraint // solver's interations constraint
{
    virtual bool operator()(const type x, const vector<type> &y, unsigned long long i) const = 0;
};

template <template <typename type> class vector, typename type>
struct СounterConstraint final : public IConstraint<vector, type> // constraint on the amount of iterations
{
    СounterConstraint(unsigned long long N) : N(N) {}
    bool operator()(const type x, const vector<type> &y, unsigned long long i) const override
    {
        return (i < N);
    }

private:
    const unsigned long long N;
};

template <template <typename type> class vector, typename type>
struct AnalyticalDeviationConstraint final : public IConstraint<vector, type> // constraint on reletive deviation of chosen coordinates
{
    AnalyticalDeviationConstraint(
        const IAnalyticalProblem<vector, type> &problem, const vector<type> &y0, const std::__1::slice &comparison_mask, const type reletive_deviation_limit)
        : problem(problem), comparison_mask(comparison_mask),
          deviation_limit2(norm2<type>(y0[comparison_mask]) * reletive_deviation_limit * reletive_deviation_limit) {}

    type current_deviation2(const type x, const vector<type> &y) const
    {
        return norm2<type>((problem.AnalyticalValue(x) - y)[comparison_mask]);
    }

    bool operator()(const type x, const vector<type> &y, unsigned long long i) const override
    {
        return current_deviation2(x, y) < deviation_limit2;
    }

private:
    const IAnalyticalProblem<vector, type> &problem;
    const std::slice comparison_mask;
    const type deviation_limit2;
};

template <template <typename type> class vector, typename type>
struct InvariantDeviationConstraint final : public IConstraint<vector, type> // constraint on reletive deviation of the integral of motion
{
    InvariantDeviationConstraint(
        const IHaveInvariantProblem<vector, type> &problem, const vector<type> &y0, const type reletive_deviation_limit)
        : problem(problem), invariant(problem.Invariant(y0)), deviation_limit(invariant * reletive_deviation_limit) {}

    bool operator()(const type x, const vector<type> &y, unsigned long long i) const override
    {
        return abs(problem.Invariant(y) - invariant) < deviation_limit;
    }

protected:
    const IHaveInvariantProblem<vector, type> &problem;
    const type invariant;
    const type deviation_limit;
};

// parse:

template <template <typename type> class vector, typename type>
const std::shared_ptr<IConstraint<vector, type>> parse_constraint(const nlohmann::json &run, const Problem<vector, type> &problem)
{
    auto cons_j = run["constraint"];
    std::string cons_s = cons_j["type"];
    if (cons_s == "counter")
        return std::shared_ptr<IConstraint<vector, type>>{new СounterConstraint<vector, type>((unsigned long long int)cons_j["N"])};
    else if (cons_s == "analytical")
    {
        auto mask = cons_j["comparison_mask"];
        return std::shared_ptr<IConstraint<vector, type>>{new AnalyticalDeviationConstraint<vector, type>(
            *dynamic_cast<IAnalyticalProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem)),
            problem.y0, std::slice(mask["start"], mask["size"], mask["stride"]), (type)cons_j["reletive_deviation_limit"])};
    }
    else if (cons_s == "invariant")
        return std::shared_ptr<IConstraint<vector, type>>{new InvariantDeviationConstraint<vector, type>(
            *dynamic_cast<IHaveInvariantProblem<vector, type> *>(const_cast<Problem<vector, type> *>(&problem)),
            problem.y0, (type)cons_j["reletive_deviation_limit"])};

    throw std::runtime_error("invalid configuration json");
}

#pragma endregion // TODO: bound constraint, multiconstraint

#pragma region Printer

// general printer params
struct PrinterConfig
{
    std::ostream *stream = nullptr;
    std::string el_sep, zone_sep, row_sep, run_sep;
    bool do_log;
};


// Class for print uotput data into a stream
template <template <typename type> class vector, typename type>
struct Printer final
{
    Printer(const PrinterConfig &conf, const Problem<vector, type> &problem) : conf(conf),
        // A(dynamic_cast<const IAnalyticalProblem<vector, type>*>(&problem)),///////////////////
        // I(dynamic_cast<const IHaveInvariantProblem<vector, type>*>(&problem)) {}
        A(dynamic_cast<IAnalyticalProblem<vector, type>*>(const_cast<Problem<vector, type>*>(&problem))),
        I(dynamic_cast<IHaveInvariantProblem<vector, type>*>(const_cast<Problem<vector, type>*>(&problem))){}

    // printing current
    void print(const type x, const vector<type> &y) const
    {
        if (!conf.do_log)
            return;
        print(x);
        print(y);
        if (I != nullptr)
        {
            print();
            print(I->Invariant(y));
        }
        if (A != nullptr)
        {
            print();
            print(A->AnalyticalValue(x));
        }
        *conf.stream << conf.row_sep;
    }

    // print run's ending and results
    void stop(const clock_t time, const unsigned long long int n, const type x, const vector<type> &y, const vector<type> &y0) const
    {
        *conf.stream << conf.run_sep << " time: " << time << ", iteraitions: " << n;
        if (I != nullptr)
        {
            *conf.stream << ", dI: ";
            *conf.stream << I->Invariant(y) - I->Invariant(y0);
        }
        if (A != nullptr)
        {
            *conf.stream << ", dy: ";
            print(A->AnalyticalValue(x) - y);
        }
        *conf.stream << conf.row_sep;
    }

private:
    // printing zone-separator
    void print() const
    {
        *conf.stream << conf.zone_sep << conf.el_sep;
    }
    // printing a value
    void print(const type x) const
    {
        *conf.stream << x << conf.el_sep;
    }
    // printing a vector
    void print(const vector<type> &y) const
    {
        for (int i = 0; i < y.size(); i++)
            *conf.stream << y[i] << conf.el_sep;
    }

    const IAnalyticalProblem<vector, type> *A;
    const IHaveInvariantProblem<vector, type> *I;
    const PrinterConfig &conf;
};

#pragma endregion

#pragma region Solver

struct ISolver
{
    virtual void run() = 0;
};

template <template <typename type> class vector, typename type>
struct DisposableSolver final : ISolver // Iterative solver of Cauchy problem
{
    // problem - Cauchy problem to solve
    // delta - x stride of one iteration
    // method - function returned y_{n+1} vector 
    DisposableSolver(const Problem<vector, type> &problem, const IConstraint<vector, type> &cons, 
                     const Printer<vector, type>& printer, type delta,
                     void (*method)(vector<type> &y, const type x, const type delta, const Problem<vector, type> &problem))
        : problem(problem), cons(cons), delta(delta), y(problem.y0), x(0), method(method), printer(printer) {}

    // do iterate
    void next()
    {
        method(y, x, delta, problem);
        x += delta;
    }

    // iterates till cons, prints logs of each iteration and results including elapsed time
    void run() override
    {
        clock_t start_time = clock();
        unsigned long long i = 0;
        for (; cons(x, y, i); i++)
        {
            next();
            printer.print(x, y);
        }
        printer.stop(clock() - start_time, i, x, y, problem.y0);
    }

protected:
    type x;
    vector<type> y;

    const Problem<vector, type> &problem;
    const IConstraint<vector, type> &cons;
    void (*method)(vector<type> &y, const type x, const type delta, const Problem<vector, type> &problem);
    const type delta;
    const Printer<vector, type> printer;
};

// methods:

template <template <typename type> class vector, typename type>
void euler(vector<type> &y, const type x, const type delta, const Problem<vector, type> &f)
{
    y += delta * f(x, y);
}

template <template <typename type> class vector, typename type>
void heun(vector<type> &y, const type x, const type delta, const Problem<vector, type> &f)
{
    vector<type> k = y + delta * f(x, y);
    y += delta * (f(x, y) + f(x + delta, k)) / 2;
}

template <template <typename type> class vector, typename type>
void runge_kutta(vector<type> &y, const type x, const type delta, const Problem<vector, type> &f)
{
    vector<type> k1 = f(x, y);
    vector<type> k2 = f(x + delta / 2, y + delta / 2 * k1);
    vector<type> k3 = f(x + delta / 2, y + delta / 2 * k2);
    vector<type> k4 = f(x + delta, y + delta * k3);
    y += delta / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
}

// parse:

template <template <typename type> class vector, typename type>
std::shared_ptr<DisposableSolver<vector, type>> parse_disposable_solver(const nlohmann::json &run, const Problem<vector, type> &problem, const IConstraint<vector, type> &cons, const PrinterConfig& conf)
{
    Printer<vector, type> printer(conf, problem);
    std::string method = run["method"];
    type delta = run["delta"];
    if (method == "euler")
        return std::shared_ptr<DisposableSolver<vector, type>>{new DisposableSolver<vector, type>(problem, cons, printer, delta, euler<vector, type>)};
    else if (method == "heun")
        return std::shared_ptr<DisposableSolver<vector, type>>{new DisposableSolver<vector, type>(problem, cons, printer, delta, heun<vector, type>)};
    else if (method == "runge_kutta")
        return std::shared_ptr<DisposableSolver<vector, type>>{new DisposableSolver<vector, type>(problem, cons, printer, delta, runge_kutta<vector, type>)};

    throw std::runtime_error("invalid configuration json");
}

#pragma endregion

#pragma region Manager

struct SolvingManager final
{
    SolvingManager() = default;

    void from_json(nlohmann::json config)
    {
        auto head = config["head"];
        printer_config.el_sep = head["el_sep"],
        printer_config.zone_sep = head["zone_sep"],
        printer_config.row_sep = head["row_sep"],
        printer_config.run_sep = head["run_sep"];
        printer_config.do_log = head["do_log"];

        std::string stream_s = head["stream"];
        if (stream_s == "std")
            printer_config.stream = &std::cout; // TODO: printer

        for (auto run : config["runs"])
            parse_type(run);
    }

    void run_all()
    {
        for (auto solver : solvers)
            solver->run();
        solvers.clear();
    }

private:
    void parse_type(const nlohmann::json &run)
    {
        std::string type = run["type"];
        if (type == "float")
            parse_vector<float>(run);
        else if (type == "double")
            parse_vector<double>(run);
        else if (type == "long double")
            parse_vector<long double>(run);
    }

    template <typename type>
    void parse_vector(const nlohmann::json &run)
    {
        std::string vector = run["vector"];
        if (vector == "naive")
            parse_solver<NaiveVector, type>(run);
        else if (vector == "kahan")
            parse_solver<KahanVector, type>(run);
    }

    template <template <typename type> class vector, typename type>
    void parse_solver(const nlohmann::json &run)
    {
        std::shared_ptr<Problem<vector, type>> problem = parse_problem<vector, type>(run);
        std::shared_ptr<IConstraint<vector, type>> cons = parse_constraint<vector, type>(run, *problem);
        // std::cout << typeid(cons).name() << std::endl;
        solvers.push_back(parse_disposable_solver<vector, type>(run, *problem, *cons, printer_config));
    }

    PrinterConfig printer_config{};
    std::vector<std::shared_ptr<ISolver>> solvers;
};

#pragma endregion

int main()
{
    SolvingManager manager;
    try
    {
        std::ifstream f("/Users/samedi/Documents/факультатив/study_modelling/term_1/solver/config.json");
        nlohmann::json config = nlohmann::json::parse(f);
        manager.from_json(config);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Configuration failed. " << e.what() << '\n';
        return -1;
    }
    //manager.run_all();
}
