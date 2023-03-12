#include <iostream>
#include <fstream>
#include "json.hpp"
#include <vector>
#include <unordered_map>

using namespace std;

struct DotType
{
    string name;
    double mass;
    double radius;
};

struct Dot
{
    DotType type;
    double x;
    double vx;
};

struct Box
{
    double l, r;
};

class Gas
{
    Box box;
    vector<Dot> dots;
    size_t N;
    unordered_map<string, size_t> Ns;

    struct CollisionInfo
    {
        double time = INFINITY;
        size_t left_index;
    };

    // finds the nearest in time collidion with the least x (min index or the most left one)
    CollisionInfo collision_search_l()
    {
        CollisionInfo result;
        double tmptime;

        //  handle wall-dot
        if (dots[0].vx != 0)
        {
            tmptime = (dots[0].x - box.l - dots[0].type.radius) / (-dots[0].vx);
            if (signbit(tmptime) == 0)
                {
                    result.time = tmptime;
                    result.left_index = -1;
                }
        }
        
        // handel dot-dot
        for (int i = 0; i < N-1; i++)
        if (dots[i].vx - dots[i+1].vx != 0)
        {
            tmptime = (dots[i+1].x - dots[i].x - dots[i+1].type.radius - dots[i].type.radius) 
                    / (dots[i].vx - dots[i+1].vx);
            if (signbit(tmptime) == 0 && tmptime < result.time)
            {
                result.time = tmptime;
                result.left_index = i;
            }
        }
    
        // hadle dot-wall
        if (dots[0].vx != 0)
        {
            tmptime = (box.r - dots[N-1].x - dots[N-1].type.radius) / (dots[N-1].vx);
            if (signbit(tmptime) == 0 && tmptime < result.time)
                {
                    result.time = tmptime;
                    result.left_index = N-1;
                }
        }
        
        return result;
    }

    // exchanges velosities
    void collide(CollisionInfo collision)
    {
        size_t i = collision.left_index;
        if (i == -1)
            dots[0].vx *= -1;
        else if (i == N-1)
            dots[N-1].vx *= -1;
        else
        {
            double u = (dots[i].type.mass * dots[i].vx + dots[i+1].type.mass * dots[i+1].vx) 
                    / (dots[i].type.mass + dots[i+1].type.mass);
            dots[i].vx = 2*u - dots[i].vx;
            dots[i+1].vx = 2*u - dots[i+1].vx;
        }
    }

public:

    Gas(Box box, vector<Dot> dots): box(box)
    {
        sort(dots.begin(), dots.end(), [](Dot d1, Dot d2){return d1.x < d2.x;});
        this->dots = dots;
        N = dots.size();
        for (int i = 0; i < N; i++)
            Ns[dots[i].type.name] = 0;
        for (int i = 0; i < N; i++)
            Ns[dots[i].type.name]++;
    }

    // moves particles untill the next collision occures, including walls, and handle the collision
    // returns the time passed from start to next collision
    double move_untill_collide()
    {
        CollisionInfo collision = collision_search_l();
        if (collision.time == INFINITY)
            return -1;
        if (collision.time != 0)
            for (int i = 0; i < N; i++)
                dots[i].x += dots[i].vx * collision.time;
        collide(collision);
        return collision.time;
    }

    void print_state()
    {
        unordered_map<string, double> energy_summator;
        for (int i = 0; i < N; i++)
        {
            cout << dots[i].x << "," << dots[i].vx << ",";
            energy_summator[dots[i].type.name] += dots[i].type.mass * dots[i].vx * dots[i].vx / 2;
        }
        cout << " energy: ";
        for (auto energy : energy_summator)
            cout << energy.second / Ns[energy.first] << ",";
    }

    // checks if particles stays sorted, so collisions have been found and handled correct in some way
    bool check_if_valid()
    {
        for (int i = 0; i < N-1; i++)
            if (dots[i].x > dots[i+1].x)
                return false;
        if (dots[0].x < box.l)
            return false;
        if (dots[N-1].x > box.r)
            return false;
        return true;
    }
};

Gas load_gas(string path)
{
    ifstream f(path);
    nlohmann::json config = nlohmann::json::parse(f);

    unordered_map<string, DotType> types;
    for (auto type_j: config["types"])
    {
        DotType type;
        type.mass = type_j["mass"];
        type.radius = type_j["radius"];
        type.name = type_j["name"];
        types[type.name] = type;
    }

    Box box;
    box.l = config["box"]["l"];
    box.r = config["box"]["r"];

    vector<Dot> dots;
    for (auto dot_j : config["dots"])
    {
        Dot dot;
        dot.type = types[dot_j["type"]];
        dot.x = dot_j["x"];
        dot.vx = dot_j["vx"];
        dots.push_back(dot);
    }

    return Gas(box, dots);
}

int main(int argc, char** argv)
{
    Gas gas = load_gas("/Users/samedi/Documents/прога/study_modelling/term_2/gas/config.json");
    
    /*switch (*argv[1])
    {
    case 'n':
        for 
        break;
    case 'v'
        break;
    default:
        cout << "Stop condition is undefined. Please use " 
             << '"' << "n <natural number>" << '"' 
             << " if you want to limit the maximum amount of collisions, or " 
             << '"' << "t <real value>" << '"'
             << " if you want ti limit the maximum evaluation time.";
        break;
    }*/

    size_t N;
    if (argc > 1)
        N = atoi(argv[1]);
    else
        N = 0;

    cout << "collisions: " << 0 << " time: " << 0 << " state: ";
    gas.print_state();
    cout << "\n";
    for (int i = 0; i < N; i++)
    {
        auto time = gas.move_untill_collide();
        if (time == -1)
        {
            cout << "no more collisions" << endl;
            break;
        }
        cout << "collisions: " << i+1 << " time: " << time << " state: ";
        gas.print_state();
        cout << "\n";
    }
    cout << "control: " << (gas.check_if_valid());
}