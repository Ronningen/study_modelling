#include <iostream>
#include <fstream>
#include <vector>
#include <valarray>
#include "json.hpp"

using namespace std;
using namespace nlohmann;

using vec = valarray<double>;

extern "C" void dgesv( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
#define N 2
#define NRHS 1
#define LDA N
#define LDB N

//v1x*v2x+v1y+v2y
double dot(vec v1, vec v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1];
}

struct Box
{
    double l,r,b,t;
    Box(json j):l(j["l"]),r(j["r"]),b(j["b"]),t(j["t"]){}
};

struct RayPart
{
    int from;
    vec x, v;
    RayPart(json j):x({j["x0"], j["y0"]}), v({j["vx"], j["vy"]}), from(-1){}
    RayPart(vec x, vec v, int from): x(x), v(v), from(from){}
};

struct Ray
{
    int tag;
    vector<RayPart> parts;
    Ray(json j):parts({RayPart(j)}), tag(j["tag"]){}
};

struct IMirror
{
    int tag;
    IMirror(int tag):tag(tag){}
    vec virtual Trace(Ray, bool*) = 0;
    RayPart virtual Reflect(Ray,vec) = 0;
};

struct Line : IMirror
{
    vec x1, x2, n;
    Line(json j):x1({j["x1"],j["y1"]}), x2({j["x2"],j["y2"]}), IMirror(j["tag"])
    {
        auto v = x2 - x1;
        n = vec({v[1],v[0]}) / sqrt(dot(v,v));
    }
    vec Trace(Ray ray, bool* has_result) final
    {
        auto p1 = ray.parts.back().v;
        auto p2 = x2 - x1;
        auto dx = x1 - ray.parts.back().x;
        
        int n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
        int ipiv[N];
        double a[LDA*N] = {
            p1[0],p2[0],
            p1[1],p2[1]
        };
        double b[LDB*NRHS] = {
            dx[0], dx[1]
        };
        dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
        double t = -b[1];

        /*TODO: strange values...*/ cout << b[0] << " " << t << " ";

        *has_result = !(info > 0 or t < 0 or t > 1 or b[0] <= 0);
        if (*has_result)
            return x1 + p2*t;
        else
            return {};
    }
    RayPart Reflect(Ray ray, vec at) final
    {
        vec v = ray.parts.back().v;
        v -= 2*dot(v,n)*n;
        return RayPart(at, v, tag);
    }
};

//TODO: Arc
struct Arc : IMirror
{
    vec x, r;
    double a, r2;
    Arc(json j):x({j["x"],j["y"]}), r({j["rx"],j["ry"]}), a(j["a"]), IMirror(j["tag"])
    {
        r2 = dot(r,r);
    }
    vec Trace(Ray ray, bool* has_result) final
    {
        auto v = ray.parts.back().v;
        auto x0 = ray.parts.back().x;
        auto a = dot(v, v);
        auto b = dot(v, x0-x);
        auto c = dot(x0-x, x0-x) - r2;
        auto D = b*b - a*c;
        if (D <= 0)
            return {};
        else
            return x0 + v*(-b-sqrt(D))/a;
    }
    RayPart Reflect(Ray ray, vec at) final
    {
        vec v = ray.parts.back().v;
        vec n = (at - x) / sqrt(r2);
        v -= 2*dot(v,n)*n;
        return RayPart(at, v, tag);
    }
};


int main(int argc, char **argv)
{
    // parse
    ifstream f("term_2/reflector/config.json");
    json config = json::parse(f);

    vector<IMirror*> mirrors{};
    for (auto line : config["lines"])
        mirrors.push_back(new Line(line));
    for (auto arc : config["arcs"])
        mirrors.push_back(new Line(arc));

    vector<Ray*> rays{};
    for (auto ray : config["rays"])
        rays.push_back(new Ray(ray));

    // trace
    bool tracing = true;
    bool has_result;
    while (tracing)
    {
        tracing = false;
        for (auto ray : rays)
        {
            vec intersection{{0,0}};
            double min_distance2 = -1;
            IMirror* nearest_mirror = nullptr;

            for (auto mirror_ : mirrors)
            {
                vec x = mirror_->Trace(*ray, &has_result);
                if (has_result)
                {
                    auto l = intersection - x;
                    double d2 = dot(l,l);
                    if (d2 < min_distance2 || min_distance2 < 0)
                    {
                        tracing = true;
                        min_distance2 = d2;
                        intersection = x;
                        nearest_mirror = mirror_;
                    }
                }
            }

            if (nearest_mirror)
                ray->parts.push_back(nearest_mirror->Reflect(*ray, intersection));
        }
    }

    // print
    for (auto ray : rays)
    {
        cout << endl << " ray ";
        for (auto part : ray->parts)
            cout << part.x[0] << " " << part.x[1] << " ";
        cout << ray->parts.back().v[0]*10000 << " " << ray->parts.back().v[1]*10000;
    }

    for (auto ray : rays)
        delete ray;
}