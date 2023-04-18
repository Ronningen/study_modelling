#include <iostream>
#include <vector>
#include <valarray>
#include "json.hpp"

using namespace std;
using namespace nlohmann;

using vec = valarray<double>;

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
    RayPart(json j):x({j["x0"], j["y0"]}), v({j["rx"], j["ry"]}), from(-1){}
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
    vec virtual Trace(Ray) = 0;
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
    vec Trace(Ray ray) final
    {
        auto v = ray.parts.back().v;
        auto x0 = ray.parts.back().x;
        auto u = x2 - x1;
        auto t = ((x1[0]-x0[0])*v[1] - (x1[1]-x0[1])*v[0]) / (u[0]*v[1] - u[1]*v[0]);
        if (t < 0 or t > 1)
            return {};
        else
            return x1 + u*t;
    }
    RayPart Reflect(Ray ray, vec at) final
    {
        vec v = ray.parts.back().v;
        v -= 2*dot(v,n)*n;
        return RayPart(at, v, tag);
    }
};

struct Arc : IMirror
{
    vec x, r;
    double a, r2;
    Arc(json j):x({j["x"],j["y"]}), r({j["rx"],j["ry"]}), a(j["a"]), IMirror(j["tag"])
    {
        r2 = dot(r,r);
    }
    vec Trace(Ray ray) final
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

}