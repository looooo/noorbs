#include <iostream>
#include <vector>

#include "unwrap.h"

typedef std::array<long, 3> triangle;
typedef std::array<double, 3> point;

int main(int argc, char **argv) {
    // cerate vertices
    point p0 {-1, -1, 0};
    point p1 {1, -1, 0};
    point p2 {1, 1, 0};
    point p3 {-1, 1, 0};
    point p4 {0, 0, 1};

    std::vector<point> vertices = {p0, p1, p2, p3, p4};

    // create triangles
    std::vector<triangle> triangles;
    triangles.push_back(triangle{ {0, 1, 4} });
    triangles.push_back(triangle{ {1, 2, 4} });
    triangles.push_back(triangle{ {2, 3, 4} });
    triangles.push_back(triangle{ {3, 0, 4} });


    nurbs::LscmRelax flattener(vertices, triangles, std::vector<long>());
    flattener.lscm();
    flattener.relax(0.001);
    std::cout << flattener.flat_vertices << std::endl;
}
