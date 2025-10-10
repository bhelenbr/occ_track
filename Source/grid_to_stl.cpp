//
//  grid_to_stl.cpp
//  occ_track
//
//  Created by Brian Helenbrook on 10/10/25.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include "occ_track.h"


static Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

static Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static Vec3 normalize(const Vec3& v) {
    double len = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (len < 1e-12) return {0,0,0};
    return {v.x/len, v.y/len, v.z/len};
}

void writeSTL(const std::vector<std::vector<Vec3>>& grid,
              const std::string& filename)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    out << "solid grid\n";
    int nx = static_cast<int>(grid.size());
    int ny = static_cast<int>(grid[0].size());

    for (int i = 0; i < nx - 1; ++i) {
        for (int j = 0; j < ny - 1; ++j) {
            Vec3 p1 = grid[i][j];
            Vec3 p2 = grid[i+1][j];
            Vec3 p3 = grid[i][j+1];
            Vec3 p4 = grid[i+1][j+1];

            // Triangle 1: p1, p2, p3
            {
                Vec3 n = normalize(cross(p2 - p1, p3 - p1));
                out << "  facet normal " << n.x << " " << n.y << " " << n.z << "\n";
                out << "    outer loop\n";
                out << "      vertex " << p1.x << " " << p1.y << " " << p1.z << "\n";
                out << "      vertex " << p2.x << " " << p2.y << " " << p2.z << "\n";
                out << "      vertex " << p3.x << " " << p3.y << " " << p3.z << "\n";
                out << "    endloop\n";
                out << "  endfacet\n";
            }

            // Triangle 2: p2, p4, p3
            {
                Vec3 n = normalize(cross(p4 - p2, p3 - p2));
                out << "  facet normal " << n.x << " " << n.y << " " << n.z << "\n";
                out << "    outer loop\n";
                out << "      vertex " << p2.x << " " << p2.y << " " << p2.z << "\n";
                out << "      vertex " << p4.x << " " << p4.y << " " << p4.z << "\n";
                out << "      vertex " << p3.x << " " << p3.y << " " << p3.z << "\n";
                out << "    endloop\n";
                out << "  endfacet\n";
            }
        }
    }

    out << "endsolid grid\n";
    std::cout << "Wrote STL file: " << filename << std::endl;
}

//int main() {
//    // Example: create a simple structured grid (z = sin(x)*cos(y))
//    int nx = 50, ny = 50;
//    double xmin = -2.0, xmax = 2.0;
//    double ymin = -2.0, ymax = 2.0;
//
//    std::vector<std::vector<Vec3>> grid(nx, std::vector<Vec3>(ny));
//
//    for (int i = 0; i < nx; ++i) {
//        for (int j = 0; j < ny; ++j) {
//            double x = xmin + (xmax - xmin) * i / (nx - 1);
//            double y = ymin + (ymax - ymin) * j / (ny - 1);
//            double z = std::sin(x) * std::cos(y);
//            grid[i][j] = {x, y, z};
//        }
//    }
//
//    writeSTL(grid, "surface.stl");
//    return 0;
//}
