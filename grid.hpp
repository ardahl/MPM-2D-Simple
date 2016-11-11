#ifndef GRID_HPP_
#define GRID_HPP_

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;

template <typename T>
class Grid {
public:
    int m, n;         // number of grid cells
    Vector2d x0;      // world-space position of cell (0, 0)
    double h, hinv;        // grid cell spacing and inverse
    std::vector<T> values; // grid cell values
    Grid() {}
    Grid(int m, int n, Vector2d x0, double dx):
        m(m), n(n), x0(x0), h(dx), hinv(1.0/dx), values(m*n) {}
    // value at cell (i, j)
    T& get(int i, int j);
    T& operator()(int i, int j) {
        return get(i, j);
    }
    // assign all values
    void assign(T value);
    // interpolated value at world-space position x
    T interpolate(Vector2d x);
    // accumulate value at world-space position x onto grid
    void addInterpolated(Vector2d x, T value);
    // From Stomakhin et al. 2013 for use in interpolation
    double weight(Vector2d x, int i, int j);
    Vector2d gradWeight(Vector2d x, int i, int j);
    double N(double x);
    double dN(double x);
    // Gets the index of the lower bounds of the grid
    int lower(double x, int axis);
    int upper(double x, int axis);
};

class StaggeredGrid {
public:
    int m, n;          // number of grid cells
    Vector2d x0;       // world-space position of cell (0, 0)
    double dx;         // grid cell spacing
    Grid<double> u, v; // values at horizontal and vertical cell faces
    StaggeredGrid() {}
    StaggeredGrid(int m, int n, Vector2d x0, double dx):
        m(m), n(n), x0(x0), dx(dx),
        u(m+1, n, x0-Vector2d(dx/2,0), dx),
        v(m, n+1, x0-Vector2d(0,dx/2), dx) {}
    // assign all values
    void assign(Vector2d value);
    // interpolated value at world-space position x
    Vector2d interpolate(Vector2d x);
    // accumulate value at world-space position x onto grid
    void addInterpolated(Vector2d x, Vector2d value);
};

#endif
