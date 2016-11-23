#ifndef GRID_HPP_
#define GRID_HPP_

#include <Eigen/Dense>
#include <vector>

template <typename T>
class Grid {
public:
    int m, n;         // number of grid cells
    Eigen::Vector2d x0;      // world-space position of cell (0, 0)
    double h, hinv;        // grid cell spacing and inverse
    std::vector<T> values; // grid cell values
    Grid() {}
    Grid(int m, int n, Eigen::Vector2d x0, double dx):
        m(m), n(n), x0(x0), h(dx), hinv(1.0/dx), values(m*n) {}
    // value at cell (i, j)
    T& get(int i, int j);
    T& operator()(int i, int j) {
        return get(i, j);
    }
    //output
    friend std::ostream& operator<<(std::ostream &out, const Grid &grid) {
        for(int i = 0; i < grid.m; i++) {
            for(int j = 0; j < grid.n; j++) {
                out << grid.values[i+grid.m*j] << "    ";
            }
            out << "\n";
        }
        return out;
    }
    // assign all values
    void assign(T value);
    // interpolated value at world-space position x
    T interpolate(Eigen::Vector2d x);
    // accumulate value at world-space position x onto grid
    void addInterpolated(Eigen::Vector2d x, T value);
    // From Stomakhin et al. 2013 for use in interpolation
    double weight(Eigen::Vector2d x, int i, int j);
    Eigen::Vector2d gradWeight(Eigen::Vector2d x, int i, int j);
    double N(double x);
    double dN(double x);
    // Gets the index of the lower bounds of the grid
    int lower(double x, int axis);
    int upper(double x, int axis);
};

class StaggeredGrid {
public:
    int m, n;          // number of grid cells
    Eigen::Vector2d x0;       // world-space position of cell (0, 0)
    double dx;         // grid cell spacing
    Grid<double> u, v; // values at horizontal and vertical cell faces
    StaggeredGrid() {}
    StaggeredGrid(int m, int n, Eigen::Vector2d x0, double dx):
        m(m), n(n), x0(x0), dx(dx),
        u(m+1, n, x0-Eigen::Vector2d(dx/2,0), dx),
        v(m, n+1, x0-Eigen::Vector2d(0,dx/2), dx) {}
    // assign all values
    void assign(Eigen::Vector2d value);
    // interpolated value at world-space position x
    Eigen::Vector2d interpolate(Eigen::Vector2d x);
    // accumulate value at world-space position x onto grid
    void addInterpolated(Eigen::Vector2d x, Eigen::Vector2d value);
};

#endif
