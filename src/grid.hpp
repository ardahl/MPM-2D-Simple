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
                out << grid.values[i*grid.n+j] << "    ";
            }
            out << "\n";
        }
        return out;
    }
    // assign all values
    void assign(T value);
    // From Stomakhin et al. 2013 for use in interpolation
    double weight(Eigen::Vector2d x, int i, int j);
    Eigen::Vector2d gradWeight(Eigen::Vector2d x, int i, int j);
    double N(double x);
    double dN(double x);
    // Gets the index of the lower bounds of the grid
    int lower(Eigen::Vector2d x, int axis);
    int upper(Eigen::Vector2d x, int axis);
};

#endif
