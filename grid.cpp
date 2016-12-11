#include "grid.hpp"
#include "defines.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace Eigen;

template <typename T>
T& Grid<T>::get(int i, int j) {
    return values[i+m*j];
}

template <typename T>
void Grid<T>::assign(T value) {
    values.assign(values.size(), value);
}

//Test weight. Put positions of grid cells in world space on the grid and interpolate
//to particles and see if the particle position is recovered
//Ask if higher order interpolation should indeed recover the exact values
template <typename T>
double Grid<T>::weight(Vector2d x, int i, int j) {
    double w = N(hinv*(x(0)-i*h-h/2.0)) * N(hinv*(x(1)-j*h-h/2.0));
    return w;
}

template <typename T>
Vector2d Grid<T>::gradWeight(Vector2d x, int i, int j) {
    Vector2d grad(0.0, 0.0);
    grad(0) = hinv*dN(hinv*(x(0)-i*h))*N(hinv*(x(1)-j*h));
    grad(1) = hinv*N(hinv*(x(0)-i*h))*dN(hinv*(x(1)-j*h));
    return grad;
}

template <typename T>
double Grid<T>::N(double x) {
    double ax = std::abs(x);
    if(ax >= 0.0 && ax < 1.0) {
        return 0.5*ax*ax*ax - x*x + (2.0/3.0);
    }
    else if(ax >= 1.0 && ax < 2.0) {
        return (-1.0/6.0)*ax*ax*ax + x*x - 2.0*ax + (4.0/3.0);
    }
    return 0;
}

template <typename T>
double Grid<T>::dN(double x) {
    double ax = std::abs(x);
    if(ax >= 0.0 && ax < 1.0) {
        return 1.5*x*ax - 2.0*x;
    }
    else if(ax >= 1.0 && ax < 2.0) {
        return (-0.5)*x*ax + 2.0*x - 2.0*(x/ax);
    }
    return 0;
}

template <typename T>
int Grid<T>::lower(Vector2d x, int axis) {
    /// return std::max(0, (int)std::ceil(-2.0 + x(axis)/h));
    return std::max(0, (int)std::ceil(-2.0 + x(axis)/h) - 1);
}

template <typename T>
int Grid<T>::upper(Vector2d x, int axis) {
    int bounds = (axis == 0) ? m : n;
    /// return std::min(bounds, (int)std::floor(2.0 + x(axis)/h));
    return std::min(bounds, (int)std::floor(2.0 + x(axis)/h) + 1);
}

template class Grid<char>;
template class Grid<int>;
template class Grid<double>;
template class Grid<Vector2d>;
