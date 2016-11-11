#include "grid.hpp"
#include <cmath>

template <typename T>
T& Grid<T>::get(int i, int j) {
    return values[i+m*j];
}

template <typename T>
void Grid<T>::assign(T value) {
    values.assign(values.size(), value);
}

template <typename T>
T Grid<T>::interpolate(Vector2d x) {
    T val = T(0);
    for(int i = lower(x(0), 0); i < upper(x(0), 0); i++) {
        for(int j = lower(x(1), 1); j < upper(x(1), 1); j++) {
            val += get(i, j) * weight(x, i, j);
        }
    }
    return val;
}

//Need special case for eigen Vector because of initialization to 0
template<>
Vector2d Grid<Vector2d>::interpolate(Vector2d x) {
    Vector2d val = Vector2d::Zero();
    for(int i = lower(x(0), 0); i < upper(x(0), 0); i++) {
        for(int j = lower(x(1), 1); j < upper(x(1), 1); j++) {
            val += get(i, j) * weight(x, i, j);
        }
    }
    return val;
}

template <typename T>
double Grid<T>::weight(Vector2d x, int i, int j) {
    return N(hinv*(x(0)-i*h)) * N(hinv*(x(1)-j*h));
}

template <typename T>
Vector2d Grid<T>::gradWeight(Vector2d x, int i, int j) {
    Vector2d grad(0, 0);
    grad(0) = hinv*dN(hinv*(x(0)-i*h))*N(hinv*(x(1)-j*h));
    grad(1) = hinv*N(hinv*(x(0)-i*h))*dN(hinv*(x(1)-j*h));
    return grad;
}

template <typename T>
double Grid<T>::N(double x) {
    double ax = abs(x);
    if(ax >= 0 && ax < 1) {
        return 0.5*ax*ax*ax - ax*ax + (2.0/3.0);
    }
    else if(ax >= 1 && ax < 2) {
        return (-1.0/6.0)*ax*ax*ax + ax*ax - 2.0*ax + (4.0/3.0); 
    }
    return 0;
}

template <typename T>
double Grid<T>::dN(double x) {
    double ax = abs(x);
    if(ax >= 0.0 && ax < 1.0) {
        return 1.5*x*ax - 2.0*x;
    }
    else if(ax >= 1.0 && ax < 2.0) {
        return (-0.5)*x*ax + 2.0*x - 2.0*(x/ax); 
    }
    return 0;
}

template <typename T>
int Grid<T>::lower(double x, int axis) {
    return std::max(0.0, std::floor(-2.0 + x/h));
}

template <typename T>
int Grid<T>::upper(double x, int axis) {
    double bounds = (axis == 0) ? m : n;
    return std::min(bounds, std::ceil(2.0 + x/h));
}

template <typename T>
void Grid<T>::addInterpolated(Vector2d x, T value) {
    for(int i = lower(x(0), 0); i < upper(x(0), 0); i++) {
        for(int j = lower(x(1), 1); j < upper(x(1), 1); j++) {
            get(i, j) += weight(x, i, j) * value;
        }
    }
}

void StaggeredGrid::assign(Vector2d value) {
    u.assign(value(0));
    v.assign(value(1));
}

Vector2d StaggeredGrid::interpolate(Vector2d x) {
    return Vector2d(u.interpolate(x), v.interpolate(x));
}

void StaggeredGrid::addInterpolated(Vector2d x, Vector2d value) {
    u.addInterpolated(x, value(0));
    v.addInterpolated(x, value(1));
}

template class Grid<char>;
template class Grid<int>;
template class Grid<double>;
template class Grid<Vector2d>;
