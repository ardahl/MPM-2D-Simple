#include "grid.hpp"

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
    Vector2d index = (x - x0)/dx;
    int i = floor(index(0)), j = floor(index(1)); // integer part of index
    double fi = index(0) - i, fj = index(1) - j;  // fractional part
    if (i < 0)
        {i = 0; fi = 0;}
    else if (i >= m-1)
        {i = m-2; fi = 1;}
    if (j < 0)
        {j = 0; fj = 0;}
    else if (j >= n-1)
        {j = n-2; fj = 1;}
    return
        + (1-fi)*(1-fj)*get(i,   j)
        + (  fi)*(1-fj)*get(i+1, j)
        + (1-fi)*(  fj)*get(i,   j+1)
        + (  fi)*(  fj)*get(i+1, j+1);
}

template <typename T>
void Grid<T>::addInterpolated(Vector2d x, T value) {
    Vector2d index = (x - x0)/dx;
    int i = floor(index(0)), j = floor(index(1)); // integer part of index
    double fi = index(0) - i, fj = index(1) - j;  // fractional part
    if (i < 0)
        {i = 0; fi = 0;}
    else if (i >= m-1)
        {i = m-2; fi = 1;}
    if (j < 0)
        {j = 0; fj = 0;}
    else if (j >= n-1)
        {j = n-2; fj = 1;}
    get(i,   j)   += (1-fi)*(1-fj)*value;
    get(i+1, j)   += (  fi)*(1-fj)*value;
    get(i,   j+1) += (1-fi)*(  fj)*value;
    get(i+1, j+1) += (  fi)*(  fj)*value;
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
