//Testing interpolation
//g++ $(pkg-config --cflags eigen3) interpTest.cpp -o interpTest
#include <Eigen/Dense>
#include <cstdio>

using namespace Eigen;

inline Vector2d interpolate(Vector2d point, Vector2d* field, Vector2d origin, int res[2], double h) {
    /// wx = x - x1;
    /// wy = y - y1;
    /// v = (1-wx)*(1-wy)*v[i1][j1] + (wx)*(1-wy)*v[i1+1][j1] + (1-wx)*(wy)*v[i1][j1+1] + (wx)*(wy)*v[i1+1][j1+1];
    Vector2d x = (point - origin);
    Vector2d ij = x / h;    
    int i1 = (int)ij(0);
    int j1 = (int)ij(1);
    int i2 = i1+1;
    int j2 = j1+1;
    double x1 = h*i1;
    double y1 = h*j1;
    double wx = (x(0) - x1) / h;
    double wy = (x(1) - y1) / h;
    printf("x1, y1: %f, %f\n", x1, y1);
    printf("wx, wy: %f, %f\n", wx, wy);
    Vector2d v = (1-wx)*(1-wy)*field[i1*res[1]+j1] + 
                 (wx)  *(1-wy)*field[i2*res[1]+j1] + 
                 (1-wx)*(wy)  *field[i1*res[1]+j2] + 
                 (wx)  *(wy)  *field[i2*res[1]+j2];
    
    return v;
}

int main(int argc, char** argv) {
    int res[2] = {5, 5};
    Vector2d origin(-1, 0);
    double h = 0.5;
    Vector2d *field = new Vector2d[res[0]*res[1]];
    for(int i = 0; i < res[0]; i++) {
        for(int j = 0; j < res[1]; j++) {
            /// field[i*res[1]+j] = origin + h*Vector2d(i,j);
            field[i*res[1]+j] = Vector2d(i,j);
        }
    }
    
    for(int j = res[1]-1; j >= 0; j--) {
        for(int i = 0; i < res[0]; i++) {
            printf("(%f, %f) ", field[i*res[1]+j](0), field[i*res[1]+j](1));
        }
        printf("\n");
    }
    
    Vector2d p1(0,1);
    Vector2d p2(0.75, 0.75);
    Vector2d p3(-0.2, 0.6);
    
    Vector2d i1(0,0), i2(0,0), i3(0,0);
    i1 = interpolate(p1, field, origin, res, h);
    i2 = interpolate(p2, field, origin, res, h);
    i3 = interpolate(p3, field, origin, res, h);
    
    printf("Value at (%f, %f): (%f, %f)\n", p1(0), p1(1), i1(0), i1(1));
    printf("Value at (%f, %f): (%f, %f)\n", p2(0), p2(1), i2(0), i2(1));
    printf("Value at (%f, %f): (%f, %f)\n", p3(0), p3(1), i3(0), i3(1));
    
    return 0;
}
