#include <Eigen/Dense>
#include <ctime>
//g++ $(pkg-config --cflags eigen3) eigenTest.cpp -o eigenTest

//Create a large dynamic array of 2d vectors and a 2*N dim dynamic vector
//Test speed of accessing in a row and random values
int main(int argc, char** argv) {
    int dim = 10000;
    Eigen::Vector2d* v1 = new Eigen::Vector2d[dim];
    Eigen::VectorXd v2(2*dim);
    
    for(int i = 0; i < dim; i++) {
        //Fill vectors with test data
        v1[i] = Eigen::Vector2d(i,i);
        v2.segment(2*i,2) = Eigen::Vector2d(i,i);
    }
    
    clock_t time;
    Eigen::Vector2d tmp;
    
    //Time accessing each element and assigning it to a 2d vector
    printf("Sequential Access\n");
    //Pointer First
    time = clock();
    for(int i = 0; i < dim; i++) {
        tmp = v1[i];
        v1[i] = tmp;
    }
    time = clock() - time;
    printf("Pointer: %d\n", (double)time/CLOCKS_PER_SEC);
    //Large Array
    time = clock();
    for(int i = 0; i < dim; i++) {
        tmp = v2.segment(2*i,2);
        v2.segment(2*i,2) = tmp;
    }
    time = clock() - time;
    printf("Vector: %d\n", (double)time/CLOCKS_PER_SEC);
    
    //Time accessing random locations
    printf("\nRandom Access\n");
    
    return 0;
}
