#ifndef SPARSE_HPP_
#define SPARSE_HPP_

#include <Eigen/Sparse>
#include <vector>
using namespace Eigen;

// helper class for incremental construction of sparse matrices
class SpMatrix {
public:
    int m, n;
    std::vector< Triplet<double> > triplets;
    SpMatrix(int m, int n): m(m), n(n) {}
    void add(int i, int j, double value);
    void addBlock(int i, int j, int p, int q, MatrixXd block);
    VectorXd solve(const VectorXd &b);
};

#endif
