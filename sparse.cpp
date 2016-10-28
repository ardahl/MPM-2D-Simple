#include "sparse.hpp"

void SpMatrix::add(int i, int j, double value) {
    triplets.push_back(Triplet<double>(i,j, value));
}

VectorXd SpMatrix::solve(const VectorXd &b) {
    SparseMatrix<double> A(m, n);
    A.setFromTriplets(triplets.begin(), triplets.end());
    ConjugateGradient< SparseMatrix<double> > solver;
    solver.setTolerance(1e-6);
    solver.setMaxIterations(200);
    return solver.compute(A).solve(b);
}
