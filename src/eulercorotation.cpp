#include "eulerworld.hpp"

Corotation::Corotation(double lambda, double mu) : lambda(lambda), mu(mu) {}

//Initalize Material. Decompose F for use in stress tensor
//Keep R, S, L
void Corotation::init(const MatrixN& F) {
    svd(F, U, Fhat, V);
    R = U * V.transpose();
    S = V * Fhat * V.transpose();

    MatrixN I = MatrixN::Identity();
    MatrixN Ld = (mu * I) + (lambda * (Fhat - I).trace() - 2 * mu) * ((Fhat.trace() * I) - Fhat).inverse();

    Ld(0, 0) = Ld(0, 0) > 0 ? Ld(0, 0) : 0;
    Ld(1, 1) = Ld(1, 1) > 0 ? Ld(1, 1) : 0;
    #if DIM == 3
    Ld(2, 2) = Ld(2, 2) > 0 ? Ld(2, 2) : 0;
    #endif
    L = V * Ld * V.transpose();
}

//Compute the Singular Value Decomposition of F into U * Fhat * VT
//https://scicomp.stackexchange.com/questions/8899/robust-algorithm-for-2x2-svd
void Corotation::svd(const MatrixN& F, MatrixN& U, MatrixN& Fhat, MatrixN& V) {
    // compute the SVD using a 2x2 specific algorithm
    double A = (F(0, 0) + F(1, 1)) / 2.0;
    double B = (F(0, 0) - F(1, 1)) / 2.0;
    double C = (F(1, 0) + F(0, 1)) / 2.0;
    double D = (F(1, 0) - F(0, 1)) / 2.0;

    double Q = std::sqrt(A*A + D*D);
    double R = std::sqrt(B*B + C*C);

    double sx = Q + R;
    double sy = Q - R;

    double a1, a2;
    //Safeguard against atan2(0, 0) = NaN
    if(std::abs(C) < 1e-12 && std::abs(B) < 1e-12) {
        a1 = 0;
    }
    else {
        a1 = std::atan2(C, B);
    }
    if(std::abs(D) < 1e-12 && std::abs(A) < 1e-12) {
        a2 = 0;
    }
    else {
        a2 = std::atan2(D, A);
    }

    double theta = (a2 - a1) / 2.0;
    double phi = (a2 + a1) / 2.0;

    //Sign of sy. 1 if >= 0, -1 otherwise
    int S = (sy >= 0.0) ? 1 : -1;
    //V is rotation by theta
    V(0, 0) = V(1, 1) = std::cos(theta);
    V(1, 0) = std::sin(theta);
    V(0, 1) = -V(1, 0);

    //U is a rotation by phi
    U(0, 0) = std::cos(phi);
    U(1, 1) = S * U(0, 0);
    U(1, 0) = std::sin(phi);
    U(0, 1) = -S * U(1, 0);

    Fhat.setZero();
    Fhat(0, 0) = sx;
    Fhat(1, 1) = std::abs(sy);


    //More standardized way
    // MatrixN Fnormal3 = F.transpose() * F;
    // Eigen::JacobiSVD<MatrixN> eigenSystem(Fnormal3, Eigen::ComputeFullV);
    // VectorN eigenvalues = eigenSystem.singularValues();
    // V = eigenSystem.matrixV();
    //
    // VectorN oEig = eigenvalues;
    // MatrixN oV = V;
    //
    // // populate the diagonal
    // Fhat.setZero();
    // Fhat(0, 0) = std::sqrt(eigenvalues(0));
    // Fhat(1, 1) = std::sqrt(eigenvalues(1));
    // #if DIM == 3
    // Fhat(2, 2) = std::sqrt(eigenvalues(2));
    // #endif
    //
    // if(std::isnan(Fhat(0, 0))) Fhat(0, 0) = 0.0;
    // if(std::isnan(Fhat(1, 1))) Fhat(1, 1) = 0.0;
    // #if DIM == 3
    // if(std::isnan(Fhat(2, 2))) Fhat(2, 2) = 0.0;
    // #endif
    //
    // // if V is a reflection, multiply a column by -1
    // // to ensure it is a rotation
    // double detV = V.determinant();
    // if (detV < 0.0) {
    //     V(0, 0) *= -1.0;
    //     V(1, 0) *= -1.0;
    //     #if DIM == 3
    //     V(2, 0) *= -1.0;
    //     #endif
    // }
    //
    // // compute U
    // U.setZero();
    // U(0, 0) = (Fhat(0, 0) > 0.0f) ? 1.0f / Fhat(0, 0) : 0.0f;
    // U(1, 1) = (Fhat(1, 1) > 0.0f) ? 1.0f / Fhat(1, 1) : 0.0f;
    // #if DIM == 3
    // U(2, 2) = (Fhat(2, 2) > 0.0f) ? 1.0f / Fhat(2, 2) : 0.0f;
    // #endif
    // U = F * V * U;
    // orthogonalizeU(U, Fhat);
    //
    // // correct any reflections in U to ensure it is a rotation
    // if (F.determinant() < 0.0) removeUReflection(U, Fhat);
}

//First Piola Kirchhoff Stress Tensor
MatrixN Corotation::firstPiolaKirchhoff() {
    MatrixN E = S - MatrixN::Identity();
    MatrixN tmp = 2.0 * mu * E;
    double tr = lambda * E.trace();
    tmp(0, 0) += tr;
    tmp(1, 1) += tr;
    #if DIM == 3
    tmp(2, 2) += tr;
    #endif
    return R * tmp;
}

MatrixN Corotation::firstPKDifferential(const MatrixN& dF) {
    MatrixN dFhat = R.transpose() * dF;
    MatrixN dFsym = 0.5 * (dFhat + dFhat.transpose());
    #if DIM == 3
    MatrixN dFskew = 0.5 * (dFhat - dFhat.transpose());
    #endif

    MatrixN dPsym = (2 * mu * dFsym) + (lambda * dFsym.trace() * MatrixN::Identity());

    #if DIM == 3
    VectorN f(-dFskew(1, 2) + dFskew(2, 1), -dFskew(2, 0) + dFskew(0, 2), -dFskew(0, 1) + dFskew(1, 0));
    #else
    // VectorN f(-dFskew(1, 0) + dFskew(0, 1), dFskew(1, 0) - dFskew(0, 1));
    #endif

    double tr = lambda * dFhat.trace();
    dPsym(0, 0) += tr;
    dPsym(1, 1) += tr;
    #if DIM == 3
    dPsym(2, 2) += tr;
    #endif

    #if DIM == 3
    MatrixN dPskew = crossMatrix(L*f);
    #else
    MatrixN dPskew = MatrixN::Zero();
    #endif
    //z axis is not doing anything so the skew would be 0.
    //any part of dFskew(?,2) and dFskew(2,?) is going to be 0. f ends up with (0,0,x)
    //L is 0 in 3rd dimention, so L*f is 0. cross product matrix is thus a matrix of 0
    //Adding 0 doesn't change anything
    MatrixN deltaP = dPsym + dPskew;
    return R * deltaP;
}

#if DIM == 3
MatrixN Corotation::crossMatrix(const VectorN& vec) {
    MatrixN final;
    final(0, 0) = 0.0;
    final(1, 0) = vec[2];
    final(2, 0) = -vec[1];

    final(0, 1) = -vec[2];
    final(1, 1) = 0.0;
    final(2, 1) = vec[0];

    final(0, 2) = vec[1];
    final(1, 2) = -vec[0];
    final(2, 2) = 0.0;

    return final;
}
#endif
