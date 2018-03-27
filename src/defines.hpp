/**
 * The purpose of this file is for somewhere to put all of the global constants
 * and defines that are used. Rather than have them scattered throughout different
 * files, it's a bit easier to find if they're in one place and are properly
 * documented. Make sure everything is documented properly.
 *
 */

#ifndef DEFINES_HPP_
#define DEFINES_HPP_

#include <Eigen/Dense>
#include <Eigen/Sparse>

// Regular MPM or testing material transfer
#define MAT_TRANSFER
// Splotting material differences or material coordiantes
#define DIFF
// Velocity extrapolation by fast sweeping or simple average
#define FASTSWEEP
// Use APIC for velocity transfer
#define APICVEL

/// Currently the value of what's close enough to 0 to round down to 0 for
/// numerical calculation purposes. This may be a function of the grid later
/// in which case it will need to move. But for now it'll rest here.
const double EPS = 5e-7;

// Typedefs for Eigen
// This allows for easier swaping of framework and dimention number
#define DIM 2
#if DIM == 2
typedef Eigen::Vector2d VectorN;
typedef Eigen::Vector2i VectorI;
typedef Eigen::Matrix2d MatrixN;
typedef Eigen::Array2d ArrayN;
#define QUAD 4                     //number of quadrature points
#else
typedef Eigen::Vector3d VectorN;
typedef Eigen::Vector3i VectorI;
typedef Eigen::Matrix3d MatrixN;
typedef Eigen::Array3d ArrayN;
#define QUAD 8
#endif

typedef Eigen::VectorXd VectorX;
typedef Eigen::MatrixXd MatrixX;
typedef Eigen::ArrayXd ArrayX;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMat;
typedef Eigen::Triplet<double> ETriplet;

#endif
