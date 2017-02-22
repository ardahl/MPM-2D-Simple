/**
 * The purpose of this file is for somewhere to put all of the global constants
 * and defines that are used. Rather than have them scattered throughout different
 * files, it's a bit easier to find if they're in one place and are properly
 * documented. Make sure everything is documented properly.
 *
 */

#ifndef DEFINES_HPP_
#define DEFINES_HPP_

#include <cmath>

/// Currently the value of what's close enough to 0 to round down to 0 for
/// numerical calculation purposes. This may be a function of the grid later
/// in which case it will need to move. But for now it'll rest here.
const double EPS = 1e-5;

/// For the weight function, it is only close to 0 when 1<x<2. Solving that 
/// for x when it equals the epsilon value gives -+2 +- 6^(1/3)*EPS^(1/3). 
/// Substituting in pos/h - i for x and solving for i gives us the upper and
/// lower bounds where anything past them is less than epsilon.  
/// Cube root of 6 times cube root of our epsilon.
const double BOUNDOFFSET = 1.81712059283 * cbrt(EPS);

#endif
