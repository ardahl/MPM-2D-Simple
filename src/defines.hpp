/**
 * The purpose of this file is for somewhere to put all of the global constants
 * and defines that are used. Rather than have them scattered throughout different
 * files, it's a bit easier to find if they're in one place and are properly
 * documented. Make sure everything is documented properly.
 *
 */

#ifndef DEFINES_HPP_
#define DEFINES_HPP_

//Uncomment this to remove debug file output
/// Code for debugging should be put around a #ifndef NDEBUG guard. The define
/// below should be commented if you want the debug output, and uncommented to
/// remove debug output.
#define NDEBUG

/// Currently the value of what's close enough to 0 to round down to 0 for
/// numerical calculation purposes. This may be a function of the grid later
/// in which case it will need to move. But for now it'll rest here.
const double EPS = 1e-5;

#endif
