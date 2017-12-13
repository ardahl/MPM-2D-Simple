//OSX: g++-7 -std=c++11 -Wall -g -c partioTest.cpp -o partioTest.o
//     g++-7 -std=c++11 -Wall -g partioTest.o -o partioTest /usr/local/lib/libpartio.dylib
#include <Partio.h>
#include <iostream>

int main(int argc, char** argv) {
    Partio::ParticlesDataMutable *data = Partio::read("test.bgeo");
    std::cout << "Number of particles " << data->numParticles() << std::endl;

    return 0;
}
