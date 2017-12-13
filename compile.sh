#!/bin/bash
#Check if argument 1 is debug or release, set appropriate build type
if [ -z ${1+x} ]
then
    BUILD_TYPE=-DCMAKE_BUILD_TYPE=Release
#elif [ ${1,,} == "debug" ]
elif [ $(echo "${1}" | tr '[:upper:]' '[:lower:]') == "debug" ]
then
    BUILD_TYPE=-DCMAKE_BUILD_TYPE=Debug
fi
mkdir -p build
cd build
unameOut="$(uname -s)"
if [ $unameOut == "Darwin" ]
then
    export CC=/usr/local/bin/gcc-7
    export CXX=/usr/local/bin/g++-7
fi
cmake -Wno-dev -DCMAKE_PREFIX_PATH=/usr/local/ "$BUILD_TYPE" ..
make VERBOSE=1
cd ..
