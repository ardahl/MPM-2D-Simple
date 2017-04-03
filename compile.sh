#!/bin/bash
#Check if argument 1 is debug or release, set appropriate build type
if [ -z ${1+x} ]
then
    BUILD_TYPE=-DCMAKE_BUILD_TYPE=Release
elif [ ${1,,} == "debug" ]
then
    BUILD_TYPE=-DCMAKE_BUILD_TYPE=Debug
fi
mkdir -p build
cd build
cmake -Wno-dev -DCMAKE_PREFIX_PATH=/usr/local/ "$BUILD_TYPE" ..
make
cd ..
