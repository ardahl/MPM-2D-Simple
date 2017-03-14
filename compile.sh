#!/bin/bash
mkdir -p build
cd build
cmake -Wno-dev -DCMAKE_PREFIX_PATH=/usr/local/ ..
make
cd ..
