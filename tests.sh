#!/bin/bash

mkdir euler/none
#compile and run
echo "Running None"
./bin/euler input_files/euler_lin.json euler/euler > euler/output.txt
#move all none-folder items to the test directory
cd euler
for x in *; do
    if ! [ -d "$x" ]; then
        mv -- "$x" none/
    fi
done
cd ..

mkdir euler/linear
echo "Running Linear"
./bin/euler input_files/euler_lin.json euler/euler > euler/output.txt
cd euler
for x in *; do
    if ! [ -d "$x" ]; then
        mv -- "$x" linear/
    fi
done
cd ..

mkdir euler/rotation
echo "Running Rotation"
./bin/euler input_files/euler_rot.json euler/euler > euler/output.txt
cd euler
for x in *; do
    if ! [ -d "$x" ]; then
        mv -- "$x" rotation/
    fi
done
cd ..

mkdir euler/gravity
echo "Running Gravity"
./bin/euler input_files/euler_grav.json euler/euler > euler/output.txt
cd euler
for x in *; do
    if ! [ -d "$x" ]; then
        mv -- "$x" gravity/
    fi
done
cd ..
