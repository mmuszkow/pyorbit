#!/bin/bash

# Tested on OSX High Sierra version

clang -Xclang -fopenmp -I openmp/include -L openmp/lib -lomp -Wall -O3 find_orbit.c -o find_orbit
export DYLD_LIBRARY_PATH=/Users/mmuszkow/pyorbit/openmp/lib
for id in 25338 40069 28654 33591; do ./find_orbit n2yo.txt $id >> orbits.txt; done

