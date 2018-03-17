#!/bin/bash

# Tested on OSX High Sierra

gcc -Wall -O3 find_orbit.c -march=native -o find_orbit && \
for id in 25338 28654 40069 33591; do ./find_orbit n2yo.txt $id 5000 orbits_all_computed.txt >> orbits.txt; done

