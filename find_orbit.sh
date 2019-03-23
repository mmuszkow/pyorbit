#!/bin/bash

# Tested on OSX High Sierra

#gcc -Wall -O3 find_orbit.c -march=native -o find_orbit && \
while true; do
    for id in 17181 20253; do ./find_orbit dumped.txt $id 1500 tle_active.txt >> ${id}.txt; done
done

