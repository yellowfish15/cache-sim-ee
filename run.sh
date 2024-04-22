#!/bin/bash

# path to the Spec_Benchmark folder
SPEC_BENCHMARK_DIR="./Spec_Benchmark/Spec_Benchmark/"

# C++ simulator
CPP_PROGRAM="./sim.cpp"
# Executable file name
EXECUTABLE="main"

g++ -o "$EXECUTABLE" "$CPP_PROGRAM"

# loop through each file in the Spec_Benchmark folder
for file in "$SPEC_BENCHMARK_DIR"*
do
    if [ -f "$file" ]; then
        echo "~===================================~"
        echo "Running $CPP_PROGRAM for file: $file"
        # extract filename from full path
        filename=$(basename -- "$file")
        # run C++ program with the filename
        "./$EXECUTABLE" "$file" $1$ $2$
    fi
done
