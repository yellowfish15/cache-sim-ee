#!/bin/bash

# path to the Spec_Benchmark folder
SPEC_BENCHMARK_DIR="./Spec_Benchmark/Spec_Benchmark/"

# C++ simulator
CPP_PROGRAM="./sim.cpp"

# loop through each file in the Spec_Benchmark folder
for file in "$SPEC_BENCHMARK_DIR"*
do
    if [ -f "$file" ]; then
        echo "Running $CPP_PROGRAM for file: $file"
        # extract filename from full path
        filename=$(basename -- "$file")
        # run C++ program with the filename
        $CPP_PROGRAM "$file" 2
    fi
done
