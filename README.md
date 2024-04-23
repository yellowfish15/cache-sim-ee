
  

# Cache Simulator  

**Authors**: Yu Lim (yl38856), Leo Lei (ll36476)

This is a cache simulator written in C++ that simulates a two-level cache hierarchy with a main memory (DRAM). It provides statistics on cache hits, misses, energy consumption, and time elapsed for the Spec_Benchmark workload.

## Prerequisites

To compile and run this code, you need:

- C++ compiler (g++)
- Standard C++ Library
- Bash shell

## Usage

To run the simulator on the workload, open a terminal, navigate to the directory containing the shell script (`run.sh`).

The script takes two input arguments:

1.  **Set Associativity:** This argument determines the set associativity parameter to be passed to the C++ program. It influences the behavior of the simulation algorithm.
2.  **Number of Trials:** This argument specifies the number of trials or iterations to run the simulation for each file.


## Execution

Once the setup is complete, execute the script. It will:

1.  Compile the C++ simulator program.
2.  Iterate through each file in the Spec_Benchmark directory.
3.  For each file:
    -   Display a message indicating the file being processed.
    -   Run the simulator on the file with the specified set associativity and number of trials.


### Example Command

`./run.sh <set_associativity> <num_trials>` 

Replace `<set_associativity>` with the desired L2 set associativity level and `<num_trials>` with the desired number of trials.

## Report and Results

Our [report](https://github.com/yellowfish15/cache-sim-ee/blob/main/report.pdf) in pdf format.
