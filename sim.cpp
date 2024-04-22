#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

using namespace std;

// convert a hex string to an unsigned 32-bit integer
unsigned int hexStringToUInt(const string& hexString) {
    unsigned int result;
    stringstream ss;
    ss << hex << hexString;
    ss >> result;
    return result;
}

class DRAM {
private:
	// stats
	unsigned int misaligned = 0; // number of misaligned memory accesses
	unsigned int reads = 0; // number of DRAM reads
	unsigned int writes = 0; // number of DRAM writes
	double energy = 0; // total energy consumed
public:
    // read data from DRAM
    double read(unsigned int address, double elapsed) {
		if(address%64 != 0)
			misaligned++;
		energy += elapsed*0.8; // idle
		energy += 50*4; // active
		reads++;
        return 50; // 50 nsec
    }
	
	// add energy spent idle (state = false) or active (state = true)
	// time is in nanoseconds
	void addEnergy(double time, bool state) {
		if(state) // active reads/writes
			energy += time*4;
		else // idle
			energy += time*0.8;
	}
	
	unsigned int getMisaligned() {
		return misaligned;
	}
	
	unsigned int getReads() {
		return reads;
	}
	
	unsigned int getWrites() {
		return writes;
	}
	
	double getEnergy() {
		return energy;
	}
};

class L2 {
private:
    vector<vector<pair<bool, unsigned int>>> sets;
	const static unsigned int size = 256 * 1024; // size of L2 cache in bytes
	const static unsigned int lineSize = 64; // line size in bytes
    const unsigned int setSize;
    const unsigned int associativity;
    const unsigned int numSets;

    // address mapping function to calculate index and tag
    void addressMapping(unsigned int address, unsigned int& index, unsigned int& tag) const {
        index = (address / lineSize) % numSets; // set number
        tag = address / (lineSize * numSets);
    }
	
	// stats
	unsigned int reads = 0; // number of L2 reads
	unsigned int writes = 0; // number of L2 writes
	unsigned int hits = 0;
	unsigned int misses = 0;
	double energy = 0; // total energy consumed
	
public:
    L2(unsigned int assoc) :
        setSize(size / (lineSize * assoc)),
        associativity(assoc),
        numSets(size / lineSize) {
			sets.resize(numSets, vector<pair<bool, unsigned int>>(associativity, make_pair(false, 0)));
		}

    // Read data from cache
    double read(unsigned int address, double elapsed, DRAM& dram) {
		reads++;
        unsigned int index, tag;
        addressMapping(address, index, tag);

		double time = 5;
		energy += time*2;

		bool hit = false;
        // Check each line in the set for a match
        for (const auto& p : sets[index]) {
            if (p.first && p.second == tag) { // cache hit
                hits++;
				hit = true;
				break;
            }
        }
		if(!hit) { // cache miss
			double accessTime = dram.read(address, elapsed+time);
			time += accessTime;
			elapsed += accessTime;
			energy += 0.64; // add 640 picojoules for DRAM access
			misses++;
		}
		energy += elapsed*0.8;
        return time;
    }

    // Write data to cache
    void write(unsigned int address) {
		writes++;
        unsigned int index, tag;
        addressMapping(address, index, tag);

		energy += 5*2; // active energy

        // Update cache
		// Check each line in the set for a match
		bool hit = false;
        for (const auto& p : sets[index]) {
            if (p.first && p.second == tag) {
				// update this line
				hits++;
				hit = true;
				break;
            }
        }
        if(!hit) { // tag not in cache
			misses++;
			bool found = false;
			for (auto& p : sets[index]) {
				if (!p.first) { // replace an empty slot
					// update this line
					p.first = true;
					p.second = tag;
					found = true;
				}
			}
			if(!found) { // random eviction
				pair<bool, unsigned int>& line = sets[index][rand() % associativity];
				line.first = true;
				line.second = tag;
			}
		}
    }
	
	unsigned int getReads() {
		return reads;
	}
	
	unsigned int getWrites() {
		return writes;
	}
	
	unsigned int getHits() {
		return hits;
	}
	
	unsigned int getMisses() {
		return misses;
	}
	
	double getEnergy() {
		return energy;
	}
};

// split cache between instruction and data
class L1 {
private:
    const static unsigned int cacheSize = 32 * 1024; // size of the instruction and data caches in bytes
    const static unsigned int lineSize = 64; // Size of each cache line in bytes
	const static unsigned int numLines = cacheSize/lineSize;
	
	unsigned int indexBits; // Number of index bits for address mapping
    unsigned int offsetBits; // Number of offset bits for address mapping
    unsigned int tagBits; // Number of tag bits for address mapping
	
	// type: 0 = mem, 1 = instr
	bool valid[2][numLines]; // valid bit for each line
	unsigned int cache[2][numLines]; // tag value for each line

	// stats
	unsigned int reads = 0; // number of L1 reads
	unsigned int writes = 0; // number of L1 writes
	unsigned int hits = 0;
	unsigned int misses = 0;
	double energy = 0; // total energy consumed

    // Address mapping function to calculate index and tag
    void addressMapping(unsigned int address, unsigned int& index, unsigned int& tag) const {
        index = (address >> offsetBits) & ((1 << indexBits) - 1);
        tag = address >> (indexBits + offsetBits);
    }

public:

    L1() {
        // Calculate index bits, offset bits, and tag bits based on cache size and line size
        indexBits = log2(cacheSize / lineSize);
        offsetBits = log2(lineSize);
        tagBits = 32 - indexBits - offsetBits;

		// reset valid
		for(int i = 0; i < numLines; i++) {
			valid[0][i] = false;
			valid[1][i] = false;
		}
    }

    // Read data from cache
	// type: 0 = mem, 1 = instr
	// return amount of time elapsed
    double read(unsigned int type, unsigned int address, L2& l2, DRAM& dram) {
        unsigned int index, tag;
        addressMapping(address, index, tag);

		double time = 0.5;
		energy += 0.5; // active

        if (valid[type][index] && cache[type][index] == tag) { // cache hit
			hits++;
        } else { // cache miss
			misses++;
			// read L2
			double elapsed = l2.read(address, time, dram);
			energy += 0.5*elapsed; // idle
			energy += 0.005; // add 5 picojoules for L2 access
			time += elapsed;
			// update L1
			valid[type][index] = 1;
			cache[type][index] = tag;
        }
		reads++;
		return time;
    }

    // Write data to cache
	// assuming memory write
	// type: 0 = mem, 1 = instr
	// return amount of time elapsed
    void write(unsigned int type, unsigned int address) {
        unsigned int index, tag;
		energy += 4.5*0.5+0.5; // idle + active energy
        addressMapping(address, index, tag);
		if (!valid[type][index] || cache[type][index] != tag) {
			misses++;
			// update L1
			cache[type][index] = tag;
			valid[type][index] = true;
		} else {
			hits++;
		}
		writes++;
    }
	
	unsigned int getReads() {
		return reads;
	}
	
	unsigned int getWrites() {
		return writes;
	}
	
	unsigned int getHits() {
		return hits;
	}
	
	unsigned int getMisses() {
		return misses;
	}
	
	double getEnergy() {
		return energy;
	}
};

class System {
private:
    L1 l1;
	L2 l2;
	DRAM dram;

public:
	System(int associativity) : l2(associativity) {}

    void run(string fname) {
        ifstream inputFile(fname); // input file
        if (!inputFile) { // file doesn't exist
            cerr << "Failed to open input file!" << endl;
            exit(-1); // error
        }

		// statistics
		double elapsed = 0; // time elapsed in nanoseconds
		
        int OP;
        string address, value;
		int access = 0;
        while (inputFile >> OP >> address >> value) {
            unsigned int addr = hexStringToUInt(address);
            unsigned int val = hexStringToUInt(value);
			
			// cout << "OP: " << OP << ", Address: " << address << " -> " << addr << ", Value: " << value << " -> " << val << endl;
            switch (OP) { // trace operations
                case 0: { // memory read
					elapsed += l1.read(0, addr, l2, dram);
                    break;
                } case 1: { // memory write
					l1.write(0, addr);
					l2.write(addr);
					elapsed += 5; // 5 ns
                    break;
                } case 2: { // instruction fetch
					elapsed += l1.read(1, addr, l2, dram);
                    break;
                } case 3: // ignore
                    break;
                case 4: // flush cache
                    break;
                default:
                    exit(-1); // error unknown operation
            }
        }

        inputFile.close();
		
		// print out statistics
		cout << "--------\n";
		cout << "Test: " << fname << '\n';
		cout << "L1 Reads: " << l1.getReads() << '\n'; 
		cout << "L1 Writes: " << l1.getWrites() << '\n'; 
		cout << "L1 Cache Hits: " << l1.getHits() << '\n';
		cout << "L1 Cache Misses: " << l1.getMisses() << '\n';
		cout << "L1 Energy (mJ): " << l1.getEnergy()/1000000 << '\n';
		cout << '\n';
		
		cout << "L2 Reads: " << l2.getReads() << '\n'; 
		cout << "L2 Writes: " << l2.getWrites() << '\n'; 
		cout << "L2 Cache Hits: " << l2.getHits() << '\n';
		cout << "L2 Cache Misses: " << l2.getMisses() << '\n';
		cout << "L2 Energy (mJ): " << l2.getEnergy()/1000000 << '\n';
		cout << '\n';
		
		cout << "DRAM Reads: " << dram.getReads() << '\n'; 
		cout << "DRAM Writes: " << dram.getWrites() << '\n'; 
		cout << "DRAM Energy (mJ): " << dram.getEnergy()/1000000 << '\n';
		// cout << "Misaligned Accesses: " << dram.getMisaligned() << '\n';
		cout << '\n';
		
		cout << "Total Time Elapsed (mS): " << elapsed/1000000 << '\n';
    }
};

int main() {
	System sys(2);
    sys.run("./Spec_Benchmark/Spec_Benchmark/047.tomcatv.din");
	return 0;
}