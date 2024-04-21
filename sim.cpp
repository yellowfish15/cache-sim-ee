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
    double read(unsigned int type, unsigned int address, bool& hit) {
        unsigned int index, tag;
        addressMapping(address, index, tag);

        if (valid[type][index] && cache[type][index] == tag) { // cache hit
			hits++;
            hit = true;
        } else { // cache miss
			misses++;
            hit = false;
        }
		return 0.5;
    }

    // Write data to cache
	// assuming memory write
	// return amount of time elapsed
    double write(unsigned int address, unsigned int data) {
        unsigned int index, tag;
        addressMapping(address, index, tag);
        cache[0][index] = tag;
		valid[0][index] = true;
		return 0.5;
    }
	
	// add energy spent idle (state = false) or active (state = true)
	// time is in nanoseconds
	void addEnergy(double time, bool state) {
		if(state) // active reads/writes
			energy += time;
		else
			energy += time*0.5;
	}
	
	int getHits() {
		return hits;
	}
	
	int getMisses() {
		return misses;
	}
	
	double getEnergy() {
		return energy;
	}
};

class System {
private:
    L1 l1;

public:
	System() {}

    void run(string fname) {
        ifstream inputFile(fname); // input file
        if (!inputFile) { // file doesn't exist
            cerr << "Failed to open input file!" << endl;
            exit(-1); // error
        }

        int OP;
        string address, value;
		double elapsed = 0; // time elapsed in nanoseconds
        while (inputFile >> OP >> address >> value) {
            unsigned int addr = hexStringToUInt(address);
            unsigned int val = hexStringToUInt(value);
			// cout << "OP: " << OP << ", Address: " << address << " -> " << addr << ", Value: " << value << " -> " << val << endl;
            switch (OP) { // trace operations
                case 0: { // memory read
					bool hit = false;
					double time = l1.read(0, addr, hit);
					elapsed += time;
					l1.addEnergy(time, true);
					if(!hit) { // missed, check L2
						
					}
                    break;
                } case 1: { // memory write
                    double time = l1.write(addr, val);
					elapsed += time;
					l1.addEnergy(time, true);
                    break;
                } case 2: { // instruction fetch
					bool hit = false;
					double time = l1.read(1, addr, hit);
					elapsed += time;
					l1.addEnergy(time, true);
					if(!hit) { // missed, check L2
						
					}
                    break;
                } case 3: // ignore
                    break;
                case 4: // flush cache
                    // Not implemented yet
                    break;
                default:
                    exit(-1); // error unknown operation
            }
        }

        inputFile.close();
		
		// print out statistics
		cout << "Test: " << fname << '\n';
		cout << "L1 Cache Hits: " << l1.getHits() << '\n';
		cout << "L1 Cache Misses: " << l1.getMisses() << '\n';
		cout << "L1 Energy (nJ): " << l1.getEnergy() << '\n';
    }
};

int main() {
	System sys;
    sys.run("023.eqntott.din");
	return 0;
}