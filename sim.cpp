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
	
	// type: 0 = mem, 1 = instr
	bool valid[2][numLines]; // valid bit for each line
	unsigned int cache[2][numLines]; // tag value for each line
	
	
    unsigned int indexBits; // Number of index bits for address mapping
    unsigned int offsetBits; // Number of offset bits for address mapping
    unsigned int tagBits; // Number of tag bits for address mapping

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
    bool read(unsigned int type, unsigned int address) {
        unsigned int index, tag;
        addressMapping(address, index, tag);

        if (valid[type][index] && cache[type][index] == tag) { // cache hit
            return true;
        } else { // cache miss
            return false;
        }
    }

    // Write data to cache
	// assuming memory write
    void write(unsigned int address, unsigned int data) {
        unsigned int index, tag;
        addressMapping(address, index, tag);

        cache[0][index] = tag;
		valid[0][index] = true;
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
        while (inputFile >> OP >> address >> value) {
            unsigned int addr = hexStringToUInt(address);
            unsigned int val = hexStringToUInt(value);
			// cout << "OP: " << OP << ", Address: " << address << " -> " << addr << ", Value: " << value << " -> " << val << endl;
            switch (OP) { // trace operations
                case 0: // memory read
                    l1.read(0, addr);
                    break;
                case 1: // memory write
                    l1.write(addr, val);
                    break;
                case 2: // instruction fetch
					l1.read(1, addr);
                    break;
                case 3: // ignore
                    break;
                case 4: // flush cache
                    // Not implemented yet
                    break;
                default:
                    exit(-1); // error unknown operation
            }
        }

        inputFile.close();
    }
};

int main() {
	System sys;
    sys.run("023.eqntott.din");
	return 0;
}