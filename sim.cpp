#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
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

// trace operations
enum MemoryOperation {
    MEMORY_READ = 0,
    MEMORY_WRITE = 1,
    INSTRUCTION_FETCH = 2,
    IGNORE = 3,
    FLUSH_CACHE = 4
};

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
    }

    inputFile.close();
}

int main() {
	run("023.eqntott.din");
	return 0;
}