#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>

using namespace std;

// system variables
static double sysclock = 0;

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
	double clock = 0; // internal clock

	// stats
	unsigned int reads = 0; // number of DRAM reads
	unsigned int writes = 0; // number of DRAM writes
	double energy = 0; // total energy consumed
public:
	// sync energy with sysclock (calculate idle energy)
	void sync() {
		if(sysclock > clock) {
			energy += (sysclock-clock)*0.8;
			clock = sysclock;
		}
	}

    // read data from DRAM
    void read(unsigned int address) {
		sync(); 
		clock += 50;
		sysclock = clock;
		energy += 50*4; // active energy
		reads++;
    }

	// write data to DRAM
    void write(unsigned int address) {
		sync();
		/*
		// we don't update clock because DRAM write is asynchronous
		clock += 45;
		sysclock = clock;
		*/
		energy += 45*4; // active energy
		writes++;
    }
	
	unsigned int getReads() {
		return reads;
	}
	
	unsigned int getWrites() {
		return writes;
	}
	
	double getEnergy() {
		sync();
		return energy;
	}
};

class L2 {
private:
	double clock = 0; // internal clock

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

	// sync energy with sysclock (calculate idle energy)
	void sync() {
		if(sysclock > clock) {
			energy += (sysclock-clock)*0.8;
			clock = sysclock;
		}
	}

    // Read data from cache
    void read(unsigned int address, DRAM& dram) {
		reads++;
        unsigned int index, tag;
        addressMapping(address, index, tag);

		sync();
		clock += 5;
		sysclock = clock;
		energy += 5*2; // active energy

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
			dram.read(address);
			energy += 0.64; // add 640 picojoules for DRAM access
			misses++;

			// update L2
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
				// we need to evict this cache line to DRAM
				dram.write(address);
			}
		}
    }

    // Write data to cache
    void write(unsigned int address, DRAM& dram) {
		writes++;
        unsigned int index, tag;
        addressMapping(address, index, tag);

		sync();
		clock += 4.5;
		sysclock = clock;
		energy += 4.5*2; // active energy

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

				// we need to evict this cache line to DRAM
				dram.write(address);
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

	double getHitRate() {
		return ((double)hits)/(hits+misses);
	}
	
	double getEnergy() {
		sync();
		return energy;
	}
};

// split cache between instruction and data
class L1 {
private:
	double clock = 0; // internal clock

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
	unsigned int hits[2];
	unsigned int misses[2];
	double energy[2]; // total energy consumed

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

		hits[0] = hits[1] = 0;
		misses[0] = misses[1] = 0;
		energy[0] = energy[1] = 0;

		// reset valid
		for(int i = 0; i < numLines; i++) {
			valid[0][i] = false;
			valid[1][i] = false;
		}
    }

	// sync energy with sysclock (calculate idle energy)
	void sync() {
		if(sysclock > clock) {
			energy[0] += (sysclock-clock)*0.5;
			energy[1] += (sysclock-clock)*0.5;
			clock = sysclock;
		}
	}

    // Read data from cache
	// type: 0 = mem, 1 = instr
	// return amount of time elapsed
    void read(unsigned int type, unsigned int address, L2& l2, DRAM& dram) {
        unsigned int index, tag;
        addressMapping(address, index, tag);

		sync();
		clock += 0.5;
		sysclock = clock;
		energy[type] += 0.5*1; // active energy
		energy[1-type] += 0.5*0.5; // idle energy

        if (valid[type][index] && cache[type][index] == tag) { // cache hit
			hits[type]++;
        } else { // cache miss
			misses[type]++;
			// read L2
			l2.read(address, dram);
			energy[type] += 0.005; // add 5 picojoules for L2 access
			// update L1
			valid[type][index] = 1;
			cache[type][index] = tag;
        }
		reads++;
    }

    // Write data to cache
	// assuming memory write
	// type: 0 = mem, 1 = instr
	// return amount of time elapsed
    void write(unsigned int type, unsigned int address) {
        unsigned int index, tag;

		sync();
		clock += 0.5;
		sysclock = clock;
		energy[type] += 0.5*1; // active energy
		energy[1-type] += 0.5*0.5; // idle energy

        addressMapping(address, index, tag);
		if (!valid[type][index] || cache[type][index] != tag) {
			misses[type]++;
			// update L1
			cache[type][index] = tag;
			valid[type][index] = true;
		} else {
			hits[type]++;
		}
		writes++;
    }
	
	unsigned int getReads() {
		return reads;
	}
	
	unsigned int getWrites() {
		return writes;
	}
	
	// type: 0 = mem, 1 = instr
	unsigned int getHits(int type) {
		return hits[type];
	}
	
	// type: 0 = mem, 1 = instr
	unsigned int getMisses(int type) {
		return misses[type];
	}

	// type: 0 = mem, 1 = instr
	double getHitRate(int type) {
		return ((double)hits[type])/(hits[type]+misses[type]);
	}
	
	// type: 0 = mem, 1 = instr
	double getEnergy(int type) {
		sync();
		return energy[type];
	}
};

// contains all the data across all runs
struct Result {
	unsigned int misaligned = 0; // number of misaligned addresses
	unsigned int l1d_hits = 0;
	unsigned int l1d_miss = 0;
	double l1d_hit_rate = 0;
	double l1d_energy = 0; // energy consumption in mJ

	unsigned int l1i_hits = 0;
	unsigned int l1i_miss = 0;
	double l1i_hit_rate = 0;
	double l1i_energy = 0; // energy consumption in mJ

	unsigned int l2_hits = 0;
	unsigned int l2_miss = 0;
	double l2_hit_rate = 0;
	double l2_energy = 0; // energy consumption in mJ

	unsigned int dram_reads = 0;
	unsigned int dram_writes = 0;
	double dram_energy = 0; // energy consumption in mJ
	double time_elapsed = 0; // time elapsed in mS 
};

// contains all the data across all runs
struct Stats {
	unsigned int misaligned = 0; // number of misaligned addresses
	double l1d_hits = 0;
	double l1d_miss = 0;
	double l1d_hit_rate = 0;
	double l1d_energy = 0; // energy consumption in mJ

	double l1i_hits = 0;
	double l1i_miss = 0;
	double l1i_hit_rate = 0;
	double l1i_energy = 0; // energy consumption in mJ

	double l2_hits = 0;
	double l2_miss = 0;
	double l2_hit_rate = 0;
	double l2_energy = 0; // energy consumption in mJ

	double dram_reads = 0;
	double dram_writes = 0;
	double dram_energy = 0; // energy consumption in mJ
	double time_elapsed = 0; // time elapsed in mS 
};

void run(string fname, int associativity, struct Result& results) {
	// configure system
	L1 l1;
	L2 l2(associativity);
	DRAM dram;
	sysclock = 0; // reset the system clock
	int misalign_count = 0;

	ifstream inputFile(fname); // input file
	if (!inputFile) { // file doesn't exist
		cerr << "Failed to open input file!" << endl;
		exit(-1); // error
	}
	
	int OP;
	string address, value;
	int access = 0;
	while (inputFile >> OP >> address >> value) {
		unsigned int addr = hexStringToUInt(address);
		unsigned int val = hexStringToUInt(value);

		// misaligned address
		bool misaligned = false;
		if(addr%64 > 60) {
			misalign_count++;
			misaligned = true;
		}

		// cout << "OP: " << OP << ", Address: " << address << " -> " << addr << ", Value: " << value << " -> " << val << '\n';
		switch (OP) { // trace operations
			case 0: { // memory read
				l1.read(0, addr, l2, dram);
				if(misaligned)
					l1.read(0, addr, l2, dram);
				break;
			} case 1: { // memory write
				l1.write(0, addr);
				l2.write(addr, dram);
				if(misaligned) {
					l1.write(0, addr);
					l2.write(addr, dram);
				}
				break;
			} case 2: { // instruction fetch
				l1.read(1, addr, l2, dram);
				if(misaligned)
					l1.read(1, addr, l2, dram);
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
	
	// get results
	results.misaligned = misalign_count; // get number of misaligned addresses
	results.l1d_hits = l1.getHits(0);
	results.l1d_miss = l1.getMisses(0);
	results.l1d_hit_rate = l1.getHitRate(0);
	results.l1d_energy = l1.getEnergy(0)/1000000; // convert from nJ to mJ
	results.l1i_hits = l1.getHits(1);
	results.l1i_miss = l1.getMisses(1);
	results.l1i_hit_rate = l1.getHitRate(1);
	results.l1i_energy = l1.getEnergy(1)/1000000; // convert from nJ to mJ
	results.l2_hits = l2.getHits();
	results.l2_miss = l2.getMisses();
	results.l2_energy = l2.getEnergy()/1000000; // convert from nJ to mJ
	results.l2_hit_rate = l2.getHitRate();
	results.dram_reads = dram.getReads();
	results.dram_writes = dram.getWrites();
	results.dram_energy = dram.getEnergy()/1000000; // convert from nJ to mJ
	results.time_elapsed = sysclock/1000000; // convert from nS to mS
}

// for debugging
int main() {
	string fname = "./Spec_Benchmark/Spec_Benchmark/047.tomcatv.din";
	int associativity = 4;
	Result r;
	run(fname, associativity, r);

	const int w1 = 20; // width of fields
	const int w2 = 8; // width of value

	// print out statistics
	cout << "Test Results: " << fname << '\n';
	cout << "L2 Associativity Level: " << associativity << '\n';
	cout << '\n';

	cout << setw(w1) << "L1i Cache Hits: " << setw(w2) << r.l1i_hits << '\n';
	cout << setw(w1) << "L1i Cache Misses: " << setw(w2) << r.l1i_miss << '\n';
	cout << setw(w1) << "L1i Hit Rate: " << setw(w2) << r.l1i_hit_rate << '\n';
	cout << setw(w1) << "L1i Energy (mJ): " << setw(w2) << r.l1i_energy << '\n';
	cout << '\n';

	cout << setw(w1) << "L1d Cache Hits: " << setw(w2) << r.l1d_hits << '\n';
	cout << setw(w1) << "L1d Cache Misses: " << setw(w2) << r.l1d_miss << '\n';
	cout << setw(w1) << "L1d Hit Rate: " << setw(w2) << r.l1d_hit_rate << '\n';
	cout << setw(w1) << "L1d Energy (mJ): " << setw(w2) << r.l1d_energy << '\n';
	cout << '\n';
	
	cout << setw(w1) << "L2 Cache Hits: " << setw(w2) << r.l2_hits << '\n';
	cout << setw(w1) << "L2 Cache Misses: " << setw(w2) << r.l2_miss << '\n';
	cout << setw(w1) << "L2 Hit Rate: " << setw(w2) << r.l2_hit_rate << '\n';
	cout << setw(w1) << "L2 Energy (mJ): " << setw(w2) << r.l2_energy << '\n';
	cout << '\n';
	
	cout << setw(w1) << "DRAM Reads: " << setw(w2) << r.dram_reads << '\n';
	cout << setw(w1) << "DRAM Writes: " << setw(w2) << r.dram_writes << '\n';
	cout << setw(w1) << "DRAM Energy (mJ): " << setw(w2) << r.dram_energy << '\n';
	cout << '\n';

	cout << setw(w1) << "Misaligned Addr: " << setw(w2) << r.misaligned << '\n';
	cout << '\n';
	
	cout << setw(w1) << "Total Energy Consumption (mJ): " << (r.l1d_energy+r.l1i_energy+r.l2_energy+r.dram_energy) << '\n';
	cout << setw(w1) << "Total Time Elapsed (mS): " << r.time_elapsed << '\n';
}


/*
int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <filename> <value>" << endl;
        return 1;
    }
    
    string fname = argv[1];
    int associativity = stoi(argv[2]);
	int trials = stoi(argv[3]); // number of trials to run
	if(trials < 1) {
		cerr << "The number of trials to run cannot be less than 1 (input was " << trials << ")!" << endl;
        return 1;
	}
	if(associativity != 2 && associativity != 4 && associativity != 8) {
		cerr << "L2 Set Associativitiy must be 2, 4, or 8 (input was " << associativity << ")!" << endl;
        return 1;
	}

	Result results[trials]; // list of all results across all trials
	for(int i = 0; i < trials; i++) {
    	run(fname, associativity, results[i]);
	}

	if(trials > 1) { // calculate distribution statistics
		Stats mean; // Calculate mean
		for (const auto& result : results) {
			mean.misaligned = result.misaligned;
			mean.l1d_hits += result.l1d_hits;
			mean.l1d_miss += result.l1d_miss;
			mean.l1d_hit_rate += result.l1d_hit_rate;
			mean.l1d_energy += result.l1d_energy;
			mean.l1i_hits += result.l1i_hits;
			mean.l1i_miss += result.l1i_miss;
			mean.l1i_hit_rate += result.l1i_hit_rate;
			mean.l1i_energy += result.l1i_energy;
			mean.l2_hits += result.l2_hits;
			mean.l2_miss += result.l2_miss;
			mean.l2_hit_rate += result.l2_hit_rate;
			mean.l2_energy += result.l2_energy;
			mean.dram_reads += result.dram_reads;
			mean.dram_writes += result.dram_writes;
			mean.dram_energy += result.dram_energy;
			mean.time_elapsed += result.time_elapsed;
		}
		mean.l1d_hits /= trials;
		mean.l1d_miss /= trials;
		mean.l1d_hit_rate /= trials;
		mean.l1d_energy /= trials;
		mean.l1i_hits /= trials;
		mean.l1i_miss /= trials;
		mean.l1i_hit_rate /= trials;
		mean.l1i_energy /= trials;
		mean.l2_hits /= trials;
		mean.l2_miss /= trials;
		mean.l2_hit_rate /= trials;
		mean.l2_energy /= trials;
		mean.dram_reads /= trials;
		mean.dram_writes /= trials;
		mean.dram_energy /= trials;
		mean.time_elapsed /= trials;

		// Calculate standard deviation
		Stats stddev;
		for (const auto& result : results) {
			stddev.l1d_hits += pow(result.l1d_hits - mean.l1d_hits, 2);
			stddev.l1d_miss += pow(result.l1d_miss - mean.l1d_miss, 2);
			stddev.l1d_hit_rate += pow(result.l1d_hit_rate - mean.l1d_hit_rate, 2);
			stddev.l1d_energy += pow(result.l1d_energy - mean.l1d_energy, 2);
			stddev.l1i_hits += pow(result.l1i_hits - mean.l1i_hits, 2);
			stddev.l1i_miss += pow(result.l1i_miss - mean.l1i_miss, 2);
			stddev.l1i_hit_rate += pow(result.l1i_hit_rate - mean.l1i_hit_rate, 2);
			stddev.l1i_energy += pow(result.l1i_energy - mean.l1i_energy, 2);
			stddev.l2_hits += pow(result.l2_hits - mean.l2_hits, 2);
			stddev.l2_miss += pow(result.l2_miss - mean.l2_miss, 2);
			stddev.l2_hit_rate += pow(result.l2_hit_rate - mean.l2_hit_rate, 2);
			stddev.l2_energy += pow(result.l2_energy - mean.l2_energy, 2);
			stddev.dram_reads += pow(result.dram_reads - mean.dram_reads, 2);
			stddev.dram_writes += pow(result.dram_writes - mean.dram_writes, 2);
			stddev.dram_energy += pow(result.dram_energy - mean.dram_energy, 2);
			stddev.time_elapsed += pow(result.time_elapsed - mean.time_elapsed, 2);
		}
		stddev.l1d_hits = sqrt(stddev.l1d_hits / trials);
		stddev.l1d_miss = sqrt(stddev.l1d_miss / trials);
		stddev.l1d_hit_rate = sqrt(stddev.l1d_hit_rate / trials);
		stddev.l1d_energy = sqrt(stddev.l1d_energy / trials);
		stddev.l1i_hits = sqrt(stddev.l1i_hits / trials);
		stddev.l1i_miss = sqrt(stddev.l1i_miss / trials);
		stddev.l1i_hit_rate = sqrt(stddev.l1i_hit_rate / trials);
		stddev.l1i_energy = sqrt(stddev.l1i_energy / trials);
		stddev.l2_hits = sqrt(stddev.l2_hits / trials);
		stddev.l2_miss = sqrt(stddev.l2_miss / trials);
		stddev.l2_hit_rate = sqrt(stddev.l2_hit_rate / trials);
		stddev.l2_energy = sqrt(stddev.l2_energy / trials);
		stddev.dram_reads = sqrt(stddev.dram_reads / trials);
		stddev.dram_writes = sqrt(stddev.dram_writes / trials);
		stddev.dram_energy = sqrt(stddev.dram_energy / trials);
		stddev.time_elapsed = sqrt(stddev.time_elapsed / trials);

		// print out statistics
		cout << "Test Results: " << '\n';
		cout << "L2 Associativity Level: " << associativity << '\n';
		cout << '\n';

		const int w1 = 18; // width of fields
		const int w2 = 14; // width of mean
		const int w3 = 12; // width of standard deviation
		const int prec = 3;

		cout << setw(w1) << "Metric:" << setw(w2) << "Mean:" << setw(w3) << "StdDev:" << '\n';
		// Print statistics
		cout << fixed << '\n';
		cout << setw(w1) << "L1d Cache Hits:" << setw(w2) << setprecision(prec) << mean.l1d_hits << setw(w3) << setprecision(prec) << stddev.l1d_hits << '\n';
		cout << setw(w1) << "L1d Cache Misses:" << setw(w2) << setprecision(prec) << mean.l1d_miss << setw(w3) << setprecision(prec) << stddev.l1d_miss << '\n';
		cout << setw(w1) << "L1d Hit Rate:" << setw(w2) << setprecision(prec) << mean.l1d_hit_rate << setw(w3) << setprecision(prec) << stddev.l1d_hit_rate << '\n';
		cout << setw(w1) << "L1d Energy (mJ):" << setw(w2) << setprecision(prec) << mean.l1d_energy << setw(w3) << setprecision(prec) << stddev.l1d_energy << '\n';
		cout << '\n';

		cout << setw(w1) << "L1i Cache Hits:" << setw(w2) << setprecision(prec) << mean.l1i_hits << setw(w3) << setprecision(prec) << stddev.l1i_hits << '\n';
		cout << setw(w1) << "L1i Cache Misses:" << setw(w2) << setprecision(prec) << mean.l1i_miss << setw(w3) << setprecision(prec) << stddev.l1i_miss << '\n';
		cout << setw(w1) << "L1i Hit Rate:" << setw(w2) << setprecision(prec) << mean.l1i_hit_rate << setw(w3) << setprecision(prec) << stddev.l1i_hit_rate << '\n';
		cout << setw(w1) << "L1i Energy (mJ):" << setw(w2) << setprecision(prec) << mean.l1i_energy << setw(w3) << setprecision(prec) << stddev.l1i_energy << '\n';
		cout << '\n';

		cout << setw(w1) << "L2 Cache Hits:" << setw(w2) << setprecision(prec) << mean.l2_hits << setw(w3) << setprecision(prec) << stddev.l2_hits << '\n';
		cout << setw(w1) << "L2 Cache Misses:" << setw(w2) << setprecision(prec) << mean.l2_miss << setw(w3) << setprecision(prec) << stddev.l2_miss << '\n';
		cout << setw(w1) << "L2 Hit Rate:" << setw(w2) << setprecision(prec) << mean.l2_hit_rate << setw(w3) << setprecision(prec) << stddev.l2_hit_rate << '\n';
		cout << setw(w1) << "L2 Energy (mJ):" << setw(w2) << setprecision(prec) << mean.l2_energy << setw(w3) << setprecision(prec) << stddev.l2_energy << '\n';
		cout << '\n';

		cout << setw(w1) << "DRAM Reads:" << setw(w2) << setprecision(prec) << mean.dram_reads << setw(w3) << setprecision(prec) << stddev.dram_reads << '\n';
		cout << setw(w1) << "DRAM Writes:" << setw(w2) << setprecision(prec) << mean.dram_writes << setw(w3) << setprecision(prec) << stddev.dram_writes << '\n';
		cout << setw(w1) << "DRAM Energy (mJ):" << setw(w2) << setprecision(prec) << mean.dram_energy << setw(w3) << setprecision(prec) << stddev.dram_energy << '\n';
		cout << '\n';

		cout << "Misaligned Addr: " << mean.misaligned << '\n';
		cout << '\n';

		cout << "Mean Total Energy Consumption (mJ): " << (mean.l1d_energy+mean.l1i_energy+mean.l2_energy+mean.dram_energy) << '\n';
		cout << "Total Time Elapsed (mS): " << mean.time_elapsed << '\n';
	} else if (trials == 1) {
		const int w1 = 20; // width of fields
		const int w2 = 8; // width of value

		// print out statistics
		cout << "Test Results (" << trials << " trial): " << fname << '\n';
		cout << "L2 Associativity Level: " << associativity << '\n';
		cout << '\n';

		cout << setw(w1) << "L1d Cache Hits: " << setw(w2) << results[0].l1d_hits << '\n';
		cout << setw(w1) << "L1d Cache Misses: " << setw(w2) << results[0].l1d_miss << '\n';
		cout << setw(w1) << "L1d Hit Rate: " << setw(w2) << results[0].l1d_hit_rate << '\n';
		cout << setw(w1) << "L1d Energy (mJ): " << setw(w2) << results[0].l1d_energy << '\n';
		cout << '\n';

		cout << setw(w1) << "L1i Cache Hits: " << setw(w2) << results[0].l1i_hits << '\n';
		cout << setw(w1) << "L1i Cache Misses: " << setw(w2) << results[0].l1i_miss << '\n';
		cout << setw(w1) << "L1i Hit Rate: " << setw(w2) << results[0].l1i_hit_rate << '\n';
		cout << setw(w1) << "L1i Energy (mJ): " << setw(w2) << results[0].l1i_energy << '\n';
		cout << '\n';
		
		cout << setw(w1) << "L2 Cache Hits: " << setw(w2) << results[0].l2_hits << '\n';
		cout << setw(w1) << "L2 Cache Misses: " << setw(w2) << results[0].l2_miss << '\n';
		cout << setw(w1) << "L2 Hit Rate: " << setw(w2) << results[0].l2_hit_rate << '\n';
		cout << setw(w1) << "L2 Energy (mJ): " << setw(w2) << results[0].l2_energy << '\n';
		cout << '\n';
		
		cout << setw(w1) << "DRAM Reads: " << setw(w2) << results[0].dram_reads << '\n';
		cout << setw(w1) << "DRAM Writes: " << setw(w2) << results[0].dram_writes << '\n';
		cout << setw(w1) << "DRAM Energy (mJ): " << setw(w2) << results[0].dram_energy << '\n';
		cout << '\n';

		cout << setw(w1) << "Misaligned Addr: " << setw(w2) << results[0].misaligned << '\n';
		cout << '\n';
		
		cout << setw(w1) << "Total Energy Consumption (mJ): " << (results[0].l1d_energy+results[0].l1i_energy+results[0].l2_energy+results[0].dram_energy) << '\n';
		cout << setw(w1) << "Total Time Elapsed (mS): " << results[0].time_elapsed << '\n';
	}
    return 0;
}
*/