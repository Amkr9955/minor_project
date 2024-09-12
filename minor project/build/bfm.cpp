#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>

using namespace std;

const int LATTICE_SIZE = 100;
const int CHAIN_LENGTH = 25;  
const double MIN_BOND_LENGTH_SQUARED = 4.0;   
const double MAX_BOND_LENGTH_SQUARED = 16.0;  


const double allowed_bond_lengths_squared[] = {4.0, 5.0, 8.0, 9.0, 10.0, 13.0};

// Fixed X positions for the two chains
const int X_FIXED_1 = 0;  // X position for first chain
const int X_FIXED_2 = 2;  // X position for second chain
const int Z_MIN = 2;      // Z position must be >= 2

// Structure to hold coordinates of a monomer in 3D
struct Monomer {
    int x, y, z;
};

void initialize_chains_from_file(vector<Monomer>& chain1, vector<Monomer>& chain2, const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    string line;

    // Read coordinates for chain 1
    for (int i = 0; i < CHAIN_LENGTH; ) {
        getline(infile, line);  // Read each line
        if (line.empty() || line[0] == '#') {
            continue; // Skip comment lines or empty lines
        }
        stringstream ss(line);
        ss >> chain1[i].x >> chain1[i].y >> chain1[i].z;
        i++;
    }

    // Read coordinates for chain 2
    for (int i = 0; i < CHAIN_LENGTH; ) {
        getline(infile, line);  // Read each line
        if (line.empty() || line[0] == '#') {
            continue; // Skip comment lines or empty lines
        }
        stringstream ss(line);
        ss >> chain2[i].x >> chain2[i].y >> chain2[i].z;
        i++;
    }

    infile.close();
}

// Calculate squared Euclidean distance between two monomers in 3D
double distance_squared(const Monomer& m1, const Monomer& m2) {
    return pow(m1.x - m2.x, 2) + pow(m1.y - m2.y, 2) + pow(m1.z - m2.z, 2);
}

// Check if the bond length between two monomers is allowed
bool is_valid_bond_length(double dist_sq) {
    for (double allowed_dist_sq : allowed_bond_lengths_squared) {
        if (fabs(dist_sq - allowed_dist_sq) < 1e-6) {  // Check if distance is one of the allowed values
            return true;
        }
    }
    return false;
}

// Check if a proposed move is valid (self-avoidance and bond length constraints)
bool is_valid_move(const vector<Monomer>& chain1, const vector<Monomer>& chain2, int chain_index, int monomer_index, const Monomer& new_pos) {
    const vector<Monomer>& current_chain = (chain_index == 0) ? chain1 : chain2;
    const vector<Monomer>& other_chain = (chain_index == 0) ? chain2 : chain1;

    if (new_pos.z < 2) {
        return false;  // Invalid move, z cannot be less than 2
    }

    // Check bond length constraints within the current chain
    if (monomer_index > 0 && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index - 1]))) return false; 
    if (monomer_index < CHAIN_LENGTH - 1 && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index + 1]))) return false;

    // Check self-avoidance within the current chain
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        if (i != monomer_index && distance_squared(new_pos, current_chain[i]) < MIN_BOND_LENGTH_SQUARED) {
            return false; // Overlap detected
        }
    }

    // Check distance constraints with the other chain
    if (distance_squared(new_pos, other_chain[monomer_index]) > MAX_BOND_LENGTH_SQUARED) {  // Max distance allowed is 4
            return false; // bond of 2 monomers of same index in different chain breaks 
        }
 

    return true;
}

// Randomly move a monomer in the YZ plane (X is fixed)
void monte_carlo_move(vector<Monomer>& chain1, vector<Monomer>& chain2) {
    int chain_index = rand() % 2;  // Select a random chain (0 for chain1, 1 for chain2)
    int monomer_index = rand() % CHAIN_LENGTH;  // Select a random monomer in the chosen chain

    Monomer new_pos;
    if (chain_index == 0) {
        new_pos = chain1[monomer_index];
    } else {
        new_pos = chain2[monomer_index];
    }

    int direction = rand() % 4;  // Random direction (0: +y, 1: -y, 2: +z, 3: -z)

    // Update Y or Z position based on the direction
    switch (direction) {
        case 0: new_pos.y += 1; break;  // Move +y
        case 1: new_pos.y -=1; break;  // Move -y,
        case 2: new_pos.z += 1; break;  // Move +z
        case 3: new_pos.z -=1; break;  // Move -z,
    }

    // Check if the move is valid, otherwise revert to original position
    if (is_valid_move(chain1, chain2, chain_index, monomer_index, new_pos)) {
        if (chain_index == 0) {
            chain1[monomer_index] = new_pos;
        } else {
            chain2[monomer_index] = new_pos;
        }
    }
}

// Save chains' conformation to a file
void save_chains(const vector<Monomer>& chain1, const vector<Monomer>& chain2, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    outfile << "Chain 1 (X = 0):" << endl;
    for (const auto& monomer : chain1) {
        outfile << monomer.x << " " << monomer.y << " " << monomer.z << endl;
    }

    outfile << "\nChain 2 (X = 2):" << endl;
    for (const auto& monomer : chain2) {
        outfile << monomer.x << " " << monomer.y << " " << monomer.z << endl;
    }

    outfile.close();
}

void save_chains_to_csv(const vector<Monomer>& chain1, const vector<Monomer>& chain2, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    // Write the CSV header
    outfile << "Chain,Monomer,X,Y,Z" << endl;

    // Save chain 1 coordinates with correct formatting
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        outfile << "Chain1," << i << "," << chain1[i].x << "," << chain1[i].y << "," << chain1[i].z << endl;
    }

    // Save chain 2 coordinates with correct formatting
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        outfile << "Chain2," << i << "," << chain2[i].x << "," << chain2[i].y << "," << chain2[i].z << endl;
    }

    outfile.close();
}

// Main simulation loop
int main() {
    srand(time(0));

    vector<Monomer> chain1(CHAIN_LENGTH);  // Chain 1
    vector<Monomer> chain2(CHAIN_LENGTH);  // Chain 2

    // Initialize chains by reading from the input file
    string input_file = "polymer_input.txt";
    initialize_chains_from_file(chain1, chain2, input_file);

    cout << "Initial Chains' Conformation loaded from " << input_file << endl;
    save_chains(chain1, chain2, "initial_chains.txt");  // Save initial conformation
    save_chains_to_csv(chain1, chain2, "initial_chains.csv");
    

    // Perform Monte Carlo simulation
    int monte_carlo_steps = 1000000;
    for (int step = 0; step < monte_carlo_steps; step++) {
        monte_carlo_move(chain1, chain2);
    }

    // Save final conformation to output file
    save_chains(chain1, chain2, "final_chains.txt");
    save_chains_to_csv(chain1, chain2, "final_chains.csv");

    cout << "Final Chains' Conformation saved to final_chains.txt" << endl;

    return 0;
}
