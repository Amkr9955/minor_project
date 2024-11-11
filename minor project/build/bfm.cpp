#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>


using namespace std;

// Constants
const int LATTICE_SIZE = 100;
const int CHAIN_LENGTH = 40;  
const double MIN_BOND_LENGTH_SQUARED = 4.0;   
const double MAX_BOND_LENGTH_SQUARED = 16.0;  
const double WALL_DISTANCE = 2.0;
const double WALL_DISTANCE_SQUARED = WALL_DISTANCE * WALL_DISTANCE; 

// Temperature range
const double T_START = 0.1;
const double T_END = 3.0;
const double T_STEP = 0.1;
const int k = 1; // Boltzmann constant (set to 1 for reduced units)

// Allowed bond lengths for 4 site lattice model
const double allowed_bond_lengths_squared[] = {4.0, 5.0, 8.0, 9.0, 10.0, 13.0};

// Chains fixed values
const int X_FIXED_1 = 0;  // X position for first chain
const int X_FIXED_2 = 2;  // X position for second chain
const int Z_MIN = 2;      // Z position must be >=2

// Data type of chains
struct Monomer {
    int x, y, z;
};

// Function to initialize chains from polymer input text file
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
        if (line.empty() || line[0] == '#') 
            continue;  // Skip comment lines or empty lines

        stringstream ss(line);
        ss >> chain1[i].x >> chain1[i].y >> chain1[i].z;
        i++;
    }
    // Read coordinates for chain 2
    for (int i = 0; i < CHAIN_LENGTH; ) {
        getline(infile, line);  // Read each line
        if (line.empty() || line[0] == '#') 
            continue;   // Skip comment lines or empty lines

        stringstream ss(line);
        ss >> chain2[i].x >> chain2[i].y >> chain2[i].z;
        i++;
    }

    infile.close();
}

// Function to calculate squared distance between two monomers in 3D
double distance_squared(const Monomer& m1, const Monomer& m2) {
    return pow(m1.x - m2.x, 2) + pow(m1.y - m2.y, 2) + pow(m1.z - m2.z, 2);
}

// Function to check if a bond length is allowed
bool is_valid_bond_length(double dist_sq) {
    for (double allowed_dist_sq : allowed_bond_lengths_squared) {
        if (fabs(dist_sq - allowed_dist_sq) < 1e-6) { // Check if distance is one of the allowed values, 1e-6 for overflow control
            return true;
        }
    }
    return false;
}

// Function to check if a proposed move is valid (self-avoidance and bond length constraints)
bool is_valid_move(const vector<Monomer>& chain1, const vector<Monomer>& chain2, int chain_index, int monomer_index, const Monomer& new_pos) {
    const vector<Monomer>& current_chain = (chain_index == 0) ? chain1 : chain2;
    const vector<Monomer>& other_chain = (chain_index == 0) ? chain2 : chain1;

    if (new_pos.z < Z_MIN) {
        return false;  // z cannot be less than Z_MIN
    }

    // Check bond length constraints within the current chain
    if (monomer_index > 0 && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index - 1]))) return false;
    if (monomer_index < CHAIN_LENGTH - 1 && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index + 1]))) return false;
    
    // Check self-avoidance within the current chain
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        if (i != monomer_index && distance_squared(new_pos, current_chain[i]) < MIN_BOND_LENGTH_SQUARED) {
            return false;  // Overlap detected
        }
    }
     /*
    if (distance_squared(new_pos, other_chain[monomer_index]) < MAX_BOND_LENGTH_SQUARED) {
        return true;    
    } else {
        return false;  // Bond of 2 monomers of same index in different chain breaks
    }
    */

    return true;
}

// Function to calculate the energy change for a monomer move
int calculate_energy_change(const Monomer& current_monomer, const Monomer& other_chain_monomer, const Monomer& new_pos) {
    int energy_change = 0;

    // Chain to wall interaction
    if (current_monomer.z > WALL_DISTANCE && new_pos.z == static_cast<int>(WALL_DISTANCE)) {
        energy_change -= 1;  // Moving towards the wall
    } else if (current_monomer.z == static_cast<int>(WALL_DISTANCE) && new_pos.z > WALL_DISTANCE) {
        energy_change += 1;  // Moving away from the wall
    }

    // Monomer interaction with the same index from the other chain
    double dist_sq = distance_squared(current_monomer, other_chain_monomer);
    double new_dist_sq = distance_squared(new_pos, other_chain_monomer);

    if (dist_sq > WALL_DISTANCE_SQUARED && new_dist_sq == WALL_DISTANCE_SQUARED) {
        energy_change -= 1;  // Getting closer to the other chain's monomer
    } else if (dist_sq == WALL_DISTANCE_SQUARED && new_dist_sq > WALL_DISTANCE_SQUARED) {
        energy_change += 1;  // Moving farther from the other chain's monomer
    }

    return energy_change;
}

// Function to calculate the total energy of the system 
int calculate_total_energy(const vector<Monomer>& chain1, const vector<Monomer>& chain2) {
    int total_energy = 0;

    // Calculate energy for chain1
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        // Chain to wall interaction
        if (chain1[i].z > WALL_DISTANCE) {
            // No interaction
        } else if (chain1[i].z == static_cast<int>(WALL_DISTANCE)) {
            total_energy -= 1;
        }

        // Interaction with other chain's monomer
        if (distance_squared(chain1[i], chain2[i]) == WALL_DISTANCE_SQUARED) {
            total_energy -= 1;
        }
    }

    // Calculate energy for chain2
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        // Chain to wall interaction
        if (chain2[i].z > WALL_DISTANCE) {
            // No interaction
        } else if (chain2[i].z == static_cast<int>(WALL_DISTANCE)) {
            total_energy -= 1;
        }

        // Interaction with other chain's monomer
        if (distance_squared(chain2[i], chain1[i]) == WALL_DISTANCE_SQUARED) {
            total_energy -= 1;
        }
    }

    return total_energy;
}

// Metropolis Monte Carlo move
void monte_carlo_move(vector<Monomer>& chain1, vector<Monomer>& chain2, int& total_energy, vector<int>& energy_history,double BETA) {
    // Consider initial configuration Ci
    // Calculate energy E(Ci) = Ei (already tracked as total_energy)

    // Select random chain and monomer (excluding last monomer)
    int chain_index = rand() % 2;  // 0 for chain1, 1 for chain2
    int monomer_index = rand() % (CHAIN_LENGTH - 1);  // Exclude last monomer

    // Get references to the selected chain and monomer
    vector<Monomer>& selected_chain = (chain_index == 0) ? chain1 : chain2;
    Monomer current_monomer = selected_chain[monomer_index];
    Monomer other_chain_monomer = (chain_index == 0) ? chain2[monomer_index] : chain1[monomer_index];

    // Attempt to move the monomer in a random direction by 1 unit
    Monomer trial_pos = current_monomer;  // Initialize trial position
    int direction = rand() % 4; // Random direction (0: +y, 1: -y, 2: +z, 3: -z)
    
    // Update Y or Z position based on the direction
    switch (direction) {
        case 0: trial_pos.y += 1; break; // Move +y
        case 1: trial_pos.y -= 1; break; // Move -y
        case 2: trial_pos.z += 1; break; // Move +z
        case 3: trial_pos.z -= 1; break; // Move -z
    }

    // Check if the move satisfies SAW and bond length conditions
    if (!is_valid_move(chain1, chain2, chain_index, monomer_index, trial_pos)) {
        // If not satisfied, reject the move and proceed to next step
        return;
    }

    // Calculate energy E(Ct) = Et
    int delta_energy = calculate_energy_change(current_monomer, other_chain_monomer, trial_pos);

    // Decide whether to accept the trial move
    bool accept_move = false;
    if (delta_energy < 0 || delta_energy==0) {
        // Accept the move unconditionally
        accept_move = true;
    } else {
        // Accept the move with probability exp(-beta * delta_energy)
        double probability = exp(-BETA * delta_energy);
        double random_number = static_cast<double>(rand()) / RAND_MAX;  // Random number between 0 and 1
        if (random_number < probability) {
            accept_move = true;
        }
    }

    if (accept_move) {
        // Accept the move
        selected_chain[monomer_index] = trial_pos;
        total_energy += delta_energy;

        // Set Ci+1 = Ct
        // Record energy after this move
        if(BETA==0.5){
        energy_history.push_back(total_energy);
        }
    }
    // If rejected, do not shift the monomer (Ci remains the same)
    // No action needed as the monomer position and total_energy remain unchanged
}

// average energy and distance calculations
void average(vector<Monomer>& chain1, vector<Monomer>& chain2,int& total_energy, double& E, double& E2,double& d1,double& d2,double BETA) {
    // Consider initial configuration Ci
    // Calculate energy E(Ci) = Ei (already tracked as total_energy)

    // Select random chain and monomer (excluding last monomer)
    int chain_index = rand() % 2;  // 0 for chain1, 1 for chain2
    int monomer_index = rand() % (CHAIN_LENGTH - 1);  // Exclude last monomer

    // Get references to the selected chain and monomer
    vector<Monomer>& selected_chain = (chain_index == 0) ? chain1 : chain2;
    Monomer current_monomer = selected_chain[monomer_index];
    Monomer other_chain_monomer = (chain_index == 0) ? chain2[monomer_index] : chain1[monomer_index];

    // Attempt to move the monomer in a random direction by 1 unit
    Monomer trial_pos = current_monomer;  // Initialize trial position
    int direction = rand() % 4; // Random direction (0: +y, 1: -y, 2: +z, 3: -z)
    
    // Update Y or Z position based on the direction
    switch (direction) {
        case 0: trial_pos.y += 1; break; // Move +y
        case 1: trial_pos.y -= 1; break; // Move -y
        case 2: trial_pos.z += 1; break; // Move +z
        case 3: trial_pos.z -= 1; break; // Move -z
    }
    
    // Check if the move satisfies SAW and bond length conditions
    if (!is_valid_move(chain1, chain2, chain_index, monomer_index, trial_pos)) {
         E += static_cast<double>(total_energy);
    E2 += static_cast<double>(total_energy * total_energy);
    d1 += static_cast<double>(chain1[0].z - 0);
    d2 += static_cast<double>(chain2[0].z - 0);// If not satisfied, reject the move and proceed to next step
        return;
    }

    // Calculate energy E(Ct) = Et

    int delta_energy = calculate_energy_change(current_monomer, other_chain_monomer, trial_pos);

    // Decide whether to accept the trial move
    bool accept_move = false;
    if (delta_energy < 0 || delta_energy==0) {
        // Accept the move unconditionally
        accept_move = true;
    } else {
        // Accept the move with probability exp(-beta * delta_energy)
        double probability = exp(-BETA * delta_energy);
        double random_number = static_cast<double>(rand()) / RAND_MAX;  // Random number between 0 and 1
        if (random_number < probability) {
            accept_move = true;
        }
    }

    if (accept_move) {
        // Accept the move
        selected_chain[monomer_index] = trial_pos;
        total_energy += delta_energy;
               // Set Ci+1 = Ct
        
    }
    // If rejected, do not shift the monomer (Ci remains the same)
    // No action needed as the monomer position and total_energy remain unchanged
    E += static_cast<double>(total_energy);
    E2 += static_cast<double>(total_energy * total_energy);
    d1 += static_cast<double>(chain1[0].z - 0);
    d2 += static_cast<double>(chain2[0].z - 0);
}


// Function to save chains' conformation to a file
void save_chains(const vector<Monomer>& chain1, const vector<Monomer>& chain2, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    outfile << "Chain 1 (X = " << X_FIXED_1 << "):" << endl;
    for (const auto& monomer : chain1) {
        outfile << monomer.x << " " << monomer.y << " " << monomer.z << endl;
    }

    outfile << "\nChain 2 (X = " << X_FIXED_2 << "):" << endl;
    for (const auto& monomer : chain2) {
        outfile << monomer.x << " " << monomer.y << " " << monomer.z << endl;
    }

    outfile.close();
}

// Function to save energy history to a file
void save_energy_history(const vector<int>& energy_history, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    for (size_t i = 0; i < energy_history.size(); i++) {
        outfile << i << "," << energy_history[i] << endl;  // Step, Energy
    }

    outfile.close();
}

// Function to save chains' conformation to a CSV file
void save_chains_to_csv(const vector<Monomer>& chain1, const vector<Monomer>& chain2, const string& filename) {
    ofstream outfile(filename);
    if (!outfile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    // CSV header
    outfile << "Chain,Monomer,X,Y,Z" << endl;

    // Save chain 1 coordinates 
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        outfile << "Chain1," << i << "," << chain1[i].x << "," << chain1[i].y << "," << chain1[i].z << endl;
    }

    // Save chain 2 coordinates 
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        outfile << "Chain2," << i << "," << chain2[i].x << "," << chain2[i].y << "," << chain2[i].z << endl;
    }

    outfile.close();
}

// Function to save results to a text file
void save_results_to_txt(const vector<double>& T_values, const vector<double>& E_, const vector<double>& Esquared_,
                         const vector<double>& d1_, const vector<double>& d2_, const vector<double>& Cv_,
                         const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "Temperature\tE\tEsquared\td1\td2\tCv\n";
    for (size_t i = 0; i < T_values.size(); ++i) {
        file << T_values[i] << "\t" << E_[i] << "\t" << Esquared_[i]
             << "\t" << d1_[i] << "\t" << d2_[i] << "\t" << Cv_[i] << "\n";
    }
    file.close();
}

// Function to save results to a CSV file
void save_results_to_csv(const vector<double>& T_values, const vector<double>& E_, const vector<double>& Esquared_,
                         const vector<double>& d1_, const vector<double>& d2_, const vector<double>& Cv_,
                         const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "Temperature,E,Esquared,d1,d2,Cv\n";
    for (size_t i = 0; i < T_values.size(); ++i) {
        file << T_values[i] << "," << E_[i] << "," << Esquared_[i]
             << "," << d1_[i] << "," << d2_[i] << "," << Cv_[i] << "\n";
    }
    file.close();
}

// Main simulation loop
int main() {
    srand(static_cast<unsigned int>(time(0)));

    // Initialize chains
    vector<Monomer> chain1(CHAIN_LENGTH); // Chain 1
    vector<Monomer> chain2(CHAIN_LENGTH); // Chain 2
    
    // Initialize chains by reading from the input file
    string input_file = "polymer_input.txt";
    initialize_chains_from_file(chain1, chain2, input_file);

    // Initialize total energy
    int total_energy = -3 * CHAIN_LENGTH;
    vector<int> energy_history;  // To store the energy after each accepted Monte Carlo move
    
    vector<double> E_,Esquared_,d1_,d2_,Cv_,T_values;
    

    double E=0.00,Esquared=0.00,d1=0.00,d2=0.00;

    cout << "Initial Energy: " << total_energy << endl;

    cout << "Initial Chains' Conformation loaded from " << input_file << endl;
    save_chains(chain1, chain2, "initial_chains.txt");  // Save initial conformation
    save_chains_to_csv(chain1, chain2, "initial_chains.csv");
    
    double T;
    for (T = T_START; T <= T_END+0.0001; T += T_STEP) {
        double BETA = 1.0 / (k * T);
        T_values.push_back(T);
    
    // Perform Monte Carlo simulation using the Metropolis algorithm
    int monte_carlo_steps1 = 50000;
    for (int step = 0; step < monte_carlo_steps1; step++) {
        monte_carlo_move(chain1, chain2, total_energy, energy_history,BETA);
        
    }
    
    //for only 1 T
    // Save final conformation and energy history
    if (BETA==0.5){
    save_chains(chain1, chain2, "equill_final_chains.txt");
    save_energy_history(energy_history, "energy_history.txt");
    
    //for only 1 T
    // Save final conformation to CSV file
    save_chains_to_csv(chain1, chain2, "equill_final_chains.csv");
    cout << "Final Chains(equillibrium)' Conformation saved to final_chains.txt and final_chains.csv" << endl;
    }

    cout << "Final Energy:(equillibrium) " << total_energy << endl;

    int total_energy1 = total_energy;
    cout << "Total Energy as Double: " << total_energy1 << endl;

    int monte_carlo_steps2 = 1000000;
    int step=0;
    for ( step = 0; step < monte_carlo_steps2; step++) {
        average(chain1, chain2,total_energy1, E,Esquared, d1,d2,BETA);
    }

    double steps=static_cast<double>(monte_carlo_steps2);
    E=static_cast<double>(E/steps);
    E_.push_back(E);
    Esquared=static_cast<double>(Esquared/steps);
    Esquared_.push_back(Esquared);
    d1=static_cast<double>(d1/steps);
    d1_.push_back(d1);
    d2=static_cast<double>(d2/steps);
    d2_.push_back(d2);

    double Cv=((Esquared-pow(E,2))/(k*T*T));
    Cv_.push_back(Cv);
    
    cout<<endl;
    cout<<"Average Energy: "<< E <<endl;
    cout<<" Average Esquared: "<< Esquared<<endl;
    cout<<" Average distance of Nth monomer of 1st chain: "<< d1<<endl;
    cout<<" Average distance of Nth monomer of 2nd chain: "<< d2<<endl;
    cout<<" Cv "<< Cv<<endl;

    if (BETA==0.5){
     // Save final conformation and energy history
    save_chains(chain1, chain2, "final_chains.txt");

    // Save final conformation to CSV file
    save_chains_to_csv(chain1, chain2, "final_chains.csv");
    }
     }
    //T loop close
     // Save results to .txt and .csv files after the temperature loop
    save_results_to_txt(T_values, E_, Esquared_, d1_, d2_, Cv_, "results.txt");
    save_results_to_csv(T_values, E_, Esquared_, d1_, d2_, Cv_, "results.csv");

    
    return 0;
}
