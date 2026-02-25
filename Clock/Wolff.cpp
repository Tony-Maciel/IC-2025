// (MUST use .FALSE. or .TRUE. in fortran input file and 0 or 1 in C++ input file) 

/* COMPILE WITH:
   FLAGS="-O3 -march=native -std=c++20 -flto"
   g++ $FLAGS -DCOMPILER_FLAGS="\"$FLAGS\"" Wolff.cpp -o Wolff

   RUN WITH:
   ./Wolff 1
*/
//NOTE: if .FALSE. is in clock_input.in, error will occur. Replace this with 0 (or 1 if .TRUE.)

#include <iostream>  // input and output
#include <vector>    // Dynamic arrays
#include <cmath>     // sin, cos, exp, pow, sqrt, atan2
#include <fstream>   // Reading and writing to files
#include <string>    // Unique file names 
#include <iomanip>   // set precision of output data
#include <algorithm> // std::fill (a little faster than for loop)
#include <cstdint>   // fixed width integer types (in prng and Wolff)
#include <numbers>   // pi constant
#include <chrono>    // getting current time
#include <ctime>     // making time from chrono more human readable
#include <sstream>   // treats strings as files 


using RealType = double;
using IntType = long long;

const RealType TWO_PI = std::numbers::pi * 2.0; // portability...

/* xoshiro256+ implementation by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * Public domain.
 */
struct Xoshiro256Plus {
    uint64_t state[4];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    // Seeding logic using SplitMix64
    void seed(uint64_t seed_val) {
        uint64_t z = seed_val;
        for (int i = 0; i < 4; i++) {
            z += 0x9e3779b97f4a7c15;
            uint64_t s = z;
            s = (s ^ (s >> 30)) * 0xbf58476d1ce4e5b9;
            s = (s ^ (s >> 27)) * 0x94d049bb133111eb;
            state[i] = s ^ (s >> 31);
        }
    }

    uint64_t next() {
        const uint64_t result = state[0] + state[3];
        const uint64_t t = state[1] << 17;

        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];

        state[2] ^= t;
        state[3] = rotl(state[3], 45);

        return result;
    }

    // Generates a random double in [0, 1)
    double next_double() {
        // 0x1.0p-53 is 1.0 / 2^53
        return (double)(next() >> 11) * 0x1.0p-53;
    }
};

// Global RNG instance
Xoshiro256Plus rng;

// Equivalent to Fortran-style modulo (no negative numbers)
inline int fast_mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

/*
 * Represents the system being simulated.
 * */
struct System {
    int dim;                         // dimension of system
    int q;                           // number of possible directions to point
    int L;                           // length of one dimension
    int Total;                       // Total number of constituents
    std::vector<int> Neighbors;      // Flattened 2D array: Index = site * (2*dim) + direction
    std::vector<int> Status;         // Integer states
    std::vector<RealType> cos_table; // precomputed lookup tables
    std::vector<RealType> sin_table; 
    std::vector<int> stack;          // data structures for Wolff algorithm
    std::vector<int> cluster;
    std::vector<uint8_t> in_cluster; // Using uint8_t is faster than vector<bool>
    std::vector<int> positions;      // Recursive neighbor setup helper
    RealType cluster_size;           // mean cluster size

    // Initializes system 
    void init_sys(int L_in, int dim_in, int q_in, bool random_start) {
        L     = L_in;
        dim   = dim_in;
        q     = q_in;
        Total = static_cast<int>(std::pow(L, dim));

        // Allocate arrays
        Neighbors.resize(Total * 2 * dim);
        Status.resize(Total);
        positions.resize(dim);
        stack.resize(2 * dim * Total);
        cluster.resize(Total);
        in_cluster.resize(Total, 0); 
        cos_table.resize(q);
        sin_table.resize(q);

        // Populate lookup tables
        for (int i = 0; i < q; ++i) {
            cos_table[i] = std::cos(RealType(i) * TWO_PI / RealType(q));
            sin_table[i] = std::sin(RealType(i) * TWO_PI / RealType(q));
        }

        if (!random_start) { // ordered start
            std::fill(Status.begin(), Status.end(), 0); 
        } else {             // random start
            for (int i = 0; i < Total; ++i) {
                Status[i] = int(rng.next_double() * q); 
            }
        }
    }

    // Auxiliary function to set up neighbors
    int f() const {
        int res = 0; 
        for (int i = 0; i < dim; ++i) {
            res += positions[i] * static_cast<int>(std::pow(L, dim - 1 - i)); 
        }
        return res;
    }

    // Recursive neighbor setup for any dimension
    void set_neighbors_recursive(int depth) {
        if (depth >= dim) {
            int temp = f();
            
            for (int j = 0; j < dim; ++j) {
                int temp2 = positions[j];
                
                // Right Neighbor
                positions[j] = (temp2 + 1) % L;
                Neighbors[temp * (2*dim) + j] = f();
                positions[j] = temp2; // Reset

                // Left Neighbor 
                positions[j] = (temp2 - 1 + L) % L;
                Neighbors[temp * (2*dim) + (j + dim)] = f();
                positions[j] = temp2; // Reset
            }
        } else {
            for (int i = 0; i < L; ++i) {
                positions[depth] = i;
                set_neighbors_recursive(depth + 1);
            }
        }
    }
  
    // Need to start recursive function somehow, so call this
    void set_neighbors() {
        set_neighbors_recursive(0);
    }

    // saves snapshot of system
    void save_config(const std::string& output_file) {
        std::ofstream io(output_file);
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                io << Status[i * L + j] << (j == L - 1 ? "" : "  ");
            }
            io << "\n";
        }
    }

    // A single cluster flip with the Wolff algorithm (optimized)
    void MCstep(const std::vector<RealType>& P_add_table) {
        int idx_cl = 0; // Cluster pointer
        int idx_st = 0; // Stack pointer
        
        int site = int(rng.next_double() * Total);           // random site
        int rand_rot = int(rng.next_double() * (q - 1)) + 1; // random rotation

        // add site to stack and cluster
        cluster[idx_cl++] = site;
        in_cluster[site] = 1; 
        stack[idx_st++] = site;

        while (idx_st > 0) {            // while still elements in stack...
            site = stack[--idx_st];     // pop one...
            int state_s = Status[site];

            int base_idx = site * (2 * dim);
            for (int i = 0; i < 2 * dim; ++i) { // for all its neighbors...
                int neighb = Neighbors[base_idx + i];
                
                if (in_cluster[neighb]) continue;  // if its already in the cluster, do nothing...

                int state_n = Status[neighb];
                int state_s_prime = (state_s + rand_rot) % q; 
                
                int k_old = fast_mod(state_s - state_n, q);
                int k_new = fast_mod(state_s_prime - state_n, q);

                if (rng.next_double() < P_add_table[k_old * q + k_new]) { // else, add to stack and cluster with this probability
                    cluster[idx_cl++] = neighb;
                    in_cluster[neighb] = 1;
                    stack[idx_st++] = neighb;
                }
            }
        }

        // rotate all sites in cluster
        for (int i = 0; i < idx_cl; ++i) {
            site = cluster[i];
            Status[site] = (Status[site] + rand_rot) % q;
            in_cluster[site] = 0; // reset
        }
        cluster_size += idx_cl; // count for mean cluster size
    }
};

/*
 * Represents the person doing the measurements on a given system
 * */
struct Measurements {
    int N, TNC, burnin, lag, length, nboot;
    RealType Ti, Tf, dT;
    
    std::vector<RealType> E, M, C, X, mc_E, mc_M, E_err, M_err, C_err, X_err;
    std::vector<RealType> corr, bincu, m_phi, U_phi, Upsilon;
    RealType E1, E2, M2, M4, M8;
    RealType Mp, Mp2, Mp4, Ex, Sx2;

    // Initializes persons notes
    void init_meas(int L, int N_in, int burnin_in, int lag_in, 
                   RealType Ti_in, RealType Tf_in, RealType dT_in, 
                   int TNC_in, int nboot_in) {
        N = N_in;
        burnin = burnin_in;
        lag = lag_in;
        Ti = Ti_in;
        Tf = Tf_in;
        dT = dT_in;
        TNC = TNC_in;
        nboot = nboot_in;
        
        length = int((Tf - Ti) / dT + 1); 
        
        E.resize(length); M.resize(length); C.resize(length); X.resize(length);
        mc_E.resize(N); mc_M.resize(N);
        E_err.resize(length); M_err.resize(length);
        C_err.resize(length); X_err.resize(length);
        corr.resize(L/2);
        bincu.resize(length); m_phi.resize(length);
        U_phi.resize(length); Upsilon.resize(length);
    }
};

// Makes the person collect measurements on the system at a given MC step
void Get(System& sys, Measurements& meas, int mcs_idx) {
    RealType Mx = 0.0, My = 0.0;
    RealType E_temp = 0.0, Exprime = 0.0, Sxprime = 0.0;

    for (int site = 0; site < sys.Total; ++site) {
        int state_s = sys.Status[site];

        // Magnetization
        Mx += sys.cos_table[state_s];
        My += sys.sin_table[state_s];

        // Right neighbor 
        int neighb = sys.Neighbors[site * (2*sys.dim) + 0];
        int state_n = sys.Status[neighb];
        int diff = fast_mod(state_s - state_n, sys.q);
        
        RealType temp = sys.cos_table[diff];
        E_temp += temp;
        Exprime += temp;
        Sxprime += sys.sin_table[diff];

        // Up neighbor 
        neighb = sys.Neighbors[site * (2*sys.dim) + 1];
        state_n = sys.Status[neighb];
        diff = fast_mod(state_s - state_n, sys.q);
        E_temp += sys.cos_table[diff];
    }

    meas.E1 = -E_temp / RealType(sys.Total); // energy per site
    meas.E2 = meas.E1 * meas.E1;
    meas.M2 = std::sqrt(Mx*Mx + My*My) / RealType(sys.Total); // magnetization per site
    meas.M4 = meas.M2 * meas.M2;
    meas.M8 = meas.M4 * meas.M4;
    meas.Mp = std::cos(sys.q * std::atan2(My, Mx)); // angular magnetization
    meas.Mp2 = meas.Mp * meas.Mp;
    meas.Mp4 = meas.Mp2 * meas.Mp2;
    meas.Ex = Exprime / RealType(sys.Total); // for spin stiffness
    meas.Sx2 = (Sxprime * Sxprime) / RealType(sys.Total); // for spin stiffness too

    meas.mc_E[mcs_idx] = meas.E1;
    meas.mc_M[mcs_idx] = meas.M2;
}

// Makes person simulate system for many MC steps, take measurements of it and save them
void Run(System& sys, Measurements& meas, RealType beta, int idx, int int_param) {
    // Create Probability Table
    std::vector<RealType> P_add_table(sys.q * sys.q);
    
    for (int i = 0; i < sys.q; ++i) {
        for (int j = 0; j < sys.q; ++j) {
            RealType energy_diff = sys.cos_table[i] - sys.cos_table[j];
            if (energy_diff > 0.0) {
                P_add_table[i * sys.q + j] = 1.0 - std::exp(-beta * energy_diff);
            } else {
                P_add_table[i * sys.q + j] = 0.0;
            }
        }
    }

    // Accumulators
    RealType Et = 0, Mt = 0, Ct = 0, Xt = 0, bc = 0;
    RealType mphi = 0, mphi2 = 0, mphi4 = 0, mExp = 0, mSx2 = 0;
    std::fill(meas.corr.begin(), meas.corr.end(), 0.0);
    sys.cluster_size = 0.0;

    // Burn-in time
    for (int i = 0; i < meas.burnin; ++i) {
        for (int j = 0; j < meas.lag; ++j) {
            sys.MCstep(P_add_table);
        }
    }

    sys.cluster_size = 0.0; // reset mean cluster size because it may be biased depending 
                            // on how the system was started (cold or hot start). What matters 
                            // is the mean cluster size when taking measurements.

    // Measurement time
    for (int i = 0; i < meas.N; ++i) {
        for (int j = 0; j < meas.lag; ++j) {
            sys.MCstep(P_add_table);
        }
        Get(sys, meas, i);

        Et += meas.E1;
        Mt += meas.M2; 
        Ct += meas.E2;
        Xt += meas.M4;
        bc += meas.M8;
        mphi += meas.Mp;
        mphi2 += meas.Mp2; 
        mphi4 += meas.Mp4;
        mExp += meas.Ex;
        mSx2 += meas.Sx2;
    }

    // Complete averages
    RealType norm = 1.0 / RealType(meas.N);
    Et *= norm; Mt *= norm; Ct *= norm; Xt *= norm;
    bc *= norm; mphi *= norm; mphi2 *= norm; mphi4 *= norm;
    mSx2 *= norm; mExp *= norm;

    meas.m_phi[idx] = mphi;
    
    meas.E[idx] = Et;
    meas.M[idx] = Mt;
    meas.C[idx] = (Ct - Et * Et) * beta * beta * sys.Total;
    meas.X[idx] = (Xt - Mt * Mt) * beta * sys.Total;
    meas.E_err[idx] = 0.0;
    meas.M_err[idx] = 0.0;
    meas.bincu[idx] = 1.0 - bc / (3.0 * Xt * Xt);
    meas.U_phi[idx] = 1.0 - mphi4 / (2.0 * mphi2 * mphi2);
    meas.C_err[idx] = 0.0;
    meas.X_err[idx] = 0.0;
    meas.Upsilon[idx] = mExp - mSx2 * beta;

    sys.cluster_size /= (RealType(meas.N) * RealType(meas.lag));
}

// Returns current date and time as a string
std::string get_current_datetime() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_time), "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

// Returns compiler info
std::string get_compiler_info() {
    std::stringstream ss;
    ss << "Compiler: " << __VERSION__;
#ifdef __OPTIMIZE__
    ss << " | Optimization: Enabled";
#else
    ss << " | Optimization: Disabled";
#endif
    return ss.str();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <int_param>\n";
        return 1;
    }

    int int_param = std::stoi(argv[1]);
    
    // Input variables
    int dim, q, init_L, n_sizes, mcs, burnin, TNC, nboot, lag, seed;
    bool random_start;
    RealType Ti, Tf, dT;

    // Read Input
    std::ifstream infile("clock_input.in");
    if (!infile) {
        std::cerr << "Error opening clock_input.in\n";
        return 1;
    }
    infile >> dim >> q >> init_L >> n_sizes >> random_start
           >> Ti >> Tf >> dT >> mcs >> burnin >> TNC >> nboot >> lag >> seed;
    infile.close();

    // Initialize Xoshiro RNG with unique seed (for that node)
    rng.seed(seed + int_param);

    System sample;
    Measurements vals;
    int L = init_L;

    for (int k = 0; k < n_sizes; ++k) {
        std::cout << "Simulating L = " << L << std::endl; 
        
        std::string output = "EMCX_" + std::to_string(L) + "_" + std::to_string(q) + "-" + std::to_string(int_param) + ".dat";
        
        // If file exists, skip it. DON'T overwrite
        std::ifstream check(output);
        if (check.good()) {
            std::cout << output << " exists already. Skipping.\n";
            check.close();
            L += init_L;
            continue;
        }
        check.close();

        lag = L * L; // definition of MC step is L^2 cluster flips with Wolff algorithm
        vals.init_meas(L, mcs, burnin, lag, Ti, Tf, dT, TNC, nboot);
        sample.init_sys(L, dim, q, random_start);
        sample.set_neighbors();

        // Write Header
        std::ofstream outfile(output);
        outfile << "# C++ Version (Xoshiro256+)\n";
        outfile << "# " << get_current_datetime() << "\n";
        outfile << "# " << get_compiler_info() << "\n";
        #ifdef COMPILER_FLAGS
        outfile << "# Build Flags: " << COMPILER_FLAGS << "\n";
        #endif
        outfile << "# dimension = " << sample.dim << "\n";
        outfile << "# q = " << sample.q << "\n";
        outfile << "# L = " << sample.L << "\n";
        outfile << "# Seed = " << seed + int_param << "\n";
        outfile << "# temperature interval (Ti, Tf, dT) = " << Ti << " " << Tf << " " << dT << "\n";
        outfile << "# random start = " << random_start << "\n";
        outfile << "# Burn in time = " << vals.burnin << "\n";
        outfile << "# Lag time = " << vals.lag << "\n";
        outfile << "# Number of MC steps = " << vals.N << "\n";
        outfile << "# Temperature / E / E error / M / M error / C / C error / X / X error / BC / Uphi / mphi / Upsilon\n";
        outfile.close();

        // Init optional correlation file (empty shell to match logic)
        // std::string corr_output = "CORR_" + std::to_string(L) + "_" + std::to_string(q) + "-" + std::to_string(int_param) + ".dat";
        // std::ofstream corr_file(corr_output);
        // corr_file.close();

        RealType current_temp = vals.Ti;
        for (int i = 0; i < vals.length; ++i) { // simulate in all these temperaturees
            std::cout << "Progress: " << i+1 << "/" << vals.length << " (T = " << current_temp << ")\n";
            
            Run(sample, vals, 1.0 / current_temp, i, int_param);

            // Append results to output file
            std::ofstream out_append(output, std::ios::app);
            out_append << std::scientific << std::setprecision(17);
            out_append << current_temp << " "
                       << vals.E[i] << " "
                       << vals.E_err[i] << " " 
                       << vals.M[i] << " "
                       << vals.M_err[i] << " "
                       << vals.C[i] << " " 
                       << vals.C_err[i] << " "
                       << vals.X[i] << " "
                       << vals.X_err[i] << " " 
                       << vals.bincu[i] << " "
                       << vals.U_phi[i] << " "
                       << vals.m_phi[i] << " " 
                       << vals.Upsilon[i] << "\n";
            out_append.close();

            current_temp += vals.dT;
        }

        L += init_L;
    }

    return 0;
}
