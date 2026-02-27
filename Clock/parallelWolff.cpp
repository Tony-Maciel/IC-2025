/* COMPILE WITH:
   g++ -fopenmp -O3 -march=native -std=c++20 parallelWolff.cpp -o parallelWolff

   RUN WITH:
   ./parallelWolff 1
*/

/*
 * This program splits the N Monte Carlo steps into the default number of threads. 
 *
 * For example, If 10^5 total MC steps were to be performed and there are 10 default threads, 
 * then (once equilibrium is reached), each thread will simulate a copy of the system for 10^4 
 * MC steps. In the end the results are brought together. 
 *
 * jump() is used by the prng to guarantee each thread doesn't generate the same sequence of random 
 * numbers once the parallel part starts. long_jump() is used once the results from each parallel part 
 * are brought together and another temperature gets simulated. This way no same sequence of random numbers 
 * gets used twice.
 *
 * */

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
#include <omp.h>     // parallelize


using RealType = double;
using IntType = long long;

const RealType TWO_PI = std::numbers::pi * 2.0; 

/* xoshiro256+ implementation by David Blackman and Sebastiano Vigna (vigna@acm.org)
 * Public domain.
 *
 * NOTE: the period of this prng is 2^256 - 1. jump() is being used for each thread 
 * at a temperature, then long_jump() is used at a new temperature to guarantee new random 
 * numbers. for a given temperature, there can be a maximum of about 10^19 threads before 
 * repeated random numbers are used , and there can be a maximum of about 10^19 different 
 * temperatures simulated before the same problem happens. In other words, this won't happen...
 *
 */
struct Xoshiro256Plus { 
    uint64_t state[4];

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    // splitmix64 for seeding as recommended by the authors of xoshiro256+
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

    double next_double() {
        return (double)(next() >> 11) * 0x1.0p-53;
    }

    // Advances 2^128 steps. Used for THREADS within a single temperature.
    void jump() {
        static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }

    // Advances 2^192 steps. Used to move to the next TEMPERATURE safely.
    void long_jump() {
        static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe69d };
        uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
        for (int i = 0; i < 4; i++) {
            for (int b = 0; b < 64; b++) {
                if (LONG_JUMP[i] & 1ULL << b) {
                    s0 ^= state[0];
                    s1 ^= state[1];
                    s2 ^= state[2];
                    s3 ^= state[3];
                }
                next();
            }
        }
        state[0] = s0;
        state[1] = s1;
        state[2] = s2;
        state[3] = s3;
    }
};

inline int fast_mod(int a, int b) {
    int r = a % b;
    return r < 0 ? r + b : r;
}

struct System {
    int dim, q, L, Total;
    std::vector<int> Neighbors, Status, stack, cluster, positions;
    std::vector<RealType> cos_table, sin_table;
    std::vector<uint8_t> in_cluster;
    RealType cluster_size;
    Xoshiro256Plus rng; // Each system has its own PRNG

    void init_sys(int L_in, int dim_in, int q_in, bool random_start) {
        L = L_in;
        dim = dim_in;
        q = q_in;
        Total = static_cast<int>(std::pow(L, dim));
        Neighbors.resize(Total * 2 * dim);
        Status.resize(Total);
        positions.resize(dim);
        stack.resize(2 * dim * Total);
        cluster.resize(Total);
        in_cluster.resize(Total, 0); 
        cos_table.resize(q);
        sin_table.resize(q);

        for (int i = 0; i < q; ++i) {
            cos_table[i] = std::cos(RealType(i) * TWO_PI / RealType(q));
            sin_table[i] = std::sin(RealType(i) * TWO_PI / RealType(q));
        }
        if (!random_start) std::fill(Status.begin(), Status.end(), 0); // ordered (cold) start
        else for (int i = 0; i < Total; ++i) Status[i] = int(rng.next_double() * q); // hot start
    }

    int f() const {
        int res = 0; 
        for (int i = 0; i < dim; ++i) res += positions[i] * static_cast<int>(std::pow(L, dim - 1 - i)); 
        return res;
    }

    void set_neighbors_recursive(int depth) {
        if (depth >= dim) {
            int temp = f();
            for (int j = 0; j < dim; ++j) {
                int temp2 = positions[j];
                positions[j] = (temp2 + 1) % L;
                Neighbors[temp * (2*dim) + j] = f();
                positions[j] = temp2;
                positions[j] = (temp2 - 1 + L) % L;
                Neighbors[temp * (2*dim) + (j + dim)] = f();
                positions[j] = temp2;
            }
        } else for (int i = 0; i < L; ++i) { positions[depth] = i; set_neighbors_recursive(depth + 1); }
    }
  
    void set_neighbors() { set_neighbors_recursive(0); }

    void MCstep(const std::vector<RealType>& P_add_table) {
        int idx_cl = 0, idx_st = 0;
        int site = int(rng.next_double() * Total);
        int rand_rot = int(rng.next_double() * (q - 1)) + 1;
        cluster[idx_cl++] = site;
        in_cluster[site] = 1;
        stack[idx_st++] = site;

        while (idx_st > 0) {
            site = stack[--idx_st];
            int state_s = Status[site];
            int base_idx = site * (2 * dim);
            for (int i = 0; i < 2 * dim; ++i) {
                int neighb = Neighbors[base_idx + i];
                if (in_cluster[neighb]) continue;
                int k_old = fast_mod(state_s - Status[neighb], q);
                int k_new = fast_mod((state_s + rand_rot) % q - Status[neighb], q);
                if (rng.next_double() < P_add_table[k_old * q + k_new]) {
                    cluster[idx_cl++] = neighb;
                    in_cluster[neighb] = 1;
                    stack[idx_st++] = neighb;
                }
            }
        }
        for (int i = 0; i < idx_cl; ++i) { Status[cluster[i]] = (Status[cluster[i]] + rand_rot) % q; in_cluster[cluster[i]] = 0; }
        cluster_size += idx_cl;
    }
};

struct Measurements {
    int N, burnin, lag, length;
    RealType Ti, Tf, dT;
    std::vector<RealType> E, M, C, X, mc_E, mc_M, E_err, M_err, C_err, X_err, bincu, m_phi, U_phi, Upsilon;
    RealType E1, E2, M2, M4, M8, Mp, Mp2, Mp4, Ex, Sx2;

    void init_meas(int L, int N_in, int burnin_in, int lag_in, RealType Ti_in, RealType Tf_in, RealType dT_in) {
        N = N_in;
        burnin = burnin_in;
        lag = lag_in;
        Ti = Ti_in;
        Tf = Tf_in;
        dT = dT_in;
        length = int((Tf - Ti) / dT + 1); 
        E.resize(length);
        M.resize(length);
        C.resize(length);
        X.resize(length);
        mc_E.resize(N);
        mc_M.resize(N);
        E_err.resize(length, 0);
        M_err.resize(length, 0);
        C_err.resize(length, 0);
        X_err.resize(length, 0);
        bincu.resize(length);
        m_phi.resize(length);
        U_phi.resize(length);
        Upsilon.resize(length);
    }
};

void Get(System& sys, Measurements& meas, int mcs_idx) {
    RealType Mx = 0.0, My = 0.0, E_temp = 0.0, Exprime = 0.0, Sxprime = 0.0;
    for (int site = 0; site < sys.Total; ++site) {
        int state_s = sys.Status[site];
        Mx += sys.cos_table[state_s];
        My += sys.sin_table[state_s];
        int n_right = sys.Neighbors[site * (2*sys.dim) + 0];
        int diff_r = fast_mod(state_s - sys.Status[n_right], sys.q);
        E_temp += sys.cos_table[diff_r];
        Exprime += sys.cos_table[diff_r];
        Sxprime += sys.sin_table[diff_r];
        int n_up = sys.Neighbors[site * (2*sys.dim) + 1];
        E_temp += sys.cos_table[fast_mod(state_s - sys.Status[n_up], sys.q)];
    }
    meas.E1 = -E_temp / RealType(sys.Total);
    meas.E2 = meas.E1 * meas.E1;
    meas.M2 = std::sqrt(Mx*Mx + My*My) / RealType(sys.Total);
    meas.M4 = meas.M2 * meas.M2;
    meas.M8 = meas.M4 * meas.M4;
    meas.Mp = std::cos(sys.q * std::atan2(My, Mx));
    meas.Mp2 = meas.Mp * meas.Mp;
    meas.Mp4 = meas.Mp2 * meas.Mp2;
    meas.Ex = Exprime / RealType(sys.Total);
    meas.Sx2 = (Sxprime * Sxprime) / RealType(sys.Total);
    meas.mc_E[mcs_idx] = meas.E1;
    meas.mc_M[mcs_idx] = meas.M2;
}

void Run(System& sys, Measurements& meas, RealType beta, int idx) {
    std::vector<RealType> P_add_table(sys.q * sys.q);
    for (int i = 0; i < sys.q; ++i) for (int j = 0; j < sys.q; ++j) {
        RealType ed = sys.cos_table[i] - sys.cos_table[j];
        P_add_table[i * sys.q + j] = (ed > 0.0) ? (1.0 - std::exp(-beta * ed)) : 0.0;
    }

    // Accumulators for global reduction
    RealType Et = 0, Mt = 0, Ct = 0, Xt = 0, bc = 0, mphi = 0, mphi2 = 0, mphi4 = 0, mExp = 0, mSx2 = 0;

    // Sequential Burn-in (Annealing preservation)
    for (int i = 0; i < meas.burnin; ++i) for (int j = 0; j < meas.lag; ++j) sys.MCstep(P_add_table);

    // Parallel Measurement Phase
    #pragma omp parallel
    {
        System l_sys = sys; // Isolated deep copy
        int tid = omp_get_thread_num();
        for(int k=0; k < tid; ++k) l_sys.rng.jump(); // Unique sequence per thread
        
        Measurements l_meas = meas; 
        RealType le = 0, lm = 0, lc = 0, lx = 0, lbc = 0, lp = 0, lp2 = 0, lp4 = 0, lex = 0, lsx = 0;

        #pragma omp for
        for (int i = 0; i < meas.N; ++i) {
            for (int j = 0; j < meas.lag; ++j) l_sys.MCstep(P_add_table);
            Get(l_sys, l_meas, i);
            le += l_meas.E1;
            lm += l_meas.M2;
            lc += l_meas.E2;
            lx += l_meas.M4;
            lbc += l_meas.M8;
            lp += l_meas.Mp;
            lp2 += l_meas.Mp2;
            lp4 += l_meas.Mp4;
            lex += l_meas.Ex;
            lsx += l_meas.Sx2;
        }

        #pragma omp critical
        {
            Et += le;
            Mt += lm;
            Ct += lc;
            Xt += lx;
            bc += lbc;
            mphi += lp;
            mphi2 += lp2;
            mphi4 += lp4;
            mExp += lex;
            mSx2 += lsx;
        }
    }

    RealType norm = 1.0 / RealType(meas.N);
    Et *= norm;
    Mt *= norm;
    Ct *= norm;
    Xt *= norm;
    bc *= norm;
    mphi *= norm;
    mphi2 *= norm;
    mphi4 *= norm;
    mExp *= norm;
    mSx2 *= norm;

    meas.E[idx] = Et;
    meas.M[idx] = Mt; 
    meas.C[idx] = (Ct - Et * Et) * beta * beta * sys.Total;
    meas.X[idx] = (Xt - Mt * Mt) * beta * sys.Total;
    meas.bincu[idx] = 1.0 - bc / (3.0 * Xt * Xt);
    meas.U_phi[idx] = 1.0 - mphi4 / (2.0 * mphi2 * mphi2);
    meas.m_phi[idx] = mphi;
    meas.Upsilon[idx] = mExp - mSx2 * beta;
}

int main(int argc, char** argv) {
    if (argc < 2) return 1;
    int int_param = std::stoi(argv[1]);
    int dim, q, init_L, n_sizes, mcs, burnin, TNC, nboot, lag, seed;
    bool random_start;
    RealType Ti, Tf, dT;

    std::ifstream infile("clock_input.in");
    infile >> dim >> q >> init_L >> n_sizes >> random_start >> Ti >> Tf >> dT >> mcs >> burnin >> TNC >> nboot >> lag >> seed;

    System sample;
    Measurements vals;
    sample.rng.seed(seed + int_param);
    int L = init_L;

    for (int k = 0; k < n_sizes; ++k) {
        std::string output = "EMCX_" + std::to_string(L) + "_" + std::to_string(q) + "-" + std::to_string(int_param) + ".dat";
        std::ofstream outfile(output);
        outfile << "# L=" << L << " q=" << q << " MCsteps=" << mcs << " burnin=" << burnin << "\n";
        outfile << "# T / E / E_err / M / M_err / C / C_err / X / X_err / BC / Uphi / mphi / Upsilon\n";
        outfile.close();

        lag = L * L; 
        vals.init_meas(L, mcs, burnin, lag, Ti, Tf, dT);
        sample.init_sys(L, dim, q, random_start);
        sample.set_neighbors();

        RealType current_temp = vals.Ti;
        for (int i = 0; i < vals.length; ++i) {
            std::cout << "L=" << L << " T=" << current_temp << std::endl;
            Run(sample, vals, 1.0 / current_temp, i);

            std::ofstream out_append(output, std::ios::app);
            out_append << std::scientific << std::setprecision(17) << current_temp << " "
                       << vals.E[i] << " " << vals.E_err[i] << " " << vals.M[i] << " " << vals.M_err[i] << " "
                       << vals.C[i] << " " << vals.C_err[i] << " " << vals.X[i] << " " << vals.X_err[i] << " " 
                       << vals.bincu[i] << " " << vals.U_phi[i] << " " << vals.m_phi[i] << " " << vals.Upsilon[i] << "\n";
            
            sample.rng.long_jump(); // Move to a completely fresh PRNG space for the next temperature 
            current_temp += vals.dT;
        }
        L += init_L;
    }
    return 0;
}
