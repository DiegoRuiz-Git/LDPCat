// Author: Diego Ruiz
// Affiliation: Alice&Bob - INRIA
// Date: 2023

#include <iostream>
#include <bitset>
#include <iostream>
#include <sstream>
#include<vector>
#include <ctime>
#include "z3++.h"

using namespace z3;
using namespace std;

int main() {

    // Define the parameters of the solver
    z3::context c;
    z3::tactic tactic1 = z3::tactic(c, "simplify");
	z3::tactic tactic2 = z3::tactic(c, "solve-eqs");
	z3::tactic tactic3 = z3::tactic(c, "aig");
	z3::tactic tactic4 = z3::tactic(c, "reduce-bv-size");
    z3::tactic tactic5 = z3::tactic(c, "bit-blast");
	z3::tactic tactic6 = z3::tactic(c, "pb2bv");
	z3::tactic tactic7 = z3::tactic(c, "sat");
    z3::tactic solver_tactic = tactic1 & tactic2 & tactic3 & tactic4 & tactic5 & /*tactic6 &*/ tactic7;
    z3::solver s = solver_tactic.mk_solver(); 

    //parameters of the problem
    int H = 7; // Height of the lattice of data qubits
    int nb_interf = 7; // nb data qubits in the bottom row
    int d = 17; // Target for the distance

    int L = 2*H - 3 + nb_interf - 1; // number of data qubits in the top row
   
    // array storing SAT variables
    // 6 booleans per shape
    // size (H-2)*6
    expr_vector Syndrom(c);
        
    // create stabilizer variables
    for (int i = 0; i < H-2; i++){
        expr_vector CurrentSyndrome(c);
        for (int j = 0; j < 6; j++) {
            std::stringstream x_name;
            x_name << "b_" << i << '_'<< j;
            Syndrom.push_back(c.bool_const(x_name.str().c_str()));
            CurrentSyndrome.push_back(c.bool_const(x_name.str().c_str()));
        }
        // constraint for at most weight-4 syndromes
        s.add(atmost(CurrentSyndrome,4));
    }

    // loop on codewords
    for (int m = 1; m < (1 << (2*nb_interf)); m++) {

        std::bitset<32> binary = std::bitset<32>(m);

        // Initialize Codeword
        // size H*L
        expr_vector Codeword(c);

        for (int i = 0; i < H-2; i++)
            Codeword.push_back(c.bool_val(false));
     
        for (int i = 0; i < nb_interf; i++)
            Codeword.push_back(c.bool_val(binary[i] == 1));

        for (int i = H-2+nb_interf; i < L; i++)
            Codeword.push_back(c.bool_val(false));

        for (int i = 0; i < H-2 ; i++)
            Codeword.push_back(c.bool_val(false));

        for (int i = 0; i < nb_interf; i++)
            Codeword.push_back(c.bool_val(binary[i+nb_interf] == 1));

        for (int i = H-2+nb_interf; i < L; i++)
            Codeword.push_back(c.bool_val(false));

        // fill codeword depending on the stabilizer

        for (int i = 2; i < H; i++)
        {
            for (int j = 0; j < L; j++) {
                Codeword.push_back(
                     (Codeword[(i - 1) * L + (j - 1 + L)%L] && Syndrom[(i - 2) * 6 + 0])^
                    ((Codeword[(i - 1) * L + j]             && Syndrom[(i - 2) * 6 + 1])^
                    ((Codeword[(i - 1) * L + (j + 1)%L]     && Syndrom[(i - 2) * 6 + 2])^
                    ((Codeword[(i - 2) * L + (j - 1 + L)%L] && Syndrom[(i - 2) * 6 + 3])^
                    ((Codeword[(i - 2) * L + j]             && Syndrom[(i - 2) * 6 + 4])^
                     (Codeword[(i - 2) * L + (j + 1L)%L]    && Syndrom[(i - 2) * 6 + 5])))))
                );
            }  
        }

        // Constaint of Hamming weight >= d
        s.add(atleast(Codeword,d));

    }

    clock_t start = clock();

    // If the SAT problem is unsolvable
    if (s.check() == unsat){
        cout << "unsat" << endl;
        clock_t end = clock();
        double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        std::cout << "Duration: " << duration << " seconds" << std::endl;
        return 0;
    }

    // the SAT problem is solvable
    cout << "sat" << endl;
    model m = s.get_model();
    std::cout << m << "\n";

    clock_t end = clock();
    double duration = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::cout << "Duration: " << duration << " seconds" << std::endl;

    // Print the found stabilizer shapes
    for (int i = H-3; i >= 0; i--){
        cout << " 1 "<< endl;
        for (int j = 0; j < 3; ++j) {
            expr result = m.eval(Syndrom[6*i + j]);
            if (result.is_true())
                cout << 1;
            else
                cout << 0;
        }
        cout << endl;
        for (int j = 3; j < 6; ++j) {
            expr result = m.eval(Syndrom[6*i + j]);
            if (result.is_true())
                cout << 1;
            else
                cout << 0;
        }
        cout << endl;
        cout << endl;
    }

    return 0;
}