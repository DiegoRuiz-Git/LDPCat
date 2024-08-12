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
    context c;
    params p(c);
    p.set("mul2concat", false);
    tactic t = 
        with(tactic(c, "simplify"), p) &
        tactic(c, "solve-eqs") &
        tactic(c, "bit-blast") &
        tactic(c, "aig") &
        tactic(c, "sat");
    solver s = t.mk_solver();

    //parameters of the problem
    int H = 8; // Height of the lattice of data qubits
    int L = 5;  // Width of the lattice of data qubits
    int d = 21; // Target for the distance
   
    // array storing SAT variables
    // 6 booleans per shape
    // size (H-2)*6
    expr_vector SyndromeRows(c);

    // create stabilizer variables
    for (int i = 0; i < H-2; ++i){
        expr_vector CurrentSyndrome(c);
        for (int j = 0; j < 6; ++j) {
            std::stringstream x_name;
            x_name << "b" << j << '_' << i;
            SyndromeRows.push_back(c.bool_const(x_name.str().c_str()));
            CurrentSyndrome.push_back(c.bool_const(x_name.str().c_str()));
        }
        // constraint for at most weight-4 syndromes
        s.add(atmost(CurrentSyndrome,3));
    }

    // loop on codewords
    for (int i = 1; i < (1 << (2 * L)); i++) {

        // array to store the codeword
        // size H*L
        expr_vector Codeword(c);

        std::bitset<32> binary = std::bitset<32>(i);

        // initialize bottom row of the codeword
        for (int m = 0; m < L; m++) {
            if (binary[m] == 1)
                Codeword.push_back(c.bool_val(true));
            else
                Codeword.push_back(c.bool_val(false));
        }
        // initialize second bottom row of the codeword
        for (int m = 0; m < L; m++) {
            if (binary[m+L] == 1)
                Codeword.push_back(c.bool_val(true));
            else
                Codeword.push_back(c.bool_val(false));
        }

	    // fill codeword depending on the stabilizer
        for (int n = 2; n < H; n++) {
            for (int m = 0; m < L; m++) {
                Codeword.push_back(
                    (Codeword[(n - 1) * L + ((m - 1 + L) % L)] && SyndromeRows[(n - 2) * 6 + 0])^
                    ((Codeword[(n - 1) * L + m] && SyndromeRows[(n - 2) * 6 + 1])^
                    ((Codeword[(n - 1) * L + ((m + 1) % L)] && SyndromeRows[(n - 2) * 6 + 2])^
                    ((Codeword[(n - 2) * L + ((m - 1 + L) % L)] && SyndromeRows[(n - 2) * 6 + 3])^
                    ((Codeword[(n - 2) * L + m] && SyndromeRows[(n - 2) * 6 + 4])^
                     (Codeword[(n - 2) * L + ((m + 1) % L)] && SyndromeRows[(n - 2) * 6 + 5]))))));
            }
        }

        // Constaint of Hamming weight >= d
        s.add(atleast(Codeword,d));
    }

    clock_t start = clock();

    // If the SAT problem is unsolvable
    if (s.check() == unsat){
        cout << "unsat" << endl;
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
    for (int i = 0; i < H-2; ++i){
        cout << " 1 "<< endl;
        for (int j = 0; j < 3; ++j) {
            expr result = m.eval(SyndromeRows[6*i + j]);
            if (result.is_true())
                cout << 1;
            else
                cout << 0;
        }
        cout << endl;
        for (int j = 3; j < 6; ++j) {
            expr result = m.eval(SyndromeRows[6*i + j]);
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