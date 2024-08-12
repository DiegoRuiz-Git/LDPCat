// Author: Diego Ruiz
// Affiliation: Alice&Bob - INRIA
// Date: 2023

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <omp.h>
#include <cstdlib>

using namespace std;

// max size of the parity check matrix
#define MAX 1200

// Compute the code parameters of a specific stabilizer shape for a given lattice size
int CodeCyclique(int l, int m, int b,vector<tuple<vector<vector<int>>,int,int,float>>& result,
		bool stabilizer[MAX][MAX],
		bool H[MAX][MAX],bool G[MAX][MAX]);

// Obtain the integer representation of the vertically mirrored stabilizer shape
int VerticalmirrorInt(int b);

// Obtain the integer representation of the horizontally mirrored stabilizer shape
int HorizontalmirrorInt(int b);

// Create the stabilizer shape from its integer representation
void CreateStabilisateur(bool stabilizer[MAX][MAX], int b);

// Create the parity check matrix from the stabilizer shape
void CreateH(bool H[MAX][MAX],bool stabilizer[MAX][MAX],int n,int l,int m);

// Put parity check matrix in normal form step 1
void HTrianglesup(bool H[MAX][MAX],int n);

// Put parity check matrix in normal form step 2
int HRemoveEmptyLine(bool H[MAX][MAX],int nb_rows, int n);

// Put parity check matrix in normal form step 3
int MakeDiagonalOne(bool H[MAX][MAX],int nb_rows,int n);

// Put parity check matrix in normal form step 4
void MakeTheIdentity(bool H[MAX][MAX],int nb_row,int n);

// Put parity check matrix in normal form step 5
void RevertTheIdentity(bool H[MAX][MAX],int nb_row,int n);

// Create the generator matrix from the parity check matrix
void CreateG(bool G[MAX][MAX], bool H[MAX][MAX],int nb_row,int n);

// Calculate the distance of the code
int CalculateDistance(bool G[MAX][MAX],int k,int n);

// Save the result
void SaveResult(vector<tuple<vector<vector<int>>,int,int,float>>&result,bool stabilizer[MAX][MAX],
		int k,int d,float gain,int l, int m);

// Export the result
void ExportResult(int l, int m,vector<tuple<vector<vector<int>>,int,int,float>>& result);

int main(int argc, char* argv[]) {

	// vector of tuple containing stabilizer, k, d, gain = kd/n
	vector< tuple< vector<vector<int>>, int, int, float > > result;

	// argument vector of the main.exe
	int l = stoi(argv[1]); // width of the lattice of data qubits
    int m = stoi(argv[2]); // height of the lattice of data qubits

	cout << "l " << l << " m " << m << endl;

	// Get the number of available threads for multiprocessing
	int num_threads = omp_get_max_threads();
	cout << "Maximum number of threads: " << num_threads << endl;

	// iterate on the stabilizer shape with multiprocessing
	#pragma omp parallel for schedule(dynamic, 1)
    for (int b = 1; b < (1<<9); b ++){

		// skip this shape if it is equivalent to another calculated shape
		if (b > VerticalmirrorInt(b) || b > HorizontalmirrorInt(b) || b > HorizontalmirrorInt(VerticalmirrorInt(b))){
			#pragma omp critical
			{
			cout << "Thread " << omp_get_thread_num() << " not handling iteration " << b << endl;
			}
		}

		else{
			auto stabilizer = new bool [MAX][MAX]; // stabilizer shape
			auto H = new bool [MAX][MAX]; // parity check matrix
			auto G = new bool [MAX][MAX]; // generator matrix
		
			#pragma omp critical
			{
			cout << "Thread " << omp_get_thread_num() << " handling iteration " << b << endl;
			}

			// Calculate the code parameters associated to this shape
			CodeCyclique(l,m,b,result,stabilizer,H,G);
			
			#pragma omp critical
			{
			cout << "Thread " << omp_get_thread_num() << " finished iteration " << b << endl;
			}
		}

	}

	//export result
	ExportResult(l,m,result);

	return 0;

}

// Obtain the integer representation of the vertically mirrored stabilizer shape
int VerticalmirrorInt(int b) {
    int mirrored = 0;
    for (int i = 0; i < 9; i++) {
        int bit = (b >> i) & 1;
        int mirroredPosition = (i / 3) * 3 + (2 - (i % 3));
        mirrored |= bit << mirroredPosition;
    }
    return mirrored;
}

// Obtain the integer representation of the horizontally mirrored stabilizer shape
int HorizontalmirrorInt(int b) {
    int mirrored = 0;
    for (int i = 0; i < 9; i++) {
        int bit = (b >> i) & 1;
        int mirroredPosition = (2 - (i / 3)) * 3 + (i % 3);
        mirrored |= bit << mirroredPosition;
    }
    return mirrored;
}

// Compute the code parameters of a specific stabilizer shape for a given lattice size
int CodeCyclique(int l, int m, int b,vector<tuple<vector<vector<int>>,int,int,float>>& result,
		bool stabilizer[MAX][MAX],
		bool H[MAX][MAX],bool G[MAX][MAX]){

	CreateStabilisateur(stabilizer,b);

	int n = m*l; // number of data qubits

	if (n>MAX){
		cout << "increase MAX to " << n << endl;
		return 1;
	}

	// Create the parity check matrix from the stabilizer shape
	CreateH(H,stabilizer,n,l,m);

	// Put parity check matrix in normal form step 1
	HTrianglesup(H,n);

	// Put parity check matrix in normal form step 2
	int nb_row = HRemoveEmptyLine(H,n,n);

	// Put parity check matrix in normal form step 3
	nb_row = MakeDiagonalOne(H,nb_row,n); // number of independant rows of the parity check matrix

	// Put parity check matrix in normal form step 4
	MakeTheIdentity(H,nb_row,n);

	// Put parity check matrix in normal form step 5
	RevertTheIdentity(H,nb_row,n);

	int k = n - nb_row; // number of logical qubits

	// Create the generator matrix from the parity check matrix
	CreateG(G,H,nb_row,n);

    #pragma omp critical
	{
	cout << "Thread " << omp_get_thread_num() << " b " << b << " k = " << k << endl;
	}

	if (k == 22){
		int d = 0; // distance of the code
		d = CalculateDistance(G,k,n);
		float gain = ((float)(d*k))/((float)n);
 		#pragma omp critical
		{
		SaveResult(result,stabilizer,k,d,gain,l,m);
		}
	}

	return 0;
}

// Create the stabilizer shape from its integer representation
void CreateStabilisateur(bool stabilizer[MAX][MAX], int b){

	// Initialize H with zeros
	for (int i = 0; i < MAX; i++) {
		for (int j = 0; j < MAX; j++)
			stabilizer[i][j] = 0;
	}

	for (int i = 0; i < 9; i ++)
		stabilizer[i/3][i%3] = (b >> i) & 1;
}

// Create the parity check matrix from the stabilizer shape
void CreateH(bool H[MAX][MAX],bool stabilizer[MAX][MAX],int n, int l,int m){

	// Initialize H with zeros
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			H[i][j] = 0;
	}

	// Check if the stabilizer has empty rows
	bool EmptyFirstLine = true;
	bool EmptySecondLine = true;
	bool EmptyLastLine = true;

	for (int j = 0; j < 3; j++) {
		if (stabilizer[0][j] != 0)
			EmptyFirstLine = false;
		if (stabilizer[1][j] != 0)
			EmptySecondLine = false;
		if (stabilizer[2][j] != 0)
			EmptyLastLine = false;
	}

	// if the stabilizer has a first or a last empty row
	// we cover an extra row so that every data qubit is 
	// part of a stabilizer
	int min_i;
	int max_i;

	if (EmptyFirstLine && EmptySecondLine)
		min_i = -2;
	else if (EmptyFirstLine)
		min_i = -1;
	else
		min_i = 0;


	if (EmptyLastLine && EmptySecondLine)
		max_i = l;
	else if (EmptyLastLine)
		max_i = l-1;
	else
		max_i = l-2;

	// Fill H
	for (int i = min_i; i < max_i; i++) {
		for (int j = 0; j < m; j++){

			// shift stabilizer
			int shift_stab[l][m];
			for (int x = 0; x < l; x++){
				for (int y = 0; y < m; y++)
					shift_stab[x][y] = stabilizer[(x-i+l)%l][(y-j+m)%m];
			}

			// flatten stabilizer
			int shift_stab_flatten[n];
			for (int x = 0; x < l; x++){
				for (int y = 0; y < m; y++)
					shift_stab_flatten[x*m+y] = shift_stab[x][y];
			}

			// Put the stabilzer in H
			for (int x = 0; x < n; x++)
				H[m*(i-min_i)+j][x] = shift_stab_flatten[x];

		}
	}
}

// Put parity check matrix in normal form step 1
// Make H an upper triangular matrix
void HTrianglesup(bool H[MAX][MAX],int n){

	for (int i = 0; i < n; i++){

		// look for a row with a 1 in the first column
		int index = i;
		while (H[index][i]==0){
			index+=1;
			if (index == n)
				break;
		}

		if (index != n){

		    // Swap the found row at the top
			for (int j = 0; j < n; j++)
				swap(H[i][j], H[index][j]);

			// sum the found row to all the other rows that 
			// have a 1 in this column
			for (int j = i+1; j < n; j++){
				if (H[j][i] == 1){
					for (int k = 0; k < n; k++)
						H[j][k] = (H[j][k] + H[i][k])%2;
				}
			}

		}
	}

}

// Put parity check matrix in normal form step 2
// remove rows with only zeros from H
int HRemoveEmptyLine(bool H[MAX][MAX],int nb_max_rows,int n){

	int newH[n][n];

	// Initialize newH with zeros
	for (int i = 0; i < nb_max_rows; i++) {
		for (int j = 0; j < n; j++)
			newH[i][j] = 0;
	}

	int nb_row = 0;

	// Copy non zero row of H in newH
	for (int i = 0; i < nb_max_rows; i++) {

		bool isZeroRow = true;
		for (int j = 0; j < n; j++) {
			if (H[i][j] != 0) {
				isZeroRow = false;
				break;
			}
		}

		if (!isZeroRow) {
			for (int j = 0; j < n; j++) {
				newH[nb_row][j] = H[i][j];
			}
			nb_row++;
		}
	}

	// copy newH into H
	for (int i = 0; i < nb_row; i++) {
		for (int j = 0; j < n; j++)
			H[i][j] = newH[i][j];
	}

	return nb_row;
}

// Put parity check matrix in normal form step 3
// Put 1 on the diagonal of the matrix
int MakeDiagonalOne(bool H[MAX][MAX],int nb_rows,int n){

	int i = 0;
	while (i < nb_rows){
		if (H[i][i] != 1){

			// Find a column with a 1 on row i
			int index = i;
			while (index<n-1 and H[i][index] != 1)
				index+=1;

			// Swap the columns
			for (int j = 0; j < nb_rows; j++)
				swap(H[j][i], H[j][index]);

			// Remove 1 in rows below
			for (int j = i+1; j < nb_rows;j++){
				if (H[j][i] == 1){
					for (int k = 0; k < n; k++)
						H[j][k] = (H[j][k] + H[i][k])%2;
				}
			}

			// Remove lines equal to zero
			nb_rows = HRemoveEmptyLine(H,nb_rows,n);

		}
		i+=1;

	}

	return nb_rows;
}

// Put parity check matrix in normal form step 4
// Put the matrix in the form (I,M) where I is the identity
void MakeTheIdentity(bool H[MAX][MAX],int nb_row,int n){
	for (int i = nb_row - 2; i> -1; i--){
		for (int j = i+1;j<nb_row;j++){
			if (H[i][j] == 1){
				for (int k = 0; k < n; k++)
					H[i][k] = (H[i][k] + H[j][k])%2;
			}
		}
	}
}

// Put parity check matrix in normal form step 5
// Put the matrix in the form (M,I) where I is the identity
void RevertTheIdentity(bool H[MAX][MAX],int nb_row,int n){
	for (int k = 0; k < nb_row; ++k) {
		for (int j = 0; j < n - 1; ++j) {
			for (int i = 0; i < nb_row; ++i) {
				std::swap(H[i][j], H[i][j + 1]);
			}
		}
	}
}

// Create the generator matrix from the parity check matrix
void CreateG(bool G[MAX][MAX], bool H[MAX][MAX],int nb_row,int n){

	// Initialize G with zeros
	for (int i = 0; i < n - nb_row; i++) {
		for (int j = 0; j < n; j++)
			G[i][j] = 0;
	}

	// Setting the identity matrix part in G
	for (int i = 0; i < n - nb_row; ++i) {
			G[i][i] = 1;
	}

	// Transposing a subset of columns from H to G
	for (int i = 0; i < nb_row; ++i) {
		for (int j = 0; j < n - nb_row; ++j) {
			G[j][n - nb_row + i] = H[i][j];
		}
	}

}

// Calculate the distance of the code
int CalculateDistance(bool G[MAX][MAX],int k,int n){

	// distance of the code initialize to the number of data qubits
	int d = n;

	/// Save G as a sparse matrix to speed up calculation

	// store the number of 1 in every row of G
	int NumberOf1InGLine[k] = {};
	for (int i = 0; i < k; i++){
		for (int j = 0; j < n; j++){
			if (G[i][j])
				NumberOf1InGLine[i] += 1;
		}
	}

	// store the locations of 1 in every row of G
	int Location1InG[k][n] = {};
	for (int i = 0; i < k; i++){
		int compt = 0;
		for (int j = 0; j < n; j++){
			if (G[i][j]){
				Location1InG[i][compt] = j;
				compt +=1;
			}
		}
	}

	// Compute Hamming weight of every codeword
	for (long long i = 1; i < (1LL << k); i++){

		int codeword[n] = {};

		// compute the codeword
		for (int j = 0; j < k; j++){
			if ((i >> j) & 1){
				for (int m = 0; m < NumberOf1InGLine[j]; m++)
					codeword[Location1InG[j][m]] = (codeword[Location1InG[j][m]] + 1)%2;
			}
		}

		// compute the Hamming weight
		int HammingWeight = 0;
		for (int m = 0; m < n; m++)
			HammingWeight += codeword[m];

		// Update distance
		if (HammingWeight<d)
			d = HammingWeight;
	}

	return d;

}

// Save the result i.e. the code parameters
void SaveResult(vector<tuple<vector<vector<int>>,int,int,float>>&result,bool stabilizer[MAX][MAX],
		int k,int d,float gain,int l, int m){

	tuple<vector<vector<int>>,int,int,float> MyTuple;

	// Save stabilizer shape
	vector<vector<int>> stab(l, std::vector<int>(m));
	for (int i = 0; i < l; ++i) {
		for (int j = 0; j < m; ++j)
			stab[i][j] = stabilizer[i][j];
	}

	get<0>(MyTuple) = stab;
	// Save number of logical qubits, distance, and kd/n
	get<1>(MyTuple) = k;
	get<2>(MyTuple) = d;
	get<3>(MyTuple) = gain;

	result.push_back(MyTuple);
}

// Export result in a txt file
void ExportResult(int l,int m,vector<tuple<vector<vector<int>>,int,int,float>>& result){
	std::ofstream outfile("./result/result"+to_string(l)+"_"+to_string(m)+".txt");

	for (const auto &entry : result) {
		// Write stabilizer in the file
		const auto& stab = std::get<0>(entry);
		for (const auto &row : stab) {
			for (const auto &elem : row) {
				outfile << elem << " ";
			}
			outfile << "\n"; // New line after each row
		}
		outfile << "---\n"; // Separator

		// Write k, d, and gain to file
		outfile << std::get<1>(entry) << "\n";
		outfile << std::get<2>(entry) << "\n";
		outfile << std::get<3>(entry) << "\n";
	}

	outfile.close();

	cout << "finish" << endl;
}

