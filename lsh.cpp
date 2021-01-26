#include "help_functions.h"
#include "calculations.h"
#include "calculations_lsh.h"

#define m 107					//a_max < m < M/2
#define NForTable 16

using namespace std;


int main(int argc, char** argv){
	string iFile, qFile, oFile;
	int i, k, L, N;
	double R, w;
	int magic_number=0, number_of_images=0;
    int n_rows=0, n_cols=0;
    int d, M, h, pos;
    unsigned int g;
    int hTableSize, probes;

	vector< vector<unsigned char> > pVec, qVec;
	vector<unsigned char> tempVec;
	vector< vector<int> > sVec;
	vector<int> aVec, tempIntVec;
	vector<distanceNode> distLsh, distTrue, distRange;
	

	read_inputLSH(&argc, argv, &iFile, &qFile, &k, &L, &oFile, &N, &R, &w);
	
	M = pow(2,floor(32/k));
	
	ifstream file (iFile);
    if (file.is_open()){
        read_data(file, &magic_number, &number_of_images, &n_rows, &n_cols, pVec, tempVec);
			
		d = n_rows * n_cols;						           // dimension
        hTableSize = number_of_images / NForTable;
        
        vector < vector< vector <hTableNode> > > lHashTables;       // vector with L hash tables
        vector< vector <hTableNode> > hashTable;       // hash table

		create_hashtables_LSH(lHashTables, hashTable, pVec, L, hTableSize, k, d, number_of_images, w, m, M);

		ifstream qfile (qFile);
	    if (qfile.is_open()){
	        read_data(qfile, &magic_number, &number_of_images, &n_rows, &n_cols, qVec, tempVec);

	    	ofstream ofile (oFile);
			if (ofile.is_open()){
				for (int i=0; i<k; i++){
					tempIntVec = get_s(w, d);                     //s_i uniform random generator
					sVec.push_back(tempIntVec);
					tempIntVec.erase(tempIntVec.begin(), tempIntVec.end());
				}

				for(int i = 0; i < number_of_images; i++){
					for (int j = 0; j < k; j++){
						aVec = calculate_a(qVec[i], sVec[j], w, d);  // calculate a for every image
						h = calculate_h(aVec, m, M, d);              // calculate h for every image
						tempIntVec.push_back(h);
					}
					g = calculate_g(tempIntVec, k);                  // calculate g for every image
					pos = g % hTableSize;                         // find the position to insert the image in the hash table
					tempIntVec.erase(tempIntVec.begin(), tempIntVec.end());

					auto t1 = chrono::high_resolution_clock::now();
					distLsh = approximate_nearest_neighbor(qVec[i], lHashTables, L, pos, d, N, g);
					auto t2 = chrono::high_resolution_clock::now();
					auto durationLsh = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();

					auto t3 = chrono::high_resolution_clock::now();
					distTrue = actual_nearest_neighbor(qVec[i], pVec, d, N);
					auto t4 = chrono::high_resolution_clock::now();
					auto durationTrue = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();

					distRange = approximate_range_search( qVec[i], lHashTables, L, pos, d, R, g);

					ofile  << "Query: " << i << endl;               // write to file
					for (int j=0; j<N; j++){
						ofile << "Nearest neighbor-" << j+1 << ": " << distLsh[j].pPos << endl;
						ofile << "distanceLSH: " << distLsh[j].dist << endl;
						ofile << "distanceTrue: " << distTrue[j].dist << endl;
					}
					ofile << "tLSH: " << durationLsh << endl;
					ofile << "tTrue: " << durationTrue << endl;
					ofile << "R-near neighbors:" << endl;
					if (distRange.size() == 0){
						ofile << "No neighbors in this range." << endl;
					}
					else{
						for(int c=0; c<distRange.size(); c++){
							ofile << distRange[c].pPos << endl;
						}
					}
					ofile << endl;
				}
			}
		}
		else{
			cout << "Could not open query file." << endl;
		}
    }
    else {
    	cout << "Could not open input file." << endl;
    }
	return 0;
}
