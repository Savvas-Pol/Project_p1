#include "help_functions.h"
#include "calculations.h"
#include "calculations_cube.h"

#define m 107					//a_max < m < M/2
#define NForTable 16

using namespace std;


int main(int argc, char** argv){
	string iFile, qFile, oFile;
	int i, k, probes, N;
	double R, w;
	int magic_number=0, number_of_images=0;
    int n_rows=0, n_cols=0;
    int d, M, h, pos, p;
    unsigned int g;
    int hTableSize;
    bool exists;

	vector< vector<unsigned char> > pVec, qVec; 
	vector<unsigned char> tempVec;
	vector< vector<int> > sVec;
	vector<int> aVec, tempIntVec;
	vector<fNode> tempfVec;
	vector< vector<fNode> > fVec, qfVec;
	vector<distanceNode> distCube, distTrue, distRange;

	hTableNode node;
	fNode fnode;
	
	read_inputCube(&argc, argv, &iFile, &qFile, &k, &M, &probes, &oFile, &N, &R, &w);
	
	srand (time(NULL));
	
	ifstream file (iFile);
    if (file.is_open()){
	    read_data(file, &magic_number, &number_of_images, &n_rows, &n_cols, pVec, tempVec);
			
		d = n_rows * n_cols;						           // dimension
	    hTableSize = pow(2,k);
	    
	    vector< vector <hTableNode> > hashTable;       // hash table
        create_hashtable_cube(hashTable, pVec, sVec, hTableSize, number_of_images, w, k, d, m, M);

		ifstream qfile (qFile);
		if (qfile.is_open()){
			read_data(qfile, &magic_number, &number_of_images, &n_rows, &n_cols, qVec, tempVec);

			ofstream ofile (oFile);
			if (ofile.is_open()){
				
				for(int i = 0; i < number_of_images; i++){
					for (int j = 0; j < k; j++){
						aVec = calculate_a(qVec[i], sVec[j], w, d);  // calculate a for every image
						h = calculate_h(aVec, m, M, d);              // calculate h for every image
						
						exists = 0;
						for (int y=0; y<tempfVec.size(); y++){
							if (h == tempfVec[y].h){
								exists = 1;
								fnode.h = h;
								fnode.f = tempfVec[y].f;
								break;
							}
						}
						if(exists == 0){
							for (int j=0; j<qfVec.size(); j++){
								for (int p=0; p<qfVec[j].size(); p++){
									if (h == qfVec[j][p].h){
										exists = 1;
										fnode.h = h;
										fnode.f = qfVec[j][p].f;
										break;
									}
								}
								if (exists){
									break;
								}
							}
						}
						
						if (exists == 0){
							fnode.h = h;
							fnode.f = get_f();
						}
						tempfVec.push_back(fnode);
					}
					
					qfVec.push_back(tempfVec);                    		  // save f*k distinct of every image
					tempfVec.erase(tempfVec.begin(), tempfVec.end());
					
					p = calculate_p(qfVec[i], k);
					
					auto t1 = chrono::high_resolution_clock::now();
					distCube = approximate_nearest_neighbor_cube(qVec[i], hashTable, p, d, N, M, probes);
					auto t2 = chrono::high_resolution_clock::now();
					auto durationCube = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();

					auto t3 = chrono::high_resolution_clock::now();
					distTrue = actual_nearest_neighbor(qVec[i], pVec, d, N);
					auto t4 = chrono::high_resolution_clock::now();
					auto durationTrue = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();

					distRange = approximate_range_search_cube(qVec[i], hashTable, p, d, R, M, probes);

					ofile  << "Query: " << i << endl;               // write to file
					for (int j=0; j<N; j++){
						ofile << "Nearest neighbor-N: " << distCube[j].pPos << endl;
						ofile << "distanceHyperCube: " << distCube[j].dist << endl;
						ofile << "distanceTrue: " << distTrue[j].dist << endl;
					}
					ofile << "tHyperCube: " << durationCube << endl;
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
