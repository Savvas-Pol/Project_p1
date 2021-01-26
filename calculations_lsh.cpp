#include "calculations_lsh.h"


unsigned int calculate_g(vector<int> hVec, int k){
	unsigned int g=0;
	
	for (int i=0; i<k; i++){
		g = hVec[i] << (k-1-i)*8 | g ;
	}
	return g;
}

void create_hashtables_LSH(vector < vector< vector <hTableNode> > > &lHashTables, vector< vector <hTableNode> > &hashTable, vector< vector<unsigned char> > pVec, int L, int hTableSize, int k, int d, int number_of_images, double w, int m, int M){
	int h, pos;
	unsigned int g;
	hTableNode node;
	
	vector<int> aVec, tempIntVec;
	vector< vector<int> > sVec, hVec;
	
	for (int l=0; l<L; l++){
		for(int y=0; y<hTableSize; y++){
			hashTable.push_back(vector<hTableNode>()); //initialize the first index with a 2D vector
		}

		for (int i=0; i<k; i++){
			tempIntVec = get_s(w, d);                     //s_i uniform random generator
			sVec.push_back(tempIntVec);
			tempIntVec.erase(tempIntVec.begin(), tempIntVec.end());
		}
		
		for (int i=0; i < number_of_images; i++){
			for (int j = 0; j < k; j++){
				aVec = calculate_a(pVec[i], sVec[j], w, d);  // calculate a for every image
				h = calculate_h(aVec, m, M, d);              // calculate h for every image
				// cout << h << endl;
				tempIntVec.push_back(h);
			}
			hVec.push_back(tempIntVec);                      // save k*h of every image
			tempIntVec.erase(tempIntVec.begin(), tempIntVec.end());
			
			g = calculate_g(hVec[i], k);                  // calculate g for every image
			pos = g % hTableSize;                         // find the position to insert the image in the hash table

			node.pPos = i;
			node.g = g;
			node.pVec = pVec[i];
			node.flag = 0;
			hashTable[pos].push_back(node);            // insert image in the hash table
		}
		lHashTables.push_back(hashTable);
		hashTable.erase(hashTable.begin(), hashTable.end());
		hVec.erase(hVec.begin(), hVec.end());
		sVec.erase(sVec.begin(), sVec.end());
	}
}

vector<distanceNode> approximate_nearest_neighbor(vector<unsigned char> qVec, vector < vector< vector <hTableNode> > > lHashTables, int L, int pos, int d, int N, unsigned int g){
	unsigned int temp;
	distanceNode node;
	vector<distanceNode> distances;
	bool flag;
	int count = 0;

	for(int i = 0; i < N; i++){
		node.pPos = -1;
		node.dist = 4294967295;                        //highest possible unsigned int
		distances.push_back(node);							
	}
	
	for(int i = 0; i < L; i++){
		for(int j = 0; j < lHashTables[i][pos].size(); j++){
			if( lHashTables[i][pos][j].g == g){													//check value of g
				temp = manhattan_dist(qVec, lHashTables[i][pos][j].pVec, d);
				if(temp < distances[N-1].dist){
					flag = 0;
					for (int c=0; c<N; c++){
						if (distances[c].pPos == lHashTables[i][pos][j].pPos){
							flag = 1;
						}
					}
					if (flag == 0){
						count++;
						distances[N-1].dist = temp;
						distances[N-1].pPos = lHashTables[i][pos][j].pPos;
						for(int c=N-2; c>=0; c--){
							if(distances[c].dist > distances[c+1].dist){
								iter_swap(distances.begin() + c, distances.begin() + c+1);
							}
							else{
								break;
							}
						}
					}
				}
			}
		}
	}
	if(count < N){																				//if N-nearest neighbors not found, then put random neighbors
		for(int i = 0; i < L; i++){
			for(int j = 0; j < lHashTables[i][pos].size(); j++){
				distances[count].dist = manhattan_dist(qVec, lHashTables[i][pos][j].pVec, d);
				distances[count].pPos = lHashTables[i][pos][j].pPos;
				count++;
				if(count >= N)
					break;
			}
			if(count >= N)
				break;
		}
	}
	return distances;
}


vector<distanceNode> approximate_range_search(vector<unsigned char> qVec, vector < vector< vector <hTableNode> > > lHashTables, int L, int pos, int d, double R, unsigned int g){
	unsigned int temp;
	vector<distanceNode> distances;
	distanceNode node;
	
	for(int i = 0; i < L; i++){
		for(int j = 0; j < lHashTables[i][pos].size(); j++){
			if( lHashTables[i][pos][j].g == g){	
				temp = manhattan_dist(qVec, lHashTables[i][pos][j].pVec, d);
				if(temp < R){
					node.pPos = j;
					node.dist = temp;
					distances.push_back(node);
				}
			}
		}
	}
	return distances;
 }
