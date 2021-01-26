#include "calculations_cube.h"

unsigned int calculate_p(vector<fNode> fVec, int k){
	unsigned int p=0;

	for (int i=0; i<k; i++){
		p = fVec[i].f << (k-1-i) | p ;
	}
	return p;
}


int get_f(){
	int f;
	
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> distribution (0, 1);
    
    f = distribution(generator);
	return f;
}


void create_hashtable_cube(vector< vector <hTableNode> > &hashTable, vector< vector<unsigned char> > pVec, vector< vector<int> > &sVec, int hTableSize, int number_of_images, double w, int k, int d, int m, int M){
	int h;
	unsigned int g, p;
	bool exists;
	hTableNode node;
	fNode fnode;

	vector<int> aVec, tempIntVec;
	vector< vector<fNode> > fVec; 
	vector<fNode> tempfVec;
	
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
				for (int j=0; j<fVec.size(); j++){
					for (int p=0; p<fVec[j].size(); p++){
						if (h == fVec[j][p].h){
							exists = 1;
							fnode.h = h;
							fnode.f = fVec[j][p].f;
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
		fVec.push_back(tempfVec);                      // save f*k distinct of every image
		tempfVec.erase(tempfVec.begin(), tempfVec.end());
		
		p = calculate_p(fVec[i], k);

		node.pPos = i;
		node.pVec = pVec[i];
		node.flag = 0;
		hashTable[p].push_back(node);            // insert image in the hash table
	}
	
}

vector<int> hamming_dist(int probes, int p){
	vector<int> hp;
	int x, count, mask;


	for (int i=0; i<probes; i++){
		mask = 1 << i;
		if (p & mask){
			x = (p & ~mask) | ((0 << i) & mask);
		}
		else{
			x = (p & ~mask) | ((1 << i) & mask);
		}
		hp.push_back(x);
	}
	return hp;
}


vector<distanceNode> approximate_nearest_neighbor_cube(vector<unsigned char> qVec, vector< vector <hTableNode> > hashTable, int p, int d, int N, int M, int probes){
	unsigned int temp;
	distanceNode node;
	vector<distanceNode> distances;
	vector<int> ph;
	int count=0;
	
	ph = hamming_dist(probes-1, p);                           // Find positions with Hamming distance 1
	ph.insert(ph.begin(), p);                                 // Put p position in position 0 to search first

	for(int i = 0; i < N; i++){
		node.pPos = -1;
		node.dist = 4294967295;                        //highest possible unsigned int
		distances.push_back(node);	
	}
	
	for(int i=0; i<ph.size(); i++){                           // For every position with Hamming distance 1 until it searches M images
		for(int j=0; j < hashTable[ph[i]].size(); j++){
			temp = manhattan_dist(qVec, hashTable[ph[i]][j].pVec, d);
			if(temp < distances[N-1].dist){
				distances[N-1].dist = temp;
				distances[N-1].pPos = hashTable[ph[i]][j].pPos;
				for(int c=N-2; c>=0; c--){
					if(distances[c].dist > distances[c+1].dist){
						iter_swap(distances.begin() + c, distances.begin() + c+1);
					}
					else{
						break;
					}
				}
			}
			count++;
			if (count >= M){
				break;
			}
		}
		if (count >= M){
			break;
		}
	}
	return distances;
 }


vector<distanceNode> approximate_range_search_cube(vector<unsigned char> qVec, vector< vector <hTableNode> > hashTable, int p, int d, double R, int M, int probes){
	unsigned int temp;
	vector<distanceNode> distances;
	vector<int> ph;
	distanceNode node;
	int count=0;
	
	ph = hamming_dist(probes-1, p);                           // Find positions with Hamming distance 1
	ph.insert(ph.begin(), p);                                 // Put p position in position 0 to search first

	for(int i=0; i<ph.size(); i++){                           // For every position with Hamming distance 1 until it searches M images
		for(int j=0; j < hashTable[ph[i]].size(); j++){
			temp = manhattan_dist(qVec, hashTable[ph[i]][j].pVec, d);
			if(temp < R){
				node.pPos = hashTable[ph[i]][j].pPos;
				node.dist = temp;
				distances.push_back(node);
			}
			count++;
			if (count >= M){
				break;
			}
		}
		if (count >= M){
			break;
		}
	}

	return distances;
}
