#include "calculations.h"


vector<int> get_s(double w, int d){
	const double range_from  = 0;
    const double range_to    = w;
    
    vector<int> sVec;
    
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_int_distribution<> distr(range_from, range_to);
 
    for (int i = 0; i < d; ++i){
        sVec.push_back(distr(generator));
    }
	return sVec;
}


vector<int> calculate_a(vector<unsigned char> pVec, vector<int> sVec, double w, int d){			//an argei na to kanoume ena-ena
	vector<int> aVec;
	float temp;
	
	for(int i = 0; i < d; i++){
		temp = ((float)pVec[i] + w - (float)sVec[i])/w;
		aVec.push_back(floor(temp));
	}
	return aVec;
}


int modular_pow(int base, int exponent, int modulator){
    int result;

    base = base % modulator;

    while(exponent > 0){
        if( (exponent%2) != 0){                       // odd number
            result = (result * base) % modulator;
        }
        exponent = exponent >> 1;
        base = (base * base) % modulator;
    }
    return result;
}


int calculate_h(vector<int> aVec, int m, int M, int d){
	int h = 0, j, x;
	j = d-1;
	int temp;
	for(int i = 0; i < d; i++){									//modulo
		
		if((aVec[j]%M) < 0)
			temp = aVec[j]%M + M;
		else
			temp = aVec[j]%M;
			
		x = (temp* modular_pow(m, i, M))%M;
		
		// if(x < 0)
			// x += M;
		
		h += x;
		j--;
	}
	//cout << endl << "==============================" << endl << endl;
	h = h % M;
	// cout << "h = " << h << endl;
	return h;
}

unsigned int manhattan_dist(vector<unsigned char> qVec, vector<unsigned char> pVec, int d){
	unsigned int dist = 0;
	
	for(int i = 0; i < d; i++){
		dist+= abs(qVec[i] - pVec[i]);
	}
	return dist;
}

vector<distanceNode> actual_nearest_neighbor(vector<unsigned char>  qVec, vector< vector< unsigned char> > pVec, int d, int N){
	unsigned int temp;							
	distanceNode node;
	vector<distanceNode> distances;

	for(int i = 0; i < N; i++){
		node.pPos = -1;
		node.dist = 4294967295;                        //highest possible unsigned int
		distances.push_back(node);							
	}                            
			
	for( int j = 0; j < pVec.size(); j++){
		temp = manhattan_dist(qVec, pVec[j], d);
		if(temp < distances[N-1].dist){
			
			distances[N-1].dist = temp;
			distances[N-1].pPos = j;
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
	return distances;
}

float get_x(float top){
	float x;
	
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_real_distribution<float> distribution (0.0, top);
    
    x = distribution(generator);
	
	return x;
}
