#ifndef __CALCULATIONS_CUBE_H__
#define __CALCULATIONS_CUBE_H__

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <bits/stdc++.h>
#include <string>
#include <list> 
#include <iterator> 
#include <random>
#include <vector>
#include <cmath>
#include <vector>
#include <chrono>
#include <algorithm>
#include <climits>

#include "calculations.h"

using namespace std;


class fNode{
public:
	int h;
	bool f;
};

unsigned int calculate_p(vector<fNode> fVec, int k);
vector<int> hamming_dist(int probes, int p);

int get_f();

void create_hashtable_cube(vector< vector <hTableNode> > &hashTable, vector< vector<unsigned char> > pVec, vector< vector<int> > &sVec, int hTableSize, int number_of_images, double w, int k, int d, int m, int M);

vector<distanceNode> approximate_nearest_neighbor_cube(vector<unsigned char> qVec, vector< vector <hTableNode> > hashTable, int p, int d, int N, int M, int probes);
vector<distanceNode> approximate_range_search_cube(vector<unsigned char> qVec, vector< vector <hTableNode> > hashTable, int p, int d, double R, int M, int probes);


#endif
