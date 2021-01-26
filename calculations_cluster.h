#ifndef __CALCULATIONS_CLUSTER_H__
#define __CALCULATIONS_CLUSTER_H__

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
#include "calculations_cube.h"
#include "help_functions.h"

using namespace std;


void k_means_init(vector< vector<unsigned char> > &centroids, int number_of_images, vector< vector<unsigned char> > pVec, int k, int d);
void lloyds_assignment(vector< vector<int> > &clusters, vector< vector<int> > temp, int number_of_images, vector< vector<unsigned char> > pVec, vector< vector<unsigned char> > centroids, int k, int d, int *changes, int first);
void update_centroids_median(vector< vector<unsigned char> > &centroids, vector <unsigned char> pDim, vector< vector<unsigned char> > pVec, vector< vector<int> > clusters, vector <unsigned char> tempC, int k, int d);

vector<distanceNode> approximate_range_search_clusterLSH(vector < vector<unsigned char> > centroids, vector < vector< vector <hTableNode> > > &lHashTables, int L, int pos, int d, double R, int cluster);
vector<distanceNode> approximate_range_search_clusterCube(vector < vector<unsigned char> > centroids, vector< vector <hTableNode> > &hashTable, int pos, int d, double R, int M, int probes, int cluster);

// void write_results(vector< vector<int> > clusters, vector< vector<unsigned char> > centroids, auto duration, int d, int k, ofstream &ofile);
void silhouette(vector< vector<int> > clusters, vector< vector<unsigned char> > centroids, vector< vector<unsigned char> > pVec, int k, int d, ofstream &ofile);

#endif
