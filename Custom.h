

#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <string.h>
#include <cmath>
using namespace std;

vector< vector<float> > genfromtxt(char* fname);

void set_custom_locations(char* fname, float **basemat_odd_x, float **basemat_odd_y, float **basemat_even_x, float **basemat_even_y, 
			  int &num_oddneu_x, int &num_oddneu_y, int &num_eveneu_x, int &num_eveneu_y);
