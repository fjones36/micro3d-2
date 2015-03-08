
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <string.h>
#include <cmath>
#include "Custom.h"
using namespace std;

vector< vector<float> > genfromtxt(char* fname)
{
	ifstream myfile(fname);
	vector< vector<float> > pts;
	char line[50];
	while( !myfile.eof() ) {
		vector<float> pt;
		myfile >> line;
		pt.push_back( atof(line) );
		myfile >> line;
		pt.push_back( atof(line) );
		pts.push_back(pt);
	}
	return pts;

}


void set_custom_locations(char* fname, float **basemat_odd_x, float **basemat_odd_y, float **basemat_even_x, float **basemat_even_y, 
			  int &num_oddneu_x, int &num_oddneu_y, int &num_eveneu_x, int &num_eveneu_y)
{
	vector<vector<float> > pts = genfromtxt(fname);

	//None fo the even odd stuff matters in this case
	num_oddneu_x = 1;
	num_oddneu_y = 1;
	num_eveneu_x = 1;
	num_eveneu_y = pts.size() - 1;

	int i, j;
	
	// ALLOCATE MEMORY
	basemat_even_x = new float*[num_eveneu_x];
	basemat_even_y = new float*[num_eveneu_x];
	for (i=0; i<num_eveneu_x; i++)
	{
	  basemat_even_x[i] = new float[num_eveneu_y];
	  basemat_even_y[i] = new float[num_eveneu_y];
	}
	basemat_odd_x = new float*[num_oddneu_x];
	basemat_odd_y = new float*[num_oddneu_x];
	for (i=0; i<num_oddneu_x; i++)
	{
	  basemat_odd_x[i] = new float[num_oddneu_y];
	  basemat_odd_y[i] = new float[num_oddneu_y];
	}

	
	for (i=0; i<num_eveneu_x; i++)
	{
	  for (j=0; j<num_eveneu_y; j++)
	  {
		basemat_even_x[i][j] = pts[j][0];
		basemat_even_y[i][j] = pts[j][1];
	  }
	}

	basemat_odd_x[0][0] = pts[pts.size()-1][0];
	basemat_odd_y[0][0] = pts[pts.size()-1][1];
}

