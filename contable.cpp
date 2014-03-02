/*
 *  contable.cpp
 *
 *  Created by SETO Ryohei on 07/01/18.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "contable.h"
#include <iostream>
#include <iomanip>

using namespace std;
#define DELETE(x) if(x){delete [] x; x = NULL;}

ConTable::~ConTable(){
	if (allocate){
		for (int i=0; i<n+2; ++i){
			DELETE(tbl[i]);
		}
		DELETE(tbl);
	}
}

void ConTable::set(int num_of_particle){
	allocate = true;
	n = num_of_particle;
	tbl = new bool* [n+2];
	for (int i=0; i < n+2; ++i){
		tbl[i] = new bool[n+2];
	}
	for (int i=0; i < n+2; ++i){
		for (int j=0; j < i; ++j){
			tbl[i][j] = false;
			tbl[j][i] = false;
		}
		tbl[i][i] = true;
	}
}

void ConTable::on_connect(int i, int j){
	tbl[i][j] = true;
	tbl[j][i] = true;
}

void ConTable::off_connect(int i, int j){
	tbl[i][j] = false;
	tbl[j][i] = false;
}

void ConTable::reset(){
	for (int i=0; i<n+2; ++i){
		for (int j=0; j<n+2; ++j){
			tbl[i][j] = false;
		}
		tbl[i][i] = true;
	}
}
