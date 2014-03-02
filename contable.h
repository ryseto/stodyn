/*
 *  contable.h
 *  small_system
 *
 *  Created by SETO Ryohei on 07/01/18.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef contable_h
#define contable_h 1
#include <cmath>
#include <cstdlib>
#include "my_utilities.h"
//#define DELETE(x) if(x){delete [] x; x = NULL;}
using namespace std;

class ConTable{
	bool allocate;
	bool **tbl;
	int n;
public:
	ConTable():allocate(false) {}
	~ConTable();
	inline bool connect(int i, int j){
		return tbl[i][j];
	}

	void set(int particleNumber);
	void reset();
	void on_connect(int i, int j);
	void off_connect(int i, int j);	
};
#endif
